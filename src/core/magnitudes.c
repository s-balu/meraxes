#ifdef CALC_MAGS
#include "meraxes.h"

enum core {MASTER};

void init_magnitudes(void) {
    int mpi_rank = run_globals.mpi_rank;
    mini_sed_params_t *mags_params = &run_globals.mags_params;

    // Initalise all relevant parameters at the master core
    if (mpi_rank == MASTER) {
        printf("#***********************************************************\n");
        printf("# Compute magnitudes\n");

        // Read target snapshots
        run_params_t *params = &run_globals.params;
        char str[STRLEN];
        char delim[] = ",";
        char *token;
        int target_snaps[MAGS_N_SNAPS];

        memcpy(str, params->TargetSnaps, sizeof(str));
        token = strtok(str, delim);
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap) {
            if (token != NULL) {
                target_snaps[i_snap] = atoi(token);
                token = strtok(NULL, delim);
            }
            else if (i_snap != MAGS_N_SNAPS - 1) {
                mlog_error("TargetSnaps does not match MAGS_N_SNAPS!");
                ABORT(EXIT_FAILURE);
            }
        }
        printf("# Target snapshots: ");
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap)
            printf("%d ", target_snaps[i_snap]);
        printf("\n");
        // Read beta filters
        double beta_bands[2*MAGS_N_BANDS];
        int n_beta = 0;

        memcpy(str, params->BetaBands, sizeof(str));
        token = strtok(str, delim);
        for(int i_band = 0; i_band < 2*MAGS_N_BANDS; ++i_band) {
            if (token != NULL) {
                beta_bands[i_band] = atof(token);
                token = strtok(NULL, delim);
                ++n_beta;
            }
            else
                break;
        }
        if (n_beta%2 == 0)
            n_beta /= 2;
        else {
            mlog_error("Wrong BetaBands!");
            ABORT(EXIT_FAILURE);
        }
        printf("# Beta filters:\n");
        for(int i_band = 0; i_band < n_beta; ++i_band)
            printf("#\t%.1f AA to %.1f\n", beta_bands[2*i_band], beta_bands[2*i_band + 1]);
        // Read rest-frame filters
        double rest_bands[2*MAGS_N_BANDS];
        int n_rest = 0;

        memcpy(str, params->RestBands, sizeof(str));
        token = strtok(str, delim);
        for(int i_band = 0; i_band < 2*MAGS_N_BANDS; ++i_band) {
            if (token != NULL) {
                rest_bands[i_band] = atof(token);
                token = strtok(NULL, delim);
                ++n_rest;
            }
            else
                break;
        }
        if (n_rest%2 == 0)
            n_rest /= 2;
        else {
            mlog_error("Wrong RestBands!");
            ABORT(EXIT_FAILURE);
        }
        printf("# Rest-frame filters:\n");
        for(int i_band = 0; i_band < n_rest; ++i_band)
            printf("#\t%.1f AA to %.1f\n", rest_bands[2*i_band], rest_bands[2*i_band + 1]);
        //
        if (n_beta + n_rest != MAGS_N_BANDS) {
            mlog_error("Number of beta and rest-frame filters do not match MAGS_N_BANDS!");
            ABORT(EXIT_FAILURE);
        }
        printf("#***********************************************************\n\n");

        // Initialise SED templates
        ////
        char *fname = params->PhotometricTablesDir;
        strcat(fname, "/sed_library.hdf5");
        ////Convert time unit to yr
        int snaplist_len = params->SnaplistLength;
        double *LTTime = malloc(snaplist_len*sizeof(double));
        double time_unit = run_globals.units.UnitTime_in_Megayears/params->Hubble_h*1e6;

        memcpy(LTTime, run_globals.LTTime, snaplist_len*sizeof(double));
        for(int i_time = 0; i_time < snaplist_len; ++i_time)
            LTTime[i_time] *= time_unit;
        ////
        init_templates_mini(mags_params, fname, LTTime, target_snaps, run_globals.ZZ,
                            beta_bands, n_beta, rest_bands, n_rest, params->BirthCloudLifetime);
    }

    // Broadcast parameters to all cores
    MPI_Comm mpi_comm = run_globals.mpi_comm;
    double *working;
    ptrdiff_t offset_inBC;
    ptrdiff_t offset_outBC;
    ptrdiff_t offset_waves;
    ptrdiff_t offset_logWaves;

    if (mpi_rank == MASTER) {
        working = mags_params->working;
        offset_inBC = mags_params->inBC - working;
        offset_outBC = mags_params->outBC - working;
        offset_waves = mags_params->centreWaves - working;
        offset_logWaves = mags_params->logWaves - working;

        mags_params->working = NULL;
        mags_params->inBC = NULL;
        mags_params->outBC = NULL;
        mags_params->centreWaves = NULL;
        mags_params->logWaves = NULL;
    }

    MPI_Bcast(mags_params, sizeof(mini_sed_params_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_inBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_outBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_waves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_logWaves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    if (mpi_rank != MASTER)
        working = (double*)malloc(mags_params->totalSize);
    MPI_Bcast(working, mags_params->totalSize, MPI_BYTE, MASTER, mpi_comm);

    mags_params->working = working;
    mags_params->inBC = working + offset_inBC;
    mags_params->outBC = working + offset_outBC;
    mags_params->centreWaves = working + offset_waves;
    mags_params->logWaves = working + offset_logWaves;
}

void cleanup_mags(void) {
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
    free(run_globals.mags_params.working);
}

void get_output_magnitudes(float *target, galaxy_t *gal, int snapshot) {
    // Check if ``snapshot`` is a target snapshot
    int iS;
    int *targetSnap = run_globals.mags_params.targetSnap;
    double *pInBCFlux = gal->inBCFlux;
    double *pOutBCFlux = gal->outBCFlux;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        if (snapshot == targetSnap[iS])
            break;
        else {
            pInBCFlux += MAGS_N_BANDS;
            pOutBCFlux += MAGS_N_BANDS;
        }
    }

    if (iS != MAGS_N_SNAPS) {
        // Convert fluxes to magnitudes
        //   -Convert the unit of SFRs
        double mags[MAGS_N_BANDS];
        double sfr_unit = -2.5*log10(run_globals.units.UnitMass_in_g
                                     /run_globals.units.UnitTime_in_s
                                     *SEC_PER_YEAR/SOLAR_MASS);
        get_magnitudes(mags, pInBCFlux, pOutBCFlux);
        for(int i_band = 0; i_band < MAGS_N_BANDS; ++i_band)
            target[i_band] = (float)(mags[i_band] + sfr_unit);
    }
}


#endif
