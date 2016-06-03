#include "meraxes.h"
#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>

void set_ReionEfficiency()
{
  if (run_globals.params.PatchyReionFlag)
  {
    // Use the params passed to Meraxes via the input file to set the HII ionising efficiency factor
    physics_params_t *params = &(run_globals.params.physics);

    // If we are using a redshift dependent escape fraction then reset
    // ReionEscapeFrac to one as we don't want to inlcude it in the
    // ReionEfficiency (it will be included in the stellar mass and SFR grids sent
    // to 21cmFAST instead).
    if (params->Flag_RedshiftDepEscFrac)
    {
      SID_log("Flag_RedshiftDepEscFrac is on => setting ReionEscapeFrac = 1.", SID_LOG_COMMENT);
      params->ReionEscapeFrac = 1.0;
    }

    // The following is based on Sobacchi & Messinger (2013) eqn 7
    // with f_* removed and f_b added since we define f_coll as M_*/M_tot rather than M_vir/M_tot,
    // and also with the inclusion of the effects of the Helium fraction.
    run_globals.params.physics.ReionEfficiency = 1.0 / run_globals.params.BaryonFrac
      * params->ReionNionPhotPerBary * params->ReionEscapeFrac / (1.0 - 0.75*run_globals.params.physics.Y_He);

    // Account for instantaneous recycling factor so that stellar mass is cumulative
    if (params->Flag_IRA)
      run_globals.params.physics.ReionEfficiency /= params->SfRecycleFraction;

    SID_log("Set value of run_globals.params.ReionEfficiency = %g", SID_LOG_COMMENT, run_globals.params.physics.ReionEfficiency);
  }
  else {
    run_globals.params.physics.ReionEfficiency = -1;
  }
}



void assign_slabs()
{
  // Allocations made in this function are free'd in `free_reionization_grids`.
  fftwf_mpi_init();

  // Assign the slab size
  int n_rank = SID.n_proc;
  int dim = run_globals.params.ReionGridDim;

  // Use fftw to find out what slab each rank should get
  ptrdiff_t local_nix, local_ix_start;
  ptrdiff_t local_n_complex = fftwf_mpi_local_size_3d(dim, dim, dim/2 + 1, SID_COMM_WORLD, &local_nix, &local_ix_start);

  // let every core know...
  ptrdiff_t **slab_nix = &run_globals.reion_grids.slab_nix;
  *slab_nix = SID_malloc(sizeof(ptrdiff_t) * n_rank);  ///< array of number of x cells of every rank
  MPI_Allgather(&local_nix, sizeof(ptrdiff_t), MPI_BYTE, *slab_nix, sizeof(ptrdiff_t), MPI_BYTE, SID_COMM_WORLD);

  ptrdiff_t **slab_ix_start = &run_globals.reion_grids.slab_ix_start;
  *slab_ix_start = SID_malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
  (*slab_ix_start)[0] = 0;
  for(int ii=1; ii<n_rank; ii++)
    (*slab_ix_start)[ii] = (*slab_ix_start)[ii-1] + (*slab_nix)[ii-1];

  ptrdiff_t **slab_n_complex = &run_globals.reion_grids.slab_n_complex;  ///< array of allocation counts for every rank
  *slab_n_complex = SID_malloc(sizeof(ptrdiff_t) * n_rank);  ///< array of allocation counts for every rank
  MPI_Allgather(&local_n_complex, sizeof(ptrdiff_t), MPI_BYTE, *slab_n_complex, sizeof(ptrdiff_t), MPI_BYTE, SID_COMM_WORLD);
}



void call_find_HII_bubbles(int snapshot, int unsampled_snapshot, int nout_gals)
{
  // Thin wrapper round find_HII_bubbles

  int total_n_out_gals = 0;

  reion_grids_t *grids = &(run_globals.reion_grids);

  SID_log("Getting ready to call find_HII_bubbles...", SID_LOG_OPEN);

  // Check to see if there are actually any galaxies at this snapshot
  SID_Allreduce(&nout_gals, &total_n_out_gals, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  if (total_n_out_gals == 0)
  {
    SID_log("No galaxies in the simulation - skipping...", SID_LOG_CLOSE);
    return;
  }

  // Construct the baryon grids
  construct_baryon_grids(snapshot, nout_gals);

  SID_log("...done", SID_LOG_CLOSE);

  // Read in the dark matter density grid
  read_dm_grid(unsampled_snapshot, 0, (float*)(grids->deltax));

  // Call find_HII_bubbles
  SID_log("Calling find_HII_bubbles...", SID_LOG_OPEN | SID_LOG_TIMER);
  grids->global_xH = find_HII_bubbles(run_globals.ZZ[snapshot]);

  SID_log("grids->global_xH = %g", SID_LOG_COMMENT, grids->global_xH);
  SID_log("...done", SID_LOG_CLOSE);
}


void malloc_reionization_grids()
{

  reion_grids_t *grids = &(run_globals.reion_grids);

  grids->galaxy_to_slab_map = NULL;

  grids->xH                 = NULL;
  grids->stars              = NULL;
  grids->stars_unfiltered   = NULL;
  grids->stars_filtered     = NULL;
  grids->deltax             = NULL;
  grids->deltax_unfiltered  = NULL;
  grids->deltax_filtered    = NULL;
  grids->sfr                = NULL;
  grids->sfr_unfiltered     = NULL;
  grids->sfr_filtered       = NULL;
  grids->z_at_ionization    = NULL;
  grids->J_21_at_ionization = NULL;
  grids->J_21               = NULL;

  grids->global_xH = 1.0;
  grids->reion_complete = false;

  if (run_globals.params.PatchyReionFlag)
  {

    assign_slabs();

    int ReionGridDim = run_globals.params.ReionGridDim;
    ptrdiff_t *slab_nix = run_globals.reion_grids.slab_nix;
    ptrdiff_t slab_n_real = slab_nix[SID.My_rank] * ReionGridDim * ReionGridDim; // TODO: NOT WORKING!!!
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[SID.My_rank];

    // create a buffer on each rank which is as large as the largest LOGICAL allocation on any single rank
    int max_cells = 0;

    for(int ii=0; ii < SID.n_proc; ii++)
      if(slab_nix[ii] > max_cells)
        max_cells = slab_nix[ii];

    max_cells *= ReionGridDim * ReionGridDim;
    grids->buffer_size = max_cells;

    grids->buffer          = fftwf_alloc_real(max_cells);
    grids->stars           = fftwf_alloc_real(slab_n_complex * 2);  // padded for in-place FFT
    grids->stars_filtered  = fftwf_alloc_complex(slab_n_complex);
    grids->deltax          = fftwf_alloc_real(slab_n_complex * 2);  // padded for in-place FFT
    grids->deltax_filtered = fftwf_alloc_complex(slab_n_complex);
    grids->sfr             = fftwf_alloc_real(slab_n_complex * 2);  // padded for in-place FFT
    grids->sfr_filtered    = fftwf_alloc_complex(slab_n_complex);
    grids->xH              = fftwf_alloc_real(slab_n_real);
    grids->z_at_ionization = fftwf_alloc_real(slab_n_real);

    if (run_globals.params.ReionUVBFlag)
    {
      grids->J_21_at_ionization = fftwf_alloc_real(slab_n_real);
      grids->J_21               = fftwf_alloc_real(slab_n_real);
      grids->Mvir_crit          = fftwf_alloc_real(slab_n_real);
    }

    SID_log("Initialising grids...", SID_LOG_COMMENT);

    for (int ii = 0; ii < slab_n_real; ii++)
    {
      grids->xH[ii] = 1.0;
      grids->z_at_ionization[ii] = -1;
    }


    for (int ii = 0; ii < slab_n_real; ii++)
      if (run_globals.params.ReionUVBFlag)
      {
        grids->J_21_at_ionization[ii] = 0.;
        grids->J_21[ii] = 0.;
        grids->Mvir_crit[ii] = 0;
      }


    for (int ii = 0; ii < slab_n_complex; ii++)
    {
      grids->stars_filtered[ii] = 0 + 0*I;
      grids->deltax_filtered[ii] = 0 + 0*I;
      grids->sfr_filtered[ii] = 0 + 0*I;
    }

    for (int ii = 0; ii < slab_n_complex*2; ii++)
    {
      grids->deltax[ii] = 0;
      grids->stars[ii] = 0;
      grids->sfr[ii] = 0;
    }

    SID_log(" ...done", SID_LOG_CLOSE);
  }
}


void free_reionization_grids()
{
  SID_log("Freeing reionization grids...", SID_LOG_OPEN);

  reion_grids_t *grids = &(run_globals.reion_grids);

  SID_free(SID_FARG run_globals.reion_grids.slab_n_complex);
  SID_free(SID_FARG run_globals.reion_grids.slab_ix_start);
  SID_free(SID_FARG run_globals.reion_grids.slab_nix);

  if (run_globals.params.ReionUVBFlag)
  {
    fftwf_free(grids->J_21);
    fftwf_free(grids->J_21_at_ionization);
  }
  fftwf_free(grids->z_at_ionization);
  fftwf_free(grids->sfr_filtered);
  fftwf_free(grids->deltax_filtered);
  fftwf_free(grids->deltax);
  fftwf_free(grids->stars_filtered);
  fftwf_free(grids->xH);

  if (run_globals.params.ReionUVBFlag)
    fftwf_free(grids->Mvir_crit);

  fftwf_free(grids->stars);
  fftwf_free(grids->sfr);
  fftwf_free(grids->buffer);

  SID_log(" ...done", SID_LOG_CLOSE);
}


int map_galaxies_to_slabs(int ngals)
{
  double box_size     = (double)(run_globals.params.BoxSize);
  int ReionGridDim         = run_globals.params.ReionGridDim;

  SID_log("Mapping galaxies to slabs...", SID_LOG_OPEN);

  // Loop through each valid galaxy and find what slab it sits in
  if (ngals > 0)
    run_globals.reion_grids.galaxy_to_slab_map = SID_malloc(sizeof(gal_to_slab_t) * ngals);
  else
    run_globals.reion_grids.galaxy_to_slab_map = NULL;

  gal_to_slab_t *galaxy_to_slab_map         = run_globals.reion_grids.galaxy_to_slab_map;
  ptrdiff_t *slab_ix_start = run_globals.reion_grids.slab_ix_start;

  galaxy_t *gal = run_globals.FirstGal;
  int gal_counter = 0;
  while (gal != NULL)
  {
    // TODO: Note that I am including ghosts here.  We will need to check the
    // validity of this.  By definition, if they are ghosts then their host
    // halo hasn't been identified at this time step and hence they haven't
    // been evolved.  Their properties (Sfr, StellarMass, etc.) will all have
    // been set when they were last identified.
    if (gal->Type < 3)
    {
      // TODO: for type 2 galaxies these positions will be set from the last
      // time they were identified.  If their host halo has moved significantly
      // since then, these positions won't reflect that and the satellites will
      // be spatially disconnected from their hosts.  We will need to fix this
      // at some point.

      // TODO: Get Greg to fix these positions to obey PBC!!
      if (gal->Pos[0] >= box_size)
        gal->Pos[0] -= box_size;
      else if (gal->Pos[0] < 0.0)
        gal->Pos[0] += box_size;
      ptrdiff_t ix = pos_to_cell(gal->Pos[0], box_size, ReionGridDim);

      assert((ix >= 0) && (ix < ReionGridDim));

      galaxy_to_slab_map[gal_counter].index = gal_counter;
      galaxy_to_slab_map[gal_counter].slab_ind = searchsorted(&ix, slab_ix_start, SID.n_proc, sizeof(ptrdiff_t), compare_ptrdiff, -1, -1);
      galaxy_to_slab_map[gal_counter++].galaxy = gal;
    }

    gal = gal->Next;
  }

  // sort the slab indices IN PLACE (n.b. compare_slab_assign is a stable comparison)
  qsort(galaxy_to_slab_map, gal_counter, sizeof(gal_to_slab_t), compare_slab_assign);

  assert(gal_counter == ngals);

  SID_log("...done.", SID_LOG_CLOSE);

  return gal_counter;
}


void assign_Mvir_crit_to_galaxies(int ngals_in_slabs)
{

  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  gal_to_slab_t *galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
  float         *Mvir_crit          = run_globals.reion_grids.Mvir_crit;
  float         *buffer             = run_globals.reion_grids.buffer;
  ptrdiff_t     *slab_nix           = run_globals.reion_grids.slab_nix;
  ptrdiff_t     *slab_ix_start      = run_globals.reion_grids.slab_ix_start;
  int           ReionGridDim             = run_globals.params.ReionGridDim;
  double        box_size            = run_globals.params.BoxSize;

  SID_log("Assigning Mvir_crit to galaxies...", SID_LOG_OPEN);

  // Work out the index of the galaxy_to_slab_map where each slab begins.
  int slab_map_offsets[SID.n_proc];
  for(int ii=0, i_gal=0; ii < SID.n_proc; ii++)
  {
    if (galaxy_to_slab_map != NULL)
    {
      while((galaxy_to_slab_map[i_gal].slab_ind < ii) && (i_gal < ngals_in_slabs))
        i_gal++;

      if(galaxy_to_slab_map[i_gal].slab_ind == ii)
        slab_map_offsets[ii] = i_gal;
      else
        slab_map_offsets[ii] = -1;
    }
    else
    {
      // if this core has no galaxies then the offsets are -1 everywhere
      slab_map_offsets[ii] = -1;
    }
  }

  // do a ring exchange of slabs between all cores
  for(int i_skip=0; i_skip < SID.n_proc; i_skip++)
  {
    int recv_from_rank = (SID.My_rank + i_skip) % SID.n_proc;
    int send_to_rank   = (SID.My_rank - i_skip + SID.n_proc) % SID.n_proc;

    bool send_flag     = false;
    bool recv_flag     = (slab_map_offsets[recv_from_rank] > -1);
      
    if (i_skip > 0)
    {
      SID_Sendrecv(&recv_flag, sizeof(bool), SID_BYTE, recv_from_rank, 6393762,
          &send_flag, sizeof(bool), SID_BYTE, send_to_rank, 6393762, SID.COMM_WORLD);

      // need to ensure sends and receives do not clash!
      if (send_to_rank > SID.My_rank)
      {
        if(send_flag)
        {
          int n_cells = slab_nix[SID.My_rank]*ReionGridDim*ReionGridDim;
          SID_Send(Mvir_crit, n_cells, SID_FLOAT, send_to_rank, 793710, SID.COMM_WORLD);
        }
        if(recv_flag)
        {
          int n_cells = slab_nix[recv_from_rank]*ReionGridDim*ReionGridDim; 
          SID_Recv(buffer, n_cells, SID_FLOAT, recv_from_rank, 793710, SID.COMM_WORLD);
        }
      }
      else
      {
        if(recv_flag)
        {
          int n_cells = slab_nix[recv_from_rank]*ReionGridDim*ReionGridDim; 
          SID_Recv(buffer, n_cells, SID_FLOAT, recv_from_rank, 793710, SID.COMM_WORLD);
        }
        if(send_flag)
        {
          int n_cells = slab_nix[SID.My_rank]*ReionGridDim*ReionGridDim;
          SID_Send(Mvir_crit, n_cells, SID_FLOAT, send_to_rank, 793710, SID.COMM_WORLD);
        }
      }
    }
    else
    {
      int n_cells = slab_nix[recv_from_rank]*ReionGridDim*ReionGridDim; 
      memcpy(buffer, Mvir_crit, sizeof(float) * n_cells);
    }

    // if this core has received a slab of Mvir_crit then assign values to the
    // galaxies which belong to this slab
    if(recv_flag)
    {
      int i_gal = slab_map_offsets[recv_from_rank];
      int ix_start = slab_ix_start[recv_from_rank];
      while((galaxy_to_slab_map[i_gal].slab_ind == recv_from_rank) && (i_gal < ngals_in_slabs))
      {
        // TODO: We should use the position of the FOF group here...
        galaxy_t *gal = galaxy_to_slab_map[i_gal].galaxy;
        int ix = pos_to_cell(gal->Pos[0], box_size, ReionGridDim) - ix_start;
        int iy = pos_to_cell(gal->Pos[1], box_size, ReionGridDim);
        int iz = pos_to_cell(gal->Pos[2], box_size, ReionGridDim);

        // Record the Mvir_crit (filtering mass) value
        gal->MvirCrit = (double)buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];

        // increment counter
        i_gal++;
      }
    }

  }

  SID_log("...done.", SID_LOG_CLOSE);

}


void construct_baryon_grids(int snapshot, int local_ngals)
{
  
  double box_size     = (double)(run_globals.params.BoxSize);
  float *stellar_grid = run_globals.reion_grids.stars;
  float *sfr_grid     = run_globals.reion_grids.sfr;
  int ReionGridDim         = run_globals.params.ReionGridDim;
  double tHubble      = hubble_time(snapshot);

  gal_to_slab_t *galaxy_to_slab_map         = run_globals.reion_grids.galaxy_to_slab_map;
  ptrdiff_t *slab_ix_start = run_globals.reion_grids.slab_ix_start;
  int local_n_complex = (int)(run_globals.reion_grids.slab_n_complex[SID.My_rank]);

  SID_log("Constructing stellar mass and sfr grids...", SID_LOG_OPEN | SID_LOG_TIMER);

  // init the grid
  for (int ii = 0; ii < local_n_complex; ii++)
  {
    *(stellar_grid + ii) = 0.0;
    *(sfr_grid + ii)     = 0.0;
  }

  // loop through each slab
  //
  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  double cell_width = box_size / (double)ReionGridDim;
  ptrdiff_t *slab_nix = run_globals.reion_grids.slab_nix;
  ptrdiff_t buffer_size = run_globals.reion_grids.buffer_size;
  float *buffer = run_globals.reion_grids.buffer;

  enum property { prop_stellar, prop_sfr };
  for(int prop = prop_stellar; prop <= prop_sfr; prop++)
  {
    int i_gal = 0;
    int skipped_gals = 0;

    for(int i_r=0; i_r < SID.n_proc; i_r++)
    {
      double min_xpos = (double)slab_ix_start[i_r] * cell_width;

      // init the buffer
      for(int ii=0; ii<buffer_size; ii++)
        buffer[ii] = 0.;

      // if this core holds no galaxies then we don't need to fill the buffer
      if(local_ngals != 0)
      {
        // fill the local buffer for this slab
        while(((i_gal - skipped_gals) < local_ngals) && (galaxy_to_slab_map[i_gal].slab_ind == i_r))
        {
          galaxy_t *gal = galaxy_to_slab_map[i_gal].galaxy;

          // Dead galaxies should not be included here and are not in the
          // local_ngals count.  They will, however, have been assigned to a
          // slab so we will need to ignore them here...
          if (gal->Type > 2)
          {
            i_gal++;
            skipped_gals++;
            continue;
          }

          assert(galaxy_to_slab_map[i_gal].index >= 0);
          assert((galaxy_to_slab_map[i_gal].slab_ind >= 0) && (galaxy_to_slab_map[i_gal].slab_ind < SID.n_proc));

          if (gal->Pos[0] >= box_size)
            gal->Pos[0] -= box_size;
          else if (gal->Pos[0] < 0.0)
            gal->Pos[0] += box_size;
          int ix = pos_to_cell(gal->Pos[0] - min_xpos, box_size, ReionGridDim);
          if (gal->Pos[1] >= box_size)
            gal->Pos[1] -= box_size;
          else if (gal->Pos[1] < 0.0)
            gal->Pos[1] += box_size;
          int iy = pos_to_cell(gal->Pos[1], box_size, ReionGridDim);
          if (gal->Pos[2] >= box_size)
            gal->Pos[2] -= box_size;
          else if (gal->Pos[2] < 0.0)
            gal->Pos[2] += box_size;
          int iz = pos_to_cell(gal->Pos[2], box_size, ReionGridDim);

          assert((ix < slab_nix[i_r]) && (ix >= 0));
          assert((iy < ReionGridDim) && (iy >= 0));
          assert((iz < ReionGridDim) && (iz >= 0));

          int ind = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

          assert((ind >=0) && (ind < slab_nix[i_r]*ReionGridDim*ReionGridDim));

          switch (prop) {
            case prop_stellar:
              buffer[ind] += gal->GrossStellarMass;
              break;

            case prop_sfr:
              buffer[ind] += gal->FescWeightedGSM;
              break;

            default:
              SID_log_error("Unrecognised property in slab creation.");
              ABORT(EXIT_FAILURE);
              break;
          }

          i_gal++;
        }
      }

      // reduce on to the correct rank
      if(SID.My_rank == i_r)
        SID_Reduce(MPI_IN_PLACE, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, SID.COMM_WORLD);
      else
        SID_Reduce(buffer, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, SID.COMM_WORLD);

      if (SID.My_rank == i_r)
      {
        // copy the buffer into the real slab
        float *slab;
        switch (prop)
        {
          case prop_stellar:
            slab = stellar_grid;
            break;
          case prop_sfr:
            slab = sfr_grid;
            break;
        }

        int slab_size = slab_nix[i_r] * ReionGridDim * ReionGridDim;
        memcpy(slab, buffer, sizeof(float)*slab_size);

				// Do one final pass and divide the sfr_grid by tHubble
        // in order to convert the stellar masses recorded into SFRs.
				switch (prop) {
					case prop_sfr:
						for(int ii=0; ii < slab_size; ii++)
						{
							if (slab[ii] > 0)
								slab[ii] /= tHubble;
							else if (slab[ii] < 0)
								slab[ii] = 0;
						}
						break;

					case prop_stellar:
						for(int ii=0; ii < slab_size; ii++)
							if (slab[ii] < 0)
                slab[ii] = 0;
            break;
				}

			}
    }
  }

  SID_log("done", SID_LOG_CLOSE);
}


static void write_grid_float(const char *name, float *data, hid_t file_id, hid_t fspace_id, hid_t memspace_id)
{
  hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, fspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id, fspace_id, plist_id, data);
  H5Pclose(plist_id);
  H5Dclose(dset_id);
}

void save_reion_grids(int snapshot)
{

  // Check if we even want to write anything...
  reion_grids_t *grids = &(run_globals.reion_grids);
  if ((grids->global_xH < REL_TOL) || (grids->global_xH > 1.0-REL_TOL))
    return;

  // If we do, well then lets!...
  int   ReionGridDim       = run_globals.params.ReionGridDim;
  int   local_nix     = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]);
  // float *ps;
  // int   ps_nbins;
  // float average_deltaT;
  // double Hubble_h = run_globals.params.Hubble_h;

  // Save tocf grids
  // ----------------------------------------------------------------------------------------------------

  SID_log("Saving tocf grids...", SID_LOG_OPEN);

  char name[STRLEN];
  sprintf(name, "%s/%s_grids.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies);

  // open the file (in parallel)
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, SID_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(name, H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);

  // create the group
  sprintf(name, "Snap%03d", snapshot);
  hid_t group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // create the filespace
  hsize_t dims[3] = {ReionGridDim, ReionGridDim, ReionGridDim};
  hid_t fspace_id = H5Screate_simple(3, dims, NULL);

  // create the memspace
  hsize_t mem_dims[3] = {local_nix, ReionGridDim, ReionGridDim};
  hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // select a hyperslab in the filespace
  hsize_t start[3] = {run_globals.reion_grids.slab_ix_start[SID.My_rank], 0, 0};
  hsize_t count[3] = {local_nix, ReionGridDim, ReionGridDim};
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // create and write the datasets
  write_grid_float("xH", grids->xH, group_id, fspace_id, memspace_id);
  write_grid_float("z_at_ionization", grids->z_at_ionization, group_id, fspace_id, memspace_id);

  if (run_globals.params.ReionUVBFlag)
  {
    write_grid_float("J_21", grids->J_21, group_id, fspace_id, memspace_id);
    write_grid_float("J_21_at_ionization", grids->J_21_at_ionization, group_id, fspace_id, memspace_id);
    write_grid_float("Mvir_crit", grids->Mvir_crit, group_id, fspace_id, memspace_id);
  }

  H5LTset_attribute_float(group_id, "xH", "global_xH", &(grids->global_xH), 1);

  // Save the escape fraction if we are using a redshift dependent escape fraction
  H5LTset_attribute_double(group_id, ".", "ReionEscapeFrac", &(run_globals.params.physics.ReionEscapeFrac), 1);

  // fftw padded grids
  float *grid = (float*)SID_calloc((int)local_nix * ReionGridDim * ReionGridDim * sizeof(float));
  for (int ii = 0; ii < local_nix; ii++)
    for (int jj = 0; jj < ReionGridDim; jj++)
      for (int kk = 0; kk < ReionGridDim; kk++)
        grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = ((float*)(grids->deltax))[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];
  write_grid_float("deltax", grids->deltax, group_id, fspace_id, memspace_id);


  // // Run delta_T_ps
  // // ----------------------------------------------------------------------------------------------------

  // SID_log("Calculating delta_T box and power spectrum...", SID_LOG_OPEN);

  // memset((void*)grid, 0, sizeof(float) * (int)dims);

  // delta_T_ps(
  //     run_globals.ZZ[snapshot],
  //     run_globals.params.numcores,
  //     grids->xH,
  //     (float*)(grids->deltax),
  //     &average_deltaT,
  //     grid,
  //     &ps,
  //     &ps_nbins);

  // H5LTmake_dataset_float(group_id, "delta_T", 1, &dims, grid);

  // dims = ps_nbins * 3;
  // H5LTmake_dataset_float(parent_group_id , "PowerSpectrum", 1               , &dims          , ps);
  // H5LTset_attribute_int(parent_group_id  , "PowerSpectrum", "nbins"         , &ps_nbins      , 1);
  // H5LTset_attribute_float(parent_group_id, "PowerSpectrum", "average_deltaT", &average_deltaT, 1);

  // free(ps);

  // SID_log("...done", SID_LOG_CLOSE);   // delta_T

  // tidy up
  SID_free(SID_FARG grid);
  H5Sclose(memspace_id);
  H5Sclose(fspace_id);
  H5Gclose(group_id);
  H5Fclose(file_id);

  SID_log("...done", SID_LOG_CLOSE);   // Saving tocf grids
}


bool check_if_reionization_complete()
{
  int complete = (int)run_globals.reion_grids.reion_complete;
  if(!complete)
  {
    complete = 1;
    float *xH = run_globals.reion_grids.xH;
    int ReionGridDim = run_globals.params.ReionGridDim;
    int slab_n_real = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]) * ReionGridDim * ReionGridDim;

    // If not all cells are ionised then reionization is still progressing...
    for (int ii=0; ii < slab_n_real; ii++)
    {
      if (xH[ii] != 0.0)
      {
        complete = 0;
        break;
      }
    }
  }
  
  SID_Allreduce(SID_IN_PLACE, &complete, 1, MPI_INT, MPI_LAND, SID.COMM_WORLD);
  run_globals.reion_grids.reion_complete = (bool)complete;
  return (bool)complete;
}
