#include "meraxes.h"
#include "parse_paramfile.h"

static void check_problem_params(run_params_t *run_params)
{
  if (run_params->NSteps != 1)
  {
    SID_log_error("The current version of the code only works if NSteps = 1. Sorry! Exiting...");
    ABORT(EXIT_FAILURE);
  }

  if (run_params->FlagReadDumpFile && run_params->FlagGenDumpFile)
  {
    SID_log_error("Both FlagReadDumpFile & FlagGenDumpFile are set!");
    ABORT(EXIT_FAILURE);
  }
}


static void inline store_params(
  entry_t       entry[PARAM_MAX_ENTRIES],
  int           n_entries,
  char        **params_tag,
  int           n_param,
  int           used_tag[PARAM_MAX_ENTRIES],
  int          *params_type,
  void        **params_addr,
  run_params_t *run_params)
{
  int level = 0;
  int tag_index;
  int temp;
  char prefix[16] = "\0";
  char key[STRLEN];

  for (int i_entry = 0; i_entry < n_entries; i_entry++)
  {
    // DEBUG
    // SID_log("Checking %s", SID_LOG_COMMENT, entry[i_entry].key);

    // reset prefix if we have descended an indentation level
    if (entry[i_entry].level < level)
      *prefix = '\0';

    strcpy(key, prefix);
    strcat(key, entry[i_entry].key);
    level = entry[i_entry].level;

    // DEBUG
    // SID_log("level = %d :: prefix = %s", SID_LOG_COMMENT, level, prefix);

    tag_index = -1;
    for (int ii = 0; ii < n_param; ii++)
    {
      if (strcmp(key, params_tag[ii]) == 0)
      {
        tag_index = ii;
        break;
      }
    }

    if (strcmp(key, "TOCF_Flag") == 0)
    {
      temp = atoi(entry[i_entry].value);
      if (used_tag[tag_index] == 0)
      {
        *((int*)params_addr[tag_index]) = atoi(entry[i_entry].value);
        used_tag[tag_index]             = 1;
      }
      if (temp != 1)
      {
        SID_log("Skipping TOCF params block...", SID_LOG_COMMENT);
        temp = entry[i_entry++].level;
        while (entry[i_entry++].level > temp)
          ;
        i_entry -= 2;
      }
      else
        sprintf(prefix, "TOCF_");
    }

    if (tag_index < 0)
    {
      SID_log_warning("%s is an unrecognised parameter (prefix='%s').", SID_LOG_COMMENT, entry[i_entry].key, prefix);
      ABORT(EXIT_FAILURE);
    }

    if (used_tag[tag_index] == 1)
      continue;

    switch (params_type[tag_index])
    {
    case PARAM_TYPE_DOUBLE:
      *((double*)params_addr[tag_index]) = atof(entry[i_entry].value);
      break;

    case PARAM_TYPE_FLOAT:
      *((float*)params_addr[tag_index]) = atof(entry[i_entry].value);
      break;

    case PARAM_TYPE_STRING:
      strcpy(params_addr[tag_index], entry[i_entry].value);
      break;

    case PARAM_TYPE_INT:
      *((int*)params_addr[tag_index]) = atoi(entry[i_entry].value);
      break;
    }
    used_tag[tag_index] = 1;
  }
}


void read_parameter_file(char *fname, int mode)
{
  // mode = 0 : for single runs
  // mode = 1 : for interactive runs (do not remalloc arrays)

  run_params_t *run_params = &(run_globals.params);

  if (SID.My_rank == 0)
  {
    int ii;
    int used_tag[PARAM_MAX_ENTRIES], required_tag[PARAM_MAX_ENTRIES];
    entry_t entry[PARAM_MAX_ENTRIES];

    int n_entries;
    hdf5_output_t *hdf5props = &(run_globals.hdf5props);
    int *params_type = hdf5props->params_type;
    void **params_addr = hdf5props->params_addr;
    char **params_tag = hdf5props->params_tag;
    int n_param = hdf5props->params_count;

    printf("\nreading parameter file:\n\n");

    // even if this is an interactive run, we have to init required_tag
    // (used_tag is initialised later)
    for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
      required_tag[ii] = 0;

    // malloc global arrays and init param properties
    if (mode == 0)
    {
      hdf5props->params_tag = SID_malloc(sizeof(char *) * PARAM_MAX_ENTRIES);
      for (int ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
        hdf5props->params_tag[ii] = SID_malloc(sizeof(char) * 128);
      params_tag = hdf5props->params_tag;
      hdf5props->params_type = SID_malloc(sizeof(int) * PARAM_MAX_ENTRIES);
      params_type = hdf5props->params_type;
      hdf5props->params_addr = SID_malloc(sizeof(void *) * PARAM_MAX_ENTRIES);
      params_addr = hdf5props->params_addr;

      // Initialise values and arrays
      n_param = 0;
      for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
      {
        required_tag[ii] = 0;
        hdf5props->params_type[ii] = PARAM_TYPE_UNUSED;
      }

      strcpy(params_tag[n_param], "DefaultsFile");
      params_addr[n_param]   = run_params->DefaultsFile;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "OutputDir");
      params_addr[n_param]   = run_params->OutputDir;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "PhotometricTablesDir");
      params_addr[n_param]   = run_params->PhotometricTablesDir;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "CoolingFuncsDir");
      params_addr[n_param]   = run_params->CoolingFuncsDir;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "SSPModel");
      params_addr[n_param]   = run_params->SSPModel;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "IMF");
      params_addr[n_param]   = run_params->IMF;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "MagSystem");
      params_addr[n_param]   = run_params->MagSystem;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "MagBands");
      params_addr[n_param]   = run_params->MagBands;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "FileNameGalaxies");
      params_addr[n_param]   = run_params->FileNameGalaxies;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "SimName");
      params_addr[n_param]   = run_params->SimName;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "SimulationDir");
      params_addr[n_param]   = run_params->SimulationDir;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "CatalogFilePrefix");
      params_addr[n_param]   = run_params->CatalogFilePrefix;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "FileWithOutputSnaps");
      params_addr[n_param]   = run_params->FileWithOutputSnaps;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "NSteps");
      params_addr[n_param]   = &(run_params->NSteps);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "BoxSize");
      params_addr[n_param]   = &(run_params->BoxSize);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "VolumeFactor");
      params_addr[n_param]   = &(run_params->VolumeFactor);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "RandomSeed");
      params_addr[n_param]   = &(run_params->RandomSeed);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "ForestIDFile");
      params_addr[n_param]        = &(run_params->ForestIDFile);
      required_tag[n_param]       = 0;
      params_type[n_param++]      = PARAM_TYPE_STRING;
      *(run_params->ForestIDFile) = '\0';

      strcpy(params_tag[n_param], "MvirCritFile");
      params_addr[n_param]        = &(run_params->MvirCritFile);
      required_tag[n_param]       = 0;
      params_type[n_param++]      = PARAM_TYPE_STRING;
      *(run_params->MvirCritFile) = '\0';

      strcpy(params_tag[n_param], "UnitVelocity_in_cm_per_s");
      params_addr[n_param]   = &(run_globals.units.UnitVelocity_in_cm_per_s);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "UnitLength_in_cm");
      params_addr[n_param]   = &(run_globals.units.UnitLength_in_cm);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "UnitMass_in_g");
      params_addr[n_param]   = &(run_globals.units.UnitMass_in_g);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "Hubble_h");
      params_addr[n_param]   = &(run_params->Hubble_h);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "BaryonFrac");
      params_addr[n_param]   = &(run_params->BaryonFrac);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "OmegaM");
      params_addr[n_param]   = &(run_params->OmegaM);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "OmegaK");
      params_addr[n_param]   = &(run_params->OmegaK);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "OmegaR");
      params_addr[n_param]   = &(run_params->OmegaR);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "OmegaLambda");
      params_addr[n_param]   = &(run_params->OmegaLambda);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "Sigma8");
      params_addr[n_param]   = &(run_params->Sigma8);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "wLambda");
      params_addr[n_param]   = &(run_params->wLambda);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SpectralIndex");
      params_addr[n_param]   = &(run_params->SpectralIndex);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "PartMass");
      params_addr[n_param]   = &(run_params->PartMass);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MergerTimeFactor");
      params_addr[n_param]   = &(run_params->physics.MergerTimeFactor);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "FlagSubhaloVirialProps");
      params_addr[n_param]   = &(run_params->FlagSubhaloVirialProps);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "FlagInteractive");
      params_addr[n_param]   = &(run_params->FlagInteractive);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "FlagGenDumpFile");
      params_addr[n_param]   = &(run_params->FlagGenDumpFile);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "FlagReadDumpFile");
      params_addr[n_param]   = &(run_params->FlagReadDumpFile);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;


      // Physics params
      
      strcpy(params_tag[n_param], "Flag_RedshiftDepEscFrac");
      params_addr[n_param]   = &(run_params->physics).Flag_RedshiftDepEscFrac;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;
      
      strcpy(params_tag[n_param], "Flag_ReionizationModifier");
      params_addr[n_param]   = &(run_params->physics).Flag_ReionizationModifier;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "Flag_BHFeedback");
      params_addr[n_param]   = &(run_params->physics).Flag_BHFeedback;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "Flag_IRA");
      params_addr[n_param]   = &(run_params->physics).Flag_IRA;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "Flag_FixDiskRadiusOnInfall");
      params_addr[n_param]   = &(run_params->physics).Flag_FixDiskRadiusOnInfall;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "Flag_FixVmaxOnInfall");
      params_addr[n_param]   = &(run_params->physics).Flag_FixVmaxOnInfall;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "Flag_ReheatToFOFGroupTemp");
      params_addr[n_param]   = &(run_params->physics).Flag_ReheatToFOFGroupTemp;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "SfDiskVelOpt");
      params_addr[n_param]   = &(run_params->physics).SfDiskVelOpt;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "SfEfficiency");
      params_addr[n_param]   = &(run_params->physics).SfEfficiency;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SfCriticalSDNorm");
      params_addr[n_param]   = &(run_params->physics).SfCriticalSDNorm;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SfRecycleFraction");
      params_addr[n_param]   = &(run_params->physics).SfRecycleFraction;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnEjectionEff");
      params_addr[n_param]   = &(run_params->physics).SnEjectionEff;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnEjectionScaling");
      params_addr[n_param]   = &(run_params->physics).SnEjectionScaling;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnEjectionNorm");
      params_addr[n_param]   = &(run_params->physics).SnEjectionNorm;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnReheatEff");
      params_addr[n_param]   = &(run_params->physics).SnReheatEff;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnReheatLimit");
      params_addr[n_param]   = &(run_params->physics).SnReheatLimit;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnReheatScaling");
      params_addr[n_param]   = &(run_params->physics).SnReheatScaling;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SnReheatNorm");
      params_addr[n_param]   = &(run_params->physics).SnReheatNorm;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReincorporationEff");
      params_addr[n_param]   = &(run_params->physics).ReincorporationEff;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MaxCoolingMassFactor");
      params_addr[n_param]   = &(run_params->physics).MaxCoolingMassFactor;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "Yield");
      params_addr[n_param]   = &(run_params->physics).Yield;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "IMFSlope");
      params_addr[n_param]   = &(run_params->physics).IMFSlope;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "EnergyPerSN");
      params_addr[n_param]   = &(run_params->physics).EnergyPerSN;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "IMFNormConst");
      params_addr[n_param]   = &(run_params->physics).IMFNormConst;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ThreshMajorMerger");
      params_addr[n_param]   = &((run_params->physics).ThreshMajorMerger);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MinMergerStellarMass");
      params_addr[n_param]   = &((run_params->physics).MinMergerStellarMass);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MinMergerRatioForBurst");
      params_addr[n_param]   = &((run_params->physics).MinMergerRatioForBurst);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MergerBurstScaling");
      params_addr[n_param]   = &((run_params->physics).MergerBurstScaling);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "MergerBurstFactor");
      params_addr[n_param]   = &((run_params->physics).MergerBurstFactor);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "RadioModeEff");
      params_addr[n_param]   = &(run_params->physics).RadioModeEff;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "BlackHoleGrowthRate");
      params_addr[n_param]   = &(run_params->physics).BlackHoleGrowthRate;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionSobacchi_Zre");
      params_addr[n_param]   = &(run_params->physics).ReionSobacchi_Zre;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionSobacchi_DeltaZre");
      params_addr[n_param]   = &(run_params->physics).ReionSobacchi_DeltaZre;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionSobacchi_DeltaZsc");
      params_addr[n_param]   = &(run_params->physics).ReionSobacchi_DeltaZsc;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionSobacchi_T0");
      params_addr[n_param]   = &(run_params->physics).ReionSobacchi_T0;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionTcool");
      params_addr[n_param]   = &(run_params->physics).ReionTcool;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionNionPhotPerBary");
      params_addr[n_param]   = &(run_params->physics).ReionNionPhotPerBary;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionEscapeFrac");
      params_addr[n_param]   = &(run_params->physics).ReionEscapeFrac;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionGnedin_z0");
      params_addr[n_param]   = &(run_params->physics).ReionGnedin_z0;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionGnedin_zr");
      params_addr[n_param]   = &(run_params->physics).ReionGnedin_zr;
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "TOCF_Flag");
      params_addr[n_param]   = &(run_params->TOCF_Flag);
      required_tag[n_param]  = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "TOCF_HII_dim");
      params_addr[n_param]   = &(tocf_params.HII_dim);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "TOCF_ion_tvir_min");
      params_addr[n_param]   = &(tocf_params.ion_tvir_min);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "TOCF_HII_eff_factor");
      params_addr[n_param]   = &(tocf_params.HII_eff_factor);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_r_bubble_min");
      params_addr[n_param]   = &(tocf_params.r_bubble_min);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_r_bubble_max");
      params_addr[n_param]   = &(tocf_params.r_bubble_max);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_delta_r_HII_factor");
      params_addr[n_param]   = &(tocf_params.delta_r_HII_factor);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_gamma_halo_bias");
      params_addr[n_param]   = &(tocf_params.gamma_halo_bias);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_m_0_sm");
      params_addr[n_param]   = &(tocf_params.m_0_sm);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_a_sm");
      params_addr[n_param]   = &(tocf_params.a_sm);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_b_sm");
      params_addr[n_param]   = &(tocf_params.b_sm);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_c_sm");
      params_addr[n_param]   = &(tocf_params.c_sm);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_d_sm");
      params_addr[n_param]   = &(tocf_params.d_sm);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;

      strcpy(params_tag[n_param], "TOCF_uvb_feedback");
      params_addr[n_param]   = &(tocf_params.uvb_feedback);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "TOCF_compute_mfp");
      params_addr[n_param]   = &(tocf_params.compute_mfp);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "TOCF_HII_filter");
      params_addr[n_param]   = &(tocf_params.HII_filter);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_INT;

      strcpy(params_tag[n_param], "TOCF_delta_k_ps");
      params_addr[n_param]   = &(tocf_params.delta_k_ps);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;
      
      strcpy(params_tag[n_param], "TOCF_alpha_uv");
      params_addr[n_param]   = &(tocf_params.alpha_uv);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;
      
      strcpy(params_tag[n_param], "TOCF_RtoM_filter");
      params_addr[n_param]   = &(tocf_params.RtoM_filter);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_INT;
      
      strcpy(params_tag[n_param], "TOCF_Y_He");
      params_addr[n_param]   = &(tocf_params.Y_He);
      required_tag[n_param]  = 0;
      params_type[n_param++] = PARAM_TYPE_FLOAT;
      
      hdf5props->params_count = n_param;
    }

    // Initialise used_tag
    for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
      used_tag[ii] = 0;

    // Parse the user parameter file first
    n_entries = parse_paramfile(fname, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr, run_params);

    // Now parse the default parameter file
    n_entries = parse_paramfile(run_params->DefaultsFile, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr, run_params);

    // Check to see if we are missing any required parameters
    for (ii = 0; ii < n_param; ii++)
    {
      if ((used_tag[ii] == 0) && (required_tag[ii] == 1))
      {
        SID_log_error("I miss a value for tag '%s' in parameter file '%s'.", params_tag[ii], fname);
        ABORT(EXIT_FAILURE);
      }
    }

    if (SID.My_rank == 0)
    {
      for (ii = 0; ii < n_param; ii++)
      {
        if (used_tag[ii] == 1)
        {
          printf("%35s\t", params_tag[ii]);

          switch (params_type[ii])
          {
          case PARAM_TYPE_DOUBLE:
            printf("%g\n", *((double*)(params_addr[ii])));
            break;

          case PARAM_TYPE_FLOAT:
            printf("%g\n", *((float*)(params_addr[ii])));
            break;

          case PARAM_TYPE_STRING:
            printf("%s\n", (char*)params_addr[ii]);
            break;

          case PARAM_TYPE_INT:
            printf("%d\n", *((int*)(params_addr[ii])));
            break;
          }

          // if (ii==(n_param-12))
          //   printf("\t\t%35s\n", "--- physics parameters ---");

          // if (ii==(n_param-4))
          //   printf("\t\t%35s\n", "--- TOCF parameters ---");
        }
      }
    }

    ii = strlen(run_params->OutputDir);
    if (ii > 0)
      if (run_params->OutputDir[ii - 1] != '/')
        strcat(run_params->OutputDir, "/");

    check_problem_params(run_params);
  }  // END if(SID.My_rank==0)

  // If running mpi then broadcast the run parameters to all cores
  SID_Bcast(run_params, sizeof(run_params_t), 0, SID.COMM_WORLD);
  SID_Bcast(&(run_globals.units), sizeof(run_units_t), 0, SID.COMM_WORLD);
  SID_Bcast(&tocf_params, sizeof(tocf_params_t), 0, SID.COMM_WORLD);
}
