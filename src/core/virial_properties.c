#include <math.h>

#include "meraxes.h"
#include "virial_properties.h"
#include <gsl/gsl_integration.h>

static float x_int_zvals[x_int_NCFVALS];
static float x_int_radvals[x_int_NCFVALS];
static float x_int_rvirvals[x_int_NCFVALS];
static float x_int_Sigmavals[x_int_NCFVALS];
static float x_int_CFvals[x_int_NCFVALS];

static inline double E_z(double z, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;

  one_plus_z = 1. + z;
  one_plus_z_sq = one_plus_z * one_plus_z;
  one_plus_z_cu = one_plus_z_sq * one_plus_z;
  result = sqrt(OmegaM * one_plus_z_cu + OmegaK * one_plus_z_sq + OmegaLambda);

  return result;
}

static inline double Omega_z(double redshift, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double Ez;
  double one_plus_z_cube;

  Ez = E_z(redshift, OmegaM, OmegaK, OmegaLambda);
  one_plus_z_cube = (1. + redshift) * (1. + redshift) * (1. + redshift);

  return OmegaM * one_plus_z_cube / (Ez * Ez);
}

static inline double Delta_vir(double redshift)
{
  // Function stolen and adapted from gbpCosmo
  double x;
  double Omega;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  Omega = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  x = Omega - 1.;

  return (18. * M_PI * M_PI + 82 * x - 39 * x * x) / Omega;
}

//! Calculates Mvir in internal units (1.e10 h^{-1}Msol), given Tvir (in K) and a redshift (z)
double Tvir_to_Mvir(double T, double z)
{
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  double mu; //!< Mean molecular weight (ionized gas)

  if (T < 9.99999e3) // Neutral IGM
    mu = 1.22;
  else // Ionised IGM
    mu = 0.59;

  double z_term = pow((1. + z) / 10., -1.5);
  double T_term = pow(T / 1.98e4, 1.5);
  double cosmo_term = pow(OmegaM / Omega_z(z, OmegaM, OmegaK, OmegaLambda) * Delta_vir(z) / 18. / (M_PI * M_PI), -0.5);
  double mol_term = pow(mu / 0.6, -1.5);

  return 0.01 * mol_term * cosmo_term * T_term * z_term;
}

double calculate_Mvir(double Mvir, int len)
{
  if ((len < 0) && (Mvir > 0))
    return Mvir;
  else
    return (double)len * run_globals.params.PartMass;
}

double hubble_at_snapshot(int snapshot)
{
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = run_globals.ZZ[snapshot] + 1;

  return Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda);
}

double hubble_time(int snapshot)
{
  return 1.0 / hubble_at_snapshot(snapshot);
}

double calculate_Rvir(double Mvir, int snapshot)
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;

  hubble_of_z_sq = pow(hubble_at_snapshot(snapshot), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  Delta = Delta_vir(run_globals.ZZ[snapshot]);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

  return cbrt(Mvir * fac);
}

/*double calculate_Rvir_2(double Mvir, double redshift) //from Mvir in 10^10 Msol/h to Rvir in comoving Mpc/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;

  hubble_of_z_sq = pow(Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  Delta = Delta_vir(redshift);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);
  //fac = 1 / (4 * M_PI / 3.0 * OmegaM * rhocrit);

  return cbrt(Mvir * fac) * zplus1; 
  //return cbrt(Mvir * fac);
}*/

double calculate_Rvir_2(double Mvir, double redshift) //from Mvir in 10^10 Msol/h to Rvir in comoving Mpc/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double little_h = run_globals.params.Hubble_h;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;
  double GG = 4.3009*1e-9;

  //rhocrit = 3 * (little_h * 100 * little_h * 100 * OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda) / (8 * M_PI * GG);
  rhocrit = 3 * little_h * little_h * 100 * 100 / (8 * M_PI * GG);

  //Delta = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  //fac = 1 / (OmegaM*zplus1*zplus1*zplus1 * 4 * M_PI / 3.0 * rhocrit);
  fac = 1 / (4 * M_PI / 3.0 * OmegaM * rhocrit);

  return cbrt(Mvir * 1e10 * fac) * little_h; 
  //return cbrt(Mvir * fac);
}

/*double calculate_Mvir_2(double Rvir, double redshift) //from Rvir in comoving Mpc/h to 10^10Msol/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;
  
  hubble_of_z_sq = pow(Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);
  
  Delta = Delta_vir(redshift);

  fac = 4.0 / 3.0 * M_PI * rhocrit * Delta;
  //fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM * Delta;

  return pow(Rvir / zplus1 , 3) * fac;
  //return pow(Rvir , 3) * fac;
}*/

double calculate_Mvir_2(double Rvir, double redshift) //from Rvir in comoving Mpc/h to 10^10Msol/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;
  double GG = 4.3009*1e-9;
  double little_h = run_globals.params.Hubble_h;
  
  //Delta = Delta_vir(redshift);
  //Delta = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);

  //rhocrit = 3 * (little_h * 100 * little_h * 100 * OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda) / (8 * M_PI * GG);
  rhocrit = 3 * little_h * little_h * 100 * 100 / (8 * M_PI * GG);

  //fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM * zplus1 * zplus1 * zplus1;
  fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM;

  return (pow(Rvir , 3) * fac) / 1e10;
  //return pow(Rvir , 3) * fac;
}

double calculate_gasMass(int snapshot, double length) //length in comoving units
{
  double hubble_of_z_sq;
  double rhocrit;
  double rhob;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaB = OmegaM * run_globals.params.BaryonFrac;
  
  hubble_of_z_sq = pow(hubble_at_snapshot(snapshot), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);
  rhob = rhocrit * OmegaB;

  return rhob * pow((length / (1.0 + run_globals.ZZ[snapshot])), 3.0);
}

double calculate_Vvir(double Mvir, double Rvir)
{
  return sqrt((run_globals.G) * Mvir / Rvir);
}

double calculate_spin_param(halo_t* halo)
{
  double angmom_mag =
    sqrt(halo->AngMom[0] * halo->AngMom[0] + halo->AngMom[1] * halo->AngMom[1] + halo->AngMom[2] * halo->AngMom[2]);
  return angmom_mag / (1.414213562 * halo->Vvir * halo->Rvir);
}

// Here you are adding functions that you need to compute the correlation function (needed for boost the probability of getting metal enriched, motivated by clustering
// As a first test these routines are copied from a previous work of Manu when he was a dumb Master student (the originals were written in Python).
// In the future it might be work to see if these can be improved. Parameters from Eisenstein & Hu 1998 or 1999 (EH98, EH99)

void initialize_interpCF_arrays()
{
  FILE* input_fileCF;
  char input_file_nameCF[500];
  char input_baseCF[] = "SpatialCF.dat";
  char modeCF[10] = "r";
  
  int i;
  
  if (run_globals.mpi_rank == 0) {
    sprintf(input_file_nameCF,"%s/%s", run_globals.params.TablesForXHeatingDir, input_baseCF); // ATM is in the same location, you might want change it later!
    input_fileCF = fopen(input_file_nameCF, modeCF);
    
    if (input_fileCF == NULL) {
        mlog("Can't open input file %s!\n", MLOG_MESG, input_file_nameCF);
        exit(1);
      }
    
    // Read in data table
      for (i = 0; i < x_int_NCFVALS; i++) {
        fscanf(input_fileCF,
               "%g %g %g",
               &x_int_zvals[i],
               &x_int_radvals[i],
               &x_int_CFvals[i]);
      }

      fclose(input_fileCF);
    }
    
  // broadcast the values to all cores
  MPI_Bcast(&x_int_zvals, sizeof(x_int_zvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_radvals, sizeof(x_int_radvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_CFvals, sizeof(x_int_CFvals), MPI_BYTE, 0, run_globals.mpi_comm);
}

void initialize_interpSigma_arrays()
{
  FILE* input_fileSigma;
  char input_file_nameSigma[500];
  char input_baseSigma[] = "SigmaVals.dat";
  char modeSigma[10] = "r";
  
  int i;
  
  if (run_globals.mpi_rank == 0) {
    sprintf(input_file_nameSigma,"%s/%s", run_globals.params.TablesForXHeatingDir, input_baseSigma); // ATM is in the same location, you might want change it later!
    input_fileSigma = fopen(input_file_nameSigma, modeSigma);
    
    if (input_fileSigma == NULL) {
        mlog("Can't open input file %s!\n", MLOG_MESG, input_file_nameSigma);
        exit(1);
      }
    
    // Read in data table
      for (i = 0; i < x_int_NCFVALS; i++) {
        fscanf(input_fileSigma,
               "%g %g %g",
               &x_int_zvals[i],
               &x_int_rvirvals[i],
               &x_int_Sigmavals[i]);
      }

      fclose(input_fileSigma);
    }
    
  // broadcast the values to all cores
  MPI_Bcast(&x_int_zvals, sizeof(x_int_zvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_rvirvals, sizeof(x_int_rvirvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_Sigmavals, sizeof(x_int_Sigmavals), MPI_BYTE, 0, run_globals.mpi_comm);
}

double read_SpatialCF(double redshift, double Radius) //Radius in cMpc/h
{
  int z_index = 0;
  int R_index = 0;
  int i = 0;
  int ii = 0;
  //mlog("Tot values %d:", MLOG_MESG, x_int_NCFVALS);
  for (i = 0; i < x_int_NCFVALS; i++) {
    if (fabs(x_int_zvals[i] - redshift) <= 0.07) {
      z_index = i;
      for (ii = z_index; ii < x_int_NCFVALS; ii++) {
        if (fabs(x_int_zvals[ii] - redshift) > 0.07 && Radius < MAX_RAD) {
          mlog("Error, you didn't find the radius value for %f!\n", MLOG_MESG, Radius);
          exit(1);
          }
        if (fabs((Radius - x_int_radvals[ii]) / Radius) < 0.09) {
          R_index = ii;
          break;
          }
        }
      if (Radius >= MAX_RAD) // If that's the case take the largest value
        R_index = ii;
      break;
      }
    }
    //mlog("Index value %d %d:", MLOG_MESG, z_index, R_index);
    //mlog("Red value %f", MLOG_MESG, x_int_zvals[z_index]);
    //mlog("Radius value %f", MLOG_MESG, x_int_radvals[R_index]);         
  return x_int_CFvals[R_index];
}

double read_Sigma(double redshift, double RvirVal) //Radius in cMpc/h
{
  int z_index = 0;
  int Rvir_index = 0;
  int i = 0;
  int ii = 0;
  //mlog("Tot values %d:", MLOG_MESG, x_int_NCFVALS);
  for (i = 0; i < x_int_NCFVALS; i++) {
    if (fabs(x_int_zvals[i] - redshift) <= 0.07) {
      z_index = i;
      for (ii = z_index; ii < x_int_NCFVALS; ii++) {
        if (fabs(x_int_zvals[ii] - redshift) > 0.07 && RvirVal < MAX_Rvir) {
          mlog("Error, you didn't find the Rvir value!\n", MLOG_MESG);
          //exit(1);
          }
        if (fabs((RvirVal - x_int_rvirvals[ii]) / RvirVal) < 0.09) {
          Rvir_index = ii;
          break;
          }
        }
      if (RvirVal >= MAX_Rvir) // If that's the case take the largest value
        Rvir_index = ii;
      break;
      }
    }
    //mlog("Index value %d %d:", MLOG_MESG, z_index, Rvir_index);
    //mlog("Red value %f", MLOG_MESG, x_int_zvals[z_index]);
    //mlog("Radius value %f", MLOG_MESG, x_int_rvirvals[Rvir_index]);         
  return x_int_Sigmavals[Rvir_index];
}

double NLBias(double Dist_Radius, double Halo_Mass, double redshift) //From your fitting function, parameters in vir_properties.h (Input in internal units
{
  double little_h = run_globals.params.Hubble_h;
  //Dist_Radius /= little_h;  //DOUBLE CHECK DIMENSION OF THE BUBBLE!! I believe this should be divided by little_h
  
  Halo_Mass = Halo_Mass * 1e10 / little_h;
  
  return (Psi_Norm * pow(Dist_Radius / 0.01, Alpha_ind) * pow(Halo_Mass / 1e6, Beta_ind) * pow(redshift / 20.0, Gamma_ind));
}
