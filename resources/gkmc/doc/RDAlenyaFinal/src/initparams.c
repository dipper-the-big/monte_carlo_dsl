/*! 
    void function InitParams - initialises various parameters
    and variables for the program ReactDiff

    Written by Manoj (Manojipr.res.in) Last modified on 28-10-2005.
*/

#include <stdio.h>
#include <math.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

void
InitParams ()
{
  /*!! Initialising variables */

  /*!! Declaring external variables */
  extern int NParts[NSpecies], NJumps[NSpecies];
  extern double *JmpFreq[NSpecies], *Em[NSpecies], Temp;

  /*!! Declaring variables initialised in this Function */

  /*!! Declaring internal variables */
  static int i, j, k, MaxHSolute, MaxH2Solute, NStore, MaxPart;
  static double rng_max;

  FILE *file_ptr;

  #ifdef DBuG1
  printf ("Initialising Parameters\n");
  #endif /* DBuG1 */

  /* Opening the log file ReactDiff.log */
  file_ptr = fopen ("../out/ReactDiff.log", "a");
  if (file_ptr == NULL)
    {
      printf ("could not open ReactDiff.log\n");
    }

  fprintf (file_ptr, "Initialising Parameters\n");

  /*!!
     Initialising the random number generator enviornment for GSL calls
     Creating a generator chosen by the Environment variable GSL_RNG_TYPE
   */
  gsl_rng_env_setup ();
  TBKL = gsl_rng_default;
  urn_BKL = gsl_rng_alloc (TBKL);
  rng_max = (double) gsl_rng_max (urn_BKL);
  rng_min = (double) gsl_rng_min (urn_BKL);
  rng_diff = rng_max - rng_min;
  Trdir2d = gsl_rng_default;
  rdir2d = gsl_rng_alloc (Trdir2d);
  fprintf (file_ptr, "Initialised Random Number Generator enviornment\n");
  printf ("Initialised Random Number Generator enviornment\n");

  /*!!
     Allocating memory for current (rx,ry,rz) and initial (rxi,ryi,rzi)
     position variables and also the boundary condition index for the
     various species
   */
  MaxHSolute = NParts[0] + 2.0e0 * NParts[1];
  fprintf (file_ptr, "MaxHSolute = %d\n", MaxHSolute);
  MaxH2Solute = NParts[0] / 2.0e0 + NParts[1];
  fprintf (file_ptr, "MaxH2Solute = %d\n", MaxH2Solute);

  MaxPart = 0;

  for (k = 0; k < NSpecies; k++)
    {
      if (NParts[k] > MaxPart) MaxPart = NParts[k];
      switch (k)
	{
	case 0: NStore = MaxHSolute;  break;
	case 1: NStore = MaxH2Solute; break;
	default:
	  printf ("Error in InitParams, value of k = %d impossible\n", k);
	  printf ("Possible error source: Check NSpecies");
	  exit (0);
	}
      if ((rx[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for rx\n");
	  exit (0);
	}
      if ((ry[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for ry\n");
	  exit (0);
	}
      if ((rz[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for rz\n");
	  exit (0);
	}
      if ((ResideTi[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for ResideTi\n");
	  exit (0);
        }
      if ((ResideTf[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for ResideTf\n");
	  exit (0);
	}
      if ((rxi[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for rx\n");
	  exit (0);
	}
      if ((ryi[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for ry\n");
	  exit (0);
	}
      if ((rzi[k] = (double *) malloc (NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for rz\n");
	  exit (0);
	}
      if ((ixcount[k] = (int *) malloc (NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to allocate space for ixcount\n");
	  exit (0);
	}
      if ((iycount[k] = (int *) malloc (NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to allocate space for iycount\n");
	  exit (0);
	}
      if ((izcount[k] = (int *) malloc (NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to allocate space for izcount\n");
	  exit (0);
	}
      if ((icurrent[k] = (int *) malloc (NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to allocate space for icurrent\n");
	  exit (0);
	}
      if ((scurrent[k] = (int *) malloc (NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to allocate space for scurrent\n");
	  exit (0);
	}
      if ((JmpRate[k] = (double *) malloc (NJumps[k] * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for JmpRate\n");
	  exit (0);
	}
      for (i = 0; i < NStore; i++)
	{
	  rx[k][i] = 0.0;
	  ry[k][i] = 0.0;
	  rz[k][i] = 0.0;
          ResideTi[k][i] = 0.0;
          ResideTf[k][i] = 0.0;
	  rxi[k][i] = 0.0;
	  ryi[k][i] = 0.0;
	  rzi[k][i] = 0.0;
	  ixcount[k][i] = 0;
	  iycount[k][i] = 0;
	  izcount[k][i] = 0;
	  icurrent[k][i] = 0;
	  scurrent[k][i] = 0;
	}
      ParticleRateSum[k] = 0.0;
      NPartsKMC[k] = 0;
      NPartsDesorb[k] = 0;
      NPartsINI[k] = NParts[k];

      for (j = 0; j < NJumps[k]; j++)
	{
	  JmpRate[k][j] = JmpFreq[k][j] *
	    exp (-Em[k][j] * eVToJoules / K_Boltzmann / Temp);
	  ParticleRateSum[k] = ParticleRateSum[k] + JmpRate[k][j];
          #ifdef DBuG2
	  printf ("JmpRate[%d][%d] = %e\n", k, j, JmpRate[k][j]);
	  printf ("ParticleRateSum[%d] = %e\n", k, ParticleRateSum[k]);
          #endif /* DBuG2 */
	}
    }

  fprintf (file_ptr, "Allocated Memory for:\n");
  fprintf (file_ptr, "                 Position Variables\n");
  fprintf (file_ptr, "                 Boundary Condition index\n");
  printf ("Allocated Memory for:\n");
  printf ("                 Position Variables\n");
  printf ("                 Boundary Condition index\n");

  fprintf (file_ptr, "Initialised all allocated memory to Zero\n");
  printf ("Initialised all allocated memory to Zero\n");

  fprintf (file_ptr, "Allocated Memory and initialised Rates:\n");

  /*!
     Initialising the time step:
     A small time step can be input and the NParts worked out from the
     incident fluxes, GammaP, and sticking coefficients, StickCoeff.
     However this has to be done at the begining of InitParams. I dont
     see any problem following the ansatz I am using now... which is to
     input a begining NParts[i] for various species and see the evolution
     to detailed balance.
  */
  CurrentTime = 0.0;
  DeltaT = 0.0;
  fprintf (file_ptr, "Initialised time step\n");
  printf ("Initialised time step\n");

  /* Initialising the KMC step counter */
  fprintf (file_ptr, "Initialised KMC step counter\n");
  printf ("Initialised KMC step counter\n");

  fprintf (file_ptr, "Done Initialising Parameters\n");

  #ifdef DBuG1
  printf ("Done Initialising Parameters\n");
  #endif /* DBuG1 */
}
