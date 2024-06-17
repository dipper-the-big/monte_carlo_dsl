/*!!
   Function BKLAlgorithm (): The main function to carry out the
   Kinetic Monte Carlo  time step.

   Last Modified by Manoj Warrier (manoj.warrier@ipp.mpg.de)
   on 05-03-2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

/*!! Function Desorb(int TheSpecie, int ThePart) described in file desorb.c */
void Desorb (int TheSpecie, int ThePart);

/*!! Function JumpOnSurface(int TheSpecie, int ThePart, TheRate)
  described in file jumponsurface.c */
void JumpOnSurface (int TheSpecie, int ThePart, int TheRate);

/*!! Function CheckRecomb(int TheSpecie, int ThePart) described
in file checkrecomb.c */
void CheckRecomb (int TheSpecie, int ThePart);

void
BKLAlgorithm ()
{
  #ifdef DBuG2
  printf ("entered BKLAlgorithm\n");
  #endif /* DBuG2 */

  /*!! Declaring varaibles */

  /*!! Declaring external variables */
  extern gsl_rng *urn_BKL;
  extern double rng_min, rng_diff;
  extern int NJumps[NSpecies];
  extern double DeltaT, *JmpRate[NSpecies], ParticleRateSum[NSpecies];
  extern int NPartsKMC[NSpecies];

  /*!! Declaring variables initialised in this function */

  /*!! Declaring local variables */
  static int i, ThePart, TheSpecie, TheRate;
  static double CumilativeRate, urn, SumOfRates[NSpecies + 1],
      WhichSpecies, WhichParticle, WhichRate, RateSum;

  /*!! Calculate the cumulative Rate */
  /*!! Get a uniform random number between [0,1] */
  urn = ((double) gsl_rng_get (urn_BKL) - rng_min) / rng_diff;
  SumOfRates[0] = 0.0;
  for (i = 0; i < NSpecies; i++)
    {
      SumOfRates[i + 1] = SumOfRates[i] +
                ParticleRateSum[i] * NPartsKMC[i];
      #ifdef DBuG2
      printf ("ParticleRateSum[i] = %e\n", ParticleRateSum[i]);
      printf ("SumOfRates[i] = %e\n", SumOfRates[i]);
      #endif /* DBuG2 */
    }
  CumilativeRate = SumOfRates[NSpecies];
  if (CumilativeRate == 0.0)
    {
      printf("Error from bklalgorithm\n");
      printf("Cumilative Rate = 0!!, No process happening\n");
      #ifdef DBuG2
      printf ("leaving BKLAlgorithm\n");
      #endif /* DBuG2 */
      exit(0);
      return;
    }

  /*!! Find the specie chosen for this KMC iteration */
  WhichSpecies = CumilativeRate * urn;
  TheSpecie = 0;
  for (i = 0; i < NSpecies; i++)
    {
      if (WhichSpecies <= SumOfRates[i + 1])
	break;
      TheSpecie++;
    }

  /*!! Find the particle chosen for this KMC iteration */
  WhichParticle = WhichSpecies - SumOfRates[TheSpecie];
  ThePart = (int) floor (WhichParticle / ParticleRateSum[TheSpecie]);

  /*!! Find the particle chosen for this KMC iteration */
  WhichRate = WhichParticle -
     ((double) ThePart * ParticleRateSum[TheSpecie]);
  #ifdef DBuG2
  printf ("WhichParticle = %e\n", WhichParticle);
  printf ("WhichRate = %e\n", WhichRate);
  #endif /* DBuG2 */
  TheRate = 0;
  RateSum = 0.0;
  for (i = 0; i < NJumps[TheSpecie]; i++)
    {
      #ifdef DBuG2
      printf ("JmpRate[TheSpecie][%d] = %e\n", i, JmpRate[TheSpecie][i]);
      #endif /* DBuG2 */
      /*!! This RateSum can be calculated in initparam and sent
           to speedup things I guess */
      RateSum = RateSum + JmpRate[TheSpecie][i];
      if (WhichRate < RateSum)
	break;
      TheRate++;
    }

  #ifdef DBuG2
  printf ("CumilativeRate = %e\n", CumilativeRate);
  printf ("TheSpecie = %d; ThePart = %d; TheRate = %d\n",
           TheSpecie, ThePart, TheRate);
  #endif /* DBuG2 */

  /*!! Carrying out the event */
  switch (TheSpecie)
    {
    case 0:			/*!! H Solute */
      switch (TheRate)
	{
	case 0:		/*!! desorption */
	  Desorb (TheSpecie, ThePart);
	  break;
	case 1:		/*!! Surface diffusion */
	  JumpOnSurface (TheSpecie, ThePart, TheRate);
          #ifdef IncludeReactions
	  CheckRecomb (TheSpecie, ThePart);
          #endif /* IncludeReactions */
	  break;
	default:
	  printf ("Error in bklalgorithm.c\n");
	  printf ("TheSpecie, TheRate = %d %d\n", TheSpecie, TheRate);
	  exit (0);
	}
      break;
    case 1:			/*!! H2 Solute */
      switch (TheRate)
	{
	case 0:		/*!! desorption */
	  Desorb (TheSpecie, ThePart);
	  break;
	default:
	  printf ("Error in bklalgorithm.c\n");
	  printf ("TheSpecie, TheRate = %d %d\n", TheSpecie, TheRate);
	  exit (0);
	}
      break;
    default:
      printf ("Error in bklalgorithm.c\n");
      printf ("TheSpecie = %d\n", TheSpecie);
      exit (0);
    }

  /* Finding and Updating time */
  urn = ((double) gsl_rng_get (urn_BKL) - rng_min) / rng_diff;
  DeltaT = -log (urn) / CumilativeRate;

  #ifdef DBuG2
  printf("DeltaT = %e\n",DeltaT);
  printf ("leaving BKLAlgorithm\n");
  #endif /* DBuG2 */
}
