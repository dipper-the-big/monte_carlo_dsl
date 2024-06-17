/*!!
   Function InitParticles (int *NP) -> Initialises NP particles on the
   material being modelled.

   Last modified by Manoj Warrier (manoj.warier@ipp.mpg.de) on 01-12-2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"


void
InitParticles (int NP[NSpecies])
{
  #ifdef DBuG2
  printf("Initialising particle positions\n");
  #endif /* DBuG2 */

  /*!! Initialising Variables */

  /*!! Declaring external Variables */
  extern gsl_rng *urn_BKL;
  extern double rng_min, rng_diff, *rx[NSpecies], *ry[NSpecies],
      *rz[NSpecies], *rxi[NSpecies], *ryi[NSpecies], *rzi[NSpecies],
      *ResideTi[NSpecies], *ResideTf[NSpecies], Lx, Ly, CelSize,
      CurrentTime;
  extern int NParts[NSpecies], NPartsDesorb[NSpecies], *ixcount[NSpecies],
      *iycount[NSpecies], *izcount[NSpecies], *icurrent[NSpecies],
      *scurrent[NSpecies]; 

  /*!! Declaring internal Variables */
  static int i, k, NIni, NFin, dN, dumi;
  static double urn;

  /*!! Initialising sources of different incident species */
  for (k = 0; k < NSpecies; k++)
    {
      if (NP[k] <= 0)
        {
          NIni = 0;
          NFin = NParts[k];
        }
      else
        {
          NIni = NParts[k] - NPartsDesorb[k] - 1 - NP[k];
          NFin = NParts[k] - NPartsDesorb[k] - 1;
          dN = NFin - NIni;
          /* Shifting the desorbed particle indices */
          for (i = NFin; i < NParts[k]; i++)
            {
              dumi = i - dN;
              #ifdef DBuG3
              printf("dumi = %d, i = %d, k = %d\n", dumi, i, k);
              #endif /* DBuG3 */
              rx[k][i] = rx[k][dumi];
              ry[k][i] = ry[k][dumi];
              rz[k][i] = rz[k][dumi];
              rxi[k][i] = rxi[k][dumi];
              ryi[k][i] = ryi[k][dumi];
              rzi[k][i] = rzi[k][dumi];
              ixcount[k][i] = ixcount[k][dumi];
              iycount[k][i] = iycount[k][dumi];
              izcount[k][i] = izcount[k][dumi];
              ResideTi[k][i] = ResideTi[k][dumi];
              ResideTf[k][i] = ResideTf[k][dumi];
              icurrent[k][i] = icurrent[k][dumi];
              #ifdef DBuG3
              printf("dumi = %d, i = %d, k = %d\n", dumi, i, k);
              #endif /* DBuG3 */
            }
        }

      for (i = NIni; i < NFin; i++)
        {
	  /* X Position */
          urn = ((double) gsl_rng_get (urn_BKL) - rng_min) / rng_diff;
          rx[k][i] = urn * Lx;
          rxi[k][i] = rx[k][i];
          ixcount[k][i] = 0;

	  /* Y Position */
          urn = ((double) gsl_rng_get (urn_BKL) - rng_min) / rng_diff;
          ry[k][i] = urn * Ly;
          ryi[k][i] = ry[k][i];
          iycount[k][i] = 0;

	  /* Z Position */
          rz[k][i] = CelSize / 2.0;
          rzi[k][i] = rz[k][i];
          izcount[k][i] = 0;

          if (NP[k] < 0) ResideTi[k][i] = 0;
	  else ResideTi[k][i] = CurrentTime;
          icurrent[k][i] = i;
          scurrent[k][i] = k;
        }
    }

#ifdef DBuG2
  printf("Done initialising particle positions\n");
#endif /* DBuG2 */
}
