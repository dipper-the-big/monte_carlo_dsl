/*!! void MakeAJump(int k, int i, double dx, double dy, double dz) :
     Makes the $i^{th}$ particle of the $k^{th}$ specie do a random
     3D jump by a distance equal to JumpDist[k][i]; dx, dy and dz
     provide the random 3d components.

     Written by Manoj Warrier (Manoj.Warrier@ipp.mpg.de) on
     10 - 04 - 2004
*/
#include <stdio.h>
#include "GlobalVars.h"

void
MakeAJump (int TheSpecie, int ThePart, double dx, double dy, double dz)
{
#ifdef DEBUG2
  printf("Entering makeajump\n");
  printf("dx, dy, dz = %e %e %e\n",dx, dy, dz);
  printf("TheSpecie, ThePart = %d %d\n",TheSpecie, ThePart);
#endif /*DEBUG2*/
  /*!! Declaring external variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies];

  /*!! Declaring internal variables */

  /*!! Making a Jump */
  rx[TheSpecie][ThePart] = rx[TheSpecie][ThePart] + dx;
  ry[TheSpecie][ThePart] = ry[TheSpecie][ThePart] + dy;
  rz[TheSpecie][ThePart] = rz[TheSpecie][ThePart] + dz;
#ifdef DEBUG2
  printf("Exiting makeajump\n");
#endif /*DEBUG2*/
  return;
}
