/*!
    Function JumpOnSurface (int TheSpecie, int ThePart, int TheRate):
    Carries out transport corresponding to "TheRate" subsequent to
    desorbtion/detrapping of "ThePart" particle of "TheSpecie" species.

    It jumps particle "ThePart" of specie TheSpecie by a step size equal
    to a JmpDist[TheSpecie][TheRate] in X,Y random direction chosen by a
    2-d uniform random number generator.

    Last modified by Manoj Warrier (manoj.warrier - atnospam - gmail.com)
    on 31-03-2008
*/
#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

int ApplyBC (int k, int IE);
int FindCellNum (int k, int i);
void MakeAJump (int k, int i, double dx, double dy, double dz);

void
JumpOnSurface (int TheSpecie, int ThePart, int TheRate)
{
  /*!! Declaring external variables */
  extern double *JmpDist[NSpecies];
  extern int *OccFlag;
  extern gsl_rng *rdir2d;

  /*!! Declaring internal variables */
  double dx, dy, dz;
  int CellNum;

  #ifdef DBuG2
  printf("Entering JumpOnSurface\n");
  #endif /* DBuG2 */

  /* Finding the Cell number */
  CellNum = FindCellNum (TheSpecie, ThePart);
  if (CellNum == -1)
    {
      printf ("Error in JumpOnSurface\n");
      printf ("CellNum = %d; Particle outside the material.\n", CellNum);
      exit (0);
    }
  if (OccFlag[CellNum] != 2)
    {
      printf ("Error in JumpOnSurface\n");
      printf ("Occflag[%d] = %d\n", CellNum, OccFlag[CellNum]);
      printf ("TheSpecie = %d, ThePart = %d\n", TheSpecie,ThePart);
      printf ("The particle should have been at a surface\n");
      exit (0);
    }

  /*!! Creating a random step in 3d */
  gsl_ran_dir_2d (rdir2d, &dx, &dy);
  /*!! Jumping */
  dx = dx * JmpDist[TheSpecie][TheRate];
  dy = dy * JmpDist[TheSpecie][TheRate];
  dz = 0.0;
  MakeAJump (TheSpecie, ThePart, dx, dy, dz);
  /*!! Applying BC */
  if (ApplyBC (TheSpecie, ThePart) == 1)
    {
      printf("Jumped out of the surface during trackback\n");
      printf("Error!! This should not happen (JumpOnSurface)\n");
      exit(0);
    }

  #ifdef DBuG2
  printf("Exiting JumpOnSurface\n");
  #endif
}
