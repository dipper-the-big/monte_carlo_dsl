/*
   Desorb : void function to carry out desorption of a specie.

   Variables passed: int TheSpecie
                     int ThePart

   Variables set/returned:
             Initial (rxi,ryi,rzi) and current (rx,ry,rz)
             coordinates, boundary condition index (ixcount,
             iycount,izcount) and number (NPartsKMC) of the
             chosen, interacting and produced species.

   Explanation of passed variables:
      The "ievent^th" particle of the "ichosen" specie (chosen
      by the BKL algorithm) interacts with the "intflag^th" particle
      of the "iinteract" specie to produce the "iproduced" specie.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

int ApplyBC (int k, int IE);
int FindCellNum (int k, int i);
void MakeAJump (int k, int i, double dx, double dy, double dz);

void
Desorb (int TheSpecie, int ThePart)
{
  /*!! Declaring external variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies], *rxi[NSpecies],
      *ryi[NSpecies], *rzi[NSpecies], *ResideTi[NSpecies], *ResideTf[NSpecies],
      CurrentTime, CelSize;
  extern int *ixcount[NSpecies], *iycount[NSpecies], *izcount[NSpecies],
      *icurrent[NSpecies], NPartsKMC[NSpecies], NPartsDesorb[NSpecies];

  /*!! Declaring internal variables */
  int ixcnttmp, iycnttmp, izcnttmp, icurrtmp, DsrbIndex, k, CellNum;
  double dx, dy, dz, rxtmp, rytmp, rztmp, rxitmp, ryitmp, rzitmp, ResTimitmp;


  #ifdef DBuG2
  printf("Entering Desorb\n");
  #endif /* DBuG2 */

  /* Finding the Cell number */
  CellNum = FindCellNum (TheSpecie, ThePart);
  if (CellNum == -1)
    {
      printf ("Error in Desorb\n");
      printf ("CellNum = %d; Particle outside the material.\n", CellNum);
      exit (0);
    }
  if (OccFlag[CellNum] != 2)
    {
      printf ("Error in Desorb\n");
      printf ("Occflag[%d] = %d\n", CellNum, OccFlag[CellNum]);
      printf ("TheSpecie = %d, ThePart = %d\n", TheSpecie,ThePart);
      printf ("The particle should have been at a surface\n");
      exit (0);
    }

  dx = 0.0;
  dy = 0.0;
  dz = -1.0 * CelSize;
  MakeAJump (TheSpecie, ThePart, dx, dy, dz);
  /*!! Applying BC */
  if (ApplyBC (TheSpecie, ThePart) != 1)
    {
      printf("Particle did not desorb!!\n");
      printf("Error!! This should not happen (Desorb)\n");
      exit(0);
    }
  rxtmp = rx[TheSpecie][ThePart];
  rytmp = ry[TheSpecie][ThePart];
  rztmp = rz[TheSpecie][ThePart];
  rxitmp = rxi[TheSpecie][ThePart];
  ryitmp = ryi[TheSpecie][ThePart];
  rzitmp = rzi[TheSpecie][ThePart];
  ixcnttmp = ixcount[TheSpecie][ThePart];
  iycnttmp = iycount[TheSpecie][ThePart];
  izcnttmp = izcount[TheSpecie][ThePart];
  icurrtmp = icurrent[TheSpecie][ThePart];
  ResTimitmp = ResideTi[TheSpecie][ThePart];
  for (k = ThePart; k < (NPartsKMC[TheSpecie] - 1); k++)
    {
      rx[TheSpecie][k] = rx[TheSpecie][k + 1];
      ry[TheSpecie][k] = ry[TheSpecie][k + 1];
      rz[TheSpecie][k] = rz[TheSpecie][k + 1];
      rxi[TheSpecie][k] = rxi[TheSpecie][k + 1];
      ryi[TheSpecie][k] = ryi[TheSpecie][k + 1];
      rzi[TheSpecie][k] = rzi[TheSpecie][k + 1];
      ixcount[TheSpecie][k] = ixcount[TheSpecie][k+1];
      iycount[TheSpecie][k] = iycount[TheSpecie][k+1];
      izcount[TheSpecie][k] = izcount[TheSpecie][k+1];
      ResideTi[TheSpecie][k] = ResideTi[TheSpecie][k+1];
      ResideTf[TheSpecie][k] = ResideTf[TheSpecie][k+1];
      icurrent[TheSpecie][k] = icurrent[TheSpecie][k+1];
    }
  DsrbIndex = NPartsKMC[TheSpecie];
  rx[TheSpecie][DsrbIndex] = rxtmp;
  ry[TheSpecie][DsrbIndex] = rytmp;
  rz[TheSpecie][DsrbIndex] = rztmp;
  rxi[TheSpecie][DsrbIndex] = rxitmp;
  ryi[TheSpecie][DsrbIndex] = ryitmp;
  rzi[TheSpecie][DsrbIndex] = rzitmp;
  ixcount[TheSpecie][DsrbIndex] = ixcnttmp;
  iycount[TheSpecie][DsrbIndex] = iycnttmp;
  izcount[TheSpecie][DsrbIndex] = izcnttmp;
  ResideTi[TheSpecie][DsrbIndex] = ResTimitmp;
  ResideTf[TheSpecie][DsrbIndex] = CurrentTime;
  icurrent[TheSpecie][DsrbIndex] = icurrtmp;
  icurrent[TheSpecie][ThePart] = DsrbIndex;
  NPartsKMC[TheSpecie] = NPartsKMC[TheSpecie] - 1;
  NPartsDesorb[TheSpecie] = NPartsDesorb[TheSpecie] + 1;
  return;
}
