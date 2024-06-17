/*
   GetDiff (int TheSpecie, int ThePart) : void function to find the mean square
   displacement per unit time of ThePart^th particle of TheSpecie.

   The final position of particle ThePart of specie TheSPecie at time ResideTf
   is given by (rx,ry,rz) with the initial positions at time ResideTi given by
   (rxi,ryi,rzi).

   Last modified by Manoj Warrier (manoj@ipr.res.in) on 11-06-2008.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

void
GetDiff (int TheSpecie, int ThePart)
{
  /*!! Declaring external variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies], *rxi[NSpecies],
      *ryi[NSpecies], *rzi[NSpecies], *ResideTi[NSpecies], *ResideTf[NSpecies],
      Lx, Ly, Lz, CurrentTime, dsqrbydt, dxsqrbydt, dysqrbydt, dzsqrbydt;
  extern int *ixcount[NSpecies], *iycount[NSpecies], *izcount[NSpecies],
      *icurrent[NSpecies];

  /*!! Declaring variables initialised in this function */

  /*!! Declaring interanal variables */
  static int i, k, OrigI;
  static double dx, dy, dz;

  k = TheSpecie;
  i = ThePart;
  OrigI = icurrent[k][i];
  ResideTf[k][OrigI] = CurrentTime;
  dx = (double) ixcount[k][OrigI] * Lx + rx[k][OrigI] - rxi[k][OrigI];
  dy = (double) iycount[k][OrigI] * Ly + ry[k][OrigI] - ryi[k][OrigI];
  dz = (double) izcount[k][OrigI] * Lz + rz[k][OrigI] - rzi[k][OrigI];
  dt = ResideTf[k][OrigI] - ResideTi[k][OrigI];
  dxsqrbydt = dx * dx / dt;
  dysqrbydt = dy * dy / dt;
  dzsqrbydt = dz * dz / dt;
  dsqrbydt = dxsqrbydt + dysqrbydt + dzsqrbydt;
  #ifdef DBuG1
  printf("TheSpecie = %d, ThePart = %d\n", k, i);
  printf("ixcount = %d, iycount = %d, izcount = %d\n",
    ixcount[k][i], iycount[k][i], izcount[k][i]);
  printf("Lx = %e, Ly = %e, Lz = %e\n", Lx, Ly, Lz);
  printf("dx = %e, dy = %e, dz = %e\n", dx, dy, dz);
  printf("dxsqrbydt = %e, dysqrbydt = %e, dzsqrbydt = %e, dsqrbydt = %e\n",
    dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt);
  printf("Final Time = %e, Initial Time = %e\n",
    ResideTf[k][OrigI], ResideTi[k][OrigI]);
  #endif /*DBuG1*/

  #ifdef DBuG1
  printf("Exiting GetDiff\n");
  #endif /* DBuG1 */
}
