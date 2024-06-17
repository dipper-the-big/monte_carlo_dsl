/*!
   InteractFunc (int ievent, int ichosen, int intflag, int iinteract,
     int iproduced): void function to carry out the interactions
     of the chosen event.

   Variables set/returned:
     Initial (rxi,ryi,rzi) and current (rx,ry,rz) coordinates, boundary
     condition index (ixcount, iycount,izcount) and number (NPartsKMC) of
     the chosen, interacting and produced species.

   Explanation of passed variables:
     The "ievent^th" particle of the "ichosen" specie (chosen
     by the BKL algorithm) interacts with the "intflag^th" particle
     of the "iinteract" specie to produce the "iproduced" specie.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

void GetDiff (int TheSpecie, int ThePart);

void
InteractFunc (int ievent, int ichosen, int intflag, int iinteract,
	      int iproduced)
{
  /*!! Declaring external variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies], *rxi[NSpecies],
      *ryi[NSpecies], *rzi[NSpecies], *ResideTi[NSpecies], *ResideTf[NSpecies],
      CurrentTime, dt, dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt;
  extern int *ixcount[NSpecies], *iycount[NSpecies], *izcount[NSpecies],
      *icurrent[NSpecies], *scurrent[NSpecies], NPartsKMC[NSpecies];

  /*!! Declaring internal variables */
  int k, idum;
  FILE *HAtomDiff;

  #ifdef DBuG2
  printf("Entered InteractFunc\n");
  #endif /* DBuG2 */

  /*!! Opening Onput file for writing out H Diffusion data */
  HAtomDiff = fopen ("../out/HAtomDiff.out", "a");
  if (HAtomDiff == NULL)
    {
      printf ("could not write ../out/HAtomDiff.out\n");
      exit(0);
    }

  /* Finding the Diffusion Co-efficient */
  dt = 0.0;
  dxsqrbydt = 0.0; 
  dysqrbydt = 0.0; 
  dzsqrbydt = 0.0; 
  dsqrbydt = 0.0; 
  GetDiff (ichosen, ievent);
  fprintf (HAtomDiff, "%e %e %e %e %e\n",
    dt, dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt);
  #ifdef DBuG2
  printf ("InteractFunc %e %e %e %e %e\n",
    dt, dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt);
  #endif /* DBuG2 */
  dt = 0.0;
  dxsqrbydt = 0.0; 
  dysqrbydt = 0.0; 
  dzsqrbydt = 0.0; 
  dsqrbydt = 0.0; 
  GetDiff (iinteract, intflag);
  fprintf (HAtomDiff, "%e %e %e %e %e\n",
    dt, dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt);
  #ifdef DBuG2
  printf ("InteractFunc %e %e %e %e %e\n",
    dt, dxsqrbydt, dysqrbydt, dzsqrbydt, dsqrbydt);
  #endif /* DBuG2 */

  fclose (HAtomDiff);

  /* Creating the produced specie */
  rx[iproduced][NPartsKMC[iproduced]] = rx[iinteract][intflag];
  ry[iproduced][NPartsKMC[iproduced]] = ry[iinteract][intflag];
  rz[iproduced][NPartsKMC[iproduced]] = rz[iinteract][intflag];
  rxi[iproduced][NPartsKMC[iproduced]] = rx[iinteract][intflag];
  ryi[iproduced][NPartsKMC[iproduced]] = ry[iinteract][intflag];
  rzi[iproduced][NPartsKMC[iproduced]] = rz[iinteract][intflag];
  ixcount[iproduced][NPartsKMC[iproduced]] = ixcount[iinteract][intflag];
  iycount[iproduced][NPartsKMC[iproduced]] = iycount[iinteract][intflag];
  izcount[iproduced][NPartsKMC[iproduced]] = izcount[iinteract][intflag];
  icurrent[iinteract][intflag] = NPartsKMC[iproduced];
  scurrent[iinteract][intflag] = iproduced;
  icurrent[ichosen][ievent] = NPartsKMC[iproduced];
  scurrent[ichosen][ievent] = iproduced;
  ResideTi[iproduced][NPartsKMC[iproduced]] = CurrentTime;
  ResideTf[iinteract][intflag] = CurrentTime;

  /* incrementing the number of produced specie */
  NPartsKMC[iproduced] = NPartsKMC[iproduced] + 1;

  /* Correcting the redistribution error when ichosen == iinteract */
  if ( (ichosen == iinteract) && (ievent < intflag) )
    {
      idum = ievent;
      ievent = intflag;
      intflag = idum;
    }

  /* redistributing chosen species */
  for (k = ievent; k < (NPartsKMC[ichosen] - 1); k++)
    {
      rx[ichosen][k] = rx[ichosen][k + 1];
      ry[ichosen][k] = ry[ichosen][k + 1];
      rz[ichosen][k] = rz[ichosen][k + 1];
      rxi[ichosen][k] = rxi[ichosen][k + 1];
      ryi[ichosen][k] = ryi[ichosen][k + 1];
      rzi[ichosen][k] = rzi[ichosen][k + 1];
      ixcount[ichosen][k] = ixcount [ichosen][k+1];
      iycount[ichosen][k] = iycount [ichosen][k+1];
      izcount[ichosen][k] = izcount [ichosen][k+1];
      icurrent[ichosen][k+1] = k;
    }

  /* Updating number of chosen species */
  NPartsKMC[ichosen] = NPartsKMC[ichosen] - 1;

  /* redistributing interacting atoms */
  for (k = intflag; k < (NPartsKMC[iinteract] - 1); k++)
    {
      rx[iinteract][k] = rx[iinteract][k + 1];
      ry[iinteract][k] = ry[iinteract][k + 1];
      rz[iinteract][k] = rz[iinteract][k + 1];
      rxi[iinteract][k] = rxi[iinteract][k + 1];
      ryi[iinteract][k] = ryi[iinteract][k + 1];
      rzi[iinteract][k] = rzi[iinteract][k + 1];
      ixcount[iinteract][k] = ixcount[iinteract][k+1];
      iycount[iinteract][k] = iycount[iinteract][k+1];
      izcount[iinteract][k] = izcount[iinteract][k+1];
      icurrent[iinteract][k+1] = k;
    }

  /* Updating number of interacting atoms */
  NPartsKMC[iinteract] = NPartsKMC[iinteract] - 1;

  #ifdef DBuG2
  printf("NPartsKMC[%d] = %d\n", ichosen, NPartsKMC[ichosen]);
  printf("NPartsKMC[%d] = %d\n", iinteract, NPartsKMC[iinteract]);
  printf("NPartsKMC[%d] = %d\n", iproduced, NPartsKMC[iproduced]);
  #endif /* DBuG2 */

  #ifdef DBuG2
  printf("Exiting InteractFunc\n");
  #endif /* DBuG2 */
}
