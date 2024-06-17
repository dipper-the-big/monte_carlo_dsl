/*!! void FindInteract(int ievent, int ichosen, int iinteract)

   FindInteract : void function to find the index of the
                  interacting specie that lies closest to the
                  chosen specie for a given interaction.

   Variables passed:
                  int ievent -> ievent^th particle of "ichosen" was
                                chosen by the BKL algorithm.
                  int ichosen -> index of the chosen specie.
                  int iinteract -> index of interacting specie.

   Variables set and returned:
          double MinDist -> Minimum distance between the closest
                            interacting species (returned via GlobalVars.h).
          int InteractSpecie -> Index of the closest interacting specie.
          int MinDistIndex -> The index of the closest particle of
                            InteractSpecie.
Last modified by Manoj Warrier (manoj.warrier@ipp.mpg.de) on
   11-03-2004.
*/

#include <stdio.h>
#include <math.h>
#include "GlobalVars.h"

void
FindInteract (int ievent, int ichosen, int iinteract)
{
  /*!! Declaring External variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies];
  extern int NPartsKMC[NSpecies];
  extern double MinDist;

  /*!! Declaring internal variables */
  static int i;
  static double dx, dy, dz, dist;

  /* loop over interacting specie to find the closest lying
     interacting specie with respect to the chosen specie
   */
  if (ichosen == iinteract)	/* Same specie recombination */
    {
      for (i = 0; i < NPartsKMC[iinteract]; i++)
	{
	  if (i == ievent) i++;
	  dx = rx[iinteract][i] - rx[ichosen][ievent];
	  dy = ry[iinteract][i] - ry[ichosen][ievent];
	  dz = rz[iinteract][i] - rz[ichosen][ievent];
	  dist = sqrt (dx * dx + dy * dy + dz * dz);
	  if (dist < MinDist)
	    {
	      MinDist = dist;
	      MinDistIndex = i;
	      InteractSpecie = ichosen;
	    }
	}
    }
  else				/* different specie recombination */
    {
      for (i = 0; i < NPartsKMC[iinteract]; i++)
	{
	  dx = rx[iinteract][i] - rx[ichosen][ievent];
	  dy = ry[iinteract][i] - ry[ichosen][ievent];
	  dz = rz[iinteract][i] - rz[ichosen][ievent];
	  dist = sqrt (dx * dx + dy * dy + dz * dz);
	  if (dist < MinDist)
	    {
	      MinDist = dist;
	      MinDistIndex = i;
	      InteractSpecie = iinteract;
	    }
	}
    }
}
