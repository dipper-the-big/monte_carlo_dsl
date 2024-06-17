/*
   StatPostProc (): Void function which does statistical post processing
   of parameters of interest after the BKL algorithm is carried out by
   the KMCDriver.

   Last modified on 07-05-2003
   by Manoj Warrier (manoj.warrier@ipp.mpg.de)
*/
#include <stdio.h>
#include "GlobalVars.h"

void DoStats (char *VarDescrip, double *StatVar, int Num);

void
StatPostProc ()
{
  /*!! Declaring external variables */
  extern int NTrials;
  extern double *FracKMC[NSpecies], 
      *FracDesorb[NSpecies], *DiffCoeffx[NSpecies], *DiffCoeffy[NSpecies],
       *DiffCoeffz[NSpecies], *DiffCoeff[NSpecies], *Zeit;

  /*!! Declaring internal variables */
  static int k;

  for (k = 0; k < NSpecies; k++)
    {
      printf ("Specie Number %d\n", k);
      if (NParts[k] == 0) goto nodiff;
      DoStats ("Diffusion Coeff X", DiffCoeffx[k], NTrials);
      DoStats ("Diffusion Coeff Y", DiffCoeffy[k], NTrials);
      DoStats ("Diffusion Coeff Z", DiffCoeffz[k], NTrials);
      DoStats ("Diffusion Coeff", DiffCoeff[k], NTrials);
      nodiff:
      DoStats ("Frac of atoms left", Frac[k], NTrials);
      DoStats ("FracKMC", FracKMC[k], NTrials);
      DoStats ("FracDesorb", FracDesorb[k], NTrials);
    }
  DoStats ("Time Simulated", Zeit, NTrials);
}
