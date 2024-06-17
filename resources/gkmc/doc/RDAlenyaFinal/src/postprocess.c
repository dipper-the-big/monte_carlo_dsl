/*
   PostProcess (): Void function which does statistical post processing
   of parameters of interest after the BKL algorithm is carried out by
   the KMCDriver.

   Last modified on 07-05-2003
   by Manoj Warrier (manoj.warrier@ipp.mpg.de)
*/

#include <stdio.h>
#include "GlobalVars.h"

/*!! Function DiffCoeffEval (int k, double timei, double timef)
     described in file diffcoeffeval.c
*/
void DiffCoeffEval (int TheSpecie);

void
PostProcess (M)
{
  /*!! Declaring external variables */
  extern int NPartsKMC[NSpecies], NPartsDesorb[NSpecies];
  extern int NParts[NSpecies];
  extern double *FracKMC[NSpecies], *FracDesorb[NSpecies];
  #ifndef IncludeReactions
    extern double *DiffCoeff[NSpecies], *DiffCoeffx[NSpecies],
      *DiffCoeffy[NSpecies], *DiffCoeffz[NSpecies], dsqrbydtsum,
      dxsqrbydtsum, dysqrbydtsum, dzsqrbydtsum;
  #endif /* IncludeReactions not defined */

  /*!! Declaring local variables */
  int k;

  for (k = 0; k < NSpecies; k++)
    {
      printf("Specie number = %d\n",k);
      if (NParts[k] != 0)
        {
	  /* If reactions are included, interactfunc.c calls GetDiff and
	     diffusion coeff outputs are made in ../out/HAtomDiff.out */
	  #ifndef IncludeReactions
          DiffCoeffEval(k);
	  /* Diffusion coeff = mean sqr displacement / (2*dimension*dt) */
	  DiffCoeff[k][M] = dsqrbydtsum / 4.0 / NParts[k];
	  DiffCoeffx[k][M] = dxsqrbydtsum / 2.0 / NParts[k];
	  DiffCoeffy[k][M] = dysqrbydtsum / 2.0 / NParts[k];
	  DiffCoeffz[k][M] = dzsqrbydtsum / 2.0 / NParts[k];
          printf("DiffCoeff[%d][%d] = %e\n", k, M, DiffCoeff[k][M]);
          printf("DiffCoeffx[%d][%d] = %e\n", k, M, DiffCoeffx[k][M]);
          printf("DiffCoeffy[%d][%d] = %e\n", k, M, DiffCoeffy[k][M]);
          printf("DiffCoeffz[%d][%d] = %e\n", k, M, DiffCoeffz[k][M]);
          #endif /* IncludeReactions not defined */
          FracKMC[k][M] =
	    (double) NPartsKMC[k] / (double) NParts[k];
          FracDesorb[k][M] =
	    (double) NPartsDesorb[k] / (double) NParts[k];
          printf("FracKMC[%d][%d] = %e\n", k,M,FracKMC[k][M]);
          printf("FracDesorb[%d][%d] = %e\n", k,M,FracDesorb[k][M]);
          printf("NPartsKMC[%d] = %d\n", k,NPartsKMC[k]);
          printf("NPartsDesorb[%d] = %d\n", k,NPartsDesorb[k]);
        }
      else
        {
          printf("Nparts[%d] = %d\n", k, NParts[k]);
        }
    }
}
