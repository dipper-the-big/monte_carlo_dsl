/*!
   AllocStatVarMem (): Void function to allocate memory to the variables
   for which statistics post processing is done.

   Last modified on 07-05-2003 by
   Manoj Warrier (manoj.warrier@ipp.mpg.de)
*/

#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

void
AllocStatVarMem ()
{
  /*!! Declaring Variables */
  /*!! Declaring external variables */
  extern int NTrials;

  /*!! Declaring internal variables */
  static int i, k;

  /*!! Allocating memory to variables used in statistical analysis
   */
  if ((Zeit = (double *) malloc (NTrials * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for Zeit\n");
      exit (1);
    }

  for (k = 0; k < NSpecies; k++)
    {
      if ((Frac[k] = (double *) malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for Frac\n");
	  exit (1);
	}
      if ((FracKMC[k] = (double *) malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for FracKMC\n");
	  exit (1);
	}
      if ((FracDesorb[k] = (double *) malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for FracDesorb\n");
	  exit (1);
	}
      if ((DiffCoeffx[k] = (double *)
	   malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for DiffCoeffx\n");
	  exit (1);
	}
      if ((DiffCoeffy[k] = (double *)
	   malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for DiffCoeffy\n");
	  exit (1);
	}
      if ((DiffCoeffz[k] = (double *)
	   malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for DiffCoeffz\n");
	  exit (1);
	}
      if ((DiffCoeff[k] = (double *)
	   malloc (NTrials * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for DiffCoeff\n");
	  exit (1);
	}
    }

  /*!! Initialising the memory allocated to the variables used in
     statistical analysis to 0.0
   */
  for (i = 0; i < NTrials; i++)
    {
      Zeit[i] = 0.0;
      for (k = 0; k < NSpecies; k++)
	{
	  Frac[k][i] = 0.0;
	  FracKMC[k][i] = 0.0;
	  FracDesorb[k][i] = 0.0;
	  DiffCoeffx[k][i] = 0.0;
	  DiffCoeffy[k][i] = 0.0;
	  DiffCoeffz[k][i] = 0.0;
	  DiffCoeff[k][i] = 0.0;
	}
    }
}
