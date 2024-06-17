/* Void Function AllocArrays(): Allocates memory for various
   arrays used in SiAsgt.

   Written by Manoj Warrier (manoj.warrier@ipp.mpg.de) on 15-10-2002.

*/
#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

void
AllocArrays ()
{
  printf ("Allocating memory\n");

/* Allocating memory for number of particles of various species */
  if ((NParts = (int *) malloc (NSpecies * sizeof (int))) == NULL)
    {
      printf ("Not enough memory to allocate space for NParts\n");
      exit (1);
    }

/* Allocating memory for jump frequencies */
  if ((JmpFreq = (double *) malloc (NSpecies * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for JmpFreq\n");
      exit (1);
    }

/* Allocating memory for jump Rates */
  if ((JmpRate = (double *) malloc (NSpecies * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for JmpRate\n");
      exit (1);
    }

/* Allocating memory for Jump distances */
  if ((JmpDist = (double *) malloc (NSpecies * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for JmpDist\n");
      exit (1);
    }

/* Allocating memory for migration activation energies */
  if ((Em = (double *) malloc (NSpecies * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for Em\n");
      exit (1);
    }

/* Allocating memory for standard deviations of initial distributions */
  if ((StdDev = (double *) malloc (NSpecies * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for StdDev\n");
      exit (1);
    }

/* Allocating memory for Interaction distances */
  if ((InteractDist = (double *) malloc (NIntracts * sizeof (double))) ==
      NULL)
    {
      printf ("Not enough memory to allocate space for StdDev\n");
      exit (1);
    }

  printf ("Allocated memory\n");
}
