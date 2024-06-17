/*  PrintInputs - C void function which Prints Inputs for SiAsgt
	onto a log file.

    Written by Manoj Warrier (manoj.warrier@ipp.mpg.de) on 15-10-2002

*/

#include <stdio.h>
#include "GlobalVars.h"

void
PrintInputs ()
{
  /* Declaring variables */
  extern double *JmpFreq;
  extern double *Em;
  extern double *JmpDist;
  extern double *InteractDist;
  extern double *StdDev;
  extern double Temp;
  extern double Lx;
  extern double Ly;
  extern double Lz;
  extern int *NParts;
  extern int MaxSteps;
  int nloop;
  FILE *file_ptr;

  printf ("Reading Inputs\n");

  /* Opening the Input file mdlearn.in */
  file_ptr = fopen ("SiAsgt.log", "a");
  if (file_ptr == NULL)
    {
      printf ("could not read SiAsgt.log\n");
    }

  /* Printing the variables read in */
  fprintf (file_ptr, "Number of Species %d\n", NSpecies);
  fprintf (file_ptr, "Number of Interactions %d\n", NIntracts);
  for (nloop = 0; nloop < NSpecies; nloop++)
    {
      fprintf (file_ptr, "Species number %d\n", nloop + 1);
      fprintf (file_ptr, "Number of particles %d\n", NParts[nloop]);
      fprintf (file_ptr, "Jump Frequency %e\n", JmpFreq[nloop]);
      fprintf (file_ptr, "Jump Distance %e\n", JmpDist[nloop]);
      fprintf (file_ptr, "Migration Activation Energy %e\n", Em[nloop]);
      fprintf (file_ptr, "Initial Standard Deviation %e\n", StdDev[nloop]);
    }
  for (nloop = 0; nloop < NIntracts; nloop++)
    {
      fprintf (file_ptr, "Interaction number %d\n", nloop + 1);
      fprintf (file_ptr, "Interaction DIstance %e\n", InteractDist[nloop]);
    }
  fprintf (file_ptr, "Max steps before stoping simulation %d\n", MaxSteps);
  fprintf (file_ptr, "Number of trials %d\n", NTrials);
  fprintf (file_ptr, "Temperature of the Material %e\n", Temp);
  fprintf (file_ptr, "Box Dimensions along X,Y and Z %e %e %e\n", Lx, Ly, Lz);
  fflush (file_ptr);
  fclose (file_ptr);

  printf ("Done Reading Inputs\n");
}
