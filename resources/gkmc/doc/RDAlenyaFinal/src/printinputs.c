/*!
    void function PrintInputs() - Prints Inputs for DiG
	onto a log file DiG.log. The inputs are evident from
        what is being printed. Else check the file readinputs.c
        for details on the inputs.

    Written by Manoj Warrier (manoj.warrier@ipp.mpg.de) on
    17-02-2004.

*/
#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

void
PrintInputs ()
{
  /* Declaring variables */

  /* External variables */
  extern int NParts[NSpecies], NJumps[NSpecies], MaxSteps, NTrials, NPosition,
      NTransient;
  extern double *JmpFreq[NSpecies], *JmpDist[NSpecies], *Em[NSpecies],
      InteractDist[NRecomb], GammaP[NSpecies], StickCoef[NSpecies],
      Temp;

  FILE *file_ptr;

  /* Variables used only in this function */
  static int i, j;

  #ifdef DBuG1
  printf ("Printing Inputs\n");
  #endif /* DBuG1 */

  /* Opening the log file ReactDiff.log */
  file_ptr = fopen ("../out/RDAlenya.log", "w");
  if (file_ptr == NULL)
    {
      printf ("could not open ../out/RDAlenya.log\n");
    }

  /* Printing the variables read in */
  for (i = 0; i < NSpecies; i++)
    {
      fprintf (file_ptr, "Particle Flux (/m^2/s) %e\n", GammaP[i]);
      fprintf (file_ptr, "Sticking Coefficient %e\n", StickCoef[i]);
      fprintf (file_ptr, "Species number %d\n", i);
      fprintf (file_ptr, "Number of particles %d\n", NParts[i]);
      fprintf (file_ptr, "Number of varieties of jumps %d\n", NJumps[i]);
      for (j = 0; j < NJumps[i]; j++)
	{
	  fprintf (file_ptr, "Jump Frequency %e\n", JmpFreq[i][j]);
	  fprintf (file_ptr, "Jump Distance %e\n", JmpDist[i][j]);
	  fprintf (file_ptr, "Migration Energy %e\n", Em[i][j]);
	}
    }
  for (i = 0; i < NRecomb; i++)
    {
      fprintf (file_ptr, "Interaction number %d\n", i + 1);
      printf ("Interaction number %d\n", i + 1);
      fprintf (file_ptr, "Interaction Distance %e\n", InteractDist[i]);
      printf ("Interaction Distance %e\n", InteractDist[i]);
    }
  fprintf (file_ptr, "Max steps before stoping simulation %d\n", MaxSteps);
  fprintf (file_ptr, "Number of trials %d\n", NTrials);
  fprintf (file_ptr, "Temperature of the Material %e\n", Temp);
  fprintf (file_ptr, "No of iterations for position output %d\n", NPosition);
  fprintf (file_ptr, "No of iterations for transient output %d\n", NTransient);

  fflush (file_ptr);

  #ifdef DBuG1
  printf ("Done Printing Inputs\n");
  #endif /* DBuG1 */
}
