/*!
    void function ReadInputs() - reads Inputs for the program DiG.
    Parameters used which are set in "SimulationParams.h":
	NSpecies -> Number of species.
        NRecomb  -> Total number of recombination reactions.
    Assigned/Input Parameters Description:
        GammaP[i] -> Incident particle flux of i^th species.
        StickCoef[i] - > Sticking co-efficient of i^th species.
	NParts[i] -> Number of particles of type i (i=1-n_species).
	NJumps[i] -> Number of varieties of jumps of particles of
                     type i (i=1-n_species).
	JmpFreq[i][j] -> Jump frequency of the j^th interaction of
                         specie i (1/sec).
	JmpDist[i][j] -> Jump distance of the j^th interaction of
                         specie i (meters).
	Em[i][j] -> Migration activation energy of the j^th interaction
                    of specie i (eV).
	InteractDist[n] -> Minimum distance below which the n^th 
                           interaction recombines (Angstroms).
	MaxSteps -> Maximum number of steps without interaction
	            to stop simulation.
	NTrials -> Number of trials for statistics.
	Temp -> Temperature of the material (Kelvins).
        NPosition -> Number of iterations after which Positions of the various
            particles must be output.
        NTransient -> Number of iterations after which transient quantities of
            interest must be output.
    Last modified on 28-10-2005 by Manoj Warrier (manoj@ipr.res.in)
*/
#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

void
ReadInputs ()
{
  /*!
     Declaring variables
  */

  /*!! Declaring external variables */

  /*!! Variables initialised in this file and used globally */

  /*!! Variables used only in this function */
  static int i, j;
  FILE *file_ptr;

  #ifdef DBuG1
  printf ("Reading Inputs\n");
  #endif /* DBuG1 */

  /*!! Opening the Input file RDAlenya.inp */
  file_ptr = fopen ("../inp/RDAlenya.inp", "r");
  if (file_ptr == NULL)
    {
      printf ("could not read ../inp/RDAlenya.inp\n");
      exit(0);
    }

  /* Reading in the variables */
  for (i = 0; i < NSpecies; i++)
    {
      fscanf (file_ptr, "%lf\n", &GammaP[i]);
      printf ("GammaP[%d] = %e\n", i, GammaP[i]);
      fscanf (file_ptr, "%lf\n", &StickCoef[i]);
      printf ("StickCoef[%d] = %e\n", i, StickCoef[i]);
      fscanf (file_ptr, "%d\n", &NParts[i]);
      fscanf (file_ptr, "%d\n", &NJumps[i]);
      if ((JmpFreq[i] =
	   (double *) malloc (NJumps[i] * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for JmpFreq\n");
	  exit (0);
	}
      if ((JmpDist[i] =
	   (double *) malloc (NJumps[i] * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for JmpDist\n");
	  exit (0);
	}
      if ((Em[i] =
	   (double *) malloc (NJumps[i] * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to allocate space for Em\n");
	  exit (0);
	}
      for (j = 0; j < NJumps[i]; j++)
	{
	  fscanf (file_ptr, "%lf\n", &JmpFreq[i][j]);
	  fscanf (file_ptr, "%lf\n", &JmpDist[i][j]);
	  fscanf (file_ptr, "%lf\n", &Em[i][j]);
	}
    }

  for (i = 0; i < NRecomb; i++)
    {
      fscanf (file_ptr, "%lf\n", &InteractDist[i]);
    }

  fscanf (file_ptr, "%d\n", &MaxSteps);
  fscanf (file_ptr, "%d\n", &NTrials);
  fscanf (file_ptr, "%lf\n", &Temp);
  fscanf (file_ptr, "%d\n", &NPosition);
  fscanf (file_ptr, "%d\n", &NTransient);

  #ifdef DBuG1
  printf ("Done Reading Inputs\n");
  #endif /* DBuG1 */
}
