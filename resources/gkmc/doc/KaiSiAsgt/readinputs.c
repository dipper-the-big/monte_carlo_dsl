/*  ReadInputs - C void function which reads Inputs for the program SiAsgt.

    Written by Manoj Warrier (manoj.warrier@ipp.mpg.de) on 14-10-2002

    Input Parameters Description:
	NSpecies -> Number of species.
	NParts[i] -> Number of particles of type i (i=1-n_species).
	Lx, Ly, Lz -> Material dimensions (Angstroms).
	JmpFreq[i] -> Jump frequency of specie i (1/femto-sec).
	JmpDist[i] -> Jump distance of specie i (Angstroms).
	Em[i] -> Migration activation energy of specie i (eV).
        StdDev -> Standard deviation in initial position of specie.
	InteractDist[i][j] -> Interaction distance for specie i and j
		(Angstroms) for the species to recombine.
	Temp -> Temperature of the material (Kelvins).
	MaxSteps -> Maximum number of steps without interaction
		to stop simulation.
	NTrials -> Number of trials for statistics.
*/

#include <stdio.h>
#include "GlobalVars.h"

void
ReadInputs ()
{

/* Declaring variables */
  int nloop;
  FILE *file_ptr;

  printf ("Reading Inputs\n");

/* Opening the Input file mdlearn.in */
  file_ptr = fopen ("SiAsgt.inp", "r");
  if (file_ptr == NULL)
    {
      printf ("could not read SiAsgt.inp\n");
    }

/* Reading in the variables
   For each of the specie, input NParts, jmpFreq, JmpDist, Em and
   StdDev */

  for (nloop = 0; nloop < NSpecies; nloop++)
    {
      fscanf (file_ptr, "%d\n", &NParts[nloop]);
      fscanf (file_ptr, "%lf\n", &JmpFreq[nloop]);
      fscanf (file_ptr, "%lf\n", &JmpDist[nloop]);
      fscanf (file_ptr, "%lf\n", &Em[nloop]);
      fscanf (file_ptr, "%lf\n", &StdDev[nloop]);
    }
/* If Vacancy and Interstitial come closer that InteractDist, they recombine */
  for (nloop = 0; nloop < NIntracts; nloop++)
    {
      fscanf (file_ptr, "%lf\n", &InteractDist[nloop]);
    }
  fscanf (file_ptr, "%d\n", &MaxSteps);
  fscanf (file_ptr, "%d\n", &NTrials);
  fscanf (file_ptr, "%lf\n", &Temp);
  fscanf (file_ptr, "%lf\n", &Lx);
  fscanf (file_ptr, "%lf\n", &Ly);
  fscanf (file_ptr, "%lf\n", &Lz);

  printf ("Done Reading Inputs\n");
}
