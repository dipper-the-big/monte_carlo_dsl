#define NSpecies 3
#define NIntracts 1

/* You must have as many *rx, *ry, *rz (position variables)
   as the number of species NSpecies. I could have made
   *rx, *ry and *rz a 2 D array, but that would use surplus
   memory if the number of particles of different species
   are not equal ...

   Also remember to allocate memory for the position
   variables in initmaterial.c after U add species.
*/

double *JmpFreq, *Em, *JmpDist, *InteractDist,
  *StdDev, *JmpRate, Temp, Lx, Ly, Lz,
  *rxi, *ryi, *rzi, *rxv, *ryv, *rzv, *Ri,
  *rxr, *ryr, *rzr, xmin, xmax, ymin, ymax, zmin, zmax;

int *NParts, MaxSteps, NTrials;
