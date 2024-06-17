/*  InitParams - C void function which initialises parameters
    for the program SiAsgt.

    Initially not much to initialise because we have not yet
    introduced the porous structure and we are also not working
    with scaled coordinates.

    Written by Manoj (manoj.warrier@ipp.mpg.de) on 16-10-2002
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

void
InitParams ()
{
  extern int *NParts;
  extern double Lx;
  extern double Ly;
  extern double Lz;
  int NTotal, k;

  printf ("Initialising Parameters\n");

  /* Creating a generator chosen by the Environment
     variable GSL_RNG_TYPE */
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r_urn = gsl_rng_alloc (T);
  r_gauss = gsl_rng_alloc (T);
  r_rdir3d = gsl_rng_alloc (T);
  rng_max = (double) gsl_rng_max (r_urn);
  rng_min = (double) gsl_rng_min (r_urn);
  rng_diff = rng_max - rng_min;

  /* Allocating memory for position variables */
  if ((rxi = (double *) malloc (NParts[0] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rxi\n");
      exit (1);
    }
  if ((ryi = (double *) malloc (NParts[0] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for ryi\n");
      exit (1);
    }
  if ((rzi = (double *) malloc (NParts[0] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rzi\n");
      exit (1);
    }
  if ((rxv = (double *) malloc (NParts[1] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rxv\n");
      exit (1);
    }
  if ((ryv = (double *) malloc (NParts[1] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for ryv\n");
      exit (1);
    }
  if ((rzv = (double *) malloc (NParts[1] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rzv\n");
      exit (1);
    }
  /* Creating enough space to store recombined interstitial-vacancy
     positions. Note NParts[2]=0 after assignment of memory.
   */
  NParts[2] = GSL_MIN_INT (NParts[1], NParts[2]);
  if ((rxr = (double *) malloc (NParts[2] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rxr\n");
      exit (1);
    }
  if ((ryr = (double *) malloc (NParts[2] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for ryr\n");
      exit (1);
    }
  if ((rzr = (double *) malloc (NParts[2] * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for rzr\n");
      exit (1);
    }
  NParts[2] = 0;

  /* Initialising material dimensions */
  xmax = Lx / 2.0;
  xmin = 0.0 - xmax;
  ymax = Ly / 2.0;
  ymin = 0.0 - ymax;
  zmax = Lz / 2.0;
  zmin = 0.0 - zmax;

  /* Initialising the array for rates of interactions */
  NTotal = 0;
  for (k = 0; k < NSpecies; k++)
    {
      NTotal = NTotal + NParts[k];
    }

  if ((Ri = (double *) malloc ((NTotal+1) * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for Ri\n");
      exit (1);
    }

  printf ("Done Initialising Parameters\n");
}
