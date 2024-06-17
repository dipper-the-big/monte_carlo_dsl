/* Function InitMaterial () -> Initialises the material being modelled.
   i.e.  it sets the size of the region being modelled, and also the
   distribution of various particles that participate in the KMC.

   Probably this is the place to introduce a crystalline target structure
   and a porous geometry too. We will cross that bridge when we come to
   it.

   Written by Manoj (manoj.warier@ipp.mpg.de) on 16-10-2002.
*/

#include <stdio.h>

#include "GlobalVars.h"
#include "GlobalGSL.h"

void
InitMaterial ()
{
  /* Initialising Variables */
  extern int *NParts;
  extern gsl_rng *r_rdir3d;
  extern double *rxi;
  extern double *ryi;
  extern double *rzi;
  extern double *rxv;
  extern double *ryv;
  extern double *rzv;
  extern double *StdDev;
  int i, j, k;
  double posr, sigmai, sigmav, rx, ry, rz, dx, dy, dz, dist;

//  printf ("Initialising the material\n");

  sigmai = StdDev[0];
  sigmav = StdDev[1];

  for (i = 0; i < NParts[0]; i++)
    {
      posr = gsl_ran_gaussian (r_gauss, sigmai);
      gsl_ran_dir_3d (r_rdir3d, &rx, &ry, &rz);
      rxi[i] = posr * rx;
      ryi[i] = posr * ry;
      rzi[i] = posr * rz;
    }

  for (i = 0; i < NParts[1]; i++)
    {
      posr = gsl_ran_gaussian (r_gauss, sigmav);
      gsl_ran_dir_3d (r_rdir3d, &rx, &ry, &rz);
      rxv[i] = posr * rx;
      ryv[i] = posr * ry;
      rzv[i] = posr * rz;
    }

/*
  // Taking care of recombining particles 
  for (i = 0; i < NParts[0]; i++)
    {
      for (j = 0; j < NParts[1]; j++)
        {
          dx = ( rxi[i] - rxv[j] );
          dy = ( ryi[i] - ryv[j] );
          dz = ( rzi[i] - rzv[j] );
          dist = dx * dx + dy * dy + dz * dz;
          if (dist <= InteractDist[0])
            {
              NParts[0] = NParts[0] - 1;
              NParts[1] = NParts[1] - 1;

              for (k = i; k < (NParts[0]-1); k++)
                {
                  rxi[k] = rxi[k+1];
                  ryi[k] = ryi[k+1];
                  rzi[k] = rzi[k+1];
                }

              for (k = j; k < (NParts[1]-1); k++)
                {
                  rxv[k] = rxv[k+1];
                  ryv[k] = ryv[k+1];
                  rzv[k] = rzv[k+1];
                }
            }
        }
    }
*/

//  printf ("Finished Initialising the material\n");
}
