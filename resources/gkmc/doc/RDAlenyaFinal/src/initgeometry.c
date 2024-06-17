/*!
    Function InitGeometry(): Function to create a surface.

    Written by Manoj (manoj.warrier - atnospam - gmail.com).
    Last modified on 27-03-2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

void
InitGeometry ()
{
  /*!! This function needs no external variables */

  /*!! Declaring internal variables */
  int i, j, k, m, TotalCells;
  double XPos, YPos, ZPos;

  FILE *inpfile;
  FILE *outfile;

  /*
     n (void)
     Reading in inputs
   */
  if ((inpfile = fopen ("../inp/geometry.inp", "r")) == NULL)
    {
      printf ("Could not open file ../inp/geometry.inp\n");
      exit (1);
    }
  fscanf (inpfile, "%lf\n", &CelSize);
  fscanf (inpfile, "%lf\n", &Lx);
  fscanf (inpfile, "%lf\n", &Ly);
  fscanf (inpfile, "%lf\n", &Lz);
  printf ("Read Geometry Inputs\n");
  #ifdef DBuG2
  printf ("CelSize = %lf\n", CelSize);
  printf ("Lx = %lf\n", Lx);
  printf ("Ly = %lf\n", Ly);
  printf ("Lz = %lf\n", Lz);
  #endif /* DBuG2 */
  /*
     Initialising various parameters in the simulation
   */
  nx = (int) rint (Lx / CelSize);
  ny = (int) rint (Ly / CelSize);
  nz = (int) rint (Lz / CelSize);
  TotalCells = nx * ny * nz;
  printf ("Initialised parameters\n");
  printf ("nx, ny, nz = %d %d %d\n", nx, ny, nz);

  /*
     Allocating memory for flag showing occupation of cell.
   */
  if ((OccFlag = (int *) malloc ((nx * ny * nz) * sizeof (int))) == NULL)
    {
      printf ("Not enough memory to allocate space for OccFlag\n");
      exit (0);
    }

  /*
     Occupied flag initialization
     The plasma facing surface is at z=0,
       OccFlag[m] = 1 => a bulk cell
       OccFlag[m] = 2 => a surface cell
     The Original ReactDiff, in addition to these also had OccFlag[m] = 0,
     to represent voids. Note that this can easily be adapted to create voids
     and study reaction-diffusion on a rough surface.
   */
  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
	  for (k = 0; k < nz; k++)
	    {
	      m = nx * ny * k + nx * j + i;
              if (k == 0) OccFlag[m] = 2;
              else OccFlag[m] = 1;
              #ifdef DBuG2
              printf ("i = %d, j = %d, k = %d, OccFlag[%d] = %d\n",
		i, j, k, m, OccFlag[m]);
              #endif /* DBuG2 */
            }
	}
    }

//  Print Output
  if ((outfile = fopen ("../out/Surface.dat", "w")) == NULL)
    {
      printf ("Cannot open file ../out/Surface.dat\n");
      exit (0);
    }

  for (i = 0; i < nx; i++)
    {
      XPos = ((double) i + 0.5) * CelSize;
      for (j = 0; j < ny; j++)
        {
          YPos = ((double) j + 0.5) * CelSize;
          for (k = 0; k < nz; k++)
            {
              ZPos = ((double) k + 0.5) * CelSize;
              m = nx * ny * k + nx * j + i;
              fprintf (outfile, "%e %e %e %d\n", XPos, YPos, ZPos, OccFlag[m]);
              #ifdef DBuG2
              printf ("XPos = %e, YPos = %e, ZPos = %e, OccFlag[%d] = %d\n",
		XPos, YPos, ZPos, m, OccFlag[m]);
              #endif /* DBuG2 */
            }
        }
    }
  fflush (outfile);
  fclose (outfile);

  printf ("Finished Setting up the surface\n");
  return;
}
