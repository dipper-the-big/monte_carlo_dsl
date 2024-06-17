/*!! 
    void function ReallocMem - Re-allocates memory for the position
    and transport tracking related variables when new particles are
    added to the system.

    Written by Manoj (Manoj@ipr.res.in) Last modified on 31-10-2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

void
ReallocMem (void)
{
  /*!! Initialising variables */

  /*!! Declaring external variables */
  extern int NParts[NSpecies];
  extern double  *rx[NSpecies], *ry[NSpecies], *rz[NSpecies], *rxi[NSpecies],
      *ryi[NSpecies], *rzi[NSpecies], *ResideTi[NSpecies], *ResideTf[NSpecies];
  extern int *ixcount[NSpecies], *iycount[NSpecies], *izcount[NSpecies],
      *icurrent[NSpecies], *scurrent[NSpecies];

  /*!! Declaring variables initialised in this Function */

  /*!! Declaring internal variables */
  static int k, MaxHSolute, MaxH2Solute, NStore, MaxPart;
  double *tempr[NSpecies];
  int *tempi[NSpecies];

  #ifdef DBuG2
  printf ("Reallocating Memory\n");
  #endif /* DBuG2 */

  /*!!
     Allocating memory for current (rx,ry,rz) and initial (rxi,ryi,rzi)
     position variables and also the boundary condition index for the
     various species
   */
  MaxHSolute = NParts[0] + 2.0e0 * NParts[1];
  MaxH2Solute = NParts[0] / 2.0e0 + NParts[1];

  MaxPart = 0;

  for (k = 0; k < NSpecies; k++)
    {
      if (NParts[k] >= MaxPart) MaxPart = NParts[k];
      switch (k)
	{
	case 0: NStore = MaxHSolute;  break;
	case 1: NStore = MaxH2Solute; break;
	default:
	  printf ("Error in InitParams, value of k = %d impossible\n", k);
	  printf ("Possible error source: Check NSpecies");
	  exit (0);
	}
      #ifdef DBuG1
        printf("NStore = %d\n", NStore);
      #endif /* DBuG1 */
      if ((tempr[k] = realloc (rx[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for rx\n");
	  exit (0);
	}
      rx[k] = tempr[k];
      if ((tempr[k] = realloc (ry[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for ry\n");
	  exit (0);
	}
      ry[k] = tempr[k];
      if ((tempr[k] = realloc (rz[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for rz\n");
	  exit (0);
	}
      rz[k] = tempr[k];
      if ((tempr[k] = realloc (ResideTi[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for ResideTi\n");
	  exit (0);
	}
      ResideTi[k] = tempr[k];
      if ((tempr[k] = realloc (ResideTf[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for ResideTf\n");
	  exit (0);
	}
      ResideTf[k] = tempr[k];
      if ((tempr[k] = realloc (rxi[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for rx\n");
	  exit (0);
	}
      rxi[k] = tempr[k];
      if ((tempr[k] = realloc (ryi[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for ry\n");
	  exit (0);
	}
      ryi[k] = tempr[k];
      if ((tempr[k] = realloc (rzi[k], NStore * sizeof (double))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for rz\n");
	  exit (0);
	}
      rzi[k] = tempr[k];
      if ((tempi[k] = realloc (ixcount[k], NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for ixcount\n");
	  exit (0);
	}
      ixcount[k] = tempi[k];
      if ((tempi[k] = realloc (iycount[k], NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for iycount\n");
	  exit (0);
	}
      iycount[k] = tempi[k];
      if ((tempi[k] = realloc (izcount[k], NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for izcount\n");
	  exit (0);
	}
      izcount[k] = tempi[k];
      if ((tempi[k] = realloc (icurrent[k], NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for icurrent\n");
	  exit (0);
	}
      icurrent[k] = tempi[k];
      if ((tempi[k] = realloc (scurrent[k], NStore * sizeof (int))) == NULL)
	{
	  printf ("Not enough memory to re-allocate space for scurrent\n");
	  exit (0);
	}
      scurrent[k] = tempi[k];
    }

  #ifdef DBuG2
  printf ("Done Reallocating Memory\n");
  #endif /* DBuG2 */
}
