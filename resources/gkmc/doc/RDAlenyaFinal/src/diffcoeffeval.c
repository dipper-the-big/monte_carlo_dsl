/*
   DiffCoeffEval (int k) : void function to find the mean square displacement
   per unit time of atoms/molecules with specie index k, whose final positions
   are defined at time ResideTf as (rx,ry,rz) with the initial positions at time
   ResideTi being defined by (rxi,ryi,rzi).

   This has been written in a general manner for addition of new species.
 
   Last modified by Manoj Warrier (manoj@ipr.res.in) on 17-11-2005.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

void GetDiff (int TheSpecie, int ThePart);

void
DiffCoeffEval (int k)
{
  /*!! Declaring external variables */
  extern int NParts[NSpecies];

  /*!! Declaring variables initialised in this function */

  /*!! Declaring interanal variables */
  int i;

  dsqrbydt = 0.0;
  dxsqrbydt = 0.0;
  dysqrbydt = 0.0;
  dzsqrbydt = 0.0;
  dsqrbydtsum = 0.0;
  dxsqrbydtsum = 0.0;
  dysqrbydtsum = 0.0;
  dzsqrbydtsum = 0.0;
  for (i = 0; i < NParts[k]; i++)
    {
      GetDiff (k, i);
      dsqrbydtsum = dsqrbydtsum + dsqrbydt;
      dxsqrbydtsum = dxsqrbydtsum + dxsqrbydt;
      dysqrbydtsum = dysqrbydtsum + dysqrbydt;
      dzsqrbydtsum = dzsqrbydtsum + dzsqrbydt;
    }

  #ifdef DBuG1
  printf("Exiting DiffCoeffEval\n");
  #endif /* DBuG1 */
}
