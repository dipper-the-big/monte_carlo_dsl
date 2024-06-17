/*!!
   Function PosOut (): The function that outputs positions of the various particles.

   Written by Manoj Warrier (manoj@ipr.res.in) on 03-11.2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"

/* int function FindCellNum(int k, int i) described in the file findcellnum.c */
int FindCellNum(int k, int i);

void
PosOut ()
{
  /*!! Declaring variables */

  /*!! Declaring external variables */
  extern int NParts[NSpecies], *icurrent[NSpecies];
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies];
      

  FILE *posout;

  /*!! Declaring variables set in this function */

  /*!! Declaring internal variables */
  static int i, OrigI, CellNum;

  #ifdef DBuG2
  printf ("Starting writing out positions \n");
  #endif /* DBuG2 */

  posout = fopen("../out/Pos.xyz","a");
  if (posout == NULL)
    {
      printf("Error, could not open Pos.xyz\n");
      exit(0);
    }

  for (i = 0; i < NParts[0]; i++)
    {
      OrigI = icurrent[0][i];
      CellNum = FindCellNum(0,OrigI);
      if (CellNum >= 0)
        {
          fprintf(posout,"%e %e %e %d\n",
            rx[0][OrigI],ry[0][OrigI],rz[0][OrigI],OccFlag[CellNum]+1);
        }
      else
        {
          fprintf(posout,"%e %e %e %d\n",
            rx[0][OrigI],ry[0][OrigI],rz[0][OrigI],1);
        }
    }

  fprintf(posout,"bingo\n");
  fflush(posout);
  fclose(posout);
}
