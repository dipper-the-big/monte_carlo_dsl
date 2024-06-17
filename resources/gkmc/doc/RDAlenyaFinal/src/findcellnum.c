/*!! int FindCellNum(int k, int i) : Returns the cell number of
     the $i^{th}$ particle of the $k^{th}$ specie.

     Written by Manoj Warrier (Manoj.Warrier@ipp.mpg.de) on
     10 - 04 - 2004
*/
#include <stdio.h>
#include <math.h>
#include "GlobalVars.h"

int
FindCellNum (int k, int i)
{
  /*!! Declaring external variables */
  extern double *rx[NSpecies], *ry[NSpecies], *rz[NSpecies], CelSize;
  #ifdef DBuG3
  extern double *rxi[NSpecies], *ryi[NSpecies], *rzi[NSpecies];
  #endif /* DBuG3 */
  extern int nx, ny;

  /*!! Declaring internal variables */
  int ix, iy, iz, CellNum;

  #ifdef DBuG3
  printf("Called FindCellNum\n");
  #endif /* DBuG3 */

  /* Finding the Cell number */
  ix = (int) floor(rx[k][i] / CelSize);
  iy = (int) floor(ry[k][i] / CelSize);
  iz = (int) floor(rz[k][i] / CelSize);
  if (rz[k][i] < minz) 
    {
      CellNum = -1;
      return(CellNum);
    }
  CellNum = nx * ny * iz + nx * iy + ix;
  #ifdef DBuG3
  printf("ix = %d, iy = %d, iz = %d\n", ix, iy, iz);
  printf("rx[%d][%d] = %e\n",k, i, rx[k][i]);
  printf("ry[%d][%d] = %e\n",k, i, ry[k][i]);
  printf("rz[%d][%d] = %e\n",k, i, rz[k][i]);
  printf("rxi[%d][%d] = %e\n",k, i, rxi[k][i]);
  printf("ryi[%d][%d] = %e\n",k, i, ryi[k][i]);
  printf("rzi[%d][%d] = %e\n",k, i, rzi[k][i]);
  #endif /* DBuG3 */

  #ifdef DBuG3
  printf("Returning from FindCellNum, CellNum = %d\n", CellNum);
  #endif /* DBuG3 */

  return (CellNum);
}
