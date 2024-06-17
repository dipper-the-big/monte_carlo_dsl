/*!
    ApplyBC(int k, int IE) : integer functions which applys
    periodic boundary conditions along X, Y and Z. It returns
    "1" if theparticle has left the surface and "0" otherwise.

    Last modified on 12-04-2004
    by Manoj Warrier (manoj.warrier@ipp.mpg.de)

    Inputs: K  -> Specie index
            IE -> atom index

    BCs will be applied to IE^th atom of the k^th specie
    in simple terms.

    Outputs: Void function. However the position after BCs getting
             applied will be passed via GlobalVars.h.
*/
#include <stdio.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

int
ApplyBC (int k, int IE)
{
  /*!! Declaring external variables */
extern double *rx[NSpecies];
extern double *ry[NSpecies];
extern double *rz[NSpecies];
extern int *ixcount[NSpecies];
extern int *iycount[NSpecies];
extern int *izcount[NSpecies];
extern double Lx, Ly, Lz;

  if (rx[k][IE] >= Lx)
    {
      rx[k][IE] = rx[k][IE] - Lx;
      ixcount[k][IE] = ixcount[k][IE] + 1;
    }
  if (rx[k][IE] < minx)
    {
      rx[k][IE] = rx[k][IE] + Lx;
      ixcount[k][IE] = ixcount[k][IE] - 1;
    }
  if (ry[k][IE] >= Ly)
    {
      ry[k][IE] = ry[k][IE] - Ly;
      iycount[k][IE] = iycount[k][IE] + 1;
    }
  if (ry[k][IE] < miny)
    {
      ry[k][IE] = ry[k][IE] + Ly;
      iycount[k][IE] = iycount[k][IE] - 1;
    }
  if (rz[k][IE] >= Lz)
    {
      rz[k][IE] = rz[k][IE] - Lz;
      izcount[k][IE] = izcount[k][IE] + 1;
    }
  if (rz[k][IE] < minz)
    {
      if (izcount[k][IE] == 0)
        {
          return (1);
        }
      rz[k][IE] = rz[k][IE] + Lz;
      izcount[k][IE] = izcount[k][IE] - 1;
    }
  return (0);
}
