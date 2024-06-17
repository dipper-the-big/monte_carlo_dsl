/* Program SiAsgt is an example program to study Kinetic Monte Carlo,
   specifically the Bortz-Kalos-Liebowitz (BKL) algorithm described
   in Dr. Kai Nordlund's notes on Kinetic Monte Carlo available at
   http://beam.helsinki.fi/~knordlun/mc/mc8.ps

   This program is an assignment in response to Exercise4
   http://beam.helsinki.fi/~knordlun/mc/exc4.ps

   It is written in C and uses the GNU Scientific Library available
   at http://sources.redhat.com/gsl  and aims to be written so as to
   as to make it extensible to solve more complex problems using KMC.
   See the file Readme for more information.

   Written by Manoj (manoj.warier@ipp.mpg.de) on 14-10-2002.
*/
#include <stdio.h>
#include <stdlib.h>

#include "GlobalVars.h"

void AllocArrays ();
void ReadInputs ();
void PrintInputs ();
void InitParams ();
void BKLAlgorithm ();

int
main (void)
{
  AllocArrays ();
  ReadInputs ();
  PrintInputs ();
  InitParams ();
  BKLAlgorithm ();
  exit (0);
}
