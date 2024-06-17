/*!
ReactDiff is written to study the diffusion limited reactions of hydrogen
isotopes, in Graphite by Manoj Warrier (manoj.Warrier - - gmail.com) and Abha
Rai (abra - - ipp.mpg.de). It uses Kinetic Monte Carlo (KMC), specifically,
the Bortz-Kalos-Liebowitz (BKL) algorithm [A. B. Bortz, M. H. Kalos,
J. L. Lebowitz, J. Comp. Phys., 17 (1975) 10 as elucidated in Kai Nordlund's
notes on Monte Carlo [http://beam.acclab.helsinki.fi/~knordlun/mc/mc8nc.pdf].
Duane Johnson's notes on the WWW
[http://www.mcc.uiuc.edu/summerschool/2001/duane%20johnson/kmc_ss.pdf]
The theoritical foundations of dynamical Monte-Carlo is available at
[http://www.mcc.uiuc.edu/summerschool/2001/duane%20johnson/Fichthorn.pdf]

ReactDiff is written in C (gcc) and uses the GNU Scientific Library available
at [http://www.gnu.org/software/gsl/]. GSL is available under the GNU public
license.

Units used: All units are in MKS system unless otherwise specified. This is an
adapted version for the Summer School on hydrogen surface interactions to be
held between June 8^th to 13^th 2008 at Alenya, France.

Copyright (C) <2008>  <Manoj Warrier, Abha Rai>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

/*!
   Function "void ReadInputs ()" reads in inputs from the file ReactDiff.inp
   and is described in the file readinputs.c
*/
void ReadInputs ();

/*!
   Function "void PrintInputs ()" prints out the inputs read in and is a check
   for the transfer of input global variables across functions and is described
   in the file printinputs.c
*/
void PrintInputs ();

/*!
   Function "void InitGeometry ()" creates a surface and is described in file
   initgeometry.c
*/
void InitGeometry ();

/*!
   Function "void InitParams ()" initialises the various simulation parameters,
   allocates arrays, etc and is described in the file initparams.c
*/
void InitParams ();

/*!
   Function "void KMCDriver ()" is the Kinetic Monte Carlo Driver routine and
   is described in the file kmcdriver.c
*/
void KMCDriver ();

int
main (void)
{
  /*! Reading inputs */
  ReadInputs ();

  /*! Printing Inputs */
  PrintInputs ();

  /*! Initialising Parameters */
  InitParams ();

  /*! Generating the surface */
  InitGeometry ();

  /*! Executing the BKL Algorithm */
  KMCDriver ();

  exit (0);
}
