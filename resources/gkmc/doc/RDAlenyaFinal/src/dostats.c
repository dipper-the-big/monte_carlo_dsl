/*!! DoStats (char *VarDesc, double *StatVar, int Num):
           Void function which finds the mean, variance and error
           of the passed array 'statvar' which has 'Num' dimension.
           The description of the variable is passed via the character
           variable VarDesc.

   Finding the statistical error assuming Gaussian statistics
   Method followed as described in:
   http://beam.helsinki.fi/~knordlun/mc/mc5.ps, page 8.

   Using function gsl_stats_variance, which returns
   variance = (1/(N-1)) \sum (x_i - \hat\mu),
   where \hat\mu is the mean of x_i (i=1,....,N).

   Then the error of the average is calculated using the formula
   error = \sqrt(variance/N)
*/
#include <stdio.h>
#include <math.h>
#include "GlobalGSL.h"

void
DoStats (char *VarDesc, double *StatVar, int Num)
{
  /*!! No external variables used by this function */

  /*!! Declaring local variables */
  static double MeanVal, VarVal, ErrVal;

  MeanVal = gsl_stats_mean (StatVar, 1, Num);
  VarVal = gsl_stats_variance_m (StatVar, 1, Num, MeanVal);
  ErrVal = sqrt (VarVal / Num);
  printf ("%s = %e %e\n", VarDesc, MeanVal, ErrVal);
}
