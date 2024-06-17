#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_math.h"
const gsl_rng_type *TBKL;
const gsl_rng_type *Trdir2d;
const gsl_rng_type *Trgauss;
gsl_rng *urn_BKL;
gsl_rng *rdir2d;
gsl_rng *rgauss;
double rng_min;
double rng_diff;
