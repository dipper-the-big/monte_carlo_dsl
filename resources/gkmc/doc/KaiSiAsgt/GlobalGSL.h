#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>

const gsl_rng_type *T;
gsl_rng *r_urn;
gsl_rng *r_rdir3d;
gsl_rng *r_gauss;

double rng_max, rng_min, rng_diff;
