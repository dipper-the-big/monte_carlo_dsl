/* Function BKLAlgorithm (): The main function to carry out the
   Kinetic Monte Carlo (or residence time algorithm as recommended
   in Kai's notes). See page 5 of chapter-8 in Kai's notes to
   understand the comments.

   Written by Manoj Warrier (manoj.warrier@ipp.mpg.de) on 16-10-2002
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GlobalVars.h"
#include "GlobalGSL.h"

void InitMaterial ();

void
BKLAlgorithm ()
{
  /* Declaring variables */
  extern double *Em;
  extern double *JmpDist;
  extern double *JmpFreq;
  extern double *JmpRate;
  extern double *InteractDist;
  extern double *rxi;
  extern double *ryi;
  extern double *rzi;
  extern double *rxv;
  extern double *ryv;
  extern double *rzv;
  extern double *rxr;
  extern double *ryr;
  extern double *rzr;
  extern double *Ri;
  extern double Temp;
  extern int NTrials;
  extern int *NParts;
  extern gsl_rng *r_urn;
  extern gsl_rng *r_rdir3d;
  extern double rng_min;
  extern double rng_diff;
  extern int MaxSteps;

  int i, j, k, M, IEvent, NTotal, ICount;
  double RiCum, urn, DeltaT, dist, dx, dy, dz, rx, ry, rz, *Zeit;
  double IJumps, VJumps, *FracI, *FracV, IniI, IniV, *IJmpByVJmp;
  double MeanFracI, MeanFracV, MeanIJmpByVJmp, MeanZeit, FindIEvent;
  double VarFracI, VarFracV, VarIJmpByVJmp, VarZeit;
  double ErrFracI, ErrFracV, ErrIJmpByVJmp, ErrZeit;

  printf ("Starting BKL Algorithm\n");

  /* Allocating memory for variables used in statistical error
     evaluation */
  if ((FracI = (double *) malloc (NTrials * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for FracI\n");
      exit (1);
    }
  if ((FracV = (double *) malloc (NTrials * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for FracV\n");
      exit (1);
    }
  if ((IJmpByVJmp = (double *) malloc (NTrials * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for IJmpByVJmp\n");
      exit (1);
    }
  if ((Zeit = (double *) malloc (NTrials * sizeof (double))) == NULL)
    {
      printf ("Not enough memory to allocate space for Zeit\n");
      exit (1);
    }

  /* Initialising variables */
  IniI = (double) NParts[0];
  IniV = (double) NParts[1];

  for (M = 0; M < NTrials; M++)
    {
      ICount = 0;
      IJumps = 0.0;
      VJumps = 0.0;
      FracI[M] = 0.0;
      FracV[M] = 0.0;
      IJmpByVJmp[M] = 0.0;
      NParts[0] = (int) IniI;
      NParts[1] = (int) IniV;
      NParts[2] = 0;

      InitMaterial ();

      /* Step (0) : Set time = 0 */
      Zeit[M] = 0.0;

      while (ICount <= MaxSteps)
	{
      /* Step (1) : Form a list of all possible transitions:
         Possible transitions : interstitials can jump
         vacancies can jump
         Each vacancy or interstitial jumping is a transition,
         therefore total # of transitions = NParts[0] + NParts[1]
         However since there are only 2 species, only 2 jump rates
         need be evaluated
      */
	  for (i = 0; i < NSpecies; i++)
	    {
              /* The 8.61745e-5 below converts Kelvins to eV */
	      JmpRate[i] = JmpFreq[i] * exp (-Em[i] / 8.629375e-5 / Temp);
	    }

      /* Step (2) : Calculate the cumulative function ....
         This is where the sum over all possible transitions of
         Step(1) comes. We calculate a Cumulative function Ri.
       */
	  NTotal = 0;
	  for (i = 0; i < NSpecies; i++)
	    {
	      NTotal = NTotal + NParts[i];
	    }

	  RiCum = 0.0;
	  j = 0;
          Ri[j] = 0.0e0;

	  for (k = 0; k < NSpecies; k++)
	    {
	      for (i = 0; i < NParts[k]; i++)
		{
		  RiCum = RiCum + JmpRate[k];
		  j = j + 1;
		  Ri[j] = RiCum;
		}
	    }

      /* Step (3) : Get a uniform random number between [0,1] */
	  urn = ((double) gsl_rng_get (r_urn) - rng_min) / rng_diff;

      /* Step (4) : Find out the event to carry out "IEvent" by
         finding the "IEvent" for which R_(i-1) < urn*RiCum < R_i.
         We do it using the successive bisection method.
       */
	  FindIEvent = urn * RiCum;
          for (i=0; i<NTotal; i++)
            {
              if (FindIEvent <= Ri[i+1] && FindIEvent >= Ri[i]) IEvent = i;
            }


	  /* Step (5) : Carry out event "i".
	     First we find out which type of particle does the event occur
	     to. It is not necessary in the Si assignment to do this, but
	     when the jump distances are different, and when we have different
	     position variables for the species we have to do this. We
	     propogate the particle corresponding to event "i" by a distance
	     jmp_dist in a random 3D direction from its current position.
	   */
	  gsl_ran_dir_3d (r_rdir3d, &rx, &ry, &rz);
	  if (IEvent < NParts[0])
	    {
	      j = 0;
	      rxi[IEvent] = rxi[IEvent] + JmpDist[j] * rx;
	      ryi[IEvent] = ryi[IEvent] + JmpDist[j] * ry;
	      rzi[IEvent] = rzi[IEvent] + JmpDist[j] * rz;
	      IJumps = IJumps + 1;
	    }
	  else if ((NParts[0] <= IEvent) && (IEvent < (NParts[0] + NParts[1])))
	    {
	      j = 1;
	      IEvent = IEvent - NParts[0] - 1;
	      rxv[IEvent] = rxv[IEvent] + JmpDist[j] * rx;
	      ryv[IEvent] = ryv[IEvent] + JmpDist[j] * ry;
	      rzv[IEvent] = rzv[IEvent] + JmpDist[j] * rz;
	      VJumps = VJumps + 1;
	    }
	  else
	    {
	      j = 2;
	      printf ("j=2 in BKLAlgorithm!! check error!!\n");
	      printf ("NParts[0], NParts[1], NParts[2], IEvent =%d %d %d %d\n",
                       NParts[0], NParts[1], NParts[2], IEvent);
      	      goto findtime;
              /* This is being done to tackle the special case
                 where all the interstitials and vacancies recombine.
                 ICount is set to MaxSteps+1.
              */
	    }
//          printf("%d %d %d %d %e\n",IEvent,NParts[0],NParts[1],j,RiCum);

	  ICount = ICount + 1;

	  /* Step (6) : Find all transitions and recalculate any rates
	     that might have changed due to the transition
	   */
	  if (j == 0)		/* Random choice is an interstitial */
	    {
	      for (i = 0; i < NParts[1]; i++)	/* loop over vacancies */
		{
		  dist = 0.0;	/* Calculating distance */
		  dx = rxv[i] - rxi[IEvent];
		  dy = ryv[i] - ryi[IEvent];
		  dz = rzv[i] - rzi[IEvent];
		  dist = dx * dx + dy * dy + dz * dz;
		  if (dist <= InteractDist[0]*InteractDist[0])
		    {
		      /* recombination position is the position of vacancy */
		      rxr[NParts[2]] = rxv[i];
		      ryr[NParts[2]] = ryv[i];
		      rzr[NParts[2]] = rzv[i];
		      /* incrementing the number recombined */
		      NParts[2] = NParts[2] + 1;
		      ICount = 0;
		      /* redistributing interstitials */
		      for (k = IEvent; k < (NParts[0] - 1); k++)
			{
			  rxi[k] = rxi[k + 1];
			  ryi[k] = ryi[k + 1];
			  rzi[k] = rzi[k + 1];
			}
		      /* redistributing vacancies */
		      for (k = i; k < (NParts[1] - 1); k++)
			{
			  rxv[k] = rxv[k + 1];
			  ryv[k] = ryv[k + 1];
			  rzv[k] = rzv[k + 1];
			}
		      NParts[0] = NParts[0] - 1;
		      NParts[1] = NParts[1] - 1;
		      goto findtime;
		    }
		}
	    }
	  if (j == 1)		/* Random choice is a Vacancy */
	    {
	      for (i = 0; i < NParts[0]; i++)	/* Loop over interstitials */
		{
		  dist = 0.0;	/* Calculating distance */
		  dx = rxi[i] - rxv[IEvent];
		  dy = ryi[i] - ryv[IEvent];
		  dz = rzi[i] - rzv[IEvent];
		  dist = dx * dx + dy * dy + dz * dz;
		  if (dist <= InteractDist[0]*InteractDist[0])	/* determining recombination */
		    {
		      /* recombination position is the position of vacancy */
		      rxr[NParts[2]] = rxv[IEvent];
		      ryr[NParts[2]] = ryv[IEvent];
		      rzr[NParts[2]] = rzv[IEvent];
		      /* incrementing the number recombined */
		      NParts[2] = NParts[2] + 1;
		      ICount = 0;
		      /* redistributing interstitials */
		      for (k = i; k < (NParts[0] - 1); k++)
			{
			  rxi[k] = rxi[k + 1];
			  ryi[k] = ryi[k + 1];
			  rzi[k] = rzi[k + 1];
			}
		      /* redistributing Vacancies */
		      for (k = IEvent; k < (NParts[0] - 1); k++)
			{
			  rxv[k] = rxv[k + 1];
			  ryv[k] = ryv[k + 1];
			  rzv[k] = rzv[k + 1];
			}
		      NParts[0] = NParts[0] - 1;
		      NParts[1] = NParts[1] - 1;
		      goto findtime;
		    }
		}
	    }

	findtime:
	  /* Step (7) : Get a uniform random number between [0,1] */
	  urn = ((double) gsl_rng_get (r_urn) - rng_min) / rng_diff;

	  /* Step (8) :Updating time */
	  DeltaT = -log (urn) / RiCum;
          if (RiCum == 0.0) /* Special case when all vacancies and
                               interstitials have recombined!! */
            {
              DeltaT = 0.0;
              ICount = MaxSteps +1;
            }
	  Zeit[M] = Zeit[M] + DeltaT;
	}
      FracI[M] = (double) NParts[0] / IniI;
      FracV[M] = (double) NParts[1] / IniV;
      IJmpByVJmp[M] = IJumps / VJumps;

      printf ("%d %e %e %e %e %e\n",
	      M, FracI[M], FracV[M], IJmpByVJmp[M], Zeit[M],
              JmpRate[0]/JmpRate[1]);
    }

  /* Finding the statistical error assuming Gaussian statistics
     Method followed as described in:
          http://beam.helsinki.fi/~knordlun/mc/mc5.ps, page 8.

     Using function gsl_stats_variance, which returns
          variance = (1/(N-1)) \sum ((x_i - \hat\mu)^2),
          where \hat\mu is the mean of x_i (i=1,....,N).

     Then the error of the average is calculated using the formula
          error = \sqrt(variance/N)
  */
  MeanFracI = gsl_stats_mean(FracI, 1, NTrials);
  MeanFracV = gsl_stats_mean(FracV, 1, NTrials);
  MeanIJmpByVJmp = gsl_stats_mean(IJmpByVJmp, 1, NTrials);
  MeanZeit = gsl_stats_mean(Zeit, 1, NTrials);
  VarFracI = gsl_stats_variance_m(FracI, 1, NTrials, MeanFracI);
  VarFracV = gsl_stats_variance_m(FracV, 1, NTrials, MeanFracV);
  VarIJmpByVJmp = gsl_stats_variance_m(IJmpByVJmp, 1, NTrials,
    MeanIJmpByVJmp);
  VarZeit = gsl_stats_variance_m(Zeit, 1, NTrials, MeanZeit);
  ErrFracI = sqrt(VarFracI / NTrials);
  ErrFracV = sqrt(VarFracV / NTrials);
  ErrIJmpByVJmp = sqrt(VarIJmpByVJmp / NTrials);
  ErrZeit = sqrt(VarZeit / NTrials);
  printf ("Fraction of Interstitials = %e; Error = %e\n",
    MeanFracI, ErrFracI);
  printf ("Fraction of Vacancies = %e; Error = %e\n",
    MeanFracV, ErrFracV);
  printf ("Ratio of Interstitial jumps to Vacancy jumps = %e; Error = %e\n",
    MeanIJmpByVJmp, ErrIJmpByVJmp);
  printf ("Simulation time = %e (seconds); Error = %e\n", MeanZeit, ErrZeit);
}
