/*!  Function KMCDriver (): The function that is the main driver which:
     (1) Sets up variables and proceedures for statistical analysis,
     (2) Carries out the Kinetic Monte Carlo (KMC) time stepping till
         MaxSteps specified as input in readinput.c. MaxSteps should be
         large enough for detailed balance to have set in.

   Written by Manoj Warrier (manoj@ipr.res.in) on 28-10-2005.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GlobalVars.h"

/*!! "void function AllocStatVarMem ()" described in the file
     allocstatvarmem.c */
void AllocStatVarMem ();

/*!! "void function ReallocMem ()" described in the file
     reallocmem.c */
void ReallocMem ();

/*!! Function "void InitParticles (int *NP)" which initialises the particles,
     is described in the file initparticles.c */
void InitParticles (int NP[NSpecies]);

/*!! "void function BKLAlgorithm ()" described in the file
     bklalgorithm.c */
void BKLAlgorithm ();

/*!! "void function PosOut ()" described in the file posout.c */
void PosOut ();

/*!! "void function PostProcess (M)" described in the file
     postprocess.c */
void PostProcess (int M);

/*!! "void function StatPostProc ()" described in the file
     statpostproc.c */
void StatPostProc ();

double round (double x);

void
KMCDriver ()
{
  /*!! Declaring variables */

  /*!! Declaring external variables */
  extern int MaxSteps, NTrials, NPartsKMC[NSpecies], NPartsDesorb[NSpecies],
    NParts[NSpecies], NParts[NSpecies], NPosition;
  extern double *Zeit, DeltaT, GammaP[NSpecies], StickCoef[NSpecies], Lx, Ly;

  /*!! Declaring variables set in this function */

  /*!! Declaring internal variables */
  int k, M, DeltaNParts[NSpecies], ICount, AddNew, ChkPrts;
  double PrtBal;

  FILE *PrtBalance;

  #ifdef DBuG1
  printf ("Starting KMCDriver()\n");
  #endif /* DBuG1 */

  /*!! Allocating memory for variables used in statistical error evaluation */
  AllocStatVarMem ();
  printf ("Allocated memory to statistical variables\n");

  if ((PrtBalance = fopen ("../out/ParticleBalance.log","w")) == NULL)
    {
      printf ("Could not open ../out/ParticleBalance.log\n");
      exit(0);
    }
  PrtBal = 0.0;

  /*!! Doing the main loop over NTrials for statistics */
  for (M = 0; M < NTrials; M++)
    {
      ICount = 1;
      Zeit[M] = 0.0;
      AddPDeltaT = 0.0;

      for (k = 0; k < NSpecies; k++)
        {
	  NParts[k] = NPartsINI[k];
          NPartsKMC[k] = NParts[k];
	  NPartsDesorb[k] = 0;
          DeltaNParts[k] = -1;
        }

      /*!! Initialising the particles */
      InitParticles (DeltaNParts);

      while (ICount < MaxSteps)
        {
          /* Calculating the newly arriving particles */
          AddNew = 0;
          ChkPrts = 0;
          for (k = 0; k < NSpecies; k++)
            {
              /*!
                 Calculating increase in number of particles of species k
                 in time DeltaT and reallocating memory.
              */
              DeltaNParts[k] = (int) round (GammaP[k] * StickCoef[k] *
		AddPDeltaT * Lx * Ly);
              NPartsKMC[k] = NPartsKMC[k] + DeltaNParts[k];
              NParts[k] = NPartsKMC[k] + NPartsDesorb[k];
              if (DeltaNParts[k] != 0)
		{
		  AddNew = 1;
		  AddPDeltaT = 0;
		}
              if (NPartsKMC[k] > 0) ChkPrts = 1;
              #ifdef DBuG1
                printf("DeltaNParts[%d] = %d\n",k,DeltaNParts[k]);
                printf("GammaP[%d] = %e; StickCoef[%d] = %e\n",
                      k, GammaP[k], k, StickCoef[k]);
                printf("lx = %e, ly = %e\n", Lx, Ly);
                printf("NParts[%d] = %d; NPartsKMC[%d] = %d; NPartsDesorb[%d] = %d\n",
                      k, NParts[k], k, NPartsKMC[k], k, NPartsDesorb[k]);
              #endif /* DBuG1 */
            }

          /*!! Initialising the particles */
          if (AddNew == 1)
            {
              ReallocMem ();
              InitParticles (DeltaNParts);
	    }

	  /*!! Running the KMC time step */
	  if (ChkPrts > 0) BKLAlgorithm ();
	  else
            {
              ICount = MaxSteps;
              printf("No more particles to run BKL Algorithm on\n");
	      exit(0);
            }

//          if ((ICount%NPosition)==0) PosOut();

	  /* Defining CurrentTime so that Zeit and M need not
             be passed to each function needing it */
	  Zeit[M] = Zeit[M] + DeltaT;
	  AddPDeltaT = AddPDeltaT + DeltaT;
          CurrentTime = Zeit[M];

          ICount = ICount + 1;

          if ((ICount%NPosition)==0)
            {
              PrtBal = (double) (NParts[0] + 2.0 * NParts[1] -
	        NPartsINI[0] - NPartsINI[1]) / Zeit[M] / Lx / Ly;
              fprintf(PrtBalance,"%e %d %d %d %d %d %d %e\n",
                Zeit[M], NParts[0], NPartsKMC[0], NPartsDesorb[0],
                NParts[1], NPartsKMC[1], NPartsDesorb[1], PrtBal);
            }
	  fflush (PrtBalance);
          #ifdef DBuG1
          printf("Zeit[%d] ,DeltaT = %e %e\n",M,Zeit[M],DeltaT);
          printf("ICount = %d\n",ICount);
          #endif /* DBuG1 */
	}

      /*!! Post processing of each parameter of interest to get values
         for the M^th run of the KMC.
       */
      PostProcess (M);
    }

  fclose (PrtBalance);

  StatPostProc ();

  #ifdef DBuG1
    printf("Exiting KMCDriver ()\n");
  #endif /* DBuG1 */
}
