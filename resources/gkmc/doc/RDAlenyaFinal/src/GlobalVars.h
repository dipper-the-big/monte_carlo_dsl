#define DBuG1
#undef DBuG1
#define DBuG2
#undef DBuG2
#define DBuG3
#undef DBuG3
#define DBuG4
#undef DBuG4
#define IncludeReactions
//#undef IncludeReactions

/*
Specie 0 HSolute
Specie 1 H2Solute
*/
#define NSpecies 2
#define HSolute 0
#define H2Solute 1

/*!!
  2 types of jumps for H-solute: desorption and surface diffusion,
  1 type of jump for H2: desorption.

  The various kinds of reactions that occur NRecomb = 1 (H-H Recombination)

  The various kinds of desorptions that occur NDesorb = 2 
    H and H_2 desorption)
*/
#define NRecomb 1
#define HHrecomb 0
#define NDesorb 2

/*!! DiGGeometry Variables */
#define minx 0.0
#define miny 0.0
#define minz 0.0

/* Setting the physical constants used in the simulation */
#define eVToJoules 1.6022e-19
#define K_Boltzmann 1.3807e-23
#define Sqrt3 1.73205

/*!! BKLKMC variables */
int MinDistIndex, InteractSpecie;
double MinDist;

/*!! DiGDiffcoeff variables */
double Lx, Ly, Lz, CelSize, dt, dsqrbydt, dxsqrbydt, dysqrbydt, dzsqrbydt,
       dsqrbydtsum, dxsqrbydtsum, dysqrbydtsum, dzsqrbydtsum;
int *OccFlag, nx, ny, nz;

/*!! DiGInputs Variables */
int NParts[NSpecies], NJumps[NSpecies], MaxSteps, NTrials, NPosition,
    NTransient, NPartsINI[NSpecies];
double *JmpFreq[NSpecies], *JmpDist[NSpecies], *Em[NSpecies], 
    InteractDist[NRecomb], GammaP[NSpecies], StickCoef[NSpecies], Temp;

/*!! DiGStatistics Variables */
double *Frac[NSpecies], *FracKMC[NSpecies], *FracDesorb[NSpecies],
    *DiffCoeffx[NSpecies], *DiffCoeffy[NSpecies], *DiffCoeffz[NSpecies],
    *DiffCoeff[NSpecies], *Zeit;

/*!! InitialisedParams Variables */
double DeltaT, AddPDeltaT, *rx[NSpecies], *ry[NSpecies], *rz[NSpecies],
    *rxi[NSpecies], *ryi[NSpecies], *rzi[NSpecies], *ResideTi[NSpecies],
    *ResideTf[NSpecies], *JmpRate[NSpecies], ParticleRateSum[NSpecies],
    CurrentTime;
int *ixcount[NSpecies], *iycount[NSpecies], *izcount[NSpecies],
    *icurrent[NSpecies], *scurrent[NSpecies], *ixcnttmp, *iycnttmp,
    *izcnttmp, KMCIndex[NSpecies], NPartsKMC[NSpecies],
    NPartsDesorb[NSpecies];
