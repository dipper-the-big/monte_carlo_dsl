/*!
     void CheckRecomb (int TheSpecie, int ThePart) ->
       Function to check for recombination events.

     Inputs:
       int TheSpecies -> Species for which recombination possibilities
                         are to be checked.
       int ThePart -> The specific particle of species TheSpecies for
                      which the recombination possibility is checked.

     Outputs:
       none. Note that it is a void function.

     Variables set:
       If recombination occurs,
       * resultant recombined specie if any is created (i.e. its position
         and number of particles are set).
       * TheSpecie is reset (ThePart recombines and is no longer part of
         TheSpecie, particles after ThePart are shifted down by a unit
         index).
       
*/

#include <stdio.h>
#include "GlobalVars.h"
#include "GlobalGSL.h"

/*!! Function interactFunc (int ievent, int ichosen, int intflag,
     int iinteract, int iproduced) described in file
     interactfunc.c */
void InteractFunc (int ievent, int ichosen, int intflag, int iinteract,
    int iproduced);

/*!! Function FindInteract (int ievent, int ichosen, int iinteract)
     described in file findinteract.c */
void FindInteract (int ievent, int ichosen, int iinteract);

void CheckRecomb (int TheSpecie, int ThePart)
{
  /*!! Declaring external variables */
  extern double InteractDist[NRecomb];
  extern int MinDistIndex, InteractSpecie;

  /*!! Declaring local variables */
  int iproduced;

  #ifdef DBuG2
  printf("Entering CheckRecomb\n");
  #endif

  MinDistIndex = -1;
  InteractSpecie = -1;
  switch (TheSpecie)
    {
    case 0:
      MinDist = InteractDist[HHrecomb];
      FindInteract (ThePart, TheSpecie, HSolute);
      switch (InteractSpecie)
	{
	case -1:		/*!! No recombination */
	  break;
	case 0:		/*!! H2 solute is produced */
	  iproduced = 1;
	  InteractFunc (ThePart, TheSpecie, MinDistIndex,
			InteractSpecie, iproduced);
          #ifdef DBuG3
	  printf ("Recombined, ThePart = %d; TheSpecie = %d\n",ThePart,TheSpecie);
          #endif /* DBuG3 */
	  break;
	default:
	  printf ("Error in bklalgorithm.c default\n");
	  printf ("InteractSpecie = %d\n", InteractSpecie);
	  printf ("TheSpecie, ThePart = %d %d\n", TheSpecie, ThePart);
	  exit (0);
	}
      break;
    default:
      printf ("Error in bklalgorithm.c default\n");
      printf ("InteractSpecie = %d\n", InteractSpecie);
      printf ("TheSpecie, ThePart = %d\n", TheSpecie);
      exit (0);
    }
  #ifdef DBuG2
  printf("Exiting CheckRecomb\n");
  #endif
}
