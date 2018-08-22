/*
 * findStates.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 *  Modified: AGabor 2018
 *  
 *  	originally nodes which had input edges are detected as states. This cause 
 *  	error, if a node has no incoming edge but not a stimuli node. 
 *  	According to the new implementation, we detect the nodes which are not stimulated 
 *  	and we define them as the states of the ODE sytem.
 */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


int *findStates(int nNodes, int nStimuli, int* indexStimuli )
{
	int *stateVec=malloc(nNodes*sizeof(int));
	int j;
	
	for(j=0; j<nNodes; j++)
	{
		stateVec[j]=1;
	}
	for(j=0; j<nStimuli; j++)
	{
		stateVec[indexStimuli[j]-1]=0;  //indexStim[j]-1: indexStim starts from 1
	}
	
	return(stateVec);
}
