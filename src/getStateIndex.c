#include <stdio.h>
#include <stdlib.h>
#define DEBUG 0
/*
int *getStateIndex(int **adjMatrix, int n)
{
    int* indexVec=malloc(n*sizeof(int));
    int i,j;
    int stateNumber=0;
    int first;
    for (j = 0; j <n; j++)
    {
    	indexVec[j]=-1;
    	first=1;
        for(i=0;i<n;i++)
        {
            if(adjMatrix[i][j] && first)
            {
            	indexVec[j]=stateNumber++;
            	first=0;
            }
        }
     }
    
    printf("State Index\n");
    for (j = 0; j <n; j++)
    {
    	printf("species %d \t\t index %d\n",j,indexVec[j]);
    }
    

   return(indexVec);
}
*/

int *getStateIndex(int nNodes, int nStimuli, int* indexStimuli )
{
	int *indexVec=malloc(nNodes*sizeof(int));
	int j,i;
	int stateNumber=0;
	
	for(i=0; i<nNodes; i++)
	{
		indexVec[i] = stateNumber;  // sets state index
		
		// set to -1 if it is a stimuli and not a state
		for(j=0; j<nStimuli; j++)
		{
			if(i == indexStimuli[j]-1 ) //indexStim[j]-1: indexStim starts from 1
				indexVec[i]=-1;  
		}
	
		if(indexVec[i]!=-1) ++stateNumber;
		
	}
	if(DEBUG){
		printf("State Index\n");
		for (j = 0; j < nNodes; j++)
		{
			printf("species %d \t\t index %d\n",j,indexVec[j]);
		}
	}
	
	return(indexVec);
}
