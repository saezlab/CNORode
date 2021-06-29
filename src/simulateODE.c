#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/* Sundials Header Files */

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>

#include "./include/nvector/nvector_serial.h"
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

int rhsODEF(double t, double* y, double* ydot, void *data);

#include "CNOStructure.h"

#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

#define DEBUG 0

typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data);

static int check_flag(void *flagvalue, char *funcname, int opt,int verbose);

int simulateODE
(
  	CNOStructure* data,		int exp_num, 			int verbose,
		double reltol,			double atol,			double maxStepSize,
		int maxNumSteps,		int maxErrTestFails
)
{
	int i,j,neq,counter,flag;
	realtype tout, ti, tf;
	N_Vector y;
	void *cvode_mem;
	double* yf;
  double* y0;
  double* ydotf;
  double* ydot0;
    
    
	cvode_mem = NULL;
	y = NULL;

	neq=(*data).nStates;

    y = N_VNew_Serial(neq);
	if (check_flag((void *)y, "N_VNew_Serial", 0,verbose))
    {
		//if(verbose)printf("\nSolver failed in N_VNew_Serial(neq). . .\n");
		return(0);
	}

    /* Initialize y */
	for(i=0; i<(*data).nRows; i++)
	{
		(*data).state_array[i] = 0.5;
		(*data).inhibitor_array[i]=0;
	}

	for(i=0; i<(*data).nRows; ++i)
	{
		if((*data).isState[i])
		{
			for (j = 0; j < (*data).nSignals; j++)
			{
			//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexSignals[j]==i+1)
				{
					(*data).state_array[i] = (*data).valueSignals[exp_num][j];
				}
			}

			for (j = 0; j < (*data).nInhibitors; j++)
			{
				//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexInhibitors[j]==i+1)
				{
					(*data).inhibitor_array[i] = (*data).valueInhibitors[exp_num][j];
				}
			}
		}
		else
		{
			for (j = 0; j < (*data).nStimuli; j++)
			{
				//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexStimuli[j]==i+1)
				{
					(*data).state_array[i] = (*data).valueStimuli[exp_num][j];
				}
			}
		}
		
		
	}
	
	if(DEBUG){
		printf("\nstate_array:\n\t");
		for(i=0; i<(*data).nRows; ++i)
		{
			printf("%f ",(*data).state_array[i]);
		}
		printf("\n");
	}
	
	counter=0;
	for(i=0; i<(*data).nRows; i++)
	{
		if((*data).isState[i])
		{
			Ith(y,(*data).state_index[i]) =
					(*data).state_array[i];
			(*data).sim_results[exp_num][0][(*data).state_index[i]]=
					(*data).state_array[i];
			/*
			if(verbose)
			{
				printf("species number initial value %d",i);
				printf("\t%f\n",Ith(y,(*data).state_index[i]));
			}
			*/
		}
	}

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if(check_flag((void *)cvode_mem, "CVodeCreate", 0,verbose)) {
		//if(verbose)printf("\nSolver failed in CVodeCreate(CV_BDF, CV_NEWTON) . . .\n");
		N_VDestroy_Serial(y);
		return(0);
	}

	 ti=(*data).timeSignals[0];
	 tf=(*data).timeSignals[(*data).nTimes-1];
	  flag = CVodeInit(cvode_mem,*rhsODE, ti, y);
	 if (check_flag(&flag, "CVodeMalloc", 1,verbose))
	 {
		//if(verbose)printf("\nSolver failed in CVodeInit(cvode_mem,*rhsODE, ti, y). . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
	 	return(0);
	 }

	//if(verbose)printf("Solver Memory Allocated\n");

	/* Set f_data */
	 flag = CVodeSetUserData(cvode_mem, data);
	 if(check_flag(&flag, "CVodeSetFdata", 1,verbose))
	 {
		// if(verbose)printf("\nSolver failed in flag = CVodeSetUserData(cvode_mem, data). . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		 return(0);
	 }

	 flag = CVodeSStolerances(cvode_mem,(realtype)reltol,(realtype)atol);
	  if(check_flag(&flag, "CVodeSStolerances", 1,verbose)) return(1);

	  if(!verbose)
	  {
		  flag = CVodeSetErrFile(cvode_mem, NULL);
	  }

	  flag = CVDense(cvode_mem, neq);
	  if (check_flag(&flag, "CVDense", 1,verbose)) return(1);
	//if(verbose)printf("CVDENSE Solver Initiated\n");

	 /* Set maxnumsteps */
	 flag = CVodeSetMaxNumSteps(cvode_mem, maxNumSteps);
	 if(check_flag(&flag, "CVodeSetMaxNumSteps", 1,verbose))
	 {
		//if(verbose)printf("\nSolver failed in CVodeSetMaxNumSteps(cvode_mem, maxNumSteps). . .\n");
		return(0);
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
	 }
	// if(verbose)printf("Max number of steps: %i\n", maxNumSteps);

    CVodeSetMaxStep(cvode_mem,(realtype)maxStepSize);
    CVodeSetMaxErrTestFails(cvode_mem, maxErrTestFails);
	 

	 if(DEBUG){
	 	printf("checkpoint 1: exp_num = %d\n",exp_num);
	 	printf(" i for nTimes\n");
	 }
    for (i = 1; i < (*data).nTimes; ++i)
    {
    	if(DEBUG){
    		printf("%d/%d\n",i,(*data).nTimes);
    		printf("call CVode\n");
    	}
    	
    	tout=(*data).timeSignals[i];
    	flag = CVode(cvode_mem, tout, y, &tf, CV_NORMAL);
    	
    	if(DEBUG){
    		printf("CVode passed\n");
    		printf("read out solutions:\n");
    	}
    	
    	for (j = 0; j < (*data).nStates; j++)
    	{
    		if(DEBUG){
    			printf("state %d/%d\n",j,(*data).nStates);
    			printf("call Ith(y,j)\n");
    			printf("%f\n",(double)Ith(y,j));
    			printf("%f\n",(double)Ith(y,j));
    			printf("(*data).sim_results[exp_num][i][j]= (double)Ith(y,j);\n");
    			printf("(*data).sim_results[%d][%d][%d]= (double)Ith(y,%d);\n", exp_num,i,j,j);
    		}
    		(*data).sim_results[exp_num][i][j]= (double)Ith(y,j);
    		if(DEBUG){
    			printf("assignement passed\n");
    		}
    	//	if(verbose)
    	}

    	if (check_flag(&flag, "CVode", 1,verbose))
    	{
    	//	if(verbose)fprintf(stderr,"\nSolver failed at flag = CVode(cvode_mem, tout, y, &tf, CV_NORMAL);. . .\n");
    		N_VDestroy_Serial(y);
    		CVodeFree(&cvode_mem);
    		return(0);
    	}
    	//if(verbose)printf("\n");
    }
    if(DEBUG){
    	printf("checkpoint 3\n");
    }
    yf=(double*)malloc((*data).nStates*sizeof(double));
    ydotf=(double*)malloc((*data).nStates*sizeof(double));
    if(DEBUG){
    	printf("checkpoint 4\n");
    }
    for (j = 0; j < (*data).nStates; j++){
    	yf[j]= (*data).sim_results[exp_num][(*data).nTimes-1][j];
    }
    if(DEBUG){
    	printf("checkpoint 5\n");
    }
    rhsODEF(tout,yf,ydotf, data);
    if(DEBUG){
    	printf("checkpoint 6\n");
    }
   y0=(double*)malloc((*data).nStates*sizeof(double));
   ydot0=(double*)malloc((*data).nStates*sizeof(double));
   for (j = 0; j < (*data).nStates; j++) y0[j]=(*data).sim_results[exp_num][0][j];
   rhsODEF(tout,y0,ydot0, data);

    
    for (j = 0; j < (*data).nStates; j++) (*data).ydotf[exp_num][j]=ydotf[j];
    for (j = 0; j < (*data).nStates; j++) (*data).ydot0[exp_num][j]=ydot0[j];

        

    N_VDestroy_Serial(y);
    /* Free integrator memory */
    
    free(ydotf);
    free(yf);
    
    CVodeFree(&cvode_mem);

    return(1);
}

static int check_flag(void *flagvalue, char *funcname, int opt,int verbose)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		//if(verbose)fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1); }

	/* Check if flag < 0 */
	else if (opt == 1)
	{
		errflag = (int *) flagvalue;
		if (*errflag < 0)
		{	//if(verbose)fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",funcname, *errflag);
			return(1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL)
	{
		//if(verbose)fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1);
	}

	return(0);
}


