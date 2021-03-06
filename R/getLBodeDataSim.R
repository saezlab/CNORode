#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
#' @export
getLBodeDataSim<-function
(
		cnolist,				    model,					ode_parameters=NULL,
		indices=NULL,			    timeSignals=NULL,		time=1,					
		verbose=0,					transfer_function=3,	reltol=1e-4,			
		atol=1e-3,					maxStepSize=Inf,		maxNumSteps=100000,		
		maxErrTestsFails=50
)
{

    if (class(cnolist)=="CNOlist"){
    	cnolist = compatCNOlist(cnolist)
    }
	
	adjMat=incidence2Adjacency(model)
	
	if(is.null(indices)){
		indices <- indexFinder(cnolist,model,verbose=FALSE)
	}
	
	if(is.null(ode_parameters)){
		ode_parameters = createLBodeContPars(model)
	}
	
	if(!is.null(timeSignals)){
		cnolist$timeSignals = timeSignals
	}
	
	sim_function = getLBodeSimFunction(cnolist, model, adjMat,
			indices1=indices, odeParameters1=ode_parameters$parValues, time1=time, verbose1=verbose,
			transfer_function1=transfer_function, reltol1=reltol, atol1=atol, maxStepSize1=maxStepSize,
			maxNumSteps1=maxNumSteps, maxErrTestsFails1=maxErrTestsFails)
	
	sim=sim_function(cnolist,model, odeParameters=ode_parameters$parValues)
	
	sim=lapply(sim,function(x) x[,indices$signals, drop=FALSE])
	
	return(sim);
}
