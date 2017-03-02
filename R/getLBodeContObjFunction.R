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
# $Id: getLBodeContObjFunction.R 3184 2013-01-21 13:50:31Z cokelaer $
getLBodeContObjFunction<-
  function
(
  cnolist,                model,                    ode_parameters,
  indices=NULL,           time=1,                   verbose=0,
  transfer_function=3,    reltol=1e-4,            atol=1e-3,
  maxStepSize=Inf,        maxNumSteps=100000,        maxErrTestsFails=50,
  nan_fac=1, lambda=0, boot_seed=sample(1:10000,1)
)
  {
  
  if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
  adjMatrix=incidence2Adjacency(model);
  
  if(is.null(indices))indices=indexFinder(cnolist,model,verbose=FALSE);
  
  sim_function=getLBodeSimFunction(cnolist,model,adjMatrix1=adjMatrix,
                                   indices1=indices, odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
                                   transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
                                   maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails)
  
  logic_ode_continuous_objective_function<-function
  (
    x,                            cnolist1=cnolist,        model1=model,
    adjMatrix1=adjMatrix,        indices1=indices,        ode_parameters1=ode_parameters,
    sim_function1=sim_function,        nan_fac1=nan_fac
  )
  {
    
    
    ode_parameters1$parValues[ode_parameters1$index_opt_pars]=x;
    sim=sim_function1(cnolist1,model1,ode_parameters1$parValues);
    temp=sim;
    sim=list();
    for(i in 1:length(cnolist$timeSignals)){
      sim[[i]]=temp[[i]];
    }
    
    sim<-as.vector(unlist(lapply(sim,function(x) x[,indices1$signals])));
    measured_values<-as.vector(unlist(lapply(cnolist1$valueSignals,function(x)x)));
    
    res_boot<-(sim-measured_values)^2
    set.seed(boot_seed)
    res_boot<-sample(res_boot, length(res_boot), replace=T)
    error<-sum(res_boot, na.rm = T)
    
    NaNs_sim=which(is.na(sim));
    
    # print(sample(1:5, 5, replace=T))
    
    lambda_k<-lambda
    lambda_tau<-0
    
    NApenalty<-length(NaNs_sim)*nan_fac1
    SSpenalty<-10*sum(temp[[length(cnolist$timeSignals)+1]]^2)
    SSbasalPenalty<-1000*sum(temp[[length(cnolist$timeSignals)+2]][1,]^2)
    L1reg_k<-lambda_k*(sum(abs(ode_parameters1$parValues[ode_parameters1$index_k])))
    L1reg_tau<-lambda_tau*(sum(abs(ode_parameters1$parValues[ode_parameters1$index_tau])))
    L1reg<-L1reg_k+L1reg_tau
    
    res=error+length(NaNs_sim)*nan_fac1+SSpenalty+SSbasalPenalty+L1reg_k;
    
    cat("res =", res, "\n",
        "NA penalty = ", NApenalty, "\n",
        "error =", error, "\n",
        "steady state penalty, ", SSpenalty, "\n",
        "control steady state penalty", SSbasalPenalty, "\n",
        "L1-reg: lambda =", lambda, "- P penalty (only tau)", L1reg_tau, "\n",
        "L1-reg: lambda =", lambda, "- P penalty (only k)", L1reg_k, "\n",
        "L1-reg: lambda =", lambda, "- P penalty (all)", L1reg, "\n",
        "--------\n")

    if(is.nan(res) || is.na(res))res=1e10;
    
    return(res);
  }
  return(logic_ode_continuous_objective_function);
}
