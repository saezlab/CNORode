#' ## Crossvalidation study
#' crossvalidate
#' 
#' k-fold crossvalidation for Boolean model
#' 
#' Does a k-fold cross-validation for Boolean CellNOpt models. In k-iterations a 
#' fraction of the data is eliminated from the CNOlist. The model is trained on the 
#' remaining data and then the model predicts the held-out data. Then the prediction
#' accuracy is reported for each iteration. 
#' 
#' @param CNOlist Cnolist which contains all the experiments  
#' @param model a model prepared for the training  
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV),
#'  it is not recommended for large datasets.  
#' @param foldid an optional vector of values between `1` and `nfold` identifying what fold each observation is in. If supplied, `nfold` can be missing.  
#' @param type define the way to do the crossvalidation. The default is 
#' `type="datapoint"`, which assigns the data randomly into folds. 
#' The option `type="experiment"` uses whole experiments for crossvalidation 
#' (all data corresponding to a cue combination). The `type=observable` uses the
#' subset of nodes across all experiments for crossvalidation.  
#' @param ... further arguments are passed to parEstimationLBode  
#' @seealso \link{\code{parEstimationLBode}}  

crossvalidateODE = function(CNOlist, model, nfolds=10, foldid=NULL, type='datapoint', parallel=FALSE, ode_parameters=NULL, paramsSSm=NULL, method = "essm", ...){
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  
  crossvalidate.call = match.call(expand.dots = TRUE)
  
  if(!(type%in%c('datapoint','experiment','observable'))){
    
    error("Wrong argument for type. It should either be 'datapoint', 'experiment' or 'observable'")
    
  } else {
    
    if(type=='datapoint'){
      N = prod(dim(CNOlist@signals[[1]]))
    }else if(type=='experiment'){
      N = nrow(CNOlist@signals[[1]])
    }else if(type=="observable"){
      N = ncol(CNOlist@signals[[1]])
    }
    
  }
  
  if (is.null(foldid)){
    
    foldid = sample(rep(seq(nfolds), length = N))
    
  } else {
    
    nfolds = max(foldid)
    
  }
  
  outlist = as.list(seq(nfolds))
  
  if (parallel) {
    require(doParallel)
    outlist = foreach(i = seq(nfolds), .packages = c("CNORode")) %dopar% 
    {
      whichI = foldid == i
      CNOlist.sub = CNOlist
      CNOlist.cv = CNOlist
      if(type=='datapoint'){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][whichI] = NA
          CNOlist.cv@signals[[timeIndex]][!whichI] = NA
        }
      }else if(type=='experiment'){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][whichI,] = NA
          CNOlist.cv@signals[[timeIndex]][!whichI,] = NA
        }
      }else if(type=="observable"){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][,whichI] = NA
          CNOlist.cv@signals[[timeIndex]][,!whichI] = NA
        }
      }
      
      if(is.null(ode_parameters)){
        
        ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                           LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                           default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                           opt_tau = TRUE, random = TRUE)
        
      }
      
      if(is.null(paramsSSm)){
        
        paramsSSm=defaultParametersSSm()
        paramsSSm$local_solver = "DHC"
        paramsSSm$maxtime = 600;
        paramsSSm$maxeval = Inf;
        paramsSSm$atol=1e-6;
        paramsSSm$reltol=1e-6;
        paramsSSm$nan_fac=0;
        paramsSSm$dim_refset=30;
        paramsSSm$n_diverse=1000;
        paramsSSm$maxStepSize=Inf;
        paramsSSm$maxNumSteps=10000;
        transferFun=4;
        paramsSSm$transfer_function = transferFun;
        
        paramsSSm$lambda_tau=0
        paramsSSm$lambda_k=0
        paramsSSm$bootstrap=F
        paramsSSm$SSpenalty_fac=10
        paramsSSm$SScontrolPenalty_fac=1000
        
      }
      
      outlist = list()
      res = parEstimationLBode(cnolist = CNOlist.sub, model = model, ode_parameters = ode_parameters, method = method, paramsSSm = paramsSSm)
      outlist$fit = res
    
      outlist$cvScore = res$ssm_results$fbest
      outlist
    }
    
  }
  else {
    cat("fold \t fitScore \t cvScore \t nTolModels \t bString\n")
    
    for (i in seq(nfolds)) {
      whichI = foldid == i
      CNOlist.sub = CNOlist
      CNOlist.cv = CNOlist
      if(type=='datapoint'){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][whichI] = NA
          CNOlist.cv@signals[[timeIndex]][!whichI] = NA
        }
      }else if(type=='experiment'){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][whichI,] = NA
          CNOlist.cv@signals[[timeIndex]][!whichI,] = NA
        }
      }else if(type=="observable"){
        for(timeIndex in 1:length(CNOlist@timepoints)){
          CNOlist.sub@signals[[timeIndex]][,whichI] = NA
          CNOlist.cv@signals[[timeIndex]][,!whichI] = NA
        }
      }
      
      if(is.null(ode_parameters)){
        
        ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                           LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                           default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                           opt_tau = TRUE, random = TRUE)
        
      }
      
      if(is.null(paramsSSm)){
        
        paramsSSm=defaultParametersSSm()
        paramsSSm$local_solver = "DHC"
        paramsSSm$maxtime = 600;
        paramsSSm$maxeval = Inf;
        paramsSSm$atol=1e-6;
        paramsSSm$reltol=1e-6;
        paramsSSm$nan_fac=1000;
        paramsSSm$dim_refset=30;
        paramsSSm$n_diverse=1000;
        paramsSSm$maxStepSize=Inf;
        paramsSSm$maxNumSteps=10000;
        transferFun=4;
        paramsSSm$transfer_function = transferFun;
        
        paramsSSm$lambda_tau=0
        paramsSSm$lambda_k=0
        paramsSSm$bootstrap=F
        paramsSSm$SSpenalty_fac=10
        paramsSSm$SScontrolPenalty_fac=1000
        
      }
      
      res = parEstimationLBode(cnolist = CNOlist.sub, model = model, ode_parameters = ode_parameters, method = method, paramsSSm = paramsSSm)
      
      outlist[[i]] = list()
      outlist[[i]]$fit = res
      
      outlist[[i]]$cvScore = res$ssm_results$fbest
      
    }
    
  }
  
  return(outlist)
}
