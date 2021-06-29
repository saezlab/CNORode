\name{crossvalidateODE}
\alias{crossvalidateODE}
\title{
k-fold crossvalidation for ODE models.
}
\description{
Cross-validation analysis for the logic-ode case.
}
\usage{
crossvalidateODE = function(CNOlist = cnolist, model = model, nfolds=10, foldid=NULL, type="datapoint", parallel = FALSE, 
                            ode_parameters = NULL, paramsSSm = NULL, method = "essm")
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at time 1).
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that is not a requirement of this function.
}
  \item{nfolds}{
number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets.
}
  \item{foldid}{
an optional vector of values between `1` and `nfold` identifying what fold each observation is in. If supplied, `nfold` can be missing.
}
  \item{type}{
define the way to do the crossvalidation. The default is type="datapoint"`, which assigns the data randomly into folds. The option `type="experiment"` uses whole experiments for crossvalidation (all data corresponding to a cue combination). The `type=observable` uses the subset of nodes across all experiments for crossvalidation.  
}
  \item{parallel}{
verbose parameter, indicating wheter to parallelize the cross-validation procedure or not (default set to FALSE).
}
  \item{ode_parameters}{
a list with the ODEs parameter information. Obtained with createLBodeContPars.
}
  \item{paramsSSm}{
a list of SSm parameters. default is the list returned bydefaultParametersSSm.
}


}
\details{
Does a k-fold cross-validation for logic-ode models. In k-iterations a  fraction of the data is eliminated from the CNOlist. The model is trained on the  remaining data and then the model predicts the held-out data. Then the prediction accuracy is reported for each iteration.
}
\value{
This function returns a list with elements:
\item{fit}{optimisation result}
\item{cvScore}{cross-validation scores}
}

\author{
A. Gabor, E. Gjerga
}
\examples{
\dontrun{
library(CellNOptR)
library(doParallel)

data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")
pknmodel = ToyModel
cnolist = CNOlist(CNOlistToy)

# original and preprocessed network 
plotModel(pknmodel,cnolist)
model = preprocessing(data = cnolist,model = pknmodel,compression = T,expansion = T)
plotModel(model,cnolist)
plotCNOlist(CNOlist = cnolist)

# set initial parameters 
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## Parameter Optimization
# essm
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
paramsSSm$SSpenalty_fac=0
paramsSSm$SScontrolPenalty_fac=0

# run the optimisation algorithm
opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, 
                            paramsSSm=paramsSSm)
plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars, 
                  transfer_function = 4)

# 10-fold crossvalidation using T1 data
# We use only T1 data for crossvalidation, because data in the T0 matrix is not independent.
# All rows of data in T0 describes the basal condition.

# Crossvalidation produce some text in the command window:  
registerDoParallel(cores=3)
R=crossvalidateODE(CNOlist = cnolist, model = model, type="datapoint", nfolds=3,
parallel = TRUE, ode_parameters = ode_parameters, paramsSSm = paramsSSm)
}
}

\references{
A. MacNamara, C. Terfve, D. Henriques, B.P. Bernabe, J. Saez-Rodriguez. State-time spectrum of signal transduction logic models, Phys Biol., 9(4):045003, 2012. Toy example can be found here: https://saezlab.github.io/CellNOptR/
}
