
#' @title runCNORode
#'
#' @description A one-line wrapper of the CNORode pipeline
#'
#' @param model A filename of prior knowledge network (PKN) in the SIF format
#' @param data A measurement filename in the MIDAS format
#' @param results_folder  results folder for the analysis. 
#' @param compression compress the prior knowledge network (TRUE), see \code{\link{preprocessing}}
#' @param cutNONC cut non-observable non-controllable node from PKN (TRUE), see \code{\link{preprocessing}}
#' @param expansion expand OR gates in the PKN (FALSE), see \code{\link{preprocessing}}
#' @param LB_n lower bound on parameter n, see \code{\link{createLBodeContPars}}
#' @param LB_k lower bound on parameter k, see \code{\link{createLBodeContPars}}
#' @param LB_tau lower bound on parameter tau, see \code{\link{createLBodeContPars}}
#' @param UB_n upper bound on parameter n, see \code{\link{createLBodeContPars}}
#' @param UB_k upper bound on parameter k, see \code{\link{createLBodeContPars}}
#' @param UB_tau upper bound on parameter tau, see \code{\link{createLBodeContPars}}
#' @param default_n default value of parameter n, see \code{\link{createLBodeContPars}}
#' @param default_k default value of parameter k, see \code{\link{createLBodeContPars}}
#' @param default_tau default value of parameter tau, see \code{\link{createLBodeContPars}}
#' @param opt_n should parameter n be optimised, see \code{\link{createLBodeContPars}}
#' @param opt_k should parameter k be optimised, see \code{\link{createLBodeContPars}}
#' @param opt_tau should parameter tau be optimised, see \code{\link{createLBodeContPars}}
#' @param random initial parameter vector generation (TRUE: random, FALSE: half of the LB-UB)
#' @param maxeval maximum number of funciton evaluations in the optimisation, see \code{\link{parEstimationLBodeSSm}}
#' @param maxtime maximum CPU time (in seconds) spent on optimisation before calling final refinement, see \code{\link{parEstimationLBodeSSm}}
#' @param transfer_function trandfer function types represented by the edges, see \code{\link{parEstimationLBodeSSm}}
#' @param nan_fac penalty for NA simulations, see \code{\link{parEstimationLBodeSSm}}
#' @param lambda_tau regularisation penalty for tau parameters, see \code{\link{parEstimationLBodeSSm}}
#' @param lambda_k regularisation penalty for k parameters for optimisation, see \code{\link{parEstimationLBodeSSm}}
#' 
#' @examples
#' \dontrun{
#' model = system.file("extdata", "ToyModelMMB_FeedbackAnd.sif",package="CNORode")
#' data = system.file("extdata", "ToyModelMMB_FeedbackAnd.csv", package="CNORode")
#' res = runCNORode(model,data,results_folder = "./results")
#' }
#' @import CellNOptR
#' @export

runCNORode <- function(model,
					  data,
					  compression=TRUE,
					  results_folder = "CNORode_results",
					  cutNONC=TRUE,
					  expansion=FALSE,     
					  LB_n = 1,
					  LB_k = 0.1,
					  LB_tau = 0.01,
					  UB_n = 5,
					  UB_k = 0.9,
					  UB_tau = 10,
					  default_n = 3,
					  default_k = 0.5,
					  default_tau = 1,
					  opt_n = TRUE,
					  opt_k = TRUE,
					  opt_tau = TRUE,
					  random = TRUE,
					  maxeval = 1e5,
					  maxtime = 60,
					  transfer_function = 3,
					  nan_fac = 1,
					  lambda_tau = 0,
					  lambda_k = 0
) {
	
	
	# Manual assignment of network (PKN) and experimental data (MIDAS) files
	pknmodel <- readSIF(model) # build a prior-knowledge network from SIF file
	cnolist <- CNOlist(data) # import experimental data from MIDAS file
	
	if(!dir.exists(results_folder)) dir.create(results_folder)
	
	pdf(file.path(results_folder,"1_network_original.pdf"))
	plotModel(pknmodel,CNOlist = cnolist)
	dev.off()
	
	
	# Now use CNORode specific pre-processing (expansion=FALSE by default and report)
	model <- preprocessing(data=cnolist,model=pknmodel, expansion=expansion, # model in CNORode shouldn't be expanded by pre-processing step
									   compression=compression, cutNONC=cutNONC,
									   verbose=FALSE) 
	
	pdf(file.path(results_folder,"2_network_preprocessed.pdf"))
	plotModel(model,CNOlist = cnolist)
	dev.off()
	
	# generate parameterised model from model object
	init_ode_parameters = createLBodeContPars(model, LB_n = 1, LB_k = 0.1,
										   LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3,
										   default_k = 0.5, default_tau = 1, opt_n = TRUE, opt_k = TRUE,
										   opt_tau = TRUE, random = FALSE)

	
	# fit the model to data
	fit_result = parEstimationLBodeSSm(cnolist = cnolist,
									model = model,
									ode_parameters = init_ode_parameters,
									maxeval = maxeval,
									maxtime = maxtime,
									local_solver = "DHC",
									transfer_function = transfer_function,
									nan_fac = nan_fac,
									lambda_k = lambda_k,
									lambda_tau = lambda_tau
									)
	
	# report results
	pdf(file.path(results_folder,"3_model_fittness.pdf"))
	simulatedData=plotLBodeFitness(cnolist = cnolist, 
								   model = model,
								   ode_parameters = fit_result,
								   transfer_function = transfer_function)
	
	dev.off()
	
	
	# collect results 
	output = list()
	output$model = model
	output$cnolist = cnolist
	output$init_ode_parameters = init_ode_parameters
	output$fit_results = fit_result
	output$model_siulation = simulatedData
	
	
	save(output, file = file.path(results_folder,"CNORode_PostHocResults.Rdata"))
	
	return(output)
}
