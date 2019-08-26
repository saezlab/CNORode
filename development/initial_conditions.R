load_all()

# set up initial conditions

cnodata = CNOlist("data/ToyModel_PB/ToyModelPB.csv")

pkn_dataGen= readSIF("data/ToyModel_PB/ToyModelPB_TRUE.sif")
plotModel(pkn_dataGen,cnodata)


model = preprocessing(cnodata,pkn_dataGen,expansion = F, compression = T)
plotModel(model,cnodata)

ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
								   LB_tau = 0, UB_n = 5, UB_k = 5, UB_tau = 10, default_n = 3,
								   default_k = 0.5, default_tau = 0.01, opt_n = TRUE, opt_k = TRUE,
								   opt_tau = TRUE, random = TRUE)
ode_parameters

#' createLBodeInitialConditions
#' 
#' create a matrix of initial conditions for all state variables.
#' 
#' creating an initial condition matrix gives more flexibility to the user to define
#' individual values. In earlier versions all non-measured nodes were set to 0.5
#'  or 0.1, depending on the implementation. 
#'  @param model a model list created by readSIF and optionally be preprocessing
#'  @param data a cnolist created from the MIDAS file
#'  @param initialValue in case of a scalar, all states in all experimental conditions
#'  are set to this value. In case of a matrix with nExps rows and nSpecies columns then the initial
#'  values are set accorsingly. See useMeasurements argument for further details.
#'  @param useMeasurements logical, if TRUE (default) then the initial values of the measured nodes are 
#'  set based on the measured data.
#'  @return numerical matrix with nExps rows and nSpecies columns
#'  @export 
createLBodeInitialConditions = function(model, data, initialValue=0.5, useMeasurements = TRUE){
	
	# checking inputs
	if (class(data)!="CNOlist"){
		data = CNOlist(data)
	}
	checkSignals(CNOlist = data, model = model)
	
	
	# initialise
	initialValueMatrix = matrix(0, ncol = length(model$namesSpecies),nrow = nrow(data@cues))
	colnames(initialValueMatrix) = model$namesSpecies
	
	# assignment
	if(is.null(dim(initialValue))){
		# scalar case
		initialValueMatrix[] = initialValue
	}else{
		# matrix case
		if(all(dim(initialValueMatrix) == dim(initialValue))){
			
			if(is.null(colnames(initialValue))){
				# [] will assure to keep the colnames
				initialValueMatrix[] = initialValue
			}else{
				# the user provided a names matrix, make sure that the columns are in correct order:
				colOrder = match(colnames(initialValue),colnames(initialValueMatrix))	
				initialValueMatrix[,colOrder] = initialValue
			}
			
			# [] will assure to keep the colnames
			initialValueMatrix[] = initialValue
		}else{
			stop("the provided initialValue matrix has wrong dimensions.")
		}
	}
	
	if(useMeasurements){
		initialValueMatrix[,colnames(data@signals[[1]])] = data@signals[[1]]
	}
	
	return(initialValueMatrix)
}

initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.2)
initValue2 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.2,useMeasurements = F)
initValue3 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.0,useMeasurements = T)

S1 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue1)
S2 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue2)
S3 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue3)

S1[[1]]
S2[[1]]
S3[[1]]

load_all()
plotLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters)
plotLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters, initialValueMatrix =  initValue1)
plotLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue2)
plotLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue3)
