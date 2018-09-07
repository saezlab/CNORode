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
#' createLBodeInitialConditions
#' 
#' create a matrix of initial conditions for all state variables.
#' 
#' creating an initial condition matrix gives more flexibility to the user to define
#' indicidual values. In earlier versions all non-measured nodes were set to 0.5
#'  or 0.1, depending on the implementation. 
#'  
#' @param model a model list created by readSIF and optionally be preprocessing
#' @param data a cnolist created from the MIDAS file
#' @param initialValue in case of a scalar, all states in all experimental conditions
#' are set to this value. In case of a matrix with nExps rows and nSpecies columns then the initial
#' values are set accorsingly. See useMeasurements argument for further details.
#' @param useMeasurements logical, if TRUE (default) then the initial values of the measured nodes are 
#' set based on the measured data.
#' @return numerical matrix with nExps rows and nSpecies columns
#' @export 
#' @author A.Gabor
#' @examples
#' \dontrun{
#' initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.2)
#' S1 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue1)
#' plotLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters, initialValueMatrix =  initValue1)
#' }
#' 
createLBodeInitialConditions = function(model, data, initialValue=0.0, useMeasurements = TRUE){
	
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
