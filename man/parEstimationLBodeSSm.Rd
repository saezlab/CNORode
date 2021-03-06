\name{parEstimationLBodeSSm}

\alias{parEstimationLBodeSSm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Perform parameter estimation using essR.
}

\description{
	This function uses essR to perform parameter estimation. The objective function
	is the same as the one provided by \code{\link{getLBodeContObjFunction}}.
}

\usage{
	parEstimationLBodeSSm(cnolist, model, ode_parameters = NULL, indices = NULL, 
	maxeval = Inf, maxtime = 100, ndiverse = NULL, dim_refset = NULL, local_solver = NULL,
	time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, 
	maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50, nan_fac = 1,
	lambda_tau = 0, lambda_k = 0, bootstrap = FALSE, SSpenalty_fac = 0, 
    SScontrolPenalty_fac = 0, boot_seed = sample(1:10000,1))
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{cnolist}{
		A list containing the experimental design and data.
	}
	\item{model}{
		The logic model to be simulated.
	}
	\item{ode_parameters}{
		A list with the ODEs parameter information. Obtained with \code{\link{createLBodeContPars}}.
	}
	\item{indices}{
		Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
	}
	\item{maxeval}{
		Maximum number of evaluation in the optimization procedure.
	}
	\item{maxtime}{
		Duration of the optimization procedure.
	}
	\item{ndiverse}{
		Number of diverse initial solutions.
	}
	\item{dim_refset}{
		Size of the reference set.
	}
	\item{local_solver}{
		Local solver to be used in SSm.
	}
	\item{time}{
		An integer with the index of the time point to start the simulation. Default is 1.
	}
	\item{verbose}{
		A logical value that triggers a set of comments.
	}
	\item{transfer_function}{
		The type of used transfer. Use 1 for no transfer function, 2 for Hill function and
		3 for normalized Hill function.
	}
	\item{reltol}{
		Relative Tolerance for numerical integration.
	}
  	\item{atol}{
		Absolute tolerance for numerical integration.
	}
	\item{maxStepSize}{
		The maximum step size allowed to ODE solver.
	}
	\item{maxNumSteps}{
		The maximum number of internal steps between two points being sampled before
		the solver fails.
	}
 	\item{maxErrTestsFails}{
		Specifies the maximum number of error test failures permitted in attempting
	 	one step.
	}
  	\item{nan_fac}{
		A penalty for each data point the model is not able to simulate. We recommend higher than
		0 and smaller that 1.
	}
	\item{lambda_tau}{
		penalty parameter for node parameters (tau)
	}
	\item{lambda_k}{
		penalty parameter for edge parameters (k)
	}
	\item{bootstrap}{
		Boolean, default: FALSE. If the residuals should be bootstrapped. 
	}
	\item{SSpenalty_fac}{
		Steady-state penalty: at the end of the simulation the model states should 
		reach steady state. The steady state is measured by the sum of sqares of the 
		state derivatives. 
	}
	\item{SScontrolPenalty_fac}{
		Steady-state penalty for a control experiment, the default is 0.
		The first condition should represent a control condition 
		(no stimulus or inhibition). Then the model
		simulation is penalised if it deviates from the initial conditions. This is 
		to make sure that the predicted dynamics is not due to the initial conditions, 
		but becuase of the stimuli.
	}
	\item{boot_seed}{
		random seed used for the bootsrapping. 
	}
      
}

\details{
	Check \code{\link{CellNOptR}} for details about the cnolist and the model format.
	For more details in the configuration of the ODE solver check the CVODES manual. 
}

\value{
	\item{LB_n}{
		A numeric value to be used as lower bound for all parameters of type n.
	}
	\item{LB_k}{
		A numeric value to be used as lower bound for all parameters of type k.
	}
  	\item{LB_tau}{
		A numeric value to be used as lower bound for all parameters of type tau.
	}
	\item{UB_n}{
		A numeric value to be used as upper bound for all parameters of type n.
	}
	\item{UB_k}{
		A numeric value to be used as upper bound for all parameters of type k.
	}
	\item{UB_tau}{
		A numeric value to be used as upper bound for all parameters of type tau. 
	}	

	\item{default_n}{
		The default parameter to be used for every parameter of type n. 
	}
	\item{default_k}{
		The default parameter to be used for every parameter of type k. 
	}
	\item{default_tau}{
		The default parameter to be used for every parameter of type tau. 
	}
	\item{LB_in}{
		An array with the the same length as ode_parameters$parValues with lower bounds
		for each specific parameter. 
	}
	\item{UB_in}{
		An array with the the same length as ode_parameters$parValues with upper bounds
		for each specific parameter.
	}
	\item{opt_n}{
		Add all parameter n to the index of parameters to be fitted.
	}
	\item{opt_k}{
		Add all parameter k to the index of parameters to be fitted.
	}
	\item{opt_tau}{
		Add all parameter tau to the index of parameters to be fitted.
	}
	\item{random}{
		A logical value that determines that a random solution is for the parameters 
		to be optimized.
	}
	\item{smm_results}{A list containing the information provided by the nonlinear optimization solver.}
}

\author{
	David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{createLBodeContPars}}
}

\examples{
\dontrun{
data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");

ode_parameters=createLBodeContPars(model,random=TRUE);

#Visualize intial simulation
simulatedData=plotLBodeFitness(cnolistCNORodeExample, model,ode_parameters,indices=indices)

ode_parameters=parEstimationLBodeSSm(cnolistCNORodeExample,model,ode_parameters,
indices=indices,maxtime=20,ndiverse=50,dim_refset=6);

#Visualize fitterd solution
simulatedData=plotLBodeFitness(cnolistCNORodeExample, model,indices=indices,ode_parameters=ode_parameters);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

