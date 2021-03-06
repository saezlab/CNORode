\name{getLBodeModelSim}

\alias{getLBodeModelSim}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	Simulate the logic-based ODE model
}

\description{
	This function simulates a logic-based ODE model and return a list with one matrix for each time point.
	The input species in the model are filled with NA values. If the simulation of a particular set of 
	initial conditions fails the solver will fill the experience row with NA values.
}

\usage{
	getLBodeModelSim(cnolist, model, ode_parameters = NULL, indices = NULL, timeSignals=NULL,
	time = 1,verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf,
	maxNumSteps = 1e+05, maxErrTestsFails = 50)
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
	\item{timeSignals}{
		An array containing a different timeSignals. If you use this argument, it will also
		modify the dimensions from valueSignals.
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
		The maximum number of internal steps between two points being sampled before
		the solver fails.
	}
  	\item{maxNumSteps}{
		The maximum number of internal steps between two points being sampled before
		the solver fails.
	}
	\item{maxErrTestsFails}{
		Specifies the maximum number of error test failures permitted in attempting
		one step.
	}
}

\details{
	Check \code{\link{CellNOptR}} for details about the cnolist and the model format.
	For more details in the configuration of the ODE solver check the CVODES manual. 
}

\value{
	Returns a list with simulated data with similar structure to cnolist$valueSignals.
	Contains one matrix for each time-point. Each matrix contains one row per experiment
	and one columns per model species.
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
library(CNORode)
data('ToyCNOlist',package='CNORode');
data('ToyModel',package='CNORode');
data('ToyIndices',package='CNORode');
modelSimulation=getLBodeModelSim(cnolistCNORodeExample, model,indices=indices);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

