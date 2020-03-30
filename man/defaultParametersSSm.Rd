\name{defaultParametersSSm}

\alias{defaultParametersSSm}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	Create default options to perform parameter estimation with scatter search meta-heuristic.
}

\description{
	This function returns a list with several arguments for performing parameter estimation with scatter search meta-heuristic 
	algorithm from the package essR.
}

\usage{
	defaultParametersSSm()
}

\value{
	\item{maxeval}{Inf}
	\item{maxtime}{100}
	\item{ndiverse}{NULL}
	\item{dim_refset}{NULL}
	\item{local_solver}{NULL}
	\item{verbose}{0}
	\item{transfer_function}{3}
	\item{reltol}{1e-04}
	\item{atol}{0.001}
	\item{maxStepSize}{Inf}
	\item{maxNumSteps}{1e+05}
	\item{maxErrTestsFails}{50}
	\item{nan_fac}{1}	
	\item{lambda_tau}{0}	
	\item{lambda_k}{0}	
	\item{bootstrap}{0}	
	\item{SSpenalty_fac}{0}	
	\item{SScontrolPenalty_fac}{0}	
	\item{boot_seed}{sample(1:10000,1)}	
}

\author{
	David Henriques, Thomas Cokelaer, Federica Eduati
}

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{parEstimationLBode}}
	\code{\link{parEstimationLBodeSSm}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{default}
\keyword{SSm}
\keyword{essR}
