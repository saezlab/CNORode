
#' title runCNORprob
#'
#' description A one-linear wrapper of the CNORprob pipeline
#'
#' param model A filename of prior knowledge network (PKN) in the SIF format
#' param data A measurement filename in the MIDAS format
#' param compression compre 
#' param ProbCutNONC A variable to define whether to cut non-observable non-controllable node from PKN
#' param ProbExpandOR A variable to define whether to expand OR gate in the CNOR format
#' param ORlist A list of interaction with OR gate (leave NULL is none)
#' param ProbHardConstraint Apply hard constraints on probability parameters (sum of all activatory edges = 1)
#' param ProbForce A variable to define whether to force all single activatory edge to a node to always be 1
#' param optRound_optim Number of rounds for optimisation
#' param L1Reg Weight for L1-regularisation
#' param MaxTime Time for each round of optimisation (seconds)
#' param HLbound A cut-off for high and low weights
#' param SSthresh A cut-off for states' difference to define whether steady-state is reached
#' param printCost A variable to select whether to print or not to print intermediate fitting cost (0,1)
#' param PlotIterations Number of rounds for optimisation to generate plots
#' param SaveOptResults A variable to define whether to generate reports from optimisation (TRUE,FALSE)
#' param rho Penalty weighting scaler / default = 1
#' param outer.iter Maximum major iterations / default = 400
#' param inner.iter Maximum minor iterations / default = 800
#' param delta Relative step size eval. / default = 1e-7
#' param tol Relative tol. for optim. / default = 1e-8
#' param iter Print objfunc every itereration / default = 1
#' param Analysis_edgeKO Perform edge KO analysis (T/F)
#' param Analysis_nodeKO Perform node KO analysis (T/F)
#' param Analysis_LPSA Perform local parameter sensitivity analysis (T/F)
#' param Analysis_BS Perform bootstrapping analysis (T/F)
#' param optRound_analysis Number of rounds for optmisation in each post-optimisation analysis
#' param LPSA_Increments Number of increments in LPSA analysis
#' param BS_Type Type of Bootstrapping (1=resample with replacement from residual; 2=resampling from mean & variant)
#' param BS_Round Number of rounds for bootstrapping
#' 
#' 
#' export

runCNORode = function(model,
					  data,
					  compression=FALSE,
					  cutNONC=TRUE,
					  ProbExpandOR=FALSE,     
					  ORlist=NULL,
					  ProbHardConstraint=F,      
					  ProbForce=F,
					  optRound_optim=3,        
					  L1Reg=1e-4,     
					  MaxTime=180,      
					  HLbound=0.5,      
					  SSthresh=2e-16,   
					  printCost=0,        
					  PlotIterations=1,        
					  SaveOptResults=TRUE,     
					  rho=1,        
					  outer.iter=100,       
					  inner.iter=100,       
					  delta=1e-7,     
					  tol=1e-8,     
					  trace=1,        
					  Analysis_edgeKO=T,
					  Analysis_nodeKO=T,
					  Analysis_LPSA=T,
					  Analysis_BS=T,
					  optRound_analysis=1,       
					  LPSA_Increments=2,        
					  BS_Type=1,        
					  BS_Round=5       
) {
	
	
	# Manual assignment of network (PKN) and experimental data (MIDAS) files
	pknmodel <- readSIF(model) # build a prior-knowledge network from SIF file
	CNOlist <- CNOlist(data) # import experimental data from MIDAS file
	# Now use CNORprob specific pre-processing (expansion=FALSE by default and report)
	ModDatppProb <- preprocessing_Prob(CNOlist, pknmodel, expansion=FALSE, # model in CNORprob shouldn't be expanded by pre-processing step
									   compression=ProbCompression, cutNONC=ProbCutNONC,
									   verbose=FALSE) 
	optmodel <- ModDatppProb$cutModel
	optCNOlist <- ModDatppProb$data
	CNORprob_inputs <- NULL
	CNORprob_inputs$optmodel=optmodel
	CNORprob_inputs$optCNOlist=optCNOlist
	CNORprob_inputs$ProbExpandOR=ProbExpandOR
	CNORprob_inputs$ProbHardConstraint=ProbHardConstraint
	CNORprob_inputs$ProbForce=ProbForce
	CNORprob_inputs$ORlist=ORlist
	
	
	
	estim_Result   <<- list() # Initialise global variable of results
	rsolnp_options  <- list(rho=rho,outer.iter=outer.iter,inner.iter=inner.iter,delta=delta,tol=tol,trace=trace) # Collapase optimisation options
	
	# Generate optimisation object
	estim <- CNORprob_buildModel(optCNOlist,optmodel,expandOR=CNORprob_inputs$ProbExpandOR,ORlist=CNORprob_inputs$ORlist,HardConstraint=CNORprob_inputs$ProbHardConstraint,Force=CNORprob_inputs$ProbForce,L1Reg=L1Reg,HLbound=HLbound,SSthresh=SSthresh,PlotIterations=PlotIterations,rsolnp_options=rsolnp_options)
	
	# Run Optimisation
	estim$maxtime <- MaxTime; estim$printCost <- printCost
	estim$optimOptions <- c(rho,outer.iter,inner.iter,delta,tol,trace)
	res <- CNORprob_optimise(estim,optRound_optim,SaveOptResults)
	
	# Plot results
	CNORprob_plotFit(optmodel,optCNOlist,estim,res,show=TRUE, plotPDF=TRUE, tag=NULL, plotParams=list(cex=0.8, cmap_scale=1, ymin=0))
	MappedProb <- CNORprob_mapModel(optmodel,optCNOlist,estim,res)
	pdf("Results/Figures/Optimised_CNORprob_Model.pdf")
	plotModel(MappedProb$model,MappedProb$CNOlist,round(MappedProb$bString,digits=2))
	dev.off()
	
	# Post-optimisation analyses
	estim$ProbCompression <- CNORprob_inputs$ProbCompression
	estim$ProbCutNONC <- CNORprob_inputs$ProbCutNONC
	estim$ProbExpandOR <- CNORprob_inputs$ProbExpandOR
	estim$optRound_analysis <- optRound_analysis
	estim_original <- estim
	estim_based <- estim_original; if (Analysis_edgeKO) { estim_Result  <- CNORprob_edgeKO(optmodel,optCNOlist,estim_based,res) }
	estim_based <- estim_original; if (Analysis_nodeKO) { estim_Result  <- CNORprob_nodeKO(optmodel,optCNOlist,estim_based,res) }
	estim_based <- estim_original; if (Analysis_LPSA) { estim_Result  <- CNORprob_LPSA(estim_based,res,HLbound,LPSA_Increments,Force=F) }
	estim_based <- estim_original; if (Analysis_BS) { estim_Result  <- CNORprob_BS(optmodel,optCNOlist,estim_based,res,BS_Type,BS_Round) }
	
	save(estim_Result,file="Results/CNORprob_PostHocResults.Rdata")
	
	return(estim_Result)
	
}


##

#

#	%\VignetteIndexEntry{Main vignette:Playing with networks using CNORode}
#%\VignetteKeywords{Training Signalling Pathway Maps to
	# %    Biochemical Data with Logic-Based Ordinary Differential Equations
	# %}
# 
# %\VignettePackage{CNORode}
# 
# \documentclass{article}
# \usepackage{Sweave,fullpage}
# \usepackage{url, color}
# 
# %\usepackage{cmap}
# 
# \usepackage{authblk}
# \usepackage[T1]{fontenc}
# \usepackage[utf8]{inputenc}
# 
# \usepackage{hyperref}
# 
# \hypersetup{
# 	colorlinks, linkcolor=blue
# }
# 
# \RequirePackage[T1]{fontenc}
# \RequirePackage{graphicx,ae,fancyvrb,color}
# \IfFileExists{upquote.sty}{\RequirePackage{upquote}}{}
# \setkeys{Gin}{width=0.8\textwidth}
# 

# \definecolor{gris90}{gray}{0.90}
# \definecolor{gris10}{gray}{0.1}
# \definecolor{green}{rgb}{0.6, 0.9,0.6}
# 
# 
# \DefineVerbatimEnvironment{Sinput}{Verbatim}{%
# 	fontshape=sl,
# 	frame=single,
# 	xleftmargin=2em,
# 	fillcolor=\color{gris90},
# 	%    numbers=left % prevent copy/paste entire code
# }
# \DefineVerbatimEnvironment{Soutput}{Verbatim}{
# 	fontshape=sl,
# 	frame=single,
# 	xleftmargin=2em,
# 	fillcolor=\color{green}
# }


# \DefineVerbatimEnvironment{shell}{Verbatim}{formatcom=\color{blue}}
# 
# \title{Training Signalling Pathway Maps to
# 	Biochemical Data with Logic-Based Ordinary Differential Equations
# }
# \author[1,2]{David Henriques}
# \author[1]{Thomas Cokelaer \thanks{cokelaerebi.ac.uk}}
# 
# \affil[1]{European Bioinformatics Institute, Saez-Rodriguez group, Cambridge, United Kingdom}
# \affil[2]{Instituto de Investigaciones Marinas-CSIC, Vigo, Spain}
# 
# \begin{document}
# \maketitle

# \tableofcontents

# 
# Mathematical models are used to understand protein signalling networks so as to
# provide an integrative view of pharmacological and toxicological processes at
# molecular level. \emph{CellNOptR}~\cite{CellNOptR, CellNOptR_paper} is an existing
# package (see
# 		 \url{http://bioconductor.org/packages/release/bioc/html/CellNOptR.html}) that
# provides functionalities to combine
# prior knowledge network (about protein signalling networks) and perturbation
# data to infer functional characteristics (of the signalling network).
# While \emph{CellNOptR} has demonstrated its ability to infer new functional
# characteristics, it is based on a boolean formalism where protein species are
# characterised as being fully active or inactive. In contrast, logic-based ordinary
# differential equations allow a quantitive description of a given Boolean model.
# 
# The method used here was first published by Wittmann et al. \cite{wittman} by
# the name of \emph{odefy}.
# For a detailed description of the methodology the user is adressed to \cite{wittman} and for
# a published application example to \cite{wittmanJS}.
# 
# This package implements the Odefy method and focus mainly extending the \emph{CellNOptR}
# capabilities in order to simulate and calibrate logic-based ordinary differential
# equation model. We provide direct and easy to use interface to optimization methods
# available in R such as \emph{eSSR} \cite{eSSR1} (enhanced Scatter Search Metaheuristic for R) and an R genetic
# algorithm implementation by the name of \emph{genalg} in order to perform parameter estimation.
# Additionally we were specially careful in tackling the main computanional bottlenecks
# by implementing CNORode simulation engine in the C language using the
# \emph{CVODES} library \cite{cvodes}.

# This brief tutorial shows how to use \emph{CNORode} using as a starting point a Boolean model
# and a dataset consisting in a time-series of several proteins.
# 
# \section{Installation}
# 
# \emph{CNORode} depends on \emph{CellNOptR} and \emph{genalg}, which are 2 bioconductor
# packages. Therefore, in order to install \emph{CNORode}, open a R
# session and type:
# 
# 	<<installCNOR, eval=FALSE, pgf=TRUE, echo=TRUE>>=
# 	source("http://bioconductor.org/biocLite.R")
# biocLite("CNORode")

# 
# 	It may take a few minutes to install all dependencies if you start from scratch
# (i.e, none of the R packages are installed on your system). Note also that under
# Linux system, some of these packages necessitate the R-devel package to be
# installed (e.g., under Fedora type \emph{sudo yum install R-devel}).
# 
# Additionally, for parameter estimation we recommend the use of \emph{eSSR}. This algorithm
# algorithm is part of the MEIGOR toolbox. Although this package is
# not available on BioConductor, it can be downloaded from
# \url{http://www.ebi.ac.uk/saezrodriguez/cno/downloads.html}. Once downloaded, install
# MEIGOR by typing:
# 
# 
# 	<<installMEIGOR, eval=FALSE>>=
# 	install.packages("MEIGOR_0.99.3_svn2444.tar.gz",type="source")
# 
# 
# 	Note "\~/" should point to the directory where the file is saved.
# 
# Finally, once \emph{CNORode} is installed you can load it by typing:
# 
# 	<<installCNORode2, eval=TRUE>>=
# 	library(CNORode)
# 
# 
# 	\label{sec:quickstart}
# \section{Quick Start}
# 
# In this section, we provide a quick example on how to use \emph{CNORode} to
# find the set of continuous parameters which minimize the squared difference
# between a model simulation and the experimental data.
# 
# Since here we will not be modifying the model structure as opposed to \emph{CellNOptR}
# we will use a model that already contains AND type gates. Such model can be for instance
# the result of calibrating a \textit{prior knowledge network}(PKN) with \emph{CellNOptR}.
# Please note that a PKN can also be used as Boolean model which will contain only OR type gates.
# 
# Detailed information about the model used here (ToyModelMMB\_FeedbackAnd) and additional models can be found at:
# 
# 	\url{http://www.ebi.ac.uk/~cokelaer/cellnopt/data/}
# 
# \begin{center}
# \begin{figure}[ht]
# \includegraphics[height=7cm, width=7cm]{ToyModelMMB_FeedbackAnd}
# \includegraphics[height=7cm, width=7cm]{data_ToyModelMMB_FeddbackAnd}
# \caption{The used model(left panel). A plot from the data, resulting from
# 	the \emph{plotCNOlist} function (right panel).
# 	\label{fig:}}
# \end{figure}
# \end{center}
# 
# The example used here is shipped with the CNORode. In order to load the data and model you
# should type the following commands:
# 
# 	% show data and model loading
# <<quickstart, eval=TRUE, results=hide>>=
# 	library(CNORode)
# model=readSIF(system.file("/doc/ToyModelMMB_FeedbackAnd.sif",
# 						  package="CNORode"));
# cnolist = CNOlist(system.file("/doc/ToyModelMMB_FeedbackAnd.csv",
# 							  package="CNORode"))
# 
# 
# 	The structure from the CNOlist and the Model object is exactly the same as used in the
# \emph{CellNOptR} and therefore for a detailed explanation about these structure we direct
# the reader to the \emph{CellNOptR} manual.
# 
# In order to simulate the model and perform parameter estimation we first need to create
# a list with the ODE parameters associated with each dynamic state as described in \cite{wittman}.
# Each dynamic state will have a \(\tau\) parameter, as many \(n\) and \(k\) parameters as inputs.
# Although the default is to use the normalized Hill function it also possible to use the
# standard Hill or even not to use any transfer function at all. To illustrate the shape of the equations
# associated to each dynamic state and the meaning of each parameter, let us show the differential
# of \textit{Mek}:
# 
# 
# 	\[ \dot{Mek}=\Bigg[  \Bigg( 1 - \frac{Akt^{n_1}/({k_1}^{n_1} + Akt^{n_1})}{1/({k_1}^{n_1} + 1)}\Bigg)
# 					   \cdot \Bigg(\frac{Raf^{n_2}/({k_2}^{n_2} + Raf^{n_2})}{1/({k_2}^{n_2} + 1)}\Bigg) - Mek \Bigg] \cdotp \tau_{Mek}
# 	  \]
# 
# To create a list of ODE parameters we will typically use the \emph{createLBodeContPars} function:
# 
# 	<<>>=
# 	ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0.1,
# 									   LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3,
# 									   default_k = 0.5, default_tau = 1, opt_n = TRUE, opt_k = TRUE,
# 									   opt_tau = TRUE, random = FALSE)
# 
# 
# 	This function creates a general structure where the ODE parameters are ordered according
# to the model. Some tweaks have been added in order to ease tasks we have found
# to be common, nevertheless you can edit several attributes manually. If you print the
# \emph{ode\_parameters}  list you will see the following attributes.
# 
# \clearpage
# <<>>=
# 	print(ode_parameters)
# 
# 	\clearpage
# 
# Typically before running an optimization run you will want to choose which type of parameters
# you want to optimize. The field \emph{index\_opt\_pars} defines which parameters are meant to be optimized.
# In the \emph{createLBodeContPars}, if you choose \(opt\_tau \) as \(TRUE\)  all \( \tau \) parameters will be
# added to the index\_opt\_pars array, the same idea is valid for \( n \) and \( k \) parameters.
# 
# It is also possible to choose default values for lower and upper bounds for the parameters of a given type,
# e.g. \( \tau \) ( \emph{LB\_tau} and \emph{UB\_tau}), as well as a default initial value for such parameters.
# 
# Once we have the ODE parameters structure we are ready to run a simulation or optimization process.
# To run a simulation we can use the \emph{getLBodeModel} or \emph{getLBodeDataSim}, depending on
# if we want to simulate only the signals present in the CNOlist object or all the species in the model. Additionally
# \emph{plotLBodeDataSim} or \emph{plotLBodeModelSim} will also return the values of a model simulation while
# plotting the same values. In figure \ref{fig:plotModelSimFig}, we use plotLBodeModelSim to plot all the experiments
# sampled in 5 different instants (\emph{timeSignals}).
# 
# \begin{figure}[ht]
# \begin{center}
# <<label=plotModelSim,include=TRUE,fig=TRUE>>=
# 	modelSim=plotLBodeModelSim(cnolist, model, ode_parameters,
# 							   timeSignals=seq(0,2,0.5));
# 
# 	\caption{A model simulation plotted with \emph{plotLBodeModelSim} ..
# 		\label{fig:plotModelSimFig}}
# \end{center}
# \end{figure}
# 
# \clearpage
# 
# As previously mentioned, we provide two optimization algorithms that allow parameter estimation
# Both of these algorithms have specific parameters that can be tunned on each specific problem
# (please check CNORode manual for detailed information). For instance, in order to run the
# genetic algorithm for 10 iterations and a population of size of 10, we can use the following code:
# 
# 	<<eval=TRUE, results=hide>>=
# 	initial_pars=createLBodeContPars(model, LB_n = 1, LB_k = 0.1,
# 									 LB_tau = 0.01, UB_n = 5, UB_k = 0.9, UB_tau = 10, random = TRUE)
# #Visualize initial solution
# simulatedData=plotLBodeFitness(cnolist, model,initial_pars)
# paramsGA = defaultParametersGA()
# paramsGA$maxStepSize = 1
# paramsGA$popSize = 50
# paramsGA$iters = 300
# paramsGA$transfer_function = 2
# 
# opt_pars=parEstimationLBode(cnolist,model,ode_parameters=initial_pars,
# 							paramsGA=paramsGA)
# #Visualize fitted solution
# simulatedData=plotLBodeFitness(cnolist, model,ode_parameters=opt_pars)
# 
# 
# 	\begin{center}
# \begin{figure}[ht]
# 
# <<label=plotInit,include=TRUE,fig=TRUE>>=
# 	simulatedData=plotLBodeFitness(cnolist, model,initial_pars)
# 
# 	\caption{The initial solution before optimization. Each row corresponds to an experiment with a particular combination of stimuli and inhibitors.
# 		The columns correspond to the measured values (triangles) and the simulated values (dashed blue lines) from a given signal. The background color
# 		gives an indication of squared difference where red means high error and white low error.}
# \label{fig:plotInitFit}
# 
# \end{figure}
# \end{center}
# 
# 
# \begin{figure}[ht]
# \begin{center}
# <<label=plotFinalFit_fit,include=TRUE,fig=TRUE>>=
# 	simulatedData=plotLBodeFitness(cnolist, model,ode_parameters=opt_pars)
# 
# 	\end{center}
# \caption{A solution obtained by optimization with a genetic algorithm.Each row corresponds to an experiment with a particular combination of stimuli and inhibitors.
# 	The columns correspond to the measured values (triangles) and the simulated values (dashed blue lines) from a given signal. The background color
# 	gives an indication of squared difference where red means high error and white low error.}
# \label{fig:plotFinalFit}
# \end{figure}
# 
# 
# \clearpage
# 
# In addition to eSSR and genalg its is fairly easy to use any other continuous optimization algorithm.
# In the following example we show how to generate and use an the objective function in order to use
# it with a variant of eSSR(part of MEIGOR package) that uses multiple cpus:
# 
# 	<<eval=FALSE>>=
# 	library(MEIGOR)
# f_hepato<-getLBodeContObjFunction(cnolist, model, initial_pars, indices=NULL,
# 								  time = 1, verbose = 0, transfer_function = 2, reltol = 1e-05, atol = 1e-03,
# 								  maxStepSize = Inf, maxNumSteps = 1e4, maxErrTestsFails = 50, nan_fac = 1)
# n_pars=length(initial_pars$LB);
# 
# problem<-list(f=f_hepato, x_L=initial_pars$LB[initial_pars$index_opt_pars],
# 			  x_U=initial_pars$UB[initial_pars$index_opt_pars]);
# 
# #Source a function containing the options used in the CeSSR publication
# source(system.file("benchmarks","get_paper_settings.R",package="MEIGOR"))
# #Set max time as 20 seconds per iteration
# opts<-get_paper_settings(20);
# Results<-CeSSR(problem,opts,Inf,Inf,3,TRUE,global_save_list=c('cnolist','model',
# 															  'initial_pars'))
# \end{document}
