# CNORode2017
This is a modified version of [CNORode](https://www.bioconductor.org/packages/release/bioc/html/CNORode.html) which is an add-on to [CellNOptR](https://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html). We refer to the [*vignette*](https://www.bioconductor.org/packages/release/bioc/vignettes/CNORode/inst/doc/CNORode-vignette.pdf) for learning usage of CNORode, and we will describe here below only added features. For more general information about the CellNOpt project visit: http://www.cellnopt.org/.

### How to install

*CNORode2017* requires *CellNOptR* and *MEIGOR* which are available in Bioconductor and can be installed typing:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("CellNOptR")
biocLite("MEIGOR")
```

*CNORode2017* can be installed typing:

```R
library(devtools)
install_github("saezlab/CNORode2017")
```

Note that libraries should be loaded in this exact order:

```R
library(CellNOptR)
library(MEIGOR)
library(CNORode2017)
```

This is because *CNORode* is a dependency of *MEIGOR* and has same function names as *CNORode2017*, so original *CNORode* functions will be prioritized when calling the optimisation function if *MEIGOR* (and thus *CNORode*) is loaded after *CNORode2017*.

### Basic info on logic based ODE and CNORode
For an introduction on logic based ODE we recommend reading the original publication [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/). Briefly, logic based ODE are ordinary differential equations (ODE) derived from logic rules using continuous update function (*B<sub>i</sub>*), which allows to have a continuous description of the behaviour of the species of interest both in time and in state. Each species *x<sub>i</sub>* is described by an ODE:

*d(x<sub>i</sub>)/dt=&tau;<sub>i</sub>(B<sub>i</sub>(f(x<sub>i1</sub>), f(x<sub>i2</sub>), ..., f(x<sub>iN</sub>))-x<sub>i</sub>)*

where *&tau;<sub>i</sub>* is the life-time of the species and *x<sub>i1</sub>, x<sub>i2</sub>, ... x<sub>iN</sub>* are its *N* regulators. Each regulation is described by a transfer function *f(x<sub>ij</sub>)* which can be, for example, a linear relationsip or a sigmoidal (Hill like) curve.

### New features of *CNORode2017*
Differences with respect to the original CNORode package are described below. All new features can be used in the *parEstimationLBode()* function for parameters estimation and can be passed using the *paramsSSm* argument (can be used only when using *essm* method):

```R
# sets optimisation parameters to default values
paramsSSm=defaultParametersSSm()
# passed to parameter estimation function as follows
opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
```
Using default *paramsSSm* should get to the same results as the original *CNORode*.


#### 1. New transfer function
A new transfer function *f(x<sub>ij</sub>)* was introduces to have a more straightforward interpretability of the parameters in terms of functionality of the edges. For the previously implemented transfer functions we refer to [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/) and [(Terfve et al., BMC Syst Biol, 2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3605281/). The new transfer function is characterized by two parameters (*k<sub>ij</sub>* and *n<sub>ij</sub>*) and mantains the previously introduced sigmoidal shape, but allows to describe the case of "no regulation" when parameter *k<sub>ij</sub>=0*. For fixed *n<sub>ij</sub>*, increasing values of *k<sub>ij</sub>* correspond to increasing strength of the regulation.

Previously implemented transfer functions are passed as numberes: 1 (linear), 2 (Hill) and 3 (Normalised Hill). The new transfer function has been assigned value 4 and can be passed as *paramsSSm* argument as follows:

```R
# use new transfer function (see help for other transfer functions)
paramsSSm$transfer_function=4
```

Same when running the simulation:

```R
simulatedData=plotLBodeFitness(cnolist, model, transfer_function=4, ode_parameters=opt_pars)
```

#### 2. L1 regularisation
A penalty term proportional to the L1-norm of the parameters has been added to the objective function in order to induce sparsity in the network. The balance between prioritising good fit or sparse model is regulated by a tunable term *&lambda;* which can be empirically selected testing the effect of different values on model fit and model sparsity.
A regularisation factor *&lambda;* can be assigned separately to parameters *&tau;<sub>i</sub>* and *k<sub>ij</sub>* as follows:


```R
# for parameters tau
paramsSSm$lambda_tau=0.01
# for parameters k
paramsSSm$lambda_k=0.001
```

#### 3. Steady state penalty
An additional penalty term can be considered in the objective function to penalise parameter sets for which the simulation does not reach steady state within the time range of the measured experimental data. This is especially recommended to match biological assumptions when only few time points are available in the experimental data and have been chosen considering that the biological system under investigation have reached a semi steady state. A different (typically higher) penalty can be given to the control (unperturbed) condition in order to favour that the system remains at basal condition when no perturbations are applied.

```R
# steady state penalty
paramsSSm$SSpenalty_fac=10
# additional penalty for control condition (assumed to be first row of MIDAS)
paramsSSm$SScontrolPenalty_fac=1000
```

#### 4. Bootstrap

This feature allows to perform the optimisation with random resampling (with replacement) of the experimental data. With default option for *boot_seed*, a random seed for the random sampling is selected each time the optimisation is run, but the seed can also be provided as user input.

```R
# select bootstrapping option
paramsSSm$bootstrap=TRUE
# default random assignment of the seed for random sampling
paramsSSm$boot_seed=sample(1:10000,1)
```

When using the bootsrap, the user should repeat the optimisation multiple times (recommended at least 100) to derive a bootstrapped distribution of the estimated parameters. Here below a simple code that can be used for small networks, while for medium or large networks parallelisation (preferably on cluster) is strongly recommended.

```R
# 100 bootstrap repetitions
opt_pars<-lapply(seq(1:100), function(x){
  paramsSSm$boot_seed=x
  parEstimationLBode(cnolist, model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
})
```

### Summary of new/modified arguments for *paramsSSm*

| Name | Default | Suggested | Description |
| :------------ | :---------- | :---------- | :---------- |
| ***transfer_function*** | 3 | 4 |  Defines which transfer function (*f(x<sub>ij</sub>)*) to use. Default (*transfer_function=3*) is Normalised Hill. New transfer function is *transfer_function=4*. |
| ***lambda_tau*** | 0 | 0.0001-0.1 |  Regularisation factor *&lambda;<sub>&tau;</sub>* to be used for parameters *&tau;<sub>i</sub>*. |
| ***lambda_k*** | 0 | 0.0001-0.1 | Regularisation factor *&lambda;<sub>k</sub>* to be used for parameters *k<sub>ij</sub>*. |
| ***SSpenalty_fac*** | 0 | 10 | Penalty factor for not reaching steady state within time range of experimental data |
| ***SScontrolPenalty_fac*** | 0 | 1000 | Penalty factor for not reaching steady state for control (unperturbed) condition. Control condition is assumed to be in the first row of the MIDAS or CNOlist. |
| ***bootstrap*** | FALSE |  |  If *bootstrap=TRUE* experimental data used in the optimisation are randomly sampled with replacement. If this option is selected the optimisation should be repeated multiple times (&ge; 100) to obtain a distribution of the estimated parameters. |
| ***boot_seed*** | sample(1:10000,1) |  |  Seed to be used for random sampling if *bootstrap=TRUE*. |

### Example
A working example is available [here](https://github.com/saezlab/CPT_QSPtutorial/blob/master/CellNOptR_optimisation.R) and data and network (PKN) used are available in the same repository (https://github.com/saezlab/CPT_QSPtutorial).

### Reference
The described additions to the CNORode package have been developed for the following publication (available as a preprint in [BioRxiv](http://biorxiv.org/content/early/2016/12/16/094755)):

> F. Eduati, V. Doldàn-Martelli, B. Klinger,  T. Cokelaer, A. Sieber, F. Kogera,  M. Dorel,  M. J Garnett,  N. Blüthgen,  J. Saez-Rodriguez. Dissecting cancer resistance to therapies with cell-type-specific dynamic logic models. bioRxiv 094755 (2016) doi: https://doi.org/10.1101/094755
