# CNORode2017
This is modified version of [CNORode](https://www.bioconductor.org/packages/release/bioc/html/CNORode.html) which is an add-on to [CellNOptR](https://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html). We refer to the [*vignette*](https://www.bioconductor.org/packages/release/bioc/vignettes/CNORode/inst/doc/CNORode-vignette.pdf) for learning usage of CNORode, and we will describe here below only added features. For more general information about the CellNOpt project visit: http://www.cellnopt.org/.

### How to install

*CNORode2017* requires *CellNOptR* and *MEIGOR* which are available in Bioconductior and can be installed typing:

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

This is because CNORode is a dependency of MEIGOR and has same function names as CNORode2017 so it will be prioritizes when calling the optimisation function if loaded after CNORode2017.


### New features of *CNORode2017*
Differences with respect to the original CNORode package are described below. All new features can be used in the *parEstimationLBode()* function for parameters estimation passing them using the *paramsSSm* argument (can be used only when using *essm* method):

```R
# sets optimisation parameters to default values
paramsSSm=defaultParametersSSm()
# passed to parameter estimation function as follows
opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
```


### Basic info on logic based ODE and CNORode
For an introduction on logic based ODE we recomment reading the original publication [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/). Briefly, logic based ODE are ordinary differential equations (ODE) derived from logic rules using continuous update function (*B<sub>i</sub>*), which allows to have a continuous description of the behaveour of the species of interest both in time and in state. Each species *i* is described by an ODE:

*d(x<sub>i</sub>)/dt=&tau;<sub>i</sub>(B<sub>i</sub>(f(x<sub>i1</sub>), f(x<sub>i2</sub>), ..., f(x<sub>iN</sub>))-x<sub>i</sub>)*

where *&tau;<sub>i</sub>* is the life-time of the species and *x<sub>i1</sub>, x<sub>i2</sub>, ... x<sub>iN</sub>* are its *N* regulators. Each regulation is described by a transfer function *f(x<sub>ij</sub>)* which can be, for example, a linear relationsip or a sigmoidal (Hill like) curve.


##### 1. New transfer function
A new transfer function *f(x<sub>ij</sub>)* was introduces to have a more streightforward interpretability of the parameters in terms of functionlity of the edges. For the previously implemented transfer functions we refer to [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/) and [(Terfve et al., BMC Syst Biol, 2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3605281/). The new transfer function is characterized by two parameters (*k<sub>ij</sub>* and *n<sub>ij</sub>*) and mantains the previously introduced sigmoidal shape, but allows to describe the case of "no regulation" when parameter *k<sub>ij</sub>=0*. For fixed *n<sub>ij</sub>*, increasing values of *k<sub>ij</sub>* correspond to increasing strength of the regulation.

The transfer function to be used in the *parEstimationLBode()* function for parameters estimation, is passed using the *paramsSSm* argument.

```R
# use new transfer functions (see help for other transfer functions)
paramsSSm$transfer_function=4
```


##### 2. L1 regularisation
A penalty term proportional to the L1-norm of the parameters has been added to the objective function in order to induce sparsity in the network. 

##### 3. Steady state penalty

##### 4. Bootstrap


### Summary of new/modified argumentrs for *paramsSSm*

| Name | Default | Suggested | Description |
| :------------ | :---------- | :---------- | :---------- |
| ***transfer_function*** | 3 | 4 |  Defines which transfer function (*f(x<sub>ij</sub>)*) to use. Default (*transfer_function=3*) is Normalised Hill. New transfer function is *transfer_function=4*. |
| ***lambda_tau*** | 0 | 0.0001-0.1 |  Regularisation factor *&lambda;<sub>&tau;</sub>* to be used for parameters *&tau;<sub>i</sub>*. |
| ***lambda_k*** | 0 | 0.0001-0.1 | Regularisation factor *&lambda;<sub>k</sub>* to be used for parameters *k<sub>ij</sub>*. |
| ***SSpenalty_fac*** | 0 | 10 | Penalty factor for nor reaching steady state withing time range of experimental data |
| ***SScontrolPenalty_fac*** | 0 | 1000 | Penalty factor for nor reaching steady state for control (unperturbed) condition. Control condition is assumed to be in the first row of the MIDAS or CNOlist. |
| ***bootstrap*** | FALSE |  |  If *bootrstrap=TRUE* experimental data used in the optimisation are randomly sampled with replacement. If this option is selected the optimisation should be repeated multiple times (&ge; 100) to obtain a distribution of the estimated parameters. |
| ***boot_seed*** | sample(1:10000,1) |  |  Seed to be used for random sampling if *bootstrap=TRUE*. |

### Example
A working example is available [here](https://github.com/saezlab/CPT_QSPtutorial/blob/master/CellNOptR_optimisation.R) and data and network (PKN) used are avilable in the same reposotory (https://github.com/saezlab/CPT_QSPtutorial).

### Reference
The described additions to the CNORode package has been developed for the following publication (available as a preprint in [BioRxiv](http://biorxiv.org/content/early/2016/12/16/094755)):

> F. Eduati, V. Doldàn-Martelli, B. Klinger,  T. Cokelaer, A. Sieber, F. Kogera,  M. Dorel,  M. J Garnett,  N. Blüthgen,  J. Saez-Rodriguez. Dissecting cancer resistance to therapies with cell-type-specific dynamic logic models. bioRxiv 094755 (2016) doi: https://doi.org/10.1101/094755
