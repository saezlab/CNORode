# CNORode2017
This is modified version of [CNORode](https://www.bioconductor.org/packages/release/bioc/html/CNORode.html) which is an add-on to [CellNOptR](https://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html). For more general information about the CellNOpt project visit: http://www.cellnopt.org/

### Basic info on logic based ODE and CNORode
For an introduction on logic based ODE we recomment reading the original publication [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/). Briefly, logic based ODE are ordinary differential equations (ODE) derived from logic rules using continuous update function (*B<sub>i</sub>*), which allows to have a continuous description of the behaveour of the species of interest both in time and in state. Each species *i* is described by an ODE:

*d(x<sub>i</sub>)/dt=&tau;<sub>i</sub>(B<sub>i</sub>(f(x<sub>i1</sub>), f(x<sub>i2</sub>), ..., f(x<sub>iN</sub>))-x<sub>i</sub>)*

where *&tau;<sub>i</sub>* is the life-time of the species and *x<sub>i1</sub>, x<sub>i2</sub>, ... x<sub>iN</sub>* are its *N* regulators. Each regulation is described by a transfer function *f(x<sub>ij</sub>)* which can be, for example, a linear relationsip or a sigmoidal (Hill like) curve.

### New features of CNORode2017
Differences with respect to the original CNORode package are described below.

###### New transfer function
A new transfer function *f(x<sub>ij</sub>)* was introduces to have a more streightforward interpretability of the parameters in terms of functionlity of the edges. For the previously implemented transfer functions we refer to [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/) and [(Terfve et al., BMC Syst Biol, 2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3605281/). The new transfer function is characterized by two parameters (*k<sub>ij</sub>* and *n<sub>ij</sub>*) and mantains the previously introduced sigmoidal shape, but allows to describe the case of "no regulation" when parameter *k<sub>ij</sub>=0*. For fixed *n<sub>ij</sub>*, increasing values of *k<sub>ij</sub>* correspond to increasing strength of the regulation.

The transfer function to be used in the function for parameters estimation, i.e. *parEstimationLBode()*, is passed using the *paramsSSm* argument.

```R
# sets optimisation parameters to default values
paramsSSm=defaultParametersSSm()
# use new transfer functions (see help for other transfer functions)
paramsSSm$transfer_function=4
```


###### L1 regularisation
A penalty term proportional to the L1-norm of the parameters has been added to the objective function in order to induce sparsity in the network

### Reference
The described additions to the CNORode package has been developed for the following publication (available as a preprint in [BioRxiv](http://biorxiv.org/content/early/2016/12/16/094755)):

> F. Eduati, V. Doldàn-Martelli, B. Klinger,  T. Cokelaer, A. Sieber, F. Kogera,  M. Dorel,  M. J Garnett,  N. Blüthgen,  J. Saez-Rodriguez. Dissecting cancer resistance to therapies with cell-type-specific dynamic logic models. bioRxiv 094755 (2016) doi: https://doi.org/10.1101/094755
