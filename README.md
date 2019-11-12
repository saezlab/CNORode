# CNORode: a logic based ordinary differential equation add-on for CellNOptR

This version of CNORode is a continuation of the [CNORode2017](https://github.com/saezlab/CNORode2017) package. CNORode2017 was lamor equivalent to CNORode with some added features, such as the sparsity enforcing
regularisation and parameter uncertainty analysis based on bootstrapping. 
Maintaining two packages are time consuming, therefore we decided to merge the packages and continue with CNORode2017.

The older version of CNORode can be downloaded from [release_v1.23](https://github.com/saezlab/CNORode/releases/tag/v1.23.0)

## Where to start
We refer to the [*vignette*](https://www.bioconductor.org/packages/release/bioc/vignettes/CNORode/inst/doc/CNORode-vignette.pdf)
for learning usage of the classic CNORode. 

The added features of the CNORode2017 are summarised on the main site of [CNORode2017](https://github.com/saezlab/CNORode2017). 

For more general information about the CellNOpt project visit: http://www.cellnopt.org/.


### How to install

*CNORode* requires *CellNOptR* and *MEIGOR* which are available in Bioconductor and can be installed typing:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("CellNOptR")
biocLite("MEIGOR")
```

Alternatively, you can obtain the most recent version of *CellNOptR* from our GitHub repository:

```R
library(devtools)
install_github("saezlab/CellNOptR")
```

*CNORode* can be installed from this repository by typing:

```R
library(devtools)
install_github("saezlab/CNORode")
```

Install the developement version of CNORode:
```R
library(devtools)
install_github("saezlab/CNORode",ref="crossval", build_vignettes = TRUE)
```


### Basic info on logic based ODE and CNORode
For an introduction on logic based ODE we recommend reading the original publication [(Wittmann et al., BMC Syst Biol., 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2764636/). Briefly, logic based ODE are ordinary differential equations (ODE) derived from logic rules using continuous update function (*B<sub>i</sub>*), which allows to have a continuous description of the behaviour of the species of interest both in time and in state. Each species *x<sub>i</sub>* is described by an ODE:

*d(x<sub>i</sub>)/dt=&tau;<sub>i</sub>(B<sub>i</sub>(f(x<sub>i1</sub>), f(x<sub>i2</sub>), ..., f(x<sub>iN</sub>))-x<sub>i</sub>)*

where *&tau;<sub>i</sub>* is the life-time of the species and *x<sub>i1</sub>, x<sub>i2</sub>, ... x<sub>iN</sub>* are its *N* regulators. Each regulation is described by a transfer function *f(x<sub>ij</sub>)* which can be, for example, a linear relationsip or a sigmoidal (Hill like) curve.

### New features of *CNORode2017*
Please visit [CNORode2017](https://github.com/saezlab/CNORode2017). 

