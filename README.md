# ASP-Newton

### Pre-requisites

Installing these following packages on R

```R
install.packages(c("devtools", "glmnet", "gcdnet", "Rcpp"))

library(devtools)
devtools::install_github("jasonge27/picasso")
```

You need to have a c++ compiler on your system to install and use Rcpp.

Please also install the SPAM package from http://spams-devel.gforge.inria.fr/downloads.html follow the instructions there

### Reproducing Simulations Results

```
source("benchmark.R")
```

