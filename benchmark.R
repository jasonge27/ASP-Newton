# loading and installing required packages
library(Rcpp)
library(glmnet)
library(picasso)
library(gcdnet)
library(spams)
library(ncvreg)

source("scripts.R")
sourceCpp("utils.cpp")

if (TRUE){

set.seed(111)
sim_ic <- generate_sim_lognet(n=1000, d=5000, 0.8)

sim_wc <- generate_sim_lognet(n=5000, d=1000, 0.2) 

load("datasets/madelon/madelon.RData")
madelon$X[which(is.na(madelon$X))] <- 0 # missing values
madelon$X <- madelon$X[ ,find_nonconstant_column(madelon$X)]
madelon$X <- scale(madelon$X)

load("datasets/gisette/gisette.RData")
gisette$X[which(is.na(gisette$X))] <- 0 # missing values
gisette$X <- gisette$X[ ,find_nonconstant_column(gisette$X)]
gisette$X <- scale(gisette$X)

}

load("datasets/farmads/farmads.RData")
farmdata$Y <- (farmdata$Y == -1)
farmdata$X <- as.matrix(farmdata$X)
farmdata$X <- farmdata$X[ ,find_nonconstant_column(farmdata$X)]
farmdata$X <- scale(farmdata$X)

if (TRUE){
# Reproducing Table 1
# Comparisons of Running Time Across Different Methods on Logistic Regression Tasks.
cat("===========madelon==========\n")
prec = list(picasso=5.0*1e-6, ncvreg=1e-2, glmnet=5*1e-5 )
test_lognet(madelon, prec)

cat("===========gisette==========\n")
prec = list(picasso=1.0*1e-4, ncvreg=1e-2, glmnet=5*1e-5 )
test_lognet(gisette, prec)

cat("===========farmdata==========\n")
prec = list(picasso=2.0*1e-6, ncvreg=1e-2, glmnet=5*1e-5 )
test_lognet(farmdata, prec, skip=c('gcdnet', 'fista'))
cat("===========simwc==========\n")
prec = list(picasso=2.0*1e-6, ncvreg=1e-2, glmnet=5*1e-5 )
test_lognet(sim_wc, prec)
cat("==========simic===========\n")
prec = list(picasso=2.0*1e-6, ncvreg=1e-2, glmnet=5*1e-5 )
test_lognet(sim_ic, prec)
}

# Reproducing Figure 4
# Comparisions of Convergence Rate Measured in CPU Time.
#timing_lognet(farmdata)
