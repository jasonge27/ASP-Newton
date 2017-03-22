# loading and installing required packages
library(Rcpp)
library(glmnet)
library(picasso)
library(gcdnet)
library(spams)

source("scripts.R")
sourceCpp("utils.cpp")

set.seed(111)
sim_wc <- generate_sim(n=2000, d=10000, 0.3) 
sim_ic <- generate_sim(n=2000, d=10000, 3.0)

load("madelon.RData")
madelon$X[which(is.na(madelon$X))] <- 0 # missing values
madelon$X <- madelon$X[ ,find_nonconstant_column(madelon$X)]
madelon$X <- scale(madelon$X)

load("gisette.RData")
gisette$X[which(is.na(gisette$X))] <- 0 # missing values
gisette$X <- gisette$X[ ,find_nonconstant_column(gisette$X)]
gisette$X <- scale(gisette$X)

# Reproducing Table 1
# Comparisons of Running Time Across Different Methods on Logistic Regression Tasks.
test_lognet(madelon)
test_lognet(gisette)
test_lognet(sim_wc)
test_lognet(sim_ic)

# Reproducing Figure 4
# Comparisions of Convergence Rate Measured in CPU Time.
timing_lognet(sim_ic)
