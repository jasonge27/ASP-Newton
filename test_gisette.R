library(Rcpp)
library(glmnet)
library(picasso)
library(gcdnet)
library(spams)



source("scripts.R")
sourceCpp("utils.cpp")

load("gisette_5k.RData")
n = dim(gisette$X)[1]
gisette$X[which(is.na(gisette$X))] <- 0 # missing values
gisette$X <- gisette$X[ ,find_nonconstant_column(gisette$X)]
gisette$X <- scale(gisette$X)/sqrt(n-1) * sqrt(n)

trialN = 1

test_lognet(gisette)
