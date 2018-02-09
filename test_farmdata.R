args = commandArgs(trailingOnly=TRUE)

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
load("datasets/farmads/farmads.RData")
}


nlambda = as.integer(args[1])
ratio = as.double(args[2]) 

prec = list(picasso=1.0*1e-4, ncvreg=1e-2, glmnet=1*1e-5 )
test_lognet(farmdata, prec, nlambda=nlambda, ratio=ratio, skip=c('ncvreg', 'gcdnet', 'fista'))