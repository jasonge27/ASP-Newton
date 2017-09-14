library(flare)
library(picasso)
library(Rcpp)

source("scripts.R")
sourceCpp("utils.cpp")

set.seed(111)
sim_wc <- generate_sim(n=1000, d=100, 0.3) 


ratio = 0.01
nlambda = 100
fitp<-picasso(sim_wc$X, sim_wc$Y,family="sqrtlasso", lambda.min.ratio=ratio,
                                    standardize=TRUE, verbose=TRUE, prec=1.0*1e-4, nlambda=nlambda)
print(fitp$df)

fitflare <- slim(sim_wc$X, sim_wc$Y, lambda = fitp$lambda, method='lq', q=2, prec= 1e-6)
print(fitflare$df)

browser()
nlambda = length(fitp$lambda)
n = 1000
for (i in 1:nlambda){
    loss.picasso <- sqrt(sum((sim_wc$Y - sim_wc$X %*% fitp$beta[, i] - fitp$intercept[i])^2)/n) + fitp$lambda[i]*sum(abs(fitp$beta[, i]))
    loss.flare <- sqrt(sum((sim_wc$Y - sim_wc$X %*% fitflare$beta[, i] - fitflare$intercept[i])^2)/n ) + fitp$lambda[i]*sum(abs(fitflare$beta[, i]))
    cat("picasso loss: ")
    cat(loss.picasso)
    cat("   flare loss: ")
    cat(loss.flare)
    cat('\n')
}