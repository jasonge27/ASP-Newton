lognet_KKT <- function(data, beta, intercept, lambda){
  # Y = 0/1, loss function is 
  # 1/n \sum_{i=1}^n (y_i <x_i, beta> - log(1 + exp(<x_i, beta>))) + \lambda |\beta|_1
  return(compute_lognet_KKT(data$X, c(data$Y), c(data$X%*%beta), c(beta), intercept, lambda))
}

lognet_loss <- function(data, beta, intercept, lambda){
  return(compute_lognet_loss(data$X, c(data$Y), c(data$X%*%beta), c(beta), intercept, lambda))
}

generate_sim_lognet <- function(n, d, c, seed=1024) {
  set.seed(seed)
  cor.X <- c 
  S <- matrix(cor.X,d,d) + (1-cor.X)*diag(d)
  R <- chol(S)

  X <- scale(matrix(rnorm(n*d),n,d)%*%R)*sqrt(n-1)/sqrt(n)
  attributes(X) <- NULL
  X <- matrix(X, n,d)

  s <- 20
  true_beta <- c(runif(s), rep(0, d-s)) 

  # strictly seperable
  Y <- rbinom(n=n, size=1, p = 1/(1+exp(-X%*%true_beta)))

  return(list(X=X, Y=c(Y), true_beta=true_beta))
}

generate_sim<- function(n, d, c, seed=1024) {
  set.seed(seed)
  cor.X <- c 

  S <- matrix(cor.X,d,d) + (1-cor.X)*diag(d)
  R <- chol(S)

  X <- scale(matrix(rnorm(n*d),n,d)%*%R)*sqrt(n-1)/sqrt(n)
  attributes(X) <- NULL
  X <- matrix(X, n,d)

  s <- 20
  true_beta <- c(runif(s), rep(0, d-s)) 
  Y <- X%*%true_beta+rnorm(n)*5
  return(list(X=X, Y=c(Y), true_beta = true_beta))
}

pathfista <- function(data, lambdas, tol=1e-6, max_it=100){
  nlambda <- length(lambdas)

  
  n <- dim(data$X)[1]
  d <- dim(data$X)[2]
  
  out <- list(beta = matrix(rep(0,d*nlambda), ncol=nlambda),
             b0 = rep(0,n),
             lambda = lambdas)
  
  bprev <- matrix(rep(0, d+1), ncol=1)
  for (i in 1:nlambda){
     y <- matrix(2*data$Y - 1, ncol=1)
     x <- cbind(data$X, matrix(rep(1,n), ncol=1))
     fit <- spams.fistaFlat(y, x, bprev, 
                            loss='logistic', regul = 'l1', tol=tol, max_it = max_it,
                            intercept = TRUE, # no not regularize last row of beta
                            lambda1 = lambdas[i])
     out$beta[,i] <- fit[1:d,1]
     out$b0[i] <- fit[d+1,1]
     bprev <- fit
  }
  return(out)
}

timing_lognet <- function(data, nlambda = 100, ratio = 0.01, fista_it = 20, trialN = 10, skip=c()) {
  library(glmnet)
  library(picasso)
  library(gcdnet)
  library(spams)


  fitp<-picasso(data$X, data$Y,family="binomial", lambda.min.ratio=ratio,
                                    standardize=FALSE, verbose=FALSE, prec=2.0*1e-4, nlambda=nlambda)

  regpath <- fitp$lambda
  reg <- fitp$lambda[length(regpath)]

  reg_breakpoint <- seq(1, length(regpath), 10)
  reg_breakpointN <- length(reg_breakpoint)


  cat("ASP-Newton timing:\n")
  rtime <- rep(0, reg_breakpointN) 
  loss <- rep(0, reg_breakpointN)
  for (i in 1:reg_breakpointN){
    t <- system.time(fitp<-picasso(data$X, data$Y,family="binomial", lambda = regpath[1:reg_breakpoint[i]],
                                    standardize=FALSE, verbose=FALSE, prec=2.0*1e-6))
    rtime[i] <- t[1]
    j <- reg_breakpoint[i]  # k lambdas
    loss[i] <- lognet_loss(data, fitp$beta[,j], fitp$intercept[j], reg)
  }
  cat(paste(rtime, collapse=","))
  cat("\n")
  cat(paste(loss, collapse=","))
  cat("\n")


  if (!("ncvreg" %in% skip)){
    cat("ncvreg timing:\n")
  }

   if (!("glmnet" %in% skip)){
    cat("glmnet timing:\n")
    rtime <- rep(0, reg_breakpointN) 
    loss <- rep(0, reg_breakpointN)
    
    for (i in 1:reg_breakpointN){
      t <- system.time(fit<-glmnet(data$X, data$Y, family="binomial", 
                                   lambda = regpath[1:reg_breakpoint[i]],
                                   standardize=FALSE, thresh=2*1e-6))
      rtime[i] <- t[1]
      j <- reg_breakpoint[i]
      loss[i] <- lognet_loss(data, fit$beta[,j], fit$a0[j], reg)

    }
    cat(paste(rtime, collapse=","))
    cat("\n")
    cat(paste(loss, collapse=","))
    cat("\n")

  }


  reg_breakpoint <- seq(1, 21, 2)
  reg_breakpointN <- length(reg_breakpoint)


  if (!("gcdnet" %in% skip)){
    cat("gcdnet timing:\n")
    rtime <- rep(0, reg_breakpointN) 
    loss <- rep(0, reg_breakpointN)
    for (i in 1:reg_breakpointN){
      t <- system.time(fit<-gcdnet(data$X, data$Y, method="logit", 
                                   lambda = regpath[1:reg_breakpoint[i]],
                                   standardize=FALSE, eps=2*1e-6))
      rtime[i] <- t[1]
      j <- reg_breakpoint[i]
      loss[i] <- lognet_loss(data, fit$beta[,j], fit$b0[j], reg)
      }
    cat(paste(rtime, collapse=","))
    cat("\n")
    cat(paste(loss, collapse=","))
    cat("\n")
  }



  if (!("fista" %in% skip)){
    cat("fista timing:\n")
    rtime <- rep(0, reg_breakpointN) 
    loss <- rep(0, reg_breakpointN)
    
    for (i in 1:reg_breakpointN){
      t <- system.time(fit<-pathfista(data, regpath[1:reg_breakpoint[i]], max_it=fista_it))
      rtime[i] <- t[1]
      j <- reg_breakpoint[i]
      loss[i] <- lognet_loss(data, fit$beta[,j], fit$b0[j], reg)
    }
    cat(paste(rtime, collapse=","))
    cat("\n")
    cat(paste(loss, collapse=","))
    cat("\n")
  }

}

test_lognet <- function(data, nlambda = 50, ratio=0.01, fista_it = 20, trialN = 3, skip=c()){
  p = dim(data$X)[2]
  cat("ASP-Newton timing:\n")
  picasso.rtime <- rep(0, trialN) 
  picasso.KKTerr <- rep(0, trialN)
  for (i in 1:trialN){
     t <- system.time(fitp<-picasso(data$X, data$Y,family="binomial", lambda.min.ratio=ratio,
                                    standardize=FALSE, verbose=FALSE, prec=2.0*1e-4, nlambda=nlambda))
     picasso.rtime[i] <- t[1]
     err <- rep(0, nlambda)
     for (j in 1:nlambda){
       err[j] <- lognet_KKT(data, fitp$beta[,j], fitp$intercept[j], fitp$lambda[j])
     }
     picasso.KKTerr[i] <- mean(err)
  }
  cat("mean running time: \n")
  print(mean(picasso.rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(picasso.rtime)))
  cat("mean KKT error: \n")
  print(mean(picasso.KKTerr))
  cat("standard deviation of KKT error: \n")
  print(sqrt(var(picasso.KKTerr)))
  
  if (!("ncvreg" %in% skip)){
    cat("ncverg timing:\n")
    
    tryCatch({
      rtime <- rep(0, trialN) 
      KKTerr <- rep(0, trialN)
    
      for (i in 1:trialN){
        t <- system.time(fit<- ncvreg(data$X, data$Y, family='binomial', 
              penalty='lasso',
              eps=1e-4,
              lambda=fitp$lambda))
        rtime[i] <- t[1]
        err <- rep(0, nlambda)
        for (j in 1:nlambda){
          err[j] <- lognet_KKT(data, fit$beta[2:(p+1),j], fit$beta[1,j], fit$lambda[j])
        }
        KKTerr[i] <- mean(err)
      }

      cat("mean running time: \n")
      print(mean(rtime))

      cat("standard deviation of running time: \n")
      print(sqrt(var(rtime)))

      cat("mean KKT error: \n")
      print(mean(KKTerr))

      cat("standard deviation of KKT error: \n")
      print(sqrt(var(KKTerr))) 
      },

      error=function(e){
        print("ncvreg runs into error")
        print(e)}
    )
   
  }
  
  if (!("glmnet" %in% skip)){
    cat("glmnet timing:\n")
    rtime <- rep(0, trialN) 
    KKTerr <- rep(0, trialN)
    
    for (i in 1:trialN){
      t <- system.time(fit<-glmnet(data$X, data$Y, family="binomial", 
                                   lambda = fitp$lambda,
                                   standardize=FALSE, thresh=2*1e-6))
      rtime[i] <- t[1]
      err <- rep(0, nlambda)
      for (j in 1:nlambda){
        err[j] <- lognet_KKT(data, fit$beta[,j], fit$a0[j], fit$lambda[j])
      }
      KKTerr[i] <- mean(err)
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("standard deviation of KKT error: \n")
    print(sqrt(var(KKTerr)))
  }

  if (!("gcdnet" %in% skip)){
    cat("gcdnet timing:\n")
    rtime <- rep(0, trialN) 
    KKTerr <- rep(0, trialN)
    for (i in 1:trialN){
      t <- system.time(fit<-gcdnet(data$X, data$Y, method="logit", 
                                   lambda = fitp$lambda,
                                   standardize=FALSE, eps=2*1e-6))
      rtime[i] <- t[1]
      err <- rep(0, nlambda)
      for (j in 1:nlambda){
        err[j] <- lognet_KKT(data, fit$beta[,j], fit$b0[j], fit$lambda[j])
      }
      KKTerr[i] <- mean(err)
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("standard deviation of KKT error: \n")
    print(sqrt(var(KKTerr)))
  }
  
  if (!("fista" %in% skip)){
    cat("fista timing:\n")
    rtime <- rep(0, trialN) 
    KKTerr <- rep(0, trialN)
    lambdas <- fitp$lambda
    for (i in 1:trialN){
      t <- system.time(fit<-pathfista(data, lambdas, max_it=fista_it))
      rtime[i] <- t[1]
      err <- rep(0, nlambda)
      for (j in 1:length(lambdas)){
        err[j] <- lognet_KKT(data, fit$beta[,j], fit$b0[j], lambdas[j])
      }
      KKTerr[i] <- mean(err)
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("standard deviation of KKT error: \n")
    print(sqrt(var(KKTerr)))
  }
}