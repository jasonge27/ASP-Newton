#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double compute_lognet_loss(NumericMatrix X, NumericVector y, 
                           NumericVector Xbeta, NumericVector beta, 
                           double a0, double lambda) {
  int n = X.nrow();
  int d = X.ncol();

  double L = 0.0;
  for (int i = 0; i < n; i++){
    L +=  (-y[i] * Xbeta[i]) + log(1 + exp(Xbeta[i]));
  }
  L = L/n;

  for (int i = 0; i < d; i++)
    L += lambda * fabs(beta[i]);

  return(L);
}

// [[Rcpp::export]]
double compute_lognet_KKT(NumericMatrix X, NumericVector y, 
                           NumericVector Xbeta, NumericVector beta, 
                           double a0, double lambda) {
  int n = X.nrow();
  int d = X.ncol();
  
  NumericVector grad(d);
  for (int i = 0; i < d; i++)
    grad[i] = 0;
  
  double tmp = 0;
  for (int i = 0; i < n; i++){
    tmp = -y[i] + 1/(1+ exp(-Xbeta[i]-a0));
    for (int j = 0; j < d; j++)
      grad[j] += tmp* X(i,j)/n;
  }
  
  double grad_max = 0;
  for (int i = 0; i < d; i++){
    if (fabs(beta[i]) < 1e-8 ){
      if (fabs(grad[i])<lambda){
        grad[i] = 0;
      } else{
        if (grad[i] >0 )
          grad[i] = grad[i] - lambda;
        else
          grad[i] = grad[i] + lambda;
      }
    } else {
      if (beta[i] >= 1e-8)
        grad[i] = grad[i] + lambda;
      else 
        grad[i] = grad[i] - lambda;
    }

    if (fabs(grad[i]) > grad_max)
      grad_max = fabs(grad[i]);
  }
  return(grad_max);
}


// [[Rcpp::export]]
NumericVector find_nonconstant_column(NumericMatrix m){
  int n = m.nrow();
  int d = m.ncol();

  std::vector<int> idx;
  idx.clear();
  for (int i = 0; i < d; i++){
    bool constant_column = true;
    for (int j = 1; j < n; j++)
      if (fabs(m(j, i)-m(0, i))>1e-8) {
        constant_column = false;
        break;
      }
    if (!constant_column)
      idx.push_back(i);
  }

  NumericVector idxR(idx.size());
  int i = 0;
  for (std::vector<int>::iterator it = idx.begin(); 
        it != idx.end(); it++, i++)
    idxR[i] = (*it+1);

  return(idxR);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R 
// timesTwo(42)
// */
