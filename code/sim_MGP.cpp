#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  
  int d = P.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param
  DATA_SCALAR(sig); // standard deviation of the noise
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)


  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = B * U + X * beta
  int Wdim = W.size();
  int betadim = Wdim - d;
  vector<Type> U(d);
  vector<Type> beta(betadim);
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i+d);
  PARAMETER(theta); // theta = -2log(sigma1)

  // Transformations
  vector<Type> eta = X * beta + B * U;
  Type sigma = exp(-0.5*theta);
  REPORT(sigma);

  // Log likelihood (disabled for simulation purpose)
  Type ll = 0;
  // ll = sum(dnorm(y, eta, sig, TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta) * UPU; // U part

  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  

  // Log determinant
  Type logdet = d * theta + logPdet;
  lpW += 0.5 * logdet; // P part

  REPORT(logdet);
  REPORT(lpW);

  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}
