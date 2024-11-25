## let a(.) be a given function
Compute_Prec <- function(a, c, k, region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
  a_func <- function(x) {-1/(a*(x+c))}
  a_matrix <- a_func(x)
  B1a <- as(apply(B1matrix, 2, function(x) x*a_matrix), "dgCMatrix")

  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")

  G <- t(B1a) %*% Numerical_I %*% B1a
  C <- t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}

Compute_Design <- function(x, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  Bmatrix
}

Compute_Prec_rev <- function(a, c, k, region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
  a_func <- function(x) {1/(a*(x+c))}
  a_matrix <- a_func(x)
  B1a <- as(apply(B1matrix, 2, function(x) x*a_matrix), "dgCMatrix")

  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")

  G <- t(B1a) %*% Numerical_I %*% B1a
  C <- t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}

