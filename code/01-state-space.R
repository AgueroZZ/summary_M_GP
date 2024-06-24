# case1 : alpha != 1 or 0. s > t = x:
mspline_cov_c1 <- function(s,t,c = 1, alpha = 2){
  (((3*alpha^4 + 2*alpha^3)*c^(1/alpha + 2)*s + (3*alpha^4 + 2*alpha^3)*c^(1/alpha + 3))*(c + s)^(1/alpha) - ((alpha^4 + alpha^3)*c^3 + 3*(alpha^4 + alpha^3)*c^2*t + 3*(alpha^4 + alpha^3)*c*t^2 + (alpha^4 + alpha^3)*t^3)*(c + t)^(2/alpha) + (((3*alpha^4 + 5*alpha^3 + 2*alpha^2)*c + (3*alpha^4 + 5*alpha^3 + 2*alpha^2)*s)*(c + s)^(1/alpha)*t^2 - ((3*alpha^4 + 2*alpha^3)*c^3 + (3*alpha^4 + 2*alpha^3)*c^2*s)*(c + s)^(1/alpha) + (3*alpha^4 + 2*alpha^3)*c^(1/alpha + 3) + (((3*alpha^3 + 2*alpha^2)*c^2 + (3*alpha^3 + 2*alpha^2)*c*s)*(c + s)^(1/alpha) + (3*alpha^4 + 2*alpha^3)*c^(1/alpha + 2))*t)*(c + t)^(1/alpha) - (2*alpha^4 + alpha^3)*c^(2/alpha + 3))*exp(1)^(-(log(c + s) - log(c))/alpha - (log(c + t) - log(c))/alpha - 2*log(c)/alpha)/((6*alpha^2 + 7*alpha + 2)*(alpha^2 + 2*alpha + 1))
}
mspline_deriv_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -(((3*alpha^3 + 2*alpha^2)*c^(1/alpha + 2)*s + (3*alpha^3 + 2*alpha^2)*c^(1/alpha + 3))*(c + s)^(1/alpha) - ((3*alpha^3 + 4*alpha^2 + alpha)*c^3 + 3*(3*alpha^3 + 4*alpha^2 + alpha)*c^2*t + 3*(3*alpha^3 + 4*alpha^2 + alpha)*c*t^2 + (3*alpha^3 + 4*alpha^2 + alpha)*t^3)*(c + t)^(2/alpha) - (2*((3*alpha^4 + 5*alpha^3 + 2*alpha^2)*c + (3*alpha^4 + 5*alpha^3 + 2*alpha^2)*s)*(c + s)^(1/alpha)*t^2 + ((3*alpha^3 + 2*alpha^2)*c^3 + (3*alpha^3 + 2*alpha^2)*c^2*s)*(c + s)^(1/alpha) - (3*alpha^3 + 2*alpha^2)*c^(1/alpha + 3) + (((6*alpha^4 + 13*alpha^3 + 6*alpha^2)*c^2 + (6*alpha^4 + 13*alpha^3 + 6*alpha^2)*c*s)*(c + s)^(1/alpha) - (3*alpha^3 + 2*alpha^2)*c^(1/alpha + 2))*t)*(c + t)^(1/alpha) + (2*alpha^2 + alpha)*c^(2/alpha + 3))/((((alpha^2 + 2*alpha + 1)*c + (alpha^2 + 2*alpha + 1)*s)*(c + s)^(1/alpha)*t + ((alpha^2 + 2*alpha + 1)*c^2 + (alpha^2 + 2*alpha + 1)*c*s)*(c + s)^(1/alpha))*(6*alpha^2 + 7*alpha + 2)*(c + t)^(1/alpha))
}
mspline_cross_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -(((3*alpha^3 + 2*alpha^2)*c^(1/alpha + 2)*s + (3*alpha^3 + 2*alpha^2)*c^(1/alpha + 3))*(c + s)^(1/alpha) + ((3*alpha^4 + 4*alpha^3 + alpha^2)*c^3 + 3*(3*alpha^4 + 4*alpha^3 + alpha^2)*c^2*t + 3*(3*alpha^4 + 4*alpha^3 + alpha^2)*c*t^2 + (3*alpha^4 + 4*alpha^3 + alpha^2)*t^3)*(c + t)^(2/alpha) - (2*((3*alpha^4 + 5*alpha^3 + 2*alpha^2)*c + (3*alpha^4 + 5*alpha^3 + 2*alpha^2)*s)*(c + s)^(1/alpha)*t^2 + ((3*alpha^3 + 2*alpha^2)*c^3 + (3*alpha^3 + 2*alpha^2)*c^2*s)*(c + s)^(1/alpha) + (3*alpha^4 + 2*alpha^3)*c^(1/alpha + 3) + (((6*alpha^4 + 13*alpha^3 + 6*alpha^2)*c^2 + (6*alpha^4 + 13*alpha^3 + 6*alpha^2)*c*s)*(c + s)^(1/alpha) + (3*alpha^4 + 2*alpha^3)*c^(1/alpha + 2))*t)*(c + t)^(1/alpha) - (2*alpha^3 + alpha^2)*c^(2/alpha + 3))*exp(1)^(-(log(c + s) - log(c))/alpha - log(c + t)/alpha)/(((alpha^2 + 2*alpha + 1)*c^(1/alpha)*t + (alpha^2 + 2*alpha + 1)*c^(1/alpha + 1))*(6*alpha^2 + 7*alpha + 2))
}
SS_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -(((6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c^2 + 2*(6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c*s + (6*alpha^4 + 7*alpha^3 + 2*alpha^2)*s^2)*(c + s)^(2/alpha)*t + ((4*alpha^4 + 3*alpha^3)*c^3 + 2*(3*alpha^4 + alpha^3 - alpha^2)*c^2*s - (5*alpha^3 + 4*alpha^2)*c*s^2 - 2*(alpha^4 + 2*alpha^3 + alpha^2)*s^3)*(c + s)^(2/alpha) + ((2*alpha^4 + alpha^3)*c^3 + 3*(2*alpha^4 + alpha^3)*c^2*t + 3*(2*alpha^4 + alpha^3)*c*t^2 + (2*alpha^4 + alpha^3)*t^3)*(c + t)^(2/alpha) - 2*(((3*alpha^4 + 2*alpha^3)*c + (3*alpha^4 + 2*alpha^3)*s)*(c + s)^(1/alpha)*t^2 + 2*((3*alpha^4 + 2*alpha^3)*c^2 + (3*alpha^4 + 2*alpha^3)*c*s)*(c + s)^(1/alpha)*t + ((3*alpha^4 + 2*alpha^3)*c^3 + (3*alpha^4 + 2*alpha^3)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))*exp(1)^(-2*(log(c + s) - log(c))/alpha - 2*log(c)/alpha)/((6*alpha^2 + 7*alpha + 2)*(alpha^2 + 2*alpha + 1))
}


SS_deriv_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -(((6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c^2 + 2*(6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c*s + (6*alpha^4 + 7*alpha^3 + 2*alpha^2)*s^2)*(c + s)^(2/alpha)*t - ((6*alpha^3 + 6*alpha^2 + alpha)*c^3 + (6*alpha^4 + 25*alpha^3 + 20*alpha^2 + 3*alpha)*c^2*s + (12*alpha^4 + 32*alpha^3 + 22*alpha^2 + 3*alpha)*c*s^2 + (6*alpha^4 + 13*alpha^3 + 8*alpha^2 + alpha)*s^3)*(c + s)^(2/alpha) + ((2*alpha^2 + alpha)*c^3 + 3*(2*alpha^2 + alpha)*c^2*t + 3*(2*alpha^2 + alpha)*c*t^2 + (2*alpha^2 + alpha)*t^3)*(c + t)^(2/alpha) + 2*(((3*alpha^3 + 2*alpha^2)*c + (3*alpha^3 + 2*alpha^2)*s)*(c + s)^(1/alpha)*t^2 + 2*((3*alpha^3 + 2*alpha^2)*c^2 + (3*alpha^3 + 2*alpha^2)*c*s)*(c + s)^(1/alpha)*t + ((3*alpha^3 + 2*alpha^2)*c^3 + (3*alpha^3 + 2*alpha^2)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))/(((alpha^2 + 2*alpha + 1)*c^2 + 2*(alpha^2 + 2*alpha + 1)*c*s + (alpha^2 + 2*alpha + 1)*s^2)*(6*alpha^2 + 7*alpha + 2)*(c + s)^(2/alpha))
}
SS_cross_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -(((6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c^2 + 2*(6*alpha^4 + 7*alpha^3 + 2*alpha^2)*c*s + (6*alpha^4 + 7*alpha^3 + 2*alpha^2)*s^2)*(c + s)^(2/alpha)*t + ((3*alpha^4 + alpha^3 - alpha^2)*c^3 + (3*alpha^4 - 4*alpha^3 - 5*alpha^2)*c^2*s - (3*alpha^4 + 11*alpha^3 + 7*alpha^2)*c*s^2 - 3*(alpha^4 + 2*alpha^3 + alpha^2)*s^3)*(c + s)^(2/alpha) - ((2*alpha^3 + alpha^2)*c^3 + 3*(2*alpha^3 + alpha^2)*c^2*t + 3*(2*alpha^3 + alpha^2)*c*t^2 + (2*alpha^3 + alpha^2)*t^3)*(c + t)^(2/alpha) - (((3*alpha^4 - alpha^3 - 2*alpha^2)*c + (3*alpha^4 - alpha^3 - 2*alpha^2)*s)*(c + s)^(1/alpha)*t^2 + 2*((3*alpha^4 - alpha^3 - 2*alpha^2)*c^2 + (3*alpha^4 - alpha^3 - 2*alpha^2)*c*s)*(c + s)^(1/alpha)*t + ((3*alpha^4 - alpha^3 - 2*alpha^2)*c^3 + (3*alpha^4 - alpha^3 - 2*alpha^2)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))*exp(1)^(-(log(c + s) - log(c))/alpha - log(c + s)/alpha)/(((alpha^2 + 2*alpha + 1)*c*c^(1/alpha) + (alpha^2 + 2*alpha + 1)*c^(1/alpha)*s)*(6*alpha^2 + 7*alpha + 2))
}
SS_cov_mat_c1 <- function(s,t,c = 1, alpha = 2){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[1,2] <- SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[2,2] <- SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  Matrix::forceSymmetric(M)
}
# R_trans_matrix_c1 <- function(s,x,c = 1, alpha = 2){
#   R <- matrix(nrow = 2, ncol = 2)
#   R[1,1] <- 1
#   R[1,2] <- -(alpha*s - alpha*x - (alpha*s - alpha*x)*((c + s)/(c + x))^(1/alpha))/log((c + s)/(c + x))
#   R[2,1] <- 0
#   R[2,2] <- ((c + s)/(c + x))^(1/alpha)
#   R
# }
R_trans_matrix_c1 <- function(s,x,c = 1, alpha = 2){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- ((alpha*c + alpha*x)*(c + x)^(1/alpha) + (c + s)^((alpha + 1)/alpha))/((alpha + 1)*(c + s)^(1/alpha)*c + (alpha + 1)*(c + s)^(1/alpha)*x)
  R[2,1] <- ((c + s)^(1/alpha)*(c + s) - (c + x)^((alpha + 1)/alpha))/(((alpha + 1)*c + (alpha + 1)*s)*(c + s)^(1/alpha)*x + ((alpha + 1)*c^2 + (alpha + 1)*c*s)*(c + s)^(1/alpha))
  R[1,2] <- ((alpha*c + alpha*s)*(c + s)^(1/alpha) - (alpha*c + alpha*x)*(c + x)^(1/alpha))/((alpha + 1)*(c + s)^(1/alpha))
  R[2,2] <- ((alpha*c + alpha*s)*(c + s)^(1/alpha) + (c + x)^((alpha + 1)/alpha))/(((alpha + 1)*c + (alpha + 1)*s)*(c + s)^(1/alpha))
  R
}

# case 2: alpha = 1. s > t = x:
mspline_cov_c2 <- function(s,t,c = 1){
  -1/60*(10*c*t^4 + 2*t^5 + 10*(c^2 - 2*c*s - s^2)*t^3 - 15*(2*c^2*s + c*s^2)*t^2)/((c + s)*(c + t))
}
mspline_deriv_cov_c2 <- function(s,t,c = 1){
  1/60*(40*c*t^4 + 8*t^5 + 20*(5*c^2 + 2*c*s + s^2)*t^3 + 15*(8*c^3 + 6*c^2*s + 3*c*s^2)*t^2 + 30*(2*c^4 + 2*c^3*s + c^2*s^2)*t)/(c^4 + 2*c^3*s + c^2*s^2 + (c^2 + 2*c*s + s^2)*t^2 + 2*(c^3 + 2*c^2*s + c*s^2)*t)
}
mspline_cross_cov_c2 <- function(s,t,c = 1){
  -1/60*(40*c*t^4 + 8*t^5 + 20*(3*c^2 - 2*c*s - s^2)*t^3 + 15*(2*c^3 - 6*c^2*s - 3*c*s^2)*t^2 - 30*(2*c^3*s + c^2*s^2)*t)*c/((c^3 + 2*c^2*t + c*t^2)*(c + s))
}
SS_cov_c2 <- function(s,t,c = 1){
  1/60*c^2*((20*c^2*s^3 + 25*c*s^4 + 8*s^5)/c^2 - (15*c*t^4 + 3*t^5 + 10*(2*c^2 - 2*c*s - s^2)*t^3 - 30*(2*c^2*s + c*s^2)*t^2 + 15*(4*c^2*s^2 + 4*c*s^3 + s^4)*t)/c^2)/(c + s)^2
}
SS_deriv_cov_c2 <- function(s,t,c = 1){
  1/60*(60*c^4*s + 180*c^3*s^2 + 220*c^2*s^3 + 125*c*s^4 + 28*s^5)/(c^4 + 4*c^3*s + 6*c^2*s^2 + 4*c*s^3 + s^4) - 1/60*(15*c*t^4 + 3*t^5 + 10*(4*c^2 + 2*c*s + s^2)*t^3 + 30*(2*c^3 + 2*c^2*s + c*s^2)*t^2 + 15*(4*c^4 + 8*c^3*s + 8*c^2*s^2 + 4*c*s^3 + s^4)*t)/(c^4 + 4*c^3*s + 6*c^2*s^2 + 4*c*s^3 + s^4)
}
SS_cross_cov_c2 <- function(s,t,c = 1){
  1/20*c*((10*c^3*s^2 + 20*c^2*s^3 + 15*c*s^4 + 4*s^5)/(c^3 + 2*c^2*s + c*s^2) + (10*c^3*t^2 + 10*c^2*t^3 + 5*c*t^4 + t^5 - 5*(4*c^3*s + 6*c^2*s^2 + 4*c*s^3 + s^4)*t)/(c^3 + 2*c^2*s + c*s^2))/(c + s)
}
SS_cov_mat_c2 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- SS_cov_c2(s = s, t = t, c = c)
  M[1,2] <- SS_cross_cov_c2(s = s, t = t, c = c)
  M[2,2] <- SS_deriv_cov_c2(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}
# R_trans_matrix_c2 <- function(s,x,c = 1){
#   R <- matrix(nrow = 2, ncol = 2)
#   R[1,1] <- 1
#   R[1,2] <- (s - x)^2/((c + x)*log((c + s)/(c + x)))
#   R[2,1] <- 0
#   R[2,2] <- (c + s)/(c + x)
#   R
# }
R_trans_matrix_c2 <- function(s,x,c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1/2*(2*c^2 + 2*c*s + s^2 + 2*c*x + x^2)/(c^2 + c*s + (c + s)*x)
  R[2,1] <- 1/2*(2*c*s + s^2 - 2*c*x - x^2)/(c^3 + 2*c^2*s + c*s^2 + (c^2 + 2*c*s + s^2)*x)
  R[1,2] <- 1/2*(2*c*s + s^2 - 2*c*x - x^2)/(c + s)
  R[2,2] <- 1/2*(2*c^2 + 2*c*s + s^2 + 2*c*x + x^2)/(c^2 + 2*c*s + s^2)
  R
}


# general case: automatically check alpha
mspline_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    mspline_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
mspline_deriv_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    mspline_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_deriv_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
mspline_cross_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    mspline_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    mspline_cross_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_deriv_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_deriv_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cross_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cross_cov_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
SS_cov_mat <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    SS_cov_mat_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cov_mat_c2(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
R_trans_matrix <- function(s,x,c = 1, alpha = 2){
  if(alpha != 1 & alpha != 0){
    R_trans_matrix_c1(s = s, x = x, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    R_trans_matrix_c2(s = s, x = x, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}


# Computing PSD:
PSD_compute <- function(x, h, c = 1, alpha = 2){
  sqrt(SS_cov(s = (x + h), t = (x), c = c, alpha = alpha))
}


# Computing boundary
boundary_compute <- function(x, c, alpha){
  if(alpha == 1){
    cbind(1, ln(x+c))
  }
  else{
    cbind(1, (x+c)^((alpha-1)/(alpha)))
  }
}

# Computing joint precision:
mGP_joint_prec <- function(t_vec, alpha, c){
  n <- length(t_vec)
  Blist <- list()
  AClist <- list()
  Clist <- list()

  Compute_Ci <- function(t_vec, i){
    Ci <- Matrix::forceSymmetric(solve(SS_cov_mat(t = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)))
    Ci
  }

  Compute_Ai <- function(t_vec, i, Ci) {
    Ti <- R_trans_matrix(x = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)
    t(Ti) %*% Ci %*% Ti
  }

  Compute_Bi <- function(t_vec, i, Ci) {
    Ti <- R_trans_matrix(x = t_vec[i], s = t_vec[i+1], alpha = alpha, c = c)
    -t(Ti) %*% Ci
  }

  for (i in 1:(n - 1)) {
    Clist[[i]] <- Compute_Ci(t_vec = t_vec, i = i)
  }

  for (i in 2:(n - 1)) {
    AClist[[i]] <- Compute_Ai(t_vec = t_vec, i = i, Ci = Clist[[i]]) + Clist[[i-1]]
  }

  AClist[[1]] <- Compute_Ai(t_vec = t_vec, i = 1, Ci = Clist[[1]]) + Compute_Ci(t_vec = c(0,t_vec[1]), i = 1)
  AClist[[n]] <- Compute_Ci(t_vec = t_vec, i = (n - 1))

  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(t_vec = t_vec, i = i, Ci = Clist[[i]])
  }

  Qlist <- list()
  Q <- matrix(0, nrow = 0, ncol = n*2)
  for (i in 1:(n-1)) {
    Qlist[[i]] <- cbind(matrix(0,nrow = 2, ncol = 2 * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = 2, ncol = (2 * (n-i-1))) )
    Q <- rbind(Q,Qlist[[i]])
  }
  Q <- rbind(Q,cbind(matrix(0,nrow = 2, ncol = 2 * (n-1)), AClist[[n]]))
  Q <- Matrix::forceSymmetric(Q)
  as(as.matrix(Q), "dgTMatrix")
}




