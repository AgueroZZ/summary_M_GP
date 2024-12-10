# case1 : alpha being positive but not 1 or 2 or 1/2
adj_SS_cov_c1 <- function(s,t,c = 1, alpha = 2){
  1/3*(9*alpha^3*(c - s)^(2/alpha)*c^3 + 6*(alpha^4 - 2*alpha^3 + alpha^2)*(c - s)^(2/alpha)*c^2*s - 6*(alpha^4 - 2*alpha^3 + alpha^2)*(c - s)^(2/alpha)*c*s^2 + 2*(alpha^4 - 2*alpha^3 + alpha^2)*(c - s)^(2/alpha)*s^3 - (2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c - s)^(2/alpha)*t^3 + 3*((2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c - s)^(2/alpha)*c - 2*((alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*c - (alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*s)*(c - t)^(1/alpha))*t^2 + 3*((2*alpha^4 + alpha^3)*c^3 - 2*(2*alpha^4 + alpha^3)*c^2*s + (2*alpha^4 + alpha^3)*c*s^2)*(c - t)^(2/alpha) - 6*((alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*c^3 - (alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*c^2*s)*(c - t)^(1/alpha) - 3*((2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c - s)^(2/alpha)*c^2 + ((2*alpha^4 + alpha^3)*c^2 - 2*(2*alpha^4 + alpha^3)*c*s + (2*alpha^4 + alpha^3)*s^2)*(c - t)^(2/alpha) - 4*((alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*c^2 - (alpha^4 + 2*alpha^3)*(c - s)^(1/alpha)*c*s)*(c - t)^(1/alpha))*t)/((2*alpha^4 + alpha^3 - 6*alpha^2 + alpha + 2)*(c - s)^(2/alpha))
}
adj_SS_deriv_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -2*(alpha*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c - alpha*(c - s)^(1/alpha)*(c - t)^(1/alpha)*t - (alpha*c - alpha*s)*(c - t)^(2/alpha))/((alpha - 1)*(c - s)^(2/alpha))
}
adj_SS_cross_cov_c1 <- function(s,t,c = 1, alpha = 2){
  ((alpha^3 + 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c^2 + (alpha^3 + 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*t^2 + (alpha^3 - alpha^2)*(c - s)^(2/alpha)*c^2 - 2*(alpha^3 - alpha^2)*(c - s)^(2/alpha)*c*s + (alpha^3 - alpha^2)*(c - s)^(2/alpha)*s^2 - ((2*alpha^3 + alpha^2)*c^2 - (2*alpha^3 + alpha^2)*c*s)*(c - t)^(2/alpha) - (2*(alpha^3 + 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c - ((2*alpha^3 + alpha^2)*c - (2*alpha^3 + alpha^2)*s)*(c - t)^(2/alpha))*t)/((2*alpha^3 + 3*alpha^2 - 3*alpha - 2)*(c - s)^(2/alpha))
}
adj_SS_cov_mat_c1 <- function(s,t,c = 1, alpha = 2){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[1,2] <- adj_SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  M[2,2] <- adj_SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  Matrix::forceSymmetric(M)
}
adj_R_trans_matrix_c1 <- function(s,x,c = 1, alpha = 2){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- (alpha*(c - s)^(1/alpha)*c - alpha*(c - s)^(1/alpha)*x - (alpha*c - alpha*s)*(c - x)^(1/alpha))/((alpha - 1)*(c - s)^(1/alpha))
  R[2,2] <- (c - x)^(1/alpha)/(c - s)^(1/alpha)
  R
}

# case 2: alpha = 1. s > t = x:
adj_SS_cov_c2 <- function(s,t,c = 1){
  1/12*(4*c*s^3 - 3*s^4 - 4*c*t^3 + t^4 + 6*(2*c*s - s^2)*t^2 - 4*(3*c*s^2 - 2*s^3)*t)/(c - t)
}
adj_SS_deriv_cov_c2 <- function(s,t,c = 1){
  (c*s - s^2 - (c - s)*t)/(c - t)
}
adj_SS_cross_cov_c2 <- function(s,t,c = 1){
  1/2*(c*s^2 - s^3 + (c - s)*t^2 - 2*(c*s - s^2)*t)/(c - t)
}
adj_SS_cov_mat_c2 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c2(s = s, t = t, c = c)
  M[1,2] <- adj_SS_cross_cov_c2(s = s, t = t, c = c)
  M[2,2] <- adj_SS_deriv_cov_c2(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}
adj_R_trans_matrix_c2 <- function(s,x,c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- 1/2*(2*c*s - s^2 - 2*c*x + x^2)/(c - x)
  R[2,2] <- (c - s)/(c - x)
  R
}

# case 3: alpha < 0 but not -1
adj_SS_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  1/3*(3*(2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*c^3 - 6*(2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*c^2*s + 3*(2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*c*s^2 - (2*alpha^4 - 5*alpha^3 + 2*alpha^2)*(c - t)^(2/alpha)*t^3 + 3*((2*alpha^4 - 5*alpha^3 + 2*alpha^2)*(c - t)^(2/alpha)*c - 2*((alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*c - (alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*s)*(c - t)^(1/alpha))*t^2 - (9*alpha^3*c^3 - 6*(alpha^4 + 2*alpha^3 + alpha^2)*c^2*s + 6*(alpha^4 + 2*alpha^3 + alpha^2)*c*s^2 - 2*(alpha^4 + 2*alpha^3 + alpha^2)*s^3)*(c - t)^(2/alpha) - 6*((alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*c^3 - (alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*c^2*s)*(c - t)^(1/alpha) - 3*((2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*c^2 + (2*alpha^4 - 5*alpha^3 + 2*alpha^2)*(c - t)^(2/alpha)*c^2 - 2*(2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*c*s + (2*alpha^4 - alpha^3)*(c - s)^(2/alpha)*s^2 - 4*((alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*c^2 - (alpha^4 - 2*alpha^3)*(c - s)^(1/alpha)*c*s)*(c - t)^(1/alpha))*t)/((2*alpha^4 - alpha^3 - 6*alpha^2 - alpha + 2)*(c - t)^(2/alpha))
}
adj_SS_deriv_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -2*(alpha*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c - alpha*(c - s)^(1/alpha)*(c - t)^(1/alpha)*t - alpha*(c - s)^(2/alpha)*c + alpha*(c - s)^(2/alpha)*s)/((alpha + 1)*(c - t)^(2/alpha))
}
adj_SS_cross_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  ((alpha^3 - 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c^2 + (alpha^3 - 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*t^2 - (2*alpha^3 - alpha^2)*(c - s)^(2/alpha)*c^2 + (2*alpha^3 - alpha^2)*(c - s)^(2/alpha)*c*s + ((alpha^3 + alpha^2)*c^2 - 2*(alpha^3 + alpha^2)*c*s + (alpha^3 + alpha^2)*s^2)*(c - t)^(2/alpha) - (2*(alpha^3 - 2*alpha^2)*(c - s)^(1/alpha)*(c - t)^(1/alpha)*c - (2*alpha^3 - alpha^2)*(c - s)^(2/alpha)*c + (2*alpha^3 - alpha^2)*(c - s)^(2/alpha)*s)*t)/((2*alpha^3 - 3*alpha^2 - 3*alpha + 2)*(c - t)^(2/alpha))
}
adj_SS_cov_mat_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c3(s = s, t = t, c = c, alpha = alpha)
  M[1,2] <- adj_SS_cross_cov_c3(s = s, t = t, c = c, alpha = alpha)
  M[2,2] <- adj_SS_deriv_cov_c3(s = s, t = t, c = c, alpha = alpha)
  Matrix::forceSymmetric(M)
}
adj_R_trans_matrix_c3 <- function(s,x,c = 1, alpha = -3){
  alpha <- abs(alpha)
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- -(alpha*(c - s)^(1/alpha)*c - alpha*(c - x)^(1/alpha)*c - alpha*(c - s)^(1/alpha)*s + alpha*(c - x)^(1/alpha)*x)/((alpha + 1)*(c - x)^(1/alpha))
  R[2,2] <- (c - s)^(1/alpha)/(c - x)^(1/alpha)
  R
}

# case 4: alpha = -1
adj_SS_cov_c4 <- function(s,t,c = 1){
  1/12*(4*c*s^3 - 3*s^4 - 4*c*t^3 + t^4 + 6*(2*c*s - s^2)*t^2 - 4*(3*c*s^2 - 2*s^3)*t)/(c - t)
}
adj_SS_deriv_cov_c4 <- function(s,t,c = 1){
  (c*s - s^2 - (c - s)*t)/(c - t)
}
adj_SS_cross_cov_c4 <- function(s,t,c = 1){
  1/2*(c*s^2 - s^3 + (c - s)*t^2 - 2*(c*s - s^2)*t)/(c - t)
}
adj_SS_cov_mat_c4 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c4(s = s, t = t, c = c)
  M[1,2] <- adj_SS_cross_cov_c4(s = s, t = t, c = c)
  M[2,2] <- adj_SS_deriv_cov_c4(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}
adj_R_trans_matrix_c4 <- function(s,x,c = 1){
  R <- matrix(nrow = 2, ncol = 2)
  R[1,1] <- 1
  R[2,1] <- 0
  R[1,2] <- 1/2*(2*c*s - s^2 - 2*c*x + x^2)/(c - x)
  R[2,2] <- (c - s)/(c - x)
  R
}

# case 5: alpha = 2 (root # 1)
adj_SS_cov_c5 <- function(s,t,c = 1){
  # 4/9*s^3*(log(c - s) - 1) - 4/9*c^3*log(c - s) + 16/27*c^3 - 4/3*(c*log(c - s) - c)*s^2 - 4/9*c^2*t + 4/9*c*t^2 - 4/27*t^3 - 16/27*(c^2 - c*s - (c - s)*t)*sqrt(c - s)*sqrt(c - t) + 4/3*(c^2*log(c - s) - c^2)*s + 4/9*(c^3 - 3*c^2*s + 3*c*s^2 - s^3)*log(c - t)
  2/15*(24*c^4 - 12*c^3*s - 3*c^2*s^2 + c*s^3 - 10*c*t^3 + 15*(3*c^2 - c*s)*t^2 - 24*(c^3 - 2*c^2*t + c*t^2)*sqrt(c - s)*sqrt(c - t) - 30*(2*c^3 - c^2*s)*t)/c
}
adj_SS_deriv_cov_c5 <- function(s,t,c = 1){
  # -c*log(c - s) + s*log(c - s) + (c - s)*log(c - t)
  1/2*(2*c*s - s^2 - 2*c*t + t^2)/(c - s)
}
adj_SS_cross_cov_c5 <- function(s,t,c = 1){
  # 2/9*s^2*(3*log(c - s) - 2) + 2/3*c^2*log(c - s) + 4/9*sqrt(c - s)*(c - t)^(3/2) - 4/9*c^2 - 4/9*(3*c*log(c - s) - 2*c)*s - 2/3*(c^2 - 2*c*s + s^2)*log(c - t)
  1/5*((c^2 - 2*c*s + s^2)*sqrt(c - s)*sqrt(c) - 5*(c^2 - 2*c*t + t^2)*sqrt(c - s)*sqrt(c) + 4*(c^2 - 2*c*t + t^2)*sqrt(c - t)*sqrt(c))/(sqrt(c - s)*sqrt(c))
}
adj_SS_cov_mat_c5 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c5(s = s, t = t, c = c)
  M[1,2] <- adj_SS_cross_cov_c5(s = s, t = t, c = c)
  M[2,2] <- adj_SS_deriv_cov_c5(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}


# case 6: alpha = 1/2 (root # 2)
adj_SS_cov_c6 <- function(s,t,c = 1){
  NA
}
adj_SS_deriv_cov_c6 <- function(s,t,c = 1){
  NA
}
adj_SS_cross_cov_c6 <- function(s,t,c = 1){
  NA
}
adj_SS_cov_mat_c6 <- function(s,t,c = 1){
  M <- matrix(nrow = 2, ncol = 2)
  M[1,1] <- adj_SS_cov_c6(s = s, t = t, c = c)
  M[1,2] <- adj_SS_cross_cov_c6(s = s, t = t, c = c)
  M[2,2] <- adj_SS_deriv_cov_c6(s = s, t = t, c = c)
  Matrix::forceSymmetric(M)
}

# general case: automatically check alpha
adj_SS_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha > 0){
    adj_SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    adj_SS_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    adj_SS_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    adj_SS_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == 2){
    adj_SS_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == 1/2){
    adj_SS_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
adj_SS_deriv_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha > 0){
    adj_SS_deriv_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    adj_SS_deriv_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    adj_SS_deriv_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    adj_SS_deriv_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == 2){
    adj_SS_deriv_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == 1/2){
    adj_SS_deriv_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
adj_SS_cross_cov <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha > 0 & alpha != 2 & alpha != (1/2)){
    adj_SS_cross_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    adj_SS_cross_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1){
    adj_SS_cross_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    adj_SS_cross_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == 2){
    adj_SS_cross_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == 1/2){
    adj_SS_cross_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
adj_SS_cov_mat <- function(s,t,c = 1, alpha = 2){
  if(alpha != 1 & alpha > 0 & alpha != 2 & alpha != (1/2)){
    adj_SS_cov_mat_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    adj_SS_cov_mat_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1){
    adj_SS_cov_mat_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    adj_SS_cov_mat_c4(s = s, t = t, c = c)
  }
  else if(alpha == 2){
    adj_SS_cov_mat_c5(s = s, t = t, c = c)
  }
  else if(alpha == 1/2){
    adj_SS_cov_mat_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}
adj_R_trans_matrix <- function(s,x,c = 1, alpha = 2){
  if(alpha != 1 & alpha > 0){
    adj_R_trans_matrix_c1(s = s, x = x, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    adj_R_trans_matrix_c2(s = s, x = x, c = c)
  }
  else if(alpha < 0 & alpha != -1){
    adj_R_trans_matrix_c3(s = s, x = x, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    adj_R_trans_matrix_c4(s = s, x = x, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}




