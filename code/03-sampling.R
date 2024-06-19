sampling_from_FEM <- function(x, a, k, region, boundary = TRUE, n = 1){
  Prec <- Compute_Prec(a, k, region, boundary = boundary)
  B <- Compute_Design(x, k, region, boundary = boundary)
  if(boundary){
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,((k-2))), Omega = as.matrix(Prec))
  }
  else{
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,(k)), Omega = as.matrix(Prec))
  }
  splfd <- B %*% t(coefs_samps)
  (splfd)
}

mGP_sim <- function(t = NULL, mesh_size = 0.01, max_t = 10, alpha = 2, sd = 1, c = 1){
  result <- matrix(data = 0, nrow = 1, ncol = 2)
  if (is.null(t)) {
    t <- seq(0, max_t, by = mesh_size)
  }
  for (i in 1:(length(t)-1)) {
    Trans <- R_trans_matrix(s = t[i+1],x = t[i], c = c, alpha = alpha)
    Sig <- (sd^2) * as.matrix(SS_cov_mat(s = t[i+1],t = t[i], c = c, alpha = alpha))
    result_new <- Trans%*% result[i, ] + as.numeric(LaplacesDemon::rmvn(mu = rep(0,2), Sigma = Sig))
    result <- rbind(result, matrix(data = result_new, nrow = 1))
  }
  result <- cbind(t, result)
  colnames(result) <- c("t", "func", "func_1st")
  result
}
