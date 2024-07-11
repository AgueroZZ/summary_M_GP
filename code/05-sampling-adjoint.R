adj_mGP_sim <- function(t = NULL, mesh_size = 0.01, max_t = 10, alpha = 2, sd = 1, c = 1, initial_vec = NULL){
  if(is.null(initial_vec)){
    initial_vec <- c(0,0)
  }
  result <- matrix(data = initial_vec, nrow = 1, ncol = 2)
  if (is.null(t)) {
    t <- seq(0, max_t, by = mesh_size)
  }
  for (i in 1:(length(t)-1)) {
    Trans <- adj_R_trans_matrix(s = t[i+1],x = t[i], c = c, alpha = alpha)
    Sig <- (sd^2) * as.matrix(adj_SS_cov_mat(s = t[i+1],t = t[i], c = c, alpha = alpha))
    if(sd > 0){
      result_new <- Trans %*% result[i, ] + as.numeric(LaplacesDemon::rmvn(mu = rep(0,2), Sigma = Sig))
    }
    else{
      result_new <- Trans %*% result[i, ]
    }
    result <- rbind(result, matrix(data = result_new, nrow = 1))
  }
  result <- cbind(t, result)
  colnames(result) <- c("t", "func", "func_1st")
  result
}
