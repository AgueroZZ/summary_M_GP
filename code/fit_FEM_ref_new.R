# Helper function to prepare common data
prepare_fem_data <- function(data_sim, data_train, k, a, c, ref_location = NULL) {
  # Set ref_location if not provided
  if (is.null(ref_location)) {
    ref_location <- median(data_sim$x)
  }

  # Update c and compute lambda
  c <- c + ref_location
  lambda <- (a - 1) / a

  # Fixed design matrix
  m <- function(x){
    ((x + c)^lambda - c^lambda)/(lambda * c^(lambda - 1))
  }
  # Fixed design matrix for training data
  X_train <- Matrix::sparse.model.matrix(~ 1 + m(data_train$x))

  # Prepare forward and backward design matrices for training data
  train_x_pos <- pmax(data_train$x - ref_location, 0)
  all_x_pos <- pmax(data_sim$x - ref_location, 0)

  train_x_neg <- pmax(ref_location - data_train$x, 0)
  all_x_neg <- pmax(ref_location - data_sim$x, 0)

  # Define B and penalty matrices based on non-zero regions for training data
  if (all(all_x_pos == 0)) {
    B_train <- Compute_Design(train_x_neg, k, region = range(all_x_neg))
    P <- Compute_Prec_rev(k, region = range(all_x_neg), a = a, c = c)
    B_sim <- Compute_Design(all_x_neg, k, region = range(all_x_neg))

  } else if (all(all_x_neg == 0)) {
    B_train <- Compute_Design(train_x_pos, k, region = range(all_x_pos))
    P <- Compute_Prec(k, region = range(all_x_pos), a = a, c = c)
    B_sim <- Compute_Design(all_x_pos, k, region = range(all_x_pos))
  } else {
    B_pos <- Compute_Design(train_x_pos, k, region = range(all_x_pos))
    B_neg <- Compute_Design(train_x_neg, k, region = range(all_x_neg))
    B_train <- Matrix::cbind2(B_pos, B_neg)

    B_pos <- Compute_Design(all_x_pos, k, region = range(all_x_pos))
    B_neg <- Compute_Design(all_x_neg, k, region = range(all_x_neg))
    B_sim <- Matrix::cbind2(B_pos, B_neg)

    # Penalty matrices
    P_pos <- Compute_Prec(k, region = range(all_x_pos), a = a, c = c)
    P_neg <- Compute_Prec_rev(k, region = range(all_x_neg), a = a, c = c)
    P <- Matrix::bdiag(P_pos, P_neg)
  }

  # Fixed design matrix for simulation data
  X_sim <- Matrix::sparse.model.matrix(~ 1 + m(data_sim$x))

  # Convert all matrices to "dgTMatrix" format before returning
  return(list(
    X_train = as(X_train, "dgTMatrix"),
    B_train = as(B_train, "dgTMatrix"),
    X_sim = as(X_sim, "dgTMatrix"),
    B_sim = as(B_sim, "dgTMatrix"),
    P = as(P, "dgTMatrix"),
    logPdet = as.numeric(determinant(P)$modulus),
    ref_location = ref_location
  ))
}

# Updated fit_mGP_once_FEM_ref function
fit_mGP_once_FEM_ref <- function(data_sim, data_train, u, betaprec = 0.001,
                                 k = 30, a, c, ref_location = NULL, accuracy = 0.01, boundary = TRUE) {
  # Prepare shared data
  fem_data <- prepare_fem_data(data_sim, data_train, k, a, c, ref_location)

  # Model fitting data
  tmbdat <- list(
    y = data_train$y, X = fem_data$X_train, B = fem_data$B_train, P = fem_data$P,
    logPdet = fem_data$logPdet, betaprec = betaprec, sig = sd_noise,
    u = u, alpha = 0.5
  )

  tmbparams <- list(W = numeric(ncol(fem_data$X_train) + ncol(fem_data$B_train)), theta = 0)

  # Fit model
  ff <- TMB::MakeADFun(data = tmbdat, parameters = tmbparams, DLL = "fitGP_known_sd",
                       random = "W", silent = TRUE)

  ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
  fit <- aghq::marginal_laplace_tmb(ff, k = 4, startingvalue = 0)

  # Return fit along with pre-computed data for sampling
  return(list(fit = fit, X_sim = fem_data$X_sim, B_sim = fem_data$B_sim, ref_location = fem_data$ref_location))
}

# Updated sample_model_once_FEM_ref function
sample_model_once_FEM_ref <- function(model_fit, M = 3000) {
  # Extract fit and simulation design matrices
  fit <- model_fit$fit
  X_sim <- model_fit$X_sim
  B_sim <- model_fit$B_sim

  # Sampling
  samps <- aghq::sample_marginal(quad = fit, M = M)
  beta_samps <- samps$samps[(nrow(samps$samps) - 1):nrow(samps$samps),]
  basis_weights_samps <- samps$samps[1:(nrow(samps$samps) - 2),]

  # Use pre-computed matrices for sampling on data_sim
  f_samps_combined <- X_sim %*% beta_samps + B_sim %*% basis_weights_samps

  return(f_samps_combined)
}

# Test the updated functions
ref_location_choice <- 8
a <- 2.3
mGP_FEM <- fit_mGP_once_FEM_ref(data_sim = data_sim, data_train = data_train,
                                u = 30, k = 30,
                                ref_location = ref_location_choice,
                                betaprec = 1e-7,
                                a = a, c = c)

samples_mGP_FEM <- sample_model_once_FEM_ref(model_fit = mGP_FEM, M = 3000)

mGP_FEM_summary <- data.frame(x = data_sim$x, mean = rowMeans(samples_mGP_FEM),
                              lower = apply(samples_mGP_FEM, 1, quantile, probs = 0.025),
                              upper = apply(samples_mGP_FEM, 1, quantile, probs = 0.975))

plot(data_sim$x, data_sim$y, type = "n", col = "black",
     lwd = 1, pch = 1, cex = 0.2,
     ylim = range(c(mGP_FEM_summary$lower, rev(mGP_FEM_summary$upper))), ylab = "y", xlab = "x")
points(data_train$x, data_train$y, col = "black", pch = 1, cex = 0.2)
lines(mGP_FEM_summary$x, mGP_FEM_summary$mean, col = "red", lwd = 2)
lines(mGP_FEM_summary$x, f(mGP_FEM_summary$x), col = "blue", lwd = 2, lty = 2)
matlines(x = mGP_FEM_summary$x, y = samples_mGP_FEM[,2:5], lty = 1, col = "blue", lwd = 0.5)
polygon(c(mGP_FEM_summary$x, rev(mGP_FEM_summary$x)), c(mGP_FEM_summary$lower, rev(mGP_FEM_summary$upper)), col = rgb(1, 0, 0, 0.2), border = NA)
abline(v = ref_location_choice, col = "green", lwd = 2, lty = 2)



