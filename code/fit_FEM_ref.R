fit_mGP_once_FEM_ref <- function(data_sim, data_train,
                             u, betaprec = 0.001,
                             k = 30, a, c, ref_location = NULL,
                             region, accuracy = 0.01, boundary = TRUE
                             ){

  if(is.null(ref_location)){
    ref_location <- median(data_sim$x)
  }
  if(ref_location == 0){
    stop("If ref_location is 0, please use the regular version of the function.")
  }
  if(ref_location >= max(data_sim$x)){
    stop("If ref_location is the maximum of data_sim$x, please use the regular version of the function.")
  }

  # Update c based on ref_location
  c <- c + ref_location

  # IMPORTANT:
  # It is highly suggested to reparametrize m(x) so m(0) = 0 and m'(0) = 1.
  # this helps the interpretation of the prior.
  lambda <- (a-1)/a
  m <- function(x){
    ((x + c)^lambda - c^lambda)/(lambda * c^(lambda - 1))
  }

  train_x <- data_train$x
  all_x <- data_sim$x

  # construct fixed design matrix (aka boundary condition)
  X <- as(cbind(1, m(data_train$x)), "dgCMatrix")

  # construct random design matrix
  # Forward design matrix
  train_x_positive <- train_x - ref_location
  all_x_positive <- all_x - ref_location
  train_x_positive <- ifelse(train_x_positive < 0, 0, train_x_positive)
  all_x_positive <- ifelse(all_x_positive < 0, 0, all_x_positive)
  B_pos <- as(Compute_Design(x = train_x_positive, k, region = range(all_x_positive)), "dgTMatrix")

  # Backward design matrix
  train_x_negative <- ref_location - train_x
  all_x_negative <- ref_location - all_x
  train_x_negative <- ifelse(train_x_negative < 0, 0, train_x_negative)
  all_x_negative <- ifelse(all_x_negative < 0, 0, all_x_negative)
  B_neg <- as(Compute_Design(x = train_x_negative, k, region = range(all_x_negative)), "dgTMatrix")

  # Combine into a matrix
  B <- as(cbind(B_pos, B_neg), "dgCMatrix")

  # construct penalty matrix:
  # Forwards
  P_pos <- Compute_Prec(k = k, region = range(all_x_positive), a = a, c = c)

  # Backwards
  P_neg <- Compute_Prec_rev(k = k, region = range(all_x_negative), a = a, c = c)

  # Combine into a block diagonal matrix
  P <- as(bdiag(P_pos, P_neg), "dgCMatrix")
  logPdet <- determinant(P)$modulus

  # fit the model
  tmbdat <- list(
    y = data_train$y,
    X = X,
    B = B,
    P = P,
    logPdet = as.numeric(logPdet),
    betaprec = betaprec,
    sig = sd_noise,
    u = u,
    alpha = 0.5
  )

  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B)))),
    theta = 0
  )

  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "fitGP_known_sd",
    random = "W",
    silent = TRUE
  )

  ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
  fit <- aghq::marginal_laplace_tmb(ff, k = 4, startingvalue = 0)
  return(fit)
}

sample_model_once_FEM_ref <- function(k, data_sim, model_fit, M = 3000,
                                      a, c, ref_location = NULL){

  if(is.null(ref_location)){
    ref_location <- median(data_sim$x)
  }
  if(ref_location == 0){
    stop("If ref_location is 0, please use the regular version of the function.")
  }
  if(ref_location >= max(data_sim$x)){
    stop("If ref_location is the maximum of data_sim$x, please use the regular version of the function.")
  }

  samps <- aghq::sample_marginal(quad = model_fit, M = M)
  basis_weights_samps <- samps$samps[1:(nrow(samps$samps)-2),]
  beta_samps <- samps$samps[(nrow(samps$samps)-1):nrow(samps$samps),]

  c <- c + ref_location
  lambda <- (a-1)/a

  m <- function(x) {
    ((x + c)^lambda - c^lambda)/(lambda * c^(lambda - 1))
  }
  X_refined <- as(cbind(1, m(data_sim$x)), "dgCMatrix")

  all_x <- data_sim$x
  all_x_positive <- all_x - ref_location
  all_x_positive <- ifelse(all_x_positive < 0, 0, all_x_positive)
  B_pos <- as(Compute_Design(x = all_x_positive, k, region = range(all_x_positive)), "dgTMatrix")

  all_x_negative <- ref_location - all_x
  all_x_negative <- ifelse(all_x_negative < 0, 0, all_x_negative)
  B_neg <- as(Compute_Design(x = all_x_negative, k, region = range(all_x_negative)), "dgTMatrix")

  B <- as(cbind(B_pos, B_neg), "dgCMatrix")

  f_samps_combined <- X_refined %*% beta_samps + B %*% basis_weights_samps
  f_samps_combined
}


# define data_train by subsampling data_sim
# the probability of being sampled is higher for larger data_sim$x
prob_weights <- exp(0.7*data_sim$x)
prob_weights <- prob_weights/sum(prob_weights)
data_train <- data_sim[sample(1:nrow(data_sim), 30, prob = prob_weights, replace = F)  ,]
data_train <- rbind(data_train, data.frame(x = 0.1, y = 0))
plot(data_train$x, data_train$y, type = "p", col = "black",
     lwd = 1, pch = 1, cex = 0.2,
     ylim = c(0,100), ylab = "y", xlab = "x")


ref_location_choice <- 8
a <- 2.3

mGP_FEM <- fit_mGP_once_FEM_ref(data_sim, data_train, u = 30, k = 30,
                                ref_location = ref_location_choice,
                                betaprec = 1e-7,
                                a = a, c = c)

samples_mGP_FEM <- sample_model_once_FEM_ref(k = 30, data_sim = data_sim,
                                  model_fit = mGP_FEM, M = 3000, a = a, c = c,
                                  ref_location = ref_location_choice)


# Plot the samples
mGP_FEM_summary <- data.frame(x = data_sim$x, mean = rowMeans(samples_mGP_FEM),
                              lower = apply(samples_mGP_FEM, 1, quantile, probs = 0.025),
                              upper = apply(samples_mGP_FEM, 1, quantile, probs = 0.975))

plot(data_sim$x, data_sim$y, type = "n", col = "black",
     lwd = 1, pch = 1, cex = 0.2,
     ylim = range(c(mGP_FEM_summary$lower, rev(mGP_FEM_summary$upper))), ylab = "y", xlab = "x")
points(data_train$x, data_train$y, col = "black", pch = 1, cex = 0.2)
lines(mGP_FEM_summary$x, mGP_FEM_summary$mean, col = "red", lwd = 2)
lines(mGP_FEM_summary$x, f(mGP_FEM_summary$x), col = "blue", lwd = 2, lty = 2)
matlines(x = mGP_FEM_summary$x, y = samples_mGP[,2:5], lty = 1, col = "blue", lwd = 0.5)
polygon(c(mGP_FEM_summary$x, rev(mGP_FEM_summary$x)), c(mGP_FEM_summary$lower, rev(mGP_FEM_summary$upper)), col = rgb(1, 0, 0, 0.2), border = NA)
abline(v = ref_location_choice, col = "green", lwd = 2, lty = 2)


