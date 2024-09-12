## double check PSD for tiwp2:
library(tidyverse)
library(Matrix)
source("code/01-state-space.R")
source("code/02-FEM.R")
source("code/03-sampling.R")
B = 3000
m <- function(x) {2*sqrt(c)*(sqrt(x+c) - sqrt(c))}
c <- 1
sd = 0.1

t_start <- 10
t_end <- t_start +  20

# set the original grid
t <- seq(t_start, t_end, by = 1)
# compute the transformed grid
t_trans <- m(t)
# sampling from tIWP2
samps1 <- sim_IWp_Var(t = t_trans, p = 2, sd = sd, initial_vec = c(m(t_start),1))
result1 <- cbind(c(t), samps1[,2])
for (i in 1:B) {
  samps1 <- sim_IWp_Var(t = t_trans, p = 2, sd = sd, initial_vec = c(m(t_start),1))
  result1 <- cbind(result1, samps1[,2])
}
# matplot(x = result1[,1], y = result1[,-1], type = 'l', ylab = "Post Samp", xlab = "x", main = "tIWP2", ylim = c(-2,10), lty = "dashed", col = "pink")
# lines(x = result1[,1], y = m(result1[,1]), col = "black", lty = "solid")


## PSD_fun:
theoretical_var_fun <- function(h, sd = 1, x){
  -(sd^2) * 8/3*(((3*h + 8)*x + 4*x^2 + 3*h + 4)*sqrt(h + x + 1) - (h^2 + (5*h + 8)*x + 4*x^2 + 5*h + 4)*sqrt(x + 1))/(sqrt(h + x + 1)*sqrt(x + 1))
}

var1 <- apply(result1[,-1], 1, var)
plot(result1[,1], var1, type = "o", xlab = "x", ylab = "Var", main = "PSD for tIWP2")
lines(result1[,1], (theoretical_var_fun(h = (result1[,1]-t_start), x = t_start, sd = sd)), col = "red", lty = "dashed")






t_start <- 1000
h <- 10
t <- c(t_start, t_start + h)
t_trans <- m(t)
samps1 <- sim_IWp_Var(t = t_trans, p = 2, sd = sd, initial_vec = c(m(t_start),1))
result1 <- cbind(c(t), samps1[,2])
for (i in 1:B) {
  samps1 <- sim_IWp_Var(t = t_trans, p = 2, sd = sd, initial_vec = c(m(t_start),1))
  result1 <- cbind(result1, samps1[,2])
}
var(result1[2,-1])
theoretical_var_fun(h = h, x = t_start, sd = sd)



