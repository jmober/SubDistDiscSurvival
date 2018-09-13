##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Simulations                     #####
#####                      presented in Section 3                        ##### 
##############################################################################

# The following contains R-Code to compute and draw the cumulative incidence 
# functions presented in the second part simulation study. The second part of 
# the simulation study compares (i) the proposed discrete subdistribution model, 
# and (ii) the continuous Fine & Gray model.  

# load R package and necessary files 
source("functions_II.R")
library("discSurv")
library("cmprsk")

# function to compute the (mean) true cumulative incidence function 
sim_F <- function(X, coefs){
  eta  <- c(as.matrix(X)%*%coefs)
  Ff   <- sapply(limits, function(t) cdf1(t, p, eta))
  Ff_m <- apply(Ff, 2, mean)
  return(Ff_m)
}

# function to compute the (mean) fitted discrete cumulative incidence function 
one_F <- function(X, coefs){
  eta    <- c(as.matrix(X)%*%coefs[n_limits:(n_limits+3)])
  eta    <- matrix(rep(coefs[1:(n_limits-1)], each = 500), nrow = 500) + eta 
  lambda <- 1 - exp(-exp(eta))
  S      <- t(apply(lambda, 1, function(x) cumprod(1-x)))
  Ff     <- 1-S
  Ff_m   <- apply(Ff, 2, mean)
  return(Ff_m)
}

#####
#####

# Figure 10(a). The boxplots visualize the estimates of the cumulative
# incidence function F1, as obtained from the discrete-time subdistribution hazard model 
# and the continuous-time Fine & Gray model (n = 500, q = 0.8, b = 0.85, k = 4).

# initialize 
n_limits <- 4
betas1   <- c(0.4, -0.4, 0.2, -0.2)
betas2   <- c(-0.4, 0.4, -0.2, 0.2)
n        <- 500
p        <- 0.8
b        <- 0.85
limits <- get_limits(p, n_limits)
F_true <- matrix(NA, nrow = 1000, ncol = 3)
F_disc <- matrix(NA, nrow = 1000, ncol = 3)
F_cont <- matrix(NA, nrow = 1000, ncol = 3)

# compute results 
for(i in 1:1000){
  
  # draw X
  set.seed(8617+i)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rbinom(n, 1, 0.5)
  X  <- cbind(x1, x2, x3, x4)
    
  # generate data set 
  data_disc <- all_disc_data(X, betas1, betas2, p, b, limits)
  data_disc$e1[data_disc$time == n_limits]   <- 0 
  data_disc$e2[data_disc$time == n_limits]   <- 0 
  data_disc$time[data_disc$time == n_limits] <- n_limits-1
  
  # compute true F 
  F_true[i,] <- sim_F(data_disc[ ,6:9], betas1)
  
  # fit discrete subdistribution model 
  data_disc_long_sub <- dataLongSubDist(dataSet = data_disc, timeColumn = "time", eventColumns = c("e1", "e2"), eventFocus = "e1")
  model_sub          <- glm(y ~ timeInt - 1 + x1 + x2 + x3 + x4, 
                            family = binomial(link = "cloglog"), 
                            data = data_disc_long_sub, 
                            weights = data_disc_long_sub$subDistWeights)
  # call glm() with the subDistWeights computed by dataLongSubDist() 
  # passed to the weights-argument
  est_sub    <- coef(model_sub)
  # compute discrete F 
  F_disc[i,] <- one_F(data_disc[ ,6:9], est_sub)
    
  # fit continuous Fine & Gray model 
  model_cont <- crr(data_disc$time, data_disc$status, X, failcode = 1, cencode = 0) # call crr 
  est_cont   <- summary(model_cont)$coef[ ,1]
  # compute continuous F 
  pred_cont  <- t(predict(model_cont, as.matrix(data_disc[ ,6:9])))[-1, ]
  F_cont[i,] <- apply(pred_cont, 2, mean)

}

# draw figure 
F_true <- apply(F_true, 2, mean)
par(mar=c(5, 6, 1, 1))
boxplot(F_cont[ ,1], F_disc[ ,1], F_cont[ ,2], F_disc[ ,2], F_cont[ ,3], F_disc[ ,3], col=rep(c(grey(0.5), grey(0.8)), 3), 
        pch = 21, bg = rep(c(grey(0.5), grey(0.8)),3), xaxt = "n",cex.lab = 1.8, cex.axis = 1.5, xlab = "t")
mtext(bquote(hat(F)[1]), 2, line = 3, cex = 1.8)
axis(1, c(1, 2, 3), at=c(1.5, 3.5, 5.5), cex.axis = 1.5)
lines(c(0.5, 2.5),c(F_true[1], F_true[1]), col = "red", lwd = 1.5)
lines(c(2.5, 4.5),c(F_true[2], F_true[2]), col = "red", lwd = 1.5)
lines(c(4.5, 6.5),c(F_true[3], F_true[3]), col = "red", lwd = 1.5)
legend("topleft", legend = c("continuous", "discrete"), lty = "solid", col = c(grey(0.5), grey(0.8)), bty = "n", lwd = 3, cex = 1.5)


#####
#####

# Figure 10(b). The boxplots visualize the estimates of the cumulative
# incidence function F1, as obtained from the discrete-time subdistribution hazard model 
# and the continuous-time Fine & Gray model (n = 500, q = 0.8, b = 0.85, k = 8).

# initialize 
n_limits <- 8
betas1   <- c(0.4, -0.4, 0.2, -0.2)
betas2   <- c(-0.4, 0.4, -0.2, 0.2)
n        <- 500
p        <- 0.8
b        <- 0.85
limits <- get_limits(p, n_limits)
F_true <- matrix(NA, nrow = 1000, ncol = 7)
F_disc <- matrix(NA, nrow = 1000, ncol = 7)
F_cont <- matrix(NA, nrow = 1000, ncol = 7)

# compute results 
for(i in 1:1000){
  
  # draw X
  set.seed(8617+i)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rbinom(n, 1, 0.5)
  X  <- cbind(x1, x2, x3, x4)
  
  # generate data set 
  data_disc <- all_disc_data(X, betas1, betas2, p, b, limits)
  data_disc$e1[data_disc$time == n_limits]   <- 0 
  data_disc$e2[data_disc$time == n_limits]   <- 0 
  data_disc$time[data_disc$time == n_limits] <- n_limits-1
  
  # compute true F 
  F_true[i,] <- sim_F(data_disc[ ,6:9], betas1)
  
  # fit discrete subdistribution model 
  data_disc_long_sub <- dataLongSubDist(dataSet = data_disc, timeColumn = "time", eventColumns = c("e1", "e2"), eventFocus = "e1")
  model_sub          <- glm(y ~ timeInt - 1 + x1 + x2 + x3 + x4, 
                            family = binomial(link = "cloglog"), 
                            data = data_disc_long_sub, 
                            weights = data_disc_long_sub$subDistWeights)
  # call glm() with the subDistWeights computed by dataLongSubDist() 
  # passed to the weights-argument
  est_sub    <- coef(model_sub)
  # compute discrete F 
  F_disc[i,] <- one_F(data_disc[ ,6:9], est_sub)
  
  # fit continuous Fine & Gray model 
  model_cont <- crr(data_disc$time, data_disc$status, X, failcode = 1, cencode = 0) # call crr 
  est_cont   <- summary(model_cont)$coef[ ,1]
  # compute continuous F 
  pred_cont  <- t(predict(model_cont, as.matrix(data_disc[ ,6:9])))[-1, ]
  F_cont[i,] <- apply(pred_cont, 2, mean)
  
}

# draw figure 
F_true <- apply(F_true, 2, mean)
par(mar=c(5, 6, 1, 1))
boxplot(F_cont[ ,1], F_disc[ ,1], F_cont[ ,2], F_disc[ ,2], F_cont[ ,3], F_disc[ ,3], F_cont[ ,4], F_disc[ ,4],
        F_cont[ ,5], F_disc[ ,5], F_cont[ ,6], F_disc[ ,6], F_cont[ ,7], F_disc[ ,7],
        col = rep(c(grey(0.5), grey(0.8)), 7), 
        pch = 21, bg = rep(c(grey(0.5), grey(0.8)),7), xaxt = "n",cex.lab = 1.8, cex.axis = 1.5, xlab = "t")
mtext(bquote(hat(F)[1]), 2, line = 3, cex = 1.8)
axis(1, seq(1:7), at = seq(1.5, 14.5, by = 2), cex.axis = 1.5)
for(i in 1:7){
  lines(c(-1.5+(i*2),0.5+(i*2)),c(F_true[i], F_true[i]), col = "red", lwd = 1.5)
}
legend("topleft", legend = c("continuous", "discrete"), lty = "solid", col = c(grey(0.5), grey(0.8)), bty = "n", lwd = 3, cex = 1.5)

