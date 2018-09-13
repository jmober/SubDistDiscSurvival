##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Simulations                     #####
#####                      presented in Section 3                        ##### 
##############################################################################

# The following contains R-Code to evaluate the results of the first part 
# of the simulation study. The first part of the simulation study compares 
# (i) weighted ML estimation, and (ii) unweighted ML estimation using the 
# complete censoring information. 

# The R-Code provided can be used to reprodce Table 1 and Figure 1 of the paper, 
# and Figure 2 and Figure 3 of the online appendix. 


# load R package and necessary files 
library("ggplot2")
load("res_100.rda")
load("res_300.rda")
load("res_500.rda")

# matrix containing the parameters of all scenarios 
args_org <- expand.grid(seed=c(1:1000) + 8617,
                    n=c(100, 300, 500), 
                    p=c(0.2, 0.4, 0.8), 
                    b=c(0.85, 1, 1.25))

# function to plot the estimated coefficients 

## Input
# data_list: list of results returned by one_sim() 
# args:      matrix of parameters 
# n:         number of observations 
# ylim:      y-limit 
# model:     type of model 

one_plot <- function(data_list, args, n, ylim, model = c("weighted", "unweighted")){
  
  beta <- c(0.4, -0.4, 0.2, -0.2)
  if(model == "weighted"){
    ind <- 1:4
  } 
  if(model == "unweighted"){
    ind <- 5:8 
  }
  
  args <- args[args$n == n, c(3, 4)]
  coef <- t(sapply(1:9000, function(j) data_list[[j]][[1]][ind, 1]))

  var      <- rep(c("x1", "x2", "x3", "x4"), each = 9000)
  coef     <- c(coef[, "x1"], coef[, "x2"], coef[, "x3"], coef[, "x4"])
  args     <- do.call(rbind, replicate(4, args, simplify = FALSE))
  for_plot <- data.frame(coef, var, args)
  
  ps    <- c("0.2" = "q=0.2", "0.4" = "q=0.4","0.8" = "q=0.8")
  cens  <- c("0.85" = "C: weak", "1" = "C: medium", "1.25" = "C: strong")
  
  ggplot(for_plot, aes(x = var,y = coef)) + 
    geom_boxplot(fill = "lightgray", outlier.color = "black", outlier.fill = "lightgray", outlier.shape = 21, fatten = 3) + 
    facet_grid(p~b, labeller = labeller(p = ps, b = cens)) + 
    geom_segment(aes(x = 0.45, y = beta[1], xend = 1.55, yend = beta[1]), colour = "red", size = 0.4) +
    geom_segment(aes(x = 1.45, y = beta[2], xend = 2.55, yend = beta[2]), colour = "red", size = 0.4) +
    geom_segment(aes(x = 2.45, y = beta[3], xend = 3.55, yend = beta[3]), colour = "red", size = 0.4) +
    geom_segment(aes(x = 3.45, y = beta[4], xend = 4.55, yend = beta[4]), colour = "red", size = 0.4) +
    ylim(ylim[1], ylim[2]) + 
    xlab("Covariate") +
    ylab(expression(hat(gamma))) +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
          strip.text.x = element_text(size= 15, face = "bold"),
          axis.title.y = element_text(angle = 0, vjust = 0.5, size = 15),
          axis.title.x = element_text(size = 15),
          axis.text = element_text(size = 14))
}

######
######

# Figure 1. The boxplots visualize the estimates of the parameters gamma1 that were obtained from fitting 
# a discrete subdistribution hazard model using the proposed weighted ML estimation approach (n = 500).
one_plot(data_list = res_500, args = args_org, n = 500, ylim = c(-2.1, 2.1), model = "weighted")

# Figure 2 (Online Appendix). The boxplots visualize the estimates of the parameters gamma1 that were obtained from fitting 
# a discrete subdistribution hazard model using the proposed weighted ML estimation approach (n = 100).
one_plot(data_list = res_100, args = args_org, n = 100, ylim = c(-2.1, 2.1), model = "weighted")

# Figure 3 (Online Appendix). The boxplots visualize the estimates of the parameters gamma1 that were obtained from fitting 
# a discrete subdistribution hazard model using the proposed weighted ML estimation approach (n = 300).
one_plot(data_list = res_300, args = args_org, n = 300, ylim = c(-2.1, 2.1), model = "weighted")


#####################################

# function to compute MSE(gamma), var(gamma) and E(var(gamma))

## Input 
# data_list: list of results returned by one_sim() 
# args:      matrix of parameters 
# n:         number of observations 
# p:         parameter affecting the probability of type 1 events 
# b:         parameter affecting the degree of censoring 
# model:     type of model 

one_table <- function(data_list, args, n, p, b, model = c("weighted", "unweighted")){
  
  beta <- c(0.4, -0.4, 0.2, -0.2)
  if(model == "weighted"){
    ind <- 1:4
  } 
  if(model == "unweighted"){
    ind <- 5:8 
  }

  args <- args[args$n == n, c(3,4)]
  ids  <- which(args$p == p & args$b == b)

  coef      <- t(sapply(ids, function(j) data_list[[j]][[1]][ind, 1]))
  ses       <- t(apply(coef, 1, function(x) (x-beta)^2))
  mses      <- apply(ses, 2, mean)
  vare      <- apply(coef, 2, var)
  varf      <- t(sapply(ids, function(j) data_list[[j]][[1]][ind, 2]^2))
  varf      <- apply(varf, 2, mean)
  
  to_return <- rbind(mses, vare, varf)
  colnames(to_return) <- paste0("x", 1:4)
  rownames(to_return) <- c("MSE(beta)", "Var(beta)", "E(Var)")
  
  return(to_return)
}

######
######

# Table 1. Upper part (weak) 
res_weak   <- lapply(c(0.2, 0.4, 0.8), function(j) one_table(res_500, args_org, 500, j, 0.85, "weighted"))
# Table 1. Middle part (medium)
res_medium <- lapply(c(0.2, 0.4, 0.8), function(j) one_table(res_500, args_org, 500, j, 1, "weighted"))
# Table 1. Lower part (strong)
res_strong <- lapply(c(0.2, 0.4, 0.8), function(j) one_table(res_500, args_org, 500, j, 1.25, "weighted"))


#####################################

# Figure 1 (Online Appendix). Illustration of the experimental design of the simulation study. 
# The bars display the average relative frequencies of observed events (0 = censoring event, 
# 1 = event of interest, 2 = competing event) that were obtained from 1000 simulated data sets (n = 500).

freq      <- t(sapply(1:9000, function(j) res_500[[j]][[2]]) / 500)
args      <- args_org[args_org$n == n, c(3, 4)]
for_means <- data.frame(freq, args)
all_means <- aggregate(cbind(X0, X1, X2) ~ p + b, data = for_means, FUN = "mean")

var      <- rep(c("X0", "X1", "X2"), each = 9)
coef     <- c(all_means[, "X0"], all_means[, "X1"], all_means[, "X2"])
args     <- do.call(rbind, replicate(3, all_means[, c(1, 2)], simplify = FALSE))
for_plot <- data.frame(coef, var, args)
  
ps    <- c("0.2" = "q=0.2", "0.4" = "q=0.4", "0.8" = "q=0.8")
cens  <- c("0.85" = "C: weak", "1" = "C: medium", "1.25" = "C: strong")
event <- c("X0" = "0", "X1" = "1", "X2" = "2")
  
# draw figure 
ggplot(for_plot, aes(x = var,y = coef)) + 
        geom_col(fill = rep(c("lightgray", "darkgray", "lightgray"), 9), colour = "black") + 
        facet_grid(p ~ b, labeller = labeller(p = ps, b = cens)) +
        scale_x_discrete(labels = event) +
        ylim(0, 1) + 
        xlab("Event") +
        ylab("Relative frequency") +
        theme_bw() +
        theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
              strip.text.x = element_text(size = 15, face = "bold"),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.text = element_text(size = 14))

