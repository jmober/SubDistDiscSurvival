##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Simulations                     #####
#####                      presented in Section 3                        ##### 
##############################################################################

# The following contains R-Code to conduct the second part of the 
# simulation study. The second part of the simulation study compares (i) the 
# proposed discrete subdistribution model, and (ii) the continuous Fine & Gray  
# model.  


### Function one_sim() ### 

## Input 
# seed:     integer passed to set.seed() 
# betas1:   true model parameters for type 1 events (of interest)
# betas2:   true model parameters for type 2 events (competing)
# n_limits: number of intervals (k)
# p:        parameter affecting the probability of type 1 events 
# b:        parameter affecting the degree of censoring 

## Description
# The function 
# - generates one data set according to the design described in 
#   Section 3.1 of the paper. 
# - fits the proposed discrete subdistribution model and the continuous 
#   Fine & Gray model. 
# 
# The function requires the R add-on packages "discSurv" and "cmprsk" and the functions 
# implemented in "functions_II.R".

## Output 
# List with two elements: 
# - est:  Matrix with 4 rows and 2 columns. Each row contains the 
#         parameter estimate for one covariate. The first column 
#         refers to the discrete subdistribution model, the second 
#         column to the continuous Fine & Gray model. 
# - time: Run-times of the two methods (discrete, continuous)  

one_sim <- function(seed, betas1, betas2, n_limits, p, b){
  
  # load required files 
  source("functions_II.R")
  library(discSurv)
  library(cmprsk)
  
  # compute quantiles
  limits <- get_limits(p, n_limits)
  
  # draw X 
  n <- 500
  set.seed(seed)
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

  # fit discrete subdistribution model 
  data_disc_long_sub <- dataLongSubDist(dataSet = data_disc, timeColumn = "time", eventColumns = c("e1","e2"), eventFocus = "e1")
  ptm <- proc.time()
  model_sub          <- glm(y ~ timeInt - 1 + x1 + x2 + x3 + x4, 
                            family = binomial(link="cloglog"), 
                            data = data_disc_long_sub, 
                            weights = data_disc_long_sub$subDistWeights)
  # call glm() with the subDistWeights computed by dataLongSubDist() 
  # passed to the weights-argument
  time_sub <- proc.time()-ptm # calculate run-time 
  est_sub  <- tail(coef(model_sub), 4)

  # fit continuous Fine & Gray model 
  ptm <- proc.time()
  model_cont <- crr(data_disc$time, data_disc$status, X, failcode = 1, cencode = 0) # call crr()
  time_cont <- proc.time()-ptm # calculate run-time 
  est_cont  <- summary(model_cont)$coef[ , 1]
  
  # build list with results
  est  <- cbind(est_sub, est_cont)
  time <- c(time_sub[3], time_cont[3])
  to_return <- list(est = est, time = time)
  
  return(to_return)
}

###########################

### Execution ### 

# The results of the second part of the simulation are stored in "res_II.rda". 
# The function one_sim() was called 45,000 with fixed values for betas1 and betas2 and the following values 
# for seed, n_limits, p, and b:  
args <- expand.grid(seed = c(1:1000) + 8617,
                    n_limits = c(4, 8, 16, 32, 64),
                    p = c(0.2, 0.4, 0.8), 
                    b = c(0.85, 1, 1.25))

# Note: Due to the use of different computing systems, the resulting run-times for the two models 
#       might differ to the ones reported in the paper! 

# true coefficients 
b1 <- c(0.4, -0.4, 0.2, -0.2)
b2 <- c(-0.4, 0.4, -0.2, 0.2)

## Example 1: seed=8618, n_limits=4, p=0.2, b=0.85
args[1, ]
one_sim(seed = 8618, betas1 = b1, betas2 = b2, n_limits = 4, p = 0.2, b = 0.85)

# result is stored in 
load("res_II.rda")
res[[1]]


## Example 2: seed=8620, n_limits=8, p=0.4, b=1 
args[21003, ]
one_sim(seed = 8620, betas1 = b1, betas2 = b2, n_limits = 8, p = 0.4, b = 1)

# result is stored in 
load("res_II.rda")
res[[21003]]

###########################

# The whole simulation was executed on a Linux Cluster 
# by use of the R add-on package "batchtools". 
# The original call was the following: 


# library("batchtools")
# r <- makeRegistry(file.dir="./Simulation/simulation_II")
# r <- loadRegistry(file.dir="./Simulation/simulation_II")
# 
# args <- expand.grid(seed = c(1:1000) + 8617,
#                     n_limits = c(4, 8, 16, 32, 64),
#                     p = c(0.2, 0.4, 0.8), 
#                     b = c(0.85, 1, 1.25))
# 
# batchMap(one_sim, args=args, more.args = list(betas1=c(0.4, -0.4, 0.2, -0.2), betas2=c(-0.4, 0.4, -0.2, 0.2)), reg=r)
# ch  = chunkIds(ids=1:45000, n.chunks=450, reg=r)
# submitJobs(ch, reg = r)


