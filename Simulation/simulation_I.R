##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Simulations                     #####
#####                      presented in Section 3                        ##### 
##############################################################################

# The following contains R-Code to conduct the first part of the simulation study. 
# The first part of the simulation study compares (i) weighted ML estimation, 
# and (ii) unweighted ML estimation using the complete censoring information. 


### Function one_sim() ### 

## Input 
# seed:   integer passed to set.seed() 
# betas1: true model parameters for type 1 events (of interest)
# betas2: true model parameters for type 2 events (competing)
# n:      number of observations 
# p:      parameter affecting the probability of type 1 events 
# b:      parameter affecting the degree of censoring 

## Description
# The function 
# - generates one data set according to the design described in 
#   Section 3.1 of the paper. 
# - fits the model based on weighted ML estimation and unweighted ML estimation 
#   using the complete censoring information 
# 
# The function requires the R add-on packages "discSurv" and the functions 
# implemented in "functions_I.R".

## Output 
# List with two elements: 
# - results: Matrix with 8 rows and 4 columns. Each row contains the 
#            parameter estimate, standard error, z value and p-value. 
#            Row 1-4 refer to weighted ML estimation, row 5-8 refer to 
#            unweighted ML estimation 
# - tab:     Table of the type of events in the data  

one_sim <- function(seed, betas1, betas2, n, p, b){
  
  # load required files 
  source("functions_I.R")
  library(discSurv)

  # compute quantiles 
  limits <- get_limits(p)
  
  # draw X
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rbinom(n, 1, 0.5)
  X  <- cbind(x1, x2, x3, x4)
  
  # generate data set 
  valid <- FALSE
  while(!valid){
    data_disc <- all_disc_data(X, betas1, betas2, p, b, limits)
    if(length(table(data_disc[, "status"])) == 3){
      valid <- TRUE
    }
  }
  
  # fit model by weighted ML estimation 
  data_disc_long_sub <- dataLongSubDist(dataSet = data_disc, timeColumn = "time", eventColumns = c("e1", "e2"), eventFocus = "e1")
  model_sub          <- glm(y ~ timeInt - 1 + x1 + x2 + x3 + x4, 
                            family = binomial(link = "cloglog"), 
                            data = data_disc_long_sub, 
                            weights = data_disc_long_sub$subDistWeights)
  # call glm() with the subDistWeights computed by dataLongSubDist() 
  # passed to the weights-argument
  
  # fit model by unweighted ML estimation 
  data_disc$time[data_disc$e2 == 1] <- data_disc$C[data_disc$e2 == 1] # use complete censoring information 
  data_disc_long  <- dataLong(dataSet = data_disc, timeColumn = "time", censColumn = "e1")
  model_censknown <- glm(y ~ timeInt - 1 + x1 + x2 + x3 + x4, 
                         family = binomial(link = "cloglog"), 
                         data = data_disc_long) 
  
  # build list with results
  results   <- rbind(tail(summary(model_sub)$coefficients, 4), tail(summary(model_censknown)$coefficients, 4))
  to_return <- list(results = results, tab = table(data_disc$status))
  
  return(to_return)
}

###########################

### Execution ### 

# The results of the first part of the simulation are stored in "res_100.rda", "res_300.rda" and "res_500.rda". 
# The function one_sim() was called 27,000 with fixed values for betas1 and betas2 and the following values 
# for seed, n, p, and b:  
args <- expand.grid(seed = c(1:1000) + 8617,
                    n = c(100, 300, 500), 
                    p = c(0.2, 0.4, 0.8), 
                    b = c(0.85, 1, 1.25))

# true coefficients 
b1 <- c(0.4, -0.4, 0.2, -0.2)
b2 <- c(-0.4, 0.4, -0.2, 0.2)

## Example 1: seed=8618, n=100, p=0.2, b=0.85
args[1, ]
one_sim(seed = 8618, betas1 = b1, betas2 = b2, n = 100, p = 0.2, b = 0.85)

# result is stored in 
load("res_100.rda")
res_100[[1]]


## Example 2: seed=8620, n=500, p=0.4, b=1 
args[14003, ]
one_sim(seed = 8620, betas1 = b1, betas2 = b2, n = 500, p = 0.4, b = 1)

# result is stored in 
load("res_500.rda")
res_500[[4003]]

###########################

# The whole simulation was executed on a Linux Cluster 
# by use of the R add-on package "batchtools". 
# The original call was the following: 

# library("batchtools")
# r <- makeRegistry(file.dir="./Simulation/simulation_I")
# r <- loadRegistry(file.dir="./Simulation/simulation_I")
# 
# args <- expand.grid(seed = c(1:1000) + 8617,
#                     n = c(100, 300, 500), 
#                     p = c(0.2, 0.4, 0.8), 
#                     b = c(0.85, 1, 1.25))
# 
# batchMap(one_sim, args = args, more.args = list(betas1=c(0.4, -0.4, 0.2, -0.2), betas2=c(-0.4, 0.4, -0.2, 0.2)), reg = r)
# ch  = chunkIds(ids = 1:27000, n.chunks = 100, reg = r)
# submitJobs(ch, reg = r)



