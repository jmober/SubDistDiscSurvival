##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Application                     #####
#####                      presented in Section 4                        ##### 
##############################################################################


# load data 
load("dataS.rda")
# The dataset contains a subset of n=1600 individuals of the original data
# used for the analysis presented in the paper, where the values of the 
# explanatory variable had been slightly modified for 10% of the individuals. 

# load R packages 
library("discSurv")
library("mgcv")

# prepare data 
dataS     <- dataS[,-9]    # omit LRT
dataS$age <- dataS$age-60  # center age 

# build augmented data matrix 
dataLongS <- dataLongSubDist(dataSet=dataS, timeColumn="day", eventColumns=c("e1", "e2"), eventFocus="e1")

# fit model with discrete baseline coefficients 
form1 <- as.formula(paste("y~timeInt-1+", paste(names(dataS)[c(4:16)], collapse="+"))) # formula 
mod1  <- glm(form1, family=binomial(link="cloglog"), data = dataLongS, weights=dataLongS$subDistWeights)
# call glm() with the subDistWeights computed by dataLongSubDist() 
# passed to the weights-argument 

# results analogous to the results presented in Table 2, model 1 (discrete bc)
tail(summary(mod1)$coefficients,13)[,c(1,2,4)]

# fit model with smooth baseline coefficients 
form2 <- as.formula(paste("y~s(timeInt, bs='ps', m=2) +",paste(names(dataS)[c(4:16)], collapse="+"))) # formula 
dataLongS$timeInt <- as.numeric(dataLongS$timeInt) # transform time variable 
mod2 <- gam(form2, family=binomial(link="cloglog"), data = dataLongS, weights=dataLongS$subDistWeights)
# call mgcv() with the subDistWeights computed by dataLongSubDist() 
# passed to the weights-argument 

# results analogous to the results presented in Table 2, model 2 (smooth bc)
tail(summary(mod2)$p.table,13)[,c(1,2,4)]

#####
##### 

# Figure 13 (Online Appendix). The figure shows the estimated baseline coefficients
# of model 1 (black dots) and the estimated smooth baseline function of model 2 (gray dots)

# calculate and transform baseline coefficients 
int1 <- coef(mod1)[1:60]
bh1  <- 1-exp(-exp(int1))
obs  <- sapply(1:60, function(i) which(dataLongS$timeInt == i)[1])
int2 <- predict(mod2, type = "terms")[obs, 14] + coef(mod2)[1]
bh2  <- 1-exp(-exp(int2))

# draw figure 
par(mar=c(4.5, 6.5, 1, 1))
plot(x = 1:60, y = bh1, pch = 19, 
     xlab = "Days",
     ylab = expression(paste("1-exp(-exp(", hat(gamma)[ot], "))")),
     cex = 1.8, cex.lab = 1.8, cex.axis = 1.5)
lines(x = 1:60, y = bh2, type = "b", lwd = 5, col = "grey", pch=19)

#####
#####

# Figure 3a. The figure shows the estimated cumulative incidence functions for NP acquisition 
# referring to the covariate profile of a randomly selected study participant. It refers to 
# as situation, where the risk factor elective surgery before admission would have been present 
# in this patient. 
vec1    <- vec2 <- vec3 <- dataS[62,4:16]
vec2[9] <- 1
eta_hat <- c(as.matrix(vec1)%*%coef(mod1)[61:73], as.matrix(vec2)%*%coef(mod1)[61:73])
eta_hat <- matrix(rep(mod1$coef[1:60],each=2),nrow=2)+eta_hat
lambda  <- 1-exp(-exp(eta_hat))
Shat    <- t(apply(lambda,1,function(x) cumprod(1-x)))
Fhat    <- 1-Shat

# draw figure 
par(mar=c(4.5,6.5,1,1))
plot(c(0:60),c(0,Fhat[1,]), type="s",lwd=5, ylim=c(0,0.08), col="grey", xlab="Days",
     cex.lab=1.8, cex.axis=1.5, ylab=bquote(hat(F)[1]))
lines(c(0:60),c(0,Fhat[2,]),type="s",lwd=5)
legend("topleft", c("with elective surgery","without elective surgery"), lty="solid", col=c("black","grey"), bty="n",
       lwd=3, cex=1.5)

# Figure 3b. The figure shows the estimated cumulative incidence functions for NP acquisition 
# referring to the covariate profile of a randomly selected study participant. It refers to 
# as situation, where the risk factor intubation on admission would have been present 
# in this patient. 
vec3[4] <- 1
eta_hat <- c(as.matrix(vec1)%*%coef(mod1)[61:73], as.matrix(vec3)%*%coef(mod1)[61:73])
eta_hat <- matrix(rep(mod1$coef[1:60],each=2),nrow=2)+eta_hat
lambda  <- 1-exp(-exp(eta_hat))
Shat    <- t(apply(lambda,1,function(x) cumprod(1-x)))
Fhat    <- 1-Shat

# draw figure 
par(mar=c(4.5,6.5,1,1))
plot(c(0:60),c(0,Fhat[1,]), type="s",lwd=5, ylim=c(0,0.08), col="grey", xlab="Days",
     cex.lab=1.8, cex.axis=1.5, ylab=bquote(hat(F)[1]))
lines(c(0:60),c(0,Fhat[2,]),type="s",lwd=5)
legend("topleft", c("with intubation","without intubation"), lty="solid", col=c("black","grey"), bty="n",
       lwd=3, cex=1.5)

#####
#####

# Figure 12. Check of the proportional subdistribution hazards assumption. The figure shows 
# the estimated cumulative incidence functions for NP acquisition that were obtained from fitting 
# covariate-free discrete subdistribution hazard models to subsets of the data. Solid lines refer 
# to subgroups defined by intubation on admission and/or elective surgery at admission. Dashed lines 
# refer to the respective average estimated cumulative incidence functions obtained from the 
# complementary log-log model with all covariates.

# function to fit covariate-free models 
CIFnon <- function(intub, elect){
  dataS_s     <- dataS[dataS$intubation == intub & dataS$elective == elect,]
  dataLongS_s <- dataLongSubDist(dataSet=dataS_s, timeColumn="day", eventColumns=c("e1", "e2"), eventFocus="e1")
  mod         <- glm(y~timeInt-1, family=binomial(link="cloglog"), data = dataLongS_s, weights=dataLongS_s$subDistWeights)
  coef        <- coef(mod)
  lambda      <- 1-exp(-exp(coef))
  S           <- cumprod(1-lambda)
  CIF         <- 1-S
}

# function to compute average estimated F 
CIFsub <- function(intub, elect){
  bh     <- coef(mod1)[1:60]
  coefs  <- coef(mod1)[61:73]
  X_s    <- dataS[dataS$intubation == intub & dataS$elective == elect,]
  eta    <- as.matrix(X_s[,4:16])%*%coefs
  etas   <- t(sapply(1:length(eta), function(j) eta[j]+bh))
  lambda <- 1-exp(-exp(etas))
  S      <- t(apply(lambda,1,function(x) cumprod(1-x)))
  CIF    <- 1-S
  CIFmean  <- apply(CIF,2,mean, na.rm=TRUE)
}

# fit covariate-free models 
CIF00 <- CIFnon(0,0)
CIF10 <- CIFnon(1,0)
CIF11 <- CIFnon(1,1)

# draw figure 
par(mar=c(4.5,4.5,1,1))
plot(0:60, c(0,CIF00), type="s", ylim=c(0,0.23),col = grey(0.8), xlab = "Days", ylab = "Cumulative incidence function", cex.lab=1.5, cex.axis=1.5, lwd = 2)
legend("topleft", legend=c( "with intubation & with elective surgery", "with intubation & without elective surgery", "without intubation & without elective surgery"),
       col=c("black",grey(0.5),grey(0.8)),bty="n", lwd=3)

lines(0:60, c(0,CIF10), type="s",col = grey(0.5), lwd = 2)
lines(0:60, c(0,CIF11), type="s",col = "black", lwd = 2)
lines(0:60, c(0,CIFsub(0,0)), lty="dashed", type="s",col = grey(0.8), lwd = 2)
lines(0:60, c(0,CIFsub(1,0)), lty="dashed", type="s", col = grey(0.5), lwd = 2)
lines(0:60, c(0,CIFsub(1,1)), lty="dashed", type="s",col = "black", lwd = 2)




