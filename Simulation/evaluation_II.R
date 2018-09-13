##############################################################################
##### Subdistribution Hazard Models for Competing Risks in Discrete Time #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################
#####             Content: R-Code of the Simulations                     #####
#####                      presented in Section 3                        ##### 
##############################################################################

# The following contains R-Code to evaluate the results of the second part of the 
# simulation study. The second part of the simulation study compares (i) the 
# proposed discrete subdistribution model, and (ii) the continuous Fine & Gray  
# model.  

# The R-Code provided can be used to reproduce Figures 4 to 7 and Figure 11 of 
# the online appendix. 

# load R package and necessary files 
library("ggplot2")
load("res_II.rda")

# matrix containing the parameters of all scenarios 
args_org <- expand.grid(seed=c(1:1000) + 8617,
                        n_limits=c(4, 8, 16, 32, 64),
                        p=c(0.2, 0.4, 0.8), 
                        b=c(0.85, 1, 1.25))
args_org$n_limits <- as.factor(args_org$n_limits)

# labels 
ps    <- c("0.2" = "q=0.2","0.4" = "q=0.4","0.8" = "q=0.8")
cens  <- c("0.85" = "C: weak", "1" = "C: medium", "1.25" = "C: strong")

#####
#####

# Figure 5 (Online Appendix). The boxplots visualize the differences between the discrete time
# subdistribution hazard model and the continuous-time Fine & Gray model in the estimates
# of the coefficient gamma_11 = 0.4 (n = 500).

est_x1   <- sapply(1:45000, function(j) res[[j]][[1]][1,])
diff_x1  <- apply(est_x1, 2, diff)
for_plot <- data.frame(diff_x1, args_org)

# draw figure 
ggplot(for_plot, aes(x = n_limits, y = diff_x1)) +
       geom_hline(aes(yintercept = 0), colour = "red", size = 0.4) +
       geom_boxplot(fill = "lightgray", outlier.color = "black", outlier.fill = "lightgray", outlier.shape = 21, fatten = 3) + 
       facet_grid(p~b, labeller = labeller(p = ps, b = cens)) + 
       ylab(bquote(hat(gamma)['11,cont']-hat(gamma)['11,disc'])) +
       xlab("k") +
       theme_bw() +
       theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
             strip.text.x = element_text(size = 15, face = "bold"),
             axis.title.y = element_text(size = 15),
             axis.title.x = element_text(size = 15),
             axis.text = element_text(size = 14))

#####
#####

# Figure 7 (Online Appendix). The boxplots visualize the differences between the discrete time
# subdistribution hazard model and the continuous-time Fine & Gray model in the estimates
# of the coefficient gamma_12 = - 0.4 (n = 500).

est_x2  <- sapply(1:45000, function(j) res[[j]][[1]][2, ])
diff_x2 <- apply(est_x2, 2, diff)
for_plot <- data.frame(diff_x2, args_org)

# draw figure 
ggplot(for_plot, aes(x = n_limits, y = diff_x2)) +
  geom_hline(aes(yintercept = 0), colour ="red", size = 0.4) +
  geom_boxplot(fill="lightgray", outlier.color = "black", outlier.fill = "lightgray", outlier.shape = 21, fatten = 3) + 
  facet_grid(p~b, labeller = labeller(p = ps, b = cens)) + 
  ylab(bquote(hat(gamma)['12,cont']-hat(gamma)['12,disc'])) +
  xlab("k") +
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 14))

#####
#####

# Figure 4 (Online Appendix). The boxplots visualize the estimates of the coefficient 
# gamma_11 = 0.4, as obtained from the discrete-time subdistribution hazard model and 
# the continuous time Fine & Gray model (n = 500).

args     <- do.call(rbind, replicate(2, args_org, simplify = FALSE))
est_x1   <- sapply(1:45000, function(j) res[[j]][[1]][1, ])
est_x1   <- c(est_x1[1, ], est_x1[2, ])

Method   <- rep(c("discrete","continuous"), each = 45000)
for_plot <- data.frame(est_x1, args, Method)

# draw figure 
ggplot(for_plot, aes(x = n_limits, y = est_x1))+
  scale_fill_manual(values = c(grey(0.5), grey(0.8))) +
  geom_hline(aes(yintercept = 0.4), colour = "red", size = 0.4) +
  facet_grid(p~b, labeller = labeller(p = ps, b = cens)) + 
  ylab(bquote(hat(gamma)[11])) +
  ylim(-0.3, 1.1) +
  xlab("k") +
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
        strip.text.x = element_text(size= 15, face = "bold"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.position = "top",
        legend.justification = "left")+
  geom_boxplot(aes(colour = Method))+
  geom_boxplot(aes(fill = Method), outlier.shape = 21, fatten = 3)


#####
#####

# Figure 6 (Online Appendix). The boxplots visualize the estimates of the coefficient 
# gamma_12 = -0.4, as obtained from the discrete-time subdistribution hazard model and 
# the continuous time Fine & Gray model (n = 500).

args     <- do.call(rbind, replicate(2, args_org, simplify = FALSE))
est_x2   <- sapply(1:45000, function(j) res[[j]][[1]][2, ])
est_x2   <- c(est_x2[1, ], est_x2[2, ])

Method <- rep(c("discrete","continuous"), each = 45000)
for_plot <- data.frame(est_x2, args, Method)

ggplot(for_plot, aes(x = n_limits, y = est_x2))+
  scale_fill_manual(values = c(grey(0.5), grey(0.8)))+
  geom_hline(aes(yintercept = -0.4), colour = "red", size = 0.4)+
  facet_grid(p~b, labeller=labeller(p = ps, b = cens)) + 
  ylab(bquote(hat(gamma)[12]))+
  ylim(-1.1, 0.3)+
  xlab("k")+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.position = "top",
        legend.justification = "left")+
  geom_boxplot(aes(colour = Method))+
  geom_boxplot(aes(fill = Method), outlier.shape = 21, fatten = 3)

#####
#####

# Figure 11 (Online Appendix). The boxplots illustrate the run-times that were obtained
# from fitting discrete-time subdistribution hazard models and continuous-time Fine & Gray models
# to the 1000 data sets. The following R-Code draws a panel for all 9 scenarios. 

args <- do.call(rbind, replicate(2, args_org, simplify = FALSE))
time <- sapply(1:45000, function(j) res[[j]][[2]])
time <- c(time[1, ], time[2, ])

Method <- rep(c("discrete","continuous"), each = 45000)
for_plot <- data.frame(time, args, Method)

# draw figure 
ggplot(for_plot, aes(x = n_limits,y = time)) +
  scale_fill_manual(values = c(grey(0.5), grey(0.8))) +
  ylab("Run-time (s)") +
  xlab("k") +
  geom_boxplot(aes(colour = Method)) +
  geom_boxplot(aes(fill = Method), outlier.shape = 21, fatten = 3) +
  facet_grid(p~b, labeller = labeller(p = ps, b = cens)) + 
  ylim(0, 7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.position = "top",
        legend.justification = "left",
        panel.spacing.x = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold"))










