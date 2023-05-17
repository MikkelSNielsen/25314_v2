# this is a code for solving the system for handling time and FRT 
# the system is from the article https://www.sciencedirect.com/science/article/pii/S0304380023000534
# incl. a temperature dependence on the metabolic input 
# made by Mikkel Nielsen and Camilla Jensen ¨
# last modified 16-Mar-2023 

library(deSolve)
library(ggplot2)
library(dplyr)
library(ggtext)

param <- list(
  p = 27.7, # Photosynthesis [week^-1]
  a = 0.0724, # Attack rate [week^-1]
  n_hn = 0.18, # Conversion rate of host to nitrogen [mol N mol C^-1]
  y_n = 0.0602, # Yield / conversion efficiency of host nitrogen to symbiont nitrogen [mol N mol N^-1]
  n_sn = 0.13, # Conversion of nitrogen to carbon [mol N mol C^-1]
  w = 0.00261, # Handling time of host for symbiont [week]
  d_s = 0.771, # Symbiont death rate [week^-1]
  d_h = 0.126, # Host death rate [week^-1]
  
  r_n = 0.0146, # Nitrogen uptake rate [mol N mol C^-1 week^-1] feed
  r_c = 0.00999, # Carbon uptake rate [week^-1] feed
  
  #r_n = 0, # starved 
  #r_c = 0, # starved
  
  E = 0.63*1.60217663*10**(-19), # Activation energy [J]
  k = 1.380659*10**(-23), # Boltzmann's constant [J/K]
  T = 273.15+30.11, # Temperature [K]
  T_ref = 273.15+26.11 # ref temperature [K]
)

# temperature dependence
param$b = exp(-param$E*(1/param$k*(1/param$T-1/param$T_ref))) # Arrhenius equation 
param$w_tilde = param$w/param$b  # Handling time dependent on temperature

# Initial conditions 
S0 <- 1.02*10**(-10) #[mol C]
H0 <- 6.50*10**(-6) #[mol C]
#S0 <- 0.08*10**(30)
#H0 <- 3*10**(60)
#S0 <- 6.50*10**(-6)
#H0 <- 1.02*10**(-10)

# Time 
t <- seq(0, 500, by = 1)

# test of f(S,H) when dependent on temperature
# obs H-end and S-end are for 26.11 degrees 
tempprofile <- seq(273.15+26.11, 273.15+43, length.out=length(t))
b_test <- vector("numeric", length(tempprofile))
w_test <- vector("numeric", length(tempprofile))
results <- vector("numeric", length(tempprofile))
H_end <- 0.0701689621
S_end <- 0.0005359714

for(i in 1:length(tempprofile)){
  
  b_test[i] = exp(-param$E*(1/param$k*(1/tempprofile[i]-1/param$T_ref)))
  w_test[i] = param$w/b_test[i]
  
  results[i] <- (param$a*H_end*S_end)/(S_end+param$a*H_end*w_test[i])
}

par(cex.lab=1.5, cex.axis=1.4)
plot(tempprofile-273.15, results, col="#FF6666", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Carbon")
legend("bottomright", legend=c("FRT"), col=c("#FF6666"), pch=21, lwd=c(2), cex=1.6, bty="n",text.font = 2)
title("FRT with increasing temperature")
par(cex.lab=1.5, cex.axis=1.4)
plot(tempprofile-273.15, w_test/b_test, col="magenta", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Week")
legend("topright", legend=c("Handling time w. temp.func"), col=c("magenta"), pch=21, lwd=c(2), cex=1.2, bty="n",text.font = 2)
title("Handling time w. temp.func., w. increasing temperature")


par(cex.lab=1.5, cex.axis=1.4)
plot(tempprofile-273.15, results, col="black", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Carbon")
legend("bottomright", legend=c("FRT"), col=c("black"), pch=21, lwd=c(2), cex=1.6, bty="n")
title("FRT with increasing temperature")



par(cex.lab=1.5, cex.axis=1.4)
plot(tempprofile-273.15, results, col="black", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Carbon")
legend("bottomright", legend=c("FRT"), col=c("black"), pch=21, lwd=c(2), cex=1.6, bty="n")
title("FRT with increasing temperature")

par(cex.lab=1.5, cex.axis=1.4)
plot(tempprofile-273.15, w_test/b_test, col="magenta", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Week")
legend("topright", legend=c("Handling time w. temp.func"), col=c("magenta"), pch=21, lwd=c(2), cex=1.6, bty="n")
title("Handling time w. temp.func., w. increasing temperature")

