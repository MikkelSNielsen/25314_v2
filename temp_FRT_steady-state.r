# this is a code for solving the system for FRT and steady state for different temperatures 
# the system is from the article https://www.sciencedirect.com/science/article/pii/S0304380023000534
# incl. a temperature dependence on the metabolic input 
# made by Mikkel Nielsen and Camilla Jensen ¨
# last modified 16-Mar-2023 

library(deSolve)
library(ggplot2)
library(dplyr)
library(ggtext)

param1 <- list(
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
  T = 273.15+26.11, # Temperature [K]
  T_ref = 273.15+26.11 # ref temperature [K]
)

# temperature dependence
param1$b = exp(-param1$E*(1/param1$k*(1/param1$T-1/param1$T_ref))) # Arrhenius equation 
param1$w_tilde = param1$w/param1$b  # Handling time dependent on temperature

param2 <- list(
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
param2$b = exp(-param2$E*(1/param2$k*(1/param2$T-1/param2$T_ref))) # Arrhenius equation 
param2$w_tilde = param2$w/param2$b  # Handling time dependent on temperature

param3 <- list(
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
  T = 273.15+43, # Temperature [K]
  T_ref = 273.15+26.11 # ref temperature [K]
)

# temperature dependence
param3$b = exp(-param3$E*(1/param3$k*(1/param3$T-1/param3$T_ref))) # Arrhenius equation 
param3$w_tilde = param3$w/param3$b  # Handling time dependent on temperature


# Initial conditions 
S0 <- 1.02*10**(-10) #[mol C]
H0 <- 6.50*10**(-6) #[mol C]
#S0 <- 0.08*10**(30)
#H0 <- 3*10**(60)
#S0 <- 6.50*10**(-6)
#H0 <- 1.02*10**(-10)

# Time 
t <- seq(0, 500, by = 1)

## number 1
# Dynamical equations 
HSmodel1 <- function(t, state, param1) {
  with(as.list(c(state, param1)), {
    dS1dt <- pmin(p*S1,
                  (a*H1*n_hn*S1*y_n)/(S1*n_sn+a*H1*n_sn*w_tilde)) - d_s*S1
    
    dH1dt <- -((a*H1*S1)/(S1+a*H1*w_tilde)) - d_h*H1 +
      pmin(((H1^(2/3)*r_n/n_hn)),
           H1^(2/3)*r_c+p*S1 - pmin(p*S1,
                                    (a*H1*n_hn*S1*y_n)/(S1*n_sn+a*H1*n_sn*w_tilde)))
    return(list(c(dS1dt, dH1dt)))
  })
}

# Solve the differential equations
res1 <- ode(y = c(S1 = S0, H1 = H0), times = t, func = HSmodel1, parms = param1)
S1 <- res1[, "S1"]
H1 <- res1[, "H1"]

# Plot results
# Create a data frame with the time and species abundances
df1 <- data.frame(time = t, S1 = S1, H1 = H1)

funcrespC1 <- (param1$a*H1*S1)/(S1+param1$a*H1*param1$w_tilde)

## number 2
# Dynamical equations 
HSmodel2 <- function(t, state, param2) {
  with(as.list(c(state, param2)), {
    dS2dt <- pmin(p*S2,
                  (a*H2*n_hn*S2*y_n)/(S2*n_sn+a*H2*n_sn*w_tilde)) - d_s*S2
    
    dH2dt <- -((a*H2*S2)/(S2+a*H2*w_tilde)) - d_h*H2 +
      pmin(((H2^(2/3)*r_n/n_hn)),
           H2^(2/3)*r_c+p*S2 - pmin(p*S2,
                                    (a*H2*n_hn*S2*y_n)/(S2*n_sn+a*H2*n_sn*w_tilde)))
    return(list(c(dS2dt, dH2dt)))
  })
}
# Solve the differential equations
res2 <- ode(y = c(S2 = S0, H2 = H0), times = t, func = HSmodel2, parms = param2)
S2 <- res2[, "S2"]
H2 <- res2[, "H2"]

# Plot results
# Create a data frame with the time and species abundances
df2 <- data.frame(time = t, S2 = S2, H2 = H2)

funcrespC2 <- (param2$a*H2*S2)/(S2+param2$a*H2*param2$w_tilde)

## number 3
# Dynamical equations 
HSmodel3 <- function(t, state, param3) {
  with(as.list(c(state, param3)), {
    dS3dt <- pmin(p*S3,
                  (a*H3*n_hn*S3*y_n)/(S3*n_sn+a*H3*n_sn*w_tilde)) - d_s*S3
    
    dH3dt <- -((a*H3*S3)/(S3+a*H3*w_tilde)) - d_h*H3 +
      pmin(((H3^(2/3)*r_n/n_hn)),
           H3^(2/3)*r_c+p*S3 - pmin(p*S3,
                                    (a*H3*n_hn*S3*y_n)/(S3*n_sn+a*H3*n_sn*w_tilde)))
    return(list(c(dS3dt, dH3dt)))
  })
}
# Solve the differential equations
res3 <- ode(y = c(S3 = S0, H3 = H0), times = t, func = HSmodel3, parms = param3)
S3 <- res3[, "S3"]
H3 <- res3[, "H3"]

# Plot results
# Create a data frame with the time and species abundances
df3 <- data.frame(time = t, S3 = S3, H3 = H3)

funcrespC3 <- (param3$a*H3*S3)/(S3+param3$a*H3*param3$w_tilde)


# # plot frt together for the three temps
 par(cex.lab=1.5, cex.axis=1.4)
 plot(H1,funcrespC1, col=c("#009999"), type="l", lwd=2, xlab = "Host biomass [C]", ylab = "FRT [C]")
 lines(H2,funcrespC2, col="#B266FF", type="l", lwd=2)
 lines(H3,funcrespC3, col="#FF99CC", type="l", lwd=2)
 legend("bottomright", legend=c("26.11\u00B0C","30.11\u00B0C","43\u00B0C"), col=c("#009999","#B266FF","#FF99CC"), lty=c(1), lwd=c(2), cex=1.6, bty="n",text.font = 2)
 title("FRT at different temperatures")
 
 par(cex.lab=1.5, cex.axis=1.4)
 plot(t,funcrespC1, col="#009999", type="l", lwd=2, xlab = "Time [weeks]", ylab = "FRT [C]")
 lines(t,funcrespC2, col="#B266FF", type="l", lwd=2)
 lines(t,funcrespC3, col="#FF99CC", type="l", lwd=2)
 legend("bottomright", legend=c("26.11\u00B0C","30.11\u00B0C","43\u00B0C"), col=c("#009999","#B266FF","#FF99CC"), lty=c(1), lwd=c(2), cex=1.6, bty="n",text.font = 2)
 title("FRT over time at different temperatures")

# plot with different initial values 
par(cex.lab=1.5, cex.axis=1.4)
plot(t,log(H1), col="#FF9933", type="l",lty=1, lwd=5, xlab = "Time [Weeks]", ylab = "ln(Biomass [C])", ylim = c(-16, -2), xlim = c(0, 500))
lines(t,log(S1), col="#00CC00", lty=1, lwd =5)
lines(t,log(H2), col="#CC6600", lty=2,lwd =5)
lines(t,log(S2), col="#33FF99", lty=2,lwd =5)
lines(t,log(H3), col="#663300", lty=3,lwd =5)
lines(t,log(S3), col="#006600", lty=3,lwd =5)
# Add a legend with a larger font size
legend("bottomright", legend=c("26.11\u00B0C: Host", "26.11\u00B0C: Symbiont", "30.11\u00B0C: Host", "30.11\u00B0C: Symbiont","43\u00B0C: Host", "43\u00B0C: Symbiont"), col=c("#FF9933", "#00CC00","#CC6600","#33FF99","#663300","#006600"), lty=c(1,1,2,2,3,3), lwd=c(2,2,2,2,2,2), cex=1.1, bty="n",text.font = 2)
title("Steady state solutions for different temperatures")

