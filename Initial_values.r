# this is a code for solving the system for different initital values in one plot - Done
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
S1 <- 0.08*10**(20) #much higher initial values
H1 <- 3*10**(50)
S2 <- 6.50*10**(-2) #S has a higher initial than H
H2 <- 1.02*10**(-4)

# Time 
t <- seq(0, 1000, by = 1)

# Dynamical equations - original initial values
HSmodel <- function(t, state, param) {
  with(as.list(c(state, param)), {
    dSdt <- pmin(p*S,
                 (a*H*n_hn*S*y_n)/(S*n_sn+a*H*n_sn*w_tilde)) - d_s*S
    
    dHdt <- -((a*H*S)/(S+a*H*w_tilde)) - d_h*H +
      pmin(((H^(2/3)*r_n/n_hn)),
           H^(2/3)*r_c+p*S - pmin(p*S,
                                  (a*H*n_hn*S*y_n)/(S*n_sn+a*H*n_sn*w_tilde)))
    return(list(c(dSdt, dHdt)))
  })
}


# Solve the differential equations
res <- ode(y = c(S = S0, H = H0), times = t, func = HSmodel, parms = param)
S <- res[, "S"]
H <- res[, "H"]

# Plot results
# Create a data frame with the time and species abundances
df <- data.frame(time = t, S = S, H = H)


# Solve the differential equations - 2 initial values
res1 <- ode(y = c(S = S1, H = H1), times = t, func = HSmodel, parms = param)
S1 <- res1[, "S"]
H1 <- res1[, "H"]

# Plot results
# Create a data frame with the time and species abundances
df1 <- data.frame(time = t, S1 = S1, H1 = H1)

# Solve the differential equations - 3 initial values 
res2 <- ode(y = c(S = S2, H = H2), times = t, func = HSmodel, parms = param)
S2 <- res2[, "S"]
H2 <- res2[, "H"]

# Plot results
# Create a data frame with the time and species abundances
df2 <- data.frame(time = t, S2 = S2, H2 = H2)


# plot with different initial values - log transformed
par(cex.lab=1.5, cex.axis=1.4)
plot(t,log(H), col="#FF9933", type="l", lwd=3, xlab = "Time [week]", ylab = "Biomass [C]", ylim = c(-16, 12), xlim = c(0, 1000))
lines(t,log(S), col="green",  lwd =3)
lines(t,log(H1), col="#994C00", lwd =3)
lines(t,log(S1), col="#009900", lwd =3)
lines(t,log(H2), col="#FFCC99", lwd =3)
lines(t,log(S2), col="#B2FF66", lwd =3)
# Add a legend with a larger font size
legend("topright", legend=c("IV1: Host", "IV1: Symbiont", "IV2: Host", "IV2: Symbiont","IV3: Host", "IV3: Symbiont"), col=c("orange", "green","#994C00","#009900","#FFCC99","#B2FF66"), lty=c(1,1), lwd=c(2,2), cex=0.8, bty="n",text.font = 2)
title("Steady state solutions for different initial values, temp=43\u00B0C")
