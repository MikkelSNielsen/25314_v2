# this is a code for solving the system for steady state, growth rate and growth limit - Done
# the system is from the article https://www.sciencedirect.com/science/article/pii/S0304380023000534
# incl. a temperature dependence on the metabolic input 
# made by Mikkel Nielsen and Camilla Jensen ¨
# last modified 16-Mar-2023 

library(deSolve)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggquiver)
library(fields)
library(ggquiver)
library(gridExtra)

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
  T = 273.15+43, # Temperature [K]
  T_ref = 273.15+26.11 # ref temperature [K]
)

# temperature dependence
param$b = exp(-param$E*(1/param$k*(1/param$T-1/param$T_ref))) # Arrhenius equation 
param$w_tilde = param$w/param$b  # Handling time dependent on temperature

# Initial conditions 
S0 <- 1.02*10**(-10) #[mol C]
H0 <- 6.50*10**(-6) #[mol C]

# Dynamical equations 
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

# Time 
t <- seq(0, 500, by = 1)

# Solve the differential equations
res <- ode(y = c(S = S0, H = H0), times = t, func = HSmodel, parms = param)
S <- res[, "S"]
H <- res[, "H"]

# Plot results
# Create a data frame with the time and species abundances
df <- data.frame(time = t, S = S, H = H)



# 
# par(mfrow=c(1,1))
# # Plot the dynamics of S and H over time
 ggplot(df, aes(x = time)) +
   geom_line(aes(y = log(S), color = "Symbiont"), size = 1.6) +
   geom_line(aes(y = log(H), color = "Host"), size = 1.6) +
   scale_color_manual(values = c("Symbiont" = "green", "Host" = "orange")) +
   labs(x = "Time [Weeks]", y = "Log - Biomass [Mol C]") +
   ggtitle("Steady state biomass - 43.11 \u00b0C") +
   theme(axis.text = element_text(size = 19),
         axis.title = element_text(size = 19),
         plot.title = element_text(hjust = 0.5, size = 23),
         legend.text = element_text(size = 16))
 
 #Steady state
 ggplot(df, aes(x = time)) +
   geom_line(aes(y = log(S), color = "Symbiont"), size = 1.6) +
   geom_line(aes(y = log(H), color = "Host"), size = 1.6) +
   scale_color_manual(values = c("Symbiont" = "green", "Host" = "orange")) +
   labs(x = "Time [Weeks]", y = "Log - Biomass [Mol C]") +
   ggtitle("Steady state biomass - 26.11 \u00b0C") +
   theme(axis.text = element_text(size = 19),
         axis.title = element_text(size = 19),
         plot.title = element_text(hjust = 0.5, size = 23),
         legend.text = element_text(size = 16))

 # Create plot with normal plot function of steady state
 plot(df$time, log(df$S), type = "l", col = "green",
      xlab = "Time [Weeks]", ylab = "Log - Biomass [Mol C]",
      main = "Steady state biomass - 26.11 \u00b0C", lwd = 4,
      ylim = c(-25, -3))
 lines(df$time, log(df$H), col = "orange", lwd = 4)
 legend("bottomright", legend = c("Symbiont", "Host"),
        col = c("green", "orange"), lty = 1, lwd = 4, cex = 2, bty = "n")
 
 


   #plot of the growth limitation
    par(cex.lab=1.3, cex.axis=1.2)
    plot(param$p*S, col="black", type="l", log="y", lwd=3, xlab = "Time [weeks]", ylab = "Biomass [mol C]", main = "Growth limitation - 43 \u00b0C")
    # Add the blue line
    lines((param$a*H*param$n_hn*S*param$y_n)/(S*param$n_sn+param$a*H*param$n_sn*param$w_tilde), col="blue", lwd =3)
    # Add a legend with a larger font size
    legend("bottomright", legend=c("Host", "Symbiont"), col=c("black", "blue"), lty=c(1,1), lwd=c(2,2), cex=1.6, bty="n")
    
    par(cex.lab = 1.3, cex.axis = 1.2)
    
    par(cex.lab = 1.3, cex.axis = 1.2)
  
   
    
   # Calculate the concentration of the limiting nutrient
   # Plot the biomass of the host
plot(param$p * S, col = "black", type = "l", log = "y", lwd = 3,
     xlab = "Time [weeks]", ylab = "Biomass [mol C]", main = "Symbiont growth limitation - 43 \u00b0C",
     cex.axis = 1.2, cex.lab = 1.3)

# Add the biomass of the symbiont
lines((param$a * H * param$n_hn * S * param$y_n) / (S * param$n_sn + param$a * H * param$n_sn * param$w_tilde),
      col = "blue", lwd = 3)

# Add a legend with larger font size
legend("bottomright", legend = c("Carbon", "Nitrogen"), col = c("black", "blue"),
       lty = c(1, 1), lwd = c(2, 2), cex = 1.6, bty = "n")

  

#Plot the limiting growth factor for the host
   plot((H^(2/3)*param$r_n/param$n_hn), col = "blue", type = "l", lwd = 3,
        xlab = "Time [Weeks]", ylab = "Biomass [mol C]", main = "Host growth limitation - 26.11 \u00b0C",
        ylim = c(0, 0.016), cex.axis = 1.2, cex.lab = 1.3)
   
   # Add the expression for carbon limitation
   lines(H^(2/3)*param$r_c + param$p*S - pmin(param$p*S, (param$a*H*param$n_hn*S*param$y_n) / (S*param$n_sn + param$a*H*param$n_sn*param$w_tilde)),
         col = "black", lwd = 3)
   
   # Find the time when the lines cross
   cross_time <- which(diff(sign(H^(2/3)*param$r_n/param$n_hn - H^(2/3)*param$r_c - param$p*S +
                                   (param$a*H*param$n_hn*S*param$y_n) / (S*param$n_sn + param$a*H*param$n_sn*param$w_tilde))) != 0)
   
   # Add a red dotted line at the crossing point
   abline(v = t[cross_time], lty = 3, col = "red", lwd = 3)
   
   # Add a legend with the dotted line
   legend("bottomright", legend = c("Nitrogen", "Carbon", "Crossing Point"), col = c("blue", "black", "red"),
          lty = c(1, 1, 3), lwd = c(2, 2, 3), cex = 1.6, bty = "n")
   
   


   # test of f(S,H) when dependent on temperature
   # obs H-end and S-end are for 26.11 degrees 
   tempprofile <- seq(273.15+26.11, 273.15+243, length.out=length(t))
   b_test <- vector("numeric", length(tempprofile))
   w_test <- vector("numeric", length(tempprofile))
   results <- vector("numeric", length(tempprofile))
   H_end <- 0.0701689621*10**3 
   S_end <- 0.0005359714
   
   for(i in 1:length(tempprofile)){
     
     b_test[i] = exp(-param$E*(1/param$k*(1/tempprofile[i]-1/param$T_ref)))
     w_test[i] = param$w/b_test[i]
     print(param$a*H_end*S_end/S_end)
     results[i] <- (param$a*H_end*S_end)/(S_end+param$a*H_end*w_test[i])
   }
   
   par(cex.lab=1.5, cex.axis=1.4)
   plot(tempprofile-273.15, results, col="black", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Carbon")
   legend("bottomright", legend=c("FRT"), col=c("black"), pch=21, lwd=c(2), cex=1.6, bty="n")
   title("FRT with increasing temperature")
   
   par(cex.lab=1.5, cex.axis=1.4)
   plot(tempprofile-273.15, w_test/b_test, col="magenta", pch=21, lwd=3, xlab = "Temperature [\u00B0C]", ylab = "Week")
   legend("topright", legend=c("Handling time w. temp.func"), col=c("magenta"), pch=21, lwd=c(2), cex=1.6, bty="n")
   title("Handling time w. temp.func., w. increasing temperature")
   


   # plot of the functional response incl. n and y 
   #funcresp <- (param$a*H*param$n_hn*S*param$y_n)/(S*param$n_sn+param$a*H*param$n_sn*param$w_tilde)
   funcrespC <- (param$a*H*S)/(S+param$a*H*param$w_tilde)
   
   
   par(cex.lab=1.5, cex.axis=1.4)
   plot(H,funcrespC, col="magenta", type="l", lwd=3, xlab = "Host biomass [C]", ylab = "FRT [C]")
   legend("bottomright", legend=c("Growth factor"), col=c("magenta"), lty=c(1), lwd=c(2), cex=1.6, bty="n")
   title("Growth of the host with FRT")
   
   par(cex.lab=1.5, cex.axis=1.4)
   plot(S,funcrespC, col="magenta", type="l", lwd=3, xlab = "Symbiont biomass [C]", ylab = "FRT [C]")
   legend("bottomright", legend=c("Growth factor"), col=c("magenta"), lty=c(1), lwd=c(2), cex=1.6, bty="n")
   title("Growth of the symbiont with FRT")
   
   par(cex.lab=1.5, cex.axis=1.4)
   plot(t,funcrespC, col="blue", type="l", lwd=3, xlab = "Time [weeks]", ylab = "Carbon")
   legend("bottomright", legend=c("FRT"), col=c("blue"), lty=c(1), lwd=c(2), cex=1.6, bty="n")
   title("FRT with time")
   
   #plot with different initial values 
   par(cex.lab=1.5, cex.axis=1.4)
   plot(t,H, col="orange", log = "y", type="l", lwd=3, xlab = "Host biomass [C]", ylab = "FRT [N]")
   lines(t,S, col="green", log = "y", lwd =3)
   # Add a legend with a larger font size
   legend("bottomright", legend=c("Host", "Symbiont"), col=c("orange", "green"), lty=c(1,1), lwd=c(2,2), cex=1.6, bty="n")
   title("Growth of the host with FRT")
  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    