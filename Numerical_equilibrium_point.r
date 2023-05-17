# this is a code for calculating the numerical equilibrium points for the system - Done
# the system is from the article https://www.sciencedirect.com/science/article/pii/S0304380023000534
# incl. a temperature dependence on the metabolic input 
# made by Mikkel Nielsen and Camilla Jensen ¨
# last modified 16-Mar-2023 

p = 27.7 # Photosynthesis [week^-1]
a = 0.0724 # Attack rate [week^-1]
n_hn = 0.18 # Conversion rate of host to nitrogen [mol N mol C^-1]
y_n = 0.0602 # Yield / conversion efficiency of host nitrogen to symbiont nitrogen [mol N mol N^-1]
n_sn = 0.13 # Conversion of nitrogen to carbon [mol N mol C^-1]
w = 0.00261 # Handling time of host for symbiont [week]
d_s = 0.771 # Symbiont death rate [week^-1]
d_h = 0.126 # Host death rate [week^-1]

r_n = 0.0146 # Nitrogen uptake rate [mol N mol C^-1 week^-1] feed
r_c = 0.00999 # Carbon uptake rate [week^-1] feed

E = 0.63*1.60217663*10**(-19) # Activation energy [J]
k = 1.380659*10**(-23) # Boltzmann's constant [J/K]
T = 273.15+43 # Temperature [K]
T_ref = 273.15+26.11 # ref temperature [K]

b = exp(-E*(1/k*(1/T-1/T_ref))) # Arrhenius equation 
w_tilde = w/b  # Handling time dependent on temperature


# without the temperature function
H1_star = r_n^3/d_h^3 
S1_star = 0 

H2_star = r_n^3/(d_h^3*n_hn^3)
S2_star = 0 

H3_star = -(d_s^3*n_hn^3*n_sn^3*r_c^3*y_n^3)/(a*d_s^2*n_hn*n_sn*w*y_n-a*d_s*n_hn*n_sn*p*w*y_n+a*d_s^2*n_sn^2*w-a*d_s*n_hn^2*y_n^2+a*n_hn^2*p*y_n^2-a*d_s*n_hn*n_sn*y_n-d_h*d_s*n_hn*n_sn*y_n)^3
S3_star = (a*d_s^2*n_hn^3*n_sn^2*r_c^3*y_n^3*(d_s*n_sn*w-n_hn*y_n))/(a*d_s^2*n_hn*n_sn*w*y_n-a*d_s*n_hn*n_sn*p*w*y_n+a*d_s^2*n_sn^2*w-a*d_s*n_hn^2*y_n^2+a*n_hn^2*p*y_n^2-a*d_s*n_hn*n_sn*y_n-d_h*d_s*n_hn*n_sn*y_n)^3

H4_star = -r_n^3*y_n^3/(a*d_s*n_sn*w-a*n_hn*y_n-d_h*n_hn*y_n)^3
S4_star = (a*r_n^3*y_n^3*(d_s*n_sn*w-n_hn*y_n))/((a*d_s*n_sn*w-a*n_hn*y_n-d_h*n_hn*y_n)^3*d_s*n_sn)


M = rbind(c(H1_star,H2_star,H3_star,H4_star),c(S1_star,S2_star,S3_star,S4_star))
rownames(M) = c("H star", "S star")
colnames(M) = c("Eq.1", "Eq.2", "Eq.3", "Eq.4")

# with the temperature function
H1_star_temp = r_n^3/d_h^3 
S1_star_temp = 0 

H2_star_temp = r_n^3/(d_h^3*n_hn^3)
S2_star_temp = 0 

H3_star_temp = -(d_s^3*n_hn^3*n_sn^3*r_c^3*y_n^3)/(a*d_s^2*n_hn*n_sn*w_tilde*y_n-a*d_s*n_hn*n_sn*p*w_tilde*y_n+a*d_s^2*n_sn^2*w_tilde-a*d_s*n_hn^2*y_n^2+a*n_hn^2*p*y_n^2-a*d_s*n_hn*n_sn*y_n-d_h*d_s*n_hn*n_sn*y_n)^3
S3_star_temp = (a*d_s^2*n_hn^3*n_sn^2*r_c^3*y_n^3*(d_s*n_sn*w_tilde-n_hn*y_n))/(a*d_s^2*n_hn*n_sn*w_tilde*y_n-a*d_s*n_hn*n_sn*p*w_tilde*y_n+a*d_s^2*n_sn^2*w_tilde-a*d_s*n_hn^2*y_n^2+a*n_hn^2*p*y_n^2-a*d_s*n_hn*n_sn*y_n-d_h*d_s*n_hn*n_sn*y_n)^3

H4_star_temp = -r_n^3*y_n^3/(a*d_s*n_sn*w_tilde-a*n_hn*y_n-d_h*n_hn*y_n)^3
S4_star_temp = (a*r_n^3*y_n^3*(d_s*n_sn*w_tilde-n_hn*y_n))/((a*d_s*n_sn*w_tilde-a*n_hn*y_n-d_h*n_hn*y_n)^3*d_s*n_sn)


M_temp = rbind(c(H1_star_temp,H2_star_temp,H3_star_temp,H4_star_temp),c(S1_star_temp,S2_star_temp,S3_star_temp,S4_star_temp))
rownames(M_temp) = c("H star", "S star")
colnames(M_temp) = c("Eq.1", "Eq.2", "Eq.3", "Eq.4")

