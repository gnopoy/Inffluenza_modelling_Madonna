### Code for Advanced Biomolecular Engineering Project ###
### Code is used to describe the within-host dynamics of an influenza A infection ###
### Model incorporates both latent infection and Immune Response ###
### Extension to the model is provided by inclusion of kinetic parameters for binding and fusion of viral particles ###
### Source paper is : Modelling Within-Host Dynamics of Influenza Virus Infection Including Immune Responses, by Pawelek et al. ###

### Packages to be used in this code ###
using ODE
using PyPlot

### Dictionary of model parameters, see Table 1 ####
# beta          : Infection rate
# phi           : IFN-induced antiviral efficacy
# rho           : Reversion rate from refrcatory
# delta         : Death rate
# kappa         : Killing rate of infected cells by NK cells
# p             : Viral production rate
# c             : Clearance rate of free virions
# q             : Production rate of interferon
# d             : Decay rate of interferon
# delta_A       : Time-varying death rate of infected cells during the adaptive immune response
# delta_I       : Death rate of infected cells before the adaptive immune response emerges
# I             : Infected epithelial cells
# T             : Uninfected epithelial cells that are suceptible to infection
# V             : Viral load
# R             : Epithelial cells in the refractory state
# F             : Interferon
# omega         : Drug efficacy
# t             : time

### Set up model for infection dynamics ###

function Dynamics(t, y)
  ydot=similar(y)

  ### Model Parameters values found in the paper ###
  ### Parameters are specific to each Pony, see Table 2 ###
  Beta = [8.3*10.0^(-6), 1.1*10.0^(-7), 8.5*10.0^(-5), 1.2*10.0^(-6), 2.8*10.0^(-7), 1.9*10.0^(-4)];
  Phi = [6.9*10.0^(-2), 2.2*10.0^(-2), 1.2*10.0^(0), 1.1*10.0^(-1), 5.3*10.0^(-1), 2.5*10.0^(-2)];
  Rho = [1.0*10.0^(-2), 1.0*10.0^(-1), 6.7*10.0^(-2), 5.1*10.0^(0), 3.7*10.0^(-2), 2.0*10.0^(-1)];
  Kappa = [1.6, 5.7, 11, 1.2, 3.2, 2.4];
  P = [7.7*10.0^(-5), 2.4*10.0^(-2), 1.6*10.0^(-5), 9.6*10.0^(-4), 6.8*10.0^(-3), 8.1*10.0^(-6)];
  C = [20, 14, 9.8, 19, 20, 5.8];
  Q = [6.1*10.0^(-10), 4.5*10.0^(-10), 3.3*10.0^(-10), 3.8*10.0^(-10), 1.9*10.0^(-9), 2.1*10.0^(-9)];
  D = [0.85, 2.2, 2.0, 1.9, 1.9, 2.4];
  Omega = [1.0, 0.37, 0.92, 0.23, 1.2, 2.2];

  ### Chose the Pony      ###
  i = 6;

  ### Setting the Parameters for the right Pony ###
  beta = Beta[i];           # (RNA copy)^-1*(ml of Nasal Secretion)*day^-1
  phi = Phi[i];             # (IFN fold change)^-1*day^-1
  rho = Rho[i];             # day^-1
  delta_I = 2.0;            # day
  kappa = Kappa[i]   ;      # (IFN fold change)^-1*day^-1
  p = P[i];                 # RNA copies (ml of Nasal Secretion)^-1*day^-1*cell^-1
  c = C[i];                 # day^-1
  q = Q[i];                 # (IFN fold change)^1*day^-1
  d = D[i];                 # No units
  omega = Omega[i] ;        # No units
  mu = 7.0;                 # days

  ### ydot[1]=T, ydot[2]=I, ydot[3]=R, ydot[4]=V, ydot[5]=F ###
  ### Model needs to be solved before and after onset of immune response
  if t<mu
    #delta_A = delta_I*exp(omega*(t-mu))
    delta_A=delta_I
    ydot[1]= rho*y[3]-beta*y[4]*y[1]-phi*y[5]*y[1]
    ydot[2]= beta*y[4]*y[1]-delta_A*y[2]-kappa*y[2]*y[5]
    ydot[3]= phi*y[5]*y[1]-rho*y[3]
    ydot[4]= p*y[2]-c*y[4]
    ydot[5]= q*y[2]-d*y[5]
  else
    delta_A = delta_I*exp(omega*(t-mu))
    #delta_A=delta_I
    ydot[1]= rho*y[3]-beta*y[4]*y[1]-phi*y[5]*y[1]
    ydot[2]= beta*y[4]*y[1]-delta_A*y[2]-kappa*y[2]*y[5]
    ydot[3]= phi*y[5]*y[1]-rho*y[3]
    ydot[4]= q*y[2]-d*y[5]
    ydot[5]= q*y[2]-d*y[5]
  end

  ydot
end

### Solve the system of differential Equations

# Define time span
t1 = 0.0:0.1:14.0                       # days

### Initial values for the model ###
T = 3.5*10.0^(11)  ; # Cells
V = 100.0  ;           # RNA copy *(Nasal Secretion ml)^-1
R = 0.0    ;         # Cells
F = 1.0  ;           # (IFN fold change)
I = 0.0   ;          # Cells


# Matching initial conditions
y0 = [T, I, R, V, F]           # initial conditions
t1, y1 = ODE.ode23s(Dynamics, y0, t1)   # Solver
A1 = map(v->v[1], y1)                    # Accessing the ydot[x]
A2 = map(v->v[2], y1)                    # Accessing the ydot[x]
A3 = map(v->v[3], y1)                    # Accessing the ydot[x]
A4 = map(v->v[4], y1)                    # Accessing the ydot[x]
A5 = map(v->v[5], y1)                    # Accessing the ydot[x]

figure()
plot(t1, A1, label="Target Cells")
plot(t1, A2, label="Infected Cells")
plot(t1, A3, label="Refractory Cells")
xlabel("time in days")
ylabel("Cells")
legend()

figure()
plot(t1, A4, label="Viremia")
xlabel("Time in days")
ylabel("Viral titers in RNA copies /Nasal secretion/ml")
legend()

figure()
plot(t1, A5, label="IFN")
ylabel("Inteferon fold change")
xlabel("Time in days")
legend()

### Saving Files needed for later modelling
file_path = "./viremia_Pony_6.txt"
data_array = [t1 A4]
writedlm(file_path, data_array)

file_path2 = "./IFN_Pony_6.txt"
data_array2 = [t1 A5]
writedlm(file_path2, data_array2)
