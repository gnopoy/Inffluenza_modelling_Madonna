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
  beta = 8.3*10.0^(-6);  # (RNA copy)^-1*(ml of Nasal Secretion)*day^-1
  phi = 6.9*10.0^(-2) ;  # (IFN fold change)^-1*day^-1
  rho = 1.0*10.0^(-2) ;   # day^-1
  delta_I = 2.0 ;      # day
  kappa = 1.6    ;     # (IFN fold change)^-1*day^-1
  p = 7.7*10.0^(-5);     # RNA copies (ml of Nasal Secretion)^-1*day^-1*cell^-1
  c = 20.0    ;        # day^-1
  q = 6.1*10.0^(-10.0);    # (IFN fold change)^1*day^-1
  d = 0.85    ;         # No units
  omega = 1.0 ;       # No units
  mu = 7.0;            # days
  kb = 0.20*3600.0*24.0 # binding constant (virus^-1*day^-1)
  kd =                # dissociation constant
  kf = 0.55*3600.0*24.0  # fusion rate constant


  ### ydot[1]=T, ydot[2]=I, ydot[3]=R, ydot[4]=V, ydot[5]=F ###
  ### Model needs to be solved before and after onset of immune response
  if t<mu
    #delta_A = delta_I*exp(omega*(t-mu))
    delta_A=delta_I
    ydot[1]= rho*y[3]+kd*y[6]-kb*y[4]*y[1]-phi*y[5]*y[1]
    ydot[2]= kf*y[6]-delta_A*y[2]-kappa*y[2]*y[5]
    ydot[3]= phi*y[5]*y[1]-rho*y[3]
    ydot[4]= p*y[2]-c*y[4]-kb*y[1]*y[4]+kd*y[6]
    ydot[5]= q*y[2]-d*y[5]
    ydot[6]= kb*y[1]*y[4]-kd*y[6]-kf*y[6]
  else
    delta_A = delta_I*exp(omega*(t-mu))
    #delta_A=delta_I
    ydot[1]= rho*y[3]+kd*y[6]-kb*y[4]*y[1]-phi*y[5]*y[1]
    ydot[2]= kf*y[6]-delta_A*y[2]-kappa*y[2]*y[5]
    ydot[3]= phi*y[5]*y[1]-rho*y[3]
    ydot[4]= p*y[2]-c*y[4]-kb*y[1]*y[4]+kd*y[6]
    ydot[5]= q*y[2]-d*y[5]
    ydot[6]= kb*y[1]*y[4]-kd*y[6]-kf*y[6]
  end

  ydot
end

### Solve the system of differential Equations

# Define time span
t1 = 0.0:0.1:14.0                       # days

### Initial values for the model ###
T = 3.5*10.0^(11)  ; # Cells
V = 0.1  ;           # RNA copy *(Nasal Secretion ml)^-1
R = 0.0    ;         # Cells
F = 1.0  ;           # (IFN fold change)
I = 0.0   ;          # Cells
TV= 0.0 ;           # Cells-virus bound


# Matching initial conditions
y0 = [T, I, R, V, F, TV]           # initial conditions
t1, y1 = ODE.ode23s(Dynamics, y0, t1)   # Solver
A1 = map(v->v[1], y1)                    # Accessing the T
A2 = map(v->v[2], y1)                    # Accessing the I
A3 = map(v->v[3], y1)                    # Accessing the R
A4 = map(v->v[4], y1)                    # Accessing the V
A5 = map(v->v[5], y1)                    # Accessing the F
A6 = map(v->v[6], y1)                    # Accessing the TV

figure()
plot(t1, A1, label="Target Cells")
plot(t1, A2, label="Infected Cells")
plot(t1, A3, label="Refractory Cells")
plot(t1, A6, label= "Cells with virus bound")
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
