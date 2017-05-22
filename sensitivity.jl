using PyPlot
using ODE

### Import the IFN Data from the current folder ####
# A = readdlm("./data/Extended_IFN_Pony_5_Kfup.txt")
# B = readdlm("./data/Extended_IFN_Pony_5_Kfdown.txt")
# E = readdlm("./data/Extended_IFN_Pony_5.txt")
# G = readdlm("./data/Experimental_data_Pony_3_IFN.txt")

### Import the Viremia Data from the current folder ####
C = readdlm("./data/Extended_viremia_Pony_6.txt")
D = readdlm("./data/Extended_viremia_Pony_6_kddown.txt")
Z = readdlm("./data/Extended_Viremia_Pony_6_kdup.txt")
F = readdlm("./data/Experimental_data_Pony_6_Viremia.txt")

### Make the Plots ###
# time and IFN variable
# t1 = A[:, 1];
# A1 = A[:, 2];
# t2 = B[:, 1];
# B2 = B[:, 2];
# t5 = E[:, 1];
# E5 = E[:, 2];
# t7 = G[:, 1]
# G7 = G[:, 2]

# time and Viremia variable
t3 = C[:, 1];
C3 = C[:, 2];
t4 = D[:, 1];
D4 = D[:, 2];
t6 = F[:, 1];
F6 = F[:, 2];
t8 = Z[:, 1];
Z8 = Z[:, 2];

### plots
# Viremia plots
figure()
title("Sensitivity of Viremia to Kd in Pony 6")
plot(t8, Z8, label="10*Kd")
plot(t3, C3, label="Kf")
plot(t4, D4, label="0.1*Kd")
scatter(t6, F6, label="Experimental Data")
xlabel("Time in days")
ylabel("Viral titers in RNA copies/NS seceretion.ml")
legend( loc="upper right")

# # IFN plots
# figure()
# title("Sensitivity of Interferon in Pony 5 to Kf")
# plot(t1, A1, label="Kf")
# plot(t2, B2, label="Kf*0.1")
# plot(t5, E5, label="Kf*10")
# scatter(t7, G7, label="Experimental Data")
# xlabel("Time in days")
# ylabel("Interferon fold change")
# legend()
