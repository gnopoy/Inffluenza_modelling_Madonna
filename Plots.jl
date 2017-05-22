using PyPlot
using ODE

### Import the IFN Data from the current folder ####
A = readdlm("./Extended_IFN_Pony_1.txt")
B = readdlm("./IFN_Pony_1.txt")
E = readdlm("./Experimental_data_Pony_1_IFN.txt")

### Import the Viremia Data from the current folder ####
C = readdlm("./Extended_viremia_Pony_1.txt")
D = readdlm("./viremia_Pony_1.txt")
F = readdlm("./Experimental_data_Pony_1_Viremia.txt")

### Make the Plots ###
# time and IFN variable
t1 = A[:, 1];
A1 = A[:, 2];
t2 = B[:, 1];
B2 = B[:, 2];
t5 = E[:, 1];
E5 = E[:, 2];

# time and Viremia variable
t3 = C[:, 1];
C3 = C[:, 2];
t4 = D[:, 1];
D4 = D[:, 2];
t6 = F[:, 1];
F6 = F[:, 2];


#### plots
# Viremia plots
figure()
title("Viremia in Pony 1")
plot(t3, C3, label="Extended Model")
plot(t4, D4, label="Paper's Model")
scatter(t6, F6, label="Experimental Data")
xlabel("Time in days")
ylabel("Viral titers in RNA copies/NS seceretion.ml")
legend()

# IFN plots
figure()
title("Interferon in Pony 1")
plot(t1, A1, label="Extended Model")
plot(t2, B2, label="Paper's Model")
scatter(t5, E5, label="Experimental Data")
xlabel("Time in days")
ylabel("Interferon fold change")
legend()
