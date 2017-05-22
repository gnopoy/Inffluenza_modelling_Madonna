### Root Mean Squared Calculation ####

# Import relevant Data set
A = readdlm("./data/Extended_IFN_Pony_1.txt")
B = readdlm("./data/Extended_viremia_Pony_1.txt")
C = readdlm("./data/Experimental_data_Pony_1_IFN.txt")
D = readdlm("./data/Experimental_data_Pony_1_Viremia.txt")

# Collect just the Data points
A1 = A[:, 2]
B1 = B[:, 2]
C1 = C[:, 2]
D1 = D[:, 2]

# Do the Mean Square calculation
# sum of the squares for Viremia and IFN

RMS1 = 0.0

for i=1:10
  RMS1 = RMS1 + (log10(B1[i])-log10(D1[i]))^2;
  @show RMS1
end

RMS2 = 0.0

for i=1:6
  RMS2 = RMS2 + (A1[i]-C1[i])^2;
  @show RMS2
end

RMS = sqrt((1/10)*RMS1 + (1/6)*RMS2)
