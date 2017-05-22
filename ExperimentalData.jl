#### Experimental Data for the Models
using PyPlot
# time points
t =[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] # Viremia
t1 =[0.0, 1.0, 2.0, 3.0, 4.0, 5.0] # IFN titers

### Pony 1
A1 =[630.9573445, 100000, 2511.886432, 630.9573445, 251.1886432, 158.4893192, 398.1071706, 1.0, 1.0, 1.0]
B1 =[1.0, 1.4, 2.0, 5.0, 2.0, 0.6]
# data_array1 = [t A1];
# data_array11 =[t1 B1];
# path_to_file1 ="./Experimental_data_Pony_1_Viremia.dat"
# path_to_file11 ="./Experimental_data_Pony_1_IFN.dat"
# writedlm(path_to_file1, data_array1)
# writedlm(path_to_file11, data_array11)


### Pony 2
A2 =[1.0, 63095734.4, 25118.9, 15848.9, 398107.2, 398107.2, 251.2, 251.2, 1.0, 1.0]
B2= [1.0, 0.4, 5.0, 0.6, 0.1, 0.2]
# data_array2 = [t A2];
# data_array22 = [t1 B2];
# path_to_file22 ="./Experimental_data_Pony_2_IFN.dat"
# path_to_file2 ="./Experimental_data_Pony_2_Viremia.dat"
# writedlm(path_to_file22, data_array22)
# writedlm(path_to_file2, data_array2)


### Pony 3
A3 =[1.0, 39810.7, 398.1, 501.2, 501.2, 6309.6, 316.2, 1.0, 1.0, 1.0]
B3= [1.0, 0.4, 1.2, 0.8, 0.2, 0.4]
# data_array3 = [t A2];
# data_array33 = [t1 B2];
# path_to_file33 ="./Experimental_data_Pony_3_IFN.dat"
# path_to_file3 ="./Experimental_data_Pony_3_Viremia.dat"
# writedlm(path_to_file33, data_array33)
# writedlm(path_to_file3, data_array3)

### Pony 4
A4 =[1.0, 39810.7, 398.1, 501.2, 501.2, 6309.6, 316.2, 1.0, 1.0, 1.0]
B4= [1.0, 1.0, 5.8, 2.6, 0.2, 0.4]
data_array4 = [t A4];
data_array44 = [t1 B4];
path_to_file44 ="./data/Experimental_data_Pony_4_IFN.dat"
path_to_file4 ="./data/Experimental_data_Pony_4_Viremia.dat"
writedlm(path_to_file44, data_array44)
writedlm(path_to_file4, data_array4)

### Pony 5
A5 =[3981.1, 1584893.2, 251.2, 1000.0, 1584.9, 631.0, 1.0, 1.0, 1.0, 1.0]
B5= [1.0, 0.8, 5.6, 0.4, 0.2, 0.3]
# data_array5 = [t A5];
# data_array55 = [t1 B5];
# path_to_file55 ="./Experimental_data_Pony_5_IFN.dat"
# path_to_file5 ="./Experimental_data_Pony_5_Viremia.dat"
# writedlm(path_to_file55, data_array33)
# writedlm(path_to_file5, data_array3)

### Pony 6
A6 =[1.0, 10000.0, 398.1, 1000.0, 251.2, 398.1, 1.0, 1.0, 1.0, 1.0]
B6= [1.0, 1.0, 10.8, 1.8, 0.5, 0.2]
# data_array6 = [t A6];
# data_array66 = [t1 B6];
# path_to_file66 ="./Experimental_data_Pony_6_IFN.dat"
# path_to_file6 ="./Experimental_data_Pony_6_Viremia.dat"
# writedlm(path_to_file66, data_array66)
# writedlm(path_to_file6, data_array6)

# 
# figure()
# scatter(t,A6,label="Pony_6")
# ylabel("Viral load in RNA copies (ml NS. day)^-1")
# xlabel("Time post-challenge in days")
# legend()
#
# figure()
# scatter(t1,B6,label="Pony_6")
# ylabel("IFN fold change")
# xlabel("Time post-challenge in days")
# legend()
