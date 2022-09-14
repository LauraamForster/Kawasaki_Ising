import random
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in file from KW_Ising_Part1.py 

filename='Matrices_for_density_.txt'

with open(filename) as fs:
    txtfile = []
    for line in fs:
        txtfile.append(line.strip())

size = txtfile[0].split()
box_height, box_width, air = int(size[2].split('=')[1]), int(size[3].split('=')[1]), int(size[4].split('=')[1]) #height of box (how tall top to bottom/ y coord)

total_array = np.array([])
for x in range(1, box_height+1):
    A = txtfile[x].split()
    line = []
    for y in range(box_width):
        B = float(A[y])
        line.append(B)
    # print(line)
    total_array = np.append(total_array, line)
S = np.reshape(total_array, (box_height, box_width)) 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Plot the density profile and fit to exponential

density = [] #create empty array
for x in range(box_height): 
    count = np.count_nonzero(S[x,:]==1) #counts the number of molecules in each row
    density.append(count/box_width) 
if air>101:
    ys = density[air:] #removes the air at the start of the box as it has evaporaated off
    xs = (np.arange(box_height))[air:] #makes arrays same size
else:
    ys = density
    xs = (np.arange(box_height))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#function to fit exponential
def monoExp(x, m, t, b): 
    return m * np.exp(-t * (x-air)) + b #returns fitted parameters

p0 = (100, 0.01, 100) # start with values near those expected
params, cv = scipy.optimize.curve_fit(monoExp, xs, ys, p0) #perform the fit
m, t, b = params 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the results

plt.bar(xs, ys, label="density") #bar chart of the density profile
plt.xlabel("Row Number") #axis labels
plt.ylabel("Density of Particles in Row") #axis labels
if air>101: #if simulation had evaporation, plot data with exponential fit 
    plt.plot(xs, monoExp(xs, m, t, b), '--', color='red', label="fitted exponential") #fitted exponential
else:
    plt.plot(xs, ys) #if no evaporation, plot just the density graph
plt.title("Histogram Showing Density of Particles in Each Row")
plt.show()

# inspect the equation of the fit
print(f"Y = {m} * e^(-{t} * x) + {b}")

