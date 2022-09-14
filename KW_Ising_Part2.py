import random
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

counter = 0
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# function to plot each output

def Plotter(total_array, counter):
    if counter % 10 == 0: #sets a counter so not every output is printed 
        plt.figure(figsize=(3.5, 7)) #sets figure size
        plot2= plt.imshow(total_array, cmap='gray') #plots lattice on graph in black and white
        plt.xticks([]) #removes axis labels
        plt.yticks([]) 
        plt.show()
    return counter
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in file from KW_Ising_Part1.py

filename='Matrices_for_analysis_.txt'

with open(filename) as fs:
    txtfile = []
    for line in fs:
        txtfile.append(line.strip())

# finds info from txt file such as box dimensions and initial conditions
size, IntConds, Energies = txtfile[0].split(), txtfile[1].split(), txtfile[2].split()
box_height, box_width = int(size[2].split('=')[1]), int(size[3].split('=')[1]) #height of box (how tall top to bottom/ y coord)
dil, temperature, air, run, evap_cycles = IntConds[2].split('=')[1], IntConds[3].split('=')[1], int(IntConds[4].split('=')[1]), int(IntConds[5].split('=')[1]), int(IntConds[6].split('=')[1])
meanE, meanEsq, varE = float(Energies[1].split('=')[1]), float(Energies[2].split('=')[1]), float(Energies[3].split('=')[1])


total_array = np.array([])
length = len(txtfile) 
for x in range(3, length):
    if 'number' in txtfile[x]:
        leng = len(total_array)
        if leng != 0: #if the array contains some data
            counter  = counter + 1
            print(counter)
            total_array = np.reshape(total_array, (box_height, box_width)) 
            Cluster_res = Plotter(total_array, counter)
            total_array = np.array([])
    else: 
        A = txtfile[x].split()
        line = []
        for y in range(box_width):
            B = float(A[y])
            line.append(B)
        total_array = np.append(total_array, line)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Outputs some info such as energy and peclet number for the simulation

Pe = box_height/evap_cycles #define the peclet number
print(' energy   averaged over run ',-meanE,' at kT/J ',temperature)
print(' energy^2 averaged over run ',meanEsq)
print(' variance of energy         ',varE)
print(' Peclet number              ',Pe)

