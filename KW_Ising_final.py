import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#----------------
#Define some initial conditions

Lhe = 200 #define system size height
Lwid = 100 #define system size width
T = 10.5 #temperature of system
dil = 0.05 #dilution of particles in system
run = 300*(Lwid*Lhe) #runtime
air = 2 #define amount of system which is air 
Evap_cycle = 10 #number of evaporation cycles
Evap_steps=Evap_cycle*Lwid*Lhe #no. of steps for evaporation is no. of cycles per lattice site
count, EngCount, accumE, accumEsq = 0, 0, 0.0, 0.0
Eng_arr = []
Images = []

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating the system

#Create a lattice (of defined size, with random points determined as 1 or -1)
S=np.zeros((Lhe, Lwid), dtype=int) #creates empty lattice
for i in range(Lhe):
        for j in range(Lwid):
                r=random.uniform(0, 1) #determines a random number between 0 and 1
                if r<=dil: S[i, j]=1 #if random number is less than .5 change point in lattice to 1
                else: S[i, j]=-1 #if not, make point in lattice -1
                if i < air: #the defined parts of the system which are air
                    S[0:Lwid] = 0 #are set as 0
                   
px, py =np.arange(Lhe, dtype=int)+1, np.arange(Lwid, dtype=int)+1 #create an array of same size as lattice height, integers in ascending order, add 1 to make it start from 1
px[-1]=0 #make the final value 0
py[-1]=0 #make the final value 0

mx, my=np.arange(Lhe, dtype=int)-1, np.arange(Lwid, dtype=int)-1 #create a second array of same size as height and one of width, integers in ascending order, minus 1 to make it start from -1
mx[0]=Lhe-1 #make first value lattice size-1
my[0]=Lwid-1 #make first value lattice size-1

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Define some functions to calculate energy and to simulate evaporation

def Energy(lattice): #function to calculate energy occasionally
    energy = 0 #reset energy variable to 0

    for a in range(Lhe): #for every point in the lattice
        for b in range(Lwid):
            First = lattice[a,b] #find the value of the point
            Surr = lattice[(a+1)%Lhe, b] + lattice[a,(b+1)%Lwid] + lattice[(a-1)%Lhe, b] + lattice[a,(b-1)%Lwid] #find point's neighbours and sum value of neighbours
            energy += Surr*First #multiply value of point, and its summed neighbours together, and add to energy value
    return energy/2 #return half the energy value

def Evapfunc(S, air): #function to simulate evaporation
    check_before=np.count_nonzero(S==1)  #counts the number of particles in system 
    succeed=True

    row = S[air,:] #defines a whole row, not including those which are air

    for j in range(0,Lwid): 
        n_vacancies_column=np.count_nonzero(S[air:Lhe,j]==-1) #count number of spaces without particles in that column 
        print('number vacancies in column j ',n_vacancies_column)  # prints no. of vacancies
        if(n_vacancies_column==0): #if there are no spaces
            print('evaporation step impossible, system jammed!') #prints that system is jammed
            succeed=False 
            break
        else:
            print('evaporating in column ',j) #prints which column the evaporation is on
            if(S[air,j]==1): #checks if the top value of column contains particle
                for i in range (air+1,Lhe): #checks all the rows under the air rows 
                   if S[i,j] == -1: #if there is a space in that column
                       S[i,j]=1 #the particle can move down to the nearest one
                       break
            S[air,j]=0 # converts the particle in top of column to air to complete evaporation
 
    check_after=np.count_nonzero(S==1) #counts the number of particles in the system
    print('before and after number of particles ',check_before, check_after)
    if(check_after != check_before): #checks the number of particles before and after evap is conserved
        print('PANIC, particle number not conserved!') 
    Images.append(S) #appends a copy of the matrix to a list so print outs of system are shown at the end
    return S, air, succeed, Images

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Kawasaki Dynamics Ising Model Swaps

for M in range(run): 
        count=count+1 #set a counter

        if count % Evap_steps == 0: #every so often call evaporation
            if air < Lhe-2: #if system isn't entirely air 
                S_tmp, air, evap_success, Images = Evapfunc(S, air) #call evaporation function
                if(evap_success): #if evaporation was successful
                    air = air + 1 #redefine air to represent another row turned to air
                    S=np.copy(S_tmp) #define S as new matrix

        i, j=random.choice(range(Lhe)), random.choice(range(Lwid)) #choose a random point in matrix 
        neb=random.choice([[0, 1], [0, -1], [-1, 0], [1, 0]]) #randomly chooses which neighbour to try and swap with

        ineb=i+neb[0] #defines the coordinates of the chosen neighbouring particle
        jneb=j+neb[1]

        #sets some boundary conditions, so when an edge particle is chosen, the opposing edge can be swapped with instead
        if ineb>Lhe-1: #if point is on bottom of matrix
            ineb=0 #switch neighbour to top point 
        if jneb>Lwid-1: #if point is on right hand side of matrix  
            jneb=0 #switch to left hand side
        if ineb<0: # if point is at the top of matrix
            ineb=Lhe-1 #switch to bottom
        if jneb<0: #if point is on left hand side of matrix
            jneb=Lwid-1 #switch with right hand side


        if S[i, j] != 0 and S[ineb, jneb] != 0: #If chosen particles are not air particles
            deltaE1=(S[i, j]-S[ineb, jneb])*(S[px[i], j]+S[i, py[j]]+S[mx[i], j]+S[i, my[j]]-S[ineb, jneb]) #calculate current energy of matrix
            deltaE2=(S[ineb,jneb]-S[i,j])*(S[px[ineb], jneb]+S[ineb,py[jneb]]+S[mx[ineb], jneb]+S[ineb, my[jneb]]-S[i,j]) #calculate the new energy of new matrix after swap
            deltaE=deltaE1+deltaE2 #add energies together
            W=np.exp(-deltaE/T) 
            r=random.random()
            if r<W: #Use monte carlo to determine if the swap will be ok
                    itemp=np.copy(S[i,j]) #make a copy of origional point
                    S[i,j]=S[ineb,jneb] #swap origional point with chosen neighbouring point
                    S[ineb,jneb]=itemp #set chosen neighbouring point to origional value of point
                    if count % Lwid**2 == 0: #occasionally calculate the energy
                        EngCount = EngCount + 1 #set a counter
                        Eng = Energy(S) #call energy function to calculate lattice energy
                        accumE = accumE + float(Eng) #define accumulated energy
                        accumEsq = accumEsq + float(Eng*Eng) #define accumulated energy squared
                        Eng_arr.append(Eng)

#display the energy and variance over the run
meanE=accumE/EngCount 
meanEsq=accumEsq/EngCount
varE=meanEsq-meanE**2
print(' energy   averaged over run ',-meanE,' at kT/J ',T)
print(' energy^2 averaged over run ',meanEsq)
print(' variance of energy         ',varE)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Display the evolution of the lattice over the run, and outputs for analysis

Eng_arr=np.array(Eng_arr) #create an array for graph axis
Pe = Lhe/Evap_cycle #define the peclet number
print('Pe ',Pe) #print peclet number
length = len(Images)

density = []
for x in range(Lhe):
    count = np.count_nonzero(S[x,:]==1)
    density.append(count)

yax = np.arange(Lhe)
plt.bar(yax, density)
plt.show()

f = open('./Matrices_During_Ising.txt', 'w') #open a textfile with defined name
f.write('Matrix Dimensions: ' + 'height=' + str(Lhe) + '   ' + 'width=' +str(Lwid) + "\n")
for y in range(length):
    S_curr = Images[y]
    particles = np.count_nonzero(S_curr==1)
    f.write('New Matrix no. ' + str(y+1) + ' / no. of particles ' + str(particles) + ' / matrix size ' + str(Lhe) + ' x ' + str(Lwid) +"\n")
    for i in range(Lhe):
        for j in range(Lwid):
            if S_curr[i, j] == 1: #for all the particles in the system
                f.write(str(i) + '   ' + str(j) + "\n") #write each location onto the text file in columns 
f.close() #close the text file

Time_pnts = EngCount  #define number of points to plot over
X_ax = run/(Lwid*Lhe) #Define maximum value on x axis
time = np.linspace(0, X_ax, Time_pnts) #time axis points
print('Final Lattice', S) #print the final lattice

for x in range(length): #for every lattice after evaporation, display them sequentially to see evaporation occur over the run
      plt.xticks([])
      plt.yticks([])
      plt.imshow(Images[x], cmap='gray') 
      plt.show() 

# plot2= plt.imshow(S, cmap='gray') #display only the final lattitce
#plt.show()