import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from statistics import *

counter = 0
graphmch = []
mcs = []

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to plot the output of the file to display location of clusters

def Plotter(col, row, box_width, box_height, n_big_comp, counter, listy):
    print('counter,', counter) #print what number counter is up to

    fig = plt.figure(figsize=(3.5, 7)) #sets the size of the figure pane
    ax1 = fig.add_subplot(1,1,1) #puts all subplots overlaid

    plt.scatter(col,row,alpha=0.5,linewidths=0.0) #plot the particles in the system
    
    plt.xticks([]) #removes x axis values
    plt.yticks([]) #removes y axis values
    plt.ylim([box_height+1, -1]) #sets height of box (with values in ascending order from 0 at top to height at bottom, to match matrix). (An extra point added each way so spots dont hang over edge)
    plt.xlim([-1,box_width+1]) #sets width of box (An extra point added each way so spots dont hang over edge)
    plt.tight_layout() #fits to size of pane
    #
    for i in range(0,n_big_comp):
        for ip in listy[i]: 
            plt.scatter(col[ip],row[ip],alpha=0.5,linewidths=0.0,c='red') #display the connected components in red to show them clearly
    #
    plt.xticks([]) #removes x axis values
    plt.yticks([]) #removes y axis values
    plt.ylim([box_height+1, -1]) #sets height of box (with values in ascending order from 0 at top to height at bottom, to match matrix). 
                                #(An extra point added each way so spots dont hang over edge)
    plt.xlim([-1,box_width+1]) #sets width of box (An extra point added each way so spots dont hang over edge)
    plt.tight_layout()
    plt.show()

    return counter

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the mean height and size of clusters at each interval

def meanheightcalc(row, n_big_comp, listy, graphmch, mcs): #function finds the mean point where all the clusters sit, and finds mean size of clusters
    #if clusters form towards the air interface, then value will be low to denote row number they on average sit at
    mch=[]
    clussize=[]
# determines the mean cluster height
    for i in range(0,n_big_comp): #range of clusters over threshold size
        cluster_height = 0 #empty array for each time
        for ip in listy[i]: 
            cluster_height = cluster_height + row[ip] #add height location of particles together
        clustersize = len(listy[i]) #determines number of points in each cluster
        clussize.append(clustersize) 
        meanheight = 0 #empty variable each time
        meanheight = cluster_height/clustersize #finds mean height of each cluster
        mch.append(meanheight) #adds mean height to array 

# determines the mean cluster size
    mean_clusterh2 = 0 #empty the variable each time
    for i in range(n_big_comp): #in range of clusters over threshold size
        mean_clusterh2 = mch[i] + mean_clusterh2 #add cluster means together
    numbclusts = len(clussize) #number of clusters 
    for i in range(numbclusts):
        clustersize = clussize[i] + clustersize 
    meanclussize = clustersize/numbclusts #finds the mean size
    print('mean cluster size in system is', meanclussize) #prints the mean size of the clusters at each output
    mcs.append(meanclussize)
    
    heightinsystem = mean_clusterh2/n_big_comp #divide by number of clusters
    graphmch.append(heightinsystem) 
    print('mean cluster height in system is', heightinsystem) #prints the mean height of the clusters at each output

    return heightinsystem, graphmch, mcs

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the mean size of clusters 
def plotmeansize(listy):
    sizesofclus=[]
    length = len(listy)
    for x in range(length):
        value = len(listy[x])
        sizesofclus.append(value)
    return sizesofclus


#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to find the clusters within the system and quantify them - uses NetworkX package

def Cluster(array, counter, graphmch, mcs):
    print(' ') #creates a gap to seperate print outs for each part of the system to make it easier to see what output is for what part

    row=total_array[:,0] #defined as the first value in the (x,y coords) and defines the row number (which row from top to bottom)
    col=total_array[:,1] #defined as the second value in the (x,y coords) and defines the column number (which column from left to right)

    N=len(row) #gives the number of particles in the system
    print('number of particles ',N)

    G=nx.Graph() # create networkx graph G
    
    for i in range(0,N): #adds a node for each particle in system
        G.add_node(i) # add nodes
        
    for i in range(0,N):
        for j in range(0,N):
            if(i != j): 
                # rsq=(x[i]-x[j])**2+(y[i]-y[j])**2
                if( np.abs(row[i]-row[j]) < 1.01 and np.abs(col[i]-col[j]) < 1.01 ): #add connections if there are 2 points next to each other
                    G.add_edge(i,j) # add connection one way
                    G.add_edge(j,i) # add connection the opposite way

    # print(list(nx.connected_components(G))) 
    n_connect_comp=nx.number_connected_components(G) #calculates the number of connected components
    print('total number of connected components ',n_connect_comp)

    print('sorting connected components fromn biggest to smallest')
    listy=sorted(nx.connected_components(G), key=len, reverse=True) #sorts the connected components from biggest to smallest
    
    list_clust_sizes = plotmeansize(listy) 

    biggest=len(listy[0]) #finds the biggest component's size
    print('biggest component has size ',biggest) 

    big_threshold=10 #defines a threshold of the number of connected components 
    print('counting connected components with at least ',big_threshold,' particles')

    for i in range(0,n_connect_comp): 
        if( len(listy[i])<big_threshold): #if the threshold is bigger than the coordinate for the biggest component then break
            n_big_comp=i
            break

    print(n_big_comp,' connected components with at least ',big_threshold,' particles')

    Cluster_res= Plotter(col, row, box_width, box_height, n_big_comp, counter, listy)
    meanheight, graphmch, mcs= meanheightcalc(row, n_big_comp, listy, graphmch, mcs)

    return list_clust_sizes, graphmch, mcs
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# read in file to analyse

filename='Matrices_During_Ising.txt'

with open(filename) as fs:
    txtfile = []
    for line in fs:
        txtfile.append(line.strip())

size = txtfile[0].split()
box_height = int(size[2].split('=')[1]) #height of box (how tall top to bottom/ y coord)
box_width = int(size[3].split('=')[1]) #width of box (how wide left to right/ x coord)

total_array = np.array([]) #create an empty array
length = len(txtfile) 
for i in range(length):
    if 'Matrix' in txtfile[i]: #each matrix contains a title to be read in
        leng = len(total_array)
        total_array = np.reshape(total_array, (int(leng/2), 2)) 
        total_array = total_array.astype(int)
        if leng != 0: #if the array contains some data
            counter = counter + 1
            if counter % 1 == 0:
                cluster_sizes, graphmch, mcs = Cluster(total_array, counter, graphmch, mcs) #call the Cluster Function to find connected components
        total_array = np.array([]) #empty the array 
    else:
        split = txtfile[i].split() #split the text file by each line for each coordinate
        col_coord = int(split[0]) #find the column coordinate (from left to right)
        row_coord = int(split[1]) #find the row coordinate (from top to bottom)
        point = [(col_coord, row_coord)] #define each particle as a point with coordinates
        total_array = np.append(total_array, point) #add the point to the array

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#display the histogram for the final output, displays amount of clusters at each size

sortedclusters = sorted(cluster_sizes) #sorts clusters by size
freq = {}
for item in sortedclusters: #creates a dictionary 
    if (item in freq):
        freq[item] += 1
    else:
        freq[item] = 1
#prints output of histogram to see exactly how many clusters at each cluster size
for key, value in freq.items(): 
    print ("cluster size", "% d : % d"%(key, value))

plt.hist(cluster_sizes) #plot the particles in the system  
plt.xlabel("Cluster Size")#axis label
plt.ylabel("Amount of Clusters")#axis label
plt.title("Histogram of Cluster Size")#sets graph title
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# display cluster height over time
lengthgraph = len(graphmch) 
timeaxis = np.arange(1, lengthgraph+1) #sets x axis to plot against

a, b = np.polyfit(timeaxis, graphmch, 1) #fits a straight line to data

plt.plot(timeaxis, graphmch) #plots the data
plt.plot(timeaxis, a*timeaxis+b, color='red') #plots the fit in another colour
plt.ylim([500, 0]) #reverses the axis to show the row number of box

plt.xlabel("Number of Run Cycles", fontsize=15) #axis label
plt.ylabel("Row Number", fontsize=15) #axis label
plt.xticks(fontsize=15) #makes labels larger
plt.yticks(fontsize=15) #makes labels larger
plt.title("Average Height of Clusters over Time", fontsize=20) #sets graph title
leg = plt.legend(loc="upper left", fontsize=20, frameon=False,)
# get the individual lines inside legend and set line width
for line in leg.get_lines():
    line.set_linewidth(6)
# get label texts inside legend and set font size
for text in leg.get_texts():
    text.set_fontsize('x-large')
plt.show()

print('equation of fit is=', -a, "x + ", b) #gives equation for line of best fit

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# display cluster size over time
lengthgraphcs = len(mcs)
timeaxis2 = np.arange(1, lengthgraphcs+1) #sets x axis to plot against

c, d = np.polyfit(timeaxis2, mcs, 1) #fits a straight line to data

plt.plot(timeaxis2, mcs, color='blue', label="No Evaporation") #plots the data
plt.plot(timeaxis2, a*timeaxis2+b, color='lightblue', label="Fit - No Evaporation") #plots the fit in another colour

plt.xlabel("Time", fontsize=15)
plt.ylabel("Mean Cluster Size", fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title("Average Size of Clusters over Time - Temp 10.5", fontsize=20)
leg = plt.legend(loc="upper left", fontsize=20, frameon=False,)
# get the individual lines inside legend and set line width
for line in leg.get_lines():
    line.set_linewidth(6)
# get label texts inside legend and set font size
for text in leg.get_texts():
    text.set_fontsize('x-large')
plt.show()

print('equation of fit is=', c, "x + ", d) #gives equation for line of best fit



