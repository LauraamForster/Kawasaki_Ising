import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1,1,1)

filename='Matrices_During_Ising.txt'

with open(filename) as fs:
    txtfile = []
    for line in fs:
        txtfile.append(line.strip())

size = txtfile[0].split()
box_height = int(size[2].split('=')[1]) #height of box (how tall top to bottom/ y coord)
box_width = int(size[3].split('=')[1]) #width of box (how wide left to right/ x coord)

def Cluster(array):
    print(' ')

    col=total_array[:,0] #defined as the first value in the (x,y coords) and defines the column number (which column from left to right)
    row=total_array[:,1] #defined as the second value in the (x,y coords) and defines the row number (which row from top to bottom)

    N=len(col) #gives the number of particles in the system
    print('number of particles ',N)
    
    G=nx.Graph() # create networkx graph G
    
    for i in range(0,N): #adds a node for each particle in system
        G.add_node(i) # add nodes
        
    for i in range(0,N):
        for j in range(0,N):
            if(i != j): 
                # rsq=(x[i]-x[j])**2+(y[i]-y[j])**2
                if( np.abs(col[i]-col[j]) < 1.01 and np.abs(row[i]-row[j]) < 1.01 ): #add edges if conditions met
                    G.add_edge(i,j) # add edges
                    G.add_edge(j,i) # add edges

    print(list(nx.connected_components(G))) 
    n_connect_comp=nx.number_connected_components(G) #calculates the number of connected components
    print('total number of connected components ',n_connect_comp)

    print('sorting connected components fromn biggest to smallest')
    listy=sorted(nx.connected_components(G), key=len, reverse=True) #sorts the connected components from biggest to smallest

    biggest=len(listy[0]) #finds the biggest component's size
    print('biggest component has size ',biggest) 

    big_threshold=10 #defines a threshold of the number of connected components 
    print('counting connected components with at least ',big_threshold,' particles')

    for i in range(0,n_connect_comp): 
        if( len(listy[i])<big_threshold): #if the threshold is bigger than the coordinate for the biggest component then break
            n_big_comp=i
            break

    print(n_big_comp,' connected components with at least ',big_threshold,' particles')

    plt.scatter(total_array[:,1],total_array[:,0],alpha=0.5,linewidths=0.0) #plot the particles in the system
    #
    plt.xticks([])
    plt.yticks([])
    plt.ylim([box_height+1, -1])
    plt.xlim([-1,box_width+1])
    plt.tight_layout()

    #
    for i in range(0,n_big_comp):
        for ip in listy[i]: plt.scatter(row[ip],col[ip],alpha=0.5,linewidths=0.0,c='red') #display the connected components in red to show them clearly
    #
    plt.xticks([])
    plt.yticks([])
    plt.ylim([box_height+1, -1])
    plt.xlim([-1,box_width+1])
    plt.tight_layout()
    plt.show()

total_array = np.array([]) #create an empty array
length = len(txtfile) 
for i in range(length):
    if 'Matrix' in txtfile[i]: #each matrix contains a title to be read in
        leng = len(total_array)
        total_array = np.reshape(total_array, (int(leng/2), 2)) 
        total_array = total_array.astype(int)
        if leng != 0: #if the array contains some data
            Cluster_res = Cluster(total_array) #call the Cluster Function to find connected components
        total_array = np.array([]) #empty the array 
    else:
        split = txtfile[i].split() #split the text file by each line for each coordinate
        col_coord = int(split[0]) #find the column coordinate (from left to right)
        row_coord = int(split[1]) #find the row coordinate (from top to bottom)
        point = [(col_coord, row_coord)] #define each particle as a point with coordinates
        total_array = np.append(total_array, point) #add the point to the array




