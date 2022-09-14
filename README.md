# Kawasaki_Ising
Masters Thesis Project : Computational Modelling of the Formation of Salt Crystals Through Evaporation in Droplets Containing COVID-19

4 files are contained in this repository.

KW_Ising_Part1.py
KW_Ising_Part1.py is the initial file to use. It creates an empty lattice which is populated, and runs the Ising model at the selected temperature. 
The evapoation rate can be sped up or slowed down entirely depending on the requirements.
The energy is calculated periodically throughout the simulation
At the end, intermediatary versions of the lattice are populated into 3 .txt files 
which can be run by the other python files to quantify the outputs from the Ising model

KW_Ising_Part2.py
KW_Ising_Part1.py is the second file to use. It plots the lattices at intermediate stages throughout the simulation to give a visual indication 
of how the particles are moving within the structure
The energy and peclet number are also available as outputs in this file

Density_Plot.py
Density_Plot.py is the third file to use. It plots a profile of the density of particles in each row. 
The file gives a good indication of where the particles are within the box, and shows the effects of evaporation on the system visually
A fit is provided to the graph indicating some parameters from the simulation

network.py
network.py is the final file to use. It quantifys the clusters within the system by using the NetworkX package
the file determines the nunber of particles directly adjacent in a cluster, and can calculate the number of clusters,
the location of these clusters respective to the box dimensions or other clusters,
and the size of these clusters. The file can do this for each iteration of the lattice, so it is possible to quantify how these
clusters change over the course of the simulation. 
The file outputs several graphs to this effect showing mean height of clusters and mean size of clusters over the simulation,
as well as the number of clusters at each size for the final output. 
