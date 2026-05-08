## get fdistribution of chain lengths in network ##########

### For end-linked network: T1_id=1,   T2_id=2 #########
### For side-linked network: T1_id=1,   T2_id=3 ########

import numpy as np
##import ioLAMMPS
##import ioLAMMPS
import math
import matplotlib
##matplotlib.use('Agg') 
from matplotlib import pyplot as plt
import os
import shutil
from numpy import linalg as LA
import scipy.optimize as opt
import param as p
import networkx as nx
import random
import networkx as nx
from scipy.stats import gaussian_kde

from collections import Counter


def readLAMMPS_restart(filename, vflag,frac_weak,G,Gstrong):

   f1=open(filename,"r")

   line1 = f1.readline()
   line2 = f1.readline()

   line3 = f1.readline()
   line3 = line3.strip()
   n_links = int(line3.split(" ")[0])
 
   line4 = f1.readline()
   line4 = line4.strip()
   atom_types = int(line4.split(" ")[0])

   line5 = f1.readline()
   line5 = line5.strip()
   n_chains = int(line5.split(" ")[0])

   line6 = f1.readline()\
           
   line6 = line6.strip()
   bond_types = int(line6.split(" ")[0])

   links_unsort  = np.zeros((n_links,4))
   links   = np.zeros((n_links,3), dtype = float)
   chains  = np.full((n_chains,4), -1, dtype = int)
   mass    = np.zeros(atom_types, dtype = float)

   line7 = f1.readline()
   line8 = f1.readline()
   line8 = line8.strip()
   xlo = float(line8.split(" ")[0])
   xhi = float(line8.split(" ")[1])

   line9 = f1.readline()
   line9 = line9.strip()
   ylo = float(line9.split(" ")[0])
   yhi = float(line9.split(" ")[1])

   line10 = f1.readline()
   line10 = line10.strip()
   zlo = float(line10.split(" ")[0])
   zhi = float(line10.split(" ")[1])


   for i in range (0, 3):
       f1.readline()
   
   for i in range(0, atom_types):
       line = f1.readline()
       line = line.strip()
       mass[i] = float(line.split(" ")[1])

   f1.close()


   links_unsort = np.genfromtxt(filename, usecols=(0,3,4,5), skip_header=18, max_rows=n_links)

   for i in range(0, n_links):
       index = int(links_unsort[i,0])
       links[index-1,:] = links_unsort[i,1:4]

   if(vflag==0):
      data= np.genfromtxt(filename,usecols=(1,2,3), skip_header=17+n_links+3, max_rows=n_chains)
      chains[:,0]=data[:,0]-np.ones(len(chains)) # ctype
      chains[:,1]=np.ones(len(chains)) # column of ones
      chains[:,2:4]=data[:,1:3] # cl1,cl2
   elif(vflag==1):
      data= np.genfromtxt(filename,usecols=(1,2,3), skip_header=17+2*n_links+2*3, max_rows=n_chains)
      chains[:,0]=data[:,0]
      chains[:,1]=np.ones(len(chains))
      chains[:,2]=data[:,1]
      chains[:,3]=data[:,2]
   else:
      print("Invalid Velocity Flag")


   
##   print(chains)
   directory = './'+str(int(100*frac_weak))+'/'
   filename = 'primary_loops'
   file_path = os.path.join(directory, filename)
   if not os.path.isdir(directory):
      os.mkdir(directory)  
   loop_atoms = np.genfromtxt(file_path, usecols=(1), skip_header=0)
   loop_atoms.tolist()
    
   for i in range(0,len(chains[:,0])):
      G.add_edge(chains[i,2],chains[i,3])
        
   for i in range(0,len(chains[:,0])):
        if(chains[i,0]==T2_id-1): ## add the chain only if it is part of the strong chains network
              Gstrong.add_edge(chains[i,2],chains[i,3])

   
   ##

   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms


##import networkx as nx

vflag=0
frac_weak_arr=[0.33]#,0.4,0.5]
T1_id=1
T2_id=2
frac_weak=frac_weak_arr[0]
lams_max=[4]
fig_cnt=6
R=p.N_low ## chain end-to-end distance  p.b*np.sqrt(p.N)# end to end distance

broken_data=np.genfromtxt('ite_failure.txt')
ite_broken=broken_data#[0]
ite_broken_rounded=int(p.wrt_step*int(ite_broken/p.wrt_step))
ite_arr=np.arange(0,ite_broken,max(p.wrt_step,1),dtype='int')##p.wrt_step
ite=ite_broken_rounded-p.wrt_step
step=20
##ite_arr=ite_arr[::step]
#ite_arr=[0]
print('ite_arr',ite_arr)
diameter=np.zeros(len(ite_arr))
plt.figure(1)
plt.figure(2)
cnt=-1
#ite_arr=[0]
for ite in ite_arr:
    cnt=cnt+1
    G=nx.Graph()
    Gstrong=nx.Graph()
    ##G_full=nx.MultiGraph()
    #plt.figure()
####    [xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G_full]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag,G_full)

    ##[xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G_full]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag)
    
    [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
              atom_types, bond_types, mass, loop_atoms] = readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag, frac_weak,G,Gstrong)
              
    #largest_cc = max(nx.connected_components(G_full), key=len)
    #G = G_full.subgraph(largest_cc)
    #num_edges = G.number_of_edges()
    #largest_conn_component_num_edges.append(num_edges)
    
    connected_components = list(nx.connected_components(G))

    # Identify the largest connected component
    largest_cc_nodes = max(connected_components, key=len)

    # Create a subgraph of the largest connected component
    largest_cc_subgraph = G.subgraph(largest_cc_nodes).copy()

    # Calculate the diameter of the largest connected component
    diameter[cnt] = nx.approximation.diameter(largest_cc_subgraph)  ## only of largest component
    #print(diameter[cnt])
    
    connected_components = [Gstrong.subgraph(c) for c in nx.connected_components(Gstrong)]
    #print('len(connected_components)',len(connected_components))

    # Get the number of edges in each connected component
    component_edge_counts = [nx.approximation.diameter(subgraph) for subgraph in connected_components]

    # Count the frequency of edge counts
    edge_distribution = Counter(component_edge_counts)

    # Sort by number of edges
    edges = sorted(edge_distribution.keys())
    frequencies = [edge_distribution[edge] for edge in edges]
    
    
      
    '''
    # Plot the distribution
    plt.figure()#figsize=(8, 6))
    plt.bar(edges, frequencies, color='navajowhite', edgecolor='black')
    plt.xlabel("Diameter", fontsize=12)
    plt.ylabel("Number of Components", fontsize=12)
    plt.title("Distribution of diameters of Connected Component of strong network, ite="+str(ite), fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xscale("log")
    plt.savefig('dia_dist_discrete'+str(ite))
    #unique_dia, counts_dia = np.unique(edge_distribution, return_counts=True)
    
    
    ##data=np.array(component_edge_counts).reshape(-1, 1)
    plt.figure()
    print('component_edge_counts',component_edge_counts)
    kde = gaussian_kde(component_edge_counts)
    x_vals = np.linspace(0, max(component_edge_counts), 1000)
    y_vals = kde(x_vals)#*x_vals
##               sns.kdeplot(total_mass, fill=False)#, color='skyblue')
    plt.plot(x_vals, y_vals,label='ite='+str(ite))
    plt.xscale('log')
    #plt.savefig('dia_dist_cont_kde_'+str(ite))
    '''
    
    plt.figure(1)
    freq = Counter(component_edge_counts)
    x = np.array(sorted(freq.keys()))
    y = np.array([freq[val] for val in x])
    #plt.plot(x, y, marker='o', linestyle='-',label='ite='+str(ite))  # you can customize this
    
    data=component_edge_counts
    log_min = np.log10(min(data))
    log_max = np.log10(max(data))
    bins = np.logspace(log_min, log_max, num=30)  # 50 bins on log scale
    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    plt.plot(bin_centers, counts, marker='o',label='ite='+str(ite))
    #plt.yscale('log')
    
    data=component_edge_counts
    log_min = min(data)
    log_max = max(data)
    bins = np.linspace(log_min, log_max, num=30)  # 50 bins on log scale
    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    #plt.plot(bin_centers, counts, marker='o',label='ite='+str(ite))
    
    
    
    if(ite==0):
      np.savetxt('T2_length_dist_0.txt',np.transpose(np.array([bin_centers,counts])))


    ##print('component_edge_counts',component_edge_counts)
    #kde = gaussian_kde(component_edge_counts)
    #x_vals = np.linspace(0, max(component_edge_counts), 1000)
    #y_vals = kde(x_vals)#*x_vals
##               sns.kdeplot(total_mass, fill=False)#, color='skyblue')
    #plt.plot(x_vals, y_vals,label='ite='+str(ite))
    plt.xscale('log')
    #plt.yscale('log')          
    
    
    
    
####    G=G_full.copy()
    
            
    
            #edge_connectivity.append(nx.edge_connectivity(G))

            

plt.figure(2)
plt.plot(ite_arr,diameter)
plt.ylabel('Network diameter')
plt.xlabel('Iteration')
plt.savefig('dia_largest_component_vs_ite')            
            
            
####
####plt.plot(ite_arr,largest_conn_component_num_edges,'o-')##,label='overall')
####plt.legend()
####plt.ylabel('Number of edges in largest component')
####plt.xlabel('ite')
######plt.ylim([-1,20])
####plt.savefig('largest_comp_num_edges')

plt.figure(1)
plt.legend()
plt.xlabel('Network diameter')
plt.ylabel('freq')
plt.savefig('dia_dist_linegraph_with_ite')



plt.show()
