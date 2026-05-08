### Calculate depercolation threshold without any chain preference #####
### written: Devosmita Sen #########

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

def readLAMMPS_restart(filename, vflag):

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

   G=nx.MultiGraph()  
##   chains[:,0] = N
#cnt,ctype,conn1,conn2
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
   directory = './'+str(int(100*frac_weak))+"/"
   filename = 'primary_loops'
   file_path = os.path.join(directory, filename)
   if not os.path.isdir(directory):
      os.mkdir(directory)  
   loop_atoms = np.genfromtxt(file_path, usecols=(1), skip_header=0)
   loop_atoms.tolist()
##   print(chains)

   for c in chains:
      ##if(c[1]==1): ## means a weak chain
      G.add_edge(c[2],c[3], weight=c[1])
####      print(c[0])
   #stop
   '''

   G=nx.MultiGraph()
   for i in range(n_chains): #[n1,n2]
         n1=chains[i,2]
         n2=chains[i,3]
         ctype=chains[i,0]
##         print("ctype",ctype)
         if (ctype==0):
            G.add_edge(n1,n2)
    '''

   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G

directory_path='./deperc_no_chain_pref'
os.makedirs(directory_path, exist_ok=True)

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
ite_arr=[0]#np.arange(0,ite_broken,max(p.wrt_step,1),dtype='int')##p.wrt_step
#ite=ite_broken_rounded-p.wrt_step

largest_conn_component_num_edges=[]
deperc_thesh_arr_arr=[]

list_dgel_avg_wo_largest=[]
for ite in ite_arr:
    G_full=nx.MultiGraph()
    plt.figure()
####    [xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G_full]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag,G_full)

    [xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G_full]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag)
    largest_cc = max(nx.connected_components(G_full), key=len)
    G = G_full.subgraph(largest_cc)
    num_edges = G.number_of_edges()
    largest_conn_component_num_edges.append(num_edges)
####    G=G_full.copy()
    if(ite%20==0):
            print(ite)
            
    
            #edge_connectivity.append(nx.edge_connectivity(G))

            G_orig=G.copy()

            potential_node_1_list=[]

            for i in G:
                  if(len(list(G.neighbors(i)))>0):
                     potential_node_1_list.append(i)


            density=[]
            max_connected=[]
            second_largest_connected=[]
            avg_cluster_size_without_largest=[]
            avg_cluster_size_with_largest=[]
            path_exists_array=[]
            num_paths_array=[]
            fraction_paths_connected=[]
            
            frac_cleaved_array=np.arange(0.1,0.8,0.01)#0.01)##[0.4,0.5,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.8,0.9]
            for frac_cleaved in frac_cleaved_array: 
               G_temp=G.copy()
              
               init_num_bonds=G_temp.number_of_edges()
               #print('init_num_bonds',init_num_bonds)
               #print('int(frac_cleaved*init_num_bonds)',int(frac_cleaved*init_num_bonds))
##               all_edges=
               target_weight = 0

               # Get edges with the desired weight
               edges_with_target_weight = list(G_temp.edges())#[edge for edge in G_temp.edges(data=True)]# if edge[2].get('weight') == target_weight]
               #print(edges_with_target_weight)
               #stop



               num_broken_bonds=0
               while(num_broken_bonds<int(frac_cleaved*init_num_bonds)):
####                  for node_1 in potential_node_1_list:
####                     if(num_broken_bonds<int(frac_cleaved*init_num_bonds)):
####                        a=list(nx.all_neighbors(G_temp, node_1))
####                        if(a!=[]):
####                           node_2_list=list(G_temp.neighbors(node_1))
####                           node_2=random.choice(node_2_list)
####                           G_temp.remove_edge(node_1,node_2)
####
####                           num_broken_bonds=num_broken_bonds+1

                                    # Randomly select an edge
                  if len(edges_with_target_weight)>0:
                      selected_edge = random.choice(edges_with_target_weight)
                      #idx = random.randrange(len(edges_with_target_weight))  # or np.random.randint(len(arr))
                      #selected_edge = edges_with_target_weight[idx]
                      #print('selected_edge',selected_edge)
                      #print('selected_edge[0]',selected_edge[0])
                      #edges_with_target_weight = np.delete(edges_with_target_weight, idx)
                      
                      edges_with_target_weight.remove(selected_edge)
                      #print(selected_edge)
####                      print("Selected Edge:", selected_edge)
                      G_temp.remove_edge(selected_edge[0], selected_edge[1])  ##(node_1,node_2)
                      num_broken_bonds=num_broken_bonds+1
                  else:
                      print("No edges with the specified weight found (no weak chains remaining)")
                      print('frac_cleaved',frac_cleaved)
                      break
                  
 
               print('num_broken_bonds',num_broken_bonds)
               density.append(nx.density(G_temp))
               #print('num_broken_bonds',num_broken_bonds)
               
               max_connected.append(max(len(cc) for cc in nx.connected_components(G_temp)))
               cc=list(nx.connected_components(G_temp))
               #print('number of cc',len(cc))
               
               cc.sort(key=len)
               second_largest_connected.append(len(cc[len(cc)-2]))
               all_cluster_sizes=[len(x) for x in cc[0:-1]] # without largest
               cluster_sizes, number_distribution=np.unique(all_cluster_sizes,return_counts=True)

               numerator=0
               denominator=0
               for i in range(len(cluster_sizes)):
                   numerator=numerator+(cluster_sizes[i]**2)*number_distribution[i]
                   denominator=denominator+(cluster_sizes[i])*number_distribution[i]
                   
               avg_cluster_size_without_largest.append(numerator/denominator)#np.mean([len(x) for x in cc[0:-1]]))
               
               all_cluster_sizes=[len(x) for x in cc[0:]] # with largest
               cluster_sizes, number_distribution=np.unique(all_cluster_sizes,return_counts=True)

               numerator=0
               denominator=0
               for i in range(len(cluster_sizes)):
                   numerator=numerator+(cluster_sizes[i]**2)*number_distribution[i]
                   denominator=denominator+(cluster_sizes[i])*number_distribution[i]
                   
               avg_cluster_size_with_largest.append(numerator/denominator)#np.mean([len(x) for x in cc[0:-1]]))
               

               sys_size=len(G_temp)

            
####            i=np.argmax(second_largest_connected)
####            list_dgel_second_largest.append(frac_cleaved_array[i])
####            #print('De-gel point from second largest connected',frac_cleaved_array[i])
            
            i=np.argmax(avg_cluster_size_without_largest)
            list_dgel_avg_wo_largest.append(frac_cleaved_array[i])
            print('De-gel point from average cluster size without largest',frac_cleaved_array[i])
    plt.figure(1)
    plt.plot(frac_cleaved_array,avg_cluster_size_without_largest,'o-',label=str(ite))
    plt.legend()
    plt.savefig(directory_path+'/reduced_cluster_size_plot_'+str(ite))
    np.savetxt(directory_path+'/reduced_cluster_size_ite_'+str(ite)+'.txt',np.transpose(np.array([frac_cleaved_array,avg_cluster_size_without_largest])))
    plt.figure(2)
    plt.plot(frac_cleaved_array,avg_cluster_size_with_largest,'o-',label=str(ite))
    plt.legend()
    plt.savefig(directory_path+'/cluster_size_plot_'+str(ite))
    np.savetxt(directory_path+'/cluster_size_ite_'+str(ite)+'.txt',np.transpose(np.array([frac_cleaved_array,avg_cluster_size_with_largest])))
    plt.close()        
    
    
####
####plt.plot(ite_arr,largest_conn_component_num_edges,'o-')##,label='overall')
####plt.legend()
####plt.ylabel('Number of edges in largest component')
####plt.xlabel('ite')
######plt.ylim([-1,20])
####plt.savefig('largest_comp_num_edges')

plt.figure()
plt.plot(ite_arr[::2],list_dgel_avg_wo_largest,'o-')##,label='overall')
plt.legend()
plt.ylabel('Depercolation threshold')
plt.xlabel('ite')
##plt.ylim([-1,20])
plt.savefig(directory_path+'/deperc_thresh')
np.savetxt(directory_path+'/deperc_thresh.txt',np.transpose(np.array([ite_arr[::2],list_dgel_avg_wo_largest])))



plt.show()
