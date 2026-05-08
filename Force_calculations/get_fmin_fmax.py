## get the minimum force needed to break chains and the maxinum force in chains in network at failure

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
def chain_length_dist(n_chains, chains, links, Lx, Ly, Lz):
##        print(chains[:,1:])
        chain_lengths=[]
##        dist = np.zeros((n_chains,4))
##        dist[:,0:3] = chains[:,1:]
##        dist[:,3] = -1
            
        for i in range (0, n_chains):
            if(chains[i,2] !=-1):
          
              link_1 = chains[i,2]-1
              link_2 = chains[i,3]-1
              lk = links[link_1,:] - links[link_2,:]
              
              lk[0] = lk[0] - int(round(lk[0]/Lx))*Lx
              lk[1] = lk[1] - int(round(lk[1]/Ly))*Ly
              lk[2] = lk[2] - int(round(lk[2]/Lz))*Lz
                    
##              dist[i,3] = LA.norm(lk)
              chain_lengths.append(LA.norm(lk))#dist[i,3])
##              meanr2=meanr2+(dist[i,3])**2
    ##          print(dist[i,3])
    ##          stop
    ##          print(((dist[i,3])**2)/(p.N_low*p.b_low**2))
              

        return np.array(chain_lengths)

def chain_length_dist_x(n_chains, chains, links, Lx, Ly, Lz):
##        print(chains[:,1:])
        chain_lengths=[]
##        dist = np.zeros((n_chains,4))
##        dist[:,0:3] = chains[:,1:]
##        dist[:,3] = -1
            
        for i in range (0, n_chains):
            if(chains[i,2] !=-1):
          
              link_1 = chains[i,2]-1
              link_2 = chains[i,3]-1
              lk = links[link_1,:] - links[link_2,:]
              
              lk[0] = lk[0] - int(round(lk[0]/Lx))*Lx
              lk[1] = lk[1] - int(round(lk[1]/Ly))*Ly
              lk[2] = lk[2] - int(round(lk[2]/Lz))*Lz
                    
##              dist[i,3] = LA.norm(lk)
              chain_lengths.append(abs(lk[0]))#dist[i,3])
##              meanr2=meanr2+(dist[i,3])**2
    ##          print(dist[i,3])
    ##          stop
    ##          print(((dist[i,3])**2)/(p.N_low*p.b_low**2))
              

        return np.array(chain_lengths)



def force_dist(chain_lengths):
##        print(chains[:,1:])
        force=[]
##        dist = np.zeros((n_chains,4))
##        dist[:,0:3] = chains[:,1:]
##        dist[:,3] = -1
            
        for i in chain_lengths:
            force.append(get_bondforce(i))#dist[i,3])          

        return np.array(force)
##    

def invlangevin(x):
        return x*(2.99942 - 2.57332*x + 0.654805*x**2)/(1-0.894936*x - 0.105064*x**2)

def kuhn_stretch(lam, E_b):
        def func(x, lam, E_b):
            y = lam/x
            beta = invlangevin(y)
            return E_b*np.log(x) - lam*beta/x

        if lam == 0:
           return 1
        else:
           lam_b = opt.root_scalar(func,args=(lam, E_b),bracket=[lam,lam+1],x0=lam+0.05)
           return lam_b.root

        
def get_bondforce(r):

        #K  = p.K
        r0 = 0
        Nb = p.N_low # b = 1 (lenght scale of the system)
        E_b = p.E_b_low
 
        x = (r-r0)/Nb
        if(x<0.90):
           lam_b = 1.0
           fbkT  = invlangevin(x)
##           fbond = -K*fbkT/r
        elif(x<1.4):
           lam_b = kuhn_stretch(x, E_b)
           fbkT  = invlangevin(x/lam_b)/lam_b
##           fbond = -K*fbkT/r
        else:
           lam_b = x + 0.05
           fbkT  = 325 + 400*(x-1.4)
##           stop
##           fbond = -K*fbkT/r
 
        return fbkT

def readLAMMPS_restart(filename, vflag,G):

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


##   chains[:,0] = N
#cnt,ctype,1,conn1,conn2
   if(vflag==0):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(0,1,2,3), skip_header=17+n_links+3, max_rows=n_chains)
   elif(vflag==1):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(0,1,2,3), skip_header=17+2*n_links+2*3, max_rows=n_chains)
   else:
      print("Invalid Velocity Flag")
##   print(chains)

####   for c in chains:
####      [lnk_1,lnk_2]=c[2:4]
####      G.add_edge(lnk_1,lnk_2)
   directory = './'+str(int(100*frac_weak))+'/'
   filename = 'primary_loops'
   file_path = os.path.join(directory, filename)
   #print(file_path)
######   if not os.path.isdir(directory):
######      os.mkdir(directory)  
   loop_atoms = np.genfromtxt(file_path , usecols=(1), skip_header=0)
   loop_atoms.tolist() 

   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms,G




import networkx as nx

vflag=0
frac_weak_arr=[0.33]#,0.4,0.5]
frac_weak=frac_weak_arr[0]
lams_max=[4]
fig_cnt=6
R=p.N_low ## chain end-to-end distance  p.b*np.sqrt(p.N)# end to end distance

broken_data=np.genfromtxt('ite_failure.txt')
ite_broken=broken_data#[0]
ite_broken_rounded=int(p.wrt_step*int(ite_broken/p.wrt_step))
ite_arr=np.arange(0,ite_broken,max(p.wrt_step,1),dtype='int')##p.wrt_step
ite=ite_broken_rounded-p.wrt_step

fmax=0
fmax_arr=[]
fmean_arr=[]
for ite in ite_arr:
    [xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag,frac_weak)
    Lx=xhi-xlo
    Ly=(yhi-ylo)  # to make sure that the distances are not being calculated across the PBC
    Lz=(zhi-zlo)
    dist=chain_length_dist(n_chains, chains, links, Lx, Ly, Lz)
    dist=force_dist(dist)
    fmax_this_ite=np.max(dist) ## maximum force at breaking ite
    fmax=max(fmax,fmax_this_ite)
    fmax_arr.append(fmax_this_ite)
    fmean_arr.append(np.mean(dist))
[xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite_arr[0])+".txt", vflag,frac_weak)
n_chains_0=n_chains
for ite in ite_arr[1:]:
    [xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, G]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite)+".txt", vflag,frac_weak)
    if(n_chains<n_chains_0):
        print('start of break at ite=',ite)
        [xlo_prev, xhi_prev, ylo_prev, yhi_prev, zlo_prev, zhi_prev, n_links_prev, n_chains_prev, links_prev, chains_prev, atom_types_prev, bond_types_prev, mass_prev, loop_atoms_prev, G_prev]=readLAMMPS_restart("./"+str(round(frac_weak*100))+"/restart_network_"+str(ite-p.wrt_step)+".txt", vflag,frac_weak)
        ##stop
        break
    else:
        n_chains_0=n_chains
## now ite=where breaking starts

n_chains_broken=-n_chains+n_chains_prev  ## this need not be equal to p.wrt_Step= because here- I am not doing adaptive time stepping with one break per step
chains_broken=np.zeros((n_chains_broken,4),dtype='int')
cnt=0
found_cnt=0
i_not_found_arr=[i for i in range(n_chains)]

#broken_chain_positions_lnk1=[]
#broken_chain_positions_lnk2=[]

broken_chains_lengths=[]
chains_1=chains_prev
chains_2=chains
for c in chains_1:
      lnk_1=c[2]
      lnk_2=c[3]
      
      found=False

      for i in i_not_found_arr:##range(n_chains_final):
   ##      if(i not in i_found_arr):
            if((chains_2[i,2]==lnk_1 and chains_2[i,3]==lnk_2)):# or (chains_final[i,2]==lnk_2 and chains_final[i,3]==lnk_1)): # means that chain is present
               found=True
   ##            print(lnk_1,lnk_2)
               found_cnt=found_cnt+1
               i_not_found_arr.remove(i)
               break
      if(found==False): # chain is not found in the final array
            
            
            chains_broken[cnt,:]=c
            cnt=cnt+1
   ##         break
Lx=xhi_prev-xlo_prev
Ly=(yhi_prev-ylo_prev)  # to make sure that the distances are not being calculated across the PBC
Lz=(zhi_prev-zlo_prev)
broken_chains_lengths= chain_length_dist(n_chains_broken, chains_broken, links_prev, Lx, Ly, Lz)
print('number of broken chains',n_chains_broken)
#'''
dist=chain_length_dist(n_chains, chains, links, Lx, Ly, Lz)
dist=force_dist(dist)
fmin1=np.max(dist) ## taking the mean - because the minimum force will be zero, which doesn't make sense for my calculations
#'''
dist=broken_chains_lengths#chain_length_dist(n_chains, chains, links, Lx, Ly, Lz)
dist=force_dist(dist)
fmin=np.mean(dist)

print('fmin1',fmin1,'fmin',fmin,'fmax',fmax)## fmin defined as the max force at the ite just before start of break
np.savetxt('fmin_max.txt',np.array([fmin,fmax]))
plt.plot(ite_arr,fmax_arr,'o-',label='max')
plt.plot(ite_arr,fmean_arr,'o-',label='mean')
plt.legend()

plt.ylabel('fmax')
plt.xlabel('ite')
plt.figure()
plt.plot(ite_arr,fmean_arr,'o-',label='mean')