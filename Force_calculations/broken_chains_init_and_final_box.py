###### Map broken chains in a specific iteration to the initial simulation box #########
##### by default: the time iteration chosen is the point at which the network fails, which is calculaeted usifn get_t_KMC.py code ##


### For end-linked network: T1_id=1,   T2_id=2 #########
### For side-linked network: T1_id=1,   T2_id=3 ########


import numpy as np
##

import math
import matplotlib
##matplotlib.use('Agg') 
from matplotlib import pyplot as plt
import param as p
import os
import shutil
import sys
sys.path.append('./function_files/')
import ioLAMMPS

##directory = './function_files/'
##orig_dir = os.path.dirname(directory)
##files=os.listdir(orig_dir)
##
##
####for fname in files:
####   if(fname=='ioLAMMPS.py'):
#### # copying the files to the
#### # destination directory
####       shutil.copy2(os.path.join(orig_dir,fname), directory)
##
##import function_files/ioLAMMPS.py as ioLAMMPS

def readLAMMPS_restart(filename, vflag,frac_weak):

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

   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms






vflag=0
frac_weak=0.33

broken_data=np.genfromtxt('ite_failure.txt')
ite_broken=broken_data#[1]
ite_broken_rounded=int(p.wrt_step*int(ite_broken/p.wrt_step))
ite_init=0
#transition_ite=int((p.lam_max-1)/(p.e_rate*p.del_t))
#ite_init_step=transition_ite  # init ite for just step strain sim

f1='./'+str(int(100*frac_weak))+'/'+'restart_network_'+str(ite_init)+'.txt'  # init box
f2='./'+str(int(100*frac_weak))+'/'+'restart_network_'+str(ite_broken_rounded)+'.txt'   # final box
#f3='./'+str(int(100*frac_weak))+'/'+'restart_network_'+str(ite_init_step)+'.txt'   # init box for step strain (end of continuous sim)



[xlo_final, xhi_final, ylo_final, yhi_final, zlo_final, zhi_final, n_links, n_chains_final, links_final, chains_final, atom_types, bond_types, mass, loop_atoms_final]=readLAMMPS_restart(f2, vflag,frac_weak)
[xlo_init, xhi_init, ylo_init, yhi_init, zlo_init, zhi_init, n_links, n_chains_init, links_init, chains_init, atom_types, bond_types, mass, loop_atoms]=readLAMMPS_restart(f1, vflag,frac_weak)

##[xlo_init_step, xhi_init_step, ylo_init_step, yhi_init_step, zlo_init_step, zhi_init_step, n_links_step, n_chains_init_step, links_init_step, chains_init_step, atom_types_step, bond_types_step, mass_step, loop_atoms_step]=readLAMMPS_restart(f3, vflag,frac_weak)

n_chains_broken=-n_chains_final+n_chains_init
chains_broken=np.zeros((n_chains_broken,4),dtype='int')
print('n_chains_broken- overall', n_chains_broken)
cnt=0
found_cnt=0
i_not_found_arr=[i for i in range(n_chains_final)]
for c in chains_init:
   lnk_1=c[2]
   lnk_2=c[3]
   found=False

   for i in i_not_found_arr:##range(n_chains_final):
##      if(i not in i_found_arr):
         if((chains_final[i,2]==lnk_1 and chains_final[i,3]==lnk_2)):# or (chains_final[i,2]==lnk_2 and chains_final[i,3]==lnk_1)): # means that chain is present
            found=True
##            print(lnk_1,lnk_2)
            found_cnt=found_cnt+1
            i_not_found_arr.remove(i)
            break
   if(found==False): # chain is not found in the final array
         
         chains_broken[cnt,:]=c
         cnt=cnt+1
##         break
         
         
ioLAMMPS.writeLAMMPS('./'+str(int(100*frac_weak))+'/'+'broken_chains_init_box.txt',xlo_init, xhi_init, ylo_init, yhi_init, zlo_init, zhi_init, links_init,chains_broken, atom_types, bond_types, mass, loop_atoms_final)

##ioLAMMPS.writeLAMMPS('./'+str(int(100*frac_weak))+'/'+'broken_chains_final_box.txt',xlo_final, xhi_final, ylo_final, yhi_final, zlo_final, zhi_final, links_final,chains_broken, atom_types, bond_types, mass, loop_atoms_final)
##
##ioLAMMPS.writeLAMMPS('./'+str(int(100*frac_weak))+'/'+'broken_chains_init_step_box_all.txt',xlo_init_step, xhi_init_step, ylo_init_step, yhi_init_step, zlo_init_step, zhi_init_step, links_init_step,chains_broken, atom_types, bond_types, mass, loop_atoms_final)
