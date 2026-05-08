## get force statistics of network over time


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
frac_weak_arr=[0.33]##[1.0]#,0.4,0.5]
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

fmax=0
fmax_arr=[]
fmean_arr=[]

fmax_T1=0
fmax_T1_arr=[]
fmean_T1_arr=[]

fmax_T2=0
fmax_T2_arr=[]
fmean_T2_arr=[]

order_sum_x_T1_arr=[]
order_sum_y_T1_arr=[]
order_sum_z_T1_arr=[]

order_sum_x_T2_arr=[]
order_sum_y_T2_arr=[]
order_sum_z_T2_arr=[]


order_sum_x_overall_arr=[]
order_sum_y_overall_arr=[]
order_sum_z_overall_arr=[]

n_broken_T1=[]
n_broken_T2=[]

n_chains_T1_arr=[]
n_chains_T2_arr=[]

n_T1_prev_ite=0
n_T2_prev_ite=0

scission_prob_T1=[]
scission_prob_T2=[]

scission_prob_T1_err=[]
scission_prob_T2_err=[]


scission_rate_T1=[]
scission_rate_T2=[]

scission_rate_T1_err=[]
scission_rate_T2_err=[]

r1_r2_ratio=[]

rmax_arr=[]



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
##    stop


    chains_T1=chains[np.where(chains[:,1]==T1_id)[0]]
    n_chains_T1=len(chains_T1)
    dist=chain_length_dist(n_chains_T1, chains_T1, links, Lx, Ly, Lz)
    dist=force_dist(dist)
    fmax_T1_this_ite=np.max(dist) ## maximum force at breaking ite
    fmax_T1=max(fmax_T1,fmax_T1_this_ite)
    fmax_T1_arr.append(fmax_T1_this_ite)
    fmean_T1_arr.append(np.mean(dist))

    rates_T1=np.exp(-p.U0_low+dist*p.fit_param_low)
    

    n_broken_T1.append(-n_chains_T1+n_T1_prev_ite)
    n_T1_prev_ite=n_chains_T1
    n_chains_T1_arr.append(n_chains_T1)

    

    chains_T2=chains[np.where(chains[:,1]==T2_id)[0]]
    n_chains_T2=len(chains_T2)
    dist=chain_length_dist(n_chains_T2, chains_T2, links, Lx, Ly, Lz)
    dist=force_dist(dist)
    fmax_T2_this_ite=np.max(dist) ## maximum force at breaking ite
    fmax_T2=max(fmax_T2,fmax_T2_this_ite)
    fmax_T2_arr.append(fmax_T2_this_ite)
    fmean_T2_arr.append(np.mean(dist))
    rates_T2=np.exp(-p.U0_high+dist*p.fit_param_high)

    rmax=max(max(rates_T1),max(rates_T2))
    rmax_arr.append(rmax)

    prob_T1=rates_T1/rmax
    prob_T2=rates_T2/rmax

    scission_prob_T1.append(np.mean(prob_T1))
    scission_prob_T2.append(np.mean(prob_T2))

    scission_prob_T1_err.append(np.std(prob_T1))
    scission_prob_T2_err.append(np.std(prob_T2))

    scission_rate_T1.append(np.mean(rates_T1))
    scission_rate_T2.append(np.mean(rates_T2))

    scission_rate_T1_err.append(np.std(rates_T1))
    scission_rate_T2_err.append(np.std(rates_T2))

    r1_r2_ratio.append(np.mean(rates_T1)/np.mean(rates_T2))
    
    

    n_broken_T2.append(-n_chains_T2+n_T2_prev_ite)
    n_T2_prev_ite=n_chains_T2
    n_chains_T2_arr.append(n_chains_T2)

    ## Allignment - measured by order parameter= (1/2)<3cos^2(theta)-1>

    order_sum_x_T1=0
    order_sum_y_T1=0
    order_sum_z_T1=0

    order_sum_x_T2=0
    order_sum_y_T2=0
    order_sum_z_T2=0
    for chain in chains_T1:
            link_1 = chain[2]-1
            link_2 = chain[3]-1
            lk = links[link_1,:] - links[link_2,:]

            norm=LA.norm(lk)
            if(norm>0):
                    ## wrt x axis:
                    cos_theta=lk[0]/norm
                    order_sum_x_T1=order_sum_x_T1+3*cos_theta**2-1

                    ## wrt y axis:
                    cos_theta=lk[1]/norm
                    order_sum_y_T1=order_sum_y_T1+3*cos_theta**2-1

                    ## wrt z axis:
                    cos_theta=lk[2]/norm
                    order_sum_z_T1=order_sum_z_T1+3*cos_theta**2-1

    for chain in chains_T2:
            link_1 = chain[2]-1
            link_2 = chain[3]-1
            lk = links[link_1,:] - links[link_2,:]

            norm=LA.norm(lk)
            if(norm>0):
                    ## wrt x axis:
                    cos_theta=lk[0]/norm
                    order_sum_x_T2=order_sum_x_T2+3*cos_theta**2-1

                    ## wrt y axis:
                    cos_theta=lk[1]/norm
                    order_sum_y_T2=order_sum_y_T2+3*cos_theta**2-1

                    ## wrt z axis:
                    cos_theta=lk[2]/norm
                    order_sum_z_T2=order_sum_z_T2+3*cos_theta**2-1

    order_param_T1_x=0.5*order_sum_x_T1/n_chains_T1
    order_param_T2_x=0.5*order_sum_x_T2/n_chains_T2
    order_param_overall_x=0.5*(order_sum_x_T1+order_sum_x_T2)/(n_chains_T1+n_chains_T2)
    
    order_sum_x_T1_arr.append(order_param_T1_x)
    order_sum_x_T2_arr.append(order_param_T2_x)
    order_sum_x_overall_arr.append(order_param_overall_x)

    order_param_T1_y=0.5*order_sum_y_T1/n_chains_T1
    order_param_T2_y=0.5*order_sum_y_T2/n_chains_T2
    order_param_overall_y=0.5*(order_sum_y_T1+order_sum_y_T2)/(n_chains_T1+n_chains_T2)
    
    order_sum_y_T1_arr.append(order_param_T1_y)
    order_sum_y_T2_arr.append(order_param_T2_y)
    order_sum_y_overall_arr.append(order_param_overall_y)

    order_param_T1_z=0.5*order_sum_z_T1/n_chains_T1
    order_param_T2_z=0.5*order_sum_z_T2/n_chains_T2
    order_param_overall_z=0.5*(order_sum_z_T1+order_sum_z_T2)/(n_chains_T1+n_chains_T2)
    
    order_sum_z_T1_arr.append(order_param_T1_z)
    order_sum_z_T2_arr.append(order_param_T2_z)
    order_sum_z_overall_arr.append(order_param_overall_z)
    


####print('fmin1',fmin1,'fmin',fmin,'fmax',fmax)## fmin defined as the mean force at the ite just before start of break (actually st the start of the break)
####np.savetxt('fmin_max.txt',np.array([fmin,fmax]))
plt.plot(ite_arr,fmax_arr,'o-',label='overall')
plt.plot(ite_arr,fmax_T1_arr,'o-',label='T1')
plt.plot(ite_arr,fmax_T2_arr,'o-',label='T2')

plt.legend()
plt.ylabel('f max')
plt.xlabel('ite')
plt.ylim([-1,20])
plt.savefig('fstats_max')

plt.figure()
plt.errorbar(ite_arr,scission_prob_T1,yerr=scission_prob_T1_err,fmt='o-',label='T1',capsize=3)#, alpha=0.1)
plt.errorbar(ite_arr,scission_prob_T2,yerr=scission_prob_T2_err,fmt='o-',label='T2',capsize=3)#, alpha=0.1)
plt.legend()
plt.ylabel('Scission probability')
plt.xlabel('ite')
plt.savefig('scission_prob')

plt.figure()
plt.errorbar(ite_arr,scission_rate_T1,yerr=scission_rate_T1_err,fmt='o-',label='T1',capsize=3)#, alpha=0.1)
plt.errorbar(ite_arr,scission_rate_T2,yerr=scission_rate_T2_err,fmt='o-',label='T2',capsize=3)#, alpha=0.1)
plt.legend()
plt.ylabel('Scission rate')
plt.xlabel('ite')
plt.savefig('scission_rate')


plt.figure()
plt.plot(ite_arr,rmax_arr,'o-',label='T1')
##plt.errorbar(ite_arr,scission_prob_T2,yerr=scission_prob_T2_err,fmt='o-',label='T2', alpha=0.1)
plt.legend()
plt.ylabel('r_max')
plt.xlabel('ite')
plt.savefig('r_max')

plt.figure()
plt.plot(ite_arr,r1_r2_ratio,'o-')#,label='T1')
##plt.errorbar(ite_arr,scission_prob_T2,yerr=scission_prob_T2_err,fmt='o-',label='T2', alpha=0.1)
plt.legend()
plt.ylabel('r1_r2_ratio')
plt.xlabel('ite')
plt.savefig('r1_r2_ratio')


np.savetxt('r1_r2_ratio_data.txt',np.transpose(np.array([ite_arr,r1_r2_ratio])),header='ite,r1_r2_ratio')


plt.figure()
plt.plot(ite_arr,fmean_arr,'o-',label='overall')
plt.plot(ite_arr,fmean_T1_arr,'o-',label='T1')
plt.plot(ite_arr,fmean_T2_arr,'o-',label='T2')
plt.legend()
plt.ylabel('f mean')
plt.ylim([-0.1,1.6])
plt.xlabel('ite')
##plt.title('')
plt.savefig('fstats_mean')



np.savetxt('f_max_data.txt',np.transpose(np.array([ite_arr,fmax_T1_arr,fmax_T2_arr,fmax_arr])),header='ite,T1,T2,overall')
np.savetxt('f_mean_data.txt',np.transpose(np.array([ite_arr,fmean_T1_arr,fmean_T2_arr,fmean_arr])),header='ite,T1,T2,overall')

####plt.figure()
####
####
####plt.legend()
####plt.ylabel('f')
####plt.xlabel('ite')
####plt.title('T2')
####plt.savefig('fstats_T2')

plt.figure()
plt.title('x')
plt.plot(ite_arr,order_sum_x_T1_arr,'o-',label='T1_x')
plt.plot(ite_arr,order_sum_x_T2_arr,'o-',label='T2_x')
plt.plot(ite_arr,order_sum_x_overall_arr,'o-',label='overall_x')
plt.legend()
plt.ylim([-0.05,0.9])
plt.ylabel('order parameter')
plt.xlabel('ite')
plt.savefig('order_parameter_x')


plt.figure()
plt.title('y')
plt.plot(ite_arr,order_sum_y_T1_arr,'o-',label='T1_y')
plt.plot(ite_arr,order_sum_y_T2_arr,'o-',label='T2_y')
plt.plot(ite_arr,order_sum_y_overall_arr,'o-',label='overall_y')
plt.legend()
plt.ylim([-0.55,0.1])
plt.ylabel('order_parameter')
plt.xlabel('ite')
plt.savefig('order_parameter_y')

plt.figure()
plt.title('z')
plt.plot(ite_arr,order_sum_z_T1_arr,'o-',label='T1_z')
plt.plot(ite_arr,order_sum_z_T2_arr,'o-',label='T2_z')
plt.plot(ite_arr,order_sum_z_overall_arr,'o-',label='overall_z')
plt.legend()
plt.ylim([-0.55,0.1])
plt.ylabel('order_parameter')
plt.xlabel('ite')
plt.savefig('order_parameter_z')

np.savetxt('order_param_data_x.txt',np.transpose(np.array([ite_arr,order_sum_x_T1_arr,order_sum_x_T2_arr,order_sum_x_overall_arr])),header='ite,T1,T2,overall')
np.savetxt('order_param_data_y.txt',np.transpose(np.array([ite_arr,order_sum_y_T1_arr,order_sum_y_T2_arr,order_sum_y_overall_arr])),header='ite,T1,T2,overall')
np.savetxt('order_param_data_z.txt',np.transpose(np.array([ite_arr,order_sum_z_T1_arr,order_sum_z_T2_arr,order_sum_z_overall_arr])),header='ite,T1,T2,overall')


plt.figure()
plt.ylabel('number of broken chains in this iteration')
plt.xlabel('ite')
plt.plot(ite_arr[1:],n_broken_T1[1:],'o-',label='T1')
plt.plot(ite_arr[1:],n_broken_T2[1:],'o-',label='T2')
plt.plot(ite_arr[1:],np.array(n_broken_T1[1:])+np.array(n_broken_T2[1:]),'o-',label='overall')
plt.legend()
plt.ylim([-1,35])
plt.savefig('n_broken')
np.savetxt('n_broken_data.txt',np.transpose(np.array([ite_arr[1:],n_broken_T1[1:],n_broken_T2[1:]])),header='ite,T1,T2')


plt.figure()
plt.ylabel('fraction of chains broken in this iteration')
plt.xlabel('ite')
plt.plot(ite_arr[1:],np.array(n_broken_T1[1:])/np.array(n_chains_T1_arr[1:]),'o-',label='T1')
plt.plot(ite_arr[1:],np.array(n_broken_T2[1:])/np.array(n_chains_T2_arr[1:]),'o-',label='T2')
plt.plot(ite_arr[1:],(np.array(n_broken_T1[1:])+np.array(n_broken_T2[1:]))/(np.array(n_chains_T1_arr[1:])+np.array(n_chains_T2_arr[1:])),'o-',label='overall')
plt.legend()
plt.ylim([-0.001,0.016])
plt.savefig('frac_broken')
np.savetxt('frac_broken_data.txt',np.transpose(np.array([ite_arr[1:],n_broken_T1[1:]/np.array(n_chains_T1_arr[1:]),n_broken_T2[1:]/np.array(n_chains_T1_arr[1:])])),header='ite,T1,T2')

plt.figure()
plt.ylabel('number of existing chains')
plt.xlabel('ite')
plt.plot(ite_arr[0:],n_chains_T1_arr[0:],'o-',label='T1')
plt.plot(ite_arr[0:],n_chains_T2_arr[0:],'o-',label='T2')
plt.plot(ite_arr[0:],np.array(n_chains_T1_arr[0:])+np.array(n_chains_T2_arr[0:]),'o-',label='overall')
plt.legend()
plt.savefig('n_existing_chains')
####plt.figure()
####plt.plot(ite_arr,fmean_arr,'o-',label='mean')
####plt.plot(ite_arr,fmax_arr,'o-',label='max')


plt.show()
