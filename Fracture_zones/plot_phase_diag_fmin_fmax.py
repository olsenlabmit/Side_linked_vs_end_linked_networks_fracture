import numpy as np
from  matplotlib import pyplot as plt
import os
import re

##plt.figure()
U0_2=10 ## _arr=[1]##,2,5]##[5,10,20,40]
d_2=1.0
##U0_1_arr=[50, 40, 30, 20, 15, 10, 7, 5, 4, 3, 2, 1.8, 1.5, 1.3, 1.2, 1.1, 1, 0.8, 0.5, 0.3, 0.1, 0, -0.1, -0.5, -1, -2, -5, -10]##[0,-0.1,-0.5,-1,-2,-5,-10,2,0.3,5,10,15,20,30,40,50,0.5,0.8,1.3,1.8,1.5,1,1.1,0,0.1,0.3,0.8,1.2,1.5,2,0.8,1.2,1.5,2,3,4,5,7,10]
c_dim=3

##plt.rcParams["font.family"] = "Arial"



tough_arr=[]
tough_arr_norm=[]
tough_err_arr=[]
max_ext_arr=[]
max_ext_err_arr=[]
fmax_arr=[]
fmin_arr=[]

fmin_line_arr=[]
fmax_line_arr=[]
U0_arr=[]
delta_arr=[]


## fracture regimes from theory  and  fmin line  #################

plt.figure(1,figsize=(7.5, 5.5))
plt.axvline(x=0,color='k',linestyle='--')
plt.axhline(y=0,color='k',linestyle='--')

plt.figure(2,figsize=(7.5, 5.5)) ## fmin line
plt.title('fmin line')

plt.figure(3,figsize=(7.5, 5.5)) ## fmin line
plt.title('fmax line')


plt.figure(4,figsize=(7.5, 5.5)) ## fmin line
plt.title('fmax line')

##pattern = match = re.match(r"^a_([-]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)_b_([-]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)$", folder)
##r'U0weak_(-?\d+(?:\.\d+)?)_d_(-?\d+(?:\.\d+)?)'##r'a_(-?\d+(?:\.\d+)?)_b_(-?\d+(?:\.\d+)?)'
top_folders = [f for f in os.listdir('.') if os.path.isdir(f)]

for folder in top_folders:
  match = re.match(r"^U0weak_([-]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)_d_weak_([-]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)$", folder)

  if match:
      U0_1 = float(match.group(1))
      d_1 = float(match.group(2))
      
      #print(f"In folder '{folder}': x = {U0_1}, y = {d_1}")
  else:
      
      #print(f"Folder name '{folder}' doesn't match pattern 'a_x_b_y'")
      
      continue

  try:
     
     filename=folder+"/toughness.txt"
     
     toughness=np.genfromtxt(filename)
     filename=folder+"/fmin_avg.txt"
     fmin=np.genfromtxt(filename)[0]
     filename=folder+"/fmax_avg.txt"
     fmax=np.genfromtxt(filename)[0]
     filename=folder+"/max_ext_avg.txt"
     max_ext=np.genfromtxt(filename)  
     
    

     max_ext_arr.append(max_ext[0])
     tough_arr_norm.append(toughness[0]/fmax)
     tough_arr.append(toughness[0])
     tough_err_arr.append(toughness[1])
     ##a_arr.append(((12+U0_1)/(12+U0_2))*(d_2/(d_1+small_num)))
     fmax_arr.append(fmax)
     fmin_arr.append(fmin)
     U0_arr.append(U0_1)
     delta_arr.append(d_1)
     
  except Exception as e:
      print(e) 
      continue    

#U0_arr=set(U0_arr)
#delta_arr=set(delta_arr)

#U0_arr=np.array(U0_arr)
#delta_arr=np.array(delta_arr)
#'''
combined_arr=np.transpose(np.array([U0_arr,delta_arr, fmin_arr]))
column_to_sort = combined_arr[:, 0] ## sort by U0 value
sorted_indices = column_to_sort.argsort()
sorted_arr=combined_arr[sorted_indices[::]]
U0_arr=sorted_arr[:, 0]
delta_arr=sorted_arr[:, 1]
fmin_arr=sorted_arr[:, 2]
#'''
print('U0_arr',U0_arr)
print('delta_arr',delta_arr)
print('sorted delta_arr',np.sort(delta_arr))



'''

############ U-delta phase diagram########
plt.figure(figsize=(7.5, 5.5))
print(min(tough_arr),max(tough_arr))
sc = plt.scatter(U0_arr, delta_arr, c=fmin_arr, cmap='inferno', edgecolors='k', vmin=0, vmax=max(fmin_arr))
plt.axvline(x=U0_2,color='r',linestyle='--')
plt.axhline(y=d_2,color='r',linestyle='--')    
cbar = plt.colorbar(sc)
cbar.set_label('Toughness')
plt.xlabel('U0_1',fontsize=16)#, fontname='Helvetica')
plt.ylabel('delta_1',fontsize=16)#, fontname='Helvetica')
#plt.xlim([min(U0_arr)*1.1,max(U0_arr)*1.1])
#plt.ylim([min(delta_arr)*1.1,max(delta_arr)*1.1]) 
  
plt.savefig('U-delta_phase_diag')
'''


from matplotlib.colors import LinearSegmentedColormap
# Define colors
burgundy = '#800020'    # dark red/burgundy
light_cyan = '#E0FFFF'  # light cyan/aqua
teal = '#008080'        # teal

# Create colormap passing through light cyan
custom_cmap = LinearSegmentedColormap.from_list('burgundy_cyan_teal',
                                                [burgundy,teal,light_cyan])
                                                
                                                

################## delta_U-delta_d phase diagram #####################

plt.figure(figsize=(8.7, 7.0))

delta_U0_arr=np.array(U0_arr)-U0_2
delta_d_arr=np.array(delta_arr)-d_2
plt.axvline(x=0,color='k',linestyle='--', zorder=1)
plt.axhline(y=0,color='k',linestyle='--', zorder=1) 

sc = plt.scatter(delta_U0_arr, delta_d_arr,s=145, c=fmin_arr, marker='s',  cmap=custom_cmap, vmin=0, vmax=max(fmin_arr), edgecolors='k',linewidths=0.0)

plt.xlabel('$\Delta$$U_0$',fontsize=16, fontweight='bold')#, fontname='Arial')
plt.ylabel('$\Delta$$\delta$',fontsize=16, fontweight='bold')#, fontname='Helvetica')

#plt.xlim([min(delta_U0_arr)*1.1,max(delta_U0_arr)*1.1])
#plt.ylim([min(delta_d_arr)*1.1,max(delta_d_arr)*1.1])    
plt.xticks(fontsize=14)#, fontname='Helvetica')
plt.yticks(fontsize=14)#, fontname='Helvetica')
cbar = plt.colorbar(sc)
cbar.set_label('fmin')

plt.savefig('fmin_delta_U-delta_d_phase_diag.pdf', dpi=600, transparent=True)


combined_arr=np.transpose(np.array([U0_arr,delta_arr, fmax_arr]))
column_to_sort = combined_arr[:, 0] ## sort by U0 value
sorted_indices = column_to_sort.argsort()
sorted_arr=combined_arr[sorted_indices[::]]
U0_arr=sorted_arr[:, 0]
delta_arr=sorted_arr[:, 1]
fmax_arr=sorted_arr[:, 2]
#'''
print('U0_arr',U0_arr)
print('delta_arr',delta_arr)
print('sorted delta_arr',np.sort(delta_arr))



'''

############ U-delta phase diagram########
plt.figure(figsize=(7.5, 5.5))
print(min(tough_arr),max(tough_arr))
sc = plt.scatter(U0_arr, delta_arr, c=fmin_arr, cmap='inferno', edgecolors='k', vmin=0, vmax=max(fmin_arr))
plt.axvline(x=U0_2,color='r',linestyle='--')
plt.axhline(y=d_2,color='r',linestyle='--')    
cbar = plt.colorbar(sc)
cbar.set_label('Toughness')
plt.xlabel('U0_1',fontsize=16)#, fontname='Helvetica')
plt.ylabel('delta_1',fontsize=16)#, fontname='Helvetica')
#plt.xlim([min(U0_arr)*1.1,max(U0_arr)*1.1])
#plt.ylim([min(delta_arr)*1.1,max(delta_arr)*1.1]) 
  
plt.savefig('U-delta_phase_diag')
'''

################## delta_U-delta_d phase diagram #####################

plt.figure(figsize=(8.7, 7.0))

delta_U0_arr=np.array(U0_arr)-U0_2
delta_d_arr=np.array(delta_arr)-d_2
plt.axvline(x=0,color='k',linestyle='--', zorder=1)
plt.axhline(y=0,color='k',linestyle='--', zorder=1) 

sc = plt.scatter(delta_U0_arr, delta_d_arr,s=145, c=fmax_arr, marker='s',  cmap=custom_cmap, vmin=0, vmax=max(fmax_arr), edgecolors='k',linewidths=0.0)

plt.xlabel('$\Delta$$U_0$',fontsize=16, fontweight='bold')#, fontname='Arial')
plt.ylabel('$\Delta$$\delta$',fontsize=16, fontweight='bold')#, fontname='Helvetica')

#plt.xlim([min(delta_U0_arr)*1.1,max(delta_U0_arr)*1.1])
#plt.ylim([min(delta_d_arr)*1.1,max(delta_d_arr)*1.1])    
plt.xticks(fontsize=14)#, fontname='Helvetica')
plt.yticks(fontsize=14)#, fontname='Helvetica')
cbar = plt.colorbar(sc)
cbar.set_label('fmax')

plt.savefig('fmax_delta_U-delta_d_phase_diag.pdf', dpi=600, transparent=True)


plt.show()

