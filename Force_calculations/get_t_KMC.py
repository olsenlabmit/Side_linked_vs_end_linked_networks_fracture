###### get the time iteration corresponding to network fracture ########
######## use frac_weak=1.0 for side-linked and frac_weak=0.33 for end-linked ########

import io
import matplotlib
import math
##matplotlib.use('Agg') 
from matplotlib import pyplot as plt
import numpy as np
import sys
import param as p


import numpy as np
import sys
temps = []
frac_weak=0.33  
with io.open(str(int(100*frac_weak))+"/stress", mode="r") as f:
    next(f)
    for line in f:
        temps.append(line.split())

Lx=[float(i[0]) for i in temps]
Ly=[float(i[1]) for i in temps]
Lz=[float(i[2]) for i in temps]

lam=[i[3] for i in temps]
lam=[float(i) for i in lam]


FE=[i[4] for i in temps]
FE=[float(i) for i in FE] #free energy stored in chain
deltaFE=[i[5] for i in temps]
deltaFE=[float(i) for i in deltaFE]

st0=[i[6] for i in temps]
st0=[float(i) for i in st0]
st1=[i[7] for i in temps]
st1=[float(i) for i in st1]
st2=[i[8] for i in temps]
st2=[float(i) for i in st2]
st3=[i[9] for i in temps]
st3=[float(i) for i in st3]
st4=[i[10] for i in temps]
st4=[float(i) for i in st4]
st5=[i[11] for i in temps]
st5=[float(i) for i in st5]
factor=4.11
st0=np.array(st0)*factor
st1=np.array(st1)*factor
st2=np.array(st2)*factor
st3=np.array(st3)*factor
st4=np.array(st4)*factor
st5=np.array(st5)*factor


stress=st0-0.5*(st1+st2)


#########t_KMC=[i[12] for i in temps]
#########t_KMC=[float(i) for i in t_KMC]
########
########
########index_failure=np.where(-stress[3:]<=0.0001)
######### get the t_KMC when the stress becomes less than 0
##########        strain_at_failure_i=t_KMC[index_failure[0]]
##########        strain_at_failure_i=strain_at_failure_i[1]
#########t_KMC=np.array(t_KMC)
##########        stop
#########t_KMC_at_failure_i=t_KMC[index_failure[0]+3]
##########failure_idx.append(index_failure[0][0]+3)
#########t_KMC_at_failure_i=t_KMC_at_failure_i[1]
########
########np.savetxt('ite_failure.txt', np.array([index_failure[0][1]+3]))
##########        replica_data[Run_cnt,frac_weak_cnt]=t_KMC_at_failure_i
##########t_KMC_failure.append(t_KMC_at_failure_i)
##########st6=np.array(st6)



index_failure=np.where(-stress[3:]<=0.0001)
cnt=3
while(index_failure[0][1]<10):
  cnt=cnt+1
  index_failure=np.where(-stress[cnt:]<=0.0001)
# get the t_KMC when the stress becomes less than 0
##        strain_at_failure_i=t_KMC[index_failure[0]]
##        strain_at_failure_i=strain_at_failure_i[1]
#t_KMC=np.array(t_KMC)
##        stop
#t_KMC_at_failure_i=t_KMC[index_failure[0]+3]
##failure_idx.append(index_failure[0][0]+3)
#t_KMC_at_failure_i=t_KMC_at_failure_i[1]

np.savetxt('ite_failure.txt', np.array([index_failure[0][1]+cnt]))

