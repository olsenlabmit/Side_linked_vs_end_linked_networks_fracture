c_dim=3
n_chains=7500
U0_low=1
U0_high=10
fit_param_low=0.1

##U0_low=10
##U0_high=U0_low
U0_main=U0_high
N_low=12
N_high=12
Nmain=12
##fit_param_low=1
fit_param_high=1
fit_param_main=1
frac_weak_arr=[1.0]

##c_dim=10


######U0_low = 1.0
######U0_main = 2.0
######U0_high = 5.0  # not used anywhere
######
######N_low=18
######N_high=18
######Nmain=18


######fit_param_low=1.0
######fit_param_high=1.0
######fit_param_main=1.0

##frac_weak=0.0

m=50

re_step=590


lam_max=20
del_t=0.002
e_rate=5
##n_chains=7500

b_low=1.0
b_high=1.0
b_main=1.0

K_low=1.0
K_high=1.0
K_main=1.0



E_b_low=1200.0
E_b_high=1200.0
E_b_main=1200.0

func=3

tol=0.01
max_itr=100000
write_itr = 1000
wrt_step = 10


nc=int(n_chains*2/(3*m)) # number of long polymer chains
print('nc',nc)

L=((n_chains/c_dim)*Nmain**1.5)**(1/3) #######

##########
###
lam_max=40
