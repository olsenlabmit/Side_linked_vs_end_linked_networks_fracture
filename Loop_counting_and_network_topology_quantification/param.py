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
frac_weak_arr=[0.33]

re_step=730
####U0_low = 10.0
####U0_low = 10.0
####U0_high = 10
####
####N_low=18.4
####N_high=18.4
####Nmain=18.4
##frac_weak=0.0

##c_dim=5
lam_max=30
del_t=0.002
e_rate=5
##n_chains=7500

b_low=1.0
b_high=1.0
b_main=1.0

K_low=1.0
K_high=1.0
K_main=1.0

####fit_param_low=1.0
####fit_param_high=1.0
####fit_param_main=1.0

E_b_low=1200.0
E_b_high=1200.0
E_b_main=1200.0

func=3

tol=0.01
max_itr=100000
write_itr = 100
wrt_step = 10


n_links=int(n_chains*2/func)
L=((n_chains/c_dim)*N_low**1.5)**(1/3)


########
lam_max=40
