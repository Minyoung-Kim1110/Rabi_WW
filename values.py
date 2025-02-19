import numpy as np 

# universal values 
hbar = 6.582*10**(-13) # meV s 
beta_0 = 8.2 # meV nm
c = 2.999*10**8 # m/s 
m_e = 0.511*10**9 /c**2 # meV s^2/m^2, mass of electron 
 
mu_b=5.7883818060E-2                    #bohr magneton meV/T

g=2                                     #g factor
B=0.02                                  #B field (T) range [0,1.5] considered in Ben's paper
E=10**(-2)                                   #In plane AC E field magnitude (10^-2 mv/nm or 10^4 V/m)

w_1=g*mu_b*B/2                          #w_1 has units of energy (meV)

omega=2*w_1/hbar                        #drive frequency resonant with zeeman splitting, has units of frequency
omega_t=1/hbar                          #orbital frequency




# in Si/SiGe device
# a_l=0.136*(1E-9)                        #lattice length (longitudinal)
# a_t=0.384*(1E-9)                        #lattice length (transverse)
# n_ge=0.05
# w_l=(0.01/hbar)                         #w_l due to confinement in z direction (10 meV)

E_G = 0.6 *1000 # meV 
 
m_l = 0.91 * m_e # meV s^2/m^2, effective mass of electron in z direction 
m_t = 0.19 * m_e # meV s^2/m^2, effective mass of electron in x,y direction 

a_0 = 5.431*10**(-1) # nm length of unit cell 

k_0 = 0.83* 2 * np.pi/a_0
k_1 = np.pi/1.57 # 1/nm

omega_z = 10 /hbar # 1/s 
# omega_t = 2 /hbar #1/s

l_z = np.sqrt(hbar / 2 /m_l/omega_z)*10**9 # nm 

l_t = np.sqrt(hbar/2/m_t/omega_t)*10**9 # nm 

# print(f"l_z = {l_z} nm\nl_t = {l_t} nm")

# V_0 = 30 # meV ( =E_G * 2 * n_GE )
# sigma_Delta = 0.2 # meV 


# Define constants to simulate wiggle well 
n_bar = 0.3
n_GE = 0.05
w = 1.9 # nm # width of interface 
d = a_0/4*80 # 9 # nm # well width 

lbda = 1.57 # nm # length of Ge oscillation - might change 
F_z = 5 # meV/nm e times electric field in z direction 

G_lambda = 2 * np.pi/lbda # 1/nm 
kappa = 2* k_1 - G_lambda # nm 
E_lambda = hbar**2 * G_lambda**2 / 2 /m_l *10**18 #meV



# V_0_d_E_lambda = E_G * 2 * n_GE / E_lambda 

