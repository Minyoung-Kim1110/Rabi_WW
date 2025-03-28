from src.values import *
from src.utils import *


N_t = np.floor(140/ a_t).astype(np.int64)
N_z = len(z_values) 
print(f'number of grid in x :{N_t}, z:{N_z}, total:{N_t**2*N_z}')


germanium_profile = n(z_values)
germanium_profile =np.tile(germanium_profile, (N_t, N_t, 1)) 
Ge = np.random.uniform(size = (N_t, N_t, N_z))
Ge = Ge<germanium_profile

potential_disorder = E_G * (Ge *(1-germanium_profile) - (1-Ge) *  germanium_profile)  

del germanium_profile
del Ge 
 
x = np.arange(N_t)* a_t 
y = np.arange(N_t)* a_t 
chi_0000 = lambda x, y, x0, y0 : a_t**2*(1/(2*np.pi*l_t**2)) * np.exp(- ((x-x0)**2 + (y-y0)**2)/2/l_t**2)
chi_1000= lambda x, y, x0, y0 : a_t**2*(1/(2*np.pi*l_t**2)) * np.exp(- ((x-x0)**2 + (y-y0)**2)/2/l_t**2)*(x-x0)/l_t



cut = 80 # to save memory 

well_potential = V_l()
electric_potential_ = electric_potential()

Delta_map = np.zeros((N_t-2*cut, N_t-2*cut, 4))
Rabi_vals = np.zeros((N_t-2*cut, N_t-2*cut, 2)) 

for i, x0 in enumerate(x[cut:-cut]): 
    for j, y0 in enumerate(y[cut:-cut]): 
        
        X, Y = np.meshgrid(x[i: i+2*cut], y[j: j+2*cut])
        
        V_0000 = np.sum(potential_disorder[i:i+2*cut, j:j+2*cut, :]*(chi_0000(X, Y, x0, y0).reshape(2*cut, 2*cut, 1) @ np.ones((1,N_z))), axis = (0, 1))
        V_1000 = np.sum(potential_disorder[i:i+2*cut, j:j+2*cut, :]*(chi_1000(X, Y, x0, y0).reshape(2*cut, 2*cut, 1) @ np.ones((1,N_z))), axis = (0, 1))
        
        gs = get_groundstate(z_values, well_potential=well_potential, disorder=V_0000, electric_noise=0, electric_potential=electric_potential_)
        
        Delta_0000 = expectation_value((well_potential+ V_0000) *np.exp((0-1j)*2*k_0*z_values), ground_state=gs)
        Delta_1000 = expectation_value((well_potential+ V_1000) *np.exp((0-1j)*2*k_0*z_values), ground_state=gs)
        
        Delta_map[i, j, 0]+= Delta_0000[0]
        Delta_map[i, j, 1]+= Delta_0000[1]
        
        Delta_map[i, j, 2]+= Delta_1000[0]
        Delta_map[i, j, 3]+= Delta_1000[1]
        
        
        beta =expectation_value(beta_0*np.exp((0+1j)*2*k_1*z_values), ground_state=gs)
        
        Omega_0 =calculate_Rabi_freq(Delta_0000, Delta_1000, beta, mode='simple')/10**9/B #GHz/T 
        Omega_dipole=calculate_Rabi_freq(Delta_0000, Delta_1000, beta, mode='simple_dipole')/10**9/B #GHz/T
        Rabi_vals[i, j, 0]+= Omega_0
        Rabi_vals[i, j, 1]+=Omega_dipole
        
np.save(f"./Results/{timestamp}_3D_Rabi", Rabi_vals) 
np.save(f"./Results/{timestamp}_3D_Delta", Delta_map) 
