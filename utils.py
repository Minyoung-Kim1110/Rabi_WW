# import libraries 
import numpy as np 
from values import * 

# discretize z 
MODE = 'original'
# MODE = 'shrinked'

d_z_standard = a_0/4 
Z = 20 # range of simulation, - 10 nm <z< 10 nm 
r = 2 

if MODE == 'original':
    d_z = d_z_standard
    
elif MODE == 'shrinked':
    shrink_scale = r
    d_z = d_z_standard/r # shink scale is r, might change 
    
# scale with d_z to conserve the number density 
a_t =np.sqrt(a_0**3/8/d_z)  # nm 

z_values = np.arange(-Z, Z, step = d_z)


# Define the function n(z), Ge profile 
def n(z, w = w, no_wiggle = False):
    if no_wiggle: 
        return (n_bar *(2-1/ (1 + np.exp(4 * (-z - d/2) / w))-1/(1 + np.exp(4 * (z - d/2) / w)))) #+ (2 * n_GE * np.sin(np.pi * z / lbda)**2)
    else:   
        return (n_bar *(2-1/ (1 + np.exp(4 * (-z - d/2) / w))-1/(1 + np.exp(4 * (z - d/2) / w)))) + (2 * n_GE * np.sin(np.pi * z / lbda)**2)

# Define confining potential 
def V_l(w=w, no_wiggle = False):
    return E_G * n(z_values, w, no_wiggle = no_wiggle )

def electric_potential(F_z= F_z):
    return np.arange(len(z_values))*F_z * d_z
    

# Define function that generates potential due to alloy disorders 
def generate_disorders(z_values, n, l_t, E_G, a_t, w=w, no_wiggle = False): 
    # sample disorder terms 
    V_0000 = np.zeros_like(z_values)
    V_1000 = np.zeros_like(z_values)

    Sigma_0000 = 1/(4 * np.pi * l_t**2)
    Sigma_1000 = 0
    Sigma_1100 = 1/(8 * np.pi * l_t**2)
    covariance_V_0000_1000 = lambda z: (E_G * a_t)**2 * n(z,w, no_wiggle = no_wiggle) * (1-n(z,w, no_wiggle = no_wiggle)) * np.array([[Sigma_0000, Sigma_1000],[np.conj(Sigma_1000), Sigma_1100]])
    
    for j in range(len(z_values)):
        z_j = z_values[j]
        random_variables = np.random.multivariate_normal(mean=np.zeros((2,)), cov = covariance_V_0000_1000(z_j))
        V_0000[j]+= random_variables[0]
        V_1000[j]+= random_variables[1]
    return V_0000, V_1000


# out of plane electric field noise spectral density S(f)  = A_z^2 /f 
# A_z is determined by the mesurement. 
# "Quantum logic with spin qubits crossing the surface code threshold"
A_z = 3.5* 10 **(-3) # mV/nm 

# f_max : 1 / time segements 
f_max = 1 / 10 **(-12) # Hz 
# f_min : 1/ measurement time (In Figure 3. T2 measurement is held with 8 min)
f_min = 1 / (8*60) # Hz 

# Variance of Electric field 
var_E_z = A_z**2*np.log(f_max / f_min)
sig_E_z = np.sqrt(var_E_z)

def sample_noise(shape):
    # return : electric field (times e) noise with dimension meV/nm 
    return np.random.normal(loc = 0, scale = sig_E_z, size = shape)


    
def get_groundstate(z_values, well_potential, disorder, electric_noise, electric_potential ):
    length = len(z_values)
    t = - hbar**2 / 2 / m_l / d_z**2 * 10**18

    diag_terms =-2*t * np.ones((length,)) + well_potential+disorder + electric_potential + z_values*electric_noise
    
    H = np.diag(diag_terms)
    H += np.diag(t* np.ones((length-1,)), -1)
    H += np.diag(t* np.ones((length-1,)), -1).T     

    # find ground state 
    _, V = np.linalg.eigh(H)
    ground_state = V[:, 0]
    return ground_state

def expectation_value(term, ground_state):
    val = np.dot(ground_state.conj(), term * ground_state)
    return np.abs(val), np.angle(val) 
    

def calculate_Rabi_freq(Delta_0000, Delta_1000, beta, mode = 'simple'):
    coeff = E * g * mu_b * B / (hbar*(hbar * omega_t)**2)
    if mode == 'simple':
        return coeff * beta[0]*(np.cos(Delta_0000[1]-beta[1]))
    elif mode == 'simple_dipole': # we have to check the calculation 
        return coeff * beta[0]*(np.cos(Delta_0000[1]-beta[1]) - (Delta_1000[0]/Delta_0000[0])**2 * np.sin(beta[1]-Delta_1000[1])*np.sin(Delta_0000[1]-Delta_1000[1]))