#!/usr/bin/python

from src.values import *
from src.utils import *
import numpy as np 
from scipy.spatial.distance import cdist
from scipy.linalg import cholesky

a_t_ = a_t*3
N_t = 20
N_z = len(z_values) 
print(f'number of grid in x :{N_t}, z:{N_z}, total:{N_t**2*N_z}')

x = np.arange(N_t)* a_t_
y = np.arange(N_t)* a_t_ 
X, Y = np.meshgrid(x, y)
V_0000 = np.zeros((N_t, N_t, N_z))
V_1000 = np.zeros((N_t, N_t, N_z))
points = np.column_stack([X.ravel(), Y.ravel()]) # linearlized points 

exp_cov = cdist(points, points, 'sqeuclidean')
exp_cov = np.exp(-exp_cov/4/l_t**2)

c1 = lambda z: (E_G*a_t_)**2 * n(z)*(1-n(z))/(4*np.pi*l_t**2)* exp_cov

for i,z in enumerate(z_values):
    C_joint = c1(z)
    mean = np.zeros(C_joint.shape[0])  # Zero mean
    L = cholesky(C_joint + 1e-6 * np.eye(C_joint.shape[0]), lower=True)  # Cholesky decomposition    
    random_samples = L @ np.random.randn(C_joint.shape[0])
    
    V_0000[:, :, i]+= random_samples.reshape(N_t, N_t)
    V_1000[1:-1, :, i]+= l_t * (V_0000[2:, :, i] - V_0000[:-2, :, i])/a_t_



well_potential = V_l()
electric_potential_ = electric_potential()

# Generate noises
sample_size = 200
electric_noise = sample_noise(shape = (sample_size,))

calculated_values = np.zeros((N_t, N_t, sample_size, 4, 2))

for i, x0 in enumerate(x): 
    for j, y0 in enumerate(y): 
        for k in range(sample_size):
            electric_noise_ = electric_noise[k]
        
            gs = get_groundstate(z_values, well_potential=well_potential, disorder=V_0000[i, j, :], electric_noise=electric_noise_, electric_potential=electric_potential_)
            
            Delta_0000 = expectation_value((well_potential+ V_0000[i, j, :]) *np.exp((0-1j)*2*k_0*z_values), ground_state=gs)
            Delta_1000 = expectation_value((well_potential+ V_1000[i, j, :]) *np.exp((0-1j)*2*k_0*z_values), ground_state=gs)

            beta =expectation_value(beta_0*np.exp((0+1j)*2*k_1*z_values), ground_state=gs)
            
            z_loc = expectation_value(z_values, ground_state=gs)
            z2_loc = expectation_value(z_values*z_values, ground_state=gs)
            z_center = z_loc[0]
            z_width = np.sqrt(z2_loc[0] - z_loc[0]**2)
            
            calculated_values[i, j,k, 0, 0], calculated_values[i, j,k, 0, 1] = Delta_0000 
            calculated_values[i, j, k,1, 0], calculated_values[i, j,k, 1, 1] = Delta_1000 
            calculated_values[i, j, k,2, 0], calculated_values[i, j,k, 2, 1] = beta 
            
            calculated_values[i, j, k,3, 0], calculated_values[i, j,k, 3, 1] = z_center, z_width 
         
            
np.save('./Results/Exp_vals', calculated_values)