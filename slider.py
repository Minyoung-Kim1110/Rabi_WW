import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from values import *
from utils import *

well_potential = V_l()

# Create figure and axis
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)

t_init = 0
line, = ax.plot(z_values, np.zeros_like(z_values),'-o', label="Initial")
line2, = ax.plot(z_values, well_potential/10000*3, label="Well")

ax.set_xlim(-10, 10) 
ax.set_ylim(0, 0.06)

# Create slider
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])
sample_size = 500+1

# slider = Slider(ax_slider, 't', 0, sample_size-1, valinit=t_init)
slider = Slider(ax_slider, 'index', 0, 100-1, valinit=t_init)


electric_noise = sample_noise((sample_size,))

# electric_potential_ = electric_potential()
V_0000, V_1000 = generate_disorders(z_values, n, l_t, E_G, a_t)
realization_size = 1
electric_field_sketch = np.linspace(1, 11, num = 100) # mV/nm 


# Update function
def update(val):
    j = slider.val
    j = int(j)
    electric_potential_ = electric_potential(F_z = electric_field_sketch[j])

    electric_noise_ = 0 #electric_noise[j]
    gs = get_groundstate(z_values, well_potential=well_potential, disorder=V_0000, electric_noise=electric_noise_, electric_potential=electric_potential_)    
    line.set_ydata(gs.conj() * gs)
    
    # ax.legend([r"$|\psi|^2,\;F_z$"+f"={5 + electric_noise_:.4f} mV/nm"])
    ax.legend([r"$|\psi|^2,\;F_z$"+f"={electric_field_sketch[j]:.4f} mV/nm"])
    
    fig.canvas.draw_idle()

# Connect slider to update function
slider.on_changed(update)

# plt.legend()
plt.show()
