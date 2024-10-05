"""Example 3
   ==========

   This example shows how to process the FAIR dipole measurements.

"""

import numpy as np
import pandas as pd
import cauchy
import os
import matplotlib.pyplot as plt

# %%
# Some functions
# ==============

def compute_multipoles(x, y, R, mesh, u, g, N=15, N_phi=180):
    '''Compute the multipoles at some position with coordinates x,y, given the
    boundary data.

    :param x:
        The x position.

    :param y:
        The y position.

    :param R:
        The radius.

    :param mesh:
        The boundary mesh.

    :param u:
        The potential.

    :param g:
        The normal derivative.

    :param N:
        The number of multipoles to reconstruct.

    :param N_phi:
        The number of points to evaluate along the azimuth angle.

    :return:
        The multipoles Bn and An.
    '''


    # The return values
    Bn = np.zeros((N, ))
    An = np.zeros((N, ))

    # make the points to evaluate the potential at
    phi = np.linspace(0., 2*np.pi, N_phi, endpoint=False)
    points = np.zeros((N_phi, 2))
    points[:, 0] = R*np.cos(phi) + x
    points[:, 1] = R*np.sin(phi) + y

    # Evaluate the B field
    B = cauchy.grad_representation_formula(mesh, u, g, points)

    # Compute the Br component
    Br = np.cos(phi)*B[:, 0] + np.sin(phi)*B[:, 1]

    # compute the harmonics
    for n in range(1, N+1):
        An[n-1] = 2.0/N_phi*np.sum(Br*np.cos(n*phi))
        Bn[n-1] = 2.0/N_phi*np.sum(Br*np.sin(n*phi))

    return Bn, An

# %%
# Some script parameters
# ======================

# this if condition makes sure that we load the data from the right place
# when launching the autodocs
if os.getcwd().split('\\')[-1] == 'examples':
    meas_dir =  r'..\data\fair_dipole\25A\20241004_135258_1606_D06_1606_D06\20241004_135258_1606_D06_Current_25_A'
else:
    meas_dir = r'data\fair_dipole\25A\20241004_135258_1606_D06_1606_D06\20241004_135258_1606_D06_Current_25_A'    


# %%
# Read the input data
# ===================

# we need to read the mesh information
mesh = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_mesh.csv'), header=None).values*1e-3

# and also the Neumann data (1st column)
u = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MSP_Dirichlet_data.csv'), header=None).values[:, 0]*1e-3

# this is the standard deviation of the Neumann data (2nd column)
u_std = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MSP_Dirichlet_data.csv'), header=None).values[:, 1]*1e-3

# and also the Neumann data (1st column)
g = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MSP_Neumann_data.csv'), header=None).values[:, 0]*1e-3

# this is the standard deviation of the Neumann data (2nd column)
g_std = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MSP_Neumann_data.csv'), header=None).values[:, 1]*1e-3


# %%
# Compute the multipoles along x
# ==============================

x_min = -0.1
x_max = 0.1
num_x = 21

rad = 0.05
x_eval = np.linspace(x_min, x_max, num_x)

An = np.zeros((num_x, 15))
Bn = np.zeros((num_x, 15))
An_std = np.zeros((num_x, 15))
Bn_std = np.zeros((num_x, 15))


for i, xx in enumerate(x_eval):
    Bn[i, :], An[i, :] = compute_multipoles(xx, 0., rad, mesh, u, g)
    Bn_std[i, :], An_std[i, :] = compute_multipoles(xx, 0., rad, mesh, u_std, g_std)

fig = plt.figure()
ax = fig.add_subplot(221)
ax.plot(x_eval, Bn[:, 0], label='normal dipole')
ax.plot(x_eval, An[:, 0], label='skew dipole')
ax.set_xlabel('$x$ in m')
ax.set_title('multipole component in T')
ax.legend()
ax = fig.add_subplot(222)
ax.plot(x_eval, Bn[:, 1]/Bn[:, 0]*1e4, label='normal quadrupole')
ax.plot(x_eval, An[:, 1]/Bn[:, 0]*1e4, label='skew quadrupole')
ax.set_xlabel('$x$ in m')
ax.set_title('multipole component in units')
ax.legend()
ax = fig.add_subplot(223)
ax.plot(x_eval, Bn[:, 2]/Bn[:, 0]*1e4, label='normal sextupole')
ax.plot(x_eval, An[:, 2]/Bn[:, 0]*1e4, label='skew sextupole')
ax.set_xlabel('$x$ in m')
ax.set_title('multipole component in T')
ax.legend()
ax = fig.add_subplot(224)
ax.plot(x_eval, Bn[:, 3]/Bn[:, 0]*1e4, label='normal octupole')
ax.plot(x_eval, An[:, 3]/Bn[:, 0]*1e4, label='skew octupole')
ax.set_xlabel('$x$ in m')
ax.set_title('multipole component in units')
ax.legend()

fig = plt.figure()
ax = fig.add_subplot(221)
ax.plot(x_eval, 3*Bn_std[:, 0], label='normal dipole')
ax.plot(x_eval, 3*An_std[:, 0], label='skew dipole')
ax.set_xlabel('$x$ in m')
ax.set_title('3 $\sigma$ multipole component in T')
ax.legend()
ax = fig.add_subplot(222)
ax.plot(x_eval, 3*Bn_std[:, 1]/Bn[:, 0]*1e4, label='normal quadrupole')
ax.plot(x_eval, 3*An_std[:, 1]/Bn[:, 0]*1e4, label='skew quadrupole')
ax.set_xlabel('$x$ in m')
ax.set_title('3 $\sigma$ multipole component in units')
ax.legend()
ax = fig.add_subplot(223)
ax.plot(x_eval, 3*Bn_std[:, 2]/Bn[:, 0]*1e4, label='normal sextupole')
ax.plot(x_eval, 3*An_std[:, 2]/Bn[:, 0]*1e4, label='skew sextupole')
ax.set_xlabel('$x$ in m')
ax.set_title('3 $\sigma$ multipole component in T')
ax.legend()
ax = fig.add_subplot(224)
ax.plot(x_eval, 3*Bn_std[:, 3]/Bn[:, 0]*1e4, label='normal octupole')
ax.plot(x_eval, 3*An_std[:, 3]/Bn[:, 0]*1e4, label='skew octupole')
ax.set_xlabel('$x$ in m')
ax.set_title('3 $\sigma$ multipole component in units')
ax.legend()
plt.show()
