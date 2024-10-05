"""Example 2
   ==========

   This example shows how to process the single stretched wire
   measurement results. We are working with the standard output
   of FFMM (see https://gitlab.cern.ch/te-msc-mm/ffmm).
   The data is processed using a Dirichlet to Neumann map.
   We also propagate the measurement uncertainty.
   
   See https://ieeexplore.ieee.org/document/8954937

"""

import numpy as np
import pandas as pd
import cauchy
import os
import matplotlib.pyplot as plt

# %%
# Some script parameters
# ======================

# this if condition makes sure that we load the data from the right place
# when launching the autodocs
if os.getcwd().split('\\')[-1] == 'examples':
    meas_dir = os.path.join('..\data', '20230926_104544_HCMBXW_001-BI000027_BEM_900_120_50_5')
else:
    meas_dir = os.path.join('data', '20230926_104544_HCMBXW_001-BI000027_BEM_900_120_50_5')    


# %%
# Read the input data
# ===================

# we need to read the mesh information
mesh = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_mesh.csv'), header=None).values*1e-3

# and also the Neumann data (1st column)
u = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MVP_Dirichlet_data.csv'), header=None).values[:, 0]*1e-3

# this is the standard deviation of the Neumann data (2nd column)
u_std = pd.read_csv(os.path.join(meas_dir, 'BEM_results_rep_1_MVP_Dirichlet_data.csv'), header=None).values[:, 1]*1e-3


# %%
# Launch the Dirichlet 2 Neumann map
# ==================================
g, Sigma_g, A = cauchy.dirichlet_2_neumann(mesh, u, u_std)

# %%
# Evaluate the solution and make a nice field map.
# ================================================

# make evaluation points
x_lim = (min(mesh[:, 0]), max(mesh[:, 0]))
y_lim = (min(mesh[:, 1]), max(mesh[:, 1]))

# this is the resolution
x_resol = 100
y_resol = 20

# we drop the points to not evaluate directly at the boundary
x_line = np.linspace(x_lim[0], x_lim[1], x_resol + 2)[1:-1]
y_line = np.linspace(y_lim[0], y_lim[1], y_resol + 2)[1:-1]

# make a meshgrid
X, Y = np.meshgrid(x_line, y_line)

# and then a numpy (M x 2) array
points = np.zeros((x_resol*y_resol, 2))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

# we now evaluate the interior field as well as the measurement
# uncertainty
B, var_B = cauchy.BEM_uncertainty_quantification(mesh, A, g, u_std, points, formulation='A')

# the arrays B and var_B store the mean and variance of the
# two components Bx and By in the columns. Below we show now
# how to plot the magnitude.

# we compute the field magnitude
B_mag = np.linalg.norm(B, axis=1)
# and the standard deviation of it
std_B_mag = np.sqrt(var_B[:, 0] + var_B[:, 1])

# we reshape to match the sape of the meshgrid
B_mag.shape = X.shape
std_B_mag.shape = X.shape


# %%
# Plot
# ====

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(211)
cntrf = ax.contourf(X, Y, B_mag, levels=100, cmap='gnuplot2')
plt.colorbar(cntrf)
ax.set_xlabel('$x$ in mm')
ax.set_ylabel('$y$ in mm')
ax.set_title('integral field magnitude in Tm')

ax = fig.add_subplot(212)
cntrf = ax.contourf(X, Y, 3*std_B_mag, levels=100, cmap='gnuplot2')
plt.colorbar(cntrf)
ax.set_xlabel('$x$ in mm')
ax.set_ylabel('$y$ in mm')
ax.set_title('3 standard deviations in integral field magnitude in Tm')

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.97, 
                    top=0.96, wspace=0.4,hspace=0.25)

plt.show()
