import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import os

theta_error = loadmat(os.path.join('figure 2','theta_error.mat'))['theta_error']
theta = np.linspace(np.pi/2/18,np.pi/2-np.pi/2/18,num=9,endpoint=True)

# ax = plt.subplot(111,polar=True)
# ax.bar(theta, theta_error[:,0],width=.15,bottom=0.0,alpha=.5,edgecolor='k')
# plt.show()

phi_error = loadmat(os.path.join('figure 2','phi_error.mat'))['phi_error']
phi = np.linspace(-np.pi+np.pi/18,np.pi-np.pi/18,num=18,endpoint=True)

ax = plt.subplot(111,polar=True)
ax.bar(phi, phi_error[:,0],width=.3,bottom=0.0,alpha=.5,edgecolor='k')
plt.show()


