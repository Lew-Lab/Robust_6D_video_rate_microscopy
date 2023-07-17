import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy import io
import os
from secondM2SymmConeWeighted import *
# import torch 

# import estimated object second moment information
obj_dir = 'performance test objects'
obj_fname = os.path.join(obj_dir, 'fiber gamma 1.mat')
print(obj_fname)
obj_file = io.loadmat(obj_fname)
object = obj_file['object']

m_xx = np.squeeze(object[0,:,:,:]);
m_yy = np.squeeze(object[1,:,:,:]);
m_zz = np.squeeze(object[2,:,:,:]);

s = m_xx + m_yy + m_zz; s = np.matrix.flatten(s, order = 'F')

m_xy = np.squeeze(object[3,:,:,:]);
m_xz = np.squeeze(object[4,:,:,:]);
m_yz = np.squeeze(object[5,:,:,:]);

m2 = [np.matrix.flatten(m_xx, order = 'F').T, np.matrix.flatten(m_yy, order = 'F').T, np.matrix.flatten(m_zz, order = 'F').T,
      np.matrix.flatten(m_xy, order = 'F').T, np.matrix.flatten(m_xz, order = 'F').T, np.matrix.flatten(m_yz, order = 'F').T]

m2 = np.concatenate(m2, axis = 0); m2 = np.reshape(m2, (np.shape(m2)[0]//6, 6), order = 'F')

# retrieve voxels with localizations
[z, x, y] = np.unravel_index(np.argwhere(s > 1e-2), np.shape(m_xx), order = 'F')
# print(np.argwhere(s > 1e-2))

# retrieve second moments of voxels with localizations
m2 = np.squeeze(m2[np.argwhere(s > 1e-2),:]); print(len(m2))


# import MVR basis images
B_dir = 'psf'
B_fname = os.path.join(B_dir, 'B_ZScanned_zf1000_PSFsz_61.mat')
# print(B_fname)
Bfile = scipy.io.loadmat(B_fname)

# Retrieve Bstruct and define b, B, sumNormList
Bstruct = Bfile['Bstruct'] # retrieves Bstruct from the B .mat file
BsList = Bstruct['BsList'][0,0] # retrieves BsList from Bstruct
BHaList = Bstruct['B_HaList'][0,0] # retrieves hadamard products from Bstruct
sumNormList = Bstruct['sumNormList'][0,0] # retrieves sumNormList

m1 = np.zeros((len(m2),4)); # initialize array to store estimated first moment information

for loc in range(len(m2)):
    m2Tmp = m2[loc, :] # second moments associated with localization
    # print(m2Tmp[0:3])
    signal = np.sum(m2Tmp[0:3]); # where mTmp = s*(all second moments). Sum the first three (s*mii) to recover s
    # print(signal)
    secM = m2Tmp[:] # normalized second moments associated with localization
    # print(secM)
    zSlice = int(z[loc])
    # print(loc)

    # retrieve b, B, and sumNorm for the correct z-height
    b = BsList[0, zSlice]
    B = BHaList[0, zSlice]
    sumNorm = sumNormList[0, zSlice][0,0]

    # estimate and store first moments
    estM1 = secondM2SymmConeWeighted(b, B, sumNorm, secM/signal, signal, 1e-12)
    estM1[2]= np.abs(estM1[2])
    m1[loc,:] = estM1

# Retrieve angles from first moments
mux = m1[:,0]
muy = m1[:,1]
muz = m1[:,2]
gamma = m1[:,3]

print(gamma)

# theta and phi in degrees
theta = np.degrees(np.arccos(muz))
# print(theta)
phi = np.degrees(np.arctan(muy/mux))
# print(phi)

# convert rotational mobility (gamma) to Omega
# gamma = 1-(3*Omega)/4*pi + Omega^2/(8*pi^2)
# quadratic equation: Omega^2/(8*pi^2) - (3*Omega)/4*pi + (1-gamma)
omega = np.zeros((len(gamma),2))

for loc in range(len(gamma)):
     omegaTmp = np.roots([1/(8*np.pi**2), -3/(4*np.pi), (1-gamma[loc])])
     omega[loc,:] = omegaTmp
# print(omega)

phiNormed = plt.Normalize(-90, 90)
phicolors = plt.cm.hsv(phiNormed(phi))

thetaNormed = plt.Normalize(0, 90)
thetacolors = plt.cm.cool(thetaNormed(theta))
# plot phi and theta for test object
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
a1 = ax1.scatter(x, y, z, c=phicolors, marker='o', alpha=0.4)
fig1.colorbar(plt.cm.ScalarMappable(norm=phiNormed, cmap='hsv'), ax=ax1)


ax1.set_xlabel('X voxel')
ax1.set_ylabel('Y voxel')
ax1.set_zlabel('Z voxel')


fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
a2 = ax2.scatter(x, y, z, c=thetacolors, marker='o', alpha=0.4)
fig2.colorbar(plt.cm.ScalarMappable(norm=thetaNormed, cmap='cool'), ax=ax2)

ax2.set_xlabel('X voxel')
ax2.set_ylabel('Y voxel')
ax2.set_zlabel('Z voxel')


fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
a3 = ax3.scatter(x, y, c=phicolors, marker='o', alpha=0.4)
fig3.colorbar(plt.cm.ScalarMappable(norm=phiNormed, cmap='hsv'), ax=ax3)

ax3.set_xlabel('X voxel')
ax3.set_ylabel('Y voxel')


fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
a4 = ax4.scatter(x, z, c=phicolors, marker='o', alpha=0.4)
fig4.colorbar(plt.cm.ScalarMappable(norm=phiNormed, cmap='hsv'), ax=ax4)

ax4.set_xlabel('X voxel')
ax4.set_ylabel('Z voxel')


fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
a5 = ax5.scatter(y, z, c=phicolors, marker='o', alpha=0.4)
fig5.colorbar(plt.cm.ScalarMappable(norm=phiNormed, cmap='hsv'), ax=ax5)

ax5.set_xlabel('Y voxel')
ax5.set_ylabel('Z voxel')


fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
a6 = ax6.scatter(x, y, c=thetacolors, marker='o', alpha=0.4)
fig6.colorbar(plt.cm.ScalarMappable(norm=thetaNormed, cmap='cool'), ax=ax6)

ax6.set_xlabel('X voxel')
ax6.set_ylabel('Y voxel')


fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
a7 = ax7.scatter(x, z, c=thetacolors, marker='o', alpha=0.4)
fig7.colorbar(plt.cm.ScalarMappable(norm=thetaNormed, cmap='cool'), ax=ax7)

ax7.set_xlabel('X voxel')
ax7.set_ylabel('Z voxel')


fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
a8 = ax8.scatter(y, z, c=thetacolors, marker='o', alpha=0.4)
fig8.colorbar(plt.cm.ScalarMappable(norm=thetaNormed, cmap='cool'), ax=ax8)

ax8.set_xlabel('Y voxel')
ax8.set_ylabel('Z voxel')
plt.show()
