
import numpy as np
import scipy
import os
from secondM2SymmConeWeighted import *
# import torch 

# import estimated object second moment information
obj_dir = 'objects'
obj_fname = os.path.join(obj_dir, 'sphere 2.mat')
print(obj_fname)
obj_file = scipy.io.loadmat(obj_fname)
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
print(np.argwhere(s > 1e-2))
# retrieve second moments of voxels with localizations
m2 = np.squeeze(m2[np.argwhere(s > 1e-2),:]); print(len(m2))


# import MVR basis images
B_dir = 'psf'
B_fname = os.path.join(B_dir, 'B_ZScanned_zf1000_PSFsz_61.mat')
print(B_fname)
Bfile = scipy.io.loadmat(B_fname)

# Retrieve Bstruct and define b, B, sumNormList
Bstruct = Bfile['Bstruct'] # retrieves Bstruct from the B .mat file
BsList = Bstruct['BsList'][0,0] # retrieves BsList from Bstruct ([0,0] is necessary for some reason)
BHaList = Bstruct['B_HaList'][0,0] # retrieves hadamard products from Bstruct
sumNormList = Bstruct['sumNormList'][0,0] # retrieves sumNormList

m1 = np.zeros((len(m2),4));
b = BsList[0,0] # retrieves XX, YY, ZZ, XY, XZ, YZ basis images from Bslist for specific z-height. columns are indexes through b at each z-height
B = BHaList[0,0] # retrieves hadamard products of basis images for specific z-height. columns are indexes through B at each z-height

for loc in range(len(m2)):
    m2Tmp = m2[loc, :] # second moments associated with localization
    signal = np.sum(m2Tmp[0:2]); # where mTmp = s*(all second moments). Sum the first three (s*mii) to recover s
    secM = m2Tmp[:] # normalized second moments associated with localization
    zSlice = int(z[loc])
    b = BsList[0, zSlice]
    B = BHaList[0, zSlice]
    sumNorm = sumNormList[0, zSlice]
    estM1 = secondM2SymmConeWeighted(b, B, sumNorm, secM/signal, signal, 1e-12)
    estM1[2] = np.abs(estM1[2])
    m1[loc,0:3] = estM1

    

