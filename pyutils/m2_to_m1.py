import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import io
import scipy
from scipy import optimize

def m2_to_m1(object, Bstruct): 

    m_xx = object[0,:,:,:]
    m_yy = object[1,:,:,:]
    m_zz = object[2,:,:,:]

    s_3D = m_xx + m_yy + m_zz
    s = np.matrix.flatten(s_3D, order = 'F')

    m_xy = object[3,:,:,:]
    m_xz = object[4,:,:,:]
    m_yz = object[5,:,:,:]

    m2 = [np.transpose(np.matrix.flatten(m_xx, order = 'F')), 
          np.transpose(np.matrix.flatten(m_yy, order = 'F')), np.transpose(np.matrix.flatten(m_zz, order = 'F')), 
          np.transpose(np.matrix.flatten(m_xy, order = 'F')), np.transpose(np.matrix.flatten(m_xz, order = 'F')), 
          np.transpose(np.matrix.flatten(m_yz, order = 'F'))]

    m2 = np.concatenate(m2, axis = 0)
    m2 = np.reshape(m2, (np.shape(m2)[0]//6, 6), order = 'F')

    [z, x, y] = np.unravel_index(np.argwhere(s > 1e-2), np.shape(m_xx), order = 'F')

    # retrieve second moments of voxels with localizations
    m2 = np.squeeze(m2[np.argwhere(s > 1e-2),:])

    # Retrieve Bstruct and define b, B, sumNormList
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

    return m1


def secondM2SymmConeWeighted(b, B, sumNorm, secM, signal, backg):
    # define a least squares problem structure

    bXX = b['XX'][0,0]; bYY = b['YY'][0,0]; bZZ = b['ZZ'][0,0]; bXY = b['XY'][0,0]; bXZ = b['XZ'][0,0]; bYZ = b['YZ'][0,0]

    Ix = secM[0] * bXX + secM[1] * bYY + secM[2] * bZZ + secM[3] * bXY + secM[4] * bXZ + secM[5] * bYZ
    Ix[Ix<0] = 0

    # Normalize second moment info and multiply by photon number
    Ix  = Ix / sumNorm; Ix = signal * Ix

    # Compute Fisher Information Matrix
    FIM = calFIMSecondM(B, Ix, signal, backg)


    # construct the M matrix
    M = np.array([[secM[0], secM[3], secM[4]],
        [secM[3], secM[1], secM[5]],
        [secM[4], secM[5], secM[2]]],dtype='float64')
    D,V = np.linalg.eig(M)
    V_real = np.real(V)
    D_real = np.real(D)

    idx_max_D = np.argmax(D_real)

    # initialization via SVD
    x0SVD = np.zeros((4,1))

    # Initial value based the largest eigen value
    #------------------------------------------------------------
    x0SVD[0] = np.real(V[0, idx_max_D])
    x0SVD[1] = np.real(V[1, idx_max_D])
    x0SVD[2] = np.real(V[2, idx_max_D])
    x0SVD[3] = 1.5 * np.max(D_real[idx_max_D]) - .5

    # upper and lower constraints for estimated first moments
    bound = [(-1, 1), (-1, 1), (-1, 1), (0, 1)]

    cons = {'type':'ineq', 'fun': constraint}

    # interior-point optimization
    M1estresults = optimize.minimize(objective_fun, x0SVD, bounds = bound, constraints = cons, args = (FIM, np.matrix(secM),))
    # print(M1estresults)
    estM1 = M1estresults['x']
    # print(estM1)
    # print(np.linalg.norm(estM1[0:2]))

    return estM1

def calFIMSecondM(B, I, s, backg):
    FIM = np.zeros((6, 6))

    Baa = B['aa'][0,0]; Bbb = B['bb'][0,0]; Bcc = B['cc'][0,0]; Bdd = B['dd'][0,0]; Bee = B['ee'][0,0]; Bff = B['ff'][0,0]
    Bab = B['ab'][0,0]; Bac = B['ac'][0,0]; Bad = B['ad'][0,0]; Bae = B['ae'][0,0]; Baf = B['af'][0,0]
    Bbc = B['bc'][0,0]; Bbd = B['bd'][0,0]; Bbe = B['be'][0,0]; Bbf = B['bf'][0,0]
    Bcd = B['cd'][0,0]; Bce = B['ce'][0,0]; Bcf = B['cf'][0,0]
    Bde = B['de'][0,0]; Bdf = B['df'][0,0]
    Bef = B['ef'][0,0]


    FIM[0, 0] = 0.5 * (s**2) * (sum(sum((Baa) / (I + backg))))
    FIM[1, 1] = 0.5 * (s**2) * (sum(sum((Bbb) / (I + backg))))
    FIM[2, 2] = 0.5 * (s**2) * (sum(sum((Bcc) / (I + backg))))
    FIM[3, 3] = 0.5 * (s**2) * (sum(sum((Bdd) / (I + backg))))
    FIM[4, 4] = 0.5 * (s**2) * (sum(sum((Bee) / (I + backg))))
    FIM[5, 5] = 0.5 * (s**2) * (sum(sum((Bff) / (I + backg))))

    FIM[0, 1] = (s**2) * (sum(sum((Bab) / (I + backg))))
    FIM[0, 2] = (s**2) * (sum(sum((Bac) / (I + backg))))
    FIM[0, 3] = (s**2) * (sum(sum((Bad) / (I + backg))))
    FIM[0, 4] = (s**2) * (sum(sum((Bae) / (I + backg))))
    FIM[0, 5] = (s**2) * (sum(sum((Baf) / (I + backg))))

    FIM[1, 2] = (s**2) * (sum(sum((Bbc) / (I + backg))))
    FIM[1, 3] = (s**2) * (sum(sum((Bbd) / (I + backg))))
    FIM[1, 4] = (s**2) * (sum(sum((Bbe) / (I + backg))))
    FIM[1, 5] = (s**2) * (sum(sum((Bbf) / (I + backg))))

    FIM[2, 3] = (s**2) * (sum(sum((Bcd) / (I + backg))))
    FIM[2, 4] = (s**2) * (sum(sum((Bce) / (I + backg))))
    FIM[2, 5] = (s**2) * (sum(sum((Bcf) / (I + backg))))

    FIM[3, 4] = (s**2) * (sum(sum((Bde) / (I + backg))))
    FIM[3, 5] = (s**2) * (sum(sum((Bdf) / (I + backg))))

    FIM[4, 5] = (s**2) * (sum(sum((Bef) / (I + backg))))

    FIM = FIM + np.transpose(FIM)

    return FIM

def symmCone2SecM(z):
    z = np.reshape(z, (1, 4), order='F')
    muxx = z[:,3] * z[:,0]**2 + (1 - z[:,3]) / 3
    muyy = z[:,3] * z[:,1]**2 + (1 - z[:,3]) / 3
    muzz = z[:,3] * z[:,2]**2 + (1 - z[:,3]) / 3
    muxy = z[:,3] * z[:,0] * z[:,1]
    muxz = z[:,3] * z[:,0] * z[:,2]
    muyz = z[:,3] * z[:,1] * z[:,2]
    out = np.array([[muxx[0]], [muyy[0]], [muzz[0]], [muxy[0]], [muxz[0]], [muyz[0]]])
    return out


def objective_fun(z, FIM, secM):
    '''Objective function, test lambda objective function first'''
    sC2SM = symmCone2SecM(z)
    return np.matmul(np.matmul((np.transpose(sC2SM) - secM), FIM), (sC2SM - np.transpose(secM)))[0,0]

def constraint(z):
    '''Constraint function'''
    return np.sum(z[0:3]**2, 0) - 1

if __name__ == "__main__":
    # import estimated object second moment information
    obj_dir = 'performance test objects'
    obj = scipy.io.loadmat(os.path.join(obj_dir, 'hemisphere gamma 1.mat'))['object']
    # import MVR basis images
    B_dir = 'psf'
    B_fname = os.path.join(B_dir, 'B_ZScanned_zf1000_PSFsz_61.mat')
    Bfile = scipy.io.loadmat(B_fname)
    Bstruct = Bfile['Bstruct']

    est_m1 = m2_to_m1(obj, Bstruct)

    a = 1