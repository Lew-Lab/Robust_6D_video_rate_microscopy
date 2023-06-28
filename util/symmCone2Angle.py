
import numpy as np
import scipy
import os
# import torch 

# directory for PSFs 
data_dir = 'psf'
B_fname = os.path.join(data_dir, 'B_ZScanned_zf1000_PSFsz_61.mat')
print(B_fname)
Bfile = scipy.io.loadmat(B_fname)

Bstruct = Bfile['Bstruct']
b = Bstruct['BsList']
B = Bstruct['B_HaList']
sumNormList = Bstruct['sumNormList']
print(np.size(B[0]))

def secondM2SymmConeWeighted(b, B, sumNorm, secM, signal, backg):
    # define a least squares problem structure

    Ix = secM[0] * b.XX + secM[1] * b.YY + secM[2] * (b.ZZ) + secM[3] * (b.XY) + secM[4] * (b.XZ) + secM[5] * (b.YZ)
    Ix[Ix<0] = 0

    Ix  = Ix / sumNorm; Ix = signal * Ix

    FIM = calFIMSecondM(B, Ix, signal, backg)

    # construct M matrix

    # Initial value based the largest eigen value
    #------------------------------------------------------------
    # construct the M matrix
    M = [[secM[0], secM[3], secM[4]],
        [secM[3], secM[1], secM[5]],
        [secM[4], secM[5], secM[2]]]
    D,V = np.linalg.eig(M)

    # initialization via SVD
    x0SVD = np.zeros(4,1)

    x0SVD[0] = np.real(V[0, 2])
    x0SVD[1] = np.real(V(1, 2))
    x0SVD[2] = np.real(V(2, 2))
    x0SVD[3] = 1.5 * np.real(D[2,2]) - .5

    # upper and lower constraints for estimated first moments
    lb = [[-np.ones(3, 1)], [np.zeros(1, 1)]]
    ub = np.ones(4, 1);

    # interior-point optimization
    scipy.optimize.minimize(objective_fun,x0SVD, bounds=(lb, ub), constraints=constraint)


def symmCone2SecM(z):
    z = np.reshape(z, (1, 4), order='F')
    muz_t = z[:, 2]
    muxx = z[:, 3] * z[:,0]**2 + (1 - z[:, 3]) / 3
    muyy = z[:, 3] * z[:,1]**2 + (1 - z[:,3]) / 3
    muzz = z[:, 3] * muz_t**2 + (1 - z[:,3]) / 3
    muxy = z[:, 3] * z[:,0] * z[:,1]
    muxz = z[:, 3] * z[:,0] * muz_t
    muyz = z[:, 3] * z[:,1] * muz_t
    out = [muxx, muyy, muzz, muxy, muxz, muyz].transpose()
    return out

def objective_fun(z, FIM, secM):
    '''Objective function'''
    sC2M = symmCone2SecM(z)
    return (sC2M.transpose() - secM) * FIM * (sC2M - secM.transpose())

def constraint(z):
    '''Constraint function'''
    return np.sum(z[:, 0:2]**2, 2) - 1


def calFIMSecondM(B, I, s, backg):
    FIM = np.zeros(6, 6);

    FIM[0, 0] = 0.5 * (s^2) * (sum(sum((B.aa) / (I + backg))));
    FIM[1, 1] = 0.5 * (s^2) * (sum(sum((B.bb) / (I + backg))));
    FIM[2, 2] = 0.5 * (s^2) * (sum(sum((B.cc) / (I + backg))));
    FIM[3, 3] = 0.5 * (s^2) * (sum(sum((B.dd) / (I + backg))));
    FIM[4, 4] = 0.5 * (s^2) * (sum(sum((B.ee) / (I + backg))));
    FIM[5, 5] = 0.5 * (s^2) * (sum(sum((B.ff) / (I + backg))));

    FIM[0, 1] = (s^2) * (sum(sum((B.ab) / (I + backg))));
    FIM[0, 2] = (s^2) * (sum(sum((B.ac) / (I + backg))));
    FIM[0, 3] = (s^2) * (sum(sum((B.ad) / (I + backg))));
    FIM[0, 4] = (s^2) * (sum(sum((B.ae) / (I + backg))));
    FIM[0, 5] = (s^2) * (sum(sum((B.af) / (I + backg))));

    FIM[1, 2] = (s^2) * (sum(sum((B.bc) / (I + backg))));
    FIM[1, 3] = (s^2) * (sum(sum((B.bd) / (I + backg))));
    FIM[1, 4] = (s^2) * (sum(sum((B.be) / (I + backg))));
    FIM[1, 5] = (s^2) * (sum(sum((B.bf) / (I + backg))));

    FIM[2, 3] = (s^2) * (sum(sum((B.cd) / (I + backg))));
    FIM[2, 4] = (s^2) * (sum(sum((B.ce) / (I + backg))));
    FIM[2, 5] = (s^2) * (sum(sum((B.cf) / (I + backg))));

    FIM[3, 4] = (s^2) * (sum(sum((B.de) / (I + backg))));
    FIM[3, 5] = (s^2) * (sum(sum((B.df) / (I + backg))));

    FIM[4, 5] = (s^2) * (sum(sum((B.ef) / (I + backg))));

    FIM = FIM + FIM.transpose();
