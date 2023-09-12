import os
import numpy as np
import torch
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import io
import scipy
from scipy import optimize
# from pyutils.autograd_minimize import minimize as torchminimize

def m2_to_m1(object, Bstruct, brightness_threshold, device): 

    m_xx = object[0,:,:,:]
    m_yy = object[1,:,:,:]
    m_zz = object[2,:,:,:]

    s = m_xx + m_yy + m_zz
    s = np.matrix.flatten(s, order = 'F')

    m_xy = object[3,:,:,:]
    m_xz = object[4,:,:,:]
    m_yz = object[5,:,:,:]

    m2 = [np.transpose(np.matrix.flatten(m_xx, order = 'F')), 
          np.transpose(np.matrix.flatten(m_yy, order = 'F')), np.transpose(np.matrix.flatten(m_zz, order = 'F')), 
          np.transpose(np.matrix.flatten(m_xy, order = 'F')), np.transpose(np.matrix.flatten(m_xz, order = 'F')), 
          np.transpose(np.matrix.flatten(m_yz, order = 'F'))]

    m2 = np.concatenate(m2, axis = 0)
    m2 = np.reshape(m2, (np.shape(m2)[0]//6, 6), order = 'F')

    [z, x, y] = np.unravel_index(np.argwhere(s > brightness_threshold), np.shape(m_xx), order = 'F')

    # retrieve second moments of voxels with localizations
    m2 = np.squeeze(m2[np.argwhere(s > brightness_threshold),:])
    # print(len(m2))

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
        b = BsList[0, zSlice];
        B = BHaList[0, zSlice]
        sumNorm = sumNormList[0, zSlice][0,0]

        # estimate and store first moments
        estM1 = secondM2SymmConeWeighted(b, B, sumNorm, secM/signal, signal, 1e-12, device)

        m1[loc,:] = estM1

    return s[np.argwhere(s > brightness_threshold)], m1, z, x, y


def secondM2SymmConeWeighted(b, B, sumNorm, secM, signal, backg, device):
    # define a least squares problem structure

    # print('Using ' + device)

    # construct the M matrix
    M = torch.tensor(np.array([[secM[0], secM[3], secM[4]], 
                            [secM[3], secM[1], secM[5]],
                            [secM[4], secM[5], secM[2]]]));
    # print(M);

    # Define second moment basis images
    bXX = b['XX'][0,0]; bYY = b['YY'][0,0]; bZZ = b['ZZ'][0,0]; bXY = b['XY'][0,0]; bXZ = b['XZ'][0,0]; bYZ = b['YZ'][0,0]

    # Convert to PyTorch tensors
    bXX = torch.tensor(bXX); bYY = torch.tensor(bYY); bZZ = torch.tensor(bZZ);
    bXY = torch.tensor(bXY); bXZ = torch.tensor(bXZ); bYZ = torch.tensor(bYZ);

    Ix = secM[0] * bXX + secM[1] * bYY + secM[2] * bZZ + secM[3] * bXY + secM[4] * bXZ + secM[5] * bYZ
    Ix[Ix<0] = 0

    # Normalize second moment info and multiply by photon number
    Ix  = Ix / sumNorm; Ix = signal * Ix

    # Compute Fisher Information Matrix
    FIM = calFIMSecondM(B, Ix, signal, backg)

    # Compute eigenstuffs of M and retrieve largest eigenvalues
    eigvalue, eigvector = torch.linalg.eig(M)
    eigvector_real = torch.real(eigvector)
    eigvalue_real = torch.real(eigvalue)

    idx_max_D = torch.argmax(eigvalue_real)

    # initialization via SVD
    x0SVD = torch.zeros((4,1))

    # Initial value based the largest eigenvalue
    #------------------------------------------------------------
    x0SVD[0] = torch.real(eigvector_real[0, idx_max_D])
    x0SVD[1] = torch.real(eigvector_real[1, idx_max_D])
    x0SVD[2] = torch.real(eigvector_real[2, idx_max_D])
    x0SVD[3] = 1.5 * torch.max(eigvalue_real[idx_max_D]) - .5

    # upper and lower constraints for estimated first moments
    bound = [(-1, 1), (-1, 1), (-1, 1), (0, 1)]

    cons = {'type':'ineq', 'fun': constraint}

    # interior-point optimization
    M1estresults = optimize.minimize(objective_fun, x0SVD, bounds = bound, constraints = cons, args = (FIM, np.matrix(secM),))
    # M1estresults = torchminimize(objective_fun, x0SVD.detach().numpy(), bounds = bound, constraints = cons, args = (FIM, np.matrix(secM),), torch_device=device)
    # print(M1estresults)
    estM1 = M1estresults['x']
    # print(estM1)
    # print(np.linalg.norm(estM1[0:2]))

    # If mu_z is negative, flip the first moment vector.
    if estM1[2] < 0:
        estM1[:3] = -estM1[:3]

    return estM1

def calFIMSecondM(B, I, s, backg):
    FIM = torch.zeros((6, 6), dtype=torch.double)

    Baa = B['aa'][0,0]; Bbb = B['bb'][0,0]; Bcc = B['cc'][0,0]; Bdd = B['dd'][0,0]; Bee = B['ee'][0,0]; Bff = B['ff'][0,0]
    Bab = B['ab'][0,0]; Bac = B['ac'][0,0]; Bad = B['ad'][0,0]; Bae = B['ae'][0,0]; Baf = B['af'][0,0]
    Bbc = B['bc'][0,0]; Bbd = B['bd'][0,0]; Bbe = B['be'][0,0]; Bbf = B['bf'][0,0]
    Bcd = B['cd'][0,0]; Bce = B['ce'][0,0]; Bcf = B['cf'][0,0]
    Bde = B['de'][0,0]; Bdf = B['df'][0,0]
    Bef = B['ef'][0,0]


    FIM[0, 0] = 0.5 * (s**2) * (torch.sum(torch.sum((Baa) / (I + backg))))
    FIM[1, 1] = 0.5 * (s**2) * (torch.sum(torch.sum((Bbb) / (I + backg))))
    FIM[2, 2] = 0.5 * (s**2) * (torch.sum(torch.sum((Bcc) / (I + backg))))
    FIM[3, 3] = 0.5 * (s**2) * (torch.sum(torch.sum((Bdd) / (I + backg))))
    FIM[4, 4] = 0.5 * (s**2) * (torch.sum(torch.sum((Bee) / (I + backg))))
    FIM[5, 5] = 0.5 * (s**2) * (torch.sum(torch.sum((Bff) / (I + backg))))

    FIM[0, 1] = (s**2) * (torch.sum(torch.sum((Bab) / (I + backg))))
    FIM[0, 2] = (s**2) * (torch.sum(torch.sum((Bac) / (I + backg))))
    FIM[0, 3] = (s**2) * (torch.sum(torch.sum((Bad) / (I + backg))))
    FIM[0, 4] = (s**2) * (torch.sum(torch.sum((Bae) / (I + backg))))
    FIM[0, 5] = (s**2) * (torch.sum(torch.sum((Baf) / (I + backg))))

    FIM[1, 2] = (s**2) * (torch.sum(torch.sum((Bbc) / (I + backg))))
    FIM[1, 3] = (s**2) * (torch.sum(torch.sum((Bbd) / (I + backg))))
    FIM[1, 4] = (s**2) * (torch.sum(torch.sum((Bbe) / (I + backg))))
    FIM[1, 5] = (s**2) * (torch.sum(torch.sum((Bbf) / (I + backg))))

    FIM[2, 3] = (s**2) * (torch.sum(torch.sum((Bcd) / (I + backg))))
    FIM[2, 4] = (s**2) * (torch.sum(torch.sum((Bce) / (I + backg))))
    FIM[2, 5] = (s**2) * (torch.sum(torch.sum((Bcf) / (I + backg))))

    FIM[3, 4] = (s**2) * (torch.sum(torch.sum((Bde) / (I + backg))))
    FIM[3, 5] = (s**2) * (torch.sum(torch.sum((Bdf) / (I + backg))))

    FIM[4, 5] = (s**2) * (torch.sum(torch.sum((Bef) / (I + backg))))

    FIM = FIM + torch.transpose(FIM, 0, 1)

    return FIM

def symmCone2SecM(z):
    # print(z.type())
    z = torch.reshape(torch.tensor(z, dtype=torch.double), (1, 4))
    # print(z.type())
    muxx = z[:,3] * z[:,0]**2 + (1 - z[:,3]) / 3
    muyy = z[:,3] * z[:,1]**2 + (1 - z[:,3]) / 3
    muzz = z[:,3] * z[:,2]**2 + (1 - z[:,3]) / 3
    muxy = z[:,3] * z[:,0] * z[:,1]
    muxz = z[:,3] * z[:,0] * z[:,2]
    muyz = z[:,3] * z[:,1] * z[:,2]
    out = torch.tensor([[muxx[0]], [muyy[0]], [muzz[0]], [muxy[0]], [muxz[0]], [muyz[0]]])
    return out


def objective_fun(z, FIM, secM):
    sC2SM = symmCone2SecM(z)
    return torch.matmul(torch.matmul((sC2SM.T - secM), FIM), (sC2SM - secM.T))[0,0]

def constraint(z):
    '''Constraint function'''
    return torch.sum(torch.tensor(z[0:3])**2, 0) - 1

def acc_eval(s_gt, m1_gt, z_gt, x_gt, y_gt, s_est, m1_est, z_est, x_est, y_est, total_voxel_num):
    pointer_gt = 0 
    pointer_est = 0
    acc_orientation_unweighted = np.array([])
    acc_orientation = np.array([])
    acc_gamma = np.array([])
    acc_s = np.array([])
    acc_z = np.array([])
    acc_x = np.array([])
    acc_y = np.array([])
    acc_gt_theta = np.array([])
    acc_gt_phi = np.array([])
    acc_gt_s = np.array([])
    fp = np.array([])
    fp_z = np.array([])
    fn = np.array([])
    fn_z = np.array([])

    while (pointer_gt < len(x_gt) and pointer_est < len(x_est)):

        if y_gt[pointer_gt][0] > y_est[pointer_est][0]:
            fp_z = np.append(fp_z,z_est[pointer_est][0])
            pointer_est += 1
        elif y_gt[pointer_gt][0] < y_est[pointer_est][0]:
            fn_z = np.append(fn_z,z_gt[pointer_gt][0])
            pointer_gt += 1
        elif y_gt[pointer_gt][0] == y_est[pointer_est][0]:
            if x_gt[pointer_gt][0] > x_est[pointer_est][0]:
                fp_z = np.append(fp_z,z_est[pointer_est][0])
                pointer_est += 1
            elif x_gt[pointer_gt][0] < x_est[pointer_est][0]:
                fn_z = np.append(fn_z,z_gt[pointer_gt][0])
                pointer_gt += 1
            elif x_gt[pointer_gt][0] == x_est[pointer_est][0]:
                if z_gt[pointer_gt][0] > z_est[pointer_est][0]:
                    fp_z = np.append(fp_z,z_est[pointer_est][0])
                    pointer_est += 1
                elif z_gt[pointer_gt][0] < z_est[pointer_est][0]:
                    fn_z = np.append(fn_z,z_gt[pointer_gt][0])
                    pointer_gt += 1
                elif z_gt[pointer_gt][0] == z_est[pointer_est][0]:
                    a_mu = m1_gt[pointer_gt,:3]
                    b_mu = m1_est[pointer_est,:3]
                    a_gamma = m1_gt[pointer_gt,3]
                    b_gamma = m1_est[pointer_est,3]
                    acc1 = np.degrees(np.arccos(np.sum(a_mu*b_mu)/(np.sqrt(np.sum(a_mu**2))*np.sqrt(np.sum(b_mu**2)))))
                    acc2 = np.degrees(np.arccos(np.sum(-a_mu*b_mu)/(np.sqrt(np.sum(a_mu**2))*np.sqrt(np.sum(b_mu**2)))))
                    acc_orientation = np.append(acc_orientation,np.min(np.asarray([acc1,acc2])))
                    acc_gamma = np.append(acc_gamma, np.abs(a_gamma-b_gamma))
                    acc_s = np.append(acc_s, np.abs(s_gt[pointer_gt][0]-s_est[pointer_est][0])/s_gt[pointer_gt][0])
                    acc_z = np.append(acc_z, z_gt[pointer_gt][0])
                    acc_x = np.append(acc_x, x_gt[pointer_gt][0])
                    acc_y = np.append(acc_y, y_gt[pointer_gt][0])
                    acc_gt_theta = np.append(acc_gt_theta, np.degrees(np.arccos(a_mu[2])))
                    if a_mu[1]==0 and a_mu[0] == 0:
                        acc_gt_phi = np.append(acc_gt_phi,np.nan)
                    else:
                        acc_gt_phi = np.append(acc_gt_phi,np.degrees(np.arctan2(a_mu[1],a_mu[0])))
                    acc_gt_s = np.append(acc_gt_s,s_gt[pointer_gt][0])
                    pointer_gt += 1
                    pointer_est += 1

    acc_s = np.asarray(acc_s)
    acc_orientation = np.asarray(acc_orientation)
    acc_gamma = np.asarray(acc_gamma)
    acc_z = np.asarray(acc_z)
    acc_x = np.asarray(acc_x)
    acc_y = np.asarray(acc_y)
    acc_gt_theta = np.asarray(acc_gt_theta)
    acc_gt_phi = np.asarray(acc_gt_phi)

    return acc_s, acc_orientation, acc_gamma, acc_z, acc_x, acc_y, acc_gt_s, acc_gt_theta, acc_gt_phi, fp_z, fn_z
