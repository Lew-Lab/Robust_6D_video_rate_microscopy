import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
import pyutils.GPU_module as gpu
import pyutils.CPU_module as cpu
import pyutils.m2_to_m1 as convert
import pandas as pd 
import matplotlib as mpl


### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y
object_file_list = ['flat ring gamma 1']
initial_file_list = ['flat ring initial']
# object_file_list = ['fiber gamma 1', 'fiber gamma 0.5', 'fiber gamma 0','tilted ring gamma 1', 'tilted ring gamma 0.5', 'tilted ring gamma 0'] # orientation x Z x X x Y
# initial_file_list = ['fiber initial', 'fiber initial', 'fiber initial', 'tilted ring initial', 'tilted ring initial', 'tilted ring initial']
Bstruct_file = 'B_ZScanned_zf1000_PSFsz_61'
image_file = ''
num_of_trials = 5
factor = 100

psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 350
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 5

accuracy = np.zeros((len(object_file_list),num_of_trials))

for idx_obj in range (len(object_file_list)):
    object_gt = loadmat(os.path.join('performance test objects', object_file_list[idx_obj]+'.mat'))['object']*factor
    object_initial = loadmat(os.path.join('performance test initials', initial_file_list[idx_obj]+'.mat'))['initial']*factor
    Bstruct = loadmat(os.path.join('psf',Bstruct_file+'.mat'))['Bstruct']
    object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
    object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
    model_cpu = cpu.smolm(psf, object_size)
    image_raw = model_cpu.forward(object_gt)
    s_gt, m1_gt, z_gt, x_gt, y_gt = convert.m2_to_m1(object_gt,Bstruct,factor*10)

    # convert to angle representation

    theta_gt = np.degrees(np.arccos(m1_gt[:,2]))
    phi_gt = np.degrees(np.arctan(m1_gt[:,1]/m1_gt[:,0]))

    # # plot the object
    # theta_normed = plt.Normalize(np.min(theta_gt), np.max(theta_gt))
    # theta_color = mpl.cm.cool(theta_normed(theta_gt))

    # fig1 = plt.figure()
    # ax1 = fig1.add_subplot(111, projection='3d')
    # a1 = ax1.scatter(x_gt, y_gt, z_gt, c=theta_color, marker='o', alpha=0.4)
    # fig1.colorbar(plt.cm.ScalarMappable(norm=theta_normed, cmap='cool'), ax=ax1)

    # ax1.set_xlabel('X voxel')
    # ax1.set_ylabel('Y voxel')
    # ax1.set_zlabel('Z voxel')

    # plt.title('theta gt')

    # phi_normed = plt.Normalize(np.min(phi_gt), np.max(phi_gt))
    # phi_color = mpl.cm.hsv(phi_normed(phi_gt))

    # fig2 = plt.figure()
    # ax2 = fig2.add_subplot(111, projection='3d')
    # a2 = ax2.scatter(x_gt, y_gt, z_gt, c=phi_color, marker='o', alpha=0.4)
    # fig2.colorbar(plt.cm.ScalarMappable(norm=phi_normed,cmap='hsv'), ax=ax2)

    # ax2.set_xlabel('X voxel')
    # ax2.set_ylabel('Y voxel')
    # ax2.set_zlabel('Z voxel')

    # plt.title('phi gt')

    # s_normed = plt.Normalize(np.min(s_gt), np.max(s_gt))
    # s_color = mpl.cm.hsv(s_normed(s_gt))

    # fig3 = plt.figure()
    # ax3 = fig3.add_subplot(111, projection='3d')
    # a3 = ax3.scatter(x_gt, y_gt, z_gt, c=s_color, marker='o', alpha=0.4)
    # fig3.colorbar(plt.cm.ScalarMappable(norm=s_normed,cmap='hsv'), ax=ax3)

    # ax3.set_xlabel('X voxel')
    # ax3.set_ylabel('Y voxel')
    # ax3.set_zlabel('Z voxel')

    # plt.title('s gt')

    # gamma_gt = m1_gt[:,3]
    # gamma_normed = plt.Normalize(np.min(gamma_gt), np.max(gamma_gt))
    # gamma_color = mpl.cm.hsv(gamma_normed(gamma_gt))

    # fig4 = plt.figure()
    # ax4 = fig4.add_subplot(111, projection='3d')
    # a4 = ax4.scatter(x_gt, y_gt, z_gt, c=gamma_color, marker='o', alpha=0.4)
    # fig4.colorbar(plt.cm.ScalarMappable(norm=gamma_normed,cmap='hsv'), ax=ax4)

    # ax4.set_xlabel('X voxel')
    # ax4.set_ylabel('Y voxel')
    # ax4.set_zlabel('Z voxel')

    # plt.title('gamma gt')

    plt.show()


    for trial in range (num_of_trials):
        image_noisy = np.random.poisson(image_raw)
        object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, optim_param['learning_rate'], 
                                        optim_param['max_iter'], optim_param['lambda_TV'], optim_param['lambda_TV'], 
                                        optim_param['lambda_I'], device)
        s_est, m1_est, z_est, x_est, y_est = convert.m2_to_m1(object_est,Bstruct,1000)

        

        # compare
        acc_s, acc_orientation, acc_gamma, acc_z, acc_x, acc_y= convert.orientation_acc(s_gt,m1_gt,z_gt,x_gt,
                                                                                        y_gt,s_est,m1_est,z_est,x_est,y_est)

        theta_est = np.degrees(np.arccos(m1_est[:,2]))
        phi_est = np.degrees(np.arctan(m1_est[:,1]/m1_est[:,0]))

        theta_norm = plt.Normalize(0, 90)
        theta_gt_color = mpl.cm.cool(theta_norm(theta_gt))
        theta_est_color = mpl.cm.cool(theta_norm(theta_est))

        fig1 = plt.figure()
        ax11 = fig1.add_subplot(221, projection='3d')
        a11 = ax11.scatter(x_gt, y_gt, z_gt, c=theta_gt_color, marker='o', alpha=0.4)
        ax12 = fig1.add_subplot(222, projection='3d')
        a12 = ax12.scatter(x_est, y_est, z_est, c=theta_est_color, marker='^', alpha=0.4)
        fig1.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax11)
        fig1.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax12)


        ax11.set_xlabel('X voxel')
        ax11.set_ylabel('Y voxel')
        ax11.set_zlabel('Z voxel')

        plt.title('theta')

        phi_norm = plt.Normalize(-90,90)
        phi_gt_color = mpl.cm.hsv(phi_norm(phi_gt))
        phi_est_color = mpl.cm.hsv(phi_norm(phi_est))

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        a21 = ax2.scatter(x_gt, y_gt, z_gt, c=phi_gt_color, marker='o', alpha=0.4)
        a22 = ax2.scatter(x_est, y_est, z_est, c=phi_est_color, marker='^', alpha=0.4)
        fig2.colorbar(plt.cm.ScalarMappable(norm=phi_norm, cmap='hsv'), ax=ax2)

        ax2.set_xlabel('X voxel')
        ax2.set_ylabel('Y voxel')
        ax2.set_zlabel('Z voxel')

        plt.title('phi')


        gamma_gt = m1_gt[:,3]
        gamma_est = m1_est[:,3]
        gamma_norm = plt.Normalize(np.min([np.min(gamma_gt),np.min(gamma_est)]),np.max([np.max(gamma_gt),np.max(gamma_est)]))
        gamma_norm = plt.Normalize(0,1)
        gamma_gt_color = mpl.cm.coolwarm(gamma_norm(gamma_gt))
        gamma_est_color = mpl.cm.coolwarm(gamma_norm(gamma_est))

        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111,projection='3d')
        a31 = ax3.scatter(x_gt, y_gt, z_gt, c=gamma_gt_color, marker='o', alpha=0.4)
        a32 = ax3.scatter(x_est, y_est, z_est, c=gamma_est_color, marker='^', alpha=0.4)
        fig3.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'), ax=ax3)

        ax3.set_xlabel('X voxel')
        ax3.set_ylabel('Y voxel')
        ax3.set_zlabel('Z voxel')

        plt.title('phi')

        
        # phi_normed = plt.Normalize(np.min(phi_est), np.max(phi_est))
        # phi_color = mpl.cm.hsv(phi_normed(phi_est))

        # fig4 = plt.figure()
        # ax4 = fig4.add_subplot(111, projection='3d')
        # a4 = ax4.scatter(x_est, y_est, z_est, c=phi_color, marker='o', alpha=0.4)
        # fig4.colorbar(plt.cm.ScalarMappable(norm=phi_normed,cmap='hsv'), ax=ax4)

        # ax4.set_xlabel('X voxel')
        # ax4.set_ylabel('Y voxel')
        # ax4.set_zlabel('Z voxel')

        # plt.title('phi est')

        acc_orientation_normed = plt.Normalize(np.min(acc_orientation),np.max(acc_orientation))
        acc_orientation_color = mpl.cm.inferno(acc_orientation_normed(acc_orientation))

        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111, projection='3d')
        a5 = ax5.scatter(acc_x, acc_y, acc_z, c=acc_orientation_color, marker='o', alpha=0.4)
        fig5.colorbar(plt.cm.ScalarMappable(norm=acc_orientation_normed, cmap='inferno'), ax=ax5)

        ax5.set_xlabel('X voxel')
        ax5.set_ylabel('Y voxel')
        ax5.set_zlabel('Z voxel')

        plt.title('orientation error')

        # acc_s_normed = plt.Normalize(np.min(acc_s),np.max(acc_s))
        # acc_s_color = mpl.cm.hot(acc_s_normed(acc_s))

        # fig6 = plt.figure()
        # ax6 = fig6.add_subplot(111, projection='3d')
        # a6 = ax6.scatter(acc_x, acc_y, acc_z, c=acc_s_color, marker='o', alpha=0.4)
        # fig6.colorbar(plt.cm.ScalarMappable(norm=acc_s_normed, cmap='hot'), ax=ax6)

        # ax6.set_xlabel('X voxel')
        # ax6.set_ylabel('Y voxel')
        # ax6.set_zlabel('Z voxel')

        # plt.title('s error')

        # acc_gamma_normed = plt.Normalize(np.min(acc_gamma),np.max(acc_gamma))
        # acc_gamma_color = mpl.cm.hot(acc_gamma_normed(acc_gamma))

        # fig7 = plt.figure()
        # ax7 = fig7.add_subplot(111, projection='3d')
        # a7 = ax7.scatter(acc_x, acc_y, acc_z, c=acc_gamma_color, marker='o', alpha=0.4)
        # fig7.colorbar(plt.cm.ScalarMappable(norm=acc_gamma_normed, cmap='hot'), ax=ax7)

        # ax7.set_xlabel('X voxel')
        # ax7.set_ylabel('Y voxel')
        # ax7.set_zlabel('Z voxel')

        # plt.title('gamma error')

        plt.show()

        accuracy[idx_obj,trial] = np.sum(acc_orientation)/len(acc_orientation)

pd.DataFrame(accuracy).to_csv("test.csv")