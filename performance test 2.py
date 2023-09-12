import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import time
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

st = time.time()

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available()  else 'cpu'
print('Using ' + device)
NO_ORIEN_GROUP = ['flat ring gamma 0','tilted ring gamma 0','fiber gamma 0',
                  'flat surface gamma 0','tilted surface gamma 0','hemisphere gamma 0']

### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y

### flat ring objects (comment out if not using entire list) ###
""" object_file_list = ['flat ring gamma 1','flat ring gamma 0.5','flat ring gamma 0',
                    'tilted ring gamma 1','tilted ring gamma 0.5','tilted ring gamma 0',
                    'fiber gamma 1','fiber gamma 0.5', 'fiber gamma 0']
initial_file_list = ['flat ring initial','flat ring initial','flat ring initial',
                     'tilted ring initial','tilted ring initial','tilted ring initial',
                     'fiber initial','fiber initial','fiber initial'] """

### hemisphere objects (comment out if not using entire list) ###
""" object_file_list = ['hemisphere gamma 1','hemisphere gamma 0.5','hemisphere gamma 0',
                    'flat surface gamma 1','flat surface gamma 0.5', 'flat surface gamma 0',
                    'tilted surface gamma 1','tilted surface gamma 0.5','tilted surface gamma 0']
initial_file_list = ['hemisphere initial','hemisphere initial','hemisphere initial',
                     'flat surface initial','flat surface initial','flat surface initial',
                     'tilted surface initial','tilted surface initial','tilted surface initial']
 """
### testing one object at a time ###
object_file_list = ['hemisphere gamma 1']
initial_file_list = ['hemisphere initial']

object_group_name = 'weighted 1d smax 10e4'
Bstruct_file = 'B_ZScanned_zf1000_61_2000'
image_file = ''
num_of_trials = 1
factor = 10000
ERROR_DISTRIBUTION_FLAG = True

psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.05
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 200
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 5

orien_acc_total = np.zeros((len(object_file_list),num_of_trials))
fp_total = np.zeros((len(object_file_list),num_of_trials))
fn_total = np.zeros((len(object_file_list),num_of_trials))
gamma_acc_total = np.zeros((len(object_file_list),num_of_trials))
s_acc_total = np.zeros((len(object_file_list),num_of_trials))

# for orientation error vs. gt theta/phi plot
theta_bin = np.zeros((9,1))
theta_cnt = np.zeros((9,1))
theta = np.linspace(np.pi/2/18,np.pi/2-np.pi/2/18,num=9,endpoint=True)
phi_bin = np.zeros((18,1))
phi_cnt = np.zeros((18,1))
phi = np.linspace(-np.pi+np.pi/18,np.pi-np.pi/18,num=18,endpoint=True)

# for orientation error vs. z
z_bin_orien = np.zeros((10,1))
z_cnt_orien = np.zeros((10,1))
z_bin_fp = np.zeros((10,1))
z_bin_fn = np.zeros((10,1))

for idx_obj in range (len(object_file_list)):
    object_gt = loadmat(os.path.join('performance test objects new', object_file_list[idx_obj]+'.mat'))['object']*factor
    object_initial = loadmat(os.path.join('performance test initials new', initial_file_list[idx_obj]+'.mat'))['initial']*factor
    Bstruct = loadmat(os.path.join('psf',Bstruct_file+'.mat'))['Bstruct']
    object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
    object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
    model_cpu = cpu.smolm(psf, object_size)
    image_raw = model_cpu.forward(object_gt)
    s_gt, m1_gt, z_gt, x_gt, y_gt = convert.m2_to_m1(object_gt,Bstruct,factor/10, device)
    total_voxel_num = object_gt.shape[1]*object_gt.shape[2]*object_gt.shape[3]

    # convert to angle representation
    theta_gt = np.degrees(np.arccos(m1_gt[:,2]))
    phi_gt = np.degrees(np.arctan2(m1_gt[:,1],m1_gt[:,0]))
    phi_gt[np.equal(m1_gt[:,0],0) & np.equal(m1_gt[:,1],0)] = np.nan

    for trial in range (num_of_trials):
        image_noisy = np.random.poisson(image_raw)
        object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, optim_param['learning_rate'], 
                                        optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], 
                                        optim_param['lambda_I'], device)
        s_est, m1_est, z_est, x_est, y_est = convert.m2_to_m1(object_est,Bstruct,factor/10, device)

        # compare
        acc_s, acc_orientation, acc_gamma, acc_z, acc_x, acc_y, acc_gt_s, acc_gt_theta, acc_gt_phi, fp_z, fn_z = convert.acc_eval(s_gt,m1_gt,z_gt,x_gt,y_gt,s_est,m1_est,z_est,x_est,y_est,total_voxel_num)

        theta_est = np.degrees(np.arccos(m1_est[:,2]))
        phi_est = np.degrees(np.arctan2(m1_est[:,1],m1_est[:,0]))
        phi_est[np.equal(m1_est[:,0],0) & np.equal(m1_est[:,1],0)] = np.nan

        if object_file_list[idx_obj] in NO_ORIEN_GROUP:  
            orien_acc_total[idx_obj,trial] = np.nan
        else:
            for i in range (acc_orientation.shape[0]):
                theta_idx = (acc_gt_theta[i]/10).astype('int64')
                if theta_idx == 9:
                    theta_idx = 8
                theta_bin[theta_idx] += acc_orientation[i]
                theta_cnt[theta_idx] += 1
                if np.isnan(acc_gt_phi[i]) == False:
                    phi_idx = ((acc_gt_phi[i]+180)/20).astype('int64')
                    if phi_idx == 18:
                        phi_idx = 17
                    phi_bin[phi_idx] += acc_orientation[i]
                    phi_cnt[phi_idx] += 1
                z_idx = z_gt[i]
                z_bin_orien[z_idx] += acc_orientation[i]
                z_cnt_orien[z_idx] += 1

            orien_acc_total[idx_obj,trial] = np.sum(acc_orientation)/len(acc_orientation)

        for i in range (len(fp_z)):
            z_idx = fp_z[i]
            z_bin_fp[int(z_idx)] += 1

        for i in range (len(fn_z)):
            z_idx = fn_z[i]
            z_bin_fn[int(z_idx)] += 1

        if ERROR_DISTRIBUTION_FLAG:

            theta_norm = plt.Normalize(0,90)
            phi_norm = plt.Normalize(-180,180)
            acc_norm = plt.Normalize(0,90)

            for z in range (10):

                idx_temp = np.where(np.equal(z_gt,z))
                x_temp = x_gt[idx_temp]
                y_temp = y_gt[idx_temp]
                theta_color_temp = mpl.cm.cool(theta_norm(theta_gt[idx_temp[0]]))
                phi_color_temp = mpl.cm.rainbow(phi_norm(phi_gt[idx_temp[0]]))

                fig_temp = plt.figure()
                ax11 = fig_temp.add_subplot(221)
                ax11.scatter(x_temp,y_temp,c=theta_color_temp,marker='s',s=5)
                plt.title('theta gt at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax11.set_aspect('equal', adjustable='box')
                fig_temp.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax11)

                ax12 = fig_temp.add_subplot(222)
                ax12.scatter(x_temp,y_temp,c=phi_color_temp,marker='s',alpha=1,s=5)
                plt.title('phi gt at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax12.set_aspect('equal', adjustable='box')
                fig_temp.colorbar(plt.cm.ScalarMappable(norm=phi_norm, cmap='rainbow'), ax=ax12)

                idx_temp = np.where(np.equal(z_est,z))
                x_temp = x_est[idx_temp]
                y_temp = y_est[idx_temp]
                theta_color_temp = mpl.cm.cool(theta_norm(theta_est[idx_temp[0]]))
                phi_color_temp = mpl.cm.rainbow(phi_norm(phi_est[idx_temp[0]]))

                ax13 = fig_temp.add_subplot(223)
                ax13.scatter(x_temp,y_temp,c=theta_color_temp,marker='s',alpha=1,s=5)
                plt.title('theta est at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax13.set_aspect('equal', adjustable='box')
                fig_temp.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax13)

                ax14 = fig_temp.add_subplot(224)
                ax14.scatter(x_temp,y_temp,c=phi_color_temp,marker='s',alpha=1,s=5)
                plt.title('phi est at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax14.set_aspect('equal', adjustable='box')
                fig_temp.colorbar(plt.cm.ScalarMappable(norm=phi_norm, cmap='rainbow'), ax=ax14)

                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.savefig('orientation' + "/file%02d.png" % z)
                plt.close(fig_temp)

                idx_temp_2 = np.where(np.equal(acc_z,z))
                acc_color_temp = mpl.cm.hot(acc_norm(acc_orientation[idx_temp_2[0]]))

                fig_temp = plt.figure()
                ax = fig_temp.add_subplot(111)
                ax.scatter(acc_x[idx_temp_2],acc_y[idx_temp_2],c=acc_color_temp,marker='s',alpha=1,s=50)
                plt.title('orientation accuracy at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax14.set_aspect('equal', adjustable='box')
                fig_temp.colorbar(plt.cm.ScalarMappable(norm=acc_norm, cmap='hot'), ax=ax)

                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.savefig('orientation accuracy' + "/file%02d.png" % z)
                plt.close(fig_temp)

            gamma_norm = plt.Normalize(0,1)
            gamma_gt = m1_gt[:,3]
            gamma_est = m1_est[:,3]

            for z in range (10):

                fig_gamma = plt.figure()

                idx_temp = np.where(np.equal(z_gt,z))
                x_temp = x_gt[idx_temp]
                y_temp = y_gt[idx_temp]
                gamma_color_temp = mpl.cm.coolwarm(gamma_norm(gamma_gt[idx_temp[0]]))

                ax_g1 = fig_gamma.add_subplot(131) 
                ax_g1.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
                plt.title('gamma gt at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax_g1.set_aspect('equal', adjustable='box')
                fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'), ax=ax_g1)

                idx_temp = np.where(np.equal(z_est,z))
                x_temp = x_est[idx_temp]
                y_temp = y_est[idx_temp]
                gamma_color_temp = mpl.cm.coolwarm(gamma_norm(gamma_est[idx_temp[0]]))

                ax_g2 = fig_gamma.add_subplot(132) 
                ax_g2.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
                plt.title('gamma est at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax_g2.set_aspect('equal', adjustable='box')
                fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'), ax=ax_g2)

                idx_temp = np.where(np.equal(acc_z,z))
                x_temp = acc_x[idx_temp]
                y_temp = acc_y[idx_temp]
                gamma_color_temp = mpl.cm.inferno(gamma_norm(acc_gamma[idx_temp[0]]))

                ax_g3 = fig_gamma.add_subplot(133) 
                ax_g3.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
                plt.title('gamma error at z = ' + str(z))
                plt.xlim(0,60)
                plt.ylim(0,60)
                ax_g3.set_aspect('equal', adjustable='box')
                fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='inferno'), ax=ax_g3)

                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.savefig('gamma' + "/file%02d.png" % z)
                plt.close(fig_gamma)

        fp_total[idx_obj,trial] = np.sum(fp_z)/total_voxel_num
        fn_total[idx_obj,trial] = np.sum(fn_z)/total_voxel_num
        gamma_acc_total[idx_obj,trial] = np.sum(acc_gamma)/len(acc_gamma)
        s_acc_total[idx_obj,trial] = np.sum(acc_s)/len(acc_s)

# pd.DataFrame(np.concatenate((phi_bin,phi_cnt),axis=1)).to_csv(object_group_name + ' phi vs o_acc.csv')
# pd.DataFrame(np.concatenate((theta_bin,theta_cnt),axis=1)).to_csv(object_group_name + ' theta vs o_acc.csv')
# pd.DataFrame(orien_acc_total).to_csv(object_group_name + ' orien acc summary.csv')
# pd.DataFrame(fp_total).to_csv(object_group_name + ' fp summary.csv')
# pd.DataFrame(fn_total).to_csv(object_group_name + ' fn summary.csv')
# pd.DataFrame(gamma_acc_total).to_csv(object_group_name + ' gamma summary.csv')
# pd.DataFrame(s_acc_total).to_csv(object_group_name + ' brt summary.csv')

for z in range (len(z_cnt_orien)):
    if z_cnt_orien[z] != 0:
        z_bin_orien[z] = z_bin_orien[z]/z_cnt_orien[z]

fig_z_hist = plt.figure()
ax_zhist1 = fig_z_hist.add_subplot(221)
ax_zhist1.bar(np.linspace(0,9,num=10,endpoint=True),z_bin_orien.squeeze())
ax_zhist1.set_title('z vs. orientation error')
ax_zhist1.set_xlabel('z layer')
ax_zhist1.set_ylabel('average orien. acc.')

z_voxel_cnt = object_gt.shape[2]*object_gt.shape[3]*len(object_file_list)*num_of_trials
ax_zhist2 = fig_z_hist.add_subplot(222)
ax_zhist2.bar(np.linspace(0,9,num=10,endpoint=True),z_bin_fp.squeeze()/z_voxel_cnt/2)
ax_zhist2.set_title('z vs. fp rate')
ax_zhist2.set_xlabel('z layer')
ax_zhist2.set_ylabel('average fp rate')

ax_zhist3 = fig_z_hist.add_subplot(223)
ax_zhist3.bar(np.linspace(0,9,num=10,endpoint=True),z_bin_fn.squeeze()/z_voxel_cnt/2)
ax_zhist3.set_title('z vs. fn rate')
ax_zhist3.set_xlabel('z layer')
ax_zhist3.set_ylabel('average fn rate')

plt.show()
et = time.time()
elapsed = et-st; print('Elapsed: ' + str(elapsed))