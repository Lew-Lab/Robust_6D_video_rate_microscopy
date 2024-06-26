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

plt.rcParams.update({'font.size': 14,})

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available()  else 'cpu'
print('Using ' + device)
NO_ORIEN_GROUP = ['flat ring gamma 0','tilted ring gamma 0','fiber gamma 0',
                  'flat surface gamma 0','tilted surface gamma 0','hemisphere gamma 0']

GAMMA_1_GROUP = ['flat ring gamma 1','tilted ring gamma 1','fiber gamma 1',
                  'flat surface gamma 1','tilted surface gamma 1','hemisphere gamma 1']

GAMMA_HALF_GROUP = ['flat ring gamma 0.5','tilted ring gamma 0.5','fiber gamma 0.5',
                  'flat surface gamma 0.5','tilted surface gamma 0.5','hemisphere gamma 0.5']


### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y

object_file = 'hemisphere gamma 1'
initial_file = 'fiber initial'
Bstruct_file = 'B_ZScanned_zf1000_61_2000'
num_of_trials = 1
factor = 1e5
ERROR_DISTRIBUTION_FLAG = True

psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.05
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 80
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 5

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

# for vs. gamma
gamma_error_1 = np.array([])
gamma_error_half = np.array([])
gamma_error_0 = np.array([])
gamma_est_1 = np.array([])
gamma_est_half = np.array([])
gamma_est_0 = np.array([])

object_gt = loadmat(os.path.join('performance test objects new', object_file+'.mat'))['object']*factor
object_initial = loadmat(os.path.join('performance test initials new', initial_file+'.mat'))['initial']*factor
Bstruct = loadmat(os.path.join('psf',Bstruct_file+'.mat'))['Bstruct']
object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
model_cpu = cpu.smolm(psf, object_size)
image_raw = model_cpu.forward(object_gt)
s_gt, m1_gt, z_gt, x_gt, y_gt = convert.m2_to_m1(object_gt,Bstruct,factor/10)

'''visualization of GT'''
brt_norm = plt.Normalize(factor/10,factor)
x,y,z = np.indices((32,32,11))
obj_space = np.array(np.zeros((31,31,10)),dtype='bool')
colors = np.zeros(obj_space.shape + (4,))
for i in range (z_gt.shape[0]):
    obj_space[x_gt[i]-15,y_gt[i]-15,z_gt[i]] = True
    colors[x_gt[i]-15,y_gt[i]-15,z_gt[i],:] = mpl.cm.Oranges(brt_norm(s_gt[i]))
ax = plt.figure().add_subplot(projection='3d')
colors[:,:,:,3] = .7
ax.voxels(x,y,z,obj_space,facecolors=colors)
# ax.set(xlabel='x', ylabel='y', zlabel='z')
# ax.set_axis_off()
ax.set_box_aspect([32,32,11])
plt.colorbar(plt.cm.ScalarMappable(norm=brt_norm, cmap='Oranges'))
# ax.view_init(azim=0, elev=90)

plt.show()


# total_voxel_num = object_gt.shape[1]*object_gt.shape[2]*object_gt.shape[3]

# convert to angle representation
# theta_gt = np.degrees(np.arccos(m1_gt[:,2]))
# phi_gt = np.degrees(np.arctan2(m1_gt[:,1],m1_gt[:,0]))
# phi_gt[np.equal(m1_gt[:,0],0) & np.equal(m1_gt[:,1],0)] = np.nan
# gamma_gt = m1_gt[:,3]

# theta_top = np.zeros((31,31))
# theta_top.fill(np.nan)
# phi_top = np.zeros((31,31))
# phi_top.fill(np.nan)
# gamma_top = np.zeros((31,31))
# gamma_top.fill(np.nan)
# z_track = np.zeros((31,31))
# for i in range (z_gt.shape[0]):
#     x = x_gt[i][0] - 15
#     y = y_gt[i][0] - 15
#     if (z_gt[i][0] >= z_track[x][y]):
#         theta_top[x][y] = theta_gt[i]
#         phi_top[x][y] = phi_gt[i]
#         gamma_top[x][y] = gamma_gt[3]
#         z_track[x][y] = z_gt[i]

# theta_norm = plt.Normalize(0,90)
# phi_norm = plt.Normalize(-180,180)
# gamma_norm = plt.Normalize(0,1)

# ax = plt.figure().add_subplot()
# ax.imshow(gamma_top,cmap='coolwarm',vmin=0,vmax=1)
# plt.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'))
# # plt.colorbar()
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# image_noisy = np.random.poisson(image_raw)
# object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, optim_param['learning_rate'], 
#                                 optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], 
#                                 optim_param['lambda_I'], device)
# s_est, m1_est, z_est, x_est, y_est = convert.m2_to_m1(object_est,Bstruct,factor/10)

# brt_norm = plt.Normalize(factor/10,factor)
# x,y,z = np.indices((32,32,11))
# obj_space = np.array(np.zeros((31,31,10)),dtype='bool')
# colors = np.zeros(obj_space.shape + (4,))
# for i in range (z_est.shape[0]):
#     obj_space[x_est[i]-15,y_est[i]-15,z_est[i]] = True
#     colors[x_est[i]-15,y_est[i]-15,z_est[i],:] = mpl.cm.Oranges(brt_norm(s_est[i]))
# ax = plt.figure().add_subplot(projection='3d')
# colors[:,:,:,3] = .7
# ax.voxels(x,y,z,obj_space,facecolors=colors)
# ax.set(xlabel='x', ylabel='y', zlabel='z')
# ax.set_box_aspect([32,32,11])
# plt.colorbar(plt.cm.ScalarMappable(norm=brt_norm, cmap='Oranges'))
# # ax.view_init(azim=0, elev=90)
# plt.show()

# convert to angle representation
# theta_est = np.degrees(np.arccos(m1_est[:,2]))
# phi_est = np.degrees(np.arctan2(m1_est[:,1],m1_est[:,0]))
# phi_est[np.equal(m1_est[:,0],0) & np.equal(m1_est[:,1],0)] = np.nan
# gamma_est = m1_est[:,3]

# theta_top = np.zeros((31,31))
# theta_top.fill(np.nan)
# phi_top = np.zeros((31,31))
# phi_top.fill(np.nan)
# gamma_top = np.zeros((31,31))
# gamma_top.fill(np.nan)
# z_track = np.zeros((31,31))
# for i in range (z_est.shape[0]):
#     x = x_est[i][0] - 15
#     y = y_est[i][0] - 15
#     if (z_est[i][0] >= z_track[x][y]):
#         theta_top[x][y] = theta_est[i]
#         phi_top[x][y] = phi_est[i]
#         gamma_top[x][y] = gamma_est[3]
#         z_track[x][y] = z_est[i]

# theta_norm = plt.Normalize(0,90)
# phi_norm = plt.Normalize(-180,180)
# gamma_norm = plt.Normalize(0,1)

# ax = plt.figure().add_subplot()
# ax.imshow(gamma_top,cmap='coolwarm',vmin=0,vmax=1)
# plt.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'))
# # plt.colorbar()
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# for trial in range (num_of_trials):
#     image_noisy = np.random.poisson(image_raw)
#     object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, optim_param['learning_rate'], 
#                                     optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], 
#                                     optim_param['lambda_I'], device)
#     s_est, m1_est, z_est, x_est, y_est = convert.m2_to_m1(object_est,Bstruct,factor/10)

#     # compare
#     acc_s, acc_orientation, acc_gamma, acc_z, acc_x, acc_y, acc_gt_s, acc_gt_theta, acc_gt_phi, fp_z, fn_z = convert.acc_eval(s_gt,m1_gt,z_gt,x_gt,y_gt,s_est,m1_est,z_est,x_est,y_est,total_voxel_num)

#     theta_est = np.degrees(np.arccos(m1_est[:,2]))
#     phi_est = np.degrees(np.arctan2(m1_est[:,1],m1_est[:,0]))
#     phi_est[np.equal(m1_est[:,0],0) & np.equal(m1_est[:,1],0)] = np.nan

#     if object_file_list[idx_obj] in NO_ORIEN_GROUP:  
#         orien_acc_total[idx_obj,trial] = np.nan
#         gamma_error_0 = np.append(gamma_error_0,acc_gamma)
#         gamma_est_0 = np.append(gamma_est_0,acc_gamma)
#     else:
#         if object_file_list[idx_obj] in GAMMA_1_GROUP:
#             gamma_error_1 = np.append(gamma_error_1,acc_gamma)
#             gamma_est_1 = np.append(gamma_est_1,acc_gamma)
#         elif object_file_list[idx_obj] in GAMMA_HALF_GROUP:
#             gamma_error_half = np.append(gamma_error_half,acc_gamma)
#             gamma_est_half = np.append(gamma_est_half,acc_gamma)
#         for i in range (acc_orientation.shape[0]):
#             theta_idx = (acc_gt_theta[i]/10).astype('int64')
#             if theta_idx == 9:
#                 theta_idx = 8
#             theta_bin[theta_idx] += acc_orientation[i]
#             theta_cnt[theta_idx] += 1
#             if np.isnan(acc_gt_phi[i]) == False:
#                 phi_idx = ((acc_gt_phi[i]+180)/20).astype('int64')
#                 if phi_idx == 18:
#                     phi_idx = 17
#                 phi_bin[phi_idx] += acc_orientation[i]
#                 phi_cnt[phi_idx] += 1
#             z_idx = z_gt[i]
#             z_bin_orien[z_idx] += acc_orientation[i]
#             z_cnt_orien[z_idx] += 1

#         orien_acc_total[idx_obj,trial] = np.sum(acc_orientation)/len(acc_orientation)

#     for i in range (len(fp_z)):
#         z_idx = fp_z[i]
#         z_bin_fp[int(z_idx)] += 1

#     for i in range (len(fn_z)):
#         z_idx = fn_z[i]
#         z_bin_fn[int(z_idx)] += 1

#     if ERROR_DISTRIBUTION_FLAG:

#         theta_norm = plt.Normalize(0,90)
#         phi_norm = plt.Normalize(-180,180)
#         acc_norm = plt.Normalize(0,90)

#         for z in range (10):

#             idx_temp = np.where(np.equal(z_gt,z))
#             x_temp = x_gt[idx_temp]
#             y_temp = y_gt[idx_temp]
#             theta_color_temp = mpl.cm.cool(theta_norm(theta_gt[idx_temp[0]]))
#             phi_color_temp = mpl.cm.rainbow(phi_norm(phi_gt[idx_temp[0]]))

#             fig_temp = plt.figure()
#             ax11 = fig_temp.add_subplot(221)
#             ax11.scatter(x_temp,y_temp,c=theta_color_temp,marker='s',s=5)
#             plt.title('theta gt at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax11.set_aspect('equal', adjustable='box')
#             fig_temp.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax11)

#             ax12 = fig_temp.add_subplot(222)
#             ax12.scatter(x_temp,y_temp,c=phi_color_temp,marker='s',alpha=1,s=5)
#             plt.title('phi gt at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax12.set_aspect('equal', adjustable='box')
#             fig_temp.colorbar(plt.cm.ScalarMappable(norm=phi_norm, cmap='rainbow'), ax=ax12)

#             idx_temp = np.where(np.equal(z_est,z))
#             x_temp = x_est[idx_temp]
#             y_temp = y_est[idx_temp]
#             theta_color_temp = mpl.cm.cool(theta_norm(theta_est[idx_temp[0]]))
#             phi_color_temp = mpl.cm.rainbow(phi_norm(phi_est[idx_temp[0]]))

#             ax13 = fig_temp.add_subplot(223)
#             ax13.scatter(x_temp,y_temp,c=theta_color_temp,marker='s',alpha=1,s=5)
#             plt.title('theta est at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax13.set_aspect('equal', adjustable='box')
#             fig_temp.colorbar(plt.cm.ScalarMappable(norm=theta_norm, cmap='cool'), ax=ax13)

#             ax14 = fig_temp.add_subplot(224)
#             ax14.scatter(x_temp,y_temp,c=phi_color_temp,marker='s',alpha=1,s=5)
#             plt.title('phi est at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax14.set_aspect('equal', adjustable='box')
#             fig_temp.colorbar(plt.cm.ScalarMappable(norm=phi_norm, cmap='rainbow'), ax=ax14)

#             mng = plt.get_current_fig_manager()
#             mng.full_screen_toggle()
#             plt.savefig('orientation' + "/file%02d.png" % z)
#             plt.close(fig_temp)

#             idx_temp_2 = np.where(np.equal(acc_z,z))
#             acc_color_temp = mpl.cm.hot(acc_norm(acc_orientation[idx_temp_2[0]]))

#             fig_temp = plt.figure()
#             ax = fig_temp.add_subplot(111)
#             ax.scatter(acc_x[idx_temp_2],acc_y[idx_temp_2],c=acc_color_temp,marker='s',alpha=1,s=50)
#             plt.title('orientation accuracy at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax14.set_aspect('equal', adjustable='box')
#             fig_temp.colorbar(plt.cm.ScalarMappable(norm=acc_norm, cmap='hot'), ax=ax)

#             mng = plt.get_current_fig_manager()
#             mng.full_screen_toggle()
#             plt.savefig('orientation accuracy' + "/file%02d.png" % z)
#             plt.close(fig_temp)

#         gamma_norm = plt.Normalize(0,1)
#         gamma_gt = m1_gt[:,3]
#         gamma_est = m1_est[:,3]

#         for z in range (10):

#             fig_gamma = plt.figure()

#             idx_temp = np.where(np.equal(z_gt,z))
#             x_temp = x_gt[idx_temp]
#             y_temp = y_gt[idx_temp]
#             gamma_color_temp = mpl.cm.coolwarm(gamma_norm(gamma_gt[idx_temp[0]]))

#             ax_g1 = fig_gamma.add_subplot(131) 
#             ax_g1.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
#             plt.title('gamma gt at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax_g1.set_aspect('equal', adjustable='box')
#             fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'), ax=ax_g1)

#             idx_temp = np.where(np.equal(z_est,z))
#             x_temp = x_est[idx_temp]
#             y_temp = y_est[idx_temp]
#             gamma_color_temp = mpl.cm.coolwarm(gamma_norm(gamma_est[idx_temp[0]]))

#             ax_g2 = fig_gamma.add_subplot(132) 
#             ax_g2.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
#             plt.title('gamma est at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax_g2.set_aspect('equal', adjustable='box')
#             fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='coolwarm'), ax=ax_g2)

#             idx_temp = np.where(np.equal(acc_z,z))
#             x_temp = acc_x[idx_temp]
#             y_temp = acc_y[idx_temp]
#             gamma_color_temp = mpl.cm.inferno(gamma_norm(acc_gamma[idx_temp[0]]))

#             ax_g3 = fig_gamma.add_subplot(133) 
#             ax_g3.scatter(x_temp,y_temp,c=gamma_color_temp,marker='s',s=5)
#             plt.title('gamma error at z = ' + str(z))
#             plt.xlim(0,60)
#             plt.ylim(0,60)
#             ax_g3.set_aspect('equal', adjustable='box')
#             fig_gamma.colorbar(plt.cm.ScalarMappable(norm=gamma_norm, cmap='inferno'), ax=ax_g3)

#             mng = plt.get_current_fig_manager()
#             mng.full_screen_toggle()
#             plt.savefig('gamma' + "/file%02d.png" % z)
#             plt.close(fig_gamma)