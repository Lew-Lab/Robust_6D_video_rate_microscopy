import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
import util.GPU_module as gpu
import util.CPU_module as cpu
import util.obj_func as obj
import util.plot_func as plot

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'dsf_raMVR_psfSZ_121_zstack_20' # channel number x orientation x Z x X x Y
object_file = 'sphere_obj_121_20_1.3um' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 1000
optim_param['lambda_L1'] = 0
optim_param['lambda_TV'] = 0

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### load .mat file corresponding to object data###
object = loadmat(os.path.join('object_plane', object_file+'.mat'))['sphere_obj_121_20']
object = np.transpose(object, (0,2,3,1))

z_height = np.linspace(0,1,11)

### create directories to store object slice data (actual and estimated) if not previously created ###
obj_slices_dir = os.path.join('obj_slices', object_file)
est_obj_slices_dir = os.path.join('est_obj_slices', object_file)

if not os.path.exists(obj_slices_dir):
    os.mkdir(obj_slices_dir)

### plot and save figures of actual object z-slices ###
for i in range(len(z_height)):
    filename = os.path.join(obj_slices_dir, str(np.round(z_height[i], 3))+'.png')
    plot.plot_obj_voxel(object,'z',z_height[i], filename)


object = np.transpose(object, (0,3,1,2))

### model setup ###
model_cpu = cpu.smolm(psf, object.shape)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)

plot.plot_img(img, 'image')

### deconvolution ###
initial = np.random.rand(*object.shape)

obj_est, loss = gpu.estimate(psf, initial, img, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], device)
# print(loss)
img_est = model_cpu.forward(obj_est)

obj_est = np.transpose(obj_est, (0,2,3,1))

### plot and save figures of estimated object z-slices ###
for i in range(len(z_height)):
    filename = os.path.join(est_obj_slices_dir, str(np.round(z_height[i], 3))+'.png')
    plot.plot_obj_voxel(obj_est,'z',z_height[i], filename)

plot.plot_img(img_est, 'Reconstructed image')
