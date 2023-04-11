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
psf_file = 'dsf_raMVR_psfSZ_61_zstack_20' # channel number x orientation x Z x X x Y
object_file = 'sphere_obj_121_20_1.3um' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 300
optim_param['lambda_L1'] = 0.000
optim_param['lambda_TV'] = 0.0003

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### object domain ###
object = loadmat(os.path.join('objects', object_file+'.mat'))['sphere_obj_121_20']

plot.plot_obj_voxel(object,'z',0.5)

### model setup ###
model_cpu = cpu.smolm(psf, object.shape)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)
plot.plot_img(img, 'Image of circle')

### initialization ###
# Now it's a function

initial = gpu.initialize(psf, object, img, device);

### deconvolution ###
obj_est, loss = gpu.estimate(psf, initial, img, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], device)

print(loss['total'][-1])
print(loss['lsq'][-1])
print(loss['L1'][-1])
print(loss['constraint_positive_mii'][-1])
print(loss['TV'][-1])

img_est = model_cpu.forward(obj_est)
plot.plot_obj_voxel(obj_est,'z',0.5,'Estimation of circle')
plot.plot_img(img_est, 'Reconstructed image of circle')
