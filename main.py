import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:1024"

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
psf_file = 'dsf_raMVR_61_51' # channel number x orientation x Z x X x Y
object_file = 'sphere_obj_101_51' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 1000
optim_param['lambda_L1'] = 0
optim_param['lambda_TV'] = 0

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))[psf_file]

### object domain ###
object = loadmat(os.path.join('object_plane', object_file+'.mat'))[object_file]
object = np.transpose(object, (0,2,3,1))
plot.plot_obj_voxel(object,'z',0.5)
object = np.transpose(object, (0,3,1,2))

### model setup ###
model_cpu = cpu.smolm(psf, object.shape)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)

plot.plot_img(img, 'image')

### deconvolution ###
# initial = np.random.rand(*object.shape)
initial = np.random.rand(4,object.shape[1],object.shape[2],object.shape[3])

obj_est, loss = gpu.estimate(psf, initial, img, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], device)
# print(loss)

s = obj_est[0,...]
mu_x = np.cos(obj_est[1,...])*np.sin(obj_est[2,...])
mu_y = np.cos(obj_est[1,...])*np.sin(obj_est[2,...])
mu_z = np.cos(obj_est[2,...])
wobble = obj_est[3,...]
gamma = 1-3*wobble/(3*np.pi)+wobble**2/(8*np.pi**2)
obj_6d = np.zeros((6,obj_est.shape[1],obj_est.shape[2],obj_est.shape[3]))
obj_6d[0,...] = s*gamma*mu_x**2+(1-gamma)/3
obj_6d[1,...] = s*gamma*mu_x**2+(1-gamma)/3
obj_6d[2,...] = s*gamma*mu_x**2+(1-gamma)/3
obj_6d[3,...] = s*gamma*mu_x*mu_y
obj_6d[4,...] = s*gamma*mu_x*mu_z
obj_6d[5,...] = s*gamma*mu_y*mu_z

# img_est = model_cpu.forward(obj_est)
img_est = model_cpu.forward(obj_6d)

obj_est = np.transpose(obj_6d, (0,2,3,1))
plot.plot_obj_voxel(obj_est,'z',0.5)
plot.plot_obj_voxel(obj_est,'z',0.2)

plot.plot_img(img_est, 'Reconstructed image')