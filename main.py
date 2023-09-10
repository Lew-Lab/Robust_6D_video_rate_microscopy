import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import numpy as np
from scipy.io import loadmat, savemat
import pyutils.GPU_module as gpu
import pyutils.CPU_module as cpu
import pyutils.plot_func as plot

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'MVR_zf500_pixelsz_z23' # #channel x orientation x Z x X x Y
object_file = 'simulated lipid membrane' # orientation x Z x X x Y
initial_file = '' # orientation x Z x X x Y
image_file = '' # #channel x X x Y
factor = 1
FLAG_NOISE = True

### hyperparams setup ###
optim_param = dict()
optim_param['lr'] = 0.05
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 50
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 5

optim_param_iso = dict()
optim_param_iso['lr'] = 0.01
optim_param_iso['max_iter'] = 1000
optim_param_iso['lambda_L1'] = 100
optim_param_iso['lambda_TV'] = 0
optim_param_iso['lambda_I'] = 5

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### object domain ###
object = loadmat(os.path.join('performance test objects new', object_file+'.mat'))['object']*factor
object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
# plot.video_obj(object,'GT of ' + object_file,'object GT')

### model setup ###
model_cpu = cpu.smolm(psf, object_size)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)
else:
    img = loadmat(os.path.join('image_plane', image_file+'.mat'))['image']

if FLAG_NOISE and image_file == '':
    img = np.random.poisson(img)

plot.plot_img_tight(img, 'Image of ' + object_file)

# ### initialization ###
if initial_file == '': 
    initial = gpu.initialization(psf, object_size, object_iso_size, img, optim_param_iso['lr'], 
                                 optim_param_iso['max_iter'], optim_param_iso['lambda_L1'], 
                                 optim_param_iso['lambda_TV'], optim_param_iso['lambda_I'], device)
    savemat('simulated lipid membrane initial.mat',{'initial':initial})
    # plot.video_obj(initial,'Est initial of ' + object_file, 'initial Est')
else: 
    initial = loadmat(os.path.join('performance test initials', initial_file +'.mat'))['initial']*factor
    # plot.video_obj(initial,'Est initial of ' + object_file, 'initial Est')

# ### deconvolution ###
# obj_est, loss = gpu.estimate(psf, initial, 'dipole', img, optim_param['lr'], 
#                              optim_param['max_iter'], optim_param['lambda_L1'], 
#                              optim_param['lambda_TV'], optim_param['lambda_I'], device)

# ### check the result ###
# img_est = model_cpu.forward(obj_est)
# plot.plot_img_tight(img_est, 'Reconstructed image of ' + object_file)
# plot.video_obj(obj_est,'Est ' + object_file,'object Est')