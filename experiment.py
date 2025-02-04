import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import numpy as np
from scipy.io import loadmat, savemat
import pyutils.GPU_module as gpu
import pyutils.CPU_module as cpu
import pyutils.plot_func as plot
from scipy import ndimage
<<<<<<< HEAD
=======
import pyutils.m2_to_m1 as convert

# starts with an image, performs deconvolution to estimate the signatures of an object
>>>>>>> main

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_dir = 'psf'
<<<<<<< HEAD
psf_file = 'MVR_zf800_pixelsz_z33' # #channel x orientation x Z x X x Y
image_dir = 'images'
image_file = 'OZ z1 NFP 800' # #channel x X x Y
=======
psf_file = 'MVR_zf800_pixelsz_z45' # #channel x orientation x Z x X x Y
Bstruct_file = 'B_ZScanned_zf800_101_n1000_2000'
image_dir = 'images'
image_file = 'OZdata' # #channel x X x Y
>>>>>>> main

### hyperparams setup ###
optim_param = dict()
optim_param['lr'] = 0.01
optim_param['max_iter'] = 1500
<<<<<<< HEAD
optim_param['lambda_L1'] = 5000
=======
optim_param['lambda_L1'] = 1000
>>>>>>> main
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 5

optim_param_iso = dict()
optim_param_iso['lr'] = 0.01
optim_param_iso['max_iter'] = 1500
<<<<<<< HEAD
optim_param_iso['lambda_L1'] = 5000
=======
optim_param_iso['lambda_L1'] = 1000
>>>>>>> main
optim_param_iso['lambda_TV'] = 0
optim_param_iso['lambda_I'] = 5

### file setup ###
psf = loadmat(os.path.join(psf_dir, psf_file+'.mat'))['dsf']
<<<<<<< HEAD
image = loadmat(os.path.join(image_dir, image_file+'.mat'))['image']
object_size = (6,psf.shape[2],image.shape[1],image.shape[2])
object_iso_size = (1,psf.shape[2],image.shape[1],image.shape[2])
=======


# psf = np.transpose(psf,(0,1,2,4,3))
# mxx = psf[:,1,:,:,:]
# myy = psf[:,0,:,:,:]
# mzz = psf[:,2,:,:,:]
# mxy = psf[:,4,:,:,:]
# mxz = psf[:,3,:,:,:]
# myz = psf[:,5,:,:,:]
# psf[:,0,:,:,:] = mxx
# psf[:,1,:,:,:] = myy
# psf[:,2,:,:,:] = mzz
# psf[:,3,:,:,:] = mxy
# psf[:,4,:,:,:] = mxz
# psf[:,5,:,:,:] = myz

image = loadmat(os.path.join(image_dir, image_file+'.mat'))['image']
# image = np.transpose(image,(0,2,1))
object_size = (6,psf.shape[2],image.shape[1],image.shape[2])
object_iso_size = (1,psf.shape[2],image.shape[1],image.shape[2])
Bstruct = loadmat(os.path.join('psf',Bstruct_file+'.mat'))['Bstruct']
>>>>>>> main

### model setup ###
model_cpu = cpu.smolm(psf, object_size)

<<<<<<< HEAD
plot.plot_img_tight(image, 'Image of ' + image_file)
=======
plot.plot_img_tighter(image, 'Image of ' + image_file)
>>>>>>> main

### initialization ###
initial = gpu.initialization(psf, object_size, object_iso_size, image, optim_param_iso['lr'], 
                                optim_param_iso['max_iter'], optim_param_iso['lambda_L1'], 
                                optim_param_iso['lambda_TV'], optim_param_iso['lambda_I'], device)
<<<<<<< HEAD
# plot.video_obj(initial,'Est initial of ' + image_file, 'initial Est')
=======
plot.video_obj(initial,'Est initial of ' + image_file, 'initial Est')
>>>>>>> main

### deconvolution ###
obj_est, loss = gpu.estimate(psf, initial, 'dipole', image, optim_param['lr'], 
                             optim_param['max_iter'], optim_param['lambda_L1'], 
                             optim_param['lambda_TV'], optim_param['lambda_I'], device)

<<<<<<< HEAD
### check the result ###
img_est = model_cpu.forward(obj_est)
plot.plot_img_tight(img_est, 'Reconstructed image of ' + image_file)
plot.video_obj(obj_est,'Est ' + image_file, 'object Est')
=======

### check the result ###
img_est = model_cpu.forward(obj_est)
plot.plot_img_tight(img_est, 'Reconstructed image of ' + image_file)
plot.video_obj(obj_est,'Est ' + image_file, 'object Est')
# s_est, m1_est, z_est, x_est, y_est = convert.m2_to_m1(obj_est,Bstruct,np.max(np.sum(obj_est[:3,...],axis=0))/10)

# plot.video_obj_1d(s_est, m1_est, z_est, x_est, y_est)
>>>>>>> main
