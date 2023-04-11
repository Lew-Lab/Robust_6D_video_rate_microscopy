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
psf_file = 'dsf_raMVR_psfSZ_89_zstack_1' # channel number x orientation x Z x X x Y
object_file = 'circle' # orientation x Z x X x Y
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
object = loadmat(os.path.join('objects', object_file+'.mat'))['object']
plot.plot_obj_voxel(object,'z',0.5)

### model setup ###
model_cpu = cpu.smolm(psf, object.shape)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)
plot.plot_img(img, 'Image of circle')

### initialization ###
# I will make it a function afterwards
psf_iso = np.sum(psf[:,:3,:,:,:],axis=1,keepdims=True)/3
object_iso = np.sum(object[:3,:,:,:],axis=0,keepdims=True)
plt.imshow(object_iso[0,0,...])
plt.colorbar()
plt.title('mxx+myy+mzz')
plt.show()
model_cpu_iso = cpu.smolm(psf_iso, object_iso.shape)
img_iso = model_cpu_iso.forward(object_iso)
plot.plot_img(img, '(Bxx+Byy+Bzz)/3 convolve with (mxx+myy+mzz)')
initial_iso = np.random.rand(*object_iso.shape)
obj_est_iso, loss_iso = gpu.estimate(psf_iso, initial_iso, img_iso, 0.01, 300, 0, 0, device)
img_est_iso = model_cpu_iso.forward(obj_est_iso)
plt.imshow(obj_est_iso[0,0,...])
plt.colorbar()
plt.title('Estimation of (mxx+myy+mzz)')
plt.show()
plot.plot_img(img_est_iso, 'Reconstructed image of (mxx+myy+mzz)')
initial = np.zeros(object.shape)
initial[0,:,:,:] = obj_est_iso
initial[1,:,:,:] = obj_est_iso
initial[2,:,:,:] = obj_est_iso

### deconvolution ###
obj_est, loss = gpu.estimate(psf, initial, img, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'], device)

print(loss)
img_est = model_cpu.forward(obj_est)
plot.plot_obj_voxel(obj_est,'z',0.5,'Estimation of circle')
plot.plot_img(img_est, 'Reconstructed image of circle')
