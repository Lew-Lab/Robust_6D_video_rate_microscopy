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
import matplotlib.animation as animation
import matplotlib.cm as cm

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y
object_file = 'tilted ring cubic' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.05
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 0.08
optim_param['lambda_TV'] = 0
optim_param['lambda_I'] = 1

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### object domain ###
object = loadmat(os.path.join('object_plane', object_file+'.mat'))['object']*100
# plot.plot_obj_voxel(object,'z',0.5,'GT')
# plot.video_obj(object,'tilted ring cubic','object GT')

object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])

### model setup ###
model_cpu = cpu.smolm(psf, object_size)

### image domain ###
if image_file == '':
    img = model_cpu.forward(object)
else:
    img = loadmat(os.path.join('image_plane', image_file+'.mat'))['image']

# print(np.sum(img))
# plot.plot_img_tight(img, 'Image of tilted ring')

### initialization ###
# I will make it a function afterwards
psf_iso = np.sum(psf[:,:3,:,:,:],axis=1,keepdims=True)
object_iso = np.sum(object[:3,:,:,:],axis=0,keepdims=True)
model_cpu_iso = cpu.smolm(psf_iso, object_iso_size)
img_iso = model_cpu_iso.forward(object_iso)
# plot.plot_img_tight(img_iso, 'Image of (mxx+myy+mzz)')
plot.video_obj(object_iso,'GT tilted ring initial','initial GT')

initial_iso = np.random.rand(1,psf.shape[2],psf.shape[3],psf.shape[4])
obj_est_iso, loss_iso = gpu.estimate(psf_iso, initial_iso, 's', img, 0.01, 500, 1, 0, 0, 7, device)
img_est_iso = model_cpu_iso.forward(obj_est_iso)
plot.plot_img_tight(img_est_iso, 'Reconstructed image of (mxx+myy+mzz) before scaling')
obj_est_iso = obj_est_iso*(np.sum(img)/np.sum(img_est_iso))
img_est_iso = model_cpu_iso.forward(obj_est_iso)
plot.plot_img_tight(img_est_iso, 'Reconstructed scaled image of (mxx+myy+mzz) after scaling')
plot.video_obj(obj_est_iso,'Est tilted ring initial','initial Est')
initial = np.zeros(object_size)
initial[0,:,:,:] = obj_est_iso/3
initial[1,:,:,:] = obj_est_iso/3
initial[2,:,:,:] = obj_est_iso/3

# ### deconvolution ###
# obj_est, loss = gpu.estimate(psf, initial, 'dipole', img, optim_param['learning_rate'], optim_param['max_iter'], 1, optim_param['lambda_L1'], optim_param['lambda_TV'], optim_param['lambda_I'], device)

# # print(loss)
# plt.plot(loss['lsq'])
# plt.xlabel('iteration')
# plt.ylabel('lsq')
# plt.title('Last epoch lsq:' + str(loss['lsq'][-1]))
# plt.show()

# plt.plot(loss['total'])
# plt.xlabel('iteration')
# plt.ylabel('loss')
# plt.title('Last epoch total:' + str(loss['total'][-1]))
# plt.show()

# img_est = model_cpu.forward(obj_est)
# # plot.plot_img_tight(img_est, 'Reconstructed image of tilted ring')

# # plot.video_obj(obj_est,'Est tilted ring','object Est')