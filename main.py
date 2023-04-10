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
psf_file = 'planar_pixOL' # channel number x orientation x Z x X x Y
object_file = 'circle' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 300
optim_param['lambda_L1'] = 0.035
optim_param['lambda_TV'] = 0.015

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf_raMVR_61_20']

### load .mat file corresponding to object data###
object = loadmat(os.path.join('objects', object_file+'.mat'))['sphere_obj_61_20']
object = np.transpose(object, (0,2,3,1))

z_height = np.linspace(0,0.95,20)
print(z_height)
print(len(z_height)-1)

### create directories to store object & object image slice data (actual and estimated) if not previously created ###
obj_slices_dir = os.path.join('obj_slices', object_file)
obj_img_slices_dir = os.path.join('obj_img_slices', object_file)

est_obj_slices_dir = os.path.join('est_obj_slices', object_file)
est_obj_img_slices_dir = os.path.join('est_obj_img_slices', object_file)

if not os.path.exists(obj_slices_dir):
    os.mkdir(obj_slices_dir)

if not os.path.exists(obj_img_slices_dir):
    os.mkdir(obj_img_slices_dir)

if not os.path.exists(est_obj_slices_dir):
    os.mkdir(est_obj_slices_dir)

if not os.path.exists(est_obj_img_slices_dir):
    os.mkdir(est_obj_img_slices_dir)

### plot and save figures of actual object z-slices ###
for i in range(len(z_height)-1):
    filename = os.path.join(obj_slices_dir, str(np.round(z_height[i], 3))+'.png')
    plot.plot_obj_voxel(object,'z',z_height[i], filename)


object = np.transpose(object, (0,3,1,2))

### model setup ###
model_cpu = cpu.smolm(psf, object.shape)

### Generate an object image using ground truth data and save the object z-slices. ###
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
# print(loss['total'])

### Generate an object image using estimated sphere data and save the object z-slices. ###
img_est, img_est_zstack = model_cpu.forward(obj_est)

obj_est = np.transpose(obj_est, (0,2,3,1))

### plot and save figures of estimated object z-slices ###
for i in range(len(z_height)-1):
    filename = os.path.join(est_obj_slices_dir, str(np.round(z_height[i], 3))+'.png')
    plot.plot_obj_voxel(obj_est,'z',z_height[i], filename)

# Plot estimated object image in 8 MVR channels.
filename = os.path.join(est_obj_slices_dir, 'Reconstructed object image.png')
plot.plot_img(img_est, 'Reconstructed object image', filename)

# Generate and save z-stack of estimated object image in 8 MVR channels.
plot.zstack_video(img_est_zstack,est_obj_img_slices_dir)


print(loss)
img_est = model_cpu.forward(obj_est)
plot.plot_obj_voxel(obj_est,'z',0.5,'Estimation of circle')
plot.plot_img(img_est, 'Reconstructed image of circle')
