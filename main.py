import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
import util.channel_2 as ch2
import util.obj_func as obj
import util.plot_func as plot

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'dsf_raMVR' # channel number x orientation x Z x X x Y
object_file = '' # orientation x Z x X x Y
image_file = ''

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.01
optim_param['max_iter'] = 1000
optim_param['lambda_L1'] = 0
optim_param['lambda_TV'] = 0

### PSF ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))[psf_file]

### model setup ###
model = ch2.pixOL(psf, device)

### object domain ###
object = obj.create_sphere(61,61,20,9,0.5)
object = np.transpose(object, (0,3,1,2))
# object = loadmat(os.path.join('object_plane', object_file+'.mat'))[object_file]

### image domain ###
if image_file == '':
    img = model(object)

plot.plot_img(img, 'image')

initial = np.random.rand(*object.shape)
img_est = model(initial)

### deconvolution ###
obj_est, loss = ch2.estimate(model, initial, img, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_L1'], optim_param['lambda_TV'])
# print(loss)
img_est = model(obj_est)

obj_est = np.transpose(obj_est, (0,2,3,1))

plot.plot_obj_voxel(obj_est,'z',0.5)

plot.plot_img(img_est, 'Reconstructed image')