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
object_file_list = ['tilted surface'] # orientation x Z x X x Y
initial_file_list = ['initialization tilted ring']
image_file = ''
num_of_trials = 2

psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']

### hyperparams setup ###
optim_param = dict()
optim_param['learning_rate'] = 0.05
optim_param['max_iter'] = 1500
optim_param['lambda_L1'] = 0.5
optim_param['lambda_TV'] = 0.1
optim_param['lambda_I'] = 5

accuracy = np.zeros((len(object_file_list),num_of_trials))

for idx_obj in range (len(object_file_list)):
    object_gt = loadmat(os.path.join('object_plane', object_file_list[idx_obj]+'.mat'))['object']
    object_initial = loadmat(os.path.join('object_plane', initial_file_list[idx_obj]+'.mat'))['initial']
    object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
    object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
    model_cpu = cpu.smolm(psf, object_size)
    image_raw = model_cpu.forward(object_gt)
    for trial in range (num_of_trials):
        image_noisy = np.random.poisson(image_raw)
        object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, optim_param['learning_rate'], optim_param['max_iter'], optim_param['lambda_TV'], optim_param['lambda_TV'], optim_param['lambda_I'], device)
        # compare
        acc = np.mean((object_gt-object_est)**2)
        accuracy[idx_obj,trial] = acc

print(accuracy)