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
import pandas as pd 

FLAG_TRACE_LOSS = True
FLAG_SAVE_RESULT = True

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y
object_file_raw = 'flat surface gamma 0' # orientation x Z x X x Y
initial_file = 'flat surface initial'
num_of_trials = 5

### lambda setup ###
l1_lb = 0
l1_up = 20
l1_step_num = 5
tv_lb = 0
tv_up = 0.1
tv_step_num = 6

lambda_L1 = np.linspace(l1_lb,l1_up,l1_step_num,endpoint=True)
lambda_TV = np.linspace(tv_lb,tv_up,tv_step_num,endpoint=True)

### optimization setup ###
lr = 0.01
max_iter = 1500
lambda_I = 5

### psf, object_gt, image_raw setup ###
psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']
object_gt = loadmat(os.path.join('performance test objects', object_file_raw+'.mat'))['object']*10
object_initial = loadmat(os.path.join('performance test initials', initial_file+'.mat'))['initial']*10
object_size = (6,psf.shape[2],psf.shape[3],psf.shape[4])
object_iso_size = (1,psf.shape[2],psf.shape[3],psf.shape[4])
model_cpu = cpu.smolm(psf, object_size)
image_raw = model_cpu.forward(object_gt)

### initialize result ###
accuracy = np.zeros((num_of_trials,lambda_L1.shape[0],lambda_TV.shape[0]))
loss_trace = np.zeros((num_of_trials*l1_step_num*tv_step_num,max_iter))
i = 0

### start sweep ###
for trial in range (num_of_trials):
    image_noisy = np.random.poisson(image_raw)
    for idx_L1 in range (lambda_L1.shape[0]):
        for idx_TV in range (lambda_TV.shape[0]):
            object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, lr, max_iter, lambda_L1[idx_L1], lambda_TV[idx_TV], lambda_I, device)
            # compare
            acc = np.mean((object_gt-object_est)**2)
            accuracy[trial,idx_L1,idx_TV] = acc

            loss_trace[i,...] = loss['total']
            i += 1

if FLAG_SAVE_RESULT: 
    accuracy_mean = np.mean(accuracy, axis=0)
    accuracy_var = np.var(accuracy,axis=0)

    print(accuracy_mean)

    pd.DataFrame(accuracy_mean).to_csv(object_file_raw + " mean " + str(l1_up) + ' ' + str(tv_up) + ".csv")
    pd.DataFrame(accuracy_var).to_csv(object_file_raw + " var " + str(l1_up) + ' ' + str(tv_up) + ".csv")


# trace the loss
if FLAG_TRACE_LOSS: 
    for i in range (num_of_trials*l1_step_num*tv_step_num):
        plt.plot(loss_trace[i,...], linewidth=0.5)

    plt.title('trace the total loss when lr = ' + str(lr))
    plt.xlabel('iteration')
    plt.ylabel('total loss')
    plt.show()