import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import torch
import torch.optim
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
import pyutils.GPU_module as gpu 
import pyutils.CPU_module as cpu
import pandas as pd 
from pyutils import m2_to_m1 as convert

FLAG_TRACE_LOSS = False
FLAG_SAVE_RESULT = True

### CUDA setup ###
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Using ' + device)

### document setup ###
psf_file = 'MVR_zf1000_3pixelsz_z10' # channel number x orientation x Z x X x Y
object_file_list = ['flat ring gamma 1','flat ring gamma 0.5','flat ring gamma 0',
                   'tilted ring gamma 1','tilted ring gamma 0.5','tilted ring gamma 0',
                   'fiber gamma 1','fiber gamma 0.5','fiber gamma 0',
                   'flat surface gamma 1','flat surface gamma 0.5','flat surface gamma 0',
                   'tilted surface gamma 1','tilted surface gamma 0.5','tilted surface gamma 0',
                   'hemisphere gamma 1','hemisphere gamma 0.5','hemisphere gamma 0'] # orientation x Z x X x Y
initial_file_list = ['flat ring initial','flat ring initial','flat ring initial',
                'tilted ring initial','tilted ring initial','tilted ring initial',
                'fiber initial','fiber initial','fiber initial',
                'flat surface initial','flat surface initial','flat surface initial',
                'tilted surface initial','tilted surface initial','tilted surface initial',
                'hemisphere initial','hemisphere initial','hemisphere initial']

num_of_trials = 5
factor = 10e5

### lambda setup ###
l1_lb = 5000
l1_up = 15000
l1_step_num = 11
tv_lb = 0
tv_up = 0
tv_step_num = 1

lambda_L1 = np.linspace(l1_lb,l1_up,l1_step_num,endpoint=True)
lambda_TV = np.linspace(tv_lb,tv_up,tv_step_num,endpoint=True)

### optimization setup ###
lr = 0.05
max_iter = 1500
lambda_I = 5

for idx in range (len(object_file_list)):

    object_file = object_file_list[idx]
    initial_file = initial_file_list[idx]

    ### psf, object_gt, image_raw setup ###
    psf = loadmat(os.path.join('psf', psf_file+'.mat'))['dsf']
    object_gt = loadmat(os.path.join('performance test objects new', object_file+'.mat'))['object']*factor
    object_initial = loadmat(os.path.join('performance test initials new', initial_file+'.mat'))['initial']*factor
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
                object_est, loss = gpu.estimate(psf, object_initial, 'dipole', image_noisy, 
                                                lr, max_iter, lambda_L1[idx_L1], lambda_TV[idx_TV], 
                                                lambda_I, device)
                # compare
                acc = np.mean((object_gt-object_est)**2)
                accuracy[trial,idx_L1,idx_TV] = acc

                loss_trace[i,...] = loss['total']
                i += 1

    if FLAG_SAVE_RESULT: 
        accuracy_mean = np.mean(accuracy, axis=0)
        accuracy_var = np.var(accuracy,axis=0)

        print(accuracy_mean)

        pd.DataFrame(accuracy_mean).to_csv(object_file + " mean " + str(l1_up) + ' ' + str(tv_up) + ".csv")
        pd.DataFrame(accuracy_var).to_csv(object_file + " var " + str(l1_up) + ' ' + str(tv_up) + ".csv")


    # trace the loss
    if FLAG_TRACE_LOSS: 
        for i in range (num_of_trials*l1_step_num*tv_step_num):
            plt.plot(loss_trace[i,...], linewidth=0.5)

        plt.title('trace the total loss when lr = ' + str(lr))
        plt.xlabel('iteration')
        plt.ylabel('total loss')
        plt.show()