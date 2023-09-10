import numpy as np
import torch
import torch.nn as nn
import torch.fft as fft
import numpy as np
from pyutils import CPU_module as cpu

class smolm(nn.Module):
    
    def __init__(self, psf, obj_size, device):

        super(smolm, self).__init__()

        # where the tensors are saved
        self.device = device

        self.channel = psf.shape[0]

        # storing the size of the psf
        self.psf_size = dict()
        self.psf_size['orientation'] = psf.shape[1]
        self.psf_size['depth'] = psf.shape[2]
        self.psf_size['width'] = psf.shape[3]
        self.psf_size['height'] = psf.shape[4]

        # storing the size of the object
        self.obj_size = dict()
        self.obj_size['orientation'] = obj_size[0]
        self.obj_size['depth'] = obj_size[1]
        self.obj_size['width'] = obj_size[2]
        self.obj_size['height'] = obj_size[3]

        # dimension check
        if self.psf_size['orientation'] != self.obj_size['orientation']:
            raise Exception('The number of orientation moments does not match')
        if self.psf_size['depth'] != self.obj_size['depth']:
            raise Exception('The number of z stack does not match')
        
        # padding 0 - assuming the psf size is smaller than the object
        w_diff = self.obj_size['width']-self.psf_size['width']
        h_diff = self.obj_size['height']-self.psf_size['height']

        w_left_pad = np.int16(np.floor(w_diff/2))
        w_right_pad = np.int16(np.ceil(w_diff/2))

        h_left_pad = np.int16(np.floor(h_diff/2))
        h_right_pad = np.int16(np.ceil(h_diff/2))

        psf_padded = nn.functional.pad(psf,(h_left_pad,h_right_pad,
                                            w_left_pad,w_right_pad),'constant',0)

        # store the psf and fft of the psf
        self.psf = torch.nn.Parameter(psf_padded, requires_grad=False)
        self.psf_fft = torch.nn.Parameter(fft.rfft2(fft.ifftshift
                                                    (self.psf, (3,4))), requires_grad=False)

    def forward(self, obj):
        
        # simulate the image
        shape_img = (self.obj_size['width'],self.obj_size['height'])
        obj_fft = fft.rfft2(fft.ifftshift(obj, (2,3)))
        img_fft = self.psf_fft * obj_fft
        img = fft.ifftshift(fft.irfft2(img_fft,s=shape_img), (3,4))
        img = torch.sum(img, 1, keepdim=False)
        img = torch.sum(img, 1, keepdim=False)

        return img

def initialization(psf, obj_size, obj_iso_size, img_true, lr, max_iter, lambda_l1, lambda_tv, lambda_I, device):

    ''' The initialization function estimates the shape of the object.

    Input:
    psf: #channel x 6 x Z x X x Y. Point spread function of the optical system.
    obj: 6 x Z x X x Y. Initialized object.
    img_true: #channel x X x Y. Image of the object.
    lr: learning rate. 
    max_iter: maximum iteration.
    lambda_L1: hyperparameter for L1 regularizer. 
    lambda_TV: hyperparameter for total variance regularizer.
    Lambda_I: hyperparameter for indicator function. 
    device: 'cuda' or 'cpu'. 

    Output:
    initial: 6 x Z x X x Y. Isotropic initial object. '''

    psf_iso = np.sum(psf[:,:3,:,:,:],axis=1,keepdims=True)
    model_cpu_iso = cpu.smolm(psf_iso, obj_iso_size)

    initial_iso = np.random.rand(1,psf.shape[2],psf.shape[3],psf.shape[4])
    obj_est_iso, loss_iso = estimate(psf_iso, initial_iso, 's', img_true, 
                                     lr, max_iter, lambda_l1, lambda_tv, lambda_I, device)
    img_est_iso = model_cpu_iso.forward(obj_est_iso)
    obj_est_iso = obj_est_iso*(np.sum(img_true)/np.sum(img_est_iso))
    img_est_iso = model_cpu_iso.forward(obj_est_iso)

    initial = np.zeros(obj_size)
    initial[0,:,:,:] = obj_est_iso/3
    initial[1,:,:,:] = obj_est_iso/3
    initial[2,:,:,:] = obj_est_iso/3

    return initial

def estimate(psf, obj, type, img_true, lr, max_iter, lambda_L1, lambda_TV, lambda_I, device):

    ''' The estimator function that solves the deconvolution problem.

    Input: 
    psf: #channel x 6 (or 1) x Z x X x Y. Point spread function of the optical system.
    obj: 6 (or 1) x Z x X x Y. Initialized object.
    type: 'dipole' or 's'.
    img_true: #channel x X x Y. Image of the object.
    lr: learning rate. 
    max_iter: maximum iteration.
    lambda_L1: hyperparameter for L1 regularizer. 
    lambda_TV: hyperparameter for total variance regularizer.
    Lambda_I: hyperparameter for indicator function. 
    device: 'cuda' or 'cpu'. 

    Output: 
    loss: cost values for each term 'lsq', 'l1', 'tv', 'pos', 'total'. '''

    # convert variables from numpy to tensor
    obj = torch.tensor(obj, dtype=torch.float32, device=device, requires_grad=True)
    img_true = torch.tensor(img_true, dtype=torch.float32, device=device, requires_grad=False)
    psf = torch.tensor(psf, dtype=torch.float32, device=device, requires_grad=False)

    # run the model
    model = smolm(psf, obj.shape, device)

    # initialization
    arr_loss_lsq = torch.zeros(max_iter)
    arr_loss_L1 = torch.zeros(max_iter)
    arr_loss_TV = torch.zeros(max_iter)
    arr_loss_pos = torch.zeros(max_iter)
    arr_loss_total = torch.zeros(max_iter)

    # set up the optimizer
    optimizer = torch.optim.SGD([obj], lr)

    # least square loss function
    least_sq = nn.MSELoss(reduction='sum').type(torch.float32)

    for i in range (max_iter): 

        img_est = model(obj)

        # compute the loss function
        loss_lsq = least_sq(img_est, img_true)
        loss_L1 = loss_sparsity(obj,type)
        loss_TV = loss_smoothness(obj)
        loss_pos = constraint_positive_mii(obj,device)
        loss = loss_lsq + lambda_L1*loss_L1 + lambda_TV*loss_TV + lambda_I*loss_pos

        # store the loss
        arr_loss_lsq[i] = loss_lsq
        arr_loss_L1[i] = loss_L1
        arr_loss_TV[i] = loss_TV
        arr_loss_pos[i] = loss_pos
        arr_loss_total[i] = loss

        # gradient descent
        optimizer.zero_grad()
        loss.backward()
        # nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5, norm_type=2)
        optimizer.step()

    # end estimation
    obj = obj.detach().cpu().numpy()
    loss_final = dict()
    loss_final['total'] = arr_loss_total.detach().numpy()
    loss_final['lsq'] = arr_loss_lsq.detach().numpy()
    loss_final['l1'] = arr_loss_L1.detach().numpy()
    loss_final['i'] = arr_loss_pos.detach().numpy()
    loss_final['tv'] = arr_loss_TV.detach().numpy()

    return obj, loss_final

# loss function

def loss_sparsity(obj,type):
    if type == 's':
        return torch.sum(torch.abs(obj))
    elif type == 'dipole':
        return torch.sum(torch.sqrt(torch.sum(obj**2,dim=0)))
        # return torch.sum(torch.abs(obj))

def loss_smoothness(obj):
    if obj.shape[1] == 1:
        dfdz = 0
    else:
        dfdz = torch.sum(torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])**2)
    dfdx = torch.sum(torch.abs(obj[:,:,1:,:]-obj[:,:,0:-1,:])**2)
    dfdy = torch.sum(torch.abs(obj[:,:,:,1:]-obj[:,:,:,0:-1])**2)
    tv = (dfdz+dfdx+dfdy)
    return tv
    
def constraint_positive_mii(obj,device):
    soft_indicator = torch.zeros(*obj[:3,...].shape).to(device=device)
    soft_indicator = torch.max(soft_indicator,torch.abs(obj[:3,...]))**2
    soft_indicator[torch.where(torch.greater_equal(obj[:3,...],0))] = 0
    return torch.sum(soft_indicator)
