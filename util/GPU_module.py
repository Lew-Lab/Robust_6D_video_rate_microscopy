import numpy as np
import torch
import torch.nn as nn
import torch.functional as F
import torch.fft as fft
import torchvision
import sys
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class smolm(nn.Module):
    
    def __init__(self, psf, obj_size, device):
        super(smolm, self).__init__()

        self.device = device # where the tensors saved

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
        # if self.psf_size['orientation'] != self.obj_size['orientation']:
        #     raise Exception('The number of orientation moments does not match')
        if self.psf_size['depth'] != self.obj_size['depth']:
            raise Exception('The number of z stack does not match')
        
        # padding 0 - assuming the psf size is smaller than the object
        w_diff = self.obj_size['width']-self.psf_size['width']
        h_diff = self.obj_size['height']-self.psf_size['height']

        w_left_pad = np.int16(np.floor(w_diff/2))
        w_right_pad = np.int16(np.ceil(w_diff/2))

        h_left_pad = np.int16(np.floor(h_diff/2))
        h_right_pad = np.int16(np.ceil(h_diff/2))

        psf_padded = nn.functional.pad(psf,(h_left_pad,h_right_pad,w_left_pad,w_right_pad),'constant',0)

        # store the psf and fft of the psf
        self.psf = torch.nn.Parameter(psf_padded, requires_grad=False)
        self.psf_fft = torch.nn.Parameter(fft.rfft2(fft.ifftshift(self.psf, (3,4))), requires_grad=False)

    def forward(self, obj):
        
        # simulate the image

        if (self.obj_size['orientation'] == 4):
            s = obj[0,...]
            mu_x = torch.cos(obj[1,...])*torch.sin(obj[2,...])
            mu_y = torch.cos(obj[1,...])*torch.sin(obj[2,...])
            mu_z = torch.cos(obj[2,...])
            wobble = obj[3,...]
            gamma = 1-3*wobble/(3*torch.pi)+wobble**2/(8*torch.pi**2)
            obj_6d = torch.zeros((6,self.obj_size['depth'],self.obj_size['width'],self.obj_size['height'])).to(self.device)
            obj_6d[0,...] = s*gamma*mu_x**2+(1-gamma)/3
            obj_6d[1,...] = s*gamma*mu_x**2+(1-gamma)/3
            obj_6d[2,...] = s*gamma*mu_x**2+(1-gamma)/3
            obj_6d[3,...] = s*gamma*mu_x*mu_y
            obj_6d[4,...] = s*gamma*mu_x*mu_z
            obj_6d[5,...] = s*gamma*mu_y*mu_z
        else:
            obj_6d = obj

        obj_fft = fft.rfft2(fft.ifftshift(obj_6d, (2,3)))
        img_fft = self.psf_fft * obj_fft
        img = fft.ifftshift(fft.irfft2(img_fft), (3,4))
        img = torch.sum(img, 1, keepdim=False)
        img = torch.sum(img, 1, keepdim=False)

        return img


def estimate(psf, obj, img_true, lr, max_iter, lambda_L1, lambda_TV, device):

    # convert variables from numpy to tensor
    obj = torch.tensor(obj, dtype=torch.float32, device=device, requires_grad=True).to(device)
    img_true = torch.tensor(img_true, dtype=torch.float32, device=device, requires_grad=False).to(device)
    psf = torch.tensor(psf, dtype=torch.float32, device=device, requires_grad=False).to(device)

    # run the model
    model = smolm(psf, obj.shape, device)

    # initialization
    arr_loss_lsq = torch.zeros(max_iter)
    arr_loss_L1 = torch.zeros(max_iter)
    arr_loss_TV = torch.zeros(max_iter)
    arr_loss_pos = torch.zeros(max_iter)
    arr_loss_total = torch.zeros(max_iter)

    # set up the optimizer
    optimizer = torch.optim.Adam([obj], lr)

    # least square loss function
    least_sq = nn.MSELoss(reduction='mean').type(torch.float32)

    for i in range (max_iter): 

        img_est = model(obj)

        # compute the loss function
        loss_lsq = least_sq(img_est, img_true)
        loss_L1 = loss_sparsity(obj)
        loss_TV = loss_smoothness(obj)
        loss_pos = constraint_positive_mii(obj)
        loss = loss_lsq + lambda_L1*loss_L1 + lambda_TV*loss_TV

        # store the loss
        arr_loss_lsq[i] = loss_lsq
        arr_loss_L1[i] = loss_L1
        arr_loss_TV[i] = loss_TV
        arr_loss_pos[i] = loss_pos
        arr_loss_total[i] = loss

        # gradient descent
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # end estimation
    obj = obj.detach().cpu().numpy()
    loss_final = dict()
    loss_final['total'] = arr_loss_total.detach().numpy()
    loss_final['lsq'] = arr_loss_lsq.detach().numpy()
    loss_final['L1'] = arr_loss_L1.detach().numpy()
    loss_final['constraint_positive_mii'] = arr_loss_pos.detach().numpy()
    loss_final['TV'] = arr_loss_TV.detach().numpy()

    return obj, loss_final

### regularization functions

def loss_sparsity(obj):
    return torch.mean(torch.abs(obj))

def loss_smoothness(obj):
    dfdz = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    dfdx = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    dfdy = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    tv = (torch.mean(dfdz)+torch.mean(dfdx)+torch.mean(dfdy))/3
    return tv

#### indicator functions

def constraint_positive_mii(obj):
    soft_indicator = torch.zeros(*obj[:3,...].shape)
    soft_indicator = 1/torch.exp(-obj[:3,...]**2)-1
    soft_indicator[torch.where(torch.greater_equal(obj[:3,...],0))] = 0
    return torch.sum(soft_indicator)

def constraint_positive_s(obj):
    s = obj[0,...]
    soft_indicator = 1/torch.exp(-s**2)-1
    soft_indicator[torch.where(torch.greater_equal(s,0))] = 0
    return torch.sum(soft_indicator)

def constraint_wobble(obj):
    wobble = obj[3,...]
    soft_indicator = torch.zeros(wobble)
    soft_indicator[torch.where(torch.greater_equal(wobble,2*torch.pi))] = 1/torch.exp(-(wobble-2*torch.pi)**2)-1
    soft_indicator[torch.where(torch.less_equal(wobble,0))] = 1/torch.exp(-(wobble)**2)-1
    return torch.sum(soft_indicator)

### help functions

def angle_to_second_moments(object_angle):

    obj_size = dict()
    obj_size['depth'] = object_angle.shape[1]
    obj_size['depth'] = object_angle.shape[2]
    obj_size['depth'] = object_angle.shape[3]

    s = object_angle(object_angle[0,...])
    mu_x = torch.cos(obj[1,...])*torch.sin(obj[2,...])
    mu_y = torch.cos(obj[1,...])*torch.sin(obj[2,...])
    mu_z = torch.cos(obj[2,...])
    wobble = obj[3,...]
    gamma = 1-3*wobble/(3*torch.pi)+wobble**2/(8*torch.pi**2)
    obj = torch.zeros((6,obj_size['depth'],obj_size['width'],obj_size['height']))
    obj[0,...] = s*gamma*mu_x**2+(1-gamma)/3
    obj[1,...] = s*gamma*mu_x**2+(1-gamma)/3
    obj[2,...] = s*gamma*mu_x**2+(1-gamma)/3
    obj[3,...] = s*gamma*mu_x*mu_y
    obj[4,...] = s*gamma*mu_x*mu_z
    obj[5,...] = s*gamma*mu_y*mu_z