import torch
import torch.nn as nn
import torch.fft as fft
import torchvision
import sys
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class smolm(nn.Module):
    
    def __init__(self, psf, device):
        super(smolm, self).__init__()

        self.device = device # where the tensors saved
        self.run_estimation = False # it's on when running estimation otherwise off

        self.channel = psf.shape[0]

        # storing the size of the psf
        self.psf_size = dict()
        self.psf_size['orientation'] = psf.shape[1]
        self.psf_size['depth'] = psf.shape[2]
        self.psf_size['width'] = psf.shape[3]
        self.psf_size['height'] = psf.shape[4]

        # store the psf
        self.psf = torch.nn.Parameter(torch.tensor(psf, dtype=torch.float32, device=self.device), requires_grad=False)

    def forward(self, obj):
        
        # storing the size of the object
        obj_size = dict()
        obj_size['orientation'] = obj.shape[0]
        obj_size['depth'] = obj.shape[1]
        obj_size['width'] = obj.shape[2]
        obj_size['height'] = obj.shape[3]

        # we can change this later
        if obj_size != self.psf_size:
            raise Exception('Dimension of PSF and object is different.')

            # max_depth = np.max([obj_size['depth'], self.psf_size['depth']])
            # max_width = np.max([obj_size['width'], self.psf_size['width']])
            # max_height = np.max([obj_size['height'], self.psf_size['height']])

            # # match the dimension of the psf and obj by padding 0
            # obj_padded = nn.functional.pad(obj, (0,max_depth-obj_size['depth'],max_width-obj_size['width'],max_height-obj_size['height']), 'const', 0)
            # psf_padded = nn.functional.pad(self.psf, (0,max_depth-self.psf_size['depth'],max_width-self.psf_size['width'],max_height-self.psf_size['height']), 'const', 0)
        
        if not torch.is_tensor(obj): 
            obj = torch.tensor(obj, dtype=torch.float32, device=self.device, requires_grad=True)
        
        # simulate the image
        psf_fft = fft.rfft2(fft.ifftshift(self.psf, (3,4)))
        obj_fft = fft.rfft2(fft.ifftshift(obj, (2,3)))
        img_fft = psf_fft * obj_fft
        img = fft.ifftshift(fft.irfft2(img_fft), (3,4))
        img = torch.sum(img, 1, keepdim=False)
        img = torch.sum(img, 1, keepdim=False)

        # not converting img to numpy array if we are running estimation
        if not self.run_estimation:
            img = img.detach().cpu().numpy()

        return img
    

def estimate(model, obj, img_true, lr, max_iter, lambda_L1, lambda_TV):

    model.run_estimation = True

    # initialization
    obj = torch.tensor(obj, dtype=torch.float32, device=model.device, requires_grad=True)
    img_true = torch.tensor(img_true, dtype=torch.float32, device=model.device, requires_grad=False)
    arr_loss_lsq = torch.zeros(max_iter)
    arr_loss_L1 = torch.zeros(max_iter)
    arr_loss_TV = torch.zeros(max_iter)
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
        loss = loss_lsq + lambda_L1*loss_L1 + lambda_TV*loss_TV

        # store the loss
        arr_loss_lsq[i] = loss_lsq
        arr_loss_L1[i] = loss_L1
        arr_loss_TV[i] = loss_TV
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
    loss_final['TV'] = arr_loss_TV.detach().numpy()

    model.run_estimation = False

    return obj, loss_final

def loss_sparsity(obj):
    return torch.mean(torch.abs(obj))

def loss_smoothness(obj):
    dfdz = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    dfdx = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    dfdy = torch.abs(obj[:,1:,:,:]-obj[:,0:-1,:,:])
    tv = (torch.mean(dfdz)+torch.mean(dfdx)+torch.mean(dfdy))/3
    return tv
    