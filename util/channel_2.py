import torch
import torch.nn as nn
import torch.fft as fft
import torchvision
import sys
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class pixOL(nn.Module):
    
    def __init__(self, psf, device):
        super(pixOL, self).__init__()

        self.device = device
        # setting up the PSF size
        self.channel = psf.shape[0]
        # setting up the PSF size
        self.psf_size['orientation'] = psf.shape[1]
        self.psf_size['depth'] = psf.shape[2]
        self.psf_size['width'] = psf.shape[3]
        self.psf_size['height'] = psf.shape[4]
        self.psf = torch.nn.Parameter(torch.tensor(psf, dtype=torch.float32, device=self.device), requires_grad=False)

    def forward(self, obj):
        
        obj_size = dict()
        obj_size['orientation'] = obj.shape[0]
        obj_size['depth'] = obj.shape[1]
        obj_size['width'] = obj.shape[2]
        obj_size['height'] = obj.shape[3]

        # we can change this later
        if obj_size != self.psf_size:
            raise Exception('Dimension of PSF and object is different.')
        
        obj = torch.tensor(obj, dtype=torch.float32, device=self.device)
        
        # FFT: dimension need to be checked 
        psf_fft = fft.fft2(fft.ifftshift(self.psf, (3,4)))
        obj_fft = fft.fft2(fft.ifftshift(obj, (3,4)))
        img_fft = psf_fft * obj_fft
        img = torch.real(fft.ifft2(img_fft))
        img = torch.sum(img, 2, keepdim=False)
        img = torch.sum(img, 1, keepdim=False)

        return img
    

def estimate(model, obj, img_ture, lr, max_iter, lambda_L1, lambda_TV):
    optimizer = torch.optim.Adam(obj, lr)
    
    # initialization
    obj = torch.tensor(obj, dtype=torch.float32, device=model.device)

    loss_lsq = []
    loss_L1 = []
    loss_TV = []
    loss_total = []

    for i in range (max_iter): 

        optimizer.zero_grad()
        img_est = model(obj)

        # compute the loss function
        np.append(loss_lsq, nn.MSELoss(img_est, img_ture, reduction='sum'))
        np.append(loss_L1, nn.L1Loss(obj))
        np.append(loss_TV, nn.SmoothL1Loss(obj))
        loss = loss_lsq[-1] + lambda_L1*loss_L1[-1] + lambda_TV*loss_TV[-1]
        np.append(loss_total, loss)

        # gradient descent
        loss.backward()
        optimizer.step()

    loss_final = dict()
    loss_final['total'] = loss_total[-1]
    loss_final['lsq'] = loss_lsq[-1]
    loss_final['L1'] = loss_L1[-1]
    loss_final['TV'] = loss_TV[-1]

    return obj, loss_final