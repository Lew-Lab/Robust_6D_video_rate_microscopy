import sys
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class smolm():
    
    def __init__(self, psf, obj_size):
        super(smolm, self).__init__()

        self.channel = psf.shape[0]

        # storing the size of the psf
        self.psf_size = dict()
        self.psf_size['orientation'] = psf.shape[1]
        self.psf_size['depth'] = psf.shape[2] # depth refers to z
        self.psf_size['width'] = psf.shape[3] # width refers to x
        self.psf_size['height'] = psf.shape[4] # height refers to y

        # storing the size of the object
        self.obj_size = dict()
        self.obj_size['orientation'] = obj_size[0]
        self.obj_size['depth'] = obj_size[1] # depth refers to z
        self.obj_size['width'] = obj_size[2] # width refers to x
        self.obj_size['height'] = obj_size[3] # height refers to y

        # dimension check
        if self.psf_size['orientation'] != self.obj_size['orientation']:
            raise Exception('The number of orientation moments does not match')
        if self.psf_size['depth'] != self.obj_size['depth']:
            raise Exception('The number of z stack does not match')
        
        ## PSF Padding Code is below ##
        #  padding 0 - assuming the psf size is smaller than the object
        
        w_diff = self.obj_size['width']-self.psf_size['width']
        h_diff = self.obj_size['height']-self.psf_size['height']

        w_left_pad = np.int16(np.floor(w_diff/2))
        w_right_pad = np.int16(np.ceil(w_diff/2))

        h_left_pad = np.int16(np.floor(h_diff/2))
        h_right_pad = np.int16(np.ceil(h_diff/2))

        psf_padded = np.pad(psf,((0,0),(0,0),(0,0),(w_left_pad,w_right_pad),(h_left_pad,h_right_pad)),constant_values=0)

        # store the psf and FFT of the psf
        self.psf = psf_padded
        self.psf_fft = np.fft.rfft2(np.fft.ifftshift(self.psf, (3,4)))

    def forward(self, obj):        
        # simulate the image - shift and then FFT
        obj_fft = np.fft.rfft2(np.fft.ifftshift(obj, (2,3))) # ifftshifting x and y so the result is centered
        img_fft = self.psf_fft * obj_fft # FFT(PSF) * FFT(Object)
        img = np.fft.ifftshift(np.fft.irfft2(img_fft), (3,4)) # ifftshifting x and y so the result is centered
        img_zstack = np.sum(img, axis=1) # sum the six second moments
        img = np.sum(img_zstack, axis=1) # sum along z

        return img, img_zstack