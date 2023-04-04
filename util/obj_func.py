import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from scipy.io import loadmat, savemat

def create_sphere(obj_w, obj_h, obj_d, radius, sigma, orientation=True):

    center = [np.ceil(obj_w/2), np.ceil(obj_h/2), np.ceil(obj_d/2)]
    if obj_w<2*radius+1 or obj_h<2*radius+1 or obj_d<2*radius+1:
        raise Exception ('Radius of the sphere is too large.')
    
    x, y, z = np.meshgrid(np.linspace(1,obj_w,num=obj_w), np.linspace(1,obj_h,num=obj_h), np.linspace(1,obj_d,num=obj_d), indexing='ij')
    dist = np.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
    s = np.exp(-((dist-radius)/sigma)**2)

    if orientation==False:
        return s
    
    mu_x = (x-center[0])/(dist+1e-9)
    mu_y = (y-center[1])/(dist+1e-9)
    mu_z = (z-center[2])/(dist+1e-9)

    m_xx = mu_x**2
    m_yy = mu_y**2
    m_zz = mu_z**2
    m_xy = mu_x*mu_y
    m_xz = mu_x*mu_z
    m_yz = mu_y*mu_z

    sphere = np.stack((s*m_xx,s*m_yy,s*m_zz,s*m_xy,s*m_xz,s*m_yz), axis=0)

    return sphere

if __name__ == "__main__":
    object = create_sphere(61,61,61,20,0.5)
    object = np.transpose(object, (0,3,1,2))
    object_downsample = np.zeros((6,20,61,61))
    idx = 0
    for i in range (20):
        object_downsample[:,i,:,:] = object[:,idx,:,:]
        idx += 3

    object_pad = np.pad(object_downsample,((0,0),(0,0),(30,30),(30,30)),constant_values=0)

    savemat('sphere_obj_121_20.mat',{'sphere_obj_121_20': object_pad})