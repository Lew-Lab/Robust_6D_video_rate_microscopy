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

def create_sphere_fiber(obj_w, obj_h, obj_d, radius, sigma, gamma, orientation=True):

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

    # theta = np.arccos(mu_z)
    # phi = np.sign(mu_y)*np.arccos(mu_x/(np.sqrt(mu_x**2+mu_y**2)+1e-9))
    # phi[np.isnan(phi) == True] = 0

    # theta = theta-np.pi/2

    # mu_x = np.sin(theta)*np.cos(phi)
    # mu_y = np.sin(theta)*np.sin(phi)
    # mu_z = np.cos(theta)

    m_xx = gamma*mu_x**2 + (1-gamma)/3
    m_yy = gamma*mu_y**2 + (1-gamma)/3
    m_zz = gamma*mu_z**2 + (1-gamma)/3
    m_xy = gamma*mu_x*mu_y
    m_xz = gamma*mu_x*mu_z
    m_yz = gamma*mu_y*mu_z

    sphere = np.stack((s*m_xx,s*m_yy,s*m_zz,s*m_xy,s*m_xz,s*m_yz), axis=0)

    return np.transpose(sphere, (0,3,1,2))


if __name__ == "__main__":

    # object = create_sphere_fiber(61,61,61,6,0.7,0)

    # tilted_ring_temp = np.zeros((6,61,61,61))
    # for i in range (61):
    #     tilted_ring_temp[:,i,i,:] = object[:,i,i,:]

    # tilted_ring = tilted_ring_temp[:,25:35,...]
    # savemat('tilted ring gamma 0.mat', {'object': tilted_ring})

    # flat_ring = np.zeros((6,10,61,61))
    # flat_ring[:,5,...] = object[:,30,...]
    # savemat('flat ring gamma 0.mat', {'object': flat_ring})





    # object[4,:,0:30,30] = -object[4,:,0:30,30]

    # real = np.zeros((6,10,61,61))
    # real[:,1:9,...] = object[:,30:38,...]
    # real[:,:,0:29,0:29] = 0
    # real[:,:,0:29,32:] = 0
    # real[:,:,32:,0:29] = 0
    # real[:,:,32:,32:] = 0
    # savemat('hemisphere gamma 0.mat', {'object': 100*real})

    # tilted_ring_big = np.zeros((6,10,61,61))
    # idx = 25
    # for i in range(10):
    #     if i == 0 :
    #         temp = np.sum(tilted_ring_temp[:,idx-2:idx+3,...], axis=1)
    #     else:
    #         temp = np.sum(tilted_ring_temp[:,idx-2:idx+3,...], axis=1)
    #     tilted_ring_big[:,i,...] = np.sum(tilted_ring_temp[:,idx-2:idx+3,...], axis=1)
    #     idx += 1

    # sphere = np.zeros((6,10,61,61))
    # sphere[:,0:9,...] = object[:,26:35,:,:]

    # savemat('sphere 2.mat',{'object': sphere})
    
    # and plot everything
    # ax = plt.figure().add_subplot(projection='3d')
    # ax.voxels(np.linspace(1,61,num=61), np.linspace(1,61,num=61), np.linspace(1,61,num=61), tilted_ring[0,...],
    #         linewidth=0.5)
    # ax.set(xlabel='r', ylabel='g', zlabel='b')
    # ax.set_aspect('equal')

    # plt.show()
    # object_downsample = np.zeros((6,1,61,61))

    # object_downsample[:,0,:,:] = object[:,30,:,:]
    # idx = 0
    # for i in range (23):
    #     object_downsample[:,i,:,:] = object[:,idx,:,:]
    #     idx += 4

    # savemat('test.mat',{'object': object_downsample})


    ### tilted surface ###

    # layer_xz = np.zeros((6,10,61))

    # layer_xz_mux = np.zeros((1,10,61))
    # layer_xz_muy = np.zeros((1,10,61))
    # layer_xz_muz = np.zeros((1,10,61))
    # layer_xz_s = np.zeros((1,10,61))

    # for i in range (1, 8):
    #     layer_xz_mux[0,i:i+2,25+i-1:25+i+2] = 1/np.sqrt(2)
    #     layer_xz_muy[0,i:i+2,25+i] = 0
    #     layer_xz_muz[0,i:i+2,25+i-1:25+i+2] = -1/np.sqrt(2)
    #     layer_xz_s[0,i:i+2,25+i] = 1
    #     layer_xz_s[0,i,25+i+1] = 0.1299
    #     layer_xz_s[0,i+1,25+i-1] = 0.1299

    # layer_xz_s[0,1,25] = 0.1299
    # layer_xz_s[0,8,25+8] = 0.1299

    # layer_xz_mux_dim = np.zeros((1,10,61))
    # layer_xz_muy_dim = np.zeros((1,10,61))
    # layer_xz_muz_dim = np.zeros((1,10,61))
    # layer_xz_s_dim = np.zeros((1,10,61))

    # for i in range(1,8):
    #     layer_xz_mux_dim[0,i:i+2,25+i-1:25+i+2] = 1/np.sqrt(2)
    #     layer_xz_muy_dim[0,i:i+2,25+i] = 0
    #     layer_xz_muz_dim[0,i:i+2,25+i-1:25+i+2] = -1/np.sqrt(2)
    #     layer_xz_s_dim[0,i:i+2,25+i] = 0.1299
    #     layer_xz_s_dim[0,i,25+i+1] = 0.0169
    #     layer_xz_s_dim[0,i+1,25+i-1] = 0.0169

    # mux = np.zeros((10,61,61))
    # muy = np.zeros((10,61,61))
    # muz = np.zeros((10,61,61))
    # s = np.zeros((10,61,61))

    # mux[:,:,24] = layer_xz_mux_dim
    # muy[:,:,24] = layer_xz_muy_dim
    # muz[:,:,24] = layer_xz_muz_dim
    # s[:,:,24] = layer_xz_s_dim
    # mux[:,:,35] = layer_xz_mux_dim
    # muy[:,:,35] = layer_xz_muy_dim
    # muz[:,:,35] = layer_xz_muz_dim
    # s[:,:,35] = layer_xz_s_dim

    # for i in range (25, 35):
    #     mux[:,:,i] = layer_xz_mux
    #     muz[:,:,i] = layer_xz_muz
    #     s[:,:,i] = layer_xz_s

    # gamma = 0

    # tilted_surface = np.zeros((6,10,61,61))

    # tilted_surface[0,...] = s*(gamma*mux**2 + (1-gamma)/3)
    # tilted_surface[1,...] = s*(gamma*muy**2 + (1-gamma)/3)
    # tilted_surface[2,...] = s*(gamma*muz**2 + (1-gamma)/3)
    # tilted_surface[3,...] = s*(gamma*mux*muy)
    # tilted_surface[4,...] = s*(gamma*mux*muz)
    # tilted_surface[5,...] = s*(gamma*muy*muz)

    # savemat('tilted surface gamma 0.mat', {'object': tilted_surface*100})

    ### flat  surface ###

    mu_x = np.zeros((10,61,61))
    mu_y = np.zeros((10,61,61))
    mu_z = np.zeros((10,61,61))
    s = np.zeros((10,61,61))
    mu_z[3:7,24:37,24:37] = 1
    s[4:6,25:36,25:35] = 1
    s[4:6,24,25:35] = 0.1299
    s[4:6,35,25:35] = 0.1299
    s[4:6,25:35,24] = 0.1299
    s[4:6,25:35,35] = 0.1299

    gamma = 0

    surface = np.zeros((6,10,61,61))
    surface[0,...] = s*(gamma*mu_x**2 + (1-gamma)/3)
    surface[1,...] = s*(gamma*mu_y**2 + (1-gamma)/3)
    surface[2,...] = s*(gamma*mu_z**2 + (1-gamma)/3)
    surface[3,...] = s*gamma*mu_x*mu_y
    surface[4,...] = s*gamma*mu_x*mu_z
    surface[5,...] = s*gamma*mu_y*mu_z


    savemat('flat surface gamma 0.mat', {'object': surface*100})
