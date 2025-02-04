import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from scipy import ndimage
from scipy.io import loadmat, savemat
import plot_func as plot

def create_sphere(obj_w, obj_h, obj_d, radius, sigma, gamma, orientation=True):

    center = [np.ceil(obj_w/2), np.ceil(obj_h/2), np.ceil(obj_d/2)]
    if obj_w<2*radius+1 or obj_h<2*radius+1 or obj_d<2*radius+1:
        raise Exception ('Radius of the sphere is too large.')
    
    x, y, z = np.meshgrid(np.linspace(1,obj_w,num=obj_w), 
                          np.linspace(1,obj_h,num=obj_h), np.linspace(1,obj_d,num=obj_d), indexing='ij')
    dist = np.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
    s = np.exp(-((dist-radius)/sigma)**2)

    if orientation==False:
        return s
    
    mu_x = (x-center[0])/(dist+1e-9)
    mu_y = (y-center[1])/(dist+1e-9)
    mu_z = (z-center[2])/(dist+1e-9)

    m_xx = gamma*mu_x**2 + (1-gamma)/3
    m_yy = gamma*mu_y**2 + (1-gamma)/3
    m_zz = gamma*mu_z**2 + (1-gamma)/3
    m_xy = gamma*mu_x*mu_y
    m_xz = gamma*mu_x*mu_z
    m_yz = gamma*mu_y*mu_z

    sphere = np.stack((s*m_xx,s*m_yy,s*m_zz,s*m_xy,s*m_xz,s*m_yz), axis=0)

    return np.transpose(sphere, (0,3,1,2)), np.transpose(mu_x, (2,0,1)), np.transpose(mu_y, (2,0,1)), np.transpose(mu_z, (2,0,1))

def create_sphere_fiber(obj_w, obj_h, obj_d, radius, sigma, gamma, orientation=True):

    center = [np.ceil(obj_w/2), np.ceil(obj_h/2), np.ceil(obj_d/2)]
    if obj_w<2*radius+1 or obj_h<2*radius+1 or obj_d<2*radius+1:
        raise Exception ('Radius of the sphere is too large.')
    
    x, y, z = np.meshgrid(np.linspace(1,obj_w,num=obj_w), 
                          np.linspace(1,obj_h,num=obj_h), np.linspace(1,obj_d,num=obj_d), indexing='ij')
    dist = np.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
    s = np.exp(-((dist-radius)/sigma)**2)

    if orientation==False:
        return s
    
    mu_x = (x-center[0])/(dist+1e-9)
    mu_y = (y-center[1])/(dist+1e-9)
    mu_z = (z-center[2])/(dist+1e-9)

    theta = np.arccos(mu_z)
    phi = np.sign(mu_y)*np.arccos(mu_x/(np.sqrt(mu_x**2+mu_y**2)+1e-9))
    phi[np.isnan(phi) == True] = 0

    theta = theta-np.pi/2

    mu_x = np.sin(theta)*np.cos(phi)
    mu_y = np.sin(theta)*np.sin(phi)
    mu_z = np.cos(theta)

    m_xx = gamma*mu_x**2 + (1-gamma)/3
    m_yy = gamma*mu_y**2 + (1-gamma)/3
    m_zz = gamma*mu_z**2 + (1-gamma)/3
    m_xy = gamma*mu_x*mu_y
    m_xz = gamma*mu_x*mu_z
    m_yz = gamma*mu_y*mu_z

    sphere = np.stack((s*m_xx,s*m_yy,s*m_zz,s*m_xy,s*m_xz,s*m_yz), axis=0)

    return np.transpose(sphere, (0,3,1,2))


<<<<<<< HEAD
if __name__ == "__main__":

    # gamma = [0,0.5,1]
    # smax = 1

    # ''' Flat ring '''

    # radius = 10 

    # for g in gamma: 
    #     sphere = create_sphere(61,61,61,radius,0.5,g)
    #     flat_ring = np.zeros((6,10,61,61))
    #     flat_ring[:,5,...] = sphere[:,30,...]
=======
# this checks to see if the code is being run as a script or being imported as a module.
# if run directly, the stuff under this line executes. if not, stuff above I guess?
if __name__ == "__main__":

    gamma = [0,0.5,1]
    smax = 1 # normalized
    obj_w_l = 121;
    obj_h = 10;

    # ''' Flat ring '''

    # radius = 20 

    # for g in gamma: 
    #     sphere = create_sphere(obj_w_l,obj_w_l,obj_w_l,radius,0.5,g)
    #     flat_ring = np.zeros((6,obj_h,obj_w_l,obj_w_l))
    #     flat_ring[:,np.ceil(obj_h/2),...] = sphere[:,30,...]
>>>>>>> main
    #     for m in range (6):
    #         flat_ring[m,...] = ndimage.gaussian_filter(flat_ring[m,...], sigma=0.7)
    #     flat_ring = flat_ring*smax/np.max(np.sum(flat_ring[:3,...],axis=0))
    #     savemat('flat ring gamma ' + str(g) + '.mat', {'object': flat_ring})

    # ''' Tilted ring '''

<<<<<<< HEAD
    # radius = 10
    # center = np.array([9/2,61/2,61/2])

    # for g in gamma: 
    #     sphere = create_sphere(61,61,61,radius,0.5,g)
    #     flat_ring = np.zeros((6,10,61,61))
    #     flat_ring[:,5,...] = sphere[:,30,...]
    #     tilted_ring = np.zeros((6,10,61,61))
    #     for layer in range (8): 
    #         tilted_ring[:,layer+1,2+layer*7:2+layer*7+7,:] = flat_ring[:,5,2+layer*7:2+layer*7+7,:]
=======
    # radius = 20
    # center = np.array([9/2,obj_w_l/2,obj_w_l/2])

    # for g in gamma: 
    #     sphere = create_sphere(obj_w_l,obj_w_l,obj_w_l,radius,0.5,g)
    #     flat_ring = np.zeros((6,obj_h,obj_w_l,obj_w_l))
    #     flat_ring[:,np.ceil(obj_h/2),...] = sphere[:,30,...]
    #     tilted_ring = np.zeros((6,obj_h,obj_w_l,obj_w_l))
    #     for layer in range (8): 
    #         tilted_ring[:,layer+1,2+layer*7:2+layer*7+7,:] = flat_ring[:,np.ceil(obj_h/2),2+layer*7:2+layer*7+7,:]
>>>>>>> main
    #     non_zero = np.argwhere(np.sum(tilted_ring[:3,...],axis=0)>np.max(np.sum(tilted_ring[:3,...],axis=0))/7)
    #     for idx in non_zero:
    #         mu_x = center[1]-idx[1]
    #         mu_y = center[2]-idx[2]
    #         mu_z = center[0]-idx[0]
    #         l = np.sqrt(mu_x**2+mu_y**2+mu_z**2)
    #         mu_x = mu_x/l
    #         mu_y = mu_y/l
    #         mu_z = mu_z/l
    #         mxx = g*mu_x**2+(1-g)/3
    #         myy = g*mu_y**2+(1-g)/3
    #         mzz = g*mu_z**2+(1-g)/3
    #         mxy = g*mu_x*mu_y
    #         mxz = g*mu_x*mu_z
    #         myz = g*mu_z*mu_y
    #         tilted_ring[:,idx[0],idx[1],idx[2]] = np.array([mxx,myy,mzz,mxy,mxz,myz])

    #     for m in range (6):
    #         tilted_ring[m,...] = ndimage.gaussian_filter(tilted_ring[m,...], sigma=0.7)
    #     tilted_ring = tilted_ring*smax/np.max(np.sum(tilted_ring[:3,...],axis=0))
    #     savemat('tilted ring gamma ' + str(g) + '.mat', {'object':tilted_ring})

    # ''' Flat surface '''

<<<<<<< HEAD
    # side = 20
    # for g in gamma:
    #     flat_surface = np.zeros((6,10,61,61))
    #     flat_surface[2,5,int(30-side/2):int(30+side/2),int(30-side/2):int(30+side/2)] = 1
    #     flat_surface[2,...] = ndimage.gaussian_filter(flat_surface[2,...], sigma=0.7)
    #     flat_surface = flat_surface*smax/np.max(np.sum(flat_surface[:3,...],axis=0))
    #     savemat('flat surface gamma ' + str(g) + '.mat', {'object':flat_surface})
=======
    side = 40
    mu_x = np.zeros((obj_h,obj_w_l,obj_w_l))
    mu_y = np.zeros((obj_h,obj_w_l,obj_w_l))
    mu_z = np.zeros((obj_h,obj_w_l,obj_w_l))

    # all dipoles on the flat surface are oriented orthogonally to it (i.e.) mu = [0, 0, 1]. Therefore, set all dipoles on the flat surface to have a first moment of [0, 0, 1].
    mu_z[np.ceil(obj_h/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2)] = 1

    # there are only dipoles on the surface (i.e.) s_env = 0 when object space does not contain the flat surface.
    s_env = np.zeros((obj_h,obj_w_l,obj_w_l))
    s_env[np.ceil(obj_h/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2)] = 1
    
    for g in gamma:
        flat_surface = np.zeros((6,obj_h,obj_w_l,obj_w_l))
        flat_surface[0,:,:,:] = g*mu_x**2 + (1-g)/3
        flat_surface[1,:,:,:] = g*mu_y**2 + (1-g)/3
        flat_surface[2,:,:,:] = g*mu_z**2 + (1-g)/3
        flat_surface[3,:,:,:] = g*mu_x*mu_y
        flat_surface[4,:,:,:] = g*mu_x*mu_z
        flat_surface[5,:,:,:] = g*mu_z*mu_y
        for m in range (6):
            flat_surface[m,:,:,:] = s_env*flat_surface[m,:,:,:]
            flat_surface[m,...] = ndimage.gaussian_filter(flat_surface[m,...], sigma=0.7)
        flat_surface = flat_surface*smax/np.max(np.sum(flat_surface[:3,...],axis=0))
        savemat('flat surface gamma ' + str(g) + '.mat', {'object':flat_surface})
>>>>>>> main

    # ''' Tilted surface '''

    # side = 20
    # side_flat = np.sqrt(side**2-5**2)
    # for g in gamma:
<<<<<<< HEAD
    #     flat_surface = np.zeros((6,10,61,61))
    #     flat_surface[2,5,int(30-side/2):int(30+side/2),int(30-side/2):int(30+side/2)] = 1
    #     tilted_surface = np.zeros((6,10,61,61))
    #     for layer in range (8): 
    #         tilted_surface[:,layer+1,2+layer*7:2+layer*7+7,:] = flat_surface[:,5,2+layer*7:2+layer*7+7,:]
=======
    #     flat_surface = np.zeros((6,obj_h,obj_w_l,obj_w_l))
    #     flat_surface[2,np.ceil(obj_h/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2),int(np.floor(obj_w_l/2)-side/2):int(np.floor(obj_w_l/2)+side/2)] = 1
    #     tilted_surface = np.zeros((6,obj_h,obj_w_l,obj_w_l))
    #     for layer in range (8): 
    #         tilted_surface[:,layer+1,2+layer*7:2+layer*7+7,:] = flat_surface[:,np.ceil(obj_h/2),2+layer*7:2+layer*7+7,:]
>>>>>>> main
    #     non_zero = np.argwhere(np.sum(tilted_surface[:3,...],axis=0)>np.max(np.sum(tilted_surface[:3,...],axis=0))/7)
    #     for idx in non_zero:
    #         z = 6
    #         l = np.sqrt(side_flat**2+ + z**2)
    #         mu_x = -side_flat/l
    #         mu_z = z/l
    #         mu_y = 0
    #         mxx = g*mu_x**2+(1-g)/3
    #         myy = g*mu_y**2+(1-g)/3
    #         mzz = g*mu_z**2+(1-g)/3
    #         mxy = g*mu_x*mu_y
    #         mxz = g*mu_x*mu_z
    #         myz = g*mu_z*mu_y
    #         tilted_surface[:,idx[0],idx[1],idx[2]] = np.array([mxx,myy,mzz,mxy,mxz,myz])
    #     for m in range (6):
    #         tilted_surface[m,...] = ndimage.gaussian_filter(tilted_surface[m,...], sigma=0.7)
    #     tilted_surface = tilted_surface*smax/np.max(np.sum(tilted_surface[:3,...],axis=0))
    #     savemat('tilted surface gamma ' + str(g) + '.mat', {'object':tilted_surface})

    # ''' Hemisphere '''

    # radius = 14
    # for g in gamma: 
<<<<<<< HEAD
    #     sphere = create_sphere(61,61,61,radius,0.5,g)
    #     hemisphere = np.zeros((6,10,61,61))
=======
    #     sphere = create_sphere(obj_w_l,obj_w_l,obj_w_l,radius,0.5,g)
    #     hemisphere = np.zeros((6,obj_h,obj_w_l,obj_w_l))
>>>>>>> main
    #     hemisphere[:,2:8,...] = sphere[:,39:45,...]
    #     for m in range (6):
    #         hemisphere[m,...] = ndimage.gaussian_filter(hemisphere[m,...], sigma=0.7)
    #     hemisphere = hemisphere*smax/np.max(np.sum(hemisphere[:3,...],axis=0))
    #     savemat('hemisphere gamma ' + str(g) + '.mat',{'object':hemisphere})

    # ''' Fiber '''

    # radius = 14
    # for g in gamma: 
<<<<<<< HEAD
    #     sphere = create_sphere_fiber(61,61,61,radius,0.5,g)
    #     fiber = np.zeros((6,10,61,61))
=======
    #     sphere = create_sphere_fiber(obj_w_l,obj_w_l,obj_w_l,radius,0.5,g)
    #     fiber = np.zeros((6,obj_h,obj_w_l,obj_w_l))
>>>>>>> main
    #     fiber[:,2:8,30,:] = sphere[:,39:45,30,:]
    #     fiber[:,2:8,29,:] = sphere[:,39:45,30,:]
    #     fiber[:,2:8,31,:] = sphere[:,39:45,30,:]
    #     fiber[:,2:8,:,30] = sphere[:,39:45,:,30]
    #     fiber[:,2:8,:,29] = sphere[:,39:45,:,30]
    #     fiber[:,2:8,:,31] = sphere[:,39:45,:,30]
    #     for m in range (6):
    #         fiber[m,...] = ndimage.gaussian_filter(fiber[m,...], sigma=0.7)
    #     fiber = fiber*smax/np.max(np.sum(fiber[:3,...],axis=0))
    #     savemat('fiber gamma ' + str(g) + '.mat',{'object':fiber})
    
    ''' Lipid membrane simulation '''

<<<<<<< HEAD
    radius = 2000/2/66
    sphere, mu_x, mu_y, mu_z = create_sphere(101,101,101,radius,2,0)
    sphere = sphere[:,35:58,:,:]
    mu_x = mu_x[35:58,...]
    mu_y = mu_y[35:58,...]
    mu_z = mu_z[35:58,...]
    membrane = np.zeros((6,23,101,101))
    for i in range (13):
        gamma = i/12
        membrane[0,i,...] = gamma*mu_x[i,...]**2 + (1-gamma)/3
        membrane[1,i,...] = gamma*mu_y[i,...]**2 + (1-gamma)/3
        membrane[2,i,...] = gamma*mu_z[i,...]**2 + (1-gamma)/3
        membrane[3,i,...] = gamma*mu_x[i,...]*mu_y[i,...]
        membrane[4,i,...] = gamma*mu_x[i,...]*mu_z[i,...]
        membrane[5,i,...] = gamma*mu_y[i,...]*mu_z[i,...]
    sphere[:,13:,...] = 0
    # plot.video_obj(sphere,'simulated lipid membrane','object GT')
    savemat('simulated lipid membrane.mat',{'object':sphere*3000})
=======
    # radius = 2000/2/66
    # sphere, mu_x, mu_y, mu_z = create_sphere(obj_w_l,obj_w_l,obj_w_l,radius,2,0)
    # sphere = sphere[:,35:58,:,:]
    # mu_x = mu_x[35:58,...]
    # mu_y = mu_y[35:58,...]
    # mu_z = mu_z[35:58,...]
    # membrane = np.zeros((6,23,obj_w_l,obj_w_l))
    # for i in range (13):
    #     gamma = i/12
    #     membrane[0,i,...] = gamma*mu_x[i,...]**2 + (1-gamma)/3
    #     membrane[1,i,...] = gamma*mu_y[i,...]**2 + (1-gamma)/3
    #     membrane[2,i,...] = gamma*mu_z[i,...]**2 + (1-gamma)/3
    #     membrane[3,i,...] = gamma*mu_x[i,...]*mu_y[i,...]
    #     membrane[4,i,...] = gamma*mu_x[i,...]*mu_z[i,...]
    #     membrane[5,i,...] = gamma*mu_y[i,...]*mu_z[i,...]
    # sphere[:,13:,...] = 0
    # # plot.video_obj(sphere,'simulated lipid membrane','object GT')
    # savemat('simulated lipid membrane.mat',{'object':sphere*3000})

    # ''' Line objects for convolution direction test '''
    # line = np.zeros((6,1,obj_w_l,obj_w_l))
    # line[1,:,40:60,30] = 1
    # plot.video_obj(line,'line object for convolution direction test', 'object GT')
    # savemat('line horizontal rec.mat',{'object':line})
>>>>>>> main
