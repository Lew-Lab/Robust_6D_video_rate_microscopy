import glob
import os
import subprocess
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib as mpl

def plot_obj_voxel(obj, slice_dim, slice_loc, filename):

    titles = ['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz']
    fig, axs = plt.subplots(2,3)

    if slice_dim == 'x':
        slice = obj[:,np.int16(obj.shape[1]*slice_loc+1),:,:]
        fig.suptitle('Slice at x=%.2fwidth of the object' % slice_loc)
    elif slice_dim == 'y':
        slice = obj[:,:,np.int16(obj.shape[2]*slice_loc+1),:]
        fig.suptitle('Slice at y=%.2fheight of the object' % slice_loc)
    elif slice_dim == 'z':
        slice = obj[:,:,:,np.int16(obj.shape[3]*slice_loc+1)]
        fig.suptitle('Slice at z=%.2fdepth of the object' % slice_loc)
    else: 
        raise Exception ('Dimension of the section slice needs to be specified.')

    for i in range (6):
        ax = axs[np.int16(i/3)][i-np.int16(i/3)*3]
        temp = ax.imshow(slice[i,...])
        ax.set_axis_off()
        ax.set_title(titles[i])
        fig.colorbar(temp, ax=ax, location='right')

    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    fig.savefig(filename)
    # plt.show()

    return

def plot_img(img, title_name, filename):
    
    channel = img.shape[0]
    
    if channel == 2:
        fig, axs = plt.subplots(1,2)
        temp = axs[0].imshow(img[0,:,:])
        axs[0].set_axis_off()
        fig.colorbar(temp, ax=axs[0], location='right')
        temp = axs[1].imshow(img[1,:,:])
        axs[1].set_axis_off()
        fig.colorbar(temp, ax=axs[1], location='right')
        plt.suptitle(title_name)
        plt.show()

    if channel == 8:
        fig, axs = plt.subplots(2,4)
        temp = axs[0][0].imshow(img[0,:,:])
        axs[0][0].set_axis_off()
        fig.colorbar(temp, ax=axs[0][0], location='right')
        temp = axs[0][1].imshow(img[1,:,:])
        axs[0][1].set_axis_off()
        fig.colorbar(temp, ax=axs[0][1], location='right')
        temp = axs[0][2].imshow(img[2,:,:])
        axs[0][2].set_axis_off()
        fig.colorbar(temp, ax=axs[0][2], location='right')
        temp = axs[0][3].imshow(img[3,:,:])
        axs[0][3].set_axis_off()
        fig.colorbar(temp, ax=axs[0][3], location='right')
        temp = axs[1][0].imshow(img[4,:,:])
        axs[1][0].set_axis_off()
        fig.colorbar(temp, ax=axs[1][0], location='right')
        temp = axs[1][1].imshow(img[5,:,:])
        axs[1][1].set_axis_off()
        fig.colorbar(temp, ax=axs[1][1], location='right')
        temp = axs[1][2].imshow(img[6,:,:])
        axs[1][2].set_axis_off()
        fig.colorbar(temp, ax=axs[1][2], location='right')
        temp = axs[1][3].imshow(img[7,:,:])
        axs[1][3].set_axis_off()
        fig.colorbar(temp, ax=axs[1][3], location='right')
        plt.suptitle(title_name)
        # plt.show()

    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    fig.savefig(filename)

    return

def plot_obj_dipole(obj, slice_dim, slice_loc):

    w = obj.shape[1]
    h = obj.shape[2]
    d = obj.shape[3]

    s = obj[0,...]+obj[1,...]+obj[2,...]

    mxx = obj[0,:,:,:]/s
    myy = obj[1,:,:,:]/s
    mzz = obj[2,:,:,:]/s
    mxy = obj[3,:,:,:]/s
    mxz = obj[4,:,:,:]/s
    myz = obj[5,:,:,:]/s

    cos_theta = np.sqrt(mzz)[50,:,:]
    sin_theta = np.sqrt(1-mzz)[50,:,:]
    cos_phi = mxz[50,...]/(cos_theta*sin_theta)
    sin_phi = myz[50,...]/(cos_theta*sin_theta)

    y_center = 51 + np.linspace(0,obj.shape[1],num=obj.shape[1])*100
    z_center = 51 + np.linspace(0,obj.shape[2],num=obj.shape[2])*100

    L = 40

    lines = []
    color = []

    for i in range(100):
        for j in range(100):
            start_pt = (y_center[i]-L/2*sin_theta[i,j]*sin_phi[i,j], z_center[j]-L/2*cos_theta[i,j])
            end_pt = (y_center[i]+L/2*sin_theta[i,j]*sin_phi[i,j], z_center[j]+L/2*cos_theta[i,j])
            lines.append([start_pt, end_pt])
            color.append((1,0,0))

    lc = mc.LineCollection(lines,linewidth=.5,color=color)

    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.margins(0.1)
    plt.show()

    return

def zstack_video(img_zstack, folder):
    zstack_size = len(img_zstack[0,:,0,0])
    print(zstack_size)
    for i in range(zstack_size):
        img_zslice = img_zstack[:,i,:,:]
        plot_img(img_zslice, 'Reconstructed z-slice #' +str(i+1), os.path.join(folder, 'Reconstructed z-slice #' +str(i+1)))
        # plt.savefig(folder + "/%02d.png" % i)

    #subprocess.call([
        # 'ffmpeg', '-framerate', '8', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
        # 'video_name.mp4'
    # ])
    # for file_name in glob.glob("*.png"):
        # os.remove(file_name)