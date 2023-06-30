import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib as mpl
import subprocess

def plot_obj_voxel(obj, slice_dim, slice_loc, title=''): 

    ''' Plot a sectional view of the object. '''

    titles = ['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz']
    fig, axs = plt.subplots(2,3)
    fig.suptitle(title)

    if slice_dim == 'x':
        slice = obj[:,:,np.int16(obj.shape[2]*slice_loc),:]
        fig.suptitle(title+': Slice at x=%.1f' % slice_loc)
    elif slice_dim == 'y':
        slice = obj[:,:,:,np.int16(obj.shape[3]*slice_loc)]
        fig.suptitle(title+': Slice at y=%.1f' % slice_loc)
    elif slice_dim == 'z':
        slice = obj[:,np.int16(obj.shape[1]*slice_loc),:,:]
        fig.suptitle(title+': Slice at z=%.1f' % slice_loc)
    else: 
        raise Exception ('Dimension of the section slice needs to be specified.')

    for i in range (6):
        ax = axs[np.int16(i/3)][i-np.int16(i/3)*3]
        temp = ax.imshow(np.transpose(slice[i,...]), origin='lower')
        ax.set_axis_off()
        ax.set_title(titles[i])
        fig.colorbar(temp, ax=ax, location='right')

    plt.show()

    return

def video_obj(obj, title_name, folder_name): 

    ''' Generate an stack of z layers of the 6 second moments. '''

    if obj.shape[0] == 6:
        titles = ['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz']
        stacks = obj.shape[1]
        vmin = np.zeros((6,1))
        vmax = np.zeros((6,1))

        for i in range (6):
            vmin[i] = np.min(obj[i,...])
            vmax[i] = np.max(obj[i,...])

        for z in range(stacks):
            slice = obj[:,z,...]

            fig, axs = plt.subplots(2,3)
            fig.suptitle(title_name + ' at z #%2d' %z)

            for i in range (6):
                ax = axs[np.int16(i/3)][i-np.int16(i/3)*3]
                temp = ax.imshow(np.transpose(slice[i,...]), vmin=vmin[i], vmax=vmax[i], origin='lower')
                ax.set_axis_off()
                ax.set_title(titles[i])
                fig.colorbar(temp, ax=ax, location='right')

            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            plt.savefig(folder_name + "/file%02d.png" % z)
            plt.close(fig)

    if obj.shape[0] == 1:
        stacks = obj.shape[1]
        vmin = np.min(obj)
        vmax = np.max(obj)

        for z in range(stacks):
            slice = obj[0,z,...]

            fig, axs = plt.subplots(1,1)
            fig.suptitle(title_name + ' at z #%2d' %z)

            temp = axs.imshow(np.transpose(slice), vmin=vmin, vmax=vmax, origin='lower')
            axs.set_axis_off()
            fig.colorbar(temp, ax=axs, location='right')

            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            plt.savefig(folder_name + "/file%02d.png" % z)
            plt.close(fig)

    return


def plot_img_tight(img, title_name):

    ''' Plot the image in a tight style. '''

    channel = img.shape[0]

    if channel == 8:
        image_red = np.concatenate((np.transpose(img[0,...]),np.transpose(img[1,...]),np.transpose(img[2,...]),np.transpose(img[3,...])),axis=1)
        image_blue = np.concatenate((np.transpose(img[4,...]),np.transpose(img[5,...]),np.transpose(img[6,...]),np.transpose(img[7,...])),axis=1)
        fig, axs = plt.subplots(2,1)
        temp = axs[0].imshow(image_red, cmap='gist_heat', origin='lower')
        axs[0].set_axis_off()
        fig.colorbar(temp, ax=axs[0])
        temp = axs[1].imshow(image_blue, cmap='bone', origin='lower')
        axs[1].set_axis_off()
        fig.colorbar(temp, ax=axs[1])
        plt.suptitle(title_name)
        plt.show()

    return

def plot_img(img, title_name):

    ''' Plot the image. '''
    
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
        plt.show()

    return

def plot_obj_dipole(obj, slice_dim, slice_loc):

    ''' Still working. 
    3D dipole plot. '''

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

def plot_loss(loss):

    fig, axs = plt.subplots(3,2)

    axs[0][0].remove()
    axs[0][1].remove()
    axbig = fig.add_subplot(axs[0])
    axbig.plot(loss['total'])
    axbig.set_title('total')

    ax = axs[1][0]
    ax.plot(loss['lsq'])
    ax.set_title('Least square')
    

    ax = axs[1][1]
    ax.plot(loss['L1'])
    ax.set_title('L1')

    ax = axs[2][0]
    ax.plot(loss['TV'])
    ax.set_title('TV')

    ax = axs[2][1]
    ax.plot(loss['constraint_positive_mii'])
    ax.set_title('Positive constraint')
