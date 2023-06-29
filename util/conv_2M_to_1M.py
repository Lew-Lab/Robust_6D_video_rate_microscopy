import numpy as np

def conv_2M_to_1M(object,psf):

    mxx = object(0,...)
    myy = object(1,...)
    mzz = object(2,...)
    mxy = object(3,...)
    mxz = object(4,...)
    myz = object(5,...)

    s = np.sum(object, axis=1)

    m2 = np.concatenate((np.reshape(mxx,-1),np.reshape(myy,-1),np.reshape(mzz,-1),np.reshape(mxy,-1),np.reshape(mxz,-1),np.reshape(myz,-1)),axis=0)

    idx_nonzero_3d = np.where(s>10e-2)
    idx_nonzero_1d = np.where(np.reshape(s,axis=-1)>10e-2)
    m2 = m2[:,idx_nonzero_1d]

    return 

def secondM2SymmConeWeighted(b,B,sumNorm,secM,signal,backg):
    Ix = secM[0]*b.XX+secM[1]*b.YY+secM[2]*b.ZZ+secM[3]*b.XY+secM[4]*b.XZ+secM[5]*b.YZ
    Ix[np.where(np.less(Ix,0))] = 0
    Ix = Ix/sumNorm
    Ix = signal*Ix