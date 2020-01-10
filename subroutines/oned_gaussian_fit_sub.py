import numpy as np
import h5py
from scipy.optimize import curve_fit

def oned_gaussian_fit_sub(optical_depth, fixcenter):

    def gaussian_1d(x,A,center,width, offset):
        return np.ravel(A*np.exp(-((x-center)/width)**2) + offset)

    xlength = np.shape(optical_depth)[1]
    zlength = np.shape(optical_depth)[0]

    with h5py.File('current_roi.h5','a') as f:
        if 'center' not in f.attrs:
            f.attrs.create('center',(xlength/2,zlength/2))
        center = f.attrs['center']

    xslice = 0 if center[1] < 0 else zlength-1 if center[1] >= zlength else center[1]
    zslice = 0 if center[0] < 0 else xlength-1 if center[0] >= xlength else center[0]

    if optical_depth[int(xslice)][int(zslice)] <0.3 and not fixcenter:
        for i in range(1000):
            guess = (np.random.randint(xlength-1),np.random.randint(zlength-1))
            if optical_depth[guess[1]][guess[0]] > 0.3:
                xslice = guess[1]
                zslice = guess[0]
                break

    initialguessz = (2,center[1],150,0)
    initialguessx = (2,center[0],150,0)

    z = np.arange(zlength)
    x = np.arange(xlength)

    imagex = optical_depth[int(xslice)]
    imagez = optical_depth[:,int(zslice)]

    try:
        for i in range(2):
            poptz, pcovz = curve_fit(gaussian_1d,z,imagez,p0 = initialguessz,
            bounds = ([0,-np.inf,0,0],[np.inf,np.inf,np.inf,np.inf]))
            xslice = 0 if poptz[1] < 0 else zlength-1 if poptz[1] >= zlength else poptz[1]
            imagex = optical_depth[int(xslice)]
            initialguessz = poptz
            poptx, pcovx = curve_fit(gaussian_1d,x,imagex,p0 = initialguessx,
            bounds = ([0,-np.inf,0,0],[np.inf,np.inf,np.inf,np.inf]))
            zslice = 0 if poptx[1] < 0 else xlength-1 if poptx[1] >= xlength else poptx[1]
            imagez = optical_depth[:,int(zslice)]
            initialguessx = poptx
    except RuntimeError:
        poptx = (0,0,0,0)
        poptz = (0,0,0,0)
        pass

    with h5py.File('current_roi.h5','a') as f:
        if poptx[0] > 0.2 and not fixcenter:
            f.attrs['center'] = (poptx[1],poptz[1])

    if (poptx[2]<30 or poptz[2]<30):
        poptx = (0,0,1,0)
        poptz = (0,0,1,0)

    return[x,imagex,poptx,z,imagez, poptz]
