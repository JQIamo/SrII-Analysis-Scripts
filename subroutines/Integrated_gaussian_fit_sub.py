import numpy as np
import h5py
from scipy.optimize import curve_fit

def Integrated_gaussian_fit_sub(optical_depth):
    def gaussian_1d(x,A,center,width, offset):
        return np.ravel(A*np.exp(-((x-center)/width)**2) + offset)

    xlength = np.shape(optical_depth)[1]
    zlength = np.shape(optical_depth)[0]

    with h5py.File('current_roi.h5','a') as f:
        if 'center' not in f.attrs:
            f.attrs.create('center',(xlength/2,zlength/2))
        center = f.attrs['center']


    initialguessz = (2,center[1],150,0)
    initialguessx = (2,center[0],150,0)

    z = np.arange(zlength)
    x = np.arange(xlength)

    imagex = np.sum(optical_depth, axis = 0)
    imagez = np.sum(optical_depth, axis = 1)

    try:
        poptz, pcovz = curve_fit(gaussian_1d,z,imagez,p0 = initialguessz,
        bounds = ([0,-np.inf,0,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
        poptx, pcovx = curve_fit(gaussian_1d,x,imagex,p0 = initialguessx,
        bounds = ([0,-np.inf,0,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
    except RuntimeError:
        poptx = (0,0,0,0)
        poptz = (0,0,0,0)
        pass

    with h5py.File('current_roi.h5','a') as f:
        if poptx[0] > 0.2:
            f.attrs['center'] = (poptx[1],poptz[1])

    if (poptx[2]<30 or poptz[2]<30):
        poptx = (0,0,1,0)
        poptz = (0,0,1,0)

    return[x,imagex,poptx,z,imagez,poptz]
