from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import curve_fit
import scipy.constants as constants
from matplotlib.widgets import Cursor
from Subroutines.Integrated_gaussian_fit_sub import Integrated_gaussian_fit_sub

################################################################################

ser = data(path)
run = Run(path)

sigma0 = 1.015*10**(-7) #mm^2
mass = 87.906*1.661*10**(-27) #kg
pixelsize = 3.72*10**(-3) #mm/pix with magnification


optical_depth = run.get_result_array('VerySimpleImaging','optical_depth')
TimeOfFlight = ser['TimeOfFlight']

def gaussian_1d(x,A,center,width, offset):
    return np.ravel(A*np.exp(-((x-center)/width)**2) + offset)

x,imagex,poptx,z,imagez, poptz = Integrated_gaussian_fit_sub(optical_depth)

fig = plt.figure()
fig.suptitle(poptx[0], fontsize = 70)
ax = fig.add_subplot(221)
axx = fig.add_subplot(223)
axz = fig.add_subplot(222)
ax.imshow(optical_depth,vmin = min(poptx[3],poptx[0]), vmax = max(poptx[3],poptx[0],0.01))
axx.plot(x,imagex)
axx.plot(x,gaussian_1d(x,*poptx))
axz.plot(z,imagez)
axz.plot(z,gaussian_1d(z,*poptz))

print(poptx)
print(poptz)



tempx = mass*(poptx[2]*pixelsize*10**(-3)/TimeOfFlight)**2/(2*constants.value("Boltzmann constant"))
tempz = mass*(poptz[2]*pixelsize*10**(-3)/TimeOfFlight)**2/(2*constants.value("Boltzmann constant"))
integral = (poptx[0]+poptz[0])*0.5*np.pi*((poptx[2]+poptz[2])/2)**2*pixelsize**2/sigma0

run.save_result("peakODx",poptx[0])
run.save_result("widthx",poptx[2])
run.save_result("tempx",tempx)
run.save_result("centerx",poptx[1])
run.save_result("peakODz",poptz[0])
run.save_result("widthz",poptz[2])
run.save_result("tempz",tempz)
run.save_result("centerz",poptz[1])
run.save_result("integral",integral)
run.save_result("PeakOD",(poptx[0]+poptz[0])/2)

cursor = Cursor(ax, useblit=True, color='white', linewidth=1)
cursor.set_active(False)

def setcenter(event):
    if event.key in ['C','c'] and not cursor.active:
        cursor.set_active(True)
    elif event.key in ['C','c'] and cursor.active:
        cursor.set_active(False)
        center = (int(event.xdata),int(event.ydata))
        print(center)
        with h5py.File('G:\\My Drive\\SrII\\Labscript\\AnalysisScripts\\current_roi.h5','w') as f:
            if 'center' not in f:
                f.attrs.create('center',center)
            else:
                f.attrs['center'] = center

if not getattr(routine_storage, 'keypress_connected', False):
    plt.connect('key_press_event', setcenter)
    routine_storage.keypress_connected = True
