from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy.constants as constants
from matplotlib.widgets import Cursor
from Subroutines.splice_gaussian_fit_sub import splice_gaussian_fit_sub
import SrConstants
import AnalysisSettings
from Subroutines.FitFunctions import gauss

################################################################################

ser = data(path)
run = Run(path)

camera = AnalysisSettings.Camera
sigma0 = 1.015*10**(-7)  # mm^2
mass = 87.906*1.661*10**(-27)  # kg
pixelSize = SrConstants.pixelSizeDict[camera]

optical_depth = run.get_result_array('VerySimpleImaging', 'optical_depth')
TimeOfFlight = ser['TimeOfFlight']

x, imageX, pOptX, z, imageZ, pOptZ = splice_gaussian_fit_sub(optical_depth, False)


fig = plt.figure()
fig.suptitle(pOptX[0], fontsize=70)
ax = fig.add_subplot(221)
axx = fig.add_subplot(223)
axz = fig.add_subplot(222)
ax.imshow(optical_depth, vmin=min(pOptX[3], pOptX[0]), vmax=max(pOptX[3], pOptX[0], 0.01))
axx.plot(x, imageX)
axx.plot(x, gauss(x, *pOptX))
axz.plot(z, imageZ)
axz.plot(z, gauss(z, *pOptZ))

print(pOptX)
print(pOptZ)

tempX = mass * (pOptX[2] * pixelSize * 10 ** (-3) / TimeOfFlight) ** 2 / (2 * constants.value("Boltzmann constant"))
tempZ = mass * (pOptZ[2] * pixelSize * 10 ** (-3) / TimeOfFlight) ** 2 / (2 * constants.value("Boltzmann constant"))
integral = (pOptX[0] + pOptZ[0]) * 0.5 * np.pi * pOptX[2] * pOptZ[2] * pixelSize ** 2 / sigma0
run.save_result("peakODx", pOptX[0])
run.save_result("widthx", pOptX[2])
run.save_result("tempx", tempX)
run.save_result("centerx", pOptX[1])
run.save_result("peakODz", pOptZ[0])
run.save_result("widthz", pOptZ[2])
run.save_result("tempz", tempZ)
run.save_result("centerz", pOptZ[1])
run.save_result("integral", integral)
run.save_result("PeakOD", (pOptX[0] + pOptZ[0]) / 2)

cursor = Cursor(ax, useblit=True, color='white', linewidth=1)
cursor.set_active(False)


def set_center(event):
    if event.key in ['C', 'c'] and not cursor.active:
        cursor.set_active(True)
    elif event.key in ['C', 'c'] and cursor.active:
        cursor.set_active(False)
        center = (int(event.xdata), int(event.ydata))
        print(center)
        with h5py.File('G:\\My Drive\\SrII\\Labscript\\AnalysisScripts\\current_roi.h5', 'w') as f:
            if 'center' not in f:
                f.attrs.create('center', center)
            else:
                f.attrs['center'] = center


if not getattr(routine_storage, 'keypress_connected', False):
    plt.connect('key_press_event', set_center)
    routine_storage.keypress_connected = True
