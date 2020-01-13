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


ser = data(path)
run = Run(path)

camera = AnalysisSettings.Camera
sigma0 = SrConstants.sigma0
mass = SrConstants.mass
pixelSize = SrConstants.pixelSizeDict[camera]


# Get Data

optical_depth = run.get_result_array('VerySimpleImaging', 'optical_depth')
TimeOfFlight = ser['TimeOfFlight']


# Fit OD

x, imageX, pOptX, z, imageZ, pOptZ = splice_gaussian_fit_sub(optical_depth, False)


# Make Figure

fig = plt.figure()
fig.suptitle(pOptX[0], fontsize=70)
ax_image = fig.add_subplot(221)
ax_x = fig.add_subplot(223)
ax_z = fig.add_subplot(222)
ax_image.imshow(optical_depth, vmin=min(pOptX[3], pOptX[0]), vmax=max(pOptX[3], pOptX[0], 0.01))
ax_x.plot(x, imageX)
ax_x.plot(x, gauss(x, *pOptX))
ax_z.plot(z, imageZ)
ax_z.plot(z, gauss(z, *pOptZ))


# Calculate Things of Interest

widthX = pOptX[2] * pixelSize
widthZ = pOptZ[2] * pixelSize

centerX = pOptX[1] * pixelSize
centerZ = pOptZ[1] * pixelSize

tempX = mass * (widthX / TimeOfFlight) ** 2 / (2 * constants.value("Boltzmann constant"))
tempZ = mass * (widthZ / TimeOfFlight) ** 2 / (2 * constants.value("Boltzmann constant"))

atomNumber = (pOptX[0] + pOptZ[0]) * 0.5 * np.pi * widthX * widthZ / sigma0
avgWidth = (widthX + widthZ)/2
avgPeakOD = (pOptX[0] + pOptX[2])/2
avgTemp = (tempX + tempZ)/2


# Save Values

run.save_result("peakODX", pOptX[0])
run.save_result("widthX", widthX)
run.save_result("tempX", tempX)
run.save_result("centerX", centerX)
run.save_result("peakODZ", pOptZ[0])
run.save_result("widthZ", widthZ)
run.save_result("tempZ", tempZ)
run.save_result("centerZ", centerZ)
run.save_result("atomNumber", atomNumber)
run.save_result("avgWidth", avgWidth)
run.save_result("avgPeakOD", avgPeakOD)
run.save_result("avgTemp", avgTemp)


# Set Center event

cursor = Cursor(ax_image, useblit=True, color='white', linewidth=1)
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
