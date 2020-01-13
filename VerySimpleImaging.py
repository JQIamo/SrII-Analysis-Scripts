from lyse import *
import matplotlib.pyplot as plt
import numpy as np

import AnalysisSettings
import SrConstants

ser = data(path)
run = Run(path)

camera = AnalysisSettings.Camera


# Get Images

im_atoms = run.get_image(camera, 'absorption', 'atoms').astype('float')
im_probe = run.get_image(camera, 'absorption', 'probe').astype('float')
im_background = run.get_image(camera, 'absorption', 'background').astype('float')


# Subtract Background

im_atoms = np.maximum(im_atoms - im_background, 0.000001*np.ones(im_atoms.shape))
im_probe = np.maximum(im_probe - im_background, 0.000001*np.ones(im_probe.shape))


# Calculate OD

im_absorption = im_atoms/im_probe
optical_depth = -1 * np.log(im_absorption)


# Make Figure

fig = plt.figure()
ax_atoms = fig.add_subplot(2, 2, 1)
ax_probe = fig.add_subplot(2, 2, 2)
ax_absorption = fig.add_subplot(2, 2, 3)
ax_OD = fig.add_subplot(2, 2, 4)

atoms = ax_atoms.imshow(im_atoms)
probe = ax_probe.imshow(im_probe)
absorption = ax_absorption.imshow(im_absorption, vmin=0, vmax=1)
OD = ax_OD.imshow(optical_depth, vmin=0, vmax=2)

ax_atoms.title.set_text("Atoms")
ax_probe.title.set_text("Probe")
ax_absorption.title.set_text("Absorption")
ax_OD.title.set_text("OD")

fig.colorbar(atoms, ax=ax_atoms)
fig.colorbar(probe, ax=ax_probe)
fig.colorbar(absorption, ax=ax_absorption)
fig.colorbar(OD, ax=ax_OD)


# Calculate Raw number and save

rawNum = np.sum(optical_depth)*(SrConstants.pixelSizeDict[camera])**2/SrConstants.sigma0

run.save_result_array('optical_depth', optical_depth)
run.save_result('raw_num', rawNum)
