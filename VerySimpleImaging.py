from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py

ser = data(path)
run = Run(path)
sigma0 = 1.015*10**(-7) #mm^2
camera = 'grating'
#im_clear = run.get_image('horizontal', 'absorption', 'clear').astype('float')
im_atoms = run.get_image(camera,'absorption','atoms').astype('float')
im_probe = run.get_image(camera,'absorption','probe').astype('float')
im_background = run.get_image(camera,'absorption','background').astype('float')

im_atoms = np.maximum(im_atoms - im_background, 0.000001*np.ones(im_atoms.shape))
im_probe = np.maximum(im_probe - im_background, 0.000001*np.ones(im_probe.shape))

# im_atoms_masked = np.ma.array(im_atoms, mask = True)
# im_atoms_masked.mask[0:100, 0:100] = False
# im_probe_masked = np.ma.array(im_probe, mask = True)
# im_probe_masked.mask[0:100, 0:100] = False
# multiplier = im_atoms_masked.mean()/im_probe_masked.mean()
#im_probe = im_probe*multiplier

im_absorption = im_atoms/im_probe
optical_depth = -1 * np.log(im_absorption)

fig = plt.figure()
ax_atoms = fig.add_subplot(2,2,1)
ax_probe = fig.add_subplot(2,2,2)
ax_absorption = fig.add_subplot(2,2,3)
ax_OD = fig.add_subplot(2,2,4)

atoms = ax_atoms.imshow(im_atoms)#, vmin = 0, vmax = 50000)
probe = ax_probe.imshow(im_probe)#,vmin = 0, vmax = 50000)
absorption = ax_absorption.imshow(im_absorption, vmin = 0, vmax = 1)
OD = ax_OD.imshow(optical_depth, vmin = 0, vmax = 2)

ax_atoms.title.set_text("Atoms")
ax_probe.title.set_text("Probe")
ax_absorption.title.set_text("Absorption")
ax_OD.title.set_text("OD")

fig.colorbar(atoms, ax = ax_atoms)
fig.colorbar(probe, ax = ax_probe)
fig.colorbar(absorption, ax = ax_absorption)
fig.colorbar(OD, ax = ax_OD)

rawnum = np.sum(optical_depth)*(5.226*10**(-3))**2/sigma0

#fig.suptitle(np.average(optical_depth), fontsize = 50)

run.save_result_array('optical_depth', optical_depth)
run.save_result('rawnum',rawnum)
