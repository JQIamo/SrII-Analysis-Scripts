from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py

ser = data(path)
run = Run(path)
sigma0 = 1.015*10**(-7) #mm^2

#im_clear = run.get_image('horizontal', 'absorption', 'clear').astype('float')
im_atoms = run.get_image('grating','fluorescence','atoms').astype('float')
#im_probe = run.get_image('horizontal','absorption','probe').astype('float')
im_background = run.get_image('grating','fluorescence','background').astype('float')

im_flr = np.maximum(im_atoms - im_background, 0.000001*np.ones(im_atoms.shape))
#im_probe = np.maximum(im_probe - im_background, 0.000001*np.ones(im_probe.shape))

# im_atoms_masked = np.ma.array(im_atoms, mask = True)
# im_atoms_masked.mask[0:100, 0:100] = False
# im_probe_masked = np.ma.array(im_probe, mask = True)
# im_probe_masked.mask[0:100, 0:100] = False
# multiplier = im_atoms_masked.mean()/im_probe_masked.mean()
#im_probe = im_probe*multiplier

fig = plt.figure()
ax_atoms = fig.add_subplot(2,2,1)
ax_background = fig.add_subplot(2,2,2)
ax_fluorescence = fig.add_subplot(2,2,3)

atoms = ax_atoms.imshow(im_atoms)#, vmin = 0, vmax = 50000)
ax_background.imshow(im_background)#,vmin = 0, vmax = 50000)
ax_fluorescence.imshow(im_flr)
fig.colorbar(atoms,ax = ax_atoms)

ax_atoms.title.set_text("Atoms")
ax_background.title.set_text("Background")
ax_fluorescence.title.set_text("Fluorescence")

rawnum = np.sum(im_flr)

run.save_result('rawnum',rawnum)
