from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py

ser = data(path)
run = Run(path)


thinning = 4
#im_clear = run.get_image('horizontal', 'absorption', 'clear').astype('float')
im_atoms = run.get_image('horizontal','absorption','atoms').astype('float')
im_probe = run.get_image('horizontal','absorption','probe').astype('float')
im_background = run.get_image('horizontal','absorption','background').astype('float')

im_atoms_thinned = im_atoms[::thinning,::thinning]
im_probe_thinned = im_probe[::thinning,::thinning]
im_background_thinned = im_background[::thinning,::thinning]

im_atoms_thinned = np.maximum(im_atoms_thinned - im_background_thinned, 0.000001*np.ones(im_atoms_thinned.shape))
im_probe_thinned = np.maximum(im_probe_thinned - im_background_thinned, 0.000001*np.ones(im_probe_thinned.shape))

# im_atoms_masked = np.ma.array(im_atoms, mask = True)
# im_atoms_masked.mask[0:100, 0:100] = False
# im_probe_masked = np.ma.array(im_probe, mask = True)
# im_probe_masked.mask[0:100, 0:100] = False
# multiplier = im_atoms_masked.mean()/im_probe_masked.mean()
#im_probe = im_probe*multiplier

im_absorption = im_atoms_thinned/im_probe_thinned
optical_depth = -1 * np.log(im_absorption)

fig = plt.figure()
ax_atoms = fig.add_subplot(2,2,1)
ax_probe = fig.add_subplot(2,2,2)
ax_absorption = fig.add_subplot(2,2,3)
ax_OD = fig.add_subplot(2,2,4)

ax_atoms.imshow(im_atoms_thinned, vmin = 0, vmax = 50000)
ax_probe.imshow(im_probe_thinned,vmin = 0, vmax = 50000)
ax_absorption.imshow(im_absorption, vmin = .3, vmax = 1)
ax_OD.imshow(optical_depth)

ax_atoms.title.set_text("Atoms")
ax_probe.title.set_text("Probe")
ax_absorption.title.set_text("Absorption")
ax_OD.title.set_text("OD")

#fig.suptitle(np.average(optical_depth), fontsize = 50)

run.save_result_array('optical_depth', optical_depth)
