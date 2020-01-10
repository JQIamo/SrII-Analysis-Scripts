from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py

ser = data(path)
run = Run(path)
sigma0 = 1.015*10**(-7) #mm^2
camera = 'horizontal'

image = run.get_image(camera,'absorption','clear')

fig = plt.figure()
ax_image = fig.add_subplot(111)

display = ax_image.imshow(image)
fig.colorbar(display, ax = ax_image)
