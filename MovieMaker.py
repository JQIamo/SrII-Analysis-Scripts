from lyse import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as constants

################################################################################

df = data()
shot_paths = df['filepath']
images = []

fig = plt.figure()

for i in range(len(shot_paths)):
    run = Run(shot_paths[i])
    OD = run.get_result_array('VerySimpleImaging','optical_depth')
    im = plt.imshow(OD, animated=True,vmin = 0, vmax = 1)
    images.append([im])

ani = animation.ArtistAnimation(fig, images, interval=50, blit=True,
                                repeat_delay=1000)

#ani.save('MovieMaker.gif', writer = 'imagemagick', fps = 15)
