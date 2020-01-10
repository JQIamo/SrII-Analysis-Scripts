from lyse import *
import matplotlib.pyplot as plt
import numpy as np


df = data()
shot_paths = df['filepath']

fig = plt.figure()

run = Run(shot_paths[0])
image = run.get_result_array('VerySimpleImaging','optical_depth')


average_density = np.zeros(np.shape(image))

for i in range(len(shot_paths)):
    run = Run(shot_paths[i])
    image = run.get_result_array('VerySimpleImaging','optical_depth')
    average_density += image

plt.imshow(average_density,vmin = 0, vmax = 1)
