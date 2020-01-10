from lyse import *
import matplotlib.pyplot as plt
import numpy as np

ser = data(path)
run = Run(path)


im_one = run.get_image('horizontal','absorption','clear').astype('float')

plt.imshow(im_one)
print(np.average(im_one))
run.save_result('clear_average',np.average(im_one))
