from lyse import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import curve_fit

ser = data(path)
run = Run(path)

optical_depth = run.get_result_array('VerySimpleImaging','optical_density')

def gaussian_2d(index,A,x0,y0,wx, wy, offset):
    x,y = index
    return np.ravel(A*np.exp(-((x-x0)/wx)**2-((y-y0)/wy)**2) + offset)

initialguess = (2,600,1250,500,500,1000)

x = np.arange(1448)
print(x)
y = np.arange(1928)
x,y = np.meshgrid(x,y)
index = np.vstack((x.ravel(),y.ravel()))

popt, pcov = curve_fit(gaussian_2d,index,optical_depth.ravel(),p0 = initialguess)

fit = np.zeros(optical_depth.shape)

for i in np.arange(1448):
    for j in np.arange(1928):
        fit[i][j] = gaussian_2d((i,j),*popt)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(optical_depth)
ax.contour(fit)

print(popt)
