from lyse import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as constants

################################################################################
#   Variables
################################################################################

sigma0 = 1.015*10**(-7) #mm^2
mass = 87.906*1.661*10**(-27) #kg
pixelsize = 3.72*10**(-3) #mm/pix with magnification


df = data()
################################################################################
peakODx = df["1d_gaussian_fit","peakODx"]
peakODz = df["1d_gaussian_fit","peakODz"]
widthx = df["1d_gaussian_fit","widthx"]
tempx = df["1d_gaussian_fit","tempx"]
widthz = df["1d_gaussian_fit","widthz"]
tempz = df["1d_gaussian_fit","tempz"]
width = (widthx + widthz)/2
temp = (tempx + tempz)/2
centerx = df["1d_gaussian_fit","centerx"]
centerz = df["1d_gaussian_fit","centerz"]
integral = df["1d_gaussian_fit", "integral"]

################################################################################
#   Fit Functions
################################################################################
def linear_y_offset(x,*p):
    m, y0 = p
    return m*x + y0

def linear_x_offset(x,*p):
    m, x0 = p
    return m*(x - x0)

def gauss(x, *p):
    A, center, width, offset = p
    return A*np.exp(-(x-center)**2/(2.*width**2)) + offset

def gauss_p_linear(x, *p):
    A, center, width, offset, drift = p
    return A*np.exp(-(x-center)**2/(2.*width**2)) + offset + drift*(x-center)

def double_gauss(x, *p):
    A, center, split, width, offset, drift = p
    return A*(np.exp(-(x-center+split/2)**2/(2.*width**2)) + np.exp(-(x-center-split/2)**2/(2.*width**2))) + offset + drift*(x-center)

def parabola(x, *p):
    a, b, c = p
    return a*((x-b)**2)+c

def exp_decay(x, *p):
    a, t, y0 = p
    return a*np.exp(-x*t)+y0

def temp_fit(x, *p):
    w0, Temp = p
    vrms = np.sqrt(2*constants.value("Boltzmann constant")*Temp/mass)
    return sqrt(w0**2 + (vrms*x)**2)

################################################################################
#   Daily Changes
################################################################################

independent_var_string = "run number"
misc_dependent_var_string = ("Integrated_gaussian_fit","peakODz")
num_points = 7
fit_function = gauss
p0 = [14,100,10,0]



independent_var = df[independent_var_string]
misc_dependent_var = df[misc_dependent_var_string]
################################################################################

centerdiffx = np.zeros(len(centerx))
for i in range(int(len(centerx)/2)):
    centerdiffx[2*i] = centerx[2*i]-centerx[2*i+1]
    centerdiffx[2*i+1] = centerdiffx[2*i]

centerdiffz = np.zeros(len(centerz))
for i in range(int(len(centerz)/2)):
    centerdiffz[2*i] = centerz[2*i]-centerz[2*i+1]
    centerdiffz[2*i+1] = centerdiffz[2*i]

tempdiff = np.zeros(len(temp))
for i in range(int(len(temp)/2)):
    tempdiff[2*i] = temp[2*i]-temp[2*i+1]
    tempdiff[2*i+1] = tempdiff[2*i]

centerdiff2d = np.sqrt(centerdiffx**2 + centerdiffz**2)

################################################################################
#   Making Figures
################################################################################
fig = plt.figure()
axod = fig.add_subplot(221)
axwidth = fig.add_subplot(222)
axn = fig.add_subplot(223)
axmisc = fig.add_subplot(224)
axod.plot(independent_var, (peakODx + peakODz)/2, 'bo')
axn.plot(independent_var, integral, 'bo')
axwidth.plot(independent_var, temp, 'bo')
axmisc.plot(independent_var, misc_dependent_var,'bo')
axn.title.set_text("number")
axod.title.set_text("od")
axwidth.title.set_text("temp")
axmisc.title.set_text("misc")
axn.set_xlabel(independent_var_string)
axod.set_xlabel(independent_var_string)
axwidth.set_xlabel(independent_var_string)
axmisc.set_xlabel(independent_var_string)

################################################################################
#   Fitting
################################################################################

if len(independent_var) > num_points:

    coeff, var_matrix = curve_fit(fit_function, independent_var, misc_dependent_var, p0 = p0)
    order = np.argsort(independent_var)
    fit = fit_function(independent_var, *coeff)[order]

    axmisc.plot(independent_var[order], fit, 'r-')

    print(coeff)

    tempalt = mass*(coeff[0]*pixelsize*10**(-3))**2/(2*constants.value("Boltzmann constant"))
    print(tempalt)
