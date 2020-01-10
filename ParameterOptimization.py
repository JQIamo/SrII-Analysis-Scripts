from lyse import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as constants
import AnalysisSettings
import SrConstants

################################################################################
#   Variables
################################################################################

camera = AnalysisSettings.Camera
pixelSize = SrConstants.pixelSizeDict[camera]

df = data()
################################################################################
peakODx = df["splice_gaussian_fit", "peakODx"]
peakODz = df["splice_gaussian_fit", "peakODz"]
widthX = df["splice_gaussian_fit", "widthx"]
tempX = df["splice_gaussian_fit", "tempx"]
widthZ = df["splice_gaussian_fit", "widthz"]
tempZ = df["splice_gaussian_fit", "tempz"]
width = (widthX + widthZ) / 2
temp = (tempX + tempZ) / 2
centerX = df["splice_gaussian_fit", "centerx"]
centerZ = df["splice_gaussian_fit", "centerz"]
integral = df["splice_gaussian_fit", "integral"]

################################################################################
#   Fit Functions
################################################################################


################################################################################
#   Daily Changes
################################################################################

independent_var_string = AnalysisSettings.POIndependentVar
misc_dependent_var_string = AnalysisSettings.POMiscDependentVar
num_points = 7
fit_function = AnalysisSettings.POFitFunction
p0 = [14, 100, 10, 0]

independent_var = df[independent_var_string]
misc_dependent_var = df[misc_dependent_var_string]
################################################################################

centerDiffX = np.zeros(len(centerX))
for i in range(int(len(centerX) / 2)):
    centerDiffX[2 * i] = centerX[2 * i] - centerX[2 * i + 1]
    centerDiffX[2 * i + 1] = centerDiffX[2 * i]

centerDiffZ = np.zeros(len(centerZ))
for i in range(int(len(centerZ) / 2)):
    centerDiffZ[2 * i] = centerZ[2 * i] - centerZ[2 * i + 1]
    centerDiffZ[2 * i + 1] = centerDiffZ[2 * i]

tempDiff = np.zeros(len(temp))
for i in range(int(len(temp)/2)):
    tempDiff[2 * i] = temp[2 * i] - temp[2 * i + 1]
    tempDiff[2 * i + 1] = tempDiff[2 * i]

centerDiff2d = np.sqrt(centerDiffX ** 2 + centerDiffZ ** 2)

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
axmisc.plot(independent_var, misc_dependent_var, 'bo')
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

    coeff, var_matrix = curve_fit(fit_function, independent_var, misc_dependent_var, p0=p0)
    order = np.argsort(independent_var)
    fit = fit_function(independent_var, *coeff)[order]

    axmisc.plot(independent_var[order], fit, 'r-')

    print(coeff)

    tempAlt = SrConstants.mass * (coeff[0] * pixelSize * 10 ** (-3)) ** 2 / (2 * constants.value("Boltzmann constant"))
    print(tempAlt)
