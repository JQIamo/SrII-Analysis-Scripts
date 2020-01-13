from lyse import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as constants
import AnalysisSettings
import SrConstants


camera = AnalysisSettings.Camera
pixelSize = SrConstants.pixelSizeDict[camera]

df = data()

# Get Relevant Data

avgPeakOD = df["splice_gaussian_fit", "avgPeakOD"]
avgTemp = df["splice_gaussian_fit", "avgTemp"]
atomNumber = df["splice_gaussian_fit", "atomNumber"]
centerX = df["splice_gaussian_fit", "centerX"]
centerZ = df["splice_gaussian_fit", "centerZ"]


# Set Daily Changes

independent_var_string = AnalysisSettings.POIndependentVar
misc_dependent_var_string = AnalysisSettings.POMiscDependentVar
num_points = 7
fit_function = AnalysisSettings.POFitFunction
p0 = AnalysisSettings.POFitParameters

independent_var = df[independent_var_string]
misc_dependent_var = df[misc_dependent_var_string]


# Calculate other possibly useful things

centerDiffX = np.zeros(len(centerX))
for i in range(int(len(centerX) / 2)):
    centerDiffX[2 * i] = centerX[2 * i] - centerX[2 * i + 1]
    centerDiffX[2 * i + 1] = centerDiffX[2 * i]

centerDiffZ = np.zeros(len(centerZ))
for i in range(int(len(centerZ) / 2)):
    centerDiffZ[2 * i] = centerZ[2 * i] - centerZ[2 * i + 1]
    centerDiffZ[2 * i + 1] = centerDiffZ[2 * i]

tempDiff = np.zeros(len(avgTemp))
for i in range(int(len(avgTemp)/2)):
    tempDiff[2 * i] = avgTemp[2 * i] - avgTemp[2 * i + 1]
    tempDiff[2 * i + 1] = tempDiff[2 * i]

centerDiff2d = np.sqrt(centerDiffX ** 2 + centerDiffZ ** 2)

# Make the figure

fig = plt.figure()

axOD = fig.add_subplot(221)
axTemp = fig.add_subplot(222)
axN = fig.add_subplot(223)
axMisc = fig.add_subplot(224)

axOD.plot(independent_var, avgPeakOD, 'bo')
axN.plot(independent_var, atomNumber, 'bo')
axTemp.plot(independent_var, avgTemp, 'bo')
axMisc.plot(independent_var, misc_dependent_var, 'bo')

axN.title.set_text("Number")
axOD.title.set_text("Od")
axTemp.title.set_text("Temp")
axMisc.title.set_text(misc_dependent_var_string)

axN.set_xlabel(independent_var_string)
axOD.set_xlabel(independent_var_string)
axTemp.set_xlabel(independent_var_string)
axMisc.set_xlabel(independent_var_string)


# Fit the Misc Data

if len(independent_var) > num_points:

    coeff, var_matrix = curve_fit(fit_function, independent_var, misc_dependent_var, p0=p0)
    order = np.argsort(independent_var)
    fit = fit_function(independent_var, *coeff)[order]

    axMisc.plot(independent_var[order], fit, 'r-')

    print("Coeffs: " + coeff)
