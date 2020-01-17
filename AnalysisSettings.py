"""This is a place for all the things we change on a day to day basis"""
import Subroutines.FitFunctions as Functions

# Multi-script variables
Camera = 'grating'
Isotope = 88

# Single-Script variables

# ParameterOptimization:
POIndependentVar = 'TimeOfFlight'
POMiscDependentVar1 = ("splice_gaussian_fit", "widthX")
POMiscDependentVar2 = ("splice_gaussian_fit", "widthZ")
POFitFunction = Functions.temp_fit
POFitParameters = [50,0.004]

#FIIMOTLoad settings:
FMLMaxCull = 30000000
FMLMinCull = 2000000
