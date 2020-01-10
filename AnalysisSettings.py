"""This is a place for all the things we change on a day to day basis"""
import Subroutines.FitFunctions as Functions

# Multi-script variables
Camera = 'grating'
Isotope = 88

# Single-Script variables

# ParameterOptimization:
POIndependentVar = 'run number'
POMiscDependentVar = ("Integrated_gaussian_fit", "peakODz")
POFitFunction = Functions.gaussian_1d
POFitParameters = [2, 100, 10, 0]
