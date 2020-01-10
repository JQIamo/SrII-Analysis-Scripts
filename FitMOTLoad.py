"""
@author: spe
@author: ebn
@author: dsbarker
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import
from lyse import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import h5py
#import lmfit
import scipy.constants as cts
#from analysislib.Subroutines import ImageFitter
#from analysislib.Subroutines import SeqFitter
import imageStackViewer
#from analysislib.Subroutines.alkaliAtom import alkaliAtom
from labscript_utils.connections import _ensure_str
from Subroutines.splice_gaussian_fit_sub import splice_gaussian_fit_sub
from Subroutines.Integrated_gaussian_fit_sub import Integrated_gaussian_fit_sub

# Initial setup for single shot routine:
ser = data(path)
run = Run(path)

# MOT specifer:
mot_specifier = 'Sr_Grating'
show_loading = 2
show_background = False

# Import atomic parameters:
# atom = alkaliAtom(ser['isotope'])

# Camera of interst:
# image_labels = run.get_all_image_labels()
# orientations = image_labels.keys()
# if len(orientations)>1:
#     #TODO: FIX
#     raise ValueError('More than one orientation detected, not ready to handle that.')
# else:
#     for orientation_of_interest in orientations:
#         cams_of_interest = image_labels[orientation_of_interest]
#
# if len(cams_of_interest)>1:
#     #TODO: FIX
#     raise ValueError('More than one camera detected, not ready to handle that.')
# else:
#     for cam_of_interest in cams_of_interest:
#         break
orientation_of_interest = 'grating'
cam_of_interest = 'FleaCamera_gMOT'

#Load the time axis directly from the camera exposure times:
with h5py.File(path) as f:
    # Load in the MOT video:
    im_loading_raw = []
    for ii in range(int(ser['BlueMOTLoadTime']/ser['TimeBetweenExp'])):
        im_loading_raw.append(run.get_image(orientation_of_interest, 'fluorescence', 'atoms'+str(ii)))
    im_loading_raw = np.array(im_loading_raw)
    print(im_loading_raw.shape)

    #Peter's Addition
    ############################################################################
    imagefitsx = []
    imagefitsz = []
    atomNumbers = np.array([])
    time = np.array([])
    atomNumbersCulled = np.array([])
    timeCulled = np.array([])

    for image in im_loading_raw:
        imagefitsx.append(splice_gaussian_fit_sub(image, True)[2])
        imagefitsz.append(splice_gaussian_fit_sub(image, True)[3])

    i = 0
    for i in range(len(imagefitsx)):
        atomNumber = (imagefitsx[i][0]+imagefitsz[i][0])*0.5*np.pi*((imagefitsx[i][2]+imagefitsz[i][2])/2)**2#*pixelsize**2/sigma0
        if atomNumber > 20000:
            timeCulled = np.append(timeCulled,i*ser['TimeBetweenExp'])
            atomNumbersCulled = np.append(atomNumbersCulled,atomNumber)
        time = np.append(time,i*ser['TimeBetweenExp'])
        atomNumbers = np.append(atomNumbers,atomNumber)

    def load_fit_func(t,*p):
        rate, tau = p
        return rate*tau*(1-np.exp(-t/tau))

    coeff, var_matrix = curve_fit(load_fit_func, timeCulled, atomNumbersCulled, [50000,.1])
    atomNumbersNorm = atomNumbers*ser["1d_gaussian_fit","integral"]/(coeff[0]*coeff[1])
    coeff[0] = ser["1d_gaussian_fit","integral"]/coeff[1]
    fit = load_fit_func(time, *coeff)

    atomNfig = plt.figure()
    atomNax = atomNfig.add_subplot(111)
    atomNax.plot(time,atomNumbersNorm,'bo')
    atomNax.plot(time, fit, 'r-')

    ############################################################################

    # Get the properties of all the movies:
    frame_times = f['devices'][cam_of_interest]['EXPOSURES']['t']
    exposure_times = f['devices'][cam_of_interest]['EXPOSURES']['trigger_duration']

    # Now, specifically select the properties for this movie
    ind = f['devices'][cam_of_interest]['EXPOSURES']['name'] == 'fluorescence'

    t_btw_exposures = frame_times[ind][2] - frame_times[ind][1]
    t_start = frame_times[ind][1]
    t_exposure = exposure_times[ind][0]
    t_len = t_btw_exposures*len(ind)

    # Now, make a time axies
    t = np.arange(0., t_len, t_btw_exposures)

    # Now pull in the background:
    im_background = run.get_image(orientation_of_interest, 'fluorescence', 'background')
        #im_background = np.mean(im_background, axis=0)


    # # Pull in ROI:
    # fit_roi = run.get_result_array('defineROIs', 'fit_ROI')
    # fit_center = run.get_result_array('defineROIs', 'fit_center')
    #
    # # Pull in camera info:
    # cam_attrs = f['devices'][cam_of_interest].attrs
    # p_size = cam_attrs['pixel_size']
    # mag = cam_attrs['magnification']
    # NA = cam_attrs['NA']
    # quantum_efficiency = cam_attrs['quantum_efficiency']
    # transmission = cam_attrs['transmission']
    # counts_per_photoelectron = cam_attrs['counts_per_photoelectron']

# Make the background images (subtracting background)
im_loading = np.zeros(im_loading_raw.shape, dtype='int16')
for jj in range(im_loading.shape[0]):
    im_loading[jj] = im_loading_raw[jj]-im_background#*t_exposure/t_exposure_bg

# Make an image stacker of the load:
if show_loading==1:
    plt.figure()
    plt.clf()
    fig1, ax1 = plt.subplots(1, 1)
    tracker = imageStackViewer.IndexTracker(
        fig1, ax1, np.concatenate((np.array([im_background]), im_loading_raw))
        )
elif show_loading==2:
    plt.figure()
    plt.clf()
    fig1, ax1 = plt.subplots(1, 1)
    tracker = imageStackViewer.IndexTracker(fig1, ax1, im_loading, scale=False)
elif show_loading==3:
    plt.figure()
    plt.clf()
    fig1, ax1 = plt.subplots(1, 1)
    tracker = imageStackViewer.IndexTracker(
        fig1,
        ax1,
        im_loading_raw[:,  fit_roi[1]:fit_roi[3], fit_roi[0]:fit_roi[2]],
        scale=True
    )
elif show_loading==4:
    plt.figure()
    plt.clf()
    fig1, ax1 = plt.subplots(1, 1)
    tracker = imageStackViewer.IndexTracker(
        fig1,
        ax1,
        im_loading[:,  fit_roi[1]:fit_roi[3], fit_roi[0]:fit_roi[2]],
        scale=False
    )


# if show_background:
#     plt.figure()
#     plt.imshow(im_background)
#
#
#
# # Run analysis on each MOT image during the load:
# # TODO: base the guess on reasonable things, this is after all a flr image!
# motim_load = []
# for jj in range(im_loading.shape[0]):
#     motim_load.append(
#         ImageFitter.ImageFitter(
#             #im_loading[jj, fit_roi[1]:fit_roi[3], fit_roi[0]:fit_roi[2]],
#             im_loading[jj, fit_roi[1]:fit_roi[3], fit_roi[0]:fit_roi[2]],
#             [fit_center[0]-fit_roi[0], fit_center[1]-fit_roi[1]],
#             'Gauss2D',
#             im_type='fluorescence'
#             )
#         )
#
# #motim_load[-1].showFit(fignum='Fit')
#
# # Other important parameters during the load:
# #I_load = 2*Li_3D_MOT_Fractional_Carrier_Power*\
# #         Li_3D_MOT_Power_Load/(np.pi*(laser_waist/10)**2) # mW/cm^2
# if ser['isotope']<30:
#     carrier_power = 0.66
# else:
#     carrier_power = 0.98
# intensity = 2*carrier_power*\
#             (ser[mot_specifier + '_Power_Load']*ser[mot_specifier + '_dP_dV']+
#              ser[mot_specifier + '_P0'])/\
#             (np.pi*(ser[mot_specifier + '_Waist']/10)**2)
# detuning = ser[mot_specifier + '_Detuning_Load'] # MHz
#
# # Now (perhaps) fit the images, but at least extract the raw counts:
# for motim_i in motim_load:
#     #motim_i.fitImage()
#     motim_i.detFlrNum(atom, intensity, detuning, NA, counts_per_photoelectron,
#                       quantum_efficiency, transmission, t_exposure)
#
# if 'Fit' in motim_load[0].counts.keys():
#     Nat_load = np.array([motim.counts['Fit'] for motim in motim_load])
# else:
#     Nat_load = np.array([motim.counts['Raw'] for motim in motim_load])
#
# # Now, fit the load:
# LoadOneBody = SeqFitter.SeqFitter(Nat_load, t, 'OneBodyLoad')
# #LoadTwoBody = SeqFitter.SeqFitter(Nat_load, t, 'TwoBodyLoad')
#
# # Next, plot up the time trace:
# plt.figure("Loading Curve")
# plt.clf()
# plt.plot(t,  np.array([motim.counts['Raw'] for motim in motim_load])/1e6,
#          '.',color='b',label='flr, cts')
# if 'Fit' in motim_load[0].counts.keys():
#     plt.plot(t,  np.array([motim.counts['Fit'] for motim in motim_load])/1e6,
#              '.',color='b',label='flr, fit')
# plt.plot(t, LoadOneBody.fit_result.best_fit/1e6)
# #plt.plot(t, LoadTwoBody.fit_result.best_fit/1e6)
# #plt.ylim([0, 2*10**6])
# plt.xlabel('$t$ (s)')
# plt.ylabel('$10^6 N$')
# plt.legend(loc='lower right')
#
# # save results
# run.save_result_arrays('t_load', t)
# run.save_result_arrays('N_load', Nat_load)
#
# run.save_result('N_max', Nat_load[-1])
# run.save_result('R', LoadOneBody.fitparams['R'])
# run.save_result('A', LoadOneBody.fitparams['A'])
# run.save_result('gamma',  LoadOneBody.fitparams['gamma'])
