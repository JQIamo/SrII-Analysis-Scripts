import numpy as np
import h5py
from scipy.optimize import curve_fit
from Subroutines.FitFunctions import gauss


def splice_gaussian_fit_sub(optical_depth, fix_center):

    x_length = np.shape(optical_depth)[1]
    z_length = np.shape(optical_depth)[0]

    with h5py.File('current_roi.h5', 'a') as f:
        if 'center' not in f.attrs:
            f.attrs.create('center', (x_length / 2, z_length / 2))
        center = f.attrs['center']

    x_slice = 0 if center[1] < 0 else z_length-1 if center[1] >= z_length else center[1]
    z_slice = 0 if center[0] < 0 else x_length - 1 if center[0] >= x_length else center[0]

    if optical_depth[int(x_slice)][int(z_slice)] < 0.3 and not fix_center:
        for i in range(1000):
            guess = (np.random.randint(x_length - 1), np.random.randint(z_length - 1))
            if optical_depth[guess[1]][guess[0]] > 0.3:
                x_slice = guess[1]
                z_slice = guess[0]
                break

    initial_guess_z = (2, center[1], 150, 0)
    initial_guess_x = (2, center[0], 150, 0)

    z = np.arange(z_length)
    x = np.arange(x_length)

    image_x = optical_depth[int(x_slice)]
    image_z = optical_depth[:, int(z_slice)]

    try:
        for i in range(2):
            p_opt_z, p_cov_z = curve_fit(gauss, z, image_z, p0=initial_guess_z,
                                         bounds=([0, -np.inf, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
            x_slice = 0 if p_opt_z[1] < 0 else z_length-1 if p_opt_z[1] >= z_length else p_opt_z[1]
            image_x = optical_depth[int(x_slice)]
            initial_guess_z = p_opt_z
            p_opt_x, p_cov_x = curve_fit(gauss, x, image_x, p0=initial_guess_x,
                                         bounds=([0, -np.inf, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
            z_slice = 0 if p_opt_x[1] < 0 else x_length - 1 if p_opt_x[1] >= x_length else p_opt_x[1]
            image_z = optical_depth[:, int(z_slice)]
            initial_guess_x = p_opt_x
    except RuntimeError:
        p_opt_x = (0, 0, 0, 0)
        p_opt_z = (0, 0, 0, 0)
        pass

    with h5py.File('current_roi.h5', 'a') as f:
        if p_opt_x[0] > 0.2 and not fix_center:
            f.attrs['center'] = (p_opt_x[1], p_opt_z[1])

    if p_opt_x[2] < 30 or p_opt_z[2] < 30:
        p_opt_x = (0, 0, 1, 0)
        p_opt_z = (0, 0, 1, 0)

    return[x, image_x, p_opt_x, z, image_z, p_opt_z]
