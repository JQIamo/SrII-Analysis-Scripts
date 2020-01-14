import numpy as np
import h5py
from scipy.optimize import curve_fit
from Subroutines.FitFunctions import gauss


def integrated_gaussian_fit_sub(optical_depth):

    # Find Size of Image
    x_length = np.shape(optical_depth)[1]
    z_length = np.shape(optical_depth)[0]

    z = np.arange(z_length)
    x = np.arange(x_length)

    # Get Center
    with h5py.File('current_roi.h5', 'a') as f:
        if 'center' not in f.attrs:
            f.attrs.create('center', (x_length/2, z_length/2))
        center = f.attrs['center']

    # Make Initial Guess
    initial_guess_x = (2, center[0], 150, 0)
    initial_guess_z = (2, center[1], 150, 0)

    # Integrate Image
    image_x = np.sum(optical_depth, axis=0)
    image_z = np.sum(optical_depth, axis=1)

    # Fit Images
    try:
        p_opt_z, p_cov_z = curve_fit(gauss, z, image_z, p0=initial_guess_z,
                                     bounds=([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
        p_opt_x, p_cov_x = curve_fit(gauss, x, image_x, p0=initial_guess_x,
                                     bounds=([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
    except RuntimeError:
        p_opt_x = (0, 0, 0, 0)
        p_opt_z = (0, 0, 0, 0)
        pass

    # Save New Center
    with h5py.File('current_roi.h5', 'a') as f:
        if p_opt_x[0] > 0.2:
            f.attrs['center'] = (p_opt_x[1], p_opt_z[1])

    # Fit Found Something Artificial
    if p_opt_x[2] < 20 or p_opt_z[2] < 20:
        p_opt_x = (0, 0, 1, 0)
        p_opt_z = (0, 0, 1, 0)
    # Return Everything
    return[x, image_x, p_opt_x, z, image_z, p_opt_z]
