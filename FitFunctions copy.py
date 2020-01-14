import numpy as np
import scipy.constants as constants
import SrConstants


def gaussian_1d(x, a=1, center=100, width=10, offset=0):
    return np.ravel(a * np.exp(-((x - center) / width) ** 2) + offset)


def linear_y_offset(x, *p):
    m, y0 = p
    return m*x + y0


def linear_x_offset(x, *p):
    m, x0 = p
    return m*(x - x0)


def gauss(x, *p):
    a, center, width, offset = p
    return a*np.exp(-(x-center)**2/(2.*width**2)) + offset


def gauss_p_linear(x, *p):
    a, center, width, offset, drift = p
    return a*np.exp(-(x-center)**2/(2.*width**2)) + offset + drift*(x-center)


def double_gauss(x, *p):
    a, center, split, width, offset, drift = p
    return (a * (np.exp(-(x-center+split/2)**2/(2.*width**2)) + np.exp(-(x-center-split/2)**2/(2.*width**2)))
            + offset + drift*(x-center))


def parabola(x, *p):
    a, b, c = p
    return a*((x-b)**2)+c


def exp_decay(x, *p):
    a, t, y0 = p
    return a*np.exp(-x*t)+y0


def temp_fit(x, *p):
    w0, temp = p
    v_rms = np.sqrt(2*constants.value("Boltzmann constant")*temp/SrConstants.mass)
    return np.sqrt(w0**2 + (v_rms*x)**2)
