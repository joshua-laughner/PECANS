import math
import numpy as np


def gaussian(center, sigma, x, y=None, z=None, normalized=True):
    """
    Compute a Gaussian curve. Currently only implemented for 1D

    :param center: the center point in the same coordinates as x, y, and z
    :type center: float

    :param sigma: the sigma width in the same units as x, y, and z
    :type sigma: float

    :param x: the x-coordinate to compute the Gaussian along.
    :type x: :py:class:`numpy.ndarray`

    :param y: the y-coordinate to compute the Gaussian along.
    :type y: :py:class:`numpy.ndarray`

    :param z: the z-coordinate to compute the Gaussian along.
    :type z: :py:class:`numpy.ndarray`

    :param normalized: whether or not to normalize the gaussian so that its area is 1. Default is `True`.
    :type normalized: bool

    :return: an array of size nx-by-ny-by-nz.
    :rtype: :class:`numpy.ndarray`.
    """
    if y is not None or z is not None:
        raise NotImplementedError('Gaussian not implemented for >1D at this time')
    exponent = ((x - center) ** 2) / (2 * sigma ** 2)
    prefactor = 1/math.sqrt(2 * math.pi * sigma**2) if normalized else 1
    return prefactor * np.exp(-exponent)
