import math
import numpy as np


def gaussian(center, sigma, x, y=None, z=None, normalized=True):
    if y is not None or z is not None:
        raise NotImplementedError('Gaussian not implemented for >1D at this time')
    exponent = ((x - center) ** 2) / (2 * sigma ** 2)
    prefactor = 1/math.sqrt(2 * math.pi * sigma**2) if normalized else 1
    return prefactor * np.exp(-exponent)
