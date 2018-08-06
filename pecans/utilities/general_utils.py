import math
import numpy as np


def gaussian(center_x, sigma_x, x, center_y=None, sigma_y=None, y=None,
             center_z=None, sigma_z=None, z=None, normalized=True):
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

    # Check that the inputs are sane - if any of the y or z pieces are given, all must be, and if y is not given, z
    # should not be.
    def dim_term_check(coord, center, sigma):
        check = [coord is None, center is None, sigma is None]
        return any(check) and not all(check)

    if dim_term_check(y, center_y, sigma_y):
        raise ValueError('All or none of y, center_y, and sigma_y must be given')
    if dim_term_check(z, center_z, sigma_z):
        raise ValueError('All or none of z, center_z, and sigma_z must be given')
    if y is None and z is not None:
        raise ValueError('If z is given, y must be as well')

    # Require x, y, and z to all be the same shape
    if (y is not None and y.shape != x.shape) or (z is not None and z.shape != x.shape):
        raise ValueError('x, y, and z must be the same shape')

    # Calculate the Gaussian
    def exponent_term(coord, center, sigma):
        return ((coord - center) ** 2) / (2 * sigma ** 2)

    def prefactor_term(sigma):
        return 1 / math.sqrt(2 * math.pi * sigma ** 2)

    exponent = exponent_term(x, center_x, sigma_x)

    if y is not None:
        exponent += exponent_term(y, center_y, sigma_y)
    if z is not None:
        exponent += exponent_term(z, center_z, sigma_z)

    if not normalized:
        prefactor = 1
    else:
        prefactor = prefactor_term(sigma_x)
        if y is not None:
            prefactor *= prefactor_term(sigma_y)
        if z is not None:
            prefactor *= prefactor_term(sigma_z)

    return prefactor * np.exp(-exponent)


def ensure_n_args_to_return(n_args, args_iterable, fill_val=None):
    """
    Create a tuple of arguments to return that is guaranteed to have a certain number of elements

    Since Python does tuple expansion, you can return a tuple from a function and have each element of the tuple
    expanded into individual return values. For example::

        def return_tuple():
            return ('a','b','c')

        x, y, z = return_tuple()

    would put ``'a'`` into ``x``, ``'b'`` into ``y``, and ``'c'`` into ``z``. However, this means that there must be
    exactly one or three variables to receive the output of return_tuple(), i.e. ``x, y = return_tuple()`` would fail.
    This can be a problem if a function is trying to return a dynamic number of values, since it means the call would
    have to account for that. This function will create a tuple with a guaranteed number of values.

    :param n_args: how many return values you want the tuple to contain
    :type n_args: int

    :param args_iterable: an iterable (currently only tuples and lists supported) that contains <= n_args elements to
        return.
    :type args_iterable: list or tuple

    :param fill_val: optional, the value to append to the returned tuple to fill it out. Default is None.
    :type fill_val: any

    :return: the tuple form of args_iterable with exactly n_args elements
    :rtype: tuple
    """
    if isinstance(args_iterable, list):
        args_iterable = tuple(args_iterable)
    elif not isinstance(args_iterable, tuple):
        raise TypeError('args_iterable must be a tuple or list')

    if len(args_iterable) > n_args:
        raise RuntimeError('More arguments to return ({}) than expected ({})'.format(len(args_iterable), n_args))

    return args_iterable + (fill_val,) * (n_args - len(args_iterable))
