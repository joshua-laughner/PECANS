from collections import OrderedDict
import numpy as np

from .config import BetterConfig, get_domain_size_from_config
from . import general_utils

import pdb


def compute_coordinates_from_config(config, as_vectors=True):
    """
    Calculates the x, y, and z coordinate vectors in kilometers from the origin based on the configuration.

    :param config: the BetterConfig instance representing the model configuration
    :type config: :class:`~pecans.utilities.config.BetterConfig`

    :param as_vectors: whether to return the coordinates as vectors (True, default) or full arrays (False)
    :type as_vectors: bool

    :return: the x, y, and z coordinates as 1D numpy arrays
    :rtype: three :py:class:`numpy.ndarray`
    """
    def coord_helper(n, d):
        # Since we start at the halfway point, ending at a full box will give the right number of boxes, since arange
        # excluded the stop value
        return np.arange(d / 2, n * d, d)

    dimensions = ['x', 'y', 'z']
    coords = []
    for dim in dimensions:
        d_dim = config.get('DOMAIN', 'd' + dim)
        n_dim = config.get('DOMAIN', 'n' + dim)
        if n_dim > 0:
            coords.append(coord_helper(n_dim, d_dim))
        else:
            # If the number of boxes in a dimension is 0, that signifies that that dimension is not used. This can,
            # e.g., simplify the transport equations, or change how the initial concentrations or ideal emissions are
            # computed. Setting the coordinate to None indicates that the model does not need to calculate any processes
            # occurring along that dimension
            coords.append(None)

    if not as_vectors:
        coords = coord_vecs_to_arrays(*coords)
    else:
        coords = tuple(coords)

    return coords


def coord_vecs_to_arrays(x, y=None, z=None):
    """
    Convert coordinate vectors to arrays

    In some cases it makes sense to represent the model domain coordinates as vectors, since the domain boxes follow
    cartesian coordinates, i.e. all boxes along M[0,:,:] will have the same x-coordinate. However, in some cases we
    need an individual coordinate for every box. This converts the coordinate vectors to the latter representation.

    :param x: the x-coordinates, as a vector
    :type x: :py:class:`numpy.ndarray`

    :param y: the y-coordinates, as a vector
    :type y: :py:class:`numpy.ndarray`

    :param z: the z-coordinates, as a vector
    :type z: :py:class:`numpy.ndarray`

    :return: the coordinates in array form
    :rtype: three :py:class:`numpy.ndarray`
    """

    # Check that all inputs are numpy vectors or None. Use an OrderedDict so that later we can convert the values into
    # a list to pass into meshgrid and it'll keep the proper order
    coords = OrderedDict(x=x)
    if y is not None:
        coords['y'] = y
    if z is not None:
        coords['z'] = z

    for k, v in coords.items():
        if not isinstance(v, np.ndarray):
            raise TypeError('The {} coordinate is not a numpy array'.format(k))
        if v.ndim != 1:
            raise ValueError('The {} coordinate is not a vector'.format(k))

    # Using ij indexing means that the new coordinate matrices will be such that new coord[0][:,0,0] == old coord[0],
    # new coord[1][0,:,0] == old coord[1], etc. Without this, the default flips the first two dimensions, which
    # often makes sense from a plotting point of view, but not a coordinate one.
    # Convert the list to a tuple so that tuple expansion can happen, i.e. x, y, z = coord_vecs_to_arrays()
    # will put each coordinate into the appropriate variable.
    coord_vectors = [v for v in coords.values()]
    coord_arrays = np.meshgrid(*coord_vectors, indexing='ij')

    return general_utils.ensure_n_args_to_return(3, coord_arrays)


def coord_arrays_to_vecs(X, Y=None, Z=None):
    """
    Convert coordinate arrays to vectors

    The reverse operation of :func:`coord_vecs_to_arrays`, this takes representations of the model coordinates as arrays
    (where there is an individual model coordinate for every model box) and returns the vector form.

    :param X: the x-coordinate array
    :type X: :py:class:`numpy.ndarray`

    :param Y: the y-coordinate array
    :type Y: :py:class:`numpy.ndarray`

    :param Z: the z-coordinate array
    :type Z: :py:class:`numpy.ndarray`

    :return: the coordinate vectors
    :rtype: three :py:class:`numpy.ndarray` vectors
    """
    vectors = [X[:, 0, 0]]
    if Y is not None:
        vectors.append(Y[0, :, 0])
    if Z is not None:
        vectors.append(Z[0, 0, :])

    return tuple(vectors)


def is_1D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny == 0 and nz == 0


def is_2D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny > 0 and nz == 0


def is_at_least_2D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny > 0


def is_3D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny > 0 and nz > 0


def get_domain_dimensionality(config_or_domain):
    if is_1D(config_or_domain):
        return 1
    elif is_2D(config_or_domain):
        return 2
    elif is_3D(config_or_domain):
        return 3
    else:
        raise NotImplementedError('Model is not 1D, 2D, or 3D')


def _get_domain_size(config_or_domain):
    if isinstance(config_or_domain, BetterConfig):
        return get_domain_size_from_config(config_or_domain, all_dims=True)
    else:
        # Assume it is a Domain - cannot import Domain because that creates a circular dependency
        return get_domain_size_from_config(config_or_domain._config, all_dims=True)

