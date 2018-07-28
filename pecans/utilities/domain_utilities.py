import numpy as np
from .config import BetterConfig, get_domain_size_from_config

import pdb

def compute_coordinates_from_config(config):
    """
    Calculates the x, y, and z coordinate vectors in kilometers from the origin based on the configuration.

    :param config: the BetterConfig instance representing the model configuration
    :type config: :class:`~pecans.utilities.config.BetterConfig`

    :return: the x, y, and z coordinates as 1D numpy arrays
    :rtype: three :py:class:`~numpy.ndarray`
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
            coords.append(None)

    return tuple(coords)


def is_1D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny == 0 and nz == 0


def is_2D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny > 0 and nz == 0


def is_3D(config_or_domain):
    nx, ny, nz = _get_domain_size(config_or_domain)
    return nx > 0 and ny > 0 and nz > 0


def _get_domain_size(config_or_domain):
    if isinstance(config_or_domain, BetterConfig):
        return get_domain_size_from_config(config_or_domain, all_dims=True)
    else:
        # Assume it is a Domain - cannot import Domain because that creates a circular dependency
        return get_domain_size_from_config(config_or_domain._config, all_dims=True)

