import numpy as np

from ..utilities import domain_utilities, general_utils
from ..utilities.config import ConfigurationError, get_domain_size_from_config
from . import ideal

import pdb


def setup_chemistry(config):
    mechanism = config.get('CHEMISTRY', 'mechanism')
    if mechanism == 'ideal-first-order':
        init_fxn = ideal.init_explicit_first_order_chem_solver
    else:
        raise NotImplementedError('No chemistry mechanism defined for "mechanism" value "{}"'.format(mechanism))

    try:
        return init_fxn(**config.get('CHEMISTRY', 'mechanism_opts'))
    except TypeError as err:
        # Assume that the message will be something along the lines of
        #   foo() missing 1 required positional argument: 'a'
        # We just want the list of missing arguments after the colon
        _, missing_args = err.args[0].split(':')
        raise ConfigurationError('The "{}" mechanism required the following options be given to the "mechanism_opts" '
                                 'configuration line: {}'.format(mechanism, missing_args.strip()))


def get_initial_conditions(config, specie):
    initial_cond = config.get('CHEMISTRY', 'initial_cond')
    if initial_cond == 'zero':
        domain_size = get_domain_size_from_config(config)
        return np.zeros(domain_size)

    elif initial_cond == 'gaussian':
        if len(get_domain_size_from_config(config)) > 1:
            raise NotImplementedError('Not yet configured for > 1D')

        # Since we start at the halfway point, ending at a full box will give the right number of boxes, since arange
        # excluded the stop value
        x_coord, y_coord, z_coord = domain_utilities.compute_coordinates_from_config(config)

        gaussian_opts = config.get('CHEMISTRY', 'initial_cond_opts')
        prefactor = gaussian_opts['height']
        center = gaussian_opts['center']
        sigma = gaussian_opts['width']

        return prefactor * general_utils.gaussian(center, sigma, x_coord, normalized=False)

    else:
        raise NotImplementedError('No method implemented for initial_cond == "{}"'.format(initial_cond))