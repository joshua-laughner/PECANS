from ..utilities import general_utils, domain_utilities
from ..utilities.config import ConfigurationError

import pdb

def setup_emissions(config):
    _check_grid_box_size(config, 'dy')
    _check_grid_box_size(config, 'dz')

    emis_type = config.get('EMISSIONS', 'emission_type')
    if emis_type == 'gaussian':
        get_emis_fxn = _setup_gaussian_emissions(config)
    else:
        raise NotImplementedError('No emissions set up for "{}" emissions type'.format(emis_type))

    def emissions_solver(config, concentration, specie, seconds_since_model_start):
        dz = config.get('DOMAIN', 'dz')
        emis = get_emis_fxn(specie, seconds_since_model_start)
        return concentration + emis / dz, emis

    return emissions_solver


def _setup_gaussian_emissions(config):

    # So for a Gaussian we multiply the total emissions by a normalized Gaussian so that in theory the integrated
    # emissions will be close to the total specified. (It'll probably be a little off b/c of discretization, this could
    # be improved on later.) Then we divide by area because emissions are usually given in molec./area/second. We assume
    # that the total is molec./second
    dx = config.get('DOMAIN', 'dx')
    dy = config.get('DOMAIN', 'dy')

    x, y, z = domain_utilities.compute_coordinates_from_config(config)

    if y is not None or z is not None:
        raise NotImplementedError('Gaussian emissions not set up for 2 or 3D models')

    emis_opts = config.get('EMISSIONS', 'emission_opts')
    emissions_vector = emis_opts['total'] / (dx * dy) * general_utils.gaussian(emis_opts['center'], emis_opts['width'], x=x, normalized=True)

    def return_gaussian_vector(specie, seconds_since_model_start):
        return emissions_vector

    return return_gaussian_vector


def _check_grid_box_size(config, dim):
    dx = config.get('DOMAIN', 'dx')
    d2 = config.get('DOMAIN', dim)

    # Quick check that the user didn't make dy some very small value in a 1D model.
    msg = '{} is {}. Emissions are calculated using both dx and dy to compute the grid box area, even in a 1D model.'
    if d2 == 0:
        raise ConfigurationError(msg.format(dim, 0) + ' Set dy to a nonzero value (usually == to dx is good).')
    elif d2 <= 0.1 * dx:
        print(msg.format(dim, '< 0.1 * dx') + ' If this is intended, then nothing is wrong, but usually in a 1D model, '
                                              'setting dy approximately equal to dx is the more usual approach.')