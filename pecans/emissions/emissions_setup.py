import numpy as np

from pecans.utilities import general_utils, domain_utilities
from pecans.utilities.config import ConfigurationError, get_domain_size_from_config, list_missing_subopts

import pdb


def setup_emissions(config):
    """
    Primary emissions setup method that sets up the proper emissions solver based on the configuration file

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :return: the emissions solver function. Any solver function takes as input: the config object, seconds since model
        start, and the dictionary of chemical species as keyword-value pairs. Any solver function returns the updated
        dictionary of chemical species and the current dictionary of emissions.
    :rtype: dict and dict
    """
    _check_grid_box_size(config, 'dy')
    _check_grid_box_size(config, 'dz')

    emis_type = config.get('EMISSIONS', 'emission_type')
    if emis_type == 'gaussian':
        get_emis_fxn = _setup_gaussian_emissions(config)
    elif emis_type == 'point':
        get_emis_fxn = _setup_point_emissions(config)
    else:
        raise NotImplementedError('No emissions set up for "{}" emissions type'.format(emis_type))

    def emissions_solver(config, seconds_since_model_start, **species):
        dt = config.get('DOMAIN', 'dt')
        dz = config.get('DOMAIN', 'dz')
        emis_dict = dict()
        for specie, concentration in species.items():
            emis = get_emis_fxn(specie, seconds_since_model_start)
            species[specie] = concentration + emis / dz * dt
            emis_dict[specie] = emis

        return species, emis_dict

    return emissions_solver


def _setup_gaussian_emissions(config):
    """
    Helper function that sets up a Gaussian shaped emission source based on the configuration file.

    Requires that emission_opts in the config file contains:

        * center_x, width_x, total
        * center_y, width_y if > 1D
        * center_z, width_z if > 2D

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :return: a function that, when called, returns an array of emissions in molecules cm^-2 s^-1. It accepts two inputs
        (species name and seconds since model start) but doesn't use either - just accepted for consistency with other
        expected "get emissions" functions
    :rtype: function
    """
    # So for a Gaussian we multiply the total emissions by a normalized Gaussian so that in theory the integrated
    # emissions will be close to the total specified. (It'll probably be a little off b/c of discretization, this could
    # be improved on later.) Then we divide by area because emissions are usually given in molec./area/second.
    dy = config.get('DOMAIN', 'dy')

    x, y, z = domain_utilities.compute_coordinates_from_config(config, as_vectors=False)

    emis_opts = config.get('EMISSIONS', 'emission_opts')

    # Check that the necessary options are present in the configuration
    required_opts = ['total', 'center_x', 'width_x']
    if domain_utilities.is_2D(config) or domain_utilities.is_3D(config):
        required_opts += ['center_y', 'width_y']
        use_2d_emis = True
    else:
        use_2d_emis = False

    list_missing_subopts(required_opts, config, 'EMISSIONS', 'emission_opts', raise_error=True)

    if not use_2d_emis:
        # Okay, this took far too long to figure out. Originally, I had E_tot / (dx*dy) * G_x
        # but that always gave emissions off by a factor of dx. Why? The gaussian itself has
        # units of m^{-1}, because it's normalized; essentially it's giving the probability
        # per unit length.
        #
        # What we want is \sum_x E_x dx dy = E_tot. That means E_x should have units of mol s^-1 m^-2,
        # if E_tot has units of mol s^-1. That would seem to suggest that E_x = E_tot / (dx*dy) * G_x,
        # except that G_x has units of m^-1, and it needs multiplied by dx to get from "fraction of
        # emissions per unit length" to "fraction of emissions in a box of length dx". Which means a
        # the dx cancels out, and we only really need to divide by dy to get E_x into molec m^-2 s^-1.
        #
        # We assume that the total is molec./second and we want molec./area/second
        dy = config.get('DOMAIN', 'dy')
        emissions_array = emis_opts['total'] / dy * general_utils.gaussian(emis_opts['center_x'], emis_opts['width_x'], x=x, normalized=True)
    else:
        # Similarly, in 2D, the Gaussian should be normalized to "probability per unit area". So we would
        # multiply it by dx * dy, which cancels out the dx * dy that we divide the total emissions by to
        # put it in per area.
        emissions_array = emis_opts['total'] * general_utils.gaussian(center_x=emis_opts['center_x'], sigma_x=emis_opts['width_x'], x=x,
                                                                      center_y=emis_opts['center_y'], sigma_y=emis_opts['width_y'], y=y,
                                                                      normalized=True)
        if domain_utilities.is_3D(config):
            # if 3D, we need to set the non-surface boxes to zero. x and y will have the proper 3D shape, so every layer
            # of emissions_array will be the same. Therefore, we just need to set all the levels above the surface to 0.
            emissions_array[:, :, 1:] = 0

    def return_gaussian_vector(specie, seconds_since_model_start):
        return emissions_array

    return return_gaussian_vector


def _setup_point_emissions(config):
    """
    Helper function that sets up a point emission source based on the configuration file.

    Requires that emission_opts in the config file contains center_x, center_y if > 1D, center_z if > 2D, and total.

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :return: a function that, when called, returns an array of emissions in molecules cm^-2 s^-1. It accepts two inputs
        (species name and seconds since model start) but doesn't use either - just accepted for consistency with other
        expected "get emissions" functions
    :rtype: function
    """
    dx = config.get('DOMAIN', 'dx')
    dy = config.get('DOMAIN', 'dy')

    x, y, z = domain_utilities.compute_coordinates_from_config(config)

    if y is not None or z is not None:
        raise NotImplementedError('Point emissions not set up for 2 or 3D models')
    emis_opts = config.get('EMISSIONS', 'emission_opts')

    idx = np.argmin(np.abs(x - emis_opts['center_x']))
    emis = np.zeros(get_domain_size_from_config(config))
    emis[idx] = emis_opts['total'] / (dx * dy)

    def return_emis_vector(specie, seconds_since_model_start):
        return emis

    return return_emis_vector


def _check_grid_box_size(config, dim):
    """
    Helper function that checks that the boxes are a reasonable size.

    Typically, we want to check that no box length in any dimension is 0 or much smaller than the others. In some cases
    it may not matter, but just as a safety precaution, this function helps check that.

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :param dim: which dimension (y or z) to check.
    :type dim: str

    :return: none, raises :class:`~pecans.utilities.config.ConfigurationError` if an invalid configuration is active.
    """

    dx = config.get('DOMAIN', 'dx')
    d2 = config.get('DOMAIN', dim)

    # Quick check that the user didn't make dy some very small value in a 1D model.
    msg = '{dim} is {value}. Emissions are calculated using both dx and dy to compute the grid box area, even in a 1D model.'
    if dx == 0:
        raise ConfigurationError((msg + ' Set dx to a nonzero value.').format(dim='dx', value=0))
    elif d2 == 0:
        raise ConfigurationError((msg + ' Set {dim} to a nonzero value (usually == to dx is good).').format(dim=dim, value=0))
    elif d2 <= 0.1 * dx:
        print(msg.format(dim=dim, value='< 0.1 * dx') + ' If this is intended, then nothing is wrong, but usually in a 1D model, '
                                                        'setting dy approximately equal to dx is the more usual approach.')