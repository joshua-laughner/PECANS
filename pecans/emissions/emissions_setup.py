import sys
from abc import ABC, abstractmethod
from typing import Callable, Tuple, Dict

import numpy as np

from pecans.utilities import general_utils, domain_utilities
from pecans.utilities.config import ConfigurationError, get_domain_size_from_config, list_missing_subopts
import os

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
input_dir = os.path.join(_mydir, '..', '..', 'inputs')


class EmissionsSolver:
    @abstractmethod
    def solve_emissions(self, config: dict, seconds_since_model_start: float, species_concentrations: dict):
        pass


def setup_emissions(config: dict) -> EmissionsSolver:
    """
    Primary emissions setup method that sets up the proper emissions solver based on the configuration file

    :param config: the configuration object

    :return: the emissions solver function. Any solver function takes as input: the config object, seconds since model
        start, and the dictionary of chemical species as keyword-value pairs. Any solver function returns the updated
        dictionary of chemical species and the current dictionary of emissions.
    :rtype: dict and dict
    """
    _check_grid_box_size(config, 'dy')
    _check_grid_box_size(config, 'dz')

    # Emissions can be set up globally, without reference to specific species, or with different emissions per specie.
    # These two modes are incompatible. The former looks like this in TOML:
    #   [EMISSIONS]
    #   do_emissions = true
    #   emission_type = "point"
    #   center_x = 0.0
    #   center_y = 0.0
    #   center_z = 0.0
    #   total = 1e24
    #
    # The latter would look like:
    #   [EMISSIONS]
    #   do_emission = true
    #   [[EMISSIONS.species]]
    #   emission_specie = "NO"
    #   emission_type = "point"
    #   ...
    #   [[EMISSIONS.species]]
    #   emissions_specie = "NO2"
    #   emission_type = "gaussian"
    #   ...

    if 'species' in config['EMISSIONS'] and 'emission_type' in config['EMISSIONS']:
        raise ConfigurationError('The EMISSIONS section contains both the keys "species" and "emission_type" in the '
                                 'top level. This is not allowed, you must provide either general emissions or '
                                 'per-specie emissions, you cannot have both.')
    elif 'emission_type' in config['EMISSIONS']:
        # Case 1: same emissions for all species,
        get_emis_fxn = _setup_emission_for_given_type(config, config['EMISSIONS'])
        return GenericEmissionsSolver(get_emis_fxn)
    elif 'species' in config['EMISSIONS']:
        # Case 2: per-emission species
        get_emis_fxn = dict()
        for emis_specie_section in config['EMISSIONS']['species']:
            specie = emis_specie_section['emission_specie']
            get_emis_fxn[specie] = _setup_emission_for_given_type(config, emis_specie_section)
        return PerSpeciesEmissionsSolver(get_emis_fxn)
    elif not config['EMISSION']['do_emissions']:
        # Case 3: emissions are turned off. The main code should catch that and not call this,
        # but just in case we return something that is a no-op for emissions.
        return NoEmissionsSolver()
    else:
        raise ConfigurationError('do_emissions was True, but no emissions were configured')



class NoEmissionsSolver(EmissionsSolver):
    def solve_emissions(self, config: dict, seconds_since_model_start: float, species_concentrations: dict):
        domain_shape = domain_utilities.get_domain_size_from_config(config)
        emis_dict = {name: np.zeros(domain_shape) for name in species_concentrations}
        return species_concentrations, emis_dict


class GenericEmissionsSolver(EmissionsSolver):
    def __init__(self, get_emissions_function: Callable[[str, float], np.ndarray]):
        self._get_emissions_function = get_emissions_function

    def solve_emissions(self, config: dict, seconds_since_model_start: float, species_concentrations: dict):
        dt = config['DOMAIN']['dt']
        dz = config['DOMAIN']['dz'] * 1e2  # convert from meters to centimeters
        emis_dict = dict()

        for specie_name, specie_conc in species_concentrations.items():
            emis = self._get_emissions_function(specie_name, seconds_since_model_start)
            species_concentrations[specie_name] = specie_conc + emis / dz * dt
            emis_dict[specie_name] = emis

        return species_concentrations, emis_dict


class PerSpeciesEmissionsSolver(EmissionsSolver):
    def __init__(self, get_emissions_functions: Dict[str, Callable[[str, float], np.ndarray]]):
        self._get_emissions_functions = get_emissions_functions
        self._warned_about_species = set()

    def solve_emissions(self, config: dict, seconds_since_model_start: float, species_concentrations: dict):
        dt = config['DOMAIN']['dt']
        dz = config['DOMAIN']['dz'] * 1e2  # convert from meters to centimeters
        domain_shape = domain_utilities.get_domain_size_from_config(config)
        emis_dict = {name: np.zeros(domain_shape) for name in species_concentrations}

        for specie_name, emis_fxn in self._get_emissions_functions.items():
            if specie_name not in species_concentrations:
                if not specie_name in self._warned_about_species:
                    print(f'WARNING: emissions are defined for specie {specie_name} but that is not one of the species '
                          f'being simulated. You may want to double check your configuration. Note that fixed/forced '
                          f'species cannot have emissions.', file=sys.stderr)
                    self._warned_about_species.add(specie_name)
            else:
                emis = emis_fxn(specie_name, seconds_since_model_start)
                species_concentrations[specie_name] += emis / dz * dt
                emis_dict[specie_name] = emis

        return species_concentrations, emis_dict


def _setup_emission_for_given_type(config: dict, emis_config_subsection: dict):
    emis_type = emis_config_subsection['emission_type']
    if emis_type == 'gaussian':
        return _setup_gaussian_emissions(config, emis_config_subsection)
    elif emis_type == 'point':
        return _setup_point_emissions(config, emis_config_subsection)
    elif emis_type == 'constant':
        return _setup_constant_emissions(config, emis_config_subsection)
    else:
        raise NotImplementedError('No emissions set up for "{}" emissions type'.format(emis_type))


def _setup_gaussian_emissions(config: dict, emis_config_subsection: dict):
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

    x, y, z = domain_utilities.compute_coordinates_from_config(config, as_vectors=False)

    # Check that the necessary options are present in the configuration
    if domain_utilities.is_2D(config) or domain_utilities.is_3D(config):
        use_2d_emis = True
    else:
        use_2d_emis = False

    center_x = emis_config_subsection['total']
    width_x = emis_config_subsection['width_x']
    total = emis_config_subsection['total']

    if not use_2d_emis:
        # Okay, this took far too long to figure out. Originally, I had E_tot / (dx*dy) * G_x
        # but that always gave emissions off by a factor of dx. Why? The gaussian itself has
        # units of m^{-1}, because it's normalized; essentially it's giving the probability
        # per unit length.
        #
        # What we want is \sum_x E_x dx dy = E_tot. That means E_x should have units of mol s^-1 m^-2,
        # if E_tot has units of mol s^-1. That would seem to suggest that E_x = E_tot / (dx*dy) * G_x,
        # except that G_x has units of m^-1, and it needs to be multiplied by dx to get from "fraction of
        # emissions per unit length" to "fraction of emissions in a box of length dx". Which means
        # the dx cancels out, and we only really need to divide by dy to get E_x into molec m^-2 s^-1.
        #
        # We also convert he E_x from molec m^-2 s^-1 to molec cm^-2 s^-1
        # We assume that the total is molec./second, and we want molec./area/second
        dy = config['DOMAIN']['dy']
        emissions_array = total / dy * general_utils.gaussian(center_x, width_x, x=x, normalized=True)
        emissions_array = emissions_array / 1e4

    else:
        # Similarly, in 2D, the Gaussian should be normalized to "probability per unit area". So we would
        # multiply it by dx * dy, which cancels out the dx * dy that we divide the total emissions by to
        # put it in per area.
        emissions_array = total * general_utils.gaussian(
            center_x=center_x, sigma_x=width_x, x=x, center_y=emis_config_subsection['center_y'],
            sigma_y=emis_config_subsection['width_y'], y=y, normalized=True
        )
        emissions_array = emissions_array /1e4
        if domain_utilities.is_3D(config):
            # if 3D, we need to set the non-surface boxes to zero. x and y will have the proper 3D shape, so every layer
            # of emissions_array will be the same. Therefore, we just need to set all the levels above the surface to 0.
            emissions_array[:, :, 1:] = 0

    def return_gaussian_array(specie, seconds_since_model_start):
        return emissions_array

    return return_gaussian_array


def _setup_point_emissions(config, emis_config_subsection):
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
    dx = config['DOMAIN']['dx']
    dy = config['DOMAIN']['dy']

    x, y, z = domain_utilities.compute_coordinates_from_config(config)
    idx_x = np.argmin(np.abs(x - emis_config_subsection['center_x']))
    emis_array = np.zeros(get_domain_size_from_config(config))
    point_emis_rate = emis_config_subsection['total'] / (dx * dy) / 1e4

    if y is None and z is None:
        emis_array[idx_x] = point_emis_rate
    elif z is None:
        idx_y = np.argmin(np.abs(y - emis_config_subsection['center_y']))
        emis_array[idx_x, idx_y] = point_emis_rate
    else:
        idx_y = np.argmin(np.abs(y - emis_config_subsection['center_y']))
        idx_z = np.argmin(np.abs(z - emis_config_subsection['center_z']))
        emis_array[idx_x, idx_y, idx_z] = point_emis_rate

    def return_emis_array(specie, seconds_since_model_start):
        return emis_array

    return return_emis_array


def _setup_constant_emissions(config, emis_config_subsection):
    """
    Helper function that sets up a constant emission source based on the configuration file.

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :return: a function that, when called, returns an array of emissions in molecules cm^-2 s^-1. It accepts two inputs
        (species name and seconds since model start) but doesn't use either - just accepted for consistency with other
        expected "get emissions" functions
    :rtype: function
    """
    dx = config['DOMAIN']['dx']
    dy = config['DOMAIN']['dy']

    emis = np.zeros(get_domain_size_from_config(config)) + emis_config_subsection['total'] / (dx * dy) / 1e4

    def return_emis_array(specie, seconds_since_model_start):
        return emis

    return return_emis_array


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

    dx = config['DOMAIN']['dx']
    d2 = config['DOMAIN'][dim]

    # Quick check that the user didn't make dy some very small value in a 1D model.
    msg = '{dim} is {value}. Emissions are calculated using both dx and dy to compute the grid box area, even in a 1D model.'
    if dx == 0:
        raise ConfigurationError((msg + ' Set dx to a nonzero value.').format(dim='dx', value=0))
    elif d2 == 0:
        raise ConfigurationError((msg + ' Set {dim} to a nonzero value (usually == to dx is good).').format(dim=dim, value=0))
    elif d2 <= 0.1 * dx:
        print(msg.format(dim=dim, value='< 0.1 * dx') + ' If this is intended, then nothing is wrong, but usually in a 1D model, '
                                                        'setting dy approximately equal to dx is the more usual approach.')