import os.path

import netCDF4 as ncdf
import numpy as np

from ..utilities import domain_utilities, general_utils
from ..utilities.config import ConfigurationError, get_domain_size_from_config, list_missing_subopts
from ..utilities.chem_utilities import MechanismInterface
from . import ideal
from . import chem_solvers


def setup_chemistry(config: dict) -> MechanismInterface:
    """
    Return the driver function that, when called, will calculate the change in concentrations due to chemistry.

    :param config: the configuration object. Must include the option "mechanism" in the "CHEMISTRY" section

    :return:
        1) the driver function. All driver functions must be called with dt, temperature, and number density of air
        followed by keyword-value pairs of all the chemical species in the mechanism.
        2) the tuple of species names required by the mechanism.
    :rtype: function, tuple of str
    """

    mechanism = config['CHEMISTRY']['mechanism']
    # Look up the right initialization function for the mechanism, we'll call it later in the try-except block to catch
    # cases where not enough mechanism options were provided. All init functions should use the **kwargs syntax to
    # consume extra mechanism options that do not apply to them.
    if mechanism == 'ideal_first_order':
        mech_info = ideal.init_explicit_first_order_chem_solver(config)
    elif mechanism == 'ideal_two_phases_first_order':
        mech_info = ideal.init_explicit_two_phases_first_order_chem_solver(config)
    elif 'nox' in mechanism:
        mech_info = chem_solvers.init_explicit_nox_chem_solver(config)
    else:
        raise NotImplementedError('No chemistry mechanism defined for "mechanism" value "{}"'.format(mechanism))

    const_params = get_constant_params_and_species(config)
    const_params.update(mech_info.const_params)

    forced_params = dict()  # eventually I will implement forced parameters that can change over time
    forced_params.update(mech_info.forced_params)

    # This allows us to check that e.g. temperature and number density of air are provided
    missing_params = [k for k in mech_info.required_params if k not in const_params and k not in forced_params]
    if missing_params:
        raise ConfigurationError(f'The following species/parameters are required to be specified in the constant '
                                 f'or forced parameters sections of the chemistry configuration, but were missing: '
                                 f'{", ".join(missing_params)}')

    return MechanismInterface(
        chem_solver=mech_info.chem_solver,
        species=mech_info.species,
        const_params=const_params,
        forced_params=forced_params
    )


def get_initial_conditions(config, specie):
    """
    Get the initial conditions for a given chemical specie based on how the configuration specifies it

    :param config: the configuration object
    :type config: :class:`~pecans.utilities.BetterConfig`

    :param specie: the name of the chemical specie to load
    :type specie: str

    :return: the array of initial concentrations
    :rtype: :class:`numpy.ndarray`
    """
    initial_cond = config['CHEMISTRY']['initial_cond']
    if initial_cond == 'zero' or specie not in initial_cond.keys():
        domain_size = get_domain_size_from_config(config)
        return np.zeros(domain_size)

    elif initial_cond[specie] == 'gaussian':
        x_coord, y_coord, z_coord = domain_utilities.compute_coordinates_from_config(config, as_vectors=False)
        # Will always need the x values. Append y and z as needed for 2D or 3D models
        gaussian_opts = config['CHEMISTRY']['initial_cond_opts']
        required_subopts = ['height', 'center_x', 'width_x']
        if domain_utilities.is_at_least_2D(config):
            required_subopts.extend(['center_y', 'width_y'])
        if domain_utilities.is_3D(config):
            required_subopts.extend(['center_z', 'width_z'])

        list_missing_subopts(required_subopts, config, 'CHEMISTRY', 'initial_cond_opts', raise_error=True)

        # Any options required here should be added to required_subopts before to verify that they are present and print
        # a useful error message if not
        prefactor = gaussian_opts['height']
        gaussian_kwargs = {'center_x': gaussian_opts['center_x'], 'sigma_x': gaussian_opts['width_x'], 'x': x_coord}
        if domain_utilities.is_at_least_2D(config):
            gaussian_kwargs.update(center_y=gaussian_opts['center_y'], sigma_y=gaussian_opts['width_y'], y=y_coord)
        if domain_utilities.is_3D(config):
            gaussian_kwargs.update(center_z=gaussian_opts['center_z'], sigma_z=gaussian_opts['width_z'], z=z_coord)

        return prefactor * general_utils.gaussian(normalized=False, **gaussian_kwargs)
    elif initial_cond[specie] == 'point':
        coords = domain_utilities.compute_coordinates_from_config(config)

        point_opts = config.get('CHEMISTRY', 'initial_cond_opts')
        if domain_utilities.is_1D(config):
            centers = (point_opts['center_x'],)
        elif domain_utilities.is_2D(config):
            centers = (point_opts['center_x'], point_opts['center_y'])
        elif domain_utilities.is_3D(config):
            centers = (point_opts['center_x'], point_opts['center_y'],  point_opts['center_z'])
        else:
            raise ConfigurationError('Model is not 1D, 2D, or 3D!')

        # zip is smart enough (at least in my testing) that it will stop at the end of the shortest list, so since
        # centers will be only 1, 2, or 3, indices will be the same length
        indices = tuple([np.argmax(np.abs(coord - center)) for coord, center in zip(coords, centers)])
        concentration = np.zeros(get_domain_size_from_config(config))
        concentration[indices] = point_opts['{}_concentration'.format(specie)]
        return concentration
    elif initial_cond[specie] == 'flat':
        domain_size = get_domain_size_from_config(config)
        point_opts = config['CHEMISTRY']['initial_cond_opts']
        return np.zeros(domain_size) + point_opts['{}_concentration'.format(specie)]
    else:
        raise NotImplementedError('No method implemented for initial_cond == "{}"'.format(initial_cond))


def get_constant_params_and_species(config: dict):
    # Constant parameters include chemical species and non-chemical terms like temperature and nair.
    # What makes them constant is that they do not change with time. They can vary in space; to support that,
    # we allow the user to pass "file" as the argument instead of a value, which indicates that it must be read
    # from the const_params_input_file.
    const_params = dict()
    if 'const_params' not in config['CHEMISTRY']:
        # This is not an error; it just means there are no constant parameters. It may be an error if temperature
        # and number density are not given elsewhere, but that is checked outside this function.
        return const_params

    if not isinstance(config['CHEMISTRY']['const_params'], dict):
        raise ConfigurationError('The CHEMISTRY -> const_params option must be a dictionary of species')

    domain_shape = domain_utilities.get_domain_size_from_config(config)
    for param_name, param_value in config['CHEMISTRY']['const_params'].items():
        if isinstance(param_value, (int, float)):
            const_params[param_name] = np.full(domain_shape, float(param_value))
        elif param_value == 'file':
            const_params[param_name] = _get_const_param_from_file(config, param_name, domain_shape)
        else:
            raise ConfigurationError(f'Constant parameter "{param_name}" is neither a number nor the string "file"')

    return const_params


def _get_const_param_from_file(config: dict, param_name: str, domain_shape: tuple):
    try:
        const_param_file = config['CHEMISTRY']['const_param_input_file']
    except KeyError:
        raise ConfigurationError(f'"file" was given as the value for constant parameter {param_name}, but '
                                 f'the "const_param_input_file" option was not listed in the CHEMISTRY section')

    if not os.path.exists(const_param_file):
        raise ConfigurationError(f'Configured const_param_input_file, {const_param_file}, does not exist')

    with ncdf.Dataset(const_param_file) as ds:
        if param_name not in ds.variables.keys():
            raise ConfigurationError(f'Constant parameter {param_name} is not a variable in the constant parameter '
                                     f'input file {const_param_file}')

        param_value = ds[param_name][:]

    if param_value.shape != domain_shape:
        raise ConfigurationError(f'Constant parameter {param_name} in file {const_param_file} has a different array '
                                 f'shape than the domain: domain shape = {domain_shape}, parameter shape = '
                                 f'{param_value.shape}')

    if np.any(param_value.mask):
        raise ConfigurationError(f'Constant parameter {param_name} in file {const_param_file} has fill values - '
                                 f'this is not permitted. All cells must have a valid value.')

    return param_value.filled(0.0)
