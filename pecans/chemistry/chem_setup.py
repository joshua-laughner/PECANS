import numpy as np

from ..utilities import domain_utilities, general_utils
from ..utilities.config import ConfigurationError, get_domain_size_from_config, list_missing_subopts
from . import ideal
from . import chem_solvers


def setup_chemistry(config):
    """
    Return the driver function that, when called, will calculate the change in concentrations due to chemistry.

    :param config: the configuration object. Must include the option "mechanism" in the "CHEMISTRY" section
    :type config: :class:`~pecans.utilities.BetterConfig`

    :return:
        1) the driver function. All driver functions must be called with dt, temperature, and number density of air
        followed by keyword-value pairs of all the chemical species in the mechanism.
        2) the tuple of species names required by the mechanism.
    :rtype: function, tuple of str
    """

    mechanism = config.get('CHEMISTRY', 'mechanism')
    ideal = False;
    # Look up the right initialization function for the mechanism, we'll call it later in the try-except block to catch
    # cases where not enough mechanism options were provided. All init functions should use the **kwargs syntax to
    # consume extra mechanism options that do not apply to them.
    if mechanism == 'ideal_first_order':
        init_fxn = ideal.init_explicit_first_order_chem_solver
        ideal = True
    elif mechanism == 'ideal_two_phases_first_order':
        init_fxn = ideal.init_explicit_two_phases_first_order_chem_solver
        ideal = True
    elif 'nox' in mechanism:
        init_fxn = chem_solvers.init_explicit_nox_chem_solver
    else:
        raise NotImplementedError('No chemistry mechanism defined for "mechanism" value "{}"'.format(mechanism))

    if ideal:
        try:
            return init_fxn(**config.get('CHEMISTRY', 'mechanism_opts'))
        except TypeError as err:
            # Assume that the message will be something along the lines of
            #   foo() missing 1 required positional argument: 'a'
            # We just want the list of missing arguments after the colon
            _, missing_args = err.args[0].split(':')
            raise ConfigurationError('The "{}" mechanism required the following options be given to the "mechanism_opts" '
                                     'configuration line: {}'.format(mechanism, missing_args.strip()))
    #TODO: handle no constant species case
    return init_fxn(config=config)



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
    initial_cond = config.get('CHEMISTRY', 'initial_cond')
    if initial_cond == 'zero' or specie not in initial_cond.keys():
        domain_size = get_domain_size_from_config(config)
        return np.zeros(domain_size)

    elif initial_cond[specie] == 'gaussian':
        x_coord, y_coord, z_coord = domain_utilities.compute_coordinates_from_config(config, as_vectors=False)
        # Will always need the x values. Append y and z as needed for 2D or 3D models
        gaussian_opts = config.get('CHEMISTRY', 'initial_cond_opts')
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
        point_opts = config.get('CHEMISTRY', 'initial_cond_opts')
        return np.zeros(domain_size) + point_opts['{}_concentration'.format(specie)]
    else:
        raise NotImplementedError('No method implemented for initial_cond == "{}"'.format(initial_cond))