"""
The ``ideal`` module contains the functions necessary to initialize and solve idealized chemistry cases. These are very
simplified chemical models that, rather that explicitly simulating a full chemical mechanism, instead impose some ideal
kinetic model on just a few species. For example, one is a simple first-order loss, where it is just specified that the
only chemical specie in the model has some lifetime.
"""
from typing import Sequence

# Any init methods need to return two values: the correct chem_solver function to call for each time step, which
# just have the proper signature (dt, TEMP, CAIR, species...), and a tuple of species names that should be initialized
# for that model

from ..utilities.config import ConfigurationError
from ..utilities.chem_utilities import MechanismInterface
from ..utilities import domain_utilities

# TODO: make species definable elsewhere so that this could e.g. have multiple species with different emissions


def _print_unused_mech_opts(extra_options: Sequence[str]):
    """
    Helper function that prints out unused mechanism_opt values

    :param extra_options: the list of extra keyword arguments passed into the mechanism init function

    :return: none
    """
    for key in extra_options:
        print('  mechanism_opt "{}" not needed for this mechanism'.format(key))


def init_explicit_first_order_chem_solver(config: dict) -> MechanismInterface:
    """
    Initialization function for the idealized first-order chemistry solver.

    This requires only one option in the "mechanism_opts" section, namely "lifetime_seconds". That value must be
    the first-order lifetime in seconds, i.e. how long it would take the concentration to decrease to :math:`1/e`
    of its original value. It can take a second, optional option "species_name" which changes the name of the specie
    in the output file.

    :param config: the configuration dictionary.

    :return: the information needed to run this mechanism
    """
    recognized_opts = ('lifetime_seconds', 'species_name')
    unused_options = [k for k in config['CHEMISTRY']['mechanism_opts'] if k not in recognized_opts]
    _print_unused_mech_opts(unused_options)

    lifetime_seconds = float(config['CHEMISTRY']['mechanism_opts']['lifetime_seconds'])
    species_name = config['CHEMISTRY']['mechanism_opts'].get('species_name', 'A')
    if not isinstance(species_name, str):
        raise ConfigurationError('The "species_name" mechanism_opt must be a species names, as a string')
    else:
        species = (species_name,)

    # TODO: unify setup and solver with the new interface for the compiled solvers (this might be done, just need to test this mechanism)
    def chem_solver(dt, _const_params: dict, _forced_params: dict, species_in: dict):
        dt = float(dt)
        for specie, conc in species_in.items():
            # This ideal case assumes that all species are lost with the same first order rate constant, so
            # dC/dt = k*C, and k = 1/lifetime converted to seconds.
            species_in[specie] += -1 / lifetime_seconds * conc * dt

        return species_in

    return MechanismInterface(chem_solver=chem_solver, species=species)


def init_explicit_two_phases_first_order_chem_solver(config: dict):
    """
    Initialization function for the idealized two phases first-order chemistry solver.

    :param first_lifetime_seconds: the first-order lifetime in seconds at the first phase, i.e. how long it would
        take the concentration to decrease to :math:`1/e` of its original value.
    :type first_lifetime_seconds: int or float

    :param second_lifetime_seconds: the second-order lifetime in seconds at the second phase, i.e. how long it would
        take the concentration to decrease to :math:`1/e` of its original value.
    :type second_lifetime_seconds: int or float

    :param first_phase_duration: the cutoff distance between first phase and second phase, i.e the distance
        interval it experiences the first phase lifetime in the unit of meter.
    :type first_phase_duration: int or float

    :param species_name: optional, the name that the specie in this mechanism will be referred to by. Default is "A".
        Changing this has no real effect on the mechanism, just what it is called in the output.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    req_opts = ('first_lifetime_seconds', 'second_lifetime_seconds', 'first_phase_width', 'species_name')
    unused_opts = [k for k in config['CHEMISTRY']['mechanism_opts'] if k not in req_opts]
    _print_unused_mech_opts(unused_opts)

    first_lifetime_seconds = float(config['CHEMISTRY']['mechanism_opts']['first_lifetime_seconds'])
    second_lifetime_seconds = float(config['CHEMISTRY']['mechanism_opts']['second_lifetime_seconds'])
    first_stage_width = float(config['CHEMISTRY']['mechanism_opts']['first_stage_width'])
    species_name = config['CHEMISTRY']['mechanism_opts'].get('species_name', 'A')

    # This mechanism is unusual in that it requires an emission center to know where to start the second lifetime
    try:
        emissions_center_x = config['EMISSIONS']['center_x']
    except KeyError:
        raise ConfigurationError('Running the ideal two-stage mechanism requires an emissions center to be defined, '
                                 'even if emissions are turned off (missing key EMISSIONS -> '
                                 'center_x)')

    # Then it also needs the x-coordinates so that it can determine whether each grid cell is in the first or second
    # stage
    x_coord, _, _ = domain_utilities.compute_coordinates_from_config(config, as_vectors=False)

    if not isinstance(species_name, str):
        raise ConfigurationError('The "species_name" mechanism_opt must be a species names, as a string')
    else:
        species = (species_name,)

    def chem_solver(dt, _const_params: dict, _forced_params: dict, species_in: dict):
        dt = float(dt)
        for specie, conc in species_in.items():
            # This ideal case assumes that all species are lost with the same two first order rate constants, so
            # dC/dt = k*C, and k = 1/lifetime converted to seconds for each half of the domain.
            species_in[specie] += -1 / first_lifetime_seconds * conc * dt * \
                                  (x_coord <= first_stage_width + emissions_center_x) + \
                                  -1 / second_lifetime_seconds * conc * dt * \
                                  (x_coord > first_stage_width + emissions_center_x)

        return species_in

    # Because we captured the emission center and x coordinate in the chem_solver closure, we don't need
    # to return them as additional constant parameters
    return MechanismInterface(chem_solver=chem_solver, species=species)
