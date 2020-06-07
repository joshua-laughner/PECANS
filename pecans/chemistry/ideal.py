"""
The ``ideal`` module contains the functions necessary to initialize and solve idealized chemistry cases. These are very
simplified chemical models that, rather that explicitly simulating a full chemical mechanism, instead impose some ideal
kinetic model on just a few species. For example, one is a simple first-order loss, where it is just specified that the
only chemical specie in the model has some lifetime.
"""

# Any init methods need to return two values: the correct chem_solver function to call for each time step, which
# just have the proper signature (dt, TEMP, CAIR, species...), and a tuple of species names that should be initialized
# for that model

from ..utilities.config import ConfigurationError

# TODO: make species definable elsewhere so that this could e.g. have multiple species with different emissions


def _print_unused_mech_opts(extra_kwargs):
    """
    Helper function that prints out unused mechanism_opt values

    :param extra_kwargs: the dictionary of extra keyword arguments passed into the mechanism init function
    :type extra_kwargs: dict

    :return: none
    """
    for key in extra_kwargs.keys():
        print('  mechanism_opt "{}" not needed for this mechanism'.format(key))


def init_explicit_first_order_chem_solver(lifetime_seconds, species_name='A', **kwargs):
    """
    Initialization function for the idealized first-order chemistry solver.

    :param lifetime_seconds: the first-order lifetime in seconds, i.e. how long it would take the concentration to
        decrease to :math:`1/e` of its original value.
    :type lifetime_seconds: int or float

    :param species_name: optional, the name that the specie in this mechanism will be referred to by. Default is "A".
        Changing this has no real effect on the mechanism, just what it is called in the output.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    _print_unused_mech_opts(kwargs)

    lifetime_seconds = float(lifetime_seconds)
    if not isinstance(species_name, str):
        raise ConfigurationError('The "species" mechanism_opt must be a species names, as a string')
    else:
        species = (species_name,)

    def chem_solver(dt, TEMP, CAIR, **species_in):
        dt = float(dt)
        for specie, conc in species_in.items():
            # This ideal case assumes that all species are lost with the same first order rate constant, so
            # dC/dt = k*C, and k = 1/lifetime converted to seconds.
            species_in[specie] += -1 / lifetime_seconds * conc * dt

        return species_in

    return chem_solver, ('A',)

def init_explicit_two_phases_first_order_chem_solver(first_lifetime_seconds, second_lifetime_seconds,
                                                     first_phase_width, species_name = 'A', **kwargs):
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
    _print_unused_mech_opts(kwargs)

    first_lifetime_seconds = float(first_lifetime_seconds)
    second_lifetime_seconds = float(second_lifetime_seconds)
    first_phase_width = float(first_phase_width)

    if not isinstance(species_name, str):
        raise ConfigurationError('The "species" mechanism_opt must be a species names, as a string')
    else:
        species = (species_name,)

    def chem_solver(dt, TEMP, CAIR, x_coord, E_center, **species_in):
        dt = float(dt)
        for specie, conc in species_in.items():
            # This ideal case assumes that all species are lost with the same first order rate constant, so
            # dC/dt = k*C, and k = 1/lifetime converted to seconds.
            species_in[specie] += -1 / first_lifetime_seconds * conc * dt * \
                                  (x_coord <= first_phase_width+E_center) + \
                                  -1 / second_lifetime_seconds * conc * dt * \
                                  (x_coord > first_phase_width+E_center)

        return species_in

    return chem_solver, ('A',)

