"""
The ``chem_solver`` module contains the functions necessary to initialize and solve idealized chemistry cases.
"""
import os
import sys
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp

from ..utilities.config import ConfigurationError, get_domain_size_from_config
from ..utilities.chem_utilities import MechanismInterface

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
mech_dir = os.path.join(_mydir, '..', 'Mechanisms')
input_dir = os.path.join(_mydir, '..', '..', 'inputs')


def _parse_pecan_species(spec_file):
    """
    Parse a PECANS-style species file. The PECANS format requires that each specie is defined on
    its own line. Anything following a # is considered a comment and ignored
    :param spec_file: The path to the species file, as a str
    :return: nothing. All species are added to Species.instances.
    """
    species_name = []
    with open(spec_file, 'r') as f:
        for line in f:
            try:
                i = line.index('#')
            except ValueError:
                line2 = line.strip()
            else:
                line2 = line[:i].strip()

            if len(line2) > 0:
                species_name.append(line2)
    return tuple(species_name)


def init_explicit_nox_chem_solver(config: dict):
    """
    Initialization function for the simplified nox-o3 chemistry solver.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    try:
        from ..chemderiv import rhs
    except ImportError:
        print('HINT: Your chemderiv module does not exist or does not include an `rhs` function - try rebuilding.',
              file=sys.stderr)
        raise

    # TODO: have the mechanism module include a function that just returns the species so we don't need this
    species_file = os.path.join(mech_dir, config['CHEMISTRY']['mechanism'] + '.spc')
    mech_species = _parse_pecan_species(species_file)
    if 'fixed_params' in config['CHEMISTRY']:
        temp = config['CHEMISTRY']['fixed_params']['temp']
        cair = config['CHEMISTRY']['fixed_params']['nair']
    else:
        msg = 'Temperature and air densities are not provided.'
        raise ConfigurationError(msg)

    # TODO: I see what's going on here - "const_species" are those that get set to a single value throughout
    #  the entire domain and "forced_species" are ones that have spatially-varying values. This seems a bit
    #  confusing to me - I will redo this to combine both of these into one "fixed_species" configuration option
    #  that can either specify a numeric value or the string "file" which means take the concentrations from the
    #  input file. Long term I think we should implement "forced_species" as time varying ones - but those should
    #  have more options (i.e. periodic in time, extrapolated in time, exact in time as well as allowing a single
    #  value for the whole domain per time or taking the whole domain).

    # Make sure that the driver knows not to update the constant and forced species/parameters
    const_and_forced_params = []
    if 'const_params' in config['CHEMISTRY']:
        const_and_forced_params.extend(config['CHEMISTRY']['const_params'].keys())
    if 'forced_params' in config['CHEMISTRY']:
        const_and_forced_params.extend(config['CHEMISTRY']['forced_params'].keys())

    species_out = tuple([specie for specie in mech_species if specie not in const_and_forced_params])

    def chem_solver(dt: float, const_param: dict, forced_param: dict, species_in: dict) -> dict:
        # TODO: change to solve_ive and swap the parameters for a dict (will require modifying mechgen, but
        #  the odds of getting the parameters out of order this way are just too great). Will also need modified
        #  to handle const_param and forced_param being dictionaries now.
        # TODO: have mechgen handle forced/const species correctly
        dt = float(dt)
        grid_shape, num_grid_cells = [(conc.shape, conc.size) for _, conc in species_in.items()][0]

        for i in range(num_grid_cells):
            multi_idx = np.unravel_index(i, grid_shape)
            this_param = {k: v[multi_idx] for k, v in const_param.items()}
            this_param.update({k: v[multi_idx] for k, v in forced_param.items()})
            this_conc = [species_in[name][multi_idx] for name in mech_species]

            solution = solve_ivp(rhs, (0, dt), this_conc, t_eval=dt, args=(this_param,))
            for ispec, specie in enumerate(mech_species):
                # if this_species not in const_species.keys():
                species_in[specie][ispec] = solution.y[ispec, -1]
        return species_in

    return MechanismInterface(chem_solver=chem_solver, species=species_out, required_params=('temp', 'nair'))
