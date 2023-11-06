"""
The ``chem_solver`` module contains the functions necessary to initialize and solve idealized chemistry cases.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

from ..utilities.config import ConfigurationError
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
        from ..chemderiv import rhs, mech_species
    except ImportError:
        print('HINT: Your chemderiv module does not exist or does not include an `rhs` and `mech_species` function - '
              'try rebuilding.', file=sys.stderr)
        raise

    mech_species = mech_species()

    # TODO: this and chem_setup.get_constant_params_and_species are inconsistent. Figure out how I'm going to do this,
    #  document it, and fix.
    # Make sure that the driver knows not to update the constant and forced species/parameters
    const_and_forced_params = []
    if 'const_params' in config['CHEMISTRY']:
        const_and_forced_params.extend(config['CHEMISTRY']['const_params'].keys())
    if 'forced_params' in config['CHEMISTRY']:
        const_and_forced_params.extend(config['CHEMISTRY']['forced_params'].keys())

    if 'TEMP' not in const_and_forced_params or 'CAIR' not in const_and_forced_params:
        raise ConfigurationError('Temperature (TEMP) and number density of air (CAIR) must be provided as constant or '
                                 'fixed parameters')

    species_out = tuple([specie for specie in mech_species if specie not in const_and_forced_params])
    print('Model will solve for concentrations of: {}'.format(', '.join(species_out)))
    print('The following species or parameters will be constant or fixed: {}'.format(', '.join(const_and_forced_params)))

    def chem_solver(dt: float, const_param: dict, forced_param: dict, species_in: dict) -> dict:
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
