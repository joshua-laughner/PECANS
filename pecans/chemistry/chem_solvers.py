"""
The ``chem_solver`` module contains the functions necessary to initialize and solve idealized chemistry cases.
"""
import os
from ..utilities.config import ConfigurationError
from ..chemderiv import rhs
import numpy as np
from scipy.integrate import odeint
import time
_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
mech_dir = os.path.join(_mydir, '..', 'Mechanisms')

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

def init_explicit_nox_chem_solver(config):
    """
    Initialization function for the simplified nox-o3 chemistry solver.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    species_file = os.path.join(mech_dir, config.get('CHEMISTRY', 'mechanism') + '.spc')
    species = _parse_pecan_species(species_file)
    if 'const_species' in config.section_as_dict('CHEMISTRY'):
        const_species = config.get('CHEMISTRY', 'const_species')
        species = tuple([specie for specie in species if specie not in const_species.keys()])
    def chem_solver(dt, param, **species_in):
        dt = float(dt)
        lens = [len(conc) for species, conc in species_in.items()][0]
        species = list(species_in.keys())
        init_time = time.time()
        for i in range(lens):
            this_conc = [conc[i] for species, conc in species_in.items()]
            sol = odeint(rhs, np.array(this_conc), np.array([0, dt]), args=(np.array(param),))
            for j in range(len(species)):
                # if this_species not in const_species.keys():
                species_in[species[j]][i] = sol[-1,j]
        #print("The chemistry solver took: %.6s seconds in this time step" % (time.time() - init_time))
        return species_in

    return chem_solver, species
