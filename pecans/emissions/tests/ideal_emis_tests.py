import numpy as np
import os
import sys
import unittest

from ...utilities import config
from ...core import Domain

_my_dir = os.path.abspath(os.path.dirname(__file__))


def _make_test_domain(config_file_name):
    config_obj = config.load_config_file(os.path.join(_my_dir, 'test_config_files', config_file_name))
    domain_obj = Domain(config_obj)
    return domain_obj, config_obj


def _get_box_volume(config_obj):
    dx = config_obj.get('DOMAIN', 'dx')
    dy = config_obj.get('DOMAIN', 'dy')
    dz = config_obj.get('DOMAIN', 'dz')

    return dx * dy * dz


class IdealEmisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.verbose = '-v' in sys.argv or '--verbose' in sys.argv

    def test_simple_point_emis(self):
        """
        Test that a simple emissions case produces the right concentrations without chemistry or transport

        Basically the simplest test than the emissions driver is doing the right thing, put a point source in a small
        domain with no chemistry or transport so that the concentration just builds up and check every so many time steps
        to be sure that it's correct
        """

        test_domain, test_config = _make_test_domain('simple_emissions_test.cfg')
        emis_opts = test_config.get('EMISSIONS', 'emission_opts')
        emis_rate = emis_opts['total']
        dt = test_config.get('DOMAIN', 'dt')
        box_volume = _get_box_volume(test_config)

        test_specie = 'A'
        test_concentrations = []
        check_concentrations = []
        for t in range(100):
            test_domain.step()
            if t % 10 == 0:
                test_concentrations.append(test_domain.species[test_specie][0])
                # So, since the emission rate is molecules/second, after n seconds, the concentration should be
                # E * n / V, where E is the emissions rate and V the box volume. Since there's no chemistry and no
                # transport, none of it should ever be removed.
                check_concentrations.append(emis_rate * dt * (t + 1) / box_volume)

        if self.verbose:
            print('\n  E = {} molec/s, V = {} m^3, dt = {} s'.format(emis_rate, box_volume, dt))
            print('  Expected concentrations:', check_concentrations)
            print('  Model concentrations:', test_concentrations)
            print('  Model emissions at time step {}: '.format(t), test_domain._emissions[test_specie])

        self.assertTrue(np.allclose(test_concentrations, check_concentrations), 'Simple emissions led to difference concentrations than expected')