import copy
import numpy as np
import os
import sys
import unittest

from ...utilities import config
from ...main import Domain

import pdb

_my_dir = os.path.abspath(os.path.dirname(__file__))


def _make_test_domain(config_file_name, config_only=False):
    config_obj = config.load_config_file(os.path.join(_my_dir, 'test_config_files', config_file_name))
    if not config_only:
        domain_obj = Domain(config_obj)
    else:
        domain_obj = None

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

    def test_total_gaussian_emis(self):
        """
        Verify that the sum of emissions across a domain is equal to the requested total emissions
        """
        def compute_total_domain_emissions(the_domain):
            area = the_domain._options['dx'] * the_domain._options['dy']
            return np.sum(the_domain._emissions['A']) * area

        _, base_test_config = _make_test_domain('total_gaussian_emis_test.cfg', config_only=True)
        zeroed_dimensions = {'1D': ['y', 'z'], '2D': ['z'], '3D': []}
        for dimensionality, z_dims in zeroed_dimensions.items():
            with self.subTest(dimensionality=dimensionality):
                # we can get away with a shallow copy b/c we just need to change ny/nz, which are just integers
                test_config = copy.copy(base_test_config)
                expected_total_emis = test_config.get('EMISSIONS', 'emission_opts')['total']
                for dim in z_dims:
                    test_config.set('DOMAIN', 'n'+dim, 0)
                test_domain = Domain(test_config)
                test_domain.execute()
                domain_total = compute_total_domain_emissions(test_domain)
                self.assertTrue(np.allclose(domain_total, expected_total_emis),
                                msg='Sum of {} domain emissions is different than the specified total (sum = {:.3g} * '
                                    'expected)'.format(dimensionality, domain_total / expected_total_emis))
                if self.verbose:
                    print('\n  {}: Expected total = {:.3g}, domain total = {:.3g}'.format(
                        dimensionality, expected_total_emis, domain_total)
                    )