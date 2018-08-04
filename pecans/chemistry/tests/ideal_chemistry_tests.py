from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import unittest

from ...utilities import config, domain_utilities, general_utils
from ...main import Domain

_my_dir = os.path.abspath(os.path.dirname(__file__))


def _make_test_domain(config_file_name):
    config_obj = config.load_config_file(os.path.join(_my_dir, 'test_config_files', config_file_name))
    domain_obj = Domain(config_obj)
    return domain_obj, config_obj


def _make_plot_filename(basename):
    plot_dir = os.path.join(_my_dir, 'test_output')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    return os.path.join(plot_dir, basename)


class IdealChemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.verbose = '-v' in sys.argv or '--verbose' in sys.argv

    def test_first_order_chemistry(self):
        domain, config = _make_test_domain('simple_chemistry_test.cfg')
        initial_concentration = domain.species['A'].copy()
        run_time = config.get('DOMAIN', 'run_time')
        lifetime_seconds = config.get('CHEMISTRY', 'mechanism_opts')['lifetime_seconds']
        domain.execute()

        expected_concentration = initial_concentration * np.exp(-run_time/lifetime_seconds)
        simulated_concentration = domain.species['A']

        self.assertTrue(np.allclose(expected_concentration, simulated_concentration, atol=5e-4),
                        'The expected ({}) and simulated ({}) concentrations differ by more than the tolerance'.format(
                            expected_concentration, simulated_concentration
                        ))

        if self.verbose:
            print('After {} sec with tau = {} sec, expected concentration = {}, simulated concentration = {}'.format(
                run_time, lifetime_seconds, expected_concentration, simulated_concentration
            ))