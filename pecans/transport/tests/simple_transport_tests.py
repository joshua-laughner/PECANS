import copy
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import unittest

from ...utilities import config, domain_utilities, general_utils
from ...main import Domain

_my_dir = os.path.abspath(os.path.dirname(__file__))


def _make_test_domain(config_file_name, config_only=False):
    config_obj = config.load_config_file(os.path.join(_my_dir, 'test_config_files', config_file_name))
    if not config_only:
        domain_obj = Domain(config_obj)
    else:
        domain_obj = None
    return domain_obj, config_obj


def _make_plot_filename(basename):
    plot_dir = os.path.join(_my_dir, 'test_output')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    return os.path.join(plot_dir, basename)


class IdealTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.verbose = '-v' in sys.argv or '--verbose' in sys.argv

    def plot_2d(self, data, savename):
        fig = plt.figure()
        ax = fig.add_axes(plt.axes())
        mesh = ax.pcolormesh(data)
        fig.colorbar(mesh)
        # pcolormesh seems to flip the dimensions
        plt.xlabel('Y')
        plt.ylabel('X')
        fig.savefig(_make_plot_filename(savename))

    def test_simple_1d_transport(self):
        """
        Starting from a Gaussian initial condition, do the transport and check how close to the ideal case it is.

        For a pure advection case, with wind U (in meters/second), if the initial case is a Gaussian with width
        :math:`\sigma` and center :math:`\mu`, after :math:`n` seconds, we'd expect the result to be a new Gaussian
        with center :math:`\mu + nU` and the same width. Since we must include diffusion, let us assume that we can
        add the effect of diffusion on top of that, basically that Strang operator splitting holds. Since the diffusion
        operator is:

        .. math::
            \frac{\partial c}{\partial t} = D \frac{\partial^2 c}{\partial x^2}

        if we calculate the second derivative, we get (thanks to Wolfram Alpha)

        .. math::
            G(x | \mu, \sigma) = \exp\left( -\frac{x - mu}{2\sigma^2} \right)

            \frac{\partial^2 c}{\partial x^2} G(x) = \frac{-\sigma^2 + x^2 + \mu^2 - 2x\mu}{\sigma^4} G(x | \mu, \sigma)

        so, in theory, the final result after :math:`n` seconds should be a Gaussian defined by:

        .. math::
            G'(x | \mu', \sigma) = \left(1 + \frac{-\sigma^2 + x^2 + \mu'^2 - 2x\mu'}{\sigma^4} \right) G(x | \mu', \sigma)

        where :math:`\mu'` = :math:`\mu + nU`.
        """

        _, base_config = _make_test_domain('simple_transport_test.cfg')
        zeroed_dimensions = {'1D': ['y', 'z'], '2D': ['z']}
        for dimensionality, zeroed_dims in zeroed_dimensions.items():
            with self.subTest(dimensionality=dimensionality):
                subtest_config = copy.copy(base_config)
                for dim in zeroed_dims:
                    subtest_config.set('DOMAIN', 'n'+dim, 0)
                domain = Domain(subtest_config)
                x, y, _ = domain_utilities.compute_coordinates_from_config(subtest_config)

                if dimensionality == '1D':
                    fig = plt.figure()
                    ax = fig.add_axes(plt.axes())
                    ax.plot(x/1000, domain.species['A'], 'r:', label='Model initial')
                elif dimensionality == '2D':
                    self.plot_2d(domain.species['A'], 'domain_gaussian_2d_transport_init.png')

                # Get the wind speed, diffusion coefficient, and run time from the config
                run_time = subtest_config.get('DOMAIN', 'run_time')
                wind_dict = subtest_config.get('TRANSPORT', 'wind_speeds')
                x_wind = wind_dict['x']
                if wind_dict['y'] != 0.0 or wind_dict['z'] != 0.0:
                    raise RuntimeError('y and z winds should be 0')

                diffusion_coeff_dict = subtest_config.get('TRANSPORT', 'diffusion_coeffs')
                x_diff = diffusion_coeff_dict['x']
                if diffusion_coeff_dict['y'] != 0.0 or diffusion_coeff_dict['z'] != 0.0:
                    raise RuntimeError('y and z diffusion coefficients should be 0')

                # Calculate the ideal final distribution of the species
                init_condition_opts = subtest_config.get('CHEMISTRY', 'initial_cond_opts')
                mu = init_condition_opts['center_x']
                sigma = init_condition_opts['width_x']
                height = init_condition_opts['height']
                target_gaussian = height * general_utils.gaussian(mu + x_wind * run_time, sigma, x, normalized=False)

                # Calculate and add in the effect of diffusion
                diffusion_prefactor = (-sigma**2 + x**2 + mu**2 - 2*x*mu) / (sigma**4)
                target_gaussian = (1 + diffusion_prefactor) * target_gaussian

                # Run the model, then compare the results
                domain.execute()

                simulated_gaussian = domain.species['A']
                rel_diff = np.abs(simulated_gaussian - target_gaussian)/target_gaussian

                # If both the simulated and target value are less than machine precision, we risk getting very large relative
                # differences because of tiny differences between very tiny numbers. So if both are this small, just set the
                # relative difference to 0
                eps = np.finfo(np.float).eps
                is_below_precision = np.logical_and(np.abs(simulated_gaussian) < eps, np.abs(target_gaussian) < eps)
                rel_diff[is_below_precision] = 0.0

                if dimensionality == '1D':
                    ax.plot(x/1000, target_gaussian, 'b', label='Target')
                    ax.plot(x/1000, simulated_gaussian, 'r--', label='Simulated')
                    plt.xlabel('km')
                    plt.ylabel('[A]')
                    ax.legend()
                    fig.savefig(_make_plot_filename('domain_gaussian_1d_transport.png'))

                    fig_diff = plt.figure()
                    ax = fig_diff.add_axes(plt.axes())
                    ax.plot(x/1000, rel_diff, 'r')
                    plt.xlabel('km')
                    plt.ylabel('(sim - target)/target')
                    fig_diff.savefig(_make_plot_filename('domain_gaussian_1d_transport_diff.png'))
                elif dimensionality == '2D':
                    self.plot_2d(domain.species['A'], 'domain_gaussian_2d_transport_final.png')
                    self.plot_2d(rel_diff, 'domain_gaussian_2d_transport_diff.png')

                if self.verbose:
                    print('Simulated result:', simulated_gaussian)
                    print('Target result:', target_gaussian)
                    print('Relative difference:', rel_diff)
                    print('Absolute difference:', np.abs(simulated_gaussian - target_gaussian))

                # Might make this stricter later, for now, want this to be within 1%
                max_total_diff = 0.01
                total_ref_diff = np.sum(rel_diff)
                self.assertLess(total_ref_diff, max_total_diff, 'The total difference between the simulated gaussian is greater than '
                                                                '{}% (= {}%). Right now this is due to dispersion error and large relative '
                                                                'differences caused by very small values.'.format(max_total_diff*100, total_ref_diff*100))