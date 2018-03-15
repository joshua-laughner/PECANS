from copy import copy
import numpy as np
import unittest

from ..transport_utils import check_transport_inputs, StencilPoint


class TestCheckInputs(unittest.TestCase):
    def setUp(self):
        # Create default args to check_transport_inputs. We'll override specific inputs in each test.
        self.null_kwargs = {'dt': None, 'dx': None, 'dy': None, 'u_x': None, 'u_y': None, 'u_z': None,
                       'D_x': None, 'D_y': None, 'D_z': None, 'domain_size': None, 'boundary_conditions': None}
        # These arguments should work if given to check_transport_inputs, so we can override individual ones to
        # test specific error catches
        self.scalar_kwargs = {'dt': 1.0, 'dx': 1000.0, 'dy': 1000.0, 'u_x': 10.0, 'u_y': 10.0, 'u_z': 10.0,
                       'D_x': 10.0, 'D_y': 10.0, 'D_z': 10.0, 'domain_size': (50, 50, 30), 'boundary_conditions': None}

    def test_nonscalar_dt_dx_dy_dz(self):
        # Check that a TypeError is raised if dt, dx, dy, or dz are not scalars
        for variable in ('dt', 'dx', 'dy', 'dz'):
            with self.subTest(variable=variable):
                with self.assertRaises(TypeError):
                    kwargs = copy(self.scalar_kwargs)
                    kwargs[variable] = np.zeros((5,5))
                    check_transport_inputs(**kwargs)


class TestStencilPoint(unittest.TestCase):
    def test_stencil_point_equality(self):
        pt = StencilPoint((1, 2), 1.0)
        pt_diff_coeff = StencilPoint((1, 2), 2.0)
        pt_diff_indices = StencilPoint((1, 1), 1.0)
        pt_diff = StencilPoint((0, 0), 1.0)

        self.assertEqual(pt, pt, msg='Checking that the same StencilPoint instance equals itself failed')
        self.assertEqual(pt, pt_diff_coeff, msg='Checking that two StencilPoints with the same indices are equal failed')
        self.assertNotEqual(pt, pt_diff_indices, msg='Checking that two StencilPoints with different indices are unequal failed')
        self.assertNotEqual(pt, pt_diff, msg='Checking that two StencilPoints with different indices and coefficients are unequal failed')

if __name__ == '__main__':
    unittest.main()