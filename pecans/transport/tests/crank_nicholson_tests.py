import numpy as np
import sys
import unittest

from ...utilities import io_utils
from .. import crank_nicholson as cn

class TestMatrixConstruction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.verbose = '-v' in sys.argv or '--verbose' in sys.argv

        dt = 1.0
        dx = 1000.0
        dy = 1000.0
        dz = 500.0
        u_x = 5.0
        u_y = 10.0
        u_z = 20.0
        D_x = 90.0
        D_y = 110.0
        D_z = 500.0

        cls.r_x = dt * D_x / (dx**2)
        cls.r_y = dt * D_y / (dy**2)
        cls.r_z = dt * D_z / (dz**2)
        cls.C_rx = dt * u_x / (2*dx)
        cls.C_ry = dt * u_y / (2*dy)
        cls.C_rz = dt * u_z / (2*dz)
        cls.model_1d_settings = {'dt': dt, 'dx': dx, 'dy': None, 'dz': None, 'u_x': u_x, 'u_y': None, 'u_z': None,
                                 'D_x': D_x, 'D_y': None, 'D_z': None, 'domain_size': (5,), 'boundary_conditions': None}
        cls.model_2d_settings = {'dt': dt, 'dx': dx, 'dy': dy, 'dz': None, 'u_x': u_x, 'u_y': u_y, 'u_z': None,
                                 'D_x': D_x, 'D_y': D_y, 'D_z': None, 'domain_size': (5, 5), 'boundary_conditions': None}
        if cls.verbose:
            for attr in ('r_x', 'r_y', 'r_z', 'C_rx', 'C_ry', 'C_rz'):
                print('TestMatrixConstruction: {} = {}'.format(attr, getattr(cls, attr)))

    """
    def test_1d_zero_bc_matrix(self):
        A = be.construct_transport_equation_2nd_order_centered_space(**self.model_1d_settings)
        A_check = np.array([[1 + 2*self.r_x,        self.C_rx - self.r_x, 0,                     0,                     0],
                           [-self.C_rx - self.r_x, 1 + 2*self.r_x,        self.C_rx - self.r_x,  0,                     0],
                           [0,                     -self.C_rx - self.r_x, 1 + 2*self.r_x,        self.C_rx - self.r_x,  0],
                           [0,                     0,                     -self.C_rx - self.r_x, 1 + 2*self.r_x,        self.C_rx - self.r_x],
                           [0,                     0,                     0,                     -self.C_rx - self.r_x, 1 + 2*self.r_x]],
                           dtype=np.float)

        if self.verbose:
            print('test_1d_zero_bc_matrix:')
            io_utils.pretty_print_matrix(A, name='A')
            io_utils.pretty_print_matrix(A_check, name='A_check')
        self.assertTrue(np.allclose(A, A_check), msg='The transport matrix for the 1D case differs from expected by greater than the tolerance')
    """

    def test_2d_zero_bc_matrix(self):
        """
        A = be.construct_transport_equation_2nd_order_centered_space(**self.model_2d_settings)

        # Construct the check matrix. In a 5x5 domain, there will be three zeros between the j-1 and i-1 terms and the
        # i+1 and j+1 terms.
        A_check = construct_2d_5x5_matrix(self.C_rx, self.C_ry, self.r_x, self.r_y)

        test_result = np.allclose(A, A_check)

        if self.verbose:
            print('test_2d_zero_bc_matrix:')
            io_utils.pretty_print_matrix(A, name='A      ')
            io_utils.pretty_print_matrix(A_check, name='A_check')
            if not test_result:
                io_utils.pretty_print_matrix(A - A_check, name='A-A_chk')

        self.assertTrue(test_result, msg='The transport matrix for the 2D case differs from expected by greater than the tolerance')
        """

        A_stencil = cn.construct_transport_matrix_with_stencil(**self.model_2d_settings)

        #test_result = np.allclose(A_stencil, A_check)
        if self.verbose:
            print('test_2d_zero_bc_matrix from stencil:')
            io_utils.pretty_print_matrix(A_stencil, name='A_stencil')
            #io_utils.pretty_print_matrix(A_check, name='A_check  ')
            #if not test_result:
            #    io_utils.pretty_print_matrix(A_stencil - A_check, name='A_st-A_ch')

        #self.assertTrue(test_result, msg='The transport matrix contructed with stencils for the 2D case differs from expected by greater than the tolerance')


if __name__ == '__main__':
    unittest.main()