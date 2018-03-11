import numpy as np
import sys
import unittest

from .. import backwards_euler as be
from ...utilities import io_utils

import pdb

def construct_2d_5x5_matrix(C_rx, C_ry, r_x, r_y, n=5):
    """
    Constructs the expected A matrix for a 5x5 domain given the values of the Courant and r numbers
    :param C_rx: Courant number in the x direction
    :param C_ry: Courant number in the y direction
    :param r_x: r number in the x direction
    :param r_y: y number in the y direction
    :return: a numpy array
    """

    # The form of the A matrix for a 2D domain can be thought of as having two main pieces. First, individual 1D
    # matrices representing each rows x-dimension transport populate the main block diagonal:
    #
    #            ---                                                                                               ---
    #   A_{1D} = | (1 + 2r_x + 2r_y)   (C_rx - r_x)        0                   0                   0                 |
    #            | (-C_rx - r_x)       (1 + 2r_x + 2r_y)   (C_rx - r_x)        0                   0                 |
    #            | 0                   (-C_rx - r_x)       (1 + 2r_x + 2r_y)   (C_rx - r_x)        0                 |
    #            | 0                   0                   (-C_rx - r_x)       (1 + 2r_x + 2r_y)   (C_rx - r_x)      |
    #            | 0                   0                   0                   (-C_rx - r_x)       (1 + 2r_x + 2r_y) |
    #            ---                                                                                               ---
    #
    # The first and last lines are missing a (-C_rx - r_x) and (C_rx - r_x) term, respectively, because there is no i-1
    # point in the first row and no i+1 point in the last row.
    #
    # Left this is a submatrix representing transport from the j-1 row of grid cells:
    #
    #             ---                                                                           ---
    #   A_{j-1} = | (-C_ry - r_y)   0               0               0               0             |
    #             | 0               (-C_ry - r_y)   0               0               0             |
    #             | 0               0               (-C_ry - r_y)   0               0             |
    #             | 0               0               0               (-C_ry - r_y)   0             |
    #             | 0               0               0               0               (-C_ry - r_y) |
    #             ---                                                                           ---
    #
    # And above this is a submatrix representing transport from the j+1 row of grid cells:
    #
    #             ---                                                                           ---
    #   A_{j+1} = | (C_ry - r_y)    0               0               0               0            |
    #             | 0               (C_ry - r_y)    0               0               0            |
    #             | 0               0               (C_ry - r_y)    0               0            |
    #             | 0               0               0               (C_ry - r_y)    0            |
    #             | 0               0               0               0               (C_ry - r_y) |
    #             ---                                                                           ---
    #
    # So the overall matrix is:
    #
    #       ---                                                   ---
    #   A = | A_{1D}    A_{j+1}     0           0           0       |
    #       | A_{j-1}   A_{1D}      A_{j+1}     0           0       |
    #       | 0         A_{j-1}     A_{1D}      A_{j+1}     0       |
    #       | 0         0           A_{j-1}     A_{1D}      A_{j+1} |
    #       | 0         0           0           A_{j-1}     A_{1D}  |
    #       ---                                                   ---
    #
    # where each term represents a 5x5 submatrix, either one defined above or one of all zeros. Note that in this super-
    # matrix, the first and last rows are also missing a j-1 and j+1 submatrix, respectively, just as the 1D matrix was
    # missing an i-1 and i+1 term. This again reflects the fact that the first five rows (i.e. the first row of
    # submatrices) is solving for the next timestep's concentration of the bottom (j=1) row, and so there is no row of
    # grid cells at j=0, and likewise at the j=5 row.
    #
    # To construct A, we will construct each relevant submatrix first, then insert them in the larger matrix.

    # The main, upper, and lower diagonal terms for the "1D" submatrix
    main_diag = 1 + 2*r_x + 2*r_y
    upper_diag = C_rx - r_x
    lower_diag = -C_rx - r_x
    A_1D = np.diag(main_diag * np.ones((n,))) + np.diag(upper_diag * np.ones((n-1,)), k=1) + np.diag(lower_diag * np.ones((n-1,)), k=-1)

    # The lower (j-1) and upper (j+1) block diagonal matrices
    jm1_diag = -C_ry - r_y
    A_lower = np.diag(jm1_diag * np.ones((n,)))
    jp1_diag = C_ry - r_y
    A_upper = np.diag(jp1_diag * np.ones((n,)))

    A = np.zeros((n**2, n**2), dtype=np.float, order='F')
    A[0:5, 0:5] = A_1D
    A[5:10, 5:10] = A_1D
    A[10:15, 10:15] = A_1D
    A[15:20, 15:20] = A_1D
    A[20:25, 20:25] = A_1D

    A[0:5, 5:10] = A_upper
    A[5:10, 10:15] = A_upper
    A[10:15, 15:20] = A_upper
    A[15:20, 20:25] = A_upper

    A[5:10, 0:5] = A_lower
    A[10:15, 5:10] = A_lower
    A[15:20, 10:15] = A_lower
    A[20:25, 15:20] = A_lower

    return A


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

    def test_2d_zero_bc_matrix(self):
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


if __name__ == '__main__':
    unittest.main()
