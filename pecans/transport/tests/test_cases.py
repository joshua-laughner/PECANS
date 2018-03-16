import copy
import numpy as np

dt = 60.0
dx = 1000.0
dy = 1000.0
dz = 500.0
u_x = 5.0
u_y = 0.0
u_z = 0.0
D_x = 100.0
D_y = 100.0
D_z = 500.0

model_1d_settings = {'dt': dt, 'dx': dx, 'dy': None, 'dz': None, 'u_x': u_x, 'u_y': None, 'u_z': None,
                     'D_x': D_x, 'D_y': None, 'D_z': None, 'domain_size': None, 'boundary_conditions': None}
model_2d_settings = {'dt': dt, 'dx': dx, 'dy': dy, 'dz': None, 'u_x': u_x, 'u_y': u_y, 'u_z': None,
                     'D_x': D_x, 'D_y': D_y, 'D_z': None, 'domain_size': None, 'boundary_conditions': None}


def point_1d():
    domain = np.zeros((51,), dtype=np.float)
    domain[10] = 1.0
    model_kwargs = copy.copy(model_1d_settings)
    model_kwargs['domain_size'] = domain.shape
    return domain, model_kwargs


def gaussian_1d(x0=None, sigma_x=None):
    model_kwargs = copy.copy(model_1d_settings)
    if x0 is None:
        x0 = dx * 10
    if sigma_x is None:
        sigma_x = 4*dx

    x = np.arange(51, dtype=np.float)*dx
    domain = np.exp(-0.5*((x - x0)/sigma_x)**2)
    model_kwargs['domain_size'] = domain.shape
    return domain, model_kwargs


def point_2d():
    domain = np.zeros((51, 31), dtype=np.float, order='F')
    domain[10, 15] = 1.0
    model_kwargs = copy.copy(model_2d_settings)
    model_kwargs['domain_size'] = domain.shape
    return domain, model_kwargs


def gaussian_2d(x0=None, y0=None, sigma_x=None, sigma_y=None):
    model_kwargs = copy.copy(model_2d_settings)
    if x0 is None:
        x0 = 10 * dx
    if sigma_x is None:
        sigma_x = 4 * dx
    if y0 is None:
        y0 = 15 * dy
    if sigma_y is None:
        sigma_y = 4 * dy

    y, x = np.meshgrid(np.arange(31) * dy, np.arange(51) * dx)
    domain = np.exp(-( ((x-x0)/sigma_x)**2 + ((y-y0)/sigma_y)**2 ))
    model_kwargs['domain_size'] = domain.shape
    return domain, model_kwargs
