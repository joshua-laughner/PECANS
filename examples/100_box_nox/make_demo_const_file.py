import netCDF4 as ncdf
import numpy as np


with ncdf.Dataset('100_box_nox_const.nc', 'w') as ds:
    ds.createDimension('x', 100)
    ho = ds.createVariable('HO', float, ('x',))
    ho[:] = np.linspace(0, 1e6, 100)
    ho.units = 'molec.cm^-3'

    o3 = ds.createVariable('O3', float, ('x',))
    o3[:] = np.linspace(0, 1e12, 100)
    o3.units = 'molec.cm^-3'
