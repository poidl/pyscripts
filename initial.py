#!/bin/python
# pylint: disable=C0103

import sys
from netCDF4 import Dataset as nc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils as u
import utils_nc as unc
grdname_out = 'ocean_ini.nc'
FIGNAME = 'figure.pdf'

nx = 20
ny = 64
nz = 5

zeta = np.zeros((ny, nx))
ubar = np.zeros((ny, nx))
vbar = np.zeros((ny, nx))
u = np.zeros((nz, ny, nx))
v = np.zeros((nz, ny, nx))
temp = np.zeros((nz, ny, nx))
salt = np.zeros((nz, ny, nx))
# ###################################

# write to nc file
unc.create_ini(grdname_out, nz, ny, nx)
grid = nc(grdname_out, 'a')
grid.variables['ocean_time'][:] = 0
grid.variables['zeta'][:] = zeta
grid.variables['ubar'][:] = ubar
grid.variables['vbar'][:] = vbar
grid.variables['u'][:] = u
grid.variables['v'][:] = v
grid.variables['temp'][:] = temp
grid.variables['salt'][:] = salt
grid.close()
print('done')
# ###################################

# fig = plt.figure(figsize=(8, 8))
# p1 = plt.imshow(zrot, origin='lower', interpolation='nearest')
# fig.savefig('figures/zrot.pdf')
# print('done')
