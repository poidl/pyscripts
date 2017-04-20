#!/bin/python
# pylint: disable=C0103

import sys
from netCDF4 import Dataset as nc
import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils as u
import utils_nc as unc
fname = '../../../work/data/ETOPO2v2g_f4.nc'
grdname = 'topo.nc'
grdname_out = 'ocean_grd.nc'
FIGNAME = 'figure.pdf'

# domain boundaries
# Syd:
# lonmin=150
# lonmax=152
# latmin=-35
# latmax=-32
# EAC:
# lonmin = 142
# lonmax = 164
# latmin = -39.2
# latmax = -15
# Interior Southern Californian Bight
lonmin = -121
lonmax = -117.5
latmin = 32.5
latmax = 34.7

# subsampling?
skipx = 1
skipy = 1

ff = nc(fname, 'r')

# load grid
x = ff.variables['x'][:]
y = ff.variables['y'][:]
nx = x.size
ny = y.size

# subsample mask
ixs = np.arange(0, nx, skipx)
iys = np.arange(0, ny, skipy)
bxs = np.zeros(nx, dtype=bool)
bys = np.zeros(ny, dtype=bool)
bxs[ixs] = True
bys[iys] = True

# domain mask
bx = np.logical_and(x > lonmin, x < lonmax)
by = np.logical_and(y > latmin, y < latmax)

# subsample mask within the domain
bx = np.logical_and(bx, bxs)
by = np.logical_and(by, bys)

# axes for the subsampled data
x = x[bx]
y = y[by]
nx = np.sum(bx)
ny = np.sum(by)

# extract from nc file
z = ff.variables['z'][by, bx]
z = np.reshape(z, (ny, nx))
z.mask[z > 0] = True

######################################

xx, yy = np.meshgrid(x, y)
fz = interpolate.RectBivariateSpline(y, x, z)

# center of rotation
centerx = x[int(len(x) / 2)]
centery = y[int(len(y) / 2)]

# translate, rotate, translate back
xx = xx - centerx
yy = yy - centery
xxrot, yyrot = u.rotate(xx, yy, np.pi / 3)
xxrot = xxrot + centerx
yyrot = yyrot + centery

# evaluate interpolation
zrot = fz.ev(yyrot.flatten(), xxrot.flatten())
zrot = np.reshape(zrot, (ny, nx))

# mask
zrot[zrot > 0] = np.nan
xout = np.logical_or(xxrot < x.min(), xxrot > x.max())
yout = np.logical_or(yyrot < y.min(), yyrot > y.max())
zrot[xout | yout] = np.nan

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(zrot, origin='lower', interpolation='nearest')
fig.savefig('figures/rot.pdf')
print('done')


###################################
mask = np.zeros(z.shape, dtype=bool)
mask[:, :40] = 1
mask[:, -35:] = 1

# alert if parts outside of the original domain are not cut
outside = ~mask & [xout | yout]
if np.any(outside):
    raise Exception(
        'Rotated domain outside of original domain (total of ' + str(
            sum(outside.flatten())) + ' gridpoints).')
zrot = zrot[~mask]

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(zrot, origin='lower', interpolation='nearest')
fig.savefig('figures/domain.pdf')
print('done')

###################################
fig = plt.figure(figsize=(8, 8))

pos1 = [0.1, 0.1, 0.65, 0.8]
ax1 = fig.add_axes(pos1)
p1 = ax1.imshow(z, origin='lower', interpolation='nearest')
ax1.set_aspect('auto')
ax1.set_xlabel('Gridpoints in x-direction')
ax1.set_ylabel('Gridpoints in y-direction')

pos2 = [0.88, 0.1, 0.03, 0.8]
ax2 = fig.add_axes(pos2)
plt.colorbar(p1, ax2)

ax3 = fig.add_axes(pos1, frameon=False)
ax3.set_xlim(lonmin, lonmax)
ax3.set_ylim(latmin, latmax)
ax3.grid()
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')

ax3.yaxis.set_label_position('left')
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_label_position('bottom')
ax3.xaxis.set_ticks_position('bottom')

ax1.yaxis.set_label_position('right')
ax1.yaxis.set_ticks_position('right')
ax1.xaxis.set_label_position('top')
ax1.xaxis.set_ticks_position('top')

fig.savefig('figures/' + FIGNAME)
print('done')

######################################

unc.create_grd(grdname_out, ny, nx)

grid = nc(grdname_out, 'a')
# grid.variables['xh'][:] = xh
# grid.variables['yh'][:] = yh
grid.variables['h'][:] = z
grid.variables['xl'][:] = nx
grid.variables['el'][:] = ny
grid.variables['spherical'][:] = 0
grid.close()
