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

rot_angle = np.pi / 3
zrot, mask_rotation = u.rotate(x, y, z, rot_angle)
limits = [0, -1, 40, 60]
zcut = u.cut(limits, zrot, mask_rotation)

# coriolis parameter
phi = 2 * np.pi * (y / 360)
omega = 7.2921159e-5
f = 2 * omega * np.cos(phi)
f = np.tile(f[::-1], [nx, 1]).transpose()

frot, mask_rotation = u.rotate(x, y, f, rot_angle)
fcut = u.cut(limits, frot, mask_rotation)
nycut, nxcut = zcut.shape

###################################

# write to nc file
unc.create_grd(grdname_out, nycut, nxcut)
grid = nc(grdname_out, 'a')
grid.variables['h'][:] = zcut
grid.variables['xl'][:] = nxcut
grid.variables['el'][:] = nycut
grid.variables['spherical'][:] = 0
grid.variables['f'][:] = fcut
grid.close()

###################################

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(zrot, origin='lower', interpolation='nearest')
fig.savefig('figures/zrot.pdf')
print('done')

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(zcut, origin='lower', interpolation='nearest')
fig.savefig('figures/zcut.pdf')
print('done')

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(f, origin='lower', interpolation='nearest')
fig.savefig('figures/f.pdf')
print('done')

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(frot, origin='lower', interpolation='nearest')
fig.savefig('figures/frot.pdf')
print('done')

fig = plt.figure(figsize=(8, 8))
p1 = plt.imshow(fcut, origin='lower', interpolation='nearest')
fig.savefig('figures/fcut.pdf')
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
