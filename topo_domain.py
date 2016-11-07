#!/bin/python
# pylint: disable=C0103

from netCDF4 import Dataset as nc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fname = '../../../work/data/ETOPO2v2g_f4.nc'
grdname = 'topo.nc'
FIGNAME = 'figure.pdf'

# domain boundaries
# Syd:
# lonmin=150
# lonmax=152
# latmin=-35
# latmax=-32
# EAC:
lonmin = 142
lonmax = 164
latmin = -39.2
latmax = -15

# subsampling?
skipx = 5
skipy = 5

ff = nc(fname, 'r')
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

bx = np.logical_and(bx, bxs)
by = np.logical_and(by, bys)

nxx = np.sum(bx)
nyy = np.sum(by)

z = ff.variables['z'][by, bx]
# ff.close()

z = np.reshape(z, (nyy, nxx))
z.mask[z > 0] = True

# topo_smoothing(z)


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
