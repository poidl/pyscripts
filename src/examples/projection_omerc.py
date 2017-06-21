#!/bin/python
# pylint: disable=C0103

"""
Outlines how to make a rotated grid using an oblique Mercator
projection by
1) defining a rectangular grid on a projection plane
2) inverse-projecting this grid onto a sphereoid
3) measuring the geodesic distances on the spheroid to obtain the
  grid-point spacing (and hence pm and pn)

Essentially this follows Kate's advice in the ROMS forum

https://www.myroms.org/forum/viewtopic.php?f=14&t=4540&sid=921cb348e0ba86f4a9c79c4d92645ea7

except for the use of the oblique Mercator projection. Using proj4,
it should be straightforward to adapt this to any projection.
"""

from netCDF4 import Dataset as nc
import utils as u
import utils_nc as unc
import numpy as np
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import regrid as regrid
import projections as proj

fname = '../../../../../work/data/ETOPO2v2g_f4.nc'
FIGPATH = '../../figures'

# Region covering north-west America and the North-east pacific
lonmin = -158.5
lonmax = -86
latmin = 2
latmax = 58

# subsampling?
skipx = 30
skipy = 30

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
bx = np.logical_and(x >= lonmin, x <= lonmax)
by = np.logical_and(y >= latmin, y <= latmax)

# subsample mask within the domain
bx = np.logical_and(bx, bxs)
by = np.logical_and(by, bys)

# axes for the subsampled data
lon = x[bx]
lat = y[by]
nx = np.sum(bx)
ny = np.sum(by)

# extract from nc file
z = ff.variables['z'][by, bx]
z = np.reshape(z, (ny, nx))

# rectangular grid in the projection plane
xg = 1e6 * np.arange(-1, 4.4, 0.4)
yg = 1e6 * np.arange(-1, 1.4, 0.4)
xg, yg = np.meshgrid(xg, yg)
nyg, nxg = xg.shape

# define an oblique Mercator projection by specifying a point and an
# azimuth
lonc = -131
latc = 40
angle_rot = -0.5 * np.pi / 2
k0 = 1

# currently (June 15th 2017) the +no_rot switch is broken in proj:
# https://github.com/OSGeo/proj.4/issues/523#issuecomment-306779910
# This means proj doesn't rotate the domain according to the defined
# azimuth. So we have to rotate ourselves:
xgr, ygr = u.rotate_origin(xg, yg, -angle_rot - np.pi / 2)

# inverse projection of the rotated grid
long, latg = proj.projection_omerc_v1_inv(
    xgr.flatten(), ygr.flatten(), lonc, latc, np.degrees(angle_rot), k0)

long = np.reshape(long, (nyg, nxg))
latg = np.reshape(latg, (nyg, nxg))

# get geodesic distances (in case of a sphere, these are Great Circle
# distances, but note that proj uses WGS84 per default)
long1 = long[:, :-1]
latg1 = latg[:, :-1]
long2 = long[:, 1:]
latg2 = latg[:, 1:]

dx = proj.get_dist_proj(long1.flatten(), latg1.flatten(),
                        long2.flatten(), latg2.flatten())

long1 = long[:-1, :]
latg1 = latg[:-1, :]
long2 = long[1:, :]
latg2 = latg[1:, :]

dy = proj.get_dist_proj(long1.flatten(), latg1.flatten(),
                        long2.flatten(), latg2.flatten())

dx = np.reshape(dx, (nyg, nxg - 1))
dy = np.reshape(dy, (nyg - 1, nxg))

# Sanity-check the projection. An oblique Mercator projection should
# yield constant cell sizes parallel to the central line (the line
# or segment where the projection cylinder touches the spheroid.
# So we compute the area and plot it.
# Average dx (dy) meridionally (zonally) to cell center
area = 0.25 * (
    (dx[:-1, :] + dx[1:, :]) *
    (dy[:, :-1] + dy[:, 1:])
)

# project the entire domain for plotting
lon, lat = np.meshgrid(lon, lat)

x, y = proj.projection_omerc_v1(
    lon.flatten(), lat.flatten(), lonc, latc, np.degrees(angle_rot), k0)

x = np.reshape(x, (ny, nx))
y = np.reshape(y, (ny, nx))


# regrid for a pcolormesh plot
xe = regrid.envelope(x, axis=1)
xe = regrid.envelope(xe, axis=0)
ye = regrid.envelope(y, axis=0)
ye = regrid.envelope(ye, axis=1)

# rotate the entire domain for plotting
xer, yer = u.rotate_origin(xe, ye, angle_rot + np.pi / 2)


# #####################################
# plot the system with true north pointing in y-dir at lonc, latc
fig = plt.figure(figsize=(8, 8))

width = xe.max() - xe.min()
height = ye.max() - ye.min()

pw = 0.72
pos1 = [0.12, 0.15, pw, height * (pw / width)]
ax1 = fig.add_axes(pos1)

p1 = ax1.pcolormesh(xe / 1e6, ye / 1e6, z, shading='flat')
ax1.set_aspect('equal')
ax1.set_xlabel(
    'Projected (planar) distance (1000 km).')
ax1.set_ylabel('Projected (planar) distance (1000 km).')
ax1.grid(color='k')

# sanity-check whether origin is where it's supposed to be
xc, yc = proj.projection_omerc_v1(
    lonc, latc, lonc, latc, np.degrees(angle_rot), k0)

# Draw the grid. Note that this is the non-rotated variable, since
# we rotated the entire domain before plotting
ax1.plot(xgr.flatten() / 1e6, ygr.flatten() / 1e6, 'r.')
ax1.plot(xc.flatten() / 1e6, yc.flatten() / 1e6, 'y.', markersize=15)

pos2 = [0.89, 0.1, 0.03, 0.8]
ax2 = fig.add_axes(pos2)
plt.colorbar(p1, ax2)

fig.savefig(FIGPATH + '/projection_omerc.png')
print('done')

# #####################################
# plot the system with the x-axis parallel to central line

fig = plt.figure(figsize=(8, 8))

width = xer.max() - xer.min()
height = yer.max() - yer.min()

pw = 0.72
pos1 = [0.12, 0.15, pw, height * (pw / width)]
ax1 = fig.add_axes(pos1)

p1 = ax1.pcolormesh(xer / 1e6, yer / 1e6, z, shading='flat')
ax1.set_aspect('equal')
ax1.set_xlabel(
    'Projected (planar) distance (1000 km).')
ax1.set_ylabel('Projected (planar) distance (1000 km).')
ax1.grid(color='k')

# sanity-check whether origin is where it's supposed to be
xc, yc = proj.projection_omerc_v1(
    lonc, latc, lonc, latc, np.degrees(angle_rot), k0)

# Draw the grid. Note that this is the non-rotated variable, since
# we rotated the entire domain before plotting
ax1.plot(xg.flatten() / 1e6, yg.flatten() / 1e6, 'r.')
ax1.plot(xc.flatten() / 1e6, yc.flatten() / 1e6, 'y.', markersize=15)

pos2 = [0.89, 0.1, 0.03, 0.8]
ax2 = fig.add_axes(pos2)
plt.colorbar(p1, ax2)

fig.savefig(FIGPATH + '/projection_omerc_rotated.png')
print('done')

# ###################################
# as a sanity check, plot the grid cell area and visually inspect

fig = plt.figure(figsize=(8, 8))
plt.imshow(area)
plt.colorbar()

fig.savefig(FIGPATH + '/projection_omerc_area.png')
print('done')
