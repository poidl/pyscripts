# pylint: disable=C0103

"""Useful functions"""

import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches
import numpy as np
from scipy import interpolate


import regrid as regrid


def rotate_origin(x, y, phi):
    """Rotate grids x an y by angle phi around origin of coordinate axis.

    Parameters
    ----------
    x: 2-d array
    y: 2-d array
    phi: angle in radians
    """
    xrot = x * np.cos(phi) - y * np.sin(phi)
    yrot = x * np.sin(phi) + y * np.cos(phi)
    return (xrot, yrot)


def rotate(x, y, phi, centerx, centery):
    """Rotate grids x an y by angle phi around (centerx, centery).

    Parameters
    ----------
    x: 2-d array
    y: 2-d array
    phi: angle in radians
    centerx: center x-coord
    centery: center y-coord
    """
    # translate, rotate, translate back
    x = x - centerx
    y = y - centery
    xrot, yrot = rotate_origin(x, y, phi)
    xrot = xrot + centerx
    yrot = yrot + centery

    return (xrot, yrot)

# DO NOT USE THIS: INSTEAD, INTERPOLATE ONTO POINTS OBTAINED FROM INVERSE
# PROJECTION
# def rotate_and_interpolate(x, y, z, phi, centerx, centery):
#     """Rotate grids x an y by angle phi around center (cx,cy), interpolate z onto rotated grid.

#     Parameters
#     ----------
#     x: 1-d array
#     y: 1-d array
#     z: Masked data on a rectangular mesh defined by x and y
#     phi: angle of rotation
#     centerx: center x-coord
#     centery: center y-coord

#     TODO: take account of varying lat/lon spacing
#     """

#     ny, nx = z.shape

#     dx = np.diff(x)
#     dy = np.diff(y)

#     if type(z) is not np.ma.MaskedArray:
#         raise Exception(
#             "Aborted. Array is not masked, fill values may affect interpolation.")

#     test_grid_uniform(x, y)

#     xx, yy = np.meshgrid(x, y)

#     # translate, rotate, translate back
#     xxrot, yyrot = rotate(xx, yy, phi, centerx, centery)

#     # output is initially masked everywhere
#     o = np.nan * np.ones((ny, nx))
#     zrot = np.ma.array(o, mask=np.ones((ny, nx)))

#     # suppress warning
#     zrot.unshare_mask()

#     # TODO: do this more efficiently, but in a way that accounts for masked
#     # values. Here, each rotated data point is considered individually, with
#     # its 4 surrounding gridpoints on the original grid.

#     for jrot in np.arange(ny):
#         for irot in np.arange(nx):
#             xp = xxrot[jrot, irot]
#             yp = yyrot[jrot, irot]
#             i = x <= xp
#             i = int(sum(i)) - 1
#             j = y <= yp
#             j = int(sum(j)) - 1

#             # check if rotated point is outside of original domain
#             if (i == nx - 1) | (j == ny - 1) | (i < 0) | (j < 0):
#                 continue

#             # interpolate
#             f = interpolate.interp2d(
#                 [y[j], y[j + 1]], [x[i], x[i + 1]], z[j:j + 2, i:i + 2])
#             zrot[jrot, irot] = f(yp, xp)

#     return (xxrot, yyrot, zrot)


def cut(rectangle, z):
    """Cut domain to rectangular shape (e.g. after rotating)

    Parameters
    ----------
    rectangle: tuple ((y,x), height, width) in units of indices
    z: data on a rectangular mesh
    """

    mask = np.ones(z.shape, dtype=bool)
    y1 = rectangle[0][0]
    x1 = rectangle[0][1]
    y2 = y1 + rectangle[1]
    x2 = x1 + rectangle[2]
    mask[y1:y2, x1:x2] = 0

    # alert if masked parts (e.g. fill values or nans) are not cut
    if type(z) is np.ma.MaskedArray:
        outside = ~mask & z.mask
        if np.any(outside):
            raise Exception(
                'Rotated domain outside of original domain (total of ' + str(
                    sum(outside.flatten())) + ' gridpoints).')
    z = z[~mask]
    ii = np.where(~mask)
    ny = ii[0][-1] - ii[0][0] + 1
    nx = ii[1][-1] - ii[1][0] + 1
    z = np.reshape(z, (ny, nx))
    return z


def spherical_distance(lat1, lat2, lon1, lon2):
    """Distance between the two points along a great circle of the sphere.

    Parameters
    ----------
    lat1: 2-d array of latitudes (in degrees) of point 1
    lat2: latitudes of point 2 (degrees)
    lon1: longitudes of point 1 (degrees)
    lon2: longitudes of point 2 (degrees)
    """
    # from: https://rosettacode.org/wiki/Haversine_formula#Python

    dLat = np.radians(lat2 - lat1)
    dLon = np.radians(lon2 - lon1)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    a = np.sin(dLat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dLon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return R * c

# DO NOT USE THIS
# def get_pmpn(x, y, phi, centerx, centery):
#     """Inverse of differential distances in XI at RHO-points, after optional rotation by angle phi.

#     Parameters
#     ----------
#     x: 1-d array
#     y: 1-d array
#     phi: angle of rotation
#     centerx: center x-coord
#     centery: center y-coord
#     """

#     test_grid_uniform(x, y)

#     # Must convert to double, differential lats/lons are too small for using
#     # f32
#     dx = np.double(np.diff(x)[0])
#     cs = dx * np.ones(len(x))
#     cs = np.cumsum(cs) - dx
#     x = np.double(x[0]) + cs
#     dy = np.double(np.diff(y)[0])
#     cs = dy * np.ones(len(y))
#     cs = np.cumsum(cs) - dy
#     y = np.double(y[0]) + cs

#     nx = len(x)
#     ny = len(y)

#     lx = regrid.envelope_1d(x)
#     ly = regrid.envelope_1d(y)

#     ddx1 = (slice(None), slice(None, -1))
#     ddx2 = (slice(None), slice(1, None))
#     ddy1 = (slice(None, -1), slice(None))
#     ddy2 = (slice(1, None), slice(None),)

#     # pm
#     xx, yy = np.meshgrid(lx, y)
#     xrot, yrot = rotate(xx, yy, phi, centerx, centery)

#     i1 = ddx1
#     i2 = ddx2
#     dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
#         i1].flatten(), xrot[i2].flatten())
#     dist = np.reshape(dist, (ny, nx))

#     pm = 1 / dist

#     # pn
#     xx, yy = np.meshgrid(x, ly)
#     xrot, yrot = rotate(xx, yy, phi, centerx, centery)

#     i1 = ddy1
#     i2 = ddy2
#     dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
#         i1].flatten(), xrot[i2].flatten())
#     dist = np.reshape(dist, (ny, nx))

#     pn = 1 / dist

#     # dmde: ETA-derivative of inverse metric factor pm, d(1/pm)/d(ETA).
#     xx, yy = np.meshgrid(lx, ly)
#     xrot, yrot = rotate(xx, yy, phi, centerx, centery)

#     i1 = ddx1
#     i2 = ddx2
#     dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
#         i1].flatten(), xrot[i2].flatten())
#     dist = np.reshape(dist, (ny + 1, nx))

#     i1 = ddy1
#     i2 = ddy2
#     dmde = dist[i2] - dist[i1]

#     # dndx: XI-derivative of inverse metric factor pn, d(1/pn)/d(XI).
#     i1 = ddy1
#     i2 = ddy2
#     dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
#         i1].flatten(), xrot[i2].flatten())
#     dist = np.reshape(dist, (ny, nx + 1))

#     i1 = ddx1
#     i2 = ddx2
#     dndx = dist[i2] - dist[i1]

#     return (pm, pn, dmde, dndx)


def test_grid_uniform(x, y):
    if (x.ndim != 1) | (y.ndim != 1):
        raise Exception(
            "Aborted. Input arrays must be one-dimensional")

    dx = np.diff(x)
    dy = np.diff(y)
    # testing for equidistant lat/lon grid fails below if the increments are
    # not rounded (float32)
    dx0 = np.around(dx[0], decimals=5)
    dx1 = np.around(dx[-1], decimals=5)
    dy0 = np.around(dy[0], decimals=5)
    dy1 = np.around(dy[-1], decimals=5)
    if (dx0 != dx1) | (dy0 != dy1) | (dx0 != dy0) | (dx1 != dy1):
        raise Exception(
            "Aborted. Not smart enough to handle non-uniform grids.")


def myplot2d(z, name):
    fig = plt.figure(figsize=(8, 8))
    p1 = plt.imshow(z, origin='lower', interpolation='nearest')
    plt.colorbar()
    fig.savefig('figures/' + name + '.pdf')
    print('done')


def myplot_add_path(vertices, axes):
    p = path.Path(vertices)
    patch = patches.PathPatch(p, facecolor='none', lw=1)
    axes.add_patch(patch)


def hanning_smoother(h):
    """Copied from Romstools. TODO: find out whether this is Hann or Hamming?"""
    [M, L] = np.array(h.shape)
    L = L - 1
    M = M - 1
    Mm = M - 1
    Mmm = M - 2
    Lm = L - 1
    Lmm = L - 2

    h[1:Mm, 1:Lm] = 0.125 * (h[0:Mmm, 1:Lm] + h[2:M, 1:Lm] +
                             h[1:Mm, 0:Lmm] + h[1:Mm, 2:L] +
                             4 * h[1:Mm, 1:Lm])

    h[0, :] = h[1, :]
    h[M, :] = h[Mm, :]
    h[:, 0] = h[:, 1]
    h[:, L] = h[:, Lm]
    return h
