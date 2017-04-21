# pylint: disable=C0103

"""Useful functions"""

import numpy as np
from scipy import interpolate


def rotate_meshgrid(x, y, phi):
    """Rotate grids x an y by angle phi.

    Parameters
    ----------
    x: 2-d array
    y: 2-d array
    """
    xrot = x * np.cos(phi) - y * np.sin(phi)
    yrot = x * np.sin(phi) + y * np.cos(phi)
    return (xrot, yrot)


def rotate(x, y, z, phi):
    """Rotate grids x an y by angle phi, interpolate z onto rotated grid.

    Parameters
    ----------
    x: 1-d array
    y: 1-d array
    z: data on a rectangular mesh defined by x and y
    phi: angle of rotation
    """
    xx, yy = np.meshgrid(x, y)
    fz = interpolate.RectBivariateSpline(y, x, z)

    # center of rotation
    centerx = x[int(len(x) / 2)]
    centery = y[int(len(y) / 2)]

    # translate, rotate, translate back
    xx = xx - centerx
    yy = yy - centery
    xxrot, yyrot = rotate_meshgrid(xx, yy, phi)
    xxrot = xxrot + centerx
    yyrot = yyrot + centery

    # evaluate interpolation
    zrot = fz.ev(yyrot.flatten(), xxrot.flatten())
    nx = len(x)
    ny = len(y)
    zrot = np.reshape(zrot, (ny, nx))

    # mask
    if type(z) is np.ma.MaskedArray:
        zrot[zrot > 0] = np.nan
    xout = np.logical_or(xxrot < x.min(), xxrot > x.max())
    yout = np.logical_or(yyrot < y.min(), yyrot > y.max())
    mask = xout | yout
    zrot[mask] = np.nan

    # return (xxrot, yyrot, zrot)
    return (zrot, mask)


def cut(limits, z, mask_rotation):
    """Cut domain to rectangular shape (e.g. after rotating)

    Parameters
    ----------
    limits: boundary indices of the domain (included) 
    z: data on a rectangular mesh
    mask_rotation: true if points were left undefined after rotating
    """
    mask = np.ones(z.shape, dtype=bool)
    y1 = limits[0]
    y2 = limits[1]
    x1 = limits[2]
    x2 = limits[3]
    mask[y1:y2, x1:x2] = 0

    # alert if parts outside of the original domain are not cut
    outside = ~mask & mask_rotation
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
