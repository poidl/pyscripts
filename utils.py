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
    z: Masked data on a rectangular mesh defined by x and y
    phi: angle of rotation

    TODO: take account of varying lat/lon spacing
    """

    ny, nx = z.shape
    dx = np.diff(x)
    dy = np.diff(y)

    if type(z) is not np.ma.MaskedArray:
        raise Exception(
            "Aborted. Array is not masked, fill values may affect interpolation.")

    if (len(x) != nx) | (len(y) != ny) | (dx[0] != dx[-1]) | (dy[0] != dy[-1]) | (dx[0] != dy[0]):
        raise Exception(
            "Aborted. Not smart enough to handle irregular grids.")

    xx, yy = np.meshgrid(x, y)

    # center of rotation
    centerx = x[int(len(x) / 2)]
    centery = y[int(len(y) / 2)]

    # translate, rotate, translate back
    xx = xx - centerx
    yy = yy - centery
    xxrot, yyrot = rotate_meshgrid(xx, yy, phi)
    xxrot = xxrot + centerx
    yyrot = yyrot + centery

    # output is initially masked everywhere
    o = np.nan * np.ones((ny, nx))
    zrot = np.ma.array(o, mask=np.ones((ny, nx)))

    # suppress warning
    zrot.unshare_mask()

    # TODO: do this more efficiently, but in a way that accounts for masked
    # values. Here, each rotated data point is considered individually, with
    # its 4 surrounding gridpoints on the original grid.

    for jrot in np.arange(ny):
        for irot in np.arange(nx):
            xp = xxrot[jrot, irot]
            yp = yyrot[jrot, irot]
            i = x < xp
            i = int(sum(i)) - 1
            j = y < yp
            j = int(sum(j)) - 1

            # check if rotated point is outside of original domain
            if (i == nx - 1) | (j == ny - 1) | (i < 0) | (j < 0):
                continue

            # interpolate
            f = interpolate.interp2d(
                [y[j], y[j + 1]], [x[i], x[i + 1]], z[j:j + 2, i:i + 2])
            zrot[jrot, irot] = f(yp, xp)

    return (zrot)


def cut(limits, z):
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
