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


def smooth_hanning(h):

    h[1:-1, 1:-1] = 0.125 * (h[0:-2, 1:-1] + h[2:, 1:-1] +
                             h[1:-1, 0:-2] + h[1:-1, 2:] +
                             4 * h[1:-1, 1:-1])

    h[0, :] = h[1, :]
    h[-1, :] = h[-2, :]
    h[1:-1, 0] = h[1:-1, 1]
    h[1:-1, -1] = h[1:-1, -1]
    return h


def rxold(h1, h2):
    return (h1 - h2) / (h1 + h2)


def iter(h):
    for j in range(h.shape[0]):
        for i in range(h.shape[1] - 1):
            r = rxold(h[j, i], h[j, i + 1])
            rt = 0.2
            if r > 0.2:
                # rt = r
                h[j, i + 1] = ((1 - rt) / (1 + rt)) * h[j, i]
    return h


def rx(h2, h1):
    """Slope parameter, Beckmann Haidvogel (1993)"""
    # but note that depth variation is of opposite sign than bathymetry slope
    # rx is positive for downward (negative) slope
    return (h2 - h1) / (h2 + h1)


def lower_topo(h):
    for i in range(h.shape[1] - 1):
        r = rx(h[:, i + 1], h[:, i])
        rt = -0.2
        # r is negative if slope is positive
        j = r < - 0.2
        h[j, i + 1] = ((1 + rt) / (1 - rt)) * h[j, i]
    return h


def raise_topo(h):
    for i in range(h.shape[1] - 1):
        r = rx(h[:, i + 1], h[:, i])
        rt = 0.2
        # r is positive if slope is negative
        j = r > 0.2
        h[j, i + 1] = ((1 + rt) / (1 - rt)) * h[j, i]
    return h


def get_rx(h):
    """Magnitude of Beckman & Haidvogel number

    Parameters
    ----------
    h: topography, 2-d array
    """
    assert h.ndim == 2, "Aborted. Input array must be two-dimensional."

    hxm = h[:, 1:] - h[:, :-1]
    hxp = h[:, 1:] + h[:, :-1]
    hym = h[1:, :] - h[:-1, :]
    hyp = h[1:, :] + h[:-1, :]
    rx0x = np.abs(hxm) / hxp
    rx0y = np.abs(hym) / hyp
    return rx0x, rx0y


def smooth_martinho(h):
    """Bathymetry smoother, loosely following Martinho and Batteen (2006).
    This version is without volume conservation. I don't have access to
    Ocean Modelling, so I can't download the supplementary data in which
    they seem to describe their algorithm in detail.

    Parameters
    ----------
    h: topography, 2-d array
    """
    heast = lower_topo(h.copy())
    heast = iter(h.copy())
    hnorth = lower_topo(heast.copy().transpose()[::-1, :])
    hwest = lower_topo(hnorth.copy().transpose()[::-1, :])
    hsouth = lower_topo(hwest.copy().transpose()[::-1, :])
    hnew = hsouth.copy().transpose()[::-1, :]
    return hnew
