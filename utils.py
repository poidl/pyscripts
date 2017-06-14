# pylint: disable=C0103

"""Useful functions"""

import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches
import numpy as np
from scipy import interpolate
import subprocess

import regrid as regrid

# mean radius of the earth
R = 6371008.8


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


def rotate_and_interpolate(x, y, z, phi, centerx, centery):
    """Rotate grids x an y by angle phi around center (cx,cy), interpolate z onto rotated grid.

    Parameters
    ----------
    x: 1-d array
    y: 1-d array
    z: Masked data on a rectangular mesh defined by x and y
    phi: angle of rotation
    centerx: center x-coord
    centery: center y-coord

    TODO: take account of varying lat/lon spacing
    """

    ny, nx = z.shape

    dx = np.diff(x)
    dy = np.diff(y)

    if type(z) is not np.ma.MaskedArray:
        raise Exception(
            "Aborted. Array is not masked, fill values may affect interpolation.")

    test_grid_uniform(x, y)

    xx, yy = np.meshgrid(x, y)

    # translate, rotate, translate back
    xxrot, yyrot = rotate(xx, yy, phi, centerx, centery)

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
            i = x <= xp
            i = int(sum(i)) - 1
            j = y <= yp
            j = int(sum(j)) - 1

            # check if rotated point is outside of original domain
            if (i == nx - 1) | (j == ny - 1) | (i < 0) | (j < 0):
                continue

            # interpolate
            f = interpolate.interp2d(
                [y[j], y[j + 1]], [x[i], x[i + 1]], z[j:j + 2, i:i + 2])
            zrot[jrot, irot] = f(yp, xp)

    return (xxrot, yyrot, zrot)


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


def get_pmpn(x, y, phi, centerx, centery):
    """Inverse of differential distances in XI at RHO-points, after optional rotation by angle phi.

    Parameters
    ----------
    x: 1-d array
    y: 1-d array
    phi: angle of rotation
    centerx: center x-coord
    centery: center y-coord
    """

    test_grid_uniform(x, y)

    # Must convert to double, differential lats/lons are too small for using
    # f32
    dx = np.double(np.diff(x)[0])
    cs = dx * np.ones(len(x))
    cs = np.cumsum(cs) - dx
    x = np.double(x[0]) + cs
    dy = np.double(np.diff(y)[0])
    cs = dy * np.ones(len(y))
    cs = np.cumsum(cs) - dy
    y = np.double(y[0]) + cs

    nx = len(x)
    ny = len(y)

    lx = regrid.envelope_1d(x)
    ly = regrid.envelope_1d(y)

    ddx1 = (slice(None), slice(None, -1))
    ddx2 = (slice(None), slice(1, None))
    ddy1 = (slice(None, -1), slice(None))
    ddy2 = (slice(1, None), slice(None),)

    # pm
    xx, yy = np.meshgrid(lx, y)
    xrot, yrot = rotate(xx, yy, phi, centerx, centery)

    i1 = ddx1
    i2 = ddx2
    dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
        i1].flatten(), xrot[i2].flatten())
    dist = np.reshape(dist, (ny, nx))

    pm = 1 / dist

    # pn
    xx, yy = np.meshgrid(x, ly)
    xrot, yrot = rotate(xx, yy, phi, centerx, centery)

    i1 = ddy1
    i2 = ddy2
    dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
        i1].flatten(), xrot[i2].flatten())
    dist = np.reshape(dist, (ny, nx))

    pn = 1 / dist

    # dmde: ETA-derivative of inverse metric factor pm, d(1/pm)/d(ETA).
    xx, yy = np.meshgrid(lx, ly)
    xrot, yrot = rotate(xx, yy, phi, centerx, centery)

    i1 = ddx1
    i2 = ddx2
    dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
        i1].flatten(), xrot[i2].flatten())
    dist = np.reshape(dist, (ny + 1, nx))

    i1 = ddy1
    i2 = ddy2
    dmde = dist[i2] - dist[i1]

    # dndx: XI-derivative of inverse metric factor pn, d(1/pn)/d(XI).
    i1 = ddy1
    i2 = ddy2
    dist = spherical_distance(yrot[i1].flatten(), yrot[i2].flatten(), xrot[
        i1].flatten(), xrot[i2].flatten())
    dist = np.reshape(dist, (ny, nx + 1))

    i1 = ddx1
    i2 = ddx2
    dndx = dist[i2] - dist[i1]

    return (pm, pn, dmde, dndx)


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


def projection_mercator(lon, lat):
    """Mercator projection (forward)

    Parameters
    ----------
    lon: 1-d array of longitudes (degrees)
    lat: 1-d array of latitudes (degrees)
    """

    assert ((lat.ndim == 1) & (lon.ndim == 1)
            ), "Aborted. Input arrays must be one-dimensional."
    lon = np.radians(lon)
    lat = np.radians(lat)
    x = R * lon
    y = R * np.log(np.tan(np.pi / 4 + lat / 2))
    return x, y


def projection_omerc(lon, lat, lon1, lat1, lon2, lat2, beta, k0=1):
    """
    Oblique mercator projection of a sphere (forward). Optionally secant.

    Parameters
    ----------
    lon: 1-d array of longitudes (degrees)
    lat: 1-d array of latitudes (degrees)
    lon1: Lon. of 1st point defining central line
    lat1: Lat. of 1st point defining central line
    lon2: Lon. of 2nd point defining central line
    lat2: Lon. of 2nd point defining central line
    k0=1: Scale factor for secant projections
    """
    # The oblique Mercator projection is cylindrical (oblique) and
    # conformal. It has constant scale across an
    # oblique great circle (or across the two bounding circles of a
    # "secant"). This is handy for regional and coast modelling
    # if the coastline is oriented neither zonally nor meridionally.

    # Snyder, John Parr. Map projections--A working manual. Vol. 1395. US
    # Government Printing Office, 1987.
    # See their page 69 for the spherical formulas used here.

    assert ((lat.ndim == 1) & (lon.ndim == 1)
            ), "Aborted. Input arrays must be one-dimensional."
    lon = np.radians(lon)
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    lat = np.radians(lat)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    # rename variables to correspond to Snyder's notation
    phi = lat
    phi1 = lat1
    phi2 = lat2
    lam = lon
    lam1 = lon1
    lam2 = lon2

    # Compute pole of the oblique transformation

    # # Following Snyder, use arctan2 here ...
    # lamp = np.arctan2(
    #     (
    #         + np.cos(phi1) * np.sin(phi2) * np.cos(lam1)
    #         - np.sin(phi1) * np.cos(phi2) * np.cos(lam2)
    #     ), (
    #         + np.sin(phi1) * np.cos(phi2) * np.sin(lam2)
    #         - np.cos(phi1) * np.sin(phi2) * np.sin(lam1)
    #     )
    # )
    # # ... but use arctan (not arctan2) here
    # phip = np.arctan(
    #     - np.cos(lamp - lam1) / np.tan(phi1)
    # )

    phip = np.arcsin(
        np.cos(phi1) * np.sin(beta)
    )
    lamp = np.arctan(
        -np.cos(beta) / (
            - np.sin(phi1) * np.sin(beta)
        )
    ) + lam1

    # assert phip > 0, "phip not lager than zero"

    # use positive (northern) value for the pole
    (phip, lamp) = (phip, lamp) if phip >= 0 else (-phip, lamp + np.pi)
    # phip = phip if phip >= 0 else -phip

    # Calculate planar coordinates for point (phi, lam). Note that the
    # X-axis lies along the central line, x increasing easterly

    # origin of planar coordinates
    # phi0 = 0
    lam0 = lamp + np.pi / 2
    lam_, phi_ = np.meshgrid(lam, phi)
    # phi_, lam_ = np.meshgrid(phi, lam)
    A = (
        + np.sin(phip) * np.sin(phi_)
        - np.cos(phip) * np.cos(phi_) * np.sin(lam_ - lam0)
    )
    x = R * k0 * np.arctan2(
        (
            np.tan(phi_) * np.cos(phip) + np.sin(phip) * np.sin(lam_ - lam0)
        ), (
            np.cos(lam_ - lam0)
        )
    )
    y = (R / 2) * k0 * np.log(
        (1 + A) / (1 - A)
    )
    # return location of poles in degrees
    latp = np.degrees(phip)
    lonp = np.degrees(lamp)
    # return coordinates of the origin of planar coordinate system
    lato = 0
    lono = np.mod(np.degrees(lam0), 360)
    return x, y, lonp, latp, lono, lato


def projection_omerc_proj4_v1(lon, lat, lon1, lat1, alpha, k0=1):
    """
    Forward blique mercator projection of a spheroid. Optionally secant.

    Parameters
    ----------
    lon: 1-d array of longitudes (degrees)
    lat: 1-d array of latitudes (degrees)
    lon1: Lon. of point defining central line
    lat1: Lat. of point defining central line
    alpha: Azimuth of the central line in degrees clockwise from north
    k0=1: Scale factor for secant projections
    """

    if not np.isscalar(lon):
        assert ((lat.ndim == 1) & (lon.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."

    proj_cmd = (
        'proj '
        '+proj=omerc '
        '+alpha={0:.5f} '
        '+lat_0={1:.5f} +lonc={2:.5f} '
        'proj_input > proj_output'
    ).format(alpha, lat1, lon1)

    return call_proj(lon, lat, proj_cmd)


def projection_omerc_proj4_v1_inv(x, y, lon1, lat1, alpha, k0=1):
    """
    Inverse oblique mercator projection of a spheroid. Optionally secant.

    Parameters
    ----------
    x: 1-d array of longitudes (degrees)
    y: 1-d array of latitudes (degrees)
    lon1: Lon. of point defining central line
    lat1: Lat. of point defining central line
    alpha: Azimuth of the central line in degrees clockwise from north
    k0=1: Scale factor for secant projections
    """

    if not np.isscalar(x):
        assert ((x.ndim == 1) & (y.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."

    proj_cmd = (
        'proj -I -f \'% .4f\' '
        '+proj=omerc '
        '+alpha={0:.5f} '
        '+lat_0={1:.5f} +lonc={2:.5f} '
        'proj_input > proj_output'
    ).format(alpha, lat1, lon1)

    return call_proj(x, y, proj_cmd)


def projection_omerc_proj4_v2(lon, lat, lon1, lat1, lon2, lat2, k0=1):
    """
    Oblique mercator projection of a sphere (forward). Optionally secant.

    Parameters
    ----------
    lon: 1-d array of longitudes (degrees)
    lat: 1-d array of latitudes (degrees)
    lon1: Lon. of 1st point defining central line
    lat1: Lat. of 1st point defining central line
    lon2: Lon. of 2nd point defining central line
    lat2: Lon. of 2nd point defining central line
    k0=1: Scale factor for secant projections
    """

    if not np.isscalar(lon):
        assert ((lat.ndim == 1) & (lon.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."

    proj_cmd = (
        'proj '
        '+proj=omerc '
        '+lon_1={0:.5f} +lat_1={1:.5f} +lon_2={2:.5f} +lat_2={3:.5f} '
        'proj_input > proj_output'
    ).format(lon1, lat1, lon2, lat2)

    return call_proj(lon, lat, proj_cmd)


def call_proj(x, y, cmd):
    """System call to proj"""

    f = open('./proj_input', 'w')
    if np.isscalar(x):
        f.write('{0:.2f} {1:.2f}\n'.format(x, y))
    else:
        for i, _ in enumerate(x):
            f.write('{0:.2f} {1:.2f}\n'.format(x[i], y[i]))

    f.close()
    print('Running \'' + cmd + '\'')
    subprocess.run(
        cmd,
        shell=True,
        check=True
    )

    res = np.loadtxt('proj_output')

    if np.isscalar(x):
        xout = res[0]
        yout = res[1]
    else:
        xout = res[:, 0]
        yout = res[:, 1]

    return xout, yout


def get_dist_proj(lon1, lat1, lon2, lat2):
    """Geodesic computations with geod (proj4)

    Parameters
    ----------
    lon1: 1-d array
    lat1: 1-d array
    lon2: 1-d array
    lat2: 1-d array
    """

    assert ((lat1.ndim == 1) &
            (lon1.ndim == 1) &
            (lat2.ndim == 1) &
            (lon2.ndim == 1)
            ), "Aborted. Input arrays must be one-dimensional."

    geod_cmd = (
        'geod -I -f \'% .4f\' '
        '+ellps=WGS84 '
        'geod_input > geod_output'
    )

    return call_geod(lon1, lat1, lon2, lat2, geod_cmd)


def call_geod(x1, y1, x2, y2, cmd):
    """System call to geod"""

    f = open('./geod_input', 'w')
    if np.isscalar(x1):
        f.write('{0:.2f} {1:.2f} {2:.2f} {3:.2f}\n'.format(x1, y1, x2, y2))
    else:
        for i, _ in enumerate(x1):
            f.write('{0:.2f} {1:.2f} {2:.2f} {3:.2f}\n'.format(
                y1[i], x1[i], y2[i], x2[i]))

    f.close()
    print('Running \'' + cmd + '\'')
    subprocess.run(
        cmd,
        shell=True,
        check=True
    )

    res = np.loadtxt('geod_output')

    if np.isscalar(x1):
        xout = res[0]
        yout = res[1]
    else:
        xout = res[:, 0]
        yout = res[:, 1]

    return xout, yout
