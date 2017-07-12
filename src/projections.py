"""Projections."""

import numpy as np
import subprocess
import os

# mean radius of the earth
R = 6371008.8


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


def projection_omerc_spherical(lon, lat, lon1, lat1, lon2, lat2, beta, k0=1):
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


def projection_omerc_v1(lon, lat, lon1, lat1, alpha, k0=1):
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
        '+alpha={0:.8f} '
        '+lat_0={1:.14f} +lonc={2:.14f} +k={3:.14f} '
        'proj_input > proj_output'
    ).format(alpha, lat1, lon1, k0)

    return call_proj(lon, lat, proj_cmd)


def projection_omerc_v1_inv(x, y, lon1, lat1, alpha, k0=1):
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
        'proj -I -f \'% .30f\' '
        '+proj=omerc '
        '+alpha={0:.8f} '
        '+lat_0={1:.14f} +lonc={2:.14f} +k={3:.14f} '
        'proj_input > proj_output'
    ).format(alpha, lat1, lon1, k0)

    return call_proj(x, y, proj_cmd)


def projection_omerc_v2(lon, lat, lon1, lat1, lon2, lat2, k0=1):
    """
    Oblique mercator projection of a spheroid (forward). Optionally secant.

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
        '+lon_1={0:.14f} +lat_1={1:.14f} +lon_2={2:.14f} +lat_2={3:.14f} '
        'proj_input > proj_output'
    ).format(lon1, lat1, lon2, lat2)

    return call_proj(lon, lat, proj_cmd)


def call_proj(x, y, cmd):
    """System call to proj"""

    f = open('./proj_input', 'w')
    if np.isscalar(x):
        f.write('{0:.14f} {1:.14f}\n'.format(x, y))
    else:
        for i, _ in enumerate(x):
            f.write('{0:.14f} {1:.14f}\n'.format(x[i], y[i]))

    f.close()
    print('Calling \'' + cmd + '\'')
    subprocess.run(
        cmd,
        shell=True,
        check=True
    )
    res = np.loadtxt('proj_output')
    # clean up
    os.remove('proj_input')
    os.remove('proj_output')

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

    allscalars = np.all(
        np.isscalar(lon1) &
        np.isscalar(lat1) &
        np.isscalar(lon2) &
        np.isscalar(lat2)
    )

    if ~allscalars:
        assert ((lat1.ndim == 1) &
                (lon1.ndim == 1) &
                (lat2.ndim == 1) &
                (lon2.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."

    geod_cmd = (
        'geod -I -f \'% .5f\' '
        '+ellps=WGS84 '
        'geod_input > geod_output'
    )

    return call_geod(lon1, lat1, lon2, lat2, geod_cmd)


def call_geod(x1, y1, x2, y2, cmd):
    """System call to geod"""

    f = open('./geod_input', 'w')
    if np.isscalar(x1):
        f.write('{0:.30f} {1:.30f} {2:.30f} {3:.30f}\n'.format(x1, y1, x2, y2))
    else:
        for i, _ in enumerate(x1):
            f.write('{0:.30f} {1:.30f} {2:.30f} {3:.30f}\n'.format(
                y1[i], x1[i], y2[i], x2[i]))

    f.close()
    print('Calling \'' + cmd + '\'')
    subprocess.run(
        cmd,
        shell=True,
        check=True
    )

    res = np.loadtxt('geod_output')

    # clean up
    os.remove('geod_input')
    os.remove('geod_output')

    if np.isscalar(x1):
        forward_azimuth = res[0]
        dist = res[2]
    else:
        forward_azimuth = res[:, 0]
        dist = res[:, 2]

    return forward_azimuth, dist


def projection_merc(x, y, lat_ts):
    """
    Mercator

    Parameters
    ----------
    x: 1-d array of longitudes (degrees)
    y: 1-d array of latitudes (degrees)
    lat_ts: latitude of true scale
    """

    proj_cmd = (
        'proj -f \'% .30f\' '
        '+proj=merc +lat_ts={0:.14f} '
        'proj_input > proj_output'
    ).format(lat_ts)

    if not np.isscalar(x):
        assert ((x.ndim == 1) & (y.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."
    return call_proj(x, y, proj_cmd)


def projection_merc_inv(x, y, lat_ts):
    """
    Inverse Mercator

    Parameters
    ----------
    x: 1-d array of longitudes (degrees)
    y: 1-d array of latitudes (degrees)
    lat_ts: latitude of true scale
    """

    proj_cmd = (
        'proj -I -f \'% .30f\' '
        '+proj=merc +lat_ts={0:.14f} '
        'proj_input > proj_output'
    ).format(lat_ts)

    if not np.isscalar(x):
        assert ((x.ndim == 1) & (y.ndim == 1)
                ), "Aborted. Input arrays must be one-dimensional."
    return call_proj(x, y, proj_cmd)
