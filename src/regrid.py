# pylint: disable=C0103

"""Helpers for regridding."""

import numpy as np


def envelope(x, axis=-1):
    """
    Regrid to envelope.

    Parameters
    ----------
    x: array of arbitrary dimension
    axis: axis along which to regrid
    """

    shape = list(x.shape)
    shape[axis] = shape[axis] + 1
    xn = np.nan * np.ones(shape)

    dx = np.diff(x, 1, axis=axis)

    s = [slice(None)]
    t = s * x.ndim
    t1 = t.copy()
    t2l = t.copy()
    t2r = t.copy()
    t3 = t.copy()
    t1[axis] = slice(0, 1)
    t2l[axis] = slice(1, -1, None)
    t2r[axis] = slice(None, -1, None)
    t3[axis] = slice(-1, None)

    xn[tuple(t2l)] = x[tuple(t2r)] + 0.5 * dx

    xn[tuple(t1)] = x[tuple(t1)] - 0.5 * dx[tuple(t1)]
    xn[tuple(t3)] = x[tuple(t3)] + 0.5 * dx[tuple(t3)]

    return xn
