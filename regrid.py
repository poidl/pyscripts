# pylint: disable=C0103

"""Helpers for regridding."""

import numpy as np


def envelope_1d(x):
    """Regrid to envelope 1d.

    Parameters
    ----------
    x: 1-d array
    """
    dx = np.diff(x)

    if np.any(np.any(dx <= 0)):
        raise Exception(
            "Aborted. Input array must be strictly monotonically increasing")

    xn = np.nan * np.ones(len(x) + 1)
    xn[1:-1] = x[:-1] + 0.5 * dx

    xn[0] = x[0] - 0.5 * dx[0]
    xn[-1] = x[-1] + 0.5 * dx[-1]

    return xn
