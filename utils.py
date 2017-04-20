# pylint: disable=C0103

"""Useful functions"""

import numpy as np


def rotate(x, y, phi):
    """Rotate grids x an y by angle phi"""
    xrot = x * np.cos(phi) - y * np.sin(phi)
    yrot = x * np.sin(phi) + y * np.cos(phi)
    return (xrot, yrot)
