import math

import numpy as np

from math_utils import numerical


def cart_to_polar_2D(xy, origin=(0, 0)) -> np.ndarray:
    (x, y) = xy
    x = x - origin[0]
    y = y - origin[1]
    r = math.sqrt((x)**2 + (y)**2)
    # if x == 0:
    if numerical.is_close(x, 0):
        # if y >= 0:
        if numerical.is_above_or_eq(y, 0):
            # quadrant II => theta=90
            theta = math.pi / 2
        else:
            # quadrant III
            theta = 3*math.pi/2
    else:
        theta = math.atan(y / x)
    theta = (180 / math.pi) * theta
    phi = quadrize((x, y), theta)
    return np.array((r, phi))


def polar_to_cart_2D(r, theta, origin=(0, 0)) -> np.ndarray:
    """THETA IS IN DEGREES> NEED TO CHANGE TO RADS"""
    x = r*math.cos(theta*math.pi/180)+origin[0]
    y = r*math.sin(theta*math.pi/180)+origin[1]
    return np.array((x, y))


def quadrize(xy, theta) -> float:
    """convert theta [180, -180] to [0, 360]"""
    # https://www.mathsisfun.com/polar-cartesian-coordinates.html
    (x, y) = xy
    # quadrant I
    # if x >= 0 and y >= 0:
    if numerical.is_above_or_eq(x, 0) and numerical.is_above_or_eq(y, 0):
        return theta
    # quadrant II
    # if x < 0 and y >= 0:
    if x < 0 and numerical.is_above_or_eq(y, 0):
        return theta + 180
    # quadrant III
    if x < 0 and y < 0:
        return theta + 180
    # quadrant IV
    # if x >= 0 and y < 0:
    if numerical.is_above_or_eq(x, 0) and y < 0:
        return theta + 360
