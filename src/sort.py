from math_utils import polar

from shapes import line, point


def polar_wrapper(inst, origin, index) -> float:
    """Function to allow one line sorting of an iterable of points and lines."""
    if isinstance(inst, line.Line):
        # is a line
        # sort wrt first point on line
        arr = inst.arrays[0]
    else:
        # is a point
        arr = inst.xy
    return polar.cart_to_polar_2D(arr, origin)[index]
