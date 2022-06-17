from __future__ import annotations

import math

from math_utils import numerical

from shapes import polygon, line, point


def cart_to_polar_2D(xy, origin=(0, 0)) -> tuple[float]:
    """2-dimensional cartesian to polar conversion."""
    if numerical.is_array_close(xy, origin):
        # when given xy is the origin
        return 0, 0
    (x, y) = xy
    x = x - origin[0]
    y = y - origin[1]
    r = math.sqrt((x)**2 + (y)**2)
    # x == 0
    if numerical.is_close(x, 0):
        # y > 0
        if numerical.is_above(y, 0):
            # quadrant II => theta=90
            theta = math.pi / 2
        # y == 0
        elif numerical.is_close(y, 0):
            theta = 0
        else:
            # quadrant III
            theta = 3*math.pi/2
    else:
        theta = math.atan(y / x)
    theta = (180 / math.pi) * theta
    return r, quadrize((x, y), theta)


def polar_to_cart_2D(r, phi, origin=(0, 0), in_rads=False) -> tuple[float]:
    """THETA IS IN DEGREES> NEED TO CHANGE TO RADS"""
    x_0, y_0 = origin[:]
    if not in_rads:
        # is in degrees, need to convert to radians
        phi = phi * (math.pi / 180)
    x = r*math.cos(phi) + x_0
    y = r*math.sin(phi) + y_0
    return x, y


def quadrize(xy, theta) -> tuple[float]:
    """convert theta [180, -180] to [0, 360]"""
    # https://www.mathsisfun.com/polar-cartesian-coordinates.html
    x, y = xy[:]
    # quadrant I
    # if x >= 0 and y >= 0:
    if numerical.is_above_or_eq(x, 0) and numerical.is_above_or_eq(y, 0):
        return theta
    # quadrant II
    # if x < 0 and y >= 0:
    if numerical.is_below(x, 0) and numerical.is_above_or_eq(y, 0):
        return theta + 180
    # quadrant III
    # if x < 0 and y < 0:
    if numerical.is_below(x, 0) and numerical.is_below(y, 0):
        return theta + 180
    # quadrant IV
    # if x >= 0 and y < 0:
    if numerical.is_above_or_eq(x, 0) and numerical.is_below(y, 0):
        return theta + 360

def get_bounding_polar_segment(origin, ngon: polygon) -> line.Segment:
    polar_arrays = tuple(seg.as_polar_arrays(origin) for seg in ngon.as_segments())
    sorted_arrays = sorted(polar_arrays, key=lambda x: x[0], reverse=True)
    sorted_arrays = sorted(sorted_arrays, key=lambda x: x[1], reverse=True)
    left, right = sorted_arrays[0], sorted_arrays[-1]
    return line.Segment(left, right, None)


def is_in_arc(point, start, end) -> bool:
    # assume points are all in polar form
    r0, phi0 = point[:]
    r1, phi1 = start[:]
    r2, phi2 = end[:]
    # only care about arc angle and not length
    print(is_in_arc(phi0, phi1, phi2))
    return numerical.is_in(phi0, phi1, phi2, inclusive=True)


def is_in_slice():
    """If a point is in a pie slice. Angle and radius are taken into account."""
    raise NotImplementedError


def is_in_arc(angle, phi_start, phi_end) -> bool:
    """If a point lies in an arc. The point's radius is ignored."""
    # https://stackoverflow.com/a/51896645
    a, s, e = angle, phi_start, phi_end
    # 4 cases
    # if s <= e and s <= a <= e
    # if s => e
    #   if a => s
    #   if a <= e
    # else point lies outside the arc
    if numerical.is_below_or_eq(s, e) and numerical.is_in(a, s, e, inclusive=True):
        return True
    if numerical.is_above_or_eq(s, e):
        if numerical.is_above_or_eq(a, s):
            return True
        if numerical.is_below_or_eq(a, e):
            return True
    return False


def get_quadrant(phi) -> int:
    """convert theta [0, 360] to [1, 2, 3, 4]"""
    phi = phi % 360.0
    # quadrant I
    # if 0 <= x <= 90
    if numerical.is_in(phi, 0, 90, inclusive=True):
        return 1
    # quadrant II
    # if 90 < x <= 180
    elif numerical.is_in(phi, 90, 180, inclusive=True):
        return 2
    # quadrant III
    # if 180 < x <= 270
    elif numerical.is_in(phi, 180, 270, inclusive=True):
        return 3
    # quadrant IV
    # if 270 < x <= 360
    elif numerical.is_in(phi, 270, 360, inclusive=True):
        return 4
    else:
        raise ValueError

def absolute_angle(phi_1, phi_2) -> float:
    """For abs(phi1-phi2)."""
    #shifted_1, shifted_2 = (phi_1 + 90) % 360, (phi_2 + 90) % 360
    #return abs(shifted_1 - shifted_2)
    quadrant_1 = get_quadrant(phi_1)
    quadrant_2 = get_quadrant(phi_2)
    if quadrant_1 > quadrant_2:
        # swap vars
        quadrant_1, quadrant_2 = quadrant_2, quadrant_1
        phi_1, phi_2 = phi_2, phi_1
    # phi2 is larger than phi1
    if quadrant_1 == quadrant_2:
        # same quadrant
        return abs(phi_1 - phi_2)
    if quadrant_2 == 4:
        # shift the angles by 90 degrees
        # 1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1
        #print(abs(((phi_1 + 90) % 360) - ((phi_2 + 90) % 360)), abs((phi_1 + 90) - ((phi_2 + 90) % 360)))
        return abs((phi_1 + 90) - ((phi_2 + 90) % 360))
    return abs(phi_1 - phi_2)

def get_major_angle(start, end) -> float:
    """Get the major angle, where major + minor = 360. Major >= minor."""
    dphi = abs(start - end)
    if numerical.is_above_or_eq(dphi, 180):
        return dphi
    return 360 - dphi

def get_arclength(rphi_start, rphi_end) -> float:
    """Get the arclength of an angle defined by start, end.
    the midpoint of r is used. r = r1 + r2 / 2
    The major angle is used. for s = r*phi"""
    r1, phi1 = rphi_start[:]
    r2, phi2 = rphi_end[:]
    r = (r1 + r2) / 2
    return r * get_major_angle(phi1, phi2)

def get_midpoint(val_1, val_2) -> float:
    return (val_1 + val_2) / 2


def maximize_angle(cart_points, origin) -> list:
    # as polar
    # (rad, phi)
    as_polar = list(tuple(cart_to_polar_2D(x.xy, origin)) for x in cart_points)
    polar_to_cart = {key: val for key, val in zip(as_polar, cart_points)}
    # sort wrt angle
    #as_polar.sort(key=lambda x: x[1])
    paired = [(a, b) for idx, a in enumerate(as_polar) for b in as_polar[idx + 1:]]
    # sort wrt distance
    paired.sort(key=lambda pair: get_midpoint(pair[0][0], pair[1][0]))
    # sort wrt angle
    max_pair = min(paired, key=lambda pair: get_major_angle(pair[0][1], pair[1][1]))
    #max_pair = max(paired, key=lambda pair: get_arclength(pair[0], pair[1]))
    return [
        polar_to_cart.get(max_pair[0]), polar_to_cart.get(max_pair[1])
    ]

