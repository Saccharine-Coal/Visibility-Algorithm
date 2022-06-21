from __future__ import annotations

from typing import Union

import math

import numpy as np

from math_utils import numerical

import timer

#@timer.timer
def get_intersection(p_1, p_2, q_1, q_2, q_infinite=True) -> Union[float, int, None]:
    """Get the point of intersection between two line segments.
    @param: p_1, p_2 is line segment 1
    @param: q_1, q_2 is line segment 2
    @param: whether line segment 2 is infinite"""
    # handle trivial intersections
    r, s = p_2 - p_1, q_2 - q_1
    r_cross_s = cross = np.cross(r, s)
    q_sub_p = q_1 - p_1
    cross_near_zero = numerical.is_close(cross, 0)
    q_sub_p_cross_r = np.cross(q_sub_p, r)
    # 4 casses
    # 1. r x s = 0 and (q - p)xr = 0
    # 2. r x s = 0 and (q - p)xr != 0
    # 3. r x s != 0 and test t and test u
    # 4. else they do not intersect
    # only care about case 1, 3
    # if (cross == 0) and (np.cross(q_sub_p, r) == 0):
    if (cross_near_zero and
            numerical.is_close(q_sub_p_cross_r, 0)):
        # dot = np.dot(r, r)
        dot = np.linalg.norm(r)**2
        # if dot == 0:
        if dot == 0:
            # inf
            return None
        t_0_num = q_sub_p_cross_r
        t_0 = t_0_num / dot
        t_1 = t_0 + (np.dot(s, r) / dot)
        # if (0 <= t_0 <= 1) or (0 <= t_1 <= 1):
        if (numerical.is_in(t_0, 0, 1, inclusive=True) or
                numerical.is_in(t_1, 0, 1, inclusive=True)):
            # lie on the same line
            return q_2
        else:
            # colinear and disjoint
            return None
    # elif cross != 0
    if not cross_near_zero:
        # t, u = parametrize(r, s, q_sub_p)[:]
        # BEGIN parametrize
        # reduce function calls place func code in here
        if numerical.is_close(r_cross_s, 0) and r_cross_s == 0:
            t = u = float("inf")
        else:
            t_num = np.cross(q_sub_p, s)
            u_num = q_sub_p_cross_r
            t = t_num / r_cross_s
            u = u_num / r_cross_s
        # END parametrize

        # if q_infinite and (0 <= t <= 1) and (0 <= u):
        if (q_infinite and numerical.is_in(t, 0, 1, inclusive=True) and
                numerical.is_above_or_eq(u, 0)):
            # segment 2 is infinite and intersection exists
            poi = p_1 + (t*r)
            return poi
        # if (0 <= t <= 1) and (0 <= u <= 1):
        if (numerical.is_in(t, 0, 1, inclusive=True) and
                numerical.is_in(u, 0, 1, inclusive=True)):
            # intersection exists
            poi = p_1 + (t*r)
            return poi
        # they do not intersect within the size of the segments
    return None


#@timer.timer
def parametrize(r, s, q_sub_p) -> tuple[float]:
    # https://stackoverflow.com/a/565282
    # assume numpy arrays
    # handle trivial intersections
    # calc t, calc denom first
    denom = np.cross(r, s)
    if numerical.is_close(denom, 0) and denom == 0:
        t = u = float("inf")
    else:
        t_num = np.cross(q_sub_p, s)
        u_num = np.cross(q_sub_p, r)
        t = t_num / denom
        # calc u, calc denom first
        u = u_num / denom
    return t, u


def is_point_in_circle(point, circle_center, circle_radius) -> bool:
    dist_from_center = numerical.distance(circle_center, point)
    # point is in circle if dist < r or dist == r
    return numerical.is_below_or_eq(dist_from_center, circle_radius)

def line_circle_collision(p_1, p_2, circle_center, circle_radius) -> bool:
    x0, y0 = p_1[:]
    x1, y1 = p_2[:]
    h, k = circle_center[:]
    r = circle_radius
    a = (x1 - x0)**2 + (y1 - y0)**2
    b = 2*((x1 - x0)*(x0 - h) + (y1 - y0)*(y0 - k))
    c = (x0 - h)**2 + (y0 - k)**2 - r**2
    pre_root = (b**2) - (4*a*c)
    if pre_root < 0:
        #print(False, line_circle_collision(p_2, p_1, circle_center, r))
        return False
    t = 2*c / (-b + math.sqrt(pre_root))
    return numerical.is_in(t, 0, 1, inclusive=True)


def is_collinear(vector_1, vector_2) -> bool:
    """Collinear if vector_1 X vector_2 == 0"""
    # cross = np.cross(vector_1, vector_2)
    return numerical.is_close(np.cross(vector_1, vector_2), 0)


def is_point_on_segment(point, start, end) -> bool:
    """Check if a point lies on a line defined by the line segment [start, end]."""
    #   https://lucidar.me/en/mathematics/check-if-a-point-belongs-on-a-line-segment/
    # handle trivial intersections
    # if numerical.is_array_close(point, start) or numerical.is_array_close(point, end):
    #    return True
    # line segment = (a, b) = (start, end)
    ab, ac = (start - end), (start - point)
    if is_collinear(ab, ac):
        # the point lies on an infinite line
        k_ac = np.dot(ab, ac)
        k_ab = np.dot(ab, ab)
        # 3 cases
        # k_ac == 0
        # k_ac == k_ab
        # 0 < k_ac < k_ab
        # merge into one condition -> 0 < k_ac < k_ab
        return numerical.is_in(k_ac, 0, k_ab, inclusive=True)
    # is not collinear thus does not lie on the line
    return False
