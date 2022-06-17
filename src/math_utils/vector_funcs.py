import numpy as np

from math_utils import numerical


def is_collinear(vector_1, vector_2) -> bool:
    """Collinear if vector_1 X vector_2 == 0"""
    # cross = np.cross(vector_1, vector_2)
    return numerical.is_close(np.cross(vector_1, vector_2), 0)


def is_point_on_segment(point, start, end) -> bool:
    """Check if a point lies on a line defined by the line segment [start, end]."""
    #   https://lucidar.me/en/mathematics/check-if-a-point-belongs-on-a-line-segment/
    # handle trivial intersections
    if numerical.is_array_close(point, start) or numerical.is_array_close(point, end):
        return True
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
        return numerical.is_in(k_ac, 0, k_ab, inclusive=False)
    # is not collinear thus does not lie on the line
    return False
