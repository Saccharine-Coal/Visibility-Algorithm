import numpy as np

from math_utils import numerical


def point_to_line(point, seg_start, seg_end, infinite=False) -> float:
    """Euclidean distance of point to segment in R^2 space..
    The segment is defined by the seg_start and seg_end.
    If infinite, then the segment is considered to be an infinite line.
    If not infinite, then the segment is finite."""
    # old code assumes line is infinite
    if infinite:
        # point = p = x0, y0
        # start = s = x1, y1
        # end = e = x2, y2
        x_0, y_0 = point[:]
        x_1, y_1 = seg_start[:]
        x_2, y_2 = seg_end[:]
        denom = numerical.distance(seg_start, seg_end)
        if numerical.is_close(denom, 0):
            # near 0
            if denom == 0:
                # division by 0 error
                # this implies that seg_start == seg_end
                # distance is between 2 points
                return numerical.distance(point, seg_start)
        # |(x2-x1)(y1-y0) - (x1-x0)(y2-y1)|
        num = abs((x_2 - x_1)*(y_1 - y_0) - (x_1 - x_0)*(y_2 - y_1))
        return num / denom
    else:
        # new code, line is finite
        # http://paulbourke.net/geometry/pointlineplane/
        # p1, p2 = seg_*
        # p3 = point
        # u = (x3 - x1)(x2 - x1) + (y3 - y1)(y2 - y1) / ||p2 - p1||
        p1, p2, p3 = seg_start, seg_end, point
        x1, y1 = p1[:]
        x2, y2 = p2[:]
        x3, y3 = p3[:]
        # reduce computation
        a, b = (x2 - x1), (y2 - y1)
        denom = np.linalg.norm(p2 - p1)
        if numerical.is_close(denom, 0):
            return True
            # p1 = p2
            # check if p1 = p2 = p3, point equality
            return numerical.is_array_close(p1, p3)
        u = (x3 - x1)*a + (y3 - y1)*b / denom
        x = x1 + u*a
        y = y1 + u*b
        return numerical.distance(point, (x, y))

