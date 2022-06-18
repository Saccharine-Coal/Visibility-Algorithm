import numpy as np


from math_utils import polar
print(polar.cart_to_polar_2D((10, 9)))

from math_utils.numerical import *

print(is_in(10, 5, 15), is_above_or_eq(10, 5), is_below_or_eq(10, 15))


from shapes import line
#seg = line.Segment((0, 0), (1, 1))
#lin = line.Line([(0, 0), (1, 1), (3, 3), (10, 10), (1000, 1000)])
#print(lin)


from math_utils import distance
dist = distance.point_to_line((0, 0), (0, 0), (0, 0))
print(dist)


from math_utils import polar
from shapes import polygon, point

poly = polygon.Polygon([(0, 0), (10, 10), (5, 5)])
x = polar.get_bounding_polar_segment((-1, -1), poly)
print(x)

p = point.Point(np.array((0, 0)), None)
print(p[:], p[1])



def complementary_slice(iterable, start, stop, inclusive=False):
    """Get a slice of an iterable that is the complement of the given slice [start:stop:step=1].
    complementary_slice + slice = full iterable
    """
    if inclusive:
        return iterable[:start+1] + iterable[stop:]
    return iterable[:start] + iterable[stop:-1]

points = [1, 2, 3, 4, 5, 6]
a, b = 1, len(points)-2
inner = points[a:b+1]
outer = points[:a+1] + points[b:]
print(complementary_slice(points, a, b))
print(complementary_slice(points, a, b, inclusive=True))
