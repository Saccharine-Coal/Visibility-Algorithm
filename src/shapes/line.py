from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from math_utils import numerical, iter_funcs, polar, intersection

from shapes import point


@dataclass
class Line:
    """Collection of n points. Points must be collinear."""
    points: tuple(point.Point)
    origin: tuple
    # parents: tuple = field(init=False,default=None)

    def __post_init__(self):
        origin = self.origin
        points = tuple(sorted(self.points, key=lambda p: polar.cart_to_polar_2D(p, origin=origin)[0]))
        #parents = tuple(p.parent for p in self.points)
        start = points[0]
        rphi = polar.cart_to_polar_2D(start, origin=origin)
        rest = points[1:]
        for r in rest:
            phi = polar.cart_to_polar_2D(r, origin=origin)[1]
            if not numerical.is_close(rphi[1], phi):
                print(tuple(polar.cart_to_polar_2D(val, start)for val in points))
                print("WARNING: given points are not collinear")
                #raise ValueError("Given points are not collinear.")
        # points are collinear
        # sort wrt to distance from start
        self.points = tuple(sorted(points, key=lambda p: polar.cart_to_polar_2D(p, start)[0]))
        self.arrays: tuple[np.ndarray] = points
        paired = iter_funcs.pair(points, closed=False)
        #self.segments: tuple[Segment] = tuple(segments)
        self.size: int = len(self.arrays)
        self.rphi = rphi
        #self.parents = parents
        #self.has_parent = all((p is not None) for p in parents)

    def __repr__(self) -> str:
        name = self.__class__.__name__
        return f"{name}(#points={len(self.points)}, parent={self.parent})"

    def is_collinear_with_origin(self) -> bool:
        """Collinear if all points lie on the same polar angle."""
        phi_vals = tuple(point.phi for point in self.points)
        return phi_vals.count(phi_vals[0]) == self.size

    def as_arrays(self, reverse=False) -> list[np.ndarray]:
        vals = list(self.arrays)
        if reverse:
            vals.reverse()
        return vals

    def as_points(self):
        return self.points
        points = []
        for seg in self.segments:
            points.extend(seg.as_points())
        return points

    def as_segments(self) -> list[Point]:
        """Get points as a collection of segment pairs.
        p1, p2, p3 ... -> [(p1, p2), (p2, p3) ...]
        """
        ...
        raise NotImplementedError

    def in_any_segment(self, point):
        raise NotImplementedError

    def get_topmost_parent(self):
        raise NotImplementedError
        return self.parent 


class Segment:
    """A line with only two points."""
    __slots__ = ("start", "end", "arrays", "parent")

    def __init__(self, start, end, parent):
        self.start: Point = point.Point(np.array(start), self)
        self.end: Point = point.Point(np.array(end), self)
        self.arrays: tuple[np.ndarray] = tuple(
            val.xy for val in (self.start, self.end)
        )
        self.parent = parent

    def __repr__(self) -> str:
        name = self.__class__.__name__
        points = [self.start, self.end]
        return f"{name}({points})"


    def __getitem__(self, val) -> point.Point:
        return (self.start, self.end)[val]

    @property
    def size(self) -> int:
        return 2

    @property
    def midpoint(self) -> np.array:
        dx_dy = self.start.xy + self.end.xy
        return dx_dy / 2

    def as_arrays(self, reverse=False) -> list[np.ndarray]:
        vals = list(self.arrays)
        if reverse:
            vals.reverse()
        return vals

    def as_points(self) -> tuple[point.Point]:
        return (self.start, self.end)

    def as_polar_arrays(self, origin, reverse=False) -> list[np.arrays]:
        return list(polar.cart_to_polar_2D(val, origin) for val in self.as_arrays(reverse=reverse))

    def collidepoint(self, point) -> bool:
        """Determine if a point is lies on the segment.
        Endpoints are included."""
        return intersection.is_point_on_segment(point, self.start.xy, self.end.xy)

    def vertex_collidepoint(self, point) -> bool:
        start = numerical.is_array_close(point, self.start)
        end = numerical.is_array_close(point, self.end)
        return (start or end)

    def get_topmost_parent(self):
        return self.parent
