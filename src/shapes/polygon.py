from __future__ import annotations

import numpy as np

from shapely import geometry

from math_utils import numerical, iter_funcs, distance, polar

from shapes import line, point


class Polygon:
    """Closed collection of segments. """
    __slots__ = ("arrays", "segments", "is_closed", "ngon")

    def __init__(self, points: list[tuple]):
        if len(points) <= 2:
            raise ValueError("Size of Polygon must be at least 3.")
        points = iter_funcs.get_unique_ordered(points)
        CLOSED = True
        self.arrays = tuple(np.array(val) for val in points)
        line_pairs = iter_funcs.pair(self.arrays, closed=CLOSED)
        segments = []
        for lp in line_pairs:
            seg = line.Segment(lp[0], lp[1], self)
            segments.append(seg)
        self.segments = tuple(segments)
        self.is_closed = CLOSED
        self.ngon = geometry.Polygon(self.arrays)
        if self.size != self._nsegments:
            raise ValueError("# of segments != # of points")

    def __repr__(self) -> str:
        return f"Polygon(#size={self.size}, closed={self.is_closed})"

    @property
    def size(self) -> int:
        return len(self.arrays)

    @property
    def _nsegments(self) -> int:
        return len(self.segments)

    def as_segments(self) -> tuple[Segment]:
        return self.segments

    def as_arrays(self) -> tuple[np.ndarray]:
        return self.arrays

    def as_points(self) -> tuple[point.Point]:
        """Get the unique points of the polygon"""
        points = []
        for seg in self.segments:
            points.append(seg.start)
        return tuple(points)

    def in_polygon(self, point) -> bool:
        '''
        sides = []
        for seg in self.as_segments():
            a, b, = seg.start.xy, seg.end.xy
            ab = b - a
            ap = point - a
            cross = np.cross(ab, ap)
            #print(cross, ab, ap)
            if numerical.is_above(cross, 0):
                sides.append(True)
                if False in sides:
                    break
            else:
                sides.append(False)
                if True in sides:
                    break
        # if True, all are True
        # if True, all are False
        if all(sides) or not any(sides):
            return True
        '''
        p = geometry.Point(point)
        return self.ngon.contains(p)

    def collidepoint(self, point) -> bool:
        if self.on_any_segment(point):
            return True
        return self.in_polygon(point)

    def on_any_segment(self, point) -> bool:
        # handle trivial intersections
        if self.is_a_vertex(point):
            return True
        for seg in self.segments:
            if seg.collidepoint(point):
                return True
        return False

    def is_a_vertex(self, point) -> bool:
        for arr in self.arrays:
            if numerical.is_array_close(point, arr):
                return True
        return False

    def get_topmost_parent(self):
        return self


    def get_arc(self, start, end, origin) -> list:

        def mean(vals) -> float:
            return sum(vals) / len(vals)

        def xy_mean(tuples) -> tuple:
            x_vals = tuple(xy[0] for xy in tuples)
            y_vals = tuple(xy[1] for xy in tuples)
            return (mean(x_vals), mean(y_vals))
        if numerical.is_array_close(start.xy, end.xy):
            return (start, end)

        if self.size == 3:
            # is a triangle
            min_point = min(self.arrays, key=lambda p: numerical.distance(p, origin))
            return (start, min_point, end)
            if False or not (self.is_a_vertex(start) and self.is_a_vertex(end)):
                # both are not vertices
                
                return (start, end)
            else:
                arrays = list(self.as_arrays())
                for i, arr in enumerate(reversed(self.as_arrays())):
                    j = self.size - 1 - i
                    if numerical.is_array_close(start.xy, arr):
                        arrays.pop(j)
                    if numerical.is_array_close(end.xy, arr):
                        arrays.pop(j)


                min_point = min(arrays, key=lambda p: numerical.distance(p, origin))
                return (start, min_point, end)
        if True: #or self.on_any_segment(start) and self.on_any_segment(end):
            #print(self.is_a_vertex(start), self.is_a_vertex(end))
            arrays = list(self.as_arrays())
            for i, arr in enumerate(reversed(self.as_arrays())):
                j = self.size - 1 - i
                if numerical.is_array_close(start.xy, arr):
                    arrays.pop(j)
                if numerical.is_array_close(end.xy, arr):
                    arrays.pop(j)

            arc = polar.in_cart_arc(start, end, arrays, origin)
            arc = sorted(arc, key=lambda p: polar.cart_to_polar_2D(p, origin)[1])
            return (start, *arc, end)
            # sort wrt dist
            #arc = sorted(arc, key=lambda p: polar.cart_to_polar_2D(p, origin)[0])
            # sort wrt angle
            return sorted(arc, key=lambda p: polar.cart_to_polar_2D(p, origin)[1])
            return (start, end)
            a = b = -1
            points = self.as_points()
            segments = self.as_segments()
            # not a vertex
            min_seg = min(self.segments, key=lambda seg: distance.point_to_line(start, *seg.as_arrays(), infinite=False))
            a = self.segments.index(min_seg)
            if numerical.is_above(numerical.distance(start, min_seg.arrays[0]), numerical.distance(start, min_seg.arrays[1])):
                a += 1
            min_seg = min(self.segments, key=lambda seg: distance.point_to_line(end, *seg.as_arrays(), infinite=False))
            b = self.segments.index(min_seg)
            if numerical.is_above(numerical.distance(end, min_seg.arrays[0]), numerical.distance(end, min_seg.arrays[1])):
                b += 1
            if a == b:
                return (start, points[a], end)
            if a > b:
                # swap indices
                a, b = b, a
            # do not want to include endpoints as they are
            # defined by start, end
            inner = points[a:b+1]
            inner = polar.in_cart_arc(start, end, inner, origin)
            outer = points[b:] + points[:a+1]
            outer = polar.in_cart_arc(start, end, outer, origin)
            if not inner or not outer:
                # inner is empty
                if outer:
                    return (start, *outer, end)
                if inner:
                    return (start, *inner, end)
                return (start, end)

            i_mean = xy_mean(inner)
            o_mean = xy_mean(outer)
            i_dist = numerical.distance(origin, i_mean)
            o_dist = numerical.distance(origin, o_mean)
            if i_dist < o_dist:
                arc = inner
            else:
                arc = outer
            return (start, *arc, end)
