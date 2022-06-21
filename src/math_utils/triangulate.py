
from __future__ import annotations

import math

import numpy as np

from shapes import line, point, polygon
from math_utils import math_funcs
from math_utils import numerical, distance, polar

import timer

def triangulate_from_origin(origin, points, ngons, boundary) -> list:

    instances = points
    counter = None
    # counter = 0, 1, 2 counter % 2 == 1 only when counter is 2
    n = len(instances)
    triangles = []
    for i, inst in enumerate(instances):
        if i < n:
            # if element is last element, next element is first element
            if i == n - 1:
                # last element
                k = 0
            else:
                k = i + 1
            next_inst = instances[k]
            tri, counter = ff(i, inst, next_inst, counter, ngons, boundary, origin)
            triangles.append(tri)

        else:
            break
    return triangles

def is_vertex_a_midpoint(point, line_pairs) -> bool:
    # point is a vertex
    count = 0
    for line_pair in line_pairs:
        if (numerical.is_array_close(point, line_pair[0]) or
                numerical.is_array_close(point, line_pair[1])):
            count += 1
            if count >= 2:
                return True
    return False


def is_point_on_any_polygon(point, ngons: list[polygon.Polygon]) -> bool:
    for poly in ngons:
        if poly.on_any_segment(point):
            return True
    return False

def get_intersecting_polygon(point, ngons: list[polygon.Polygon]) -> polygon.Polygon:
    for poly in ngons:
        if poly.on_any_segment(point) or poly.collidepoint(point):
            return poly
    return None

def get_intersecting_segment(point, segments: list[line.Segment]) -> line.Segment:
    for seg in segments:
        if seg.collidepoint(point):
            return seg
    return None




def is_point_a_polygon_vertex(point, ngons: list[polygon.Polygon]) -> bool:
    for poly in ngons:
        if poly.is_a_vertex(point):
            return True
    return False

def binary_switch(starting_state=False)-> bool:
    """A binary like switch that takes 2 iterations for the state to flip."""
    counter = 1
    while True:
        mod = counter % 2
        is_mod_zero = mod == 0
        if is_mod_zero:
            counter = 0
        counter += 1
        if starting_state:
            yield not is_mod_zero
        yield is_mod_zero


def get_states(inst, next_inst) -> tuple[bool]:
    # 4 possible states
    # (p, p) = (True, True)
    # (p, l) = (True, False)
    # (l, p) = (False, True)
    # (l, l) = (False, False)
    first_is_point = isinstance(inst, point.Point)
    second_is_point = isinstance(next_inst, point.Point)
    return first_is_point, second_is_point


def get_ordered_points(line, towards) -> list:
    points = line.as_arrays()
    if towards:
        points.reverse()
    return points

def ff(i, inst, next_inst, counter, ngons, boundary, origin):
    point, next_point = get_states(inst, next_inst)
    # 4 possible states
    # (p, p) = (True, True)
    # (p, l) = (True, False)
    # (l, p) = (False, True)
    # (l, l) = (False, False)
    #print(inst.parent, next_inst.parent)
    index = i
    gen = binary_switch
    # if i == 0 or counter is None
    counter_is_none_or_0th = i == 0 or counter is None
    # get previous state
    if counter is not None:
        next(counter)
        towards = next(counter)
    if point:
        if next_point:
            # p, p
            collection = (origin, inst.xy, next_inst.xy)
        if not next_point:
            # p, l
            if counter_is_none_or_0th:
                towards = True
                counter = gen(starting_state=towards)
                # need to determine if the point is a midpoint
                if not boundary.is_a_vertex(inst.xy) and is_point_a_polygon_vertex(inst.xy, ngons):
                    towards = False
                    counter = gen(starting_state=towards)
            else:
                if not boundary.is_a_vertex(inst.xy) and is_point_a_polygon_vertex(inst.xy, ngons):
                    # is a polygon midpoint and not a boundary vertex
                    #towards = next(counter)
                    towards = False
                    counter = gen(starting_state=towards)
                elif boundary.is_a_vertex(inst.xy):
                    towards = True
                    counter = gen(starting_state=towards)
                else:
                    # invert the state
                    towards = next(counter)
                #print(i, towards)
            next_points = get_ordered_points(next_inst, towards)
            if not towards:
                # determine whether the next instance shares a polygon touch with inst
                same_polygon = False
                index = None
                #inst_segments = math_funcs.get_all_connected_segments(inst.xy, line_pairs, recursively=True)
                instersecting_ngon = get_intersecting_polygon(inst.xy, ngons)
                if instersecting_ngon is not None:
                    for n, point in enumerate(next_points):
                        for seg in instersecting_ngon.segments:
                            if seg.collidepoint(inst.xy):
                                same_polygon = True
                                index = n
                                break
                # 2 cases
                # closest point is on the polygon of next inst
                # furthest point is on the polygon of the next inst
                if same_polygon:
                    next_points = [next_points.pop(index)]
                else:
                    next_points = [next_points[0]]
            collection = [origin, inst.xy, *next_points]
            towards = next(counter)
    elif not point:
        points = inst.as_arrays()
        if next_point:
            # l, p
            if counter_is_none_or_0th:

                # thus the line must be away
                if not boundary.is_a_vertex(next_inst.xy):
                    towards = True
                else:
                    towards = False
                counter = gen(starting_state=towards)
                points = get_ordered_points(inst, towards)
                same_polygon = False
                inst_parent = inst.as_points()[0].get_topmost_parent()
                next_parent = next_inst.get_topmost_parent()
                if inst_parent is next_parent:
                    # same polygon
                    same_polygon = True
                    for i, point in enumerate(inst.as_points()):
                        if next_inst.get_topmost_parent().collidepoint(point):
                            index = i
                            break

                # 2 cases
                # closest point is on the polygon of next inst
                # furthest point is on the polygon of the next inst
                if same_polygon:
                    points = [points.pop(index)]
                collection = [origin, *points, next_inst.xy]
            else:
                # need to determine if the point is a midpoint
                if not boundary.is_a_vertex(next_inst.xy) and is_point_a_polygon_vertex(next_inst.xy, ngons):
                    # is not a boundary point but is a midpoint
                    towards = True
                    counter = gen(starting_state=towards)
                    points = get_ordered_points(inst, towards)
                    same_polygon = False
                    index = None
                    #next_inst_segments = math_funcs.get_all_connected_segments(next_inst.xy, line_pairs, recursively=True)
                    next_ngon = get_intersecting_polygon(next_inst.xy, ngons)
                    for n, point in enumerate(points):
                        for seg in next_ngon.segments:
                            if seg.collidepoint(point):
                                same_polygon = True
                                index = n
                                break
                    # 2 cases
                    # closest point is on the polygon of next inst
                    # furthest point is on the polygon of the next inst
                    if same_polygon:
                        points = [points.pop(index)]
                    collection = [origin, *points, next_inst.xy]
                elif boundary.is_a_vertex(next_inst.xy):
                    # is a boundary point
                    towards = not towards
                    counter = gen(starting_state=towards)
                    points = get_ordered_points(inst, towards)
                    if not towards:
                        # remove extrema points
                        points = [points[0]]
                    collection = [origin, next_inst.xy, * points]
                else:
                    # not a boundary nor mid point
                    # invert the state
                    towards = not towards
                    counter = gen(starting_state=towards)
                    points = get_ordered_points(inst, towards)
                    if not towards:
                        # remove extrema points
                        points = [points[0]]
                    collection = [origin, *points, next_inst.xy]
            towards = next(counter)
        if not next_point:
            # l, l
            # line, line. Thus line must be towards and away
            towards = True
            counter = gen(starting_state=towards)
            points = get_ordered_points(inst, towards)
            # needed to use as a flag to reduce code duplication
            collection = None
            # if line, line we need to switch the state immediately
            towards = not towards # this is wrong
            next_points = get_ordered_points(next_inst, towards)
            inst_parent = inst.as_points()[0].get_topmost_parent()
            next_parent = next_inst.as_points()[0].get_topmost_parent()
            if inst_parent is next_parent:
                # they are the same polygon so we only want the closest points
                # get the points on the each line that are on the same polygon
                if inst_parent is boundary:
                    # they are being casted to the bounds
                    # arc is not needed
                    collection = [*points, origin, next_points]
                else:
                    p = inst.as_points()[0]
                    n = next_inst.as_points()[0]
                    _points = inst_parent.get_arc(p, n, origin)
                    collection = [origin, *_points]
            else:
                # lines are not immediately intersecting the same polygon
                # points on the lines may intersect same polygons at different depths

                # {int, (list #, i in list)}
                arcs = []
                id_to_index = {}
                for i, pnt in enumerate(inst.as_points()):
                    parent = pnt.parent
                    if parent is not None:
                        id_to_index.setdefault(id(parent), [])
                        id_to_index[id(parent)].append(pnt)
                for i, pnt in enumerate(next_inst.as_points()):
                    parent = pnt.parent
                    if parent is not None:
                        id_to_index.setdefault(id(parent), [])
                        id_to_index[id(parent)].append(pnt)
                for ngon_id, points in id_to_index.items():
                    if len(points) > 1:
                        # this is a polygon that both lines intersect
                        start, end = points[:]
                        for poly in ngons:
                            if id(poly) == ngon_id:
                                arc = poly.get_arc(start, end, origin)
                                arcs.extend(arc)
                collection = [origin, *inst.as_arrays(reverse=not towards), *arcs, *next_inst.as_arrays(reverse=towards)]
                if len(arcs) == 0:
                    # lines never intersect the same polygon
                    collection = [*inst.as_arrays(reverse=not towards), origin, *next_inst.as_arrays(reverse=towards)]
            towards = next(counter)
    else:
        raise ValueError
    return collection, counter


def is_between(point, seg_start, seg_end) -> bool:
    a = dist(seg_start, point)
    b = dist(point, seg_end)
    c = dist(seg_start, seg_end)
    return numerical.is_close(a + b, c)


def dist(start, end) -> float:
    return math.sqrt(sum((end-start)**2))


def is_point_on_any_segment(point, line_pairs) -> bool:
    if line_pairs is None:
        return False
    if line_pairs is not None and len(line_pairs) == 0:
        return False
    for pair in line_pairs:
        if is_between(point, *pair):
            return True
    return False


def get_segment_from_point(point, line_pairs) -> tuple:
    """Get the segment that the point intersects.
    Return None if no such segment exists."""
    for pair in line_pairs:
        if is_between(point, *pair):
            return pair
    return None


def get_connected_segments_from_point(point, line_pairs) -> list:
    seg = get_segment_from_point(point, line_pairs)
    if seg is not None:
        start, end = seg[:]
        return math_funcs.get_all_connected_segments(start, line_pairs, recursively=False)
    return None


def get_closest_poly_equiv(points, polygon_segments) -> list:
    """Get the closest point that intersects a segment in the given polygon."""
    # points: points to test
    if polygon_segments is None:
        return None
    if polygon_segments is not None and len(polygon_segments) == 0:
        return None
    same_polygon = False
    index = None
    for i, point in enumerate(points):
        if is_point_on_any_segment(point, polygon_segments):
            same_polygon = True
            index = i
            break
    # 2 cases
    # closest point is on the polygon of next inst
    # furthest point is on the polygon of the next inst
    if same_polygon:
        return [points.pop(index)]
    return None


def flatten_segments(segments):
    flattened = []
    for pair in segments:
        for val in pair:
            flattened.append(tuple(val))
    return flattened


def from_same_polygon(point_1, point_2, line_pairs) -> bool:

    p1_segs = math_funcs.get_all_connected_segments(point_1, line_pairs, recursively=True)
    p2_segs = math_funcs.get_all_connected_segments(point_1, line_pairs, recursively=True)
    if len(p1_segs) != 0 and len(p2_segs) != 0:
        # are part of polygons
        # set rep
        u_p1_segs = set(tuple(flatten_segments(p1_segs)))
        u_p2_segs = set(tuple(flatten_segments(p2_segs)))
        if u_p1_segs == u_p2_segs:
            # they are from the same polygon
            return True
    return False
