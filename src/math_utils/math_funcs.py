from __future__ import annotations
from typing import Union
import math
import copy

import numpy as np

import pygame as pg

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from math_utils import numerical, intersection, distance, iter_funcs

from shapes import line, point
import states
import timer


def triangulate_from_origin(origin, points, ngons, boundary) -> list:

    def two_switch_generator(starting_state=False)-> bool:
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

    def four_switch_generator(paired=False) -> bool:
        """A binary like switch that takes 4 iterations for the state to flip."""
        # paired -> False, True, True, False
        # paired -> True, False, False, True
        counter = 1
        while True:
            mod = counter % 4
            is_mod_zero = mod == 0
            if is_mod_zero:
                counter = 0
            counter += 1
            if paired:
                if mod == 1 or mod == 2:
                    yield False
                else:
                    yield True
            yield is_mod_zero


    def binary_four_switch_generator(starting_state=False) -> bool:
        """A binary like switch that takes 2 iterations for the state to flip."""
        # True, True, False, False
        counter = 0
        while True:
            if counter < 2:
                yield starting_state
            else:
                if counter == 3:
                    counter = -1
                yield not starting_state
            counter += 1
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
            tri, counter = states.ff(i, inst, next_inst, counter, ngons, boundary, origin)
            triangles.append(tri)

        else:
            break
    return triangles


def pair_with_collinear2(iterable, origin, collinear_indices: dict={}):
    """Arrange in iterable in pairs in the order of the original iterable. If collinear indices are given then they collinear indices are paired together."""
    # N = 3 => [1, 2, 3] => [(1, 2), (2, 3)] => Even size
    # N = 4 => [1, 2, 3, 4] => [(1, 2), (2, 3), (3, 4)] => odd size
    '''def points_from_indices(iterable, i_dict, key) -> list:
        points =  []
        for i in i_dict[key]:
            points.append(iterable[i])
        return points
    '''
    handle_collinear = False
    if len(collinear_indices) != 0:
        handle_collinear = True
        away = False
        # {1: a, 2: b, ..., n: z}
        index_dict = list_to_dict(iterable)
        col_dict = {}
        # remove indices from dict and iterate over keys while preserving indices
        for indices in collinear_indices.values():
            for i in indices:
                del index_dict[i]
        # set collinear indices to False and interate over new list
    size = len(iterable)
    if size > 2:
        pairs = []
        if handle_collinear:
            valid_indices = list(index_dict.keys())
            valid_size = len(valid_indices)
            for index, i in enumerate(valid_indices):
                # all i: val are valid
                if i < valid_size:
                    pass
                else:
                    break
                val = index_dict[i]
                if i in collinear_indices.keys():
                    # collect points from original set
                    col_indices = collinear_indices[i]
                    import lines
                    collinear_vals = tuple(iterable[j] for j in col_indices)

                    col_p = [val, *collinear_vals]

                    print(lines.Line(col_p, polar_origin=origin).is_collinear_with_origin())
                    # sort wrt distance
                    col_p.sort(key=lambda p: cart_to_polar_2D(p, origin)[0], reverse=away)
                    #away = not away
                else:
                    #point is not collinear
                    col_p = [val]
                # HANDLE NEXT ELEMENT
                next_val = iterable[k]
                if k in collinear_indices.keys():
                    # collect points from original set
                    next_col_indices = collinear_indices[k]
                    collinear_vals = tuple(iterable[j] for j in next_col_indices)
                    next_col_p = [next_val, *collinear_vals]
                    # sort wrt distance
                    next_col_p.sort(key=lambda p: cart_to_polar_2D(p, origin)[0], reverse=away)
                    #away = not away
                else:
                    #point is not collinear
                    next_col_p = [next_val]
                # change order depending on
                if len(next_col_p) > len(col_p):
                    # switch order if one is collinear and the other is not
                    tmp = next_col_p
                    next_col_p = col_p
                    col_p = tmp
                    if away:
                        pair_i = col_p + next_col_p
                    else:
                        pair_i = next_col_p + col_p
                else:
                    pair_i = col_p + next_col_p
                pairs.append(pair_i)
            # pair first with last
            pairs.append((iterable[0], iterable[-1]))
            return pairs
        else:
            return pair(iterable)
    elif size == 2:
        # smallest size
        return [(iterable[0], iterable[1])]
    else:
        # a iter of size 1 should not be allowed
        raise ValueError



def triangulate_from_origin2(origin, points) -> list:
    # get unique unordered list
    u_points = list(set(tuple(p) for p in points))
    # sort wrt angle around origin
    u_points.sort(key=lambda p: cart_to_polar_2D(p, origin)[1])
    col_i_dict = get_collinear_indices_dict(u_points, origin)
    paired_points = pair_with_collinear2(u_points, origin, col_i_dict)
    triangles = []
    for pair in paired_points:
        tri = (origin, *pair)
        triangles.append(tri)
    return triangles

def list_to_dict(values) -> dict:
    # make sure vals are sorted
    vals_dict = {}
    for i, val in enumerate(values):
        vals_dict[i] = val
    return vals_dict



def get_collinear_indices_dict(points, origin) -> dict:
    """
    @return index (int): indices list(int)
    """
    polar_points = list(cart_to_polar_2D(point, origin) for point in points)
    # {xy: [indices that are collinear with xy]}
    list_of_collinear_indices = []
    for i, polar in enumerate(polar_points):
        indices = []
        polar_copy = polar_points.copy()
        # preserve indices in copy
        polar_copy[i] = False
        for j, polar_cp in enumerate(polar_copy):
            if not (polar_cp is False):
                # use isclose for float comparison
                # lower the tolerance as we are in angle degrees
                if (numerical.is_close(polar[1], polar_cp[1])):
                    # lie on the same angle, thus are collinear
                    indices.append(j)
        if len(indices) > 0:
            # len will always be at least 1
            indices.append(i)
            indices = set(indices)
            # use set to compare equal lists with different orders
            if indices not in list_of_collinear_indices:
                list_of_collinear_indices.append(indices)
    indices_dict = {}
    for set_of_indices in list_of_collinear_indices:
        indices = list(set_of_indices)
        indices.sort()
        start, rest = indices[0], indices[1:]
        indices_dict[start] = rest
    return indices_dict



def cart_to_polar_2D(xy, origin=(0, 0)):
    (x, y) = xy
    x = x - origin[0]
    y = y - origin[1]
    r = math.sqrt((x)**2 + (y)**2)
    if x == 0:
        if y >= 0:
            # quadrant II => theta=90
            theta = math.pi / 2
        else:
            # quadrant III
            theta = 3*math.pi/2
    else:
        theta = math.atan(y / x)
    theta = (180 / math.pi) * theta
    return r, quadrize((x, y), theta)

def polar_to_cart_2D(r, theta, origin=(0, 0)):
    """THETA IS IN DEGREES> NEED TO CHANGE TO RADS"""
    return r*math.cos(theta*math.pi/180)+origin[0], r*math.sin(theta*math.pi/180)+origin[1]

def quadrize(xy, theta):
    """convert theta [180, -180] to [0, 360]"""
    # https://www.mathsisfun.com/polar-cartesian-coordinates.html
    (x, y) = xy
    # quadrant I
    if x >= 0 and y >= 0:
        return theta
    # quadrant II
    if x < 0 and y >= 0:
        return theta + 180
    # quadrant III
    if x < 0 and y < 0:
        return theta + 180
    # quadrant IV
    if x >= 0 and y < 0:
        return theta + 360


def get_slope(p_1, p_2) -> float:
    x_1, y_1 = p_1[:]
    x_2, y_2 = p_2[:]
    if x_2 - x_1 == 0:
        if y_2 - y_1 > 0:
            return float("inf")
        else:
            return -float("inf")
    return (y_2 - y_1) / (x_2 - x_1)

def is_in(val, boundary_1, boundary_2, inclusive=False) -> bool:
    if inclusive:
        return boundary_1 <= val <= boundary_2
    return boundary_1 < val < boundary_2






def parametrize(p_1, p_2, q_1, q_2) -> tuple:
    # https://stackoverflow.com/a/565282
    # assume numpy arrays

    if (tuple(p_1) in (tuple(q_1), tuple(q_2))):
        return 1, 1
    if (tuple(p_2) in (tuple(q_1), tuple(q_2))):
        return 1, 1
    r = p_2 - p_1
    s = q_2 - q_1
    # calc t, calc denom first
    t_denom = u_denom = np.cross(r, s)
    t_num = np.cross(q_1 - p_1, s)
    u_num = np.cross(q_1 - p_1, r)
    if t_denom == 0:
        t = u = float("inf")
    else:
        t = t_num / t_denom
        # calc u, calc denom first
        u = u_num / u_denom
    return t, u


def triangulate2(origin, points):
    # group collinear points by default
    # sort wrt angle
    sorted_points = sorted(points, key=lambda p: cart_to_polar_2D(p, origin)[1])
    print(get_collinear_indices(points, origin))
    #for i, pp in enumerate(sorted_points):



def triangulate(points, origin, group_collinear=False):
    """Group a list of points into triangles with a point in the triangle being the origin."""
    triangles = []
    points = list(set(tuple(p) for p in points))
    # sort wrt distance
    points.sort(key=lambda p: cart_to_polar_2D(p, origin)[0])
    # sort wrt angle
    points.sort(key=lambda p: cart_to_polar_2D(p, origin)[1])
    col_indices = get_collinear_indices(points, origin)
    paired_points = pair2(points, origin, collinear_indices=col_indices)
    # get collinear indices
    for paired in paired_points:
        new_triangle = (origin, *paired)
        triangles.append(new_triangle)
    # pair the beginning and end points
    #triangles.append([origin, points[0], points[-1]])
    #print(len(points), len(triangles))
    return triangles

def get_collinear_indices(points, origin) -> list:
    polar_points = list(cart_to_polar_2D(point, origin) for point in points)
    # {xy: [indices that are collinear with xy]}
    list_of_collinear_indices = []
    for i, polar in enumerate(polar_points):
        indices = []
        polar_copy = polar_points.copy()
        # preserve indices in copy
        polar_copy[i] = False
        for j, polar_cp in enumerate(polar_copy):
            if not (polar_cp is False):
                if (polar[1] == polar_cp[1]):
                    # lie on the same angle, thus are collinear
                    indices.append(j)
        if len(indices) > 0:
            # len will always be at least 1
            indices.append(i)
            indices = set(indices)
            # use set to compare equal lists with different orders
            if indices not in list_of_collinear_indices:
                list_of_collinear_indices.append(indices)
    ordered_indices = []
    for set_of_indices in list_of_collinear_indices:
        indices = list(set_of_indices)
        ordered_indices.sort()
        ordered_indices.append(indices)
    return ordered_indices


def area_of_triangle(point_a, point_b, point_c) -> float:
    """Get the area of a triangle.
    @return: float. If area is 0, then at least two points are collinear."""
    a_x, a_y = point_a[:]
    b_x, b_y = point_b[:]
    c_x, c_y = point_c[:]
    return (a_x*(b_y-c_y)+b_x*(c_y-a_y)+c_x*(a_y-b_y))/2

def points_to_vector(point_1, point_2) -> pg.math.Vector2:
    """Get a 2D directional vector that lies on the given line segment."""
    x_1, y_1 = point_1[:]
    x_2, y_2 = point_2[:]
    return pg.math.Vector2(x_2 - x_1, y_2 - y_1)


def pair2(iterable, origin, collinear_indices=None):
    """Arrange in iterable in pairs in the order of the original iterable. If collinear indices are given then they collinear indices are paired together."""
    # N = 3 => [1, 2, 3] => [(1, 2), (2, 3)] => Even size
    # N = 4 => [1, 2, 3, 4] => [(1, 2), (2, 3), (3, 4)] => odd size
    #iter_cp = list(iterable).copy()
    if collinear_indices is not None:
        indices = []
        col_dict = {}
        for col_indices in collinear_indices:
            key = str(iterable[col_indices[0]])
            val = [iterable[i] for i in col_indices[1:]]
            col_dict[key] = val
            indices.extend(col_indices[1:])
        for i in sorted(indices, reverse=True):
            iterable.pop(i)
    size = len(iterable)
    if size > 2:
        pairs = []
        for i, element in enumerate(iterable):
            if i < size-1:
                if collinear_indices is not None and str(element) in col_dict.keys():
                    #colinear_points = elements
                    # sort colinear points wrt distance
                    next_elem = iterable[i+1]
                    collinear_points = sorted((element, *col_dict[str(element)]), key=lambda p:
                                             cart_to_polar_2D(p, next_elem)[0])
                    pair_i = (next_elem, *collinear_points)
                else:
                    next_elem = iterable[i+1]
                    if collinear_indices is not None and str(next_elem) in col_dict.keys():
                        # sort colinear points wrt distance
                        collinear_points = sorted((element, *col_dict[str(next_elem)]), key=lambda p:
                                             cart_to_polar_2D(p, element)[0])
                        pair_i = (element, *collinear_points)
                    else:
                        pair_i = (element, next_elem)
                pairs.append(pair_i)
            else:
                break
        # pair first with last
        pairs.append((iterable[0], iterable[-1]))
        return pairs
    elif size == 2:
        # smallest size
        return [(iterable[0], iterable[1])]
    else:
        # a iter of size 1 should not be allowed
        raise ValueError

#def distance(point_1, point_2) -> float:
#    p_1, p_2 = np.array(point_1), np.array(point_2)
#    return math.sqrt(sum((p_2 - p_1)**2))


def order_collinear_points(origin, points) -> list:
    return sorted(points, key=lambda p: distance(origin, p))

def pair(iterable, collinear_indices=None, closed=True):
    """Arrange in iterable in pairs in the order of the original iterable. If collinear indices are given then they collinear indices are paired together."""
    # N = 3 => [1, 2, 3] => [(1, 2), (2, 3)] => Even size
    # N = 4 => [1, 2, 3, 4] => [(1, 2), (2, 3), (3, 4)] => odd size
    size = len(iterable)
    if size > 2:
        pairs = []
        for i, element in enumerate(iterable):
            if i < size-1:
                pair_i = (element, iterable[i+1])
                pairs.append(pair_i)
            else:
                break
        if closed:
            # pair first with last
            pairs.append((iterable[0], iterable[-1]))
        return pairs
    elif size == 2:
        # smallest size
        return [(iterable[0], iterable[1])]
    else:
        # a iter of size 1 should not be allowed
        raise ValueError

def get_duplicate_of_tuples(values: list):
    print("DOES NOT WORK!")
    unique = set()
    duplicates = []
    for val in values:
        tuple_val = tuple(val)
        if tuple_val not in unique:
            unique.add(tuple_val)
        else:
            # return first duplicate
            return val
    # no duplicates exist
    return None

def cast_ray(ray_pair: tuple[tuple], line_pairs: list, bounds: pg.Rect=None, vertices_only=True) -> Union[point.Point, line.Line]:
    poi_list = get_closest_intersection(ray_pair, line_pairs, bounds, vertices_only)
    if isinstance(poi_list, list):
        return None
    return poi_list

def get_point_line_intersections(point, segments: list[line.Segment]) -> list[int]:
    coll_i = []
    for i, seg in enumerate(segments):
        if seg.collidepoint(point):
            coll_i.append(i)
    #if len(coll_i) == 0:
    #    raise ValueError
    return coll_i

def minimize_segment_dist(ray_line, segments: list) -> list[line.Segment]:
    # minimize wrt to ray end to the segment
    return sorted(segments,
                  key=lambda seg: distance.point_to_line(ray_line[1],
                                                         seg.start.xy, seg.end.xy))
def get_closest_intersection_old(ray_pair: tuple[tuple], line_pairs: list, bounds: pg.Rect=None, vertices_only=True) -> Union[None, tuple]:
    """Get the closest vertex intersection for a given ray pair: (ray start, ray end).
    @return None or tuple from line_pairs"""
    bound_points = get_points_from_rect(bounds)
    poi_and_dist = []
    # get point of intersections with the euclidean distance from the ray origin
    for line_pair in line_pairs:
        infinite_ray = not vertices_only
        poi_with_dist = get_poi_with_dist(*ray_pair, *line_pair, infinite_ray=infinite_ray)
        if poi_with_dist is not None:
            poi_and_dist.append(poi_with_dist)
    # check if list is empty
    if len(poi_and_dist) != 0:
        # minimize w.r.t. euclidean dist
        min_poi = min(poi_and_dist, key=lambda tup: tup[1])[0]
        # determine whether this min_poi is a vertex
        if not vertices_only:
            return min_poi
        if is_point_a_vertex(min_poi, line_pairs):
            # only care if the point is a vertex
            # recast the ray recursively unless it is a boundary vertex
            # we need to check whether recasting the ray ends up in a polygon or not
            # ray = ray + normalized ray
            ray_start, ray_end = ray_pair[:]
            norm_ray = normalize(ray_end-ray_start)
            new_ray_end = ray_end + norm_ray
            # determine if new ray end is in the polygon the it intersects
            # assume all polygons are closed.
            # collect all the segments the min_poi intersects
            segments = get_all_connected_segments(min_poi, line_pairs, recursively=True)
            if is_array_in_list(min_poi, bound_points):
                # ray is casted to a boundary. Ray cannot continue past the boundary.
                return min_poi
            if in_polygon(new_ray_end, segments):
                # don't need to recast if the new ray is in the polygon
                return min_poi
            '''if is_vertex_a_midpoint(min_poi, line_pairs):
                return min_poi
                # poi is a midpoint'''
            # new ray is not in a polygon nor is it intersecting an endpoint so we recast the ray
            #print(min_poi, ray_pair[0], new_ray_end)
            trimmed_pairs = remove_vertex_endpoint(min_poi, line_pairs)
            ray_pair = (ray_pair[0], min_poi)
            return min_poi, get_closest_intersection(ray_pair, trimmed_pairs, bounds=bounds, vertices_only=False)
        # return poi if not vertices only
    # cast ray on bounds
    '''elif bounds is not None:
        poi_and_dist = []
        for line_pair in get_segments_from_rect(bounds):
            poi_with_dist = get_poi_with_dist(*ray_pair, *line_pair, infinite_ray=True)
            if poi_with_dist is not None:
                poi_and_dist.append(poi_with_dist)
        if len(poi_and_dist) == 0:
            return None
        return min(poi_and_dist, key=lambda tup: tup[1])[0]
    else:
        return None
    '''
    return None

#@timer.timer
def get_closest_intersection(ray_pair: tuple[tuple], segments: list, bounds: polygon.Polygon, vertices_only=True) -> list[Union[None, tuple]]:
    """Get the closest vertex intersection for a given ray pair: (ray start, ray end).
    @return None or tuple from line_pairs"""
    # sort segments wrt of point to segment
    segments = sorted(segments, key=lambda seg: distance.point_to_line(ray_pair[0], seg.start.xy, seg.end.xy, infinite=False))
    # we are always casting to a vertex which is a poi
    MIN_POI = 4
    possible_min_poi = []
    infinite_ray = not vertices_only
    for i, seg in enumerate(segments):
        arrays = seg.as_arrays()
        poi = intersection.get_intersection(*arrays, *ray_pair, q_infinite=infinite_ray)
        if poi is not None:
            dist = numerical.distance(poi, ray_pair[0])
            possible_min_poi.append((poi, dist, seg.get_topmost_parent()))
            if len(possible_min_poi) >= MIN_POI:
                break
    if len(possible_min_poi) == 0:
        # cast ray on bounds if the list is empty
        for seg in bounds.as_segments():
            poi = intersection.get_intersection(*seg.as_arrays(), *ray_pair, q_infinite=infinite_ray)
            if poi is not None:
                dist = numerical.distance(poi, ray_pair[0])
                possible_min_poi.append((poi, dist, seg.get_topmost_parent()))
                break
    # check if list is empty
    if len(possible_min_poi) != 0: #or len(poi_and_dist) != 0:
        # minimize w.r.t. euclidean dist
        possible_min_poi.sort(key=lambda p: p[1])
        min_poi = possible_min_poi[0][0]
        parent = possible_min_poi[0][2]
        #min_poi = min(poi_and_dist, key=lambda tup: tup[1])[0]
        # determine whether this min_poi is a vertex
        if not vertices_only:
            # will only occur on recursion
            return [min_poi]
        if parent is bounds:
            # hits bound so recursion is not needed
            return point.Point(min_poi, bounds)
        if parent.is_a_vertex(min_poi):
            # 2 cases: point or line
            # the parent of either case will be a polygon
            # only care if the point is a vertex
            # recast the ray recursively unless it is a boundary vertex
            # we need to check whether recasting the ray ends up in a polygon or not
            # new ray = min_poi + normalized ray
            ray_start, ray_end = ray_pair[:]
            norm_ray = normalize(ray_end-ray_start)
            # determine if new ray end is in the polygon the it intersects
            new_ray_end = min_poi + i*norm_ray
            if parent.collidepoint(new_ray_end):
                return [point.Point(min_poi, parent)]
            # is likely not in a polygon, so we recast the ray
            trimmed_segs = remove_vertex_endpoint(min_poi, segments)
            ray_pair = (ray_pair[0], min_poi)
            _iter = [min_poi]
            recursive_poi = get_closest_intersection(ray_pair, trimmed_segs, bounds=bounds, vertices_only=False)
            # a ray shoudd never cross the same polygon 2 times
            for pnt in recursive_poi:
                if parent.collidepoint(pnt):
                    return [point.Point(min_poi, parent)]
            _iter.extend(recursive_poi)
            # construct line object
            _line = line.Line(_iter, ray_pair[0], parent)
            return _line
    return []


def get_poi_with_dist(ray_start, ray_end, line_start, line_end, infinite_ray=False) -> Union[None, tuple]:
    poi = intersection.get_intersection(line_start, line_end, ray_start, ray_end, q_infinite=infinite_ray)
    if poi is not None:
        dist = np.linalg.norm(poi - ray_start)
        return (poi, dist)
    return None

def is_array_in_list(array, list_of_arrays) -> bool:
    if len(list_of_arrays) == 0:
        # empty list
        return False
    return any(np.array_equal(array, val) for val in list_of_arrays)

def is_point_a_vertex(point, segments) -> bool:
    for seg in segments:
        if seg.vertex_collidepoint(point):
            return True
    return False

def get_segments_from_rect(rect) -> list:
    rect_bounds = get_points_from_rect(rect, for_pairing=False)
    rect_segments = pair(rect_bounds, closed=False)
    return rect_segments

def get_points_from_rect(rect, for_pairing=False) -> list[np.array]:
    if for_pairing:
        rect_bounds = (np.array(rect.topleft), np.array(rect.topright), np.array(rect.bottomright), np.array(rect.bottomleft), np.array(rect.topleft))
    else:
        rect_bounds = (np.array(rect.topleft), np.array(rect.topright), np.array(rect.bottomright), np.array(rect.bottomleft))
    return rect_bounds

def is_vertex_a_midpoint(point, line_pairs) -> bool:
    # point is a vertex
    count = 0
    for line_pair in line_pairs:
        if np.array_equal(point, line_pair[0]) or np.array_equal(point, line_pair[1]):
            count += 1
            if count >= 2:
                return True
    return False


def normalize(vector) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm == 0:
       return vector
    return vector / norm


def get_all_connected_segments(point, line_pairs, recursively=False) -> list:
    # VERTICES ONLY
    unique_points = set([tuple(point)])
    connected_segments = []
    for i, pair in enumerate(line_pairs):
        # prevent set addition during iteration
        pairs_to_add = []
        for u_point in unique_points:
            # check if any unique point is in the line pair
            if is_array_in_list(u_point, pair):
                # collect connected pair
                pairs_to_add.append(pair)
                connected_segments.append(pair)
        # add the line pair to the unique set
        if len(pairs_to_add) != 0:
            for pta in pairs_to_add:
                for p in pta:
                    unique_points.add(tuple(p))
    if recursively:
        for i, pair in enumerate(line_pairs):
            # prevent set addition during iteration
            pairs_to_add = []
            for u_point in unique_points:
                # check if any unique point is in the line pair
                if is_array_in_list(u_point, pair):
                    # collect connected pair
                    pairs_to_add.append(pair)
                    connected_segments.append(pair)
            # add the line pair to the unique set
            if len(pairs_to_add) != 0:
                for pta in pairs_to_add:
                    for p in pta:
                        unique_points.add(tuple(p))
    #print(point, connected_segments, "\n", line_pairs)
    return connected_segments


def get_points_from_line_pairs(line_pairs, unique_only=False, for_shapely=False) -> list:
    points = []
    for lp in line_pairs:
        points.extend(lp)
    if unique_only:
        # do this to preserve order
        unique_points = []
        for point in points:
            if not is_array_in_list(point, unique_points):
                unique_points.append(point)
        return unique_points
    return points

def in_polygon(point, line_pairs) -> bool:
    # assume line_pairs is a closed polygon
    if len(line_pairs) == 0:
        return False
    if len(line_pairs) < 3:
        # line_segment
        return False
    p = Point(point)
    return Polygon(get_points_from_line_pairs(line_pairs, unique_only=True)).contains(p)

def is_closed(line_pairs) -> bool:
    # use shapely
    points = get_points_from_line_pairs(line_pairs)
    return Polygon(points).is_closed


def remove_vertex_endpoint(point, segments) -> list:
    trimmed_list = list(segments).copy()
    indices = []
    for i, seg in enumerate(segments):
        if seg.vertex_collidepoint(point):
            indices.append(i)
    for i in reversed(indices):
        trimmed_list.pop(i)
    return trimmed_list


if __name__ == "__main__":
    a = [1, 2, 3, 4]
    a = ((1, 1), (0, 0), (1, 1))
    poly = Polygon(a).is_closed
