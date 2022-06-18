from __future__ import annotations
from typing import Union

import numpy as np

import pygame as pg

from math_utils import numerical, intersection, distance

from shapes import line, point

import timer


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


def cast_ray(ray_pair: tuple[tuple], line_pairs: list, bounds: pg.Rect=None, vertices_only=True) -> Union[point.Point, line.Line]:
    poi_list = get_closest_intersection(ray_pair, line_pairs, bounds, vertices_only)
    if isinstance(poi_list, list):
        return None
    return poi_list


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


def normalize(vector) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm == 0:
       return vector
    return vector / norm


def get_points_from_rect(rect, for_pairing=False) -> list[np.array]:
    if for_pairing:
        rect_bounds = (np.array(rect.topleft), np.array(rect.topright), np.array(rect.bottomright), np.array(rect.bottomleft), np.array(rect.topleft))
    else:
        rect_bounds = (np.array(rect.topleft), np.array(rect.topright), np.array(rect.bottomright), np.array(rect.bottomleft))
    return rect_bounds


def remove_vertex_endpoint(point, segments) -> list:
    trimmed_list = list(segments).copy()
    indices = []
    for i, seg in enumerate(segments):
        if seg.vertex_collidepoint(point):
            indices.append(i)
    for i in reversed(indices):
        trimmed_list.pop(i)
    return trimmed_list
