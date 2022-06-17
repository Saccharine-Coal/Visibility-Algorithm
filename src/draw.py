from __future__ import annotations

import pygame as pg
from math_utils import math_funcs


def draw_line(surface, points, point_color=(255, 255, 255), line_color=(155, 155, 155), width=2):
    pg.draw.line(surface, line_color, points[0], points[1], width=width)
    pg.draw.circle(surface, point_color, points[0], 2*width)
    pg.draw.circle(surface, point_color, points[1], 2*width)

def draw_polygon(surface, points, point_color=(255, 255, 0), line_color=(155, 155, 0), width=1):
    pg.draw.polygon(surface, line_color, points, width=width)
    for xy in points:
        pg.draw.circle(surface, point_color, xy, 2*width)

def draw_mouse(surface):
    xy = pg.mouse.get_pos()
    pg.draw.circle(surface, (255, 0, 0), xy, 5)


def scale_to_steps(vector, steps):
    step_size = vector.length() / steps
    return vector.scale_to_length(step_size)


def draw_triangle_with_gradient(surface, points, steps=10):
    return
    origin, p_1, p_2 = points[:]
    p_1, p_2 = pg.math.Vector2(p_1), pg.math.Vector2(p_2)
    v_1 = math_funcs.points_to_vector(origin, p_1)
    v_2 = math_funcs.points_to_vector(origin, p_2).normalize()
    scaled_v_1 = scale_to_steps(v_1, steps)
    scaled_v_2 = scale_to_steps(v_2, steps)
    for i in range(0, steps):
        i_1 = p_1 + (scaled_v_1 * i)
        i_2 = P_2 + (scaled_v_2 * i)
        triangle = (origin, i_1, i_2)
        pg.draw.polygon(screen, (0, 255, 0), triangle)

