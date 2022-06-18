import math
import pygame as pg

import numpy as np


import raycaster
import draw
from shapes import polygon, line, point
from math_utils import math_funcs, polar, numerical, triangulate

import timer

# rgb(255, 255, 255)


def draw_mouse(surface, xy, DEBUG):
    pg.draw.circle(surface, (255, 200, 0), xy, 10)
    if DEBUG:
        # draw 1st quadrant
        # polar coordinates: (radius, phi) (50, 0) and (50, 90)
        r = 50
        p1 = polar.polar_to_cart_2D(r, 0, xy)
        p2 = polar.polar_to_cart_2D(r, 90, xy)
        pg.draw.circle(surface, (255, 0, 0), p1, 5)
        pg.draw.circle(surface, (0, 255, 0), p2, 5)



def darken_rgb(rgb, sub) -> tuple:

    def sub_and_limit(val, sub, ceil=255, floor=0):
        subtracted = val - sub
        return min(max(subtracted, floor), 255)
    r = sub_and_limit(rgb[0], sub)
    g = sub_and_limit(rgb[1], sub)
    b = sub_and_limit(rgb[2], sub)
    return (r, g, b)

def get_ngon(n, center, radius) -> polygon.Polygon:
    theta = np.linspace(0, 2*np.pi, n)
    circle_points = tuple((a, b) for a, b in zip(center[0] + radius*np.sin(theta), center[1] + radius*np.cos(theta)))
    return polygon.Polygon(circle_points)



def group_by_parents(instances) -> dict:
    # instances are points
    parent_to_inst = {}
    # keys are parents, values are lists
    for inst in instances:
        parent_id = id(inst.parent.parent)
        if not (parent_id in parent_to_inst.keys()):
            parent_to_inst[parent_id] = [inst]
        else:
            # parent id exists append to existing id
            # check if point is in list
            _list = parent_to_inst[parent_id]
            #print(inst)
            in_list = False
            for list_inst in _list:
                if numerical.is_array_close(inst.xy, list_inst.xy):
                    in_list = True
                    break
            if not in_list:
                _list.append(inst)
    return parent_to_inst


def arcs_from_dict(inst_dict, origin) -> list:
    arcs = []
    for instances in inst_dict.values():
        points = []
        if len(instances) < 2:
            continue
        for inst in instances:
            if isinstance(inst, line.Line):
                points.append(inst.as_points()[-1])
                raise NotImplementedError
            else:
                # is a point
                points.append((inst, inst.as_polar(origin)))
        # [
        # (cart, polar)
        # ]
        # ascending sort
        cart_points = tuple(x[0] for x in points)
        arcs.append(polar.maximize_angle(cart_points, origin))
        '''
        #points.sort(key=lambda x: x[1][0])
        points.sort(key=lambda x: x[1][1])
        left, right = points[0], points[-1]
        l_cart, l_polar = left[:]
        r_cart, r_polar = right[:]
        l_quad, r_quad = polar.get_quadrant(l_polar[1]), polar.get_quadrant(r_polar[1])
        if l_quad != r_quad:
            # are in diff quads
            # r_quad > l_quad
            if l_quad == 1 and r_quad == 4:
                # if left and right are in 1st and 4th quadrant, more needs to be done
                # shift 4 quad to 1 quad and 1 quad to 2 quad
                paired = [(a, b) for idx, a in enumerate(points) for b in points[idx + 1:]]
                # [(pair1, pair2) ...] -> [(pair1, pair2, angle) ...]
                pair_with_angle = []
                for pair in paired:
                    p1, p2 = pair[:]
                    ang = polar.get_major_angle(p1[1][1], p2[1][1])
                    dangle = 360-ang
                    #if dangle < 180:
                    #    dangle = 360 - dangle
                    #print(dangle)
                    pair_with_angle.append((*pair, dangle))

                max_pair = max(pair_with_angle, key=lambda x: x[2])
                left, right = max_pair[0], max_pair[1]
                l_cart, l_polar = left[:]
                r_cart, r_polar = right[:]
        if not numerical.is_array_close(l_cart, r_cart):
            points = tuple(p[0] for p in points)
            print(l_cart, r_cart, polar.maximize_angle(points, origin))
            arcs.append((l_cart, r_cart))
        '''
    # arcs are cartesian
    return arcs



#@timer.timer
def cast_rays(ray_pairs, segments, caster, bounding_ngon):
    intersections = []
    # sort wrt angle
    ray_pairs.sort(key=lambda pair: polar.cart_to_polar_2D(pair[0], caster.xy)[1])
    # sort wrt dist
    ray_pairs.sort(key=lambda pair: polar.cart_to_polar_2D(pair[1], caster.xy)[0])
    for ray_pair in ray_pairs:
        intersection_point = math_funcs.cast_ray(ray_pair, segments, bounding_ngon)
        if intersection_point is not None:
            intersections.append(intersection_point)
    return intersections


def draw_triangles(screen, DEBUG, triangles, hollow=False):
    for i, tri in enumerate(triangles):
        if DEBUG:
            mod = i % 6
            if mod == 0:
                pc = (255, 0, 0)
            if mod == 1:
                pc = (0, 255, 0)
            if mod == 2:
                pc = (0, 0, 255)
            if mod == 3:
                # rgb(255, 255, 0)
                pc = (255, 255, 0)
            if mod == 4:
                # rgb(255, 0, 255)
                pc = (255, 0, 255)
            if mod == 5:
                pc = (0, 255, 255)
            if not hollow:
                pg.draw.polygon(screen, darken_rgb(pc, 150), tri, width=0)
            pg.draw.lines(screen, pc, True, tri, width=4)
            for pnt in tri:
                pg.draw.circle(screen, darken_rgb(pc, 100), pnt, 15)
        else:
            # draw triangles nicely
            YELLOW = (255, 255, 0)
            if not hollow:
                pg.draw.polygon(screen, YELLOW, tri, width=0)
                pg.draw.polygon(screen, YELLOW, tri, width=4)
            else:
                pg.draw.polygon(screen, YELLOW, tri, width=4)

def main(DEBUG=False, PRETTY=True):
    # PRETTY is fancy graphics and no debug utils
    pg.init()

    # Set up the drawing window
    screen = pg.display.set_mode([500, 500])

    # Run until the user asks to quit
    running = True
    clock = pg.time.Clock()
    # game states
    center = (250, 250)
    x, y = 300, 300
    x3, y3 = 100, 100
    r = 25
    r3 = 10
    #wall_1 = polygon.Polygon([(100+0, y+(5*r)), (100+(4.755*r), 100+(1.545*r))])#, (x+(2.939*r), y+(-4.045*r)), (x+(-2.939*r), y+(-4.045*r)), (x+(-4.755*r), y+(1.545*r))])

    wall_2 = polygon.Polygon([(x+0, y+(5*r)), (x+(4.755*r), y+(1.545*r)), (x+(2.939*r), y+(-4.045*r)), (x+(-2.939*r), y+(-4.045*r)), (x+(-4.755*r), y+(1.545*r))])

    wall_3 = polygon.Polygon([(x3+0, y3+(5*r3)), (x3+(4.755*r3), y3+(1.545*r3)), (x3+(2.939*r3), y3+(-4.045*r3)), (x3+(-2.939*r3), y3+(-4.045*r3)), (x3+(-4.755*r3), y3+(1.545*r3))])
    r = 50
    theta = np.linspace(0, 2*np.pi, 4)
    circle_points = tuple((a, b) for a, b in zip(300 + r*np.sin(theta), 100 + r*np.cos(theta)))
    circle = polygon.Polygon(circle_points)
    circle_2 = get_ngon(6, (300, 300), 115)
    circle_3 = get_ngon(9, (100, 300), 50)
    circle_4 = get_ngon(5, (100, 425), 30)
    circle_5 = get_ngon(11, (400, 400), 30)
    walls = [wall_3, circle, circle_2, circle_3, circle_4, circle_5]
    # bounding rect
    bounding_rect = pg.Rect((0, 0,), (450, 450))
    bounding_rect.center = screen.get_rect().center
    boundary_ngon = polygon.Polygon(math_funcs.get_points_from_rect(bounding_rect, for_pairing=False))
    # ray caster origin
    caster = raycaster.MouseCaster()
    # get all vertices here, then construct ray pairs in game loop
    vertices = []
    boundary_points = []
    segments = []
    # polygons
    for wall in walls:
        for seg in wall.segments:
            vertices.extend(seg.as_points())
            segments.append(seg)
    # bounding rect
    for seg in boundary_ngon.segments:
        boundary_points.extend(seg.as_arrays())
        #vertices.extend(seg.as_points())
        #segments.append(seg)
    id_to_list = group_by_parents(vertices)
    caster.update()
    # convert list to tuples
    vertices = tuple(vertices)
    segments = tuple(segments)
    boundary_points = tuple(boundary_points)
    # reduce draw calls, since polygons are static
    polygon_surface = pg.Surface(screen.get_rect().size)
    polygon_surface.set_colorkey((0, 0, 0))
    # visibility surface
    vis_surface = pg.Surface(screen.get_rect().size, pg.SRCALPHA)
    #vis_surface.set_colorkey((0, 0, 0)) # set color key
    # radial surface
    rad_surface = pg.Surface(screen.get_rect().size, pg.SRCALPHA)   # per-pixel alpha
    # boundary rect
    pg.draw.rect(polygon_surface, (255, 255, 255), bounding_rect, width=4)
    pg.draw.rect(polygon_surface, (255, 255, 255), bounding_rect.inflate(8, 8), width=4)
    for wall in walls:
        draw.draw_polygon(polygon_surface, wall.as_arrays(), line_color=(255, 255, 255), point_color=(200, 0, 0), width=3)
    previous_pos = [0, 0]
    # GAME LOOP
    intersections = []
    arcs = []
    draw_vis = True
    while running:
        clock.tick(30)
        print(f'fps={clock.get_fps():.0f}')
        caster.update()
        new_pos = caster.xy
        update = False
        if not numerical.is_array_close(previous_pos, new_pos):
            update = draw_vis = True
            # check if caster is in any polygons
            # stuff breaks, so it's easier to prevent
            # raycasting when this occurs
            for ngon in walls:
                if ngon.collidepoint(caster.xy):
                    update = draw_vis = False
                    break
            # check if outside bounding ngon
            if not boundary_ngon.collidepoint(caster.xy):
                update = draw_vis = False
        previous_pos = new_pos
        # Did the user click the window close button?
        for event in pg.event.get():
            if event.type == pg.QUIT:
                running = False
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_ESCAPE:
                    running = False
        if update:
            arcs = arcs_from_dict(id_to_list, caster.xy)
            # these arcs are segments
            ray_pairs = []
            for arc in arcs:
                for p in arc:
                    ray_pair = (caster.xy, p.xy)
                    ray_pairs.append(ray_pair)
            ray_pairs.extend(tuple((caster.xy, p) for p in boundary_ngon.as_arrays()))
            # THIS IS THE LARGEST BOTTLENECK
            intersections = cast_rays(ray_pairs, segments, caster, boundary_ngon)
        screen.fill((0, 0, 0))

        if len(intersections) > 1 and draw_vis:
            # sort wrt dist
            intersections.sort(key=lambda p: polar.polar_wrapper(p, caster.xy, 0))
            # sort wrt angle
            intersections.sort(key=lambda p: polar.polar_wrapper(p, caster.xy, 1))
            triangles = triangulate.triangulate_from_origin(caster.xy, intersections, walls, boundary_ngon)
            if PRETTY:
                # clean vis surface
                vis_surface.fill((0, 0, 0, 0))
                draw_triangles(vis_surface, DEBUG, triangles, hollow=False)
                rad_surface.fill((0, 0, 0, 0))
                # draw expanding circles
                radius = 5
                yellow = (255, 255, 0)
                START, STOP, STEP = 60, -1, -1
                RANGE = abs(abs(START) - abs(STEP)) + 1
                for i in range(START, STOP, STEP):
                    val = 255 - (255*(i / (RANGE)))
                    alpha = min(val, 255)
                    color = (*yellow, alpha)
                    pg.draw.circle(rad_surface, color, caster.xy, radius*i)
                vis_surface.blit(rad_surface, (0, 0), special_flags=pg.BLEND_RGBA_MIN)
                screen.blit(vis_surface, (0, 0))
            else:
                draw_triangles(screen, DEBUG, triangles, hollow=False)
        if (DEBUG and draw_vis):
            COLOR = (0, 255, 0)
            for pnt in intersections:
                if isinstance(pnt, line.Line):
                    for i, p in enumerate(pnt.as_arrays()):
                        color = darken_rgb(COLOR, i*75)
                        pg.draw.circle(screen, color, p, 10)
                        pg.draw.line(screen, COLOR, caster.xy, p)
                else:
                    pg.draw.circle(screen, (0, 255, 0), pnt.xy, 10)
                    pg.draw.line(screen, COLOR, caster.xy, pnt.xy)
            for arc in arcs:
                pg.draw.line(screen, (0, 255, 255), arc[0], arc[1], width=10)

        screen.blit(polygon_surface, (0, 0))
        draw_mouse(screen, caster.xy, DEBUG)
        pg.display.flip()

    # Done! Time to quit.
    pg.quit()
    return 0


main(DEBUG=False, PRETTY=True)
