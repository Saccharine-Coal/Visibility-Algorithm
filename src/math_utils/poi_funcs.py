from __future__ import annotations


def get_scalar_params(p_1, p_2, p_3, p_4) -> tuple[float, float]:
    # https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment

    near_zero = 0.000001
    # line_1 = (P_1, P_2)
    # line_2 = (P_2, P_3)
    x_1, y_1 = p_1[:]
    x_2, y_2 = p_2[:]
    x_3, y_3 = p_3[:]
    x_4, y_4 = p_4[:]
    # t = t_num / t_denom
    t_num = (x_1 - x_3)*(y_3 - y_4) - (y_1 - y_3)*(x_3 - x_4)
    if t_num == 0:
        t = 0
    else:
        t_denom = (x_1 - x_2)*(y_3 - y_4) - (y_1 - y_2)*(x_3 - x_4)
        if t_denom == 0:
            t = t_num / near_zero
        else:
            t = t_num / t_denom
    # u = u_num / u_denom
    u_num = (x_1 - x_3)*(y_1 - y_2) - (y_1 - y_3)*(x_1 - x_2)
    if u_num == 0:
        u = 0
    else:
        u_denom = (x_1 - x_2)*(y_3 - y_4) - (y_1 - y_2)*(x_3 - x_4)
        if u_denom == 0:
            u = u_num / near_zero
        else:
            u = u_num / u_denom
    return t, u

def check_collision(t, u) -> bool:
    # assume t is finite [0, 1]
    # assume u is infinite [0, infinity]
    collision: bool = False
    if u >= 0:
        if 0 <= t <= 1:
            collision: bool = True
    return collision

def get_poi(p_1, p_2, t):
    x_1, y_1 = p_1[:]
    x_2, y_2 = p_2[:]
    poi_x = x_1 + t*(x_2 - x_1)
    poi_y = y_1 + t*(y_2 - y_1)
    return poi_x, poi_y
