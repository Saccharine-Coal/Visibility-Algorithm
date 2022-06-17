from __future__ import annotations

import pygame as pg

import numpy as np



class MouseCaster:
    """This is a ray caster whose position is bound to the mouse position."""

    @property
    def x(self) -> int:
        return self.xy[0]

    @property
    def y(self) -> int:
        return self.xy[1]

    def update(self) -> None:
        self.xy: np.ndarray = np.array(pg.mouse.get_pos())

