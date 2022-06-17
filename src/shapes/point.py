from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shapely import geometry

from math_utils import math_funcs, numerical, iter_funcs, coord_conversions

'''
@dataclass
class PolarPoint:
    __slots__ = ("rphi")
    # (r, phi)
    rphi: np.ndarray

    @property
    def r(self) -> float:
        return self.rphi[0]

    @property
    def phi(self) -> float:
        return self.rphi[1]
'''


@dataclass
class Point:
    __slots__ = ("xy", "parent")
    xy: np.ndarray
    parent: object


    def __repr__(self) -> str:
        name = self.__class__.__name__
        has_parent = not (self.parent is None)
        return f"name(xy={self.xy}, has_parent={has_parent})"

    def __post_init__(self):
        # type and value checking
        if not isinstance(self.xy, np.ndarray):
            raise TypeError("Not an numpy ndarray")
        if len(self.xy) != 2:
            raise TypeError("Size of xy is not 2")

    def __iter__(self):
        return iter(self.xy)

    def __len__(self) -> int:
        return self.size

    def __delitem__(self, key):
        raise NotImplementedError("Slotted class.")
        #self.__delattr__(key)

    def __getitem__(self, val):
        return self.xy[val]

    def __setitem__(self, key, value):
        raise NotImplementedError("Slotted class.")
        #self.__setattr__(key, value)

    @property
    def x(self) -> float:
        return self._cart_xy[0]

    @property
    def y(self) -> float:
        return self._cart_xy[1]

    @property
    def size(self) -> int:
        return 2

    def as_polar(self, polar_origin) -> np.ndarray:
        """Get polar coordinates from cartesian xy relative to the polar_origin."""
        # (r, phi)
        return coord_conversions.cart_to_polar_2D(self.xy, origin=polar_origin)

    def get_topmost_parent(self):
        if self.parent is not None:
            return self.parent.get_topmost_parent()
        return self.parent
