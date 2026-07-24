from _typeshed import Incomplete

class _YStrip:
    y0: Incomplete
    y1: Incomplete
    slope: Incomplete
    fillColor: Incomplete
    fillColorShaded: Incomplete
    def __init__(self, y0, y1, slope, fillColor, fillColorShaded, shading: float = 0.1) -> None: ...

def mod_2pi(radians): ...

class _Segment:
    a: Incomplete
    b: Incomplete
    x0: Incomplete
    x1: Incomplete
    y0: Incomplete
    y1: Incomplete
    series: Incomplete
    i: Incomplete
    s: Incomplete
    def __init__(self, s, i, data) -> None: ...
    def intersect(self, o, I): ...

def find_intersections(data, small: int = 0): ...
