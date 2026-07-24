from collections.abc import Iterable, Sequence
from typing_extensions import TypeAlias

from Xlib.protocol import rq

# Aliases used in other modules
_RGB3IntIterable: TypeAlias = Iterable[int]  # noqa: Y047
_Rectangle4IntSequence: TypeAlias = Sequence[int]  # noqa: Y047
_Segment4IntSequence: TypeAlias = Sequence[int]  # noqa: Y047
_Arc6IntSequence: TypeAlias = Sequence[int]  # noqa: Y047

# TODO: Complete all classes using WindowValues and GCValues
# Currently *object is used to represent the ValueList instead of the possible attribute types
def WindowValues(arg: str) -> rq.ValueList: ...
def GCValues(arg: str) -> rq.ValueList: ...

TimeCoord: rq.Struct
Host: rq.Struct
CharInfo: rq.Struct
FontProp: rq.Struct
ColorItem: rq.Struct
RGB: rq.Struct
Point: rq.Struct
Segment: rq.Struct
Rectangle: rq.Struct
Arc: rq.Struct
