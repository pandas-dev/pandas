from _typeshed import Incomplete
from abc import ABC
from collections.abc import Iterable
from typing import Final, Literal

from .drawing import DeviceCMYK, DeviceGray, DeviceRGB
from .fpdf import FPDF
from .syntax import Name, PDFObject

class Pattern(PDFObject):
    type: Name
    pattern_type: int
    def __init__(self, shading: LinearGradient | RadialGradient) -> None: ...
    @property
    def shading(self) -> str: ...

class Type2Function(PDFObject):
    function_type: Final = 2
    domain: str
    c0: str
    c1: str
    n: int
    def __init__(self, color_1, color_2) -> None: ...

class Type3Function(PDFObject):
    function_type: Final = 3
    domain: str
    bounds: str
    encode: str
    n: int

    def __init__(self, functions: Iterable[Incomplete], bounds: Iterable[Incomplete]) -> None: ...
    @property
    def functions(self) -> str: ...

class Shading(PDFObject):
    shading_type: Literal[2, 3]
    background: str | None
    color_space: Name
    coords: list[int]
    function: str
    extend: str
    def __init__(
        self,
        shading_type: Literal[2, 3],
        background: DeviceRGB | DeviceGray | DeviceCMYK | None,
        color_space: str,
        coords: list[int],
        function: Type2Function | Type3Function,
        extend_before: bool,
        extend_after: bool,
    ) -> None: ...

class Gradient(ABC):
    color_space: str
    colors: list[Incomplete]
    background: Incomplete | None
    extend_before: Incomplete
    extend_after: Incomplete
    bounds: Incomplete
    functions: Incomplete
    pattern: Pattern
    coords: Incomplete | None
    shading_type: int

    def __init__(self, colors, background, extend_before, extend_after, bounds): ...
    def get_shading_object(self) -> Shading: ...
    def get_pattern(self) -> Pattern: ...

class LinearGradient(Gradient):
    coords: list[str]
    shading_type: int
    def __init__(
        self,
        fpdf: FPDF,
        from_x: float,
        from_y: float,
        to_x: float,
        to_y: float,
        colors: list[Incomplete],
        background=None,
        extend_before: bool = False,
        extend_after: bool = False,
        bounds: list[int] | None = None,
    ) -> None: ...

class RadialGradient(Gradient):
    coords: list[str]
    shading_type: int
    def __init__(
        self,
        fpdf: FPDF,
        start_circle_x: float,
        start_circle_y: float,
        start_circle_radius: float,
        end_circle_x: float,
        end_circle_y: float,
        end_circle_radius: float,
        colors: list[Incomplete],
        background=None,
        extend_before: bool = False,
        extend_after: bool = False,
        bounds: list[int] | None = None,
    ): ...
