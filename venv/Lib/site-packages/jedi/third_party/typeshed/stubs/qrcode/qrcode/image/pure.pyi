from _typeshed import SupportsWrite
from collections.abc import Generator
from typing import Any, Literal
from typing_extensions import TypeAlias

from . import base

# png.Writer; no types available
_Writer: TypeAlias = Any

class PyPNGImage(base.BaseImage):
    kind: str
    allowed_kinds: tuple[Literal["PNG"]]
    # the new_image and get_image methods accept arbitrary keyword arguments to
    # accommodate subclasses with additional arguments.
    def new_image(self, **kwargs: Any) -> _Writer: ...
    def get_image(self, **kwargs: Any) -> _Writer: ...
    def drawrect(self, row: int, col: int) -> None: ...
    def save(self, stream: SupportsWrite[bytes], kind: str | None = None) -> None: ...
    def rows_iter(self) -> Generator[list[int], Any]: ...
    def border_rows_iter(self) -> Generator[list[int], Any]: ...

PymagingImage = PyPNGImage
