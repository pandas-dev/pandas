from pathlib import Path
from typing import Any, Literal

from PIL import Image

from .._types import Writeable
from . import base

class PilImage(base.BaseImage):
    kind: Literal["PNG"]
    fill_color: str
    # the new_image and get_image methods accept arbitrary keyword arguments to
    # accommodate subclasses with additional arguments.
    def new_image(self, *, back_color: str = "white", fill_color: str = "black", **kwargs: Any) -> Image.Image: ...
    def get_image(self, **kwargs: Any) -> Image.Image: ...
    def drawrect(self, row: int, col: int) -> None: ...
    # kwargs are passed on to PIL.Image.save, which also accepts arbitrary keyword arguments.
    def save(  # type: ignore[override]
        self,
        stream: str | bytes | Path | Writeable,
        format: str | None = None,
        *,
        kind: str | None = None,
        save_all: bool = ...,
        bitmap_format: Literal["bmp", "png"] = ...,
        optimize: bool = ...,
        **kwargs: Any,
    ) -> None: ...
    # attribute access is forwarded to the wrapped PIL.Image.Image instance.
    def __getattr__(self, name: str) -> Any: ...
