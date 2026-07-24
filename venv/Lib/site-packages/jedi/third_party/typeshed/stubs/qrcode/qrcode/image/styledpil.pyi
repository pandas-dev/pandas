from _typeshed import SupportsRead
from pathlib import Path
from typing import Any, Literal

from PIL import Image

from .._types import Ink, Writeable
from ..main import ModulesType
from . import base
from .styles.colormasks import QRColorMask
from .styles.moduledrawers import SquareModuleDrawer
from .styles.moduledrawers.base import QRModuleDrawer

class StyledPilImage(base.BaseImageWithDrawer):
    kind: Literal["PNG"]
    color_mask: QRColorMask
    default_drawer_class: type[SquareModuleDrawer]
    embeded_image: Image.Image
    embeded_image_resample: Image.Resampling
    paint_color: Ink
    # the class accepts arbitrary additional positional arguments to accommodate
    # subclasses with additional arguments. kwargs are forwarded to the `new_image()` call
    # via the BaseImage.__init__ method.
    def __init__(
        self,
        border: int,
        width: int,
        box_size: int,
        *args: Any,
        qrcode_modules: ModulesType | None,
        module_drawer: QRModuleDrawer | str | None = None,
        eye_drawer: QRModuleDrawer | str | None = None,
        color_mask: QRColorMask = ...,
        embeded_image_path: str | bytes | Path | SupportsRead[bytes] | None = None,
        embeded_image: Image.Image | None = None,
        embeded_image_resample: Image.Resampling = ...,
        **kwargs: Any,
    ) -> None: ...
    # the new_image method accepts arbitrary keyword arguments to accommodate
    # subclasses with additional arguments.
    def new_image(self, **kwargs: Any) -> Image.Image: ...
    def draw_embedded_image(self) -> None: ...
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
