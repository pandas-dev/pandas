import _tkinter
from _typeshed import StrPath
from typing import ClassVar

class ThemedWidget:
    pixmap_themes: ClassVar[list[str]]
    PACKAGES: ClassVar[dict[str, str]]
    tk: _tkinter.TkappType
    png_support: bool
    def __init__(self, tk_interpreter: _tkinter.TkappType, gif_override: bool = False) -> None: ...
    def set_theme(self, theme_name: str) -> None: ...
    def get_themes(self) -> list[str]: ...
    @property
    def themes(self) -> list[str]: ...
    @property
    def current_theme(self) -> str: ...
    def set_theme_advanced(
        self,
        theme_name: str,
        brightness: float = 1.0,
        saturation: float = 1.0,
        hue: float = 1.0,
        preserve_transparency: bool = True,
        output_dir: StrPath | None = None,
        advanced_name: str = "advanced",
    ) -> None: ...
