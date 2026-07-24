"""
Color-related factories. They produce Styles.
"""

from __future__ import annotations

__lazy_modules__ = {f"{__spec__.parent}.names", "functools", "operator"}

import functools
import operator
import sys
import typing
from typing import Any, Generic, TextIO, TypeVar

from .names import color_names, default_styles
from .styles import ColorNotFound, Style

if typing.TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping

    from plumbum._compat.typing import Self

__all__ = ["ColorFactory", "StyleFactory"]


def __dir__() -> list[str]:
    return list(__all__)


S = TypeVar("S", bound=Style)


class ColorFactory(Generic[S]):
    """This creates color names given fg = True/False. It usually will
    be called as part of a StyleFactory."""

    def __init__(self, fg: bool, style: type[S]):
        self._fg = fg
        self._style = style
        self.reset = style.from_color(style.color_class(fg=fg))

        # Adding the color name shortcuts for foreground colors
        for item in color_names[:16]:
            setattr(
                self, item, style.from_color(style.color_class.from_simple(item, fg=fg))
            )

    def __getattr__(self, item: str) -> S:
        """Full color names work, but do not populate __dir__."""
        try:
            return self._style.from_color(self._style.color_class(item, fg=self._fg))
        except ColorNotFound:
            raise AttributeError(item) from None

    def full(self, name: str | int) -> S:
        """Gets the style for a color, using standard name procedure: either full
        color name, html code, or number."""
        return self._style.from_color(
            self._style.color_class.from_full(name, fg=self._fg)
        )

    def simple(self, name: str | int) -> S:
        """Return the extended color scheme color for a value or name."""
        return self._style.from_color(
            self._style.color_class.from_simple(name, fg=self._fg)
        )

    def rgb(self, r: str | int, g: int | None = None, b: int | None = None) -> S:
        """Return the extended color scheme color for a value."""
        if g is None and b is None and isinstance(r, str):
            return self.hex(r)

        return self._style.from_color(self._style.color_class(r, g, b, fg=self._fg))

    def hex(self, hexcode: str) -> S:
        """Return the extended color scheme color for a value."""
        return self._style.from_color(
            self._style.color_class.from_hex(hexcode, fg=self._fg)
        )

    def ansi(self, ansiseq: str) -> S:
        """Make a style from an ansi text sequence"""
        return self._style.from_ansi(ansiseq)

    @typing.overload
    def __getitem__(self, val: slice) -> list[S]: ...

    @typing.overload
    def __getitem__(self, val: tuple[int, int, int] | str | int) -> S: ...

    def __getitem__(self, val: slice | tuple[int, int, int] | str | int) -> list[S] | S:
        """\
        Shortcut to provide way to access colors numerically or by slice.
        If end <= 16, will stay to simple ANSI version."""
        if isinstance(val, slice):
            (start, stop, stride) = val.indices(256)
            if stop <= 16:
                return [self.simple(v) for v in range(start, stop, stride)]

            return [self.full(v) for v in range(start, stop, stride)]

        if isinstance(val, tuple):
            return self.rgb(*val)

        try:
            return self.full(val)
        except ColorNotFound:
            if isinstance(val, int):
                raise ColorNotFound(f"Color number {val!r} is out of range") from None
            return self.hex(val)

    def __call__(
        self,
        val_or_r: str | int | None = None,
        g: int | None = None,
        b: int | None = None,
    ) -> S:
        """Shortcut to provide way to access colors."""
        if val_or_r is None or (isinstance(val_or_r, str) and val_or_r == ""):
            return self._style()
        if isinstance(val_or_r, self._style):
            return self._style(val_or_r)
        if isinstance(val_or_r, str) and "\033" in val_or_r:
            return self.ansi(val_or_r)
        return self._style.from_color(
            self._style.color_class(val_or_r, g, b, fg=self._fg)
        )

    def __iter__(self) -> Iterator[S]:
        """Iterates through all colors in extended colorset."""
        return (self.full(i) for i in range(256))

    def __invert__(self) -> S:
        """Allows clearing a color with ~"""
        return self.reset

    def __enter__(self) -> Self:
        """This will reset the color on leaving the with statement."""
        return self

    def __exit__(self, _type: object, _value: object, _traceback: object) -> None:
        """This resets a FG/BG color or all styles,
        due to different definition of RESET for the
        factories."""

        self.reset.now()

    def __repr__(self) -> str:
        """Simple representation of the class by name."""
        return f"<{self.__class__.__name__}>"


class StyleFactory(ColorFactory[S]):
    """Factory for styles. Holds font styles, FG and BG objects representing colors, and
    imitates the FG ColorFactory to a large degree."""

    def __init__(self, style: type[S]):
        super().__init__(True, style)

        self.fg = ColorFactory(True, style)
        self.bg = ColorFactory(False, style)

        self.do_nothing = style()
        self.reset = style(reset=True)

        for item in style.attribute_names:
            setattr(self, item, style(attributes={item: True}))

        self.load_stylesheet(default_styles)

    @property
    def use_color(self) -> int:
        """Shortcut for setting color usage on Style"""
        return self._style.use_color

    @use_color.setter
    def use_color(self, val: int) -> None:
        self._style.use_color = val

    def from_ansi(self, ansi_sequence: str) -> S:
        """Calling this is a shortcut for creating a style from an ANSI sequence."""
        return self._style.from_ansi(ansi_sequence)

    def from_ansi_string(self, ansi_string: str) -> Iterator[S | str]:
        """Calling this is a shortcut for creating a style from an ANSI string.

        .. versionadded:: 2.0
        """
        return self._style.from_ansi_string(ansi_string)

    def sequence_to_string(self, sequence: Iterable[S | str]) -> str:
        """Converts a sequence of styles and strings into a single string with ANSI codes.

        .. versionadded:: 2.0
        """
        return self._style.sequence_to_string(sequence)

    @property
    def stdout(self) -> TextIO:
        """This is a shortcut for getting stdout from a class without an instance."""
        return self._style._stdout if self._style._stdout is not None else sys.stdout

    @stdout.setter
    def stdout(self, newout: TextIO) -> None:
        self._style._stdout = newout

    def get_colors_from_string(self, color: str = "") -> S:
        """
        Sets color based on string, use `.` or space for separator,
        and numbers, fg/bg, htmlcodes, etc all accepted (as strings).
        """

        names = color.replace(".", " ").split()
        # TODO: check the logic here to see what prev can be
        prev: Any = self
        styleslist: list[S] = []
        for name in names:
            try:
                prev = getattr(prev, name)
            except AttributeError:
                try:
                    prev = prev(int(name))
                except (ColorNotFound, ValueError):
                    prev = prev(name)
            if isinstance(prev, self._style):
                styleslist.append(prev)
                prev = self

        if styleslist:
            prev = functools.reduce(operator.and_, styleslist)

        return prev if isinstance(prev, self._style) else prev.reset

    def filter(self, colored_string: str) -> str:
        """Filters out colors in a string, returning only the name."""
        if isinstance(colored_string, self._style):
            return colored_string
        return self._style.string_filter_ansi(colored_string)

    def contains_colors(self, colored_string: str) -> bool:
        """Checks to see if a string contains colors."""
        return self._style.string_contains_colors(colored_string)

    def extract(self, colored_string: str) -> S:
        """Gets colors from an ansi string, returns those colors"""
        return self._style.from_ansi(colored_string, True)

    def load_stylesheet(self, stylesheet: Mapping[str, str] | None = None) -> None:
        if stylesheet is None:
            stylesheet = default_styles
        for item in stylesheet:
            setattr(self, item, self.get_colors_from_string(stylesheet[item]))
