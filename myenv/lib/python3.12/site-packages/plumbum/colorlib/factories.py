"""
Color-related factories. They produce Styles.
"""

import functools
import operator
import sys
from typing import Any

from .names import color_names, default_styles
from .styles import ColorNotFound

__all__ = ["ColorFactory", "StyleFactory"]


class ColorFactory:
    """This creates color names given fg = True/False. It usually will
    be called as part of a StyleFactory."""

    def __init__(self, fg, style):
        self._fg = fg
        self._style = style
        self.reset = style.from_color(style.color_class(fg=fg))

        # Adding the color name shortcuts for foreground colors
        for item in color_names[:16]:
            setattr(
                self, item, style.from_color(style.color_class.from_simple(item, fg=fg))
            )

    def __getattr__(self, item):
        """Full color names work, but do not populate __dir__."""
        try:
            return self._style.from_color(self._style.color_class(item, fg=self._fg))
        except ColorNotFound:
            raise AttributeError(item) from None

    def full(self, name):
        """Gets the style for a color, using standard name procedure: either full
        color name, html code, or number."""
        return self._style.from_color(
            self._style.color_class.from_full(name, fg=self._fg)
        )

    def simple(self, name):
        """Return the extended color scheme color for a value or name."""
        return self._style.from_color(
            self._style.color_class.from_simple(name, fg=self._fg)
        )

    def rgb(self, r, g=None, b=None):
        """Return the extended color scheme color for a value."""
        if g is None and b is None:
            return self.hex(r)

        return self._style.from_color(self._style.color_class(r, g, b, fg=self._fg))

    def hex(self, hexcode):
        """Return the extended color scheme color for a value."""
        return self._style.from_color(
            self._style.color_class.from_hex(hexcode, fg=self._fg)
        )

    def ansi(self, ansiseq):
        """Make a style from an ansi text sequence"""
        return self._style.from_ansi(ansiseq)

    def __getitem__(self, val):
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
            return self.hex(val)

    def __call__(self, val_or_r=None, g=None, b=None):
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

    def __iter__(self):
        """Iterates through all colors in extended colorset."""
        return (self.full(i) for i in range(256))

    def __invert__(self):
        """Allows clearing a color with ~"""
        return self.reset

    def __enter__(self):
        """This will reset the color on leaving the with statement."""
        return self

    def __exit__(self, _type: Any, _value: Any, _traceback: Any) -> None:
        """This resets a FG/BG color or all styles,
        due to different definition of RESET for the
        factories."""

        self.reset.now()

    def __repr__(self):
        """Simple representation of the class by name."""
        return f"<{self.__class__.__name__}>"


class StyleFactory(ColorFactory):
    """Factory for styles. Holds font styles, FG and BG objects representing colors, and
    imitates the FG ColorFactory to a large degree."""

    def __init__(self, style):
        super().__init__(True, style)

        self.fg = ColorFactory(True, style)
        self.bg = ColorFactory(False, style)

        self.do_nothing = style()
        self.reset = style(reset=True)

        for item in style.attribute_names:
            setattr(self, item, style(attributes={item: True}))

        self.load_stylesheet(default_styles)

    @property
    def use_color(self):
        """Shortcut for setting color usage on Style"""
        return self._style.use_color

    @use_color.setter
    def use_color(self, val):
        self._style.use_color = val

    def from_ansi(self, ansi_sequence):
        """Calling this is a shortcut for creating a style from an ANSI sequence."""
        return self._style.from_ansi(ansi_sequence)

    @property
    def stdout(self):
        """This is a shortcut for getting stdout from a class without an instance."""
        return self._style._stdout if self._style._stdout is not None else sys.stdout

    @stdout.setter
    def stdout(self, newout):
        self._style._stdout = newout

    def get_colors_from_string(self, color=""):
        """
        Sets color based on string, use `.` or space for separator,
        and numbers, fg/bg, htmlcodes, etc all accepted (as strings).
        """

        names = color.replace(".", " ").split()
        prev = self
        styleslist = []
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

    def filter(self, colored_string):
        """Filters out colors in a string, returning only the name."""
        if isinstance(colored_string, self._style):
            return colored_string
        return self._style.string_filter_ansi(colored_string)

    def contains_colors(self, colored_string):
        """Checks to see if a string contains colors."""
        return self._style.string_contains_colors(colored_string)

    def extract(self, colored_string):
        """Gets colors from an ansi string, returns those colors"""
        return self._style.from_ansi(colored_string, True)

    def load_stylesheet(self, stylesheet=None):
        if stylesheet is None:
            stylesheet = default_styles
        for item in stylesheet:
            setattr(self, item, self.get_colors_from_string(stylesheet[item]))
