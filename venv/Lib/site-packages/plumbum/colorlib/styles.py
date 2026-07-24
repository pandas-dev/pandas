"""
This file provides two classes, `Color` and `Style`.

``Color`` is rarely used directly,
but merely provides the workhorse for finding and manipulating colors.

With the ``Style`` class, any color can be directly called or given to a with statement.
"""

from __future__ import annotations

__lazy_modules__ = {"contextlib", "copy", "platform", "re", "warnings"}

import contextlib
import dataclasses
import os
import platform
import re
import sys
import typing
import warnings
from abc import ABCMeta, abstractmethod
from copy import copy
from typing import Any, ClassVar, Literal, TextIO, TypeVar

from .names import (
    FindNearest,
    attributes_ansi,
    color_codes_simple,
    color_html,
    color_names,
    from_html,
)

if typing.TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from plumbum._compat.typing import Self

__all__ = [
    "ANSIStyle",
    "AttributeNotFound",
    "Color",
    "ColorNotFound",
    "HTMLStyle",
    "Style",
]

_lower_camel_names = [n.replace("_", "") for n in color_names]


def get_color_repr() -> int:
    """Gets best colors for current system."""
    if "NO_COLOR" in os.environ:
        return 0
    if os.environ.get("FORCE_COLOR", "") in {"0", "1", "2", "3", "4"}:
        return int(os.environ["FORCE_COLOR"])
    if (
        sys.stdout is None
        or not hasattr(sys.stdout, "isatty")
        or not sys.stdout.isatty()
    ):
        return 0

    term = os.environ.get("TERM", "")

    # Some terminals set TERM=xterm for compatibility
    if term.endswith("256color") or term == "xterm":
        return 3 if platform.system() == "Darwin" else 4
    if term.endswith("16color"):
        return 2
    if term == "screen":
        return 1
    if os.name == "nt":
        return 0

    return 3


class ColorNotFound(Exception):
    """Thrown when a color is not valid for a particular method."""


class AttributeNotFound(Exception):
    """Similar to color not found, only for attributes."""


class ResetNotSupported(Exception):
    """An exception indicating that Reset is not available
    for this Style."""


class Color:
    """\
    Loaded with ``(r, g, b, fg=fg)`` or ``(color, fg=fg)``. The second
    signature is a shortcut and will try full and hex loading.

    This class stores the idea of a color, rather than a specific implementation.
    It provides as many different tools for representations as possible, and can be subclassed
    to add more representations, though that should not be needed for most situations. ``.from_`` class methods provide quick ways to create colors given different representations.
    You will not usually interact with this class.

    Possible colors::

        reset = Color() # The reset color by default
        background_reset = Color(fg=False) # Can be a background color
        blue = Color(0,0,255) # Red, Green, Blue
        green = Color.from_full("green") # Case insensitive name, from large colorset
        red = Color.from_full(1) # Color number
        white = Color.from_html("#FFFFFF") # HTML supported
        yellow = Color.from_simple("red") # Simple colorset


    The attributes are:

    .. data:: reset

        True it this is a reset color (following attributes don't matter if True)

    .. data:: rgb

        The red/green/blue tuple for this color

    .. data:: simple

        If true will stay to 16 color mode.

    .. data:: number

        The color number given the mode, closest to rgb
        if not rgb not exact, gives position of closest name.

    .. data:: fg

        This is a foreground color if True. Background color if False.

        """

    __slots__ = ("exact", "fg", "isreset", "number", "representation", "rgb")

    def __init__(
        self,
        r_or_color: None | int | Color | str = None,
        g: int | None = None,
        b: int | None = None,
        /,
        *,
        fg: bool = True,
    ):
        """This works from color values, or tries to load non-simple ones."""

        if isinstance(r_or_color, Color):
            for item in ("fg", "isreset", "rgb", "number", "representation", "exact"):
                setattr(self, item, getattr(r_or_color, item))
            return

        self.fg = fg
        self.isreset = True  # Starts as reset color
        self.rgb = (0, 0, 0)

        # This will get set to a real int below!
        self.number = typing.cast("int", None)
        "Number of the original color, or closest color"

        self.representation = 2
        "0 for off, 1 for 8 colors, 2 for 16 colors, 3 for 256 colors, 4 for true color"

        self.exact = True
        "This is false if the named color does not match the real color"

        if (g is None) != (b is None):
            raise ColorNotFound(
                "Color requires either no g/b arguments or both g and b together"
            )

        if g is None and b is None:
            if r_or_color is None or r_or_color == "":
                return
            try:
                self._from_simple(r_or_color)
            except ColorNotFound:
                try:
                    self._from_full(r_or_color)
                except ColorNotFound:
                    if not isinstance(r_or_color, str):
                        raise ColorNotFound(
                            "Did not find color: " + repr(r_or_color)
                        ) from None
                    self._from_hex(r_or_color)

        elif isinstance(r_or_color, int) and g is not None and b is not None:
            self.rgb = (r_or_color, g, b)
            self.representation = 4
            self._init_number()
        else:
            raise ColorNotFound("Invalid parameters for a color!")

    def _init_number(self) -> None:
        """Should always be called after filling in r, g, b, and representation.
        Color will not be a reset color anymore."""
        # TODO: make this a function that returns the number, instead of mutation

        if self.representation in (0, 1):
            number = FindNearest(*self.rgb).only_basic()
        elif self.representation == 2:
            number = FindNearest(*self.rgb).only_simple()
        elif self.representation in (3, 4):
            number = FindNearest(*self.rgb).all_fast()
        else:
            raise AssertionError("Invalid representation, needs to be 0-4")

        if self.number is None:
            self.number = number

        self.isreset = False
        self.exact = self.rgb == from_html(color_html[self.number])
        if not self.exact:
            self.number = number

    @classmethod
    def from_simple(cls, color: str | int, fg: bool = True) -> Self:
        """Creates a color from simple name or color number"""
        self = cls(fg=fg)
        self._from_simple(color)
        return self

    def _from_simple(self, color: str | int) -> None:
        if isinstance(color, str):
            color = color.lower()
            color = color.replace(" ", "")
            color = color.replace("_", "")

        if color == "reset":
            return

        if isinstance(color, str) and color in _lower_camel_names[:16]:
            self.number = _lower_camel_names.index(color)
            self.rgb = from_html(color_html[self.number])

        elif isinstance(color, int) and 0 <= color < 16:
            self.number = color
            self.rgb = from_html(color_html[color])

        else:
            raise ColorNotFound("Did not find color: " + repr(color))

        self.representation = 2
        self._init_number()

    @classmethod
    def from_full(cls, color: str | int, fg: bool = True) -> Self:
        """Creates a color from full name or color number"""
        self = cls(fg=fg)
        self._from_full(color)
        return self

    def _from_full(self, color: str | int) -> None:
        if isinstance(color, str):
            color = color.lower()
            color = color.replace(" ", "")
            color = color.replace("_", "")

        if color == "reset":
            self.representation = 2
            return

        if isinstance(color, str) and color in _lower_camel_names:
            self.number = _lower_camel_names.index(color)
            self.rgb = from_html(color_html[self.number])

        elif isinstance(color, int) and 0 <= color <= 255:
            self.number = color
            self.rgb = from_html(color_html[color])

        else:
            raise ColorNotFound("Did not find color: " + repr(color))

        self.representation = 3
        self._init_number()

    @classmethod
    def from_hex(cls, color: str, fg: bool = True) -> Self:
        """Converts #123456 values to colors."""

        self = cls(fg=fg)
        self._from_hex(color)
        return self

    def _from_hex(self, color: str) -> None:
        try:
            self.rgb = from_html(color)
        except (TypeError, ValueError):
            raise ColorNotFound("Did not find htmlcode: " + repr(color)) from None

        self.representation = 4
        self._init_number()

    @property
    def name(self) -> str:
        """The (closest) name of the current color"""
        return "reset" if self.isreset else color_names[self.number]

    @property
    def name_camelcase(self) -> str:
        """The camelcase name of the color"""
        return self.name.replace("_", " ").title().replace(" ", "")

    def __repr__(self) -> str:
        """This class has a smart representation that shows name and color (if not unique)."""
        name = ["Deactivated: ", "BasicColor: ", "", "FullColor: ", "TrueColor: "][
            self.representation
        ]
        name += "Foreground" if self.fg else "Background"
        name += " " + self.name_camelcase
        name += "" if self.exact else " " + self.hex_code
        return name

    def __eq__(self, other: object) -> bool:
        """Reset colors are equal, otherwise rgb have to match."""
        if not isinstance(other, Color):
            return NotImplemented
        return other.isreset if self.isreset else self.rgb == other.rgb

    def __hash__(self) -> int:
        return hash(self.isreset or self.rgb)

    @property
    def ansi_sequence(self) -> str:
        """This is the ansi sequence as a string, ready to use."""
        return "\033[" + ";".join(map(str, self.ansi_codes)) + "m"

    @property
    def ansi_codes(
        self,
    ) -> tuple[int] | tuple[int, int, int] | tuple[int, int, int, int, int]:
        """This is the full ANSI code, can be reset, simple, 256, or full color."""
        ansi_addition = 30 if self.fg else 40

        if self.isreset:
            return (ansi_addition + 9,)
        if self.representation < 3:
            return (color_codes_simple[self.number] + ansi_addition,)
        if self.representation == 3:
            return (ansi_addition + 8, 5, self.number)

        return (ansi_addition + 8, 2, self.rgb[0], self.rgb[1], self.rgb[2])

    @property
    def hex_code(self) -> str:
        """This is the hex code of the current color, html style notation."""

        return (
            "#000000"
            if self.isreset
            else f"#{self.rgb[0]:02X}{self.rgb[1]:02X}{self.rgb[2]:02X}"
        )

    def __str__(self) -> str:
        """This just prints it's simple name"""
        return self.name

    def to_representation(self, val: int) -> Self:
        """Converts a color to any representation"""
        other = copy(self)
        other.representation = val
        if self.isreset:
            return other
        other.number = typing.cast("int", None)
        other._init_number()
        return other

    def limit_representation(self, val: int) -> Self:
        """Only converts if val is lower than representation"""

        return self if self.representation <= val else self.to_representation(val)

    def reset(self) -> Self:
        """This is a shortcut to get a reset color, keeps the same representation."""
        retval = self.__class__(fg=self.fg)
        retval.representation = self.representation
        return retval


S = TypeVar("S", "Style", str)


class Style(metaclass=ABCMeta):
    """This class allows the color changes to be called directly
    to write them to stdout, ``[]`` calls to wrap colors (or the ``.wrap`` method)
    and can be called in a with statement.
    """

    __slots__ = ("__weakref__", "attributes", "bg", "fg", "isreset")

    color_class = Color
    """The class of color to use. Never hardcode ``Color`` call when writing a Style
    method."""

    # These must be defined by subclasses
    # pylint: disable-next=declare-non-slot
    attribute_names: ClassVar[dict[str, str] | dict[str, int]]

    _stdout: ClassVar[TextIO | None] = None
    end = "\n"
    """The endline character. Override if needed in subclasses."""

    ANSI_REG = re.compile("\033\\[([\\d;]+)m")
    """The regular expression that finds ansi codes in a string."""

    use_color = 4
    """The color level. This is a default value that subclasses may override. For example, ANSIStyle uses this to control color output, while other styles may always use the maximum level (4)."""

    @property
    def stdout(self) -> TextIO:
        """\
        This property will allow custom, class level control of stdout.
        It will use current sys.stdout if set to None (default).
        Unfortunately, it only works on an instance..
        """
        # Import sys repeated here to make calling this stable in atexit function
        import sys  # pylint: disable=reimported, redefined-outer-name, import-outside-toplevel

        return (
            self.__class__._stdout if self.__class__._stdout is not None else sys.stdout
        )

    @stdout.setter
    def stdout(self, newout: TextIO) -> None:
        self.__class__._stdout = newout

    def __init__(
        self,
        attributes: Style | dict[str, bool] | None = None,
        fgcolor: Color | None = None,
        bgcolor: Color | None = None,
        reset: bool = False,
    ):
        """This is usually initialized from a factory."""
        if isinstance(attributes, Style):
            for item in ("attributes", "fg", "bg", "isreset"):
                setattr(self, item, copy(getattr(attributes, item)))
            return
        self.attributes = attributes if attributes is not None else {}
        self.fg = fgcolor
        self.bg = bgcolor
        self.isreset = reset
        invalid_attributes = set(self.attributes) - set(self.attribute_names)
        if len(invalid_attributes) > 0:
            raise AttributeNotFound(
                "Attribute(s) not valid: " + ", ".join(invalid_attributes)
            )

    @classmethod
    def from_color(cls, color: Color) -> Self:
        return cls(fgcolor=color) if color.fg else cls(bgcolor=color)

    def invert(self) -> Self:
        """This resets current color(s) and flips the value of all
        attributes present"""

        other = self.__class__()

        # Opposite of reset is reset
        if self.isreset:
            other.isreset = True
            return other

        # Flip all attributes
        for attribute in self.attributes:
            other.attributes[attribute] = not self.attributes[attribute]

        # Reset only if color present
        if self.fg:
            other.fg = self.fg.reset()

        if self.bg:
            other.bg = self.bg.reset()

        return other

    @property
    def reset(self) -> Self:
        """Shortcut to access reset as a property."""
        return self.invert()

    def __copy__(self) -> Self:
        """Copy is supported, will make dictionary and colors unique."""
        result = self.__class__()
        result.isreset = self.isreset
        result.fg = copy(self.fg)
        result.bg = copy(self.bg)
        result.attributes = copy(self.attributes)
        return result

    def __invert__(self) -> Self:
        """This allows ~color."""
        return self.invert()

    def __add__(self, other: S) -> S:
        """Adding two matching Styles results in a new style with
        the combination of both. Adding with a string results in
        the string concatenation of a style.

        Addition is non-commutative, with the rightmost Style property
        being taken if both have the same property.
        (Not safe)"""
        if isinstance(other, Style):
            result = copy(other)

            result.isreset = self.isreset or other.isreset
            for attribute in self.attributes:
                if attribute not in result.attributes:
                    result.attributes[attribute] = self.attributes[attribute]
            if not result.fg:
                result.fg = self.fg
            elif result.fg.isreset and self.fg:
                result.fg.representation = self.fg.representation
            if not result.bg:
                result.bg = self.bg
            elif result.bg.isreset and self.bg:
                result.bg.representation = self.bg.representation
            return result

        return other.__class__(self) + other

    def __radd__(self, other: str) -> str:
        """This only gets called if the string is on the left side. (Not safe)"""
        return other + other.__class__(self)

    def _clean_resets(self) -> Self:
        """This removes resets from the style"""
        other = copy(self)
        if other.fg and other.fg.isreset:
            other.fg = None
        if other.bg and other.bg.isreset:
            other.bg = None
        for attribute in list(other.attributes):
            if not other.attributes[attribute]:
                del other.attributes[attribute]
        return other

    def wrap(self, wrap_this: str) -> str:
        """Wrap a string in this style and its inverse."""
        return self + wrap_this + ~self

    def __and__(self, other: S) -> S:
        """This class supports ``color & color2`` syntax,
        and ``color & "String" syntax too.``"""
        if isinstance(other, Style):
            return self + other

        return self.wrap(other)

    def __rand__(self, other: str) -> str:
        """This class supports ``"String:" & color`` syntax."""
        return self.wrap(other)

    def __ror__(self, other: str) -> str:
        """Support for "String" | color syntax"""
        return self.wrap(other)

    def __or__(self, other: S) -> S:
        """This class supports ``color | color2`` syntax. It also supports
        ``"color | "String"`` syntax too."""
        return self.__and__(other)

    def __call__(self) -> None:
        """\
        This is a shortcut to print color immediately to the stdout. (Not safe)
        """

        self.now()

    def now(self) -> None:
        """Immediately writes color to stdout. (Not safe)"""
        # Silently handle broken pipe (e.g., when output is piped to head)
        # OSError (errno 22, EINVAL) is raised on Windows when the pipe is closed
        with contextlib.suppress(BrokenPipeError, OSError):
            self.stdout.write(str(self))

    def print(self, *printables: object, **kargs: Any) -> None:
        """\
        This acts like print; will print that argument to stdout wrapped
        in Style with the same syntax as the print function in 3.4."""

        end = kargs.get("end", self.end)
        sep = kargs.get("sep", " ")
        file = kargs.get("file", self.stdout)
        flush = kargs.get("flush", False)
        try:
            file.write(self.wrap(sep.join(map(str, printables))) + end)
            if flush:
                file.flush()
        except (BrokenPipeError, OSError):
            # Silently handle broken pipe (e.g., when output is piped to head)
            # OSError (errno 22, EINVAL) is raised on Windows when the pipe is closed
            pass

    def print_(self, *printables: object, **kargs: Any) -> None:
        """DEPRECATED: Shortcut from classic Python 2; use ``.print`` instead."""
        warnings.warn("Use .print instead", FutureWarning, stacklevel=2)
        self.print(*printables, **kargs)

    def __getitem__(self, wrapped: str) -> str:
        """The [] syntax is supported for wrapping"""
        return self.wrap(wrapped)

    def __enter__(self) -> None:
        """Context manager support"""
        try:
            self.stdout.write(str(self))
            self.stdout.flush()
        except (BrokenPipeError, OSError):
            # Silently handle broken pipe (e.g., when output is piped to head)
            # OSError (errno 22, EINVAL) is raised on Windows when the pipe is closed
            pass

    def __exit__(
        self, _type: object, _value: object, _traceback: object
    ) -> Literal[False]:
        """Runs even if exception occurred, does not catch it."""
        try:
            self.stdout.write(str(~self))
            self.stdout.flush()
        except (BrokenPipeError, OSError):
            # Silently handle broken pipe (e.g., when output is piped to head)
            # OSError (errno 22, EINVAL) is raised on Windows when the pipe is closed
            pass
        return False

    @property
    def ansi_codes(self) -> list[int]:
        """Generates the full ANSI code sequence for a Style"""

        if self.isreset:
            return [0]

        codes = []
        for attribute in self.attributes:
            if self.attributes[attribute]:
                codes.append(attributes_ansi[attribute])
            else:
                # Fixing bold inverse being 22 instead of 21 on some terminals:
                codes.append(
                    attributes_ansi[attribute] + 20
                    if attributes_ansi[attribute] != 1
                    else 22
                )

        if self.fg:
            codes.extend(self.fg.ansi_codes)

        if self.bg:
            bg = copy(self.bg)
            bg.fg = False
            codes.extend(bg.ansi_codes)

        return codes

    @property
    def ansi_sequence(self) -> str:
        """This is the string ANSI sequence."""
        codes = ";".join(str(c) for c in self.ansi_codes)
        return f"\033[{codes}m" if codes else ""

    def __repr__(self) -> str:
        name = self.__class__.__name__
        attributes = ", ".join(a for a in self.attributes if self.attributes[a])
        neg_attributes = ", ".join(
            f"-{a}" for a in self.attributes if not self.attributes[a]
        )
        colors = ", ".join(repr(c) for c in (self.fg, self.bg) if c)
        string = (
            "; ".join(s for s in (attributes, neg_attributes, colors) if s) or "empty"
        )
        if self.isreset:
            string = "reset"
        return f"<{name}: {string}>"

    def __eq__(self, other: object) -> bool:
        """Equality is true only if reset, or if attributes, fg, and bg match."""
        if isinstance(other, Style):
            if self.isreset:
                return other.isreset

            return (
                self.attributes == other.attributes
                and self.fg == other.fg
                and self.bg == other.bg
            )

        return str(self) == other

    __hash__ = None  # type: ignore[assignment]

    @abstractmethod
    def __str__(self) -> str:
        """Base Style does not implement a __str__ representation. This is the one
        required method of a subclass."""

    @classmethod
    def from_ansi(cls, ansi_string: str, filter_resets: bool = False) -> Self:
        """This generated a style from an ansi string. Will ignore resets if filter_resets is True."""
        result = cls()
        for res in cls.ANSI_REG.finditer(ansi_string):
            for group in res.groups():
                sequence = map(int, group.split(";"))
                result.add_ansi(sequence, filter_resets)
        return result

    @classmethod
    def from_ansi_string(cls, ansi_string: str) -> Iterator[str | Self]:
        """This will read in a string with interleaved codes

        .. versionadded:: 2.0
        """
        last_end = 0
        for res in cls.ANSI_REG.finditer(ansi_string):
            if res.start() > last_end:
                yield ansi_string[last_end : res.start()]
            style = cls()
            for group in res.groups():
                sequence = map(int, group.split(";"))
                style.add_ansi(sequence)
            yield style
            last_end = res.end()
        if last_end < len(ansi_string):
            yield ansi_string[last_end:]

    @classmethod
    def sequence_to_string(cls, sequence: Iterable[Style | str]) -> str:
        """This is the opposite of from_ansi_string, it takes a list of styles and strings and concatenates them.

        .. versionadded:: 2.0
        """
        return "".join(str(s) for s in sequence)

    def add_ansi(self, sequence: Iterable[int], filter_resets: bool = False) -> None:
        """Adds a sequence of ansi numbers to the class. Will ignore resets if filter_resets is True."""

        values = iter(sequence)
        try:
            while True:
                value = next(values)
                if value in {38, 48}:
                    fg = value == 38
                    value = next(values)
                    if value == 5:
                        value = next(values)
                        if fg:
                            self.fg = self.color_class.from_full(value)
                        else:
                            self.bg = self.color_class.from_full(value, fg=False)
                    elif value == 2:
                        r = next(values)
                        g = next(values)
                        b = next(values)
                        if fg:
                            self.fg = self.color_class(r, g, b)
                        else:
                            self.bg = self.color_class(r, g, b, fg=False)
                    else:
                        raise ColorNotFound("the value 5 or 2 should follow a 38 or 48")
                elif value == 0:
                    if filter_resets is False:
                        self.isreset = True
                elif value in attributes_ansi.values():
                    for name, att_value in attributes_ansi.items():
                        if value == att_value:
                            self.attributes[name] = True
                elif value == 22:
                    # ANSI 22 = "normal intensity": clears both bold (1) and dim (2).
                    # ansi_codes emits 22 for bold=False (not 21) to match terminal
                    # behaviour, so we must parse 22 as bold=False here to round-trip.
                    if filter_resets is False:
                        if "bold" in self.attribute_names:
                            self.attributes["bold"] = False
                        if "dim" in self.attribute_names:
                            self.attributes["dim"] = False
                elif value in (20 + n for n in attributes_ansi.values()):
                    if filter_resets is False:
                        for name, att_value in attributes_ansi.items():
                            if value == att_value + 20:
                                self.attributes[name] = False
                elif 30 <= value <= 37:
                    self.fg = self.color_class.from_simple(value - 30)
                elif 40 <= value <= 47:
                    self.bg = self.color_class.from_simple(value - 40, fg=False)
                elif 90 <= value <= 97:
                    self.fg = self.color_class.from_simple(value - 90 + 8)
                elif 100 <= value <= 107:
                    self.bg = self.color_class.from_simple(value - 100 + 8, fg=False)
                elif value == 39:
                    if filter_resets is False:
                        self.fg = self.color_class()
                elif value == 49:
                    if filter_resets is False:
                        self.bg = self.color_class(fg=False)
                else:
                    raise ColorNotFound(f"The code {value} is not recognised")
        except StopIteration:
            pass

    @classmethod
    def string_filter_ansi(cls, colored_string: str) -> str:
        """Filters out colors in a string, returning only the name."""
        return cls.ANSI_REG.sub("", colored_string)

    @classmethod
    def string_contains_colors(cls, colored_string: str) -> bool:
        """Checks to see if a string contains colors."""
        return len(cls.ANSI_REG.findall(colored_string)) > 0

    def to_representation(self, rep: int) -> Self:
        """This converts both colors to a specific representation"""
        other = copy(self)
        if other.fg:
            other.fg = other.fg.to_representation(rep)
        if other.bg:
            other.bg = other.bg.to_representation(rep)
        return other

    def limit_representation(self, rep: int) -> Self:
        """This only converts if true representation is higher"""

        if rep is True or rep is False:
            return self

        other = copy(self)
        if other.fg:
            other.fg = other.fg.limit_representation(rep)
        if other.bg:
            other.bg = other.bg.limit_representation(rep)
        return other

    @property
    def basic(self) -> Self:
        """The color in the 8 color representation."""
        return self.to_representation(1)

    @property
    def simple(self) -> Self:
        """The color in the 16 color representation."""
        return self.to_representation(2)

    @property
    def full(self) -> Self:
        """The color in the 256 color representation."""
        return self.to_representation(3)

    @property
    def true(self) -> Self:
        """The color in the true color representation."""
        return self.to_representation(4)


class ANSIStyle(Style):
    """This is a subclass for ANSI styles. Use it to get
    color on sys.stdout tty terminals on posix systems.

    Set ``use_color = 0 # up to 4`` if you want to control color
    for anything using this Style."""

    __slots__ = ()
    use_color = get_color_repr()

    attribute_names = attributes_ansi

    def __str__(self) -> str:
        return (
            self.limit_representation(self.use_color).ansi_sequence
            if self.use_color
            else ""
        )


@dataclasses.dataclass(init=False, eq=True)
class _HTMLTag:
    """This is a helper for HTMLStyle"""

    __slots__ = ("attributes", "closing", "name")
    name: str
    attributes: dict[str, str]
    closing: bool

    def __init__(self, name: str, closing: bool = False, /, **attributes: str) -> None:
        self.name = name
        self.attributes = attributes
        self.closing = closing

    @classmethod
    def from_attr_name(cls, tag_string: str) -> _HTMLTag:
        name, _, attr = tag_string.partition(" ")
        if not attr:
            return cls(name)

        # We don't join attributes, so only supporting one is fine for now.
        attr_name, _, attr_value = attr.partition("=")
        if attr_value.startswith('"') and attr_value.endswith('"'):
            attr_value = attr_value[1:-1]
        return cls(name, **{attr_name: attr_value})

    def __str__(self) -> str:
        if self.closing:
            return f"</{self.name}>"
        if self.attributes:
            attr_string = " ".join(f'{k}="{v}"' for k, v in self.attributes.items())
            return f"<{self.name} {attr_string}>"
        return f"<{self.name}>"

    def close(self) -> Self:
        return self.__class__(self.name, True, **self.attributes)


class HTMLStyle(Style):
    """This was meant to be a demo of subclassing Style, but
    actually can be a handy way to quickly color html text."""

    __slots__ = ()
    attribute_names: ClassVar[dict[str, str]] = {
        "bold": "b",
        "em": "em",
        "italics": "i",
        "li": "li",
        "underline": 'span style="text-decoration: underline;"',
        "code": "code",
        "ol": "ol start=0",
        "strikeout": "s",
    }
    end = "<br/>\n"

    @classmethod
    def sequence_to_string(cls, sequence: Iterable[Style | str]) -> str:
        return "".join(cls._sequence_to_sequence(sequence))

    @classmethod
    def _sequence_to_sequence(cls, sequence: Iterable[Style | str]) -> Iterator[str]:
        """Convert sequence of styles and strings to HTML, handling tag closing order"""
        current_tags: list[_HTMLTag] = []
        current_style = cls()

        for item in sequence:
            if isinstance(item, str):
                yield item
            elif not cls.use_color:
                pass
            elif item.isreset:
                for tag in reversed(current_tags):
                    yield str(tag.close())
                current_tags = []
                current_style = cls()
            else:
                current_style += item  # type: ignore[assignment]
                current_style = current_style._clean_resets()
                new_tags = list(current_style._get_tags())

                # Find common prefix between current and new tags
                common_prefix_len = 0
                for i, (old, new) in enumerate(zip(current_tags, new_tags)):
                    if old == new:
                        common_prefix_len = i + 1
                    else:
                        break

                # Close tags that are different (in reverse order)
                for tag in reversed(current_tags[common_prefix_len:]):
                    yield str(tag.close())

                # Open new tags that are different
                for tag in new_tags[common_prefix_len:]:
                    yield str(tag)

                current_tags = new_tags

        for tag in reversed(current_tags):
            yield str(tag.close())

    def _get_tags(self) -> Iterator[_HTMLTag]:
        # Don't open tags if we're closing fg/bg
        fg_resetting = self.fg and self.fg.isreset
        bg_resetting = self.bg and self.bg.isreset

        if self.bg and not self.bg.isreset:
            yield _HTMLTag("span", style=f"background-color: {self.bg.hex_code}")
        if self.fg and not self.fg.isreset:
            yield _HTMLTag("font", color=self.fg.hex_code)

        # Only open attribute tags if we're not resetting
        if not (fg_resetting or bg_resetting):
            for attr in sorted(self.attributes):
                if self.attributes[attr]:
                    yield _HTMLTag.from_attr_name(self.attribute_names[attr])

        # Closing instead
        for attr in sorted(self.attributes, reverse=True):
            if not self.attributes[attr]:
                yield _HTMLTag.from_attr_name(self.attribute_names[attr]).close()
        # Close all open attributes before closing the foreground tag
        if fg_resetting:
            yield _HTMLTag("font").close()
        if bg_resetting:
            yield _HTMLTag("span").close()

    def __str__(self) -> str:
        if not self.use_color:
            return ""
        if self.isreset:
            raise ResetNotSupported("HTML does not support global resets!")

        return "".join(str(tag) for tag in self._get_tags())

    def __eq__(self, other: object) -> bool:
        if self.isreset:
            return isinstance(other, HTMLStyle) and other.isreset
        return super().__eq__(other)

    __hash__ = None


def __dir__() -> list[str]:
    return list(__all__)
