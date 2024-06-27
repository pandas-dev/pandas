"""
This file provides two classes, `Color` and `Style`.

``Color`` is rarely used directly,
but merely provides the workhorse for finding and manipulating colors.

With the ``Style`` class, any color can be directly called or given to a with statement.
"""

import contextlib
import os
import platform
import re
import sys
from abc import ABCMeta, abstractmethod
from copy import copy
from typing import IO, Dict, Optional, Union

from .names import (
    FindNearest,
    attributes_ansi,
    color_codes_simple,
    color_html,
    color_names,
    from_html,
)

__all__ = [
    "Color",
    "Style",
    "ANSIStyle",
    "HTMLStyle",
    "ColorNotFound",
    "AttributeNotFound",
]

_lower_camel_names = [n.replace("_", "") for n in color_names]


def get_color_repr():
    """Gets best colors for current system."""
    if "NO_COLOR" in os.environ:
        return 0
    if os.environ.get("FORCE_COLOR", "") in {"0", "1", "2", "3", "4"}:
        return int(os.environ["FORCE_COLOR"])
    if not sys.stdout.isatty():
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
    Loaded with ``(r, g, b, fg)`` or ``(color, fg=fg)``. The second signature is a short cut
    and will try full and hex loading.

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

    __slots__ = ("fg", "isreset", "rgb", "number", "representation", "exact")

    def __init__(self, r_or_color=None, g=None, b=None, fg=True):
        """This works from color values, or tries to load non-simple ones."""

        if isinstance(r_or_color, type(self)):
            for item in ("fg", "isreset", "rgb", "number", "representation", "exact"):
                setattr(self, item, getattr(r_or_color, item))
            return

        self.fg = fg
        self.isreset = True  # Starts as reset color
        self.rgb = (0, 0, 0)

        self.number = None
        "Number of the original color, or closest color"

        self.representation = 4
        "0 for off, 1 for 8 colors, 2 for 16 colors, 3 for 256 colors, 4 for true color"

        self.exact = True
        "This is false if the named color does not match the real color"

        if None in (g, b):
            if not r_or_color:
                return
            try:
                self._from_simple(r_or_color)
            except ColorNotFound:
                try:
                    self._from_full(r_or_color)
                except ColorNotFound:
                    self._from_hex(r_or_color)

        elif None not in (r_or_color, g, b):
            self.rgb = (r_or_color, g, b)
            self._init_number()
        else:
            raise ColorNotFound("Invalid parameters for a color!")

    def _init_number(self):
        """Should always be called after filling in r, g, b, and representation.
        Color will not be a reset color anymore."""

        if self.representation in (0, 1):
            number = FindNearest(*self.rgb).only_basic()
        elif self.representation == 2:
            number = FindNearest(*self.rgb).only_simple()
        elif self.representation in (3, 4):
            number = FindNearest(*self.rgb).all_fast()

        if self.number is None:
            self.number = number

        self.isreset = False
        self.exact = self.rgb == from_html(color_html[self.number])
        if not self.exact:
            self.number = number

    @classmethod
    def from_simple(cls, color, fg=True):
        """Creates a color from simple name or color number"""
        self = cls(fg=fg)
        self._from_simple(color)
        return self

    def _from_simple(self, color):
        with contextlib.suppress(AttributeError):
            color = color.lower()
            color = color.replace(" ", "")
            color = color.replace("_", "")

        if color == "reset":
            return

        if color in _lower_camel_names[:16]:
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
    def from_full(cls, color, fg=True):
        """Creates a color from full name or color number"""
        self = cls(fg=fg)
        self._from_full(color)
        return self

    def _from_full(self, color):
        with contextlib.suppress(AttributeError):
            color = color.lower()
            color = color.replace(" ", "")
            color = color.replace("_", "")

        if color == "reset":
            return

        if color in _lower_camel_names:
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
    def from_hex(cls, color, fg=True):
        """Converts #123456 values to colors."""

        self = cls(fg=fg)
        self._from_hex(color)
        return self

    def _from_hex(self, color):
        try:
            self.rgb = from_html(color)
        except (TypeError, ValueError):
            raise ColorNotFound("Did not find htmlcode: " + repr(color)) from None

        self.representation = 4
        self._init_number()

    @property
    def name(self):
        """The (closest) name of the current color"""
        return "reset" if self.isreset else color_names[self.number]

    @property
    def name_camelcase(self):
        """The camelcase name of the color"""
        return self.name.replace("_", " ").title().replace(" ", "")

    def __repr__(self):
        """This class has a smart representation that shows name and color (if not unique)."""
        name = ["Deactivated:", " Basic:", "", " Full:", " True:"][self.representation]
        name += "" if self.fg else " Background"
        name += " " + self.name_camelcase
        name += "" if self.exact else " " + self.hex_code
        return name[1:]

    def __eq__(self, other):
        """Reset colors are equal, otherwise rgb have to match."""
        return other.isreset if self.isreset else self.rgb == other.rgb

    @property
    def ansi_sequence(self):
        """This is the ansi sequence as a string, ready to use."""
        return "\033[" + ";".join(map(str, self.ansi_codes)) + "m"

    @property
    def ansi_codes(self):
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
    def hex_code(self):
        """This is the hex code of the current color, html style notation."""

        return (
            "#000000"
            if self.isreset
            else f"#{self.rgb[0]:02X}{self.rgb[1]:02X}{self.rgb[2]:02X}"
        )

    def __str__(self):
        """This just prints it's simple name"""
        return self.name

    def to_representation(self, val):
        """Converts a color to any representation"""
        other = copy(self)
        other.representation = val
        if self.isreset:
            return other
        other.number = None
        other._init_number()
        return other

    def limit_representation(self, val):
        """Only converts if val is lower than representation"""

        return self if self.representation <= val else self.to_representation(val)


class Style(metaclass=ABCMeta):
    """This class allows the color changes to be called directly
    to write them to stdout, ``[]`` calls to wrap colors (or the ``.wrap`` method)
    and can be called in a with statement.
    """

    __slots__ = ("attributes", "fg", "bg", "isreset", "__weakref__")

    color_class = Color
    """The class of color to use. Never hardcode ``Color`` call when writing a Style
    method."""

    attribute_names: Union[Dict[str, str], Dict[str, int]]

    _stdout: Optional[IO] = None
    end = "\n"
    """The endline character. Override if needed in subclasses."""

    ANSI_REG = re.compile("\033\\[([\\d;]+)m")
    """The regular expression that finds ansi codes in a string."""

    @property
    def stdout(self):
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
    def stdout(self, newout):
        self.__class__._stdout = newout

    def __init__(self, attributes=None, fgcolor=None, bgcolor=None, reset=False):
        """This is usually initialized from a factory."""
        if isinstance(attributes, type(self)):
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
    def from_color(cls, color):
        return cls(fgcolor=color) if color.fg else cls(bgcolor=color)

    def invert(self):
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
            other.fg = self.fg.__class__()

        if self.bg:
            other.bg = self.bg.__class__()

        return other

    @property
    def reset(self):
        """Shortcut to access reset as a property."""
        return self.invert()

    def __copy__(self):
        """Copy is supported, will make dictionary and colors unique."""
        result = self.__class__()
        result.isreset = self.isreset
        result.fg = copy(self.fg)
        result.bg = copy(self.bg)
        result.attributes = copy(self.attributes)
        return result

    def __invert__(self):
        """This allows ~color."""
        return self.invert()

    def __add__(self, other):
        """Adding two matching Styles results in a new style with
        the combination of both. Adding with a string results in
        the string concatenation of a style.

        Addition is non-commutative, with the rightmost Style property
        being taken if both have the same property.
        (Not safe)"""
        if type(self) == type(other):
            result = copy(other)

            result.isreset = self.isreset or other.isreset
            for attribute in self.attributes:
                if attribute not in result.attributes:
                    result.attributes[attribute] = self.attributes[attribute]
            if not result.fg:
                result.fg = self.fg
            if not result.bg:
                result.bg = self.bg
            return result

        return other.__class__(self) + other

    def __radd__(self, other):
        """This only gets called if the string is on the left side. (Not safe)"""
        return other + other.__class__(self)

    def wrap(self, wrap_this):
        """Wrap a string in this style and its inverse."""
        return self + wrap_this + ~self

    def __and__(self, other):
        """This class supports ``color & color2`` syntax,
        and ``color & "String" syntax too.``"""
        if type(self) == type(other):
            return self + other

        return self.wrap(other)

    def __rand__(self, other):
        """This class supports ``"String:" & color`` syntax."""
        return self.wrap(other)

    def __ror__(self, other):
        """Support for "String" | color syntax"""
        return self.wrap(other)

    def __or__(self, other):
        """This class supports ``color | color2`` syntax. It also supports
        ``"color | "String"`` syntax too."""
        return self.__and__(other)

    def __call__(self):
        """\
        This is a shortcut to print color immediately to the stdout. (Not safe)
        """

        self.now()

    def now(self):
        """Immediately writes color to stdout. (Not safe)"""
        self.stdout.write(str(self))

    def print(self, *printables, **kargs):
        """\
        This acts like print; will print that argument to stdout wrapped
        in Style with the same syntax as the print function in 3.4."""

        end = kargs.get("end", self.end)
        sep = kargs.get("sep", " ")
        file = kargs.get("file", self.stdout)
        flush = kargs.get("flush", False)
        file.write(self.wrap(sep.join(map(str, printables))) + end)
        if flush:
            file.flush()

    print_ = print
    """DEPRECATED: Shortcut from classic Python 2"""

    def __getitem__(self, wrapped):
        """The [] syntax is supported for wrapping"""
        return self.wrap(wrapped)

    def __enter__(self):
        """Context manager support"""
        self.stdout.write(str(self))
        self.stdout.flush()

    def __exit__(self, _type, _value, _traceback):
        """Runs even if exception occurred, does not catch it."""
        self.stdout.write(str(~self))
        self.stdout.flush()
        return False

    @property
    def ansi_codes(self):
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
            self.bg.fg = False
            codes.extend(self.bg.ansi_codes)

        return codes

    @property
    def ansi_sequence(self):
        """This is the string ANSI sequence."""
        codes = ";".join(str(c) for c in self.ansi_codes)
        return f"\033[{codes}m" if codes else ""

    def __repr__(self):
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

    def __eq__(self, other):
        """Equality is true only if reset, or if attributes, fg, and bg match."""
        if type(self) == type(other):
            if self.isreset:
                return other.isreset

            return (
                self.attributes == other.attributes
                and self.fg == other.fg
                and self.bg == other.bg
            )

        return str(self) == other

    @abstractmethod
    def __str__(self):
        """Base Style does not implement a __str__ representation. This is the one
        required method of a subclass."""

    @classmethod
    def from_ansi(cls, ansi_string, filter_resets=False):
        """This generated a style from an ansi string. Will ignore resets if filter_resets is True."""
        result = cls()
        res = cls.ANSI_REG.search(ansi_string)
        for group in res.groups():
            sequence = map(int, group.split(";"))
            result.add_ansi(sequence, filter_resets)
        return result

    def add_ansi(self, sequence, filter_resets=False):
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
            return

    @classmethod
    def string_filter_ansi(cls, colored_string):
        """Filters out colors in a string, returning only the name."""
        return cls.ANSI_REG.sub("", colored_string)

    @classmethod
    def string_contains_colors(cls, colored_string):
        """Checks to see if a string contains colors."""
        return len(cls.ANSI_REG.findall(colored_string)) > 0

    def to_representation(self, rep):
        """This converts both colors to a specific representation"""
        other = copy(self)
        if other.fg:
            other.fg = other.fg.to_representation(rep)
        if other.bg:
            other.bg = other.bg.to_representation(rep)
        return other

    def limit_representation(self, rep):
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
    def basic(self):
        """The color in the 8 color representation."""
        return self.to_representation(1)

    @property
    def simple(self):
        """The color in the 16 color representation."""
        return self.to_representation(2)

    @property
    def full(self):
        """The color in the 256 color representation."""
        return self.to_representation(3)

    @property
    def true(self):
        """The color in the true color representation."""
        return self.to_representation(4)


class ANSIStyle(Style):
    """This is a subclass for ANSI styles. Use it to get
    color on sys.stdout tty terminals on posix systems.

    Set ``use_color = True/False`` if you want to control color
    for anything using this Style."""

    __slots__ = ()
    use_color = get_color_repr()

    attribute_names = attributes_ansi

    def __str__(self):
        return (
            self.limit_representation(self.use_color).ansi_sequence
            if self.use_color
            else ""
        )


class HTMLStyle(Style):
    """This was meant to be a demo of subclassing Style, but
    actually can be a handy way to quickly color html text."""

    __slots__ = ()
    attribute_names = {
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

    def __str__(self):
        if self.isreset:
            raise ResetNotSupported("HTML does not support global resets!")

        result = ""

        if self.bg and not self.bg.isreset:
            result += f'<span style="background-color: {self.bg.hex_code}">'
        if self.fg and not self.fg.isreset:
            result += f'<font color="{self.fg.hex_code}">'
        for attr in sorted(self.attributes):
            if self.attributes[attr]:
                result += "<" + self.attribute_names[attr] + ">"

        for attr in sorted(self.attributes, reverse=True):
            if not self.attributes[attr]:
                result += "</" + self.attribute_names[attr].split(" ")[0] + ">"
        if self.fg and self.fg.isreset:
            result += "</font>"
        if self.bg and self.bg.isreset:
            result += "</span>"

        return result
