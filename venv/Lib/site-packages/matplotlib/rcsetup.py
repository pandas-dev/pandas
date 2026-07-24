"""
The rcsetup module contains the validation code for customization using
Matplotlib's rc settings.

Each rc setting is assigned a function used to validate any attempted changes
to that setting.  The validation functions are defined in the rcsetup module,
and are used to construct the rcParams global object which stores the settings
and is referenced throughout Matplotlib.

The default values of the rc settings are set in the default matplotlibrc file.
Any additions or deletions to the parameter set listed here should also be
propagated to the :file:`lib/matplotlib/mpl-data/matplotlibrc` in Matplotlib's
root source directory. New rcparams also need to be added to the RcKeyType enum
in :file:`lib/matplotlib/typing.py`.
"""


import ast
from dataclasses import dataclass
from functools import lru_cache, reduce
from numbers import Real
import operator
import os
import re
from typing import Any
from collections.abc import Callable

import numpy as np

import matplotlib as mpl
from matplotlib import _api, cbook
from matplotlib.backends import backend_registry
from matplotlib.cbook import ls_mapper
from matplotlib.colors import Colormap, is_color_like
from matplotlib._fontconfig_pattern import parse_fontconfig_pattern
from matplotlib._enums import JoinStyle, CapStyle

# Don't let the original cycler collide with our validating cycler
from cycler import Cycler, concat as cconcat, cycler as ccycler


class ValidateInStrings:
    def __init__(self, key, valid, ignorecase=False, *,
                 _deprecated_since=None):
        """*valid* is a list of legal strings."""
        self.key = key
        self.ignorecase = ignorecase
        self._deprecated_since = _deprecated_since

        def func(s):
            if ignorecase:
                return s.lower()
            else:
                return s
        self.valid = {func(k): k for k in valid}

    def __call__(self, s):
        if self._deprecated_since:
            name, = (k for k, v in globals().items() if v is self)
            _api.warn_deprecated(
                self._deprecated_since, name=name, obj_type="function")
        if self.ignorecase and isinstance(s, str):
            s = s.lower()
        if s in self.valid:
            return self.valid[s]
        msg = (f"{s!r} is not a valid value for {self.key}; supported values "
               f"are {[*self.valid.values()]}")
        if (isinstance(s, str)
                and (s.startswith('"') and s.endswith('"')
                     or s.startswith("'") and s.endswith("'"))
                and s[1:-1] in self.valid):
            msg += "; remove quotes surrounding your string"
        raise ValueError(msg)

    def __repr__(self):
        return (f"{self.__class__.__name__}("
                f"key={self.key!r}, valid={[*self.valid.values()]}, "
                f"ignorecase={self.ignorecase})")

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, ValidateInStrings):
            return NotImplemented
        return (
            self.key,
            self.ignorecase,
            self._deprecated_since,
            tuple(sorted(self.valid.items()))
        ) == (
            other.key,
            other.ignorecase,
            other._deprecated_since,
            tuple(sorted(other.valid.items()))
        )

    def __hash__(self):
        return hash((
            self.key,
            self.ignorecase,
            self._deprecated_since,
            tuple(sorted(self.valid.items()))
        ))


def _single_string_color_list(s, scalar_validator):
    """
    Convert the string *s* to a list of colors interpreting it either as a
    color sequence name, or a string containing single-letter colors.
    """
    try:
        colors = mpl.color_sequences[s]
    except KeyError:
        try:
            # Sometimes, a list of colors might be a single string
            # of single-letter colornames. So give that a shot.
            colors = [scalar_validator(v.strip()) for v in s if v.strip()]
        except ValueError:
            raise ValueError(f'{s!r} is neither a color sequence name nor can '
                             'it be interpreted as a list of colors')

    return colors


@lru_cache
def _listify_validator(scalar_validator, allow_stringlist=False, *,
                       n=None, doc=None):
    def f(s):
        if isinstance(s, str):
            try:
                val = [scalar_validator(v.strip()) for v in s.split(',')
                       if v.strip()]
            except Exception:
                if allow_stringlist:
                    # Special handling for colors
                    val = _single_string_color_list(s, scalar_validator)
                else:
                    raise
        # Allow any ordered sequence type -- generators, np.ndarray, pd.Series
        # -- but not sets, whose iteration order is non-deterministic.
        elif np.iterable(s) and not isinstance(s, (set, frozenset)):
            # The condition on this list comprehension will preserve the
            # behavior of filtering out any empty strings (behavior was
            # from the original validate_stringlist()), while allowing
            # any non-string/text scalar values such as numbers and arrays.
            val = [scalar_validator(v) for v in s
                   if not isinstance(v, str) or v]
        else:
            raise ValueError(
                f"Expected str or other non-set iterable, but got {s}")
        if n is not None and len(val) != n:
            raise ValueError(
                f"Expected {n} values, but there are {len(val)} values in {s}")
        return val

    try:
        f.__name__ = f"{scalar_validator.__name__}list"
    except AttributeError:  # class instance.
        f.__name__ = f"{type(scalar_validator).__name__}List"
    f.__qualname__ = f.__qualname__.rsplit(".", 1)[0] + "." + f.__name__
    f.__doc__ = doc if doc is not None else scalar_validator.__doc__
    return f


def validate_any(s):
    return s
validate_anylist = _listify_validator(validate_any)


def _validate_date(s):
    try:
        np.datetime64(s)
        return s
    except ValueError:
        raise ValueError(
            f'{s!r} should be a string that can be parsed by numpy.datetime64')


def validate_bool(b):
    """Convert b to ``bool`` or raise."""
    if isinstance(b, str):
        b = b.lower()
    if b in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
        return True
    elif b in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
        return False
    else:
        raise ValueError(f'Cannot convert {b!r} to bool')


def validate_axisbelow(s):
    try:
        return validate_bool(s)
    except ValueError:
        if isinstance(s, str):
            if s == 'line':
                return 'line'
    raise ValueError(f'{s!r} cannot be interpreted as'
                     ' True, False, or "line"')


def validate_dpi(s):
    """Confirm s is string 'figure' or convert s to float or raise."""
    if s == 'figure':
        return s
    try:
        return float(s)
    except ValueError as e:
        raise ValueError(f'{s!r} is not string "figure" and '
                         f'could not convert {s!r} to float') from e


def _make_type_validator(cls, *, allow_none=False):
    """
    Return a validator that converts inputs to *cls* or raises (and possibly
    allows ``None`` as well).
    """

    def validator(s):
        if (allow_none and
                (s is None or cbook._str_lower_equal(s, "none"))):
            if cbook._str_lower_equal(s, "none") and s != "None":
                _api.warn_deprecated(
                    "3.11",
                    message=f"Using the capitalization {s!r} in matplotlibrc for "
                            "*None* is deprecated in %(removal)s and will lead to an "
                            "error from version 3.13 onward. Please use 'None' "
                            "instead."
                )
            return None
        if cls is str and not isinstance(s, str):
            raise ValueError(f'Could not convert {s!r} to str')
        try:
            return cls(s)
        except (TypeError, ValueError) as e:
            raise ValueError(
                f'Could not convert {s!r} to {cls.__name__}') from e

    validator.__name__ = f"validate_{cls.__name__}"
    if allow_none:
        validator.__name__ += "_or_None"
    validator.__qualname__ = (
        validator.__qualname__.rsplit(".", 1)[0] + "." + validator.__name__)
    return validator


validate_string = _make_type_validator(str)
validate_string_or_None = _make_type_validator(str, allow_none=True)
validate_stringlist = _listify_validator(
    validate_string, doc='return a list of strings')
validate_int = _make_type_validator(int)
validate_int_or_None = _make_type_validator(int, allow_none=True)
validate_intlist = _listify_validator(validate_int, n=2)
validate_float = _make_type_validator(float)
validate_float_or_None = _make_type_validator(float, allow_none=True)
validate_floatlist = _listify_validator(
    validate_float)


def _validate_marker(s):
    try:
        return validate_int(s)
    except ValueError as e:
        try:
            return validate_string(s)
        except ValueError as e:
            raise ValueError('Supported markers are [string, int]') from e


_validate_markerlist = _listify_validator(
    _validate_marker, doc='return a list of markers')


def _validate_pathlike(s):
    if isinstance(s, (str, os.PathLike)):
        # Store value as str because savefig.directory needs to distinguish
        # between "" (cwd) and "." (cwd, but gets updated by user selections).
        return os.fsdecode(s)
    else:
        return validate_string(s)


def validate_fonttype(s):
    """
    Confirm that this is a Postscript or PDF font type that we know how to
    convert to.
    """
    fonttypes = {'type3':    3,
                 'truetype': 42}
    try:
        fonttype = validate_int(s)
    except ValueError:
        try:
            return fonttypes[s.lower()]
        except KeyError as e:
            raise ValueError('Supported Postscript/PDF font types are %s'
                             % list(fonttypes)) from e
    else:
        if fonttype not in fonttypes.values():
            raise ValueError(
                'Supported Postscript/PDF font types are %s' %
                list(fonttypes.values()))
        return fonttype


_auto_backend_sentinel = object()


def validate_backend(s):
    if s is _auto_backend_sentinel or backend_registry.is_valid_backend(s):
        return s
    else:
        msg = (f"'{s}' is not a valid value for backend; supported values are "
               f"{backend_registry.list_all()}")
        raise ValueError(msg)


def _validate_toolbar(s):
    s = ValidateInStrings(
        'toolbar', ['None', 'toolbar2', 'toolmanager'], ignorecase=True)(s)
    if s == 'toolmanager':
        _api.warn_external(
            "Treat the new Tool classes introduced in v1.5 as experimental "
            "for now; the API and rcParam may change in future versions.")
    return s


def validate_color_or_inherit(s):
    """Return a valid color arg."""
    if cbook._str_equal(s, 'inherit'):
        return s
    return validate_color(s)


def validate_color_or_auto(s):
    if cbook._str_equal(s, 'auto'):
        return s
    return validate_color(s)


def _validate_color_or_edge(s):
    if cbook._str_equal(s, 'edge'):
        return s
    return validate_color(s)


def validate_color_for_prop_cycle(s):
    # N-th color cycle syntax can't go into the color cycle.
    if isinstance(s, str) and re.match("^C[0-9]$", s):
        raise ValueError(f"Cannot put cycle reference ({s!r}) in prop_cycler")
    return validate_color(s)


def _validate_color_or_linecolor(s):
    if cbook._str_equal(s, 'linecolor'):
        return s
    elif cbook._str_equal(s, 'mfc') or cbook._str_equal(s, 'markerfacecolor'):
        return 'markerfacecolor'
    elif cbook._str_equal(s, 'mec') or cbook._str_equal(s, 'markeredgecolor'):
        return 'markeredgecolor'
    elif s is None:
        return None
    elif isinstance(s, str) and len(s) == 6 or len(s) == 8:
        stmp = '#' + s
        if is_color_like(stmp):
            return stmp
        if s.lower() == 'none':
            return None
    elif is_color_like(s):
        return s

    raise ValueError(f'{s!r} does not look like a color arg')


def validate_color(s):
    """Return a valid color arg."""
    if isinstance(s, str):
        if s.lower() == 'none':
            return 'none'
        if len(s) == 6 or len(s) == 8:
            stmp = '#' + s
            if is_color_like(stmp):
                return stmp

    if is_color_like(s):
        return s

    # If it is still valid, it must be a tuple (as a string from matplotlibrc).
    try:
        color = ast.literal_eval(s)
    except (SyntaxError, ValueError):
        pass
    else:
        if is_color_like(color):
            return color

    raise ValueError(f'{s!r} does not look like a color arg')


def _validate_color_or_None(s):
    if s is None or cbook._str_equal(s, "None"):
        return None
    return validate_color(s)


validate_colorlist = _listify_validator(
    validate_color, allow_stringlist=True, doc='return a list of colorspecs')


def _validate_cmap(s):
    _api.check_isinstance((str, Colormap), cmap=s)
    return s


def validate_aspect(s):
    if s in ('auto', 'equal'):
        return s
    try:
        return float(s)
    except ValueError as e:
        raise ValueError('not a valid aspect specification') from e


def validate_fontsize_None(s):
    if s is None or s == 'None':
        return None
    else:
        return validate_fontsize(s)


def validate_fontsize(s):
    fontsizes = ['xx-small', 'x-small', 'small', 'medium', 'large',
                 'x-large', 'xx-large', 'smaller', 'larger']
    if isinstance(s, str):
        s = s.lower()
    if s in fontsizes:
        return s
    try:
        return float(s)
    except ValueError as e:
        raise ValueError("%s is not a valid font size. Valid font sizes "
                         "are %s." % (s, ", ".join(fontsizes))) from e


validate_fontsizelist = _listify_validator(validate_fontsize)


def validate_fontweight(s):
    weights = [
        'ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman',
        'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black']
    # Note: Historically, weights have been case-sensitive in Matplotlib
    if s in weights:
        return s
    try:
        return int(s)
    except (ValueError, TypeError) as e:
        raise ValueError(f'{s} is not a valid font weight.') from e


def validate_fontstretch(s):
    stretchvalues = [
        'ultra-condensed', 'extra-condensed', 'condensed', 'semi-condensed',
        'normal', 'semi-expanded', 'expanded', 'extra-expanded',
        'ultra-expanded']
    # Note: Historically, stretchvalues have been case-sensitive in Matplotlib
    if s in stretchvalues:
        return s
    try:
        return int(s)
    except (ValueError, TypeError) as e:
        raise ValueError(f'{s} is not a valid font stretch.') from e


def validate_font_properties(s):
    parse_fontconfig_pattern(s)
    return s


def _validate_mathtext_fallback(s):
    _fallback_fonts = ['cm', 'stix', 'stixsans']
    if isinstance(s, str):
        s = s.lower()
    if s is None or s == 'none':
        return None
    elif s.lower() in _fallback_fonts:
        return s
    else:
        raise ValueError(
            f"{s} is not a valid fallback font name. Valid fallback font "
            f"names are {','.join(_fallback_fonts)}. Passing 'None' will turn "
            "fallback off.")


def validate_whiskers(s):
    try:
        return _listify_validator(validate_float, n=2)(s)
    except (TypeError, ValueError):
        try:
            return float(s)
        except ValueError as e:
            raise ValueError("Not a valid whisker value [float, "
                             "(float, float)]") from e


def validate_ps_distiller(s):
    if isinstance(s, str):
        s = s.lower()
    if s in ('none', None, 'false', False):
        return None
    else:
        return ValidateInStrings('ps.usedistiller', ['ghostscript', 'xpdf'])(s)


# A validator dedicated to the named line styles, based on the items in
# ls_mapper, and a list of possible strings read from Line2D.set_linestyle
_validate_named_linestyle = ValidateInStrings(
    'linestyle',
    [*ls_mapper.keys(), *ls_mapper.values(), 'None', 'none', ' ', ''],
    ignorecase=True)


def _validate_linestyle(ls):
    """
    A validator for all possible line styles, the named ones *and*
    the on-off ink sequences.
    """
    if isinstance(ls, str):
        try:  # Look first for a valid named line style, like '--' or 'solid'.
            return _validate_named_linestyle(ls)
        except ValueError:
            pass
        try:
            ls = ast.literal_eval(ls)  # Parsing matplotlibrc.
        except (SyntaxError, ValueError):
            pass  # Will error with the ValueError at the end.

    def _is_iterable_not_string_like(x):
        # Explicitly exclude bytes/bytearrays so that they are not
        # nonsensically interpreted as sequences of numbers (codepoints).
        return np.iterable(x) and not isinstance(x, (str, bytes, bytearray))

    if _is_iterable_not_string_like(ls):
        if len(ls) == 2 and _is_iterable_not_string_like(ls[1]):
            # (offset, (on, off, on, off, ...))
            offset, onoff = ls
        else:
            # For backcompat: (on, off, on, off, ...); the offset is implicit.
            offset = 0
            onoff = ls

        if (isinstance(offset, Real)
                and len(onoff) % 2 == 0
                and all(isinstance(elem, Real) for elem in onoff)):
            return (offset, onoff)

    raise ValueError(f"linestyle {ls!r} is not a valid on-off ink sequence.")


def _validate_linestyle_or_None(s):
    if s is None or cbook._str_equal(s, "None"):
        return None

    return _validate_linestyle(s)


validate_fillstyle = ValidateInStrings(
    'markers.fillstyle', ['full', 'left', 'right', 'bottom', 'top', 'none'])


validate_fillstylelist = _listify_validator(validate_fillstyle)


def validate_markevery(s):
    """
    Validate the markevery property of a Line2D object.

    Parameters
    ----------
    s : None, int, (int, int), slice, float, (float, float), or list[int]

    Returns
    -------
    None, int, (int, int), slice, float, (float, float), or list[int]
    """
    # Validate s against type slice float int and None
    if isinstance(s, (slice, float, int, type(None))):
        return s
    # Validate s against type tuple
    if isinstance(s, tuple):
        if (len(s) == 2
                and (all(isinstance(e, int) for e in s)
                     or all(isinstance(e, float) for e in s))):
            return s
        else:
            raise TypeError(
                "'markevery' tuple must be pair of ints or of floats")
    # Validate s against type list
    if isinstance(s, list):
        if all(isinstance(e, int) for e in s):
            return s
        else:
            raise TypeError(
                "'markevery' list must have all elements of type int")
    raise TypeError("'markevery' is of an invalid type")


validate_markeverylist = _listify_validator(validate_markevery)


def validate_bbox(s):
    if isinstance(s, str):
        s = s.lower()
        if s == 'tight':
            return s
        if s == 'standard':
            return None
        raise ValueError("bbox should be 'tight' or 'standard'")
    elif s is not None:
        # Backwards compatibility. None is equivalent to 'standard'.
        raise ValueError("bbox should be 'tight' or 'standard'")
    return s


def validate_sketch(s):

    if isinstance(s, str):
        s = s.lower().strip()
        if s.startswith("(") and s.endswith(")"):
            s = s[1:-1]
    if s == 'none' or s is None:
        return None
    try:
        return tuple(_listify_validator(validate_float, n=3)(s))
    except ValueError as exc:
        raise ValueError("Expected a (scale, length, randomness) tuple") from exc


def _validate_greaterthan_minushalf(s):
    s = validate_float(s)
    if s > -0.5:
        return s
    else:
        raise RuntimeError(f'Value must be >-0.5; got {s}')


def _validate_greaterequal0_lessequal1(s):
    s = validate_float(s)
    if 0 <= s <= 1:
        return s
    else:
        raise RuntimeError(f'Value must be >=0 and <=1; got {s}')


def _validate_int_greaterequal0(s):
    s = validate_int(s)
    if s >= 0:
        return s
    else:
        raise RuntimeError(f'Value must be >=0; got {s}')


def validate_hatch(s):
    r"""
    Validate a hatch pattern.
    A hatch pattern string can have any sequence of the following
    characters: ``\ / | - + * . x o O``.
    """
    if not isinstance(s, str):
        raise ValueError("Hatch pattern must be a string")
    _api.check_isinstance(str, hatch_pattern=s)
    unknown = set(s) - {'\\', '/', '|', '-', '+', '*', '.', 'x', 'o', 'O'}
    if unknown:
        raise ValueError("Unknown hatch symbol(s): %s" % list(unknown))
    return s


validate_hatchlist = _listify_validator(validate_hatch)
validate_dashlist = _listify_validator(validate_floatlist)


def _validate_minor_tick_ndivs(n):
    """
    Validate ndiv parameter related to the minor ticks.
    It controls the number of minor ticks to be placed between
    two major ticks.
    """

    if cbook._str_lower_equal(n, 'auto'):
        return n
    try:
        n = _validate_int_greaterequal0(n)
        return n
    except (RuntimeError, ValueError):
        pass

    raise ValueError("'tick.minor.ndivs' must be 'auto' or non-negative int")


_prop_validators = {
        'color': _listify_validator(validate_color_for_prop_cycle,
                                    allow_stringlist=True),
        'linewidth': validate_floatlist,
        'linestyle': _listify_validator(_validate_linestyle),
        'facecolor': validate_colorlist,
        'edgecolor': validate_colorlist,
        'joinstyle': _listify_validator(JoinStyle),
        'capstyle': _listify_validator(CapStyle),
        'fillstyle': validate_fillstylelist,
        'markerfacecolor': validate_colorlist,
        'markersize': validate_floatlist,
        'markeredgewidth': validate_floatlist,
        'markeredgecolor': validate_colorlist,
        'markevery': validate_markeverylist,
        'alpha': validate_floatlist,
        'marker': _validate_markerlist,
        'hatch': validate_hatchlist,
        'dashes': validate_dashlist,
    }
_prop_aliases = {
        'c': 'color',
        'lw': 'linewidth',
        'ls': 'linestyle',
        'fc': 'facecolor',
        'ec': 'edgecolor',
        'mfc': 'markerfacecolor',
        'mec': 'markeredgecolor',
        'mew': 'markeredgewidth',
        'ms': 'markersize',
    }


def cycler(*args, **kwargs):
    """
    Create a `~cycler.Cycler` object much like :func:`cycler.cycler`,
    but includes input validation.

    Call signatures::

      cycler(cycler)
      cycler(label=values, label2=values2, ...)
      cycler(label, values)

    Form 1 copies a given `~cycler.Cycler` object.

    Form 2 creates a `~cycler.Cycler` which cycles over one or more
    properties simultaneously. If multiple properties are given, their
    value lists must have the same length.

    Form 3 creates a `~cycler.Cycler` for a single property. This form
    exists for compatibility with the original cycler. Its use is
    discouraged in favor of the kwarg form, i.e. ``cycler(label=values)``.

    Parameters
    ----------
    cycler : Cycler
        Copy constructor for Cycler.

    label : str
        The property key. Must be a valid `.Artist` property.
        For example, 'color' or 'linestyle'. Aliases are allowed,
        such as 'c' for 'color' and 'lw' for 'linewidth'.

    values : iterable
        Finite-length iterable of the property values. These values
        are validated and will raise a ValueError if invalid.

    Returns
    -------
    Cycler
        A new :class:`~cycler.Cycler` for the given properties.

    Examples
    --------
    Creating a cycler for a single property:

    >>> c = cycler(color=['red', 'green', 'blue'])

    Creating a cycler for simultaneously cycling over multiple properties
    (e.g. red circle, green plus, blue cross):

    >>> c = cycler(color=['red', 'green', 'blue'],
    ...            marker=['o', '+', 'x'])

    """
    if args and kwargs:
        raise TypeError("cycler() can only accept positional OR keyword "
                        "arguments -- not both.")
    elif not args and not kwargs:
        raise TypeError("cycler() must have positional OR keyword arguments")

    if len(args) == 1:
        if not isinstance(args[0], Cycler):
            raise TypeError("If only one positional argument given, it must "
                            "be a Cycler instance.")
        return validate_cycler(args[0])
    elif len(args) == 2:
        pairs = [(args[0], args[1])]
    elif len(args) > 2:
        raise _api.nargs_error('cycler', '0-2', len(args))
    else:
        pairs = kwargs.items()

    validated = []
    for prop, vals in pairs:
        norm_prop = _prop_aliases.get(prop, prop)
        validator = _prop_validators.get(norm_prop, None)
        if validator is None:
            raise TypeError("Unknown artist property: %s" % prop)
        vals = validator(vals)
        # We will normalize the property names as well to reduce
        # the amount of alias handling code elsewhere.
        validated.append((norm_prop, vals))

    return reduce(operator.add, (ccycler(k, v) for k, v in validated))


def _parse_cycler_string(s):
    """
    Parse a string representation of a cycler into a Cycler object safely,
    without using eval().

    Accepts expressions like::

        cycler('color', ['r', 'g', 'b'])
        cycler('color', 'rgb') + cycler('linewidth', [1, 2, 3])
        cycler(c='rgb', lw=[1, 2, 3])
        cycler('c', 'rgb') * cycler('linestyle', ['-', '--'])
    """
    try:
        tree = ast.parse(s, mode='eval')
    except SyntaxError as e:
        raise ValueError(f"Could not parse {s!r}: {e}") from e
    return _eval_cycler_expr(tree.body)


def _eval_cycler_expr(node):
    """Recursively evaluate an AST node to build a Cycler object."""
    if isinstance(node, ast.BinOp):
        left = _eval_cycler_expr(node.left)
        right = _eval_cycler_expr(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Mult):
            return left * right
        raise ValueError(f"Unsupported operator: {type(node.op).__name__}")
    if isinstance(node, ast.Call):
        if not (isinstance(node.func, ast.Name)
                and node.func.id in ('cycler', 'concat')):
            raise ValueError(
                "only the 'cycler()' and 'concat()' functions are allowed")
        func = cycler if node.func.id == 'cycler' else cconcat
        args = [_eval_cycler_expr(a) for a in node.args]
        kwargs = {kw.arg: _eval_cycler_expr(kw.value) for kw in node.keywords}
        return func(*args, **kwargs)
    if isinstance(node, ast.Subscript):
        sl = node.slice
        if not isinstance(sl, ast.Slice):
            raise ValueError("only slicing is supported, not indexing")
        s = slice(
            ast.literal_eval(sl.lower) if sl.lower else None,
            ast.literal_eval(sl.upper) if sl.upper else None,
            ast.literal_eval(sl.step) if sl.step else None,
        )
        value = _eval_cycler_expr(node.value)
        return value[s]
    # Allow literal values (int, strings, lists, tuples) as arguments
    # to cycler() and concat().
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError):
        raise ValueError(
            f"Unsupported expression in cycler string: {ast.dump(node)}")


# A validator dedicated to the named legend loc
_validate_named_legend_loc = ValidateInStrings(
    'legend.loc',
    [
        "best",
        "upper right", "upper left", "lower left", "lower right", "right",
        "center left", "center right", "lower center", "upper center",
        "center"],
    ignorecase=True)


def _validate_legend_loc(loc):
    """
    Confirm that loc is a type which rc.Params["legend.loc"] supports.

    .. versionadded:: 3.8

    Parameters
    ----------
    loc : str | int | (float, float) | str((float, float))
        The location of the legend.

    Returns
    -------
    loc : str | int | (float, float) or raise ValueError exception
        The location of the legend.
    """
    if isinstance(loc, str):
        try:
            return _validate_named_legend_loc(loc)
        except ValueError:
            pass
        try:
            loc = ast.literal_eval(loc)
        except (SyntaxError, ValueError):
            pass
    if isinstance(loc, int):
        if 0 <= loc <= 10:
            return loc
    if isinstance(loc, tuple):
        if len(loc) == 2 and all(isinstance(e, Real) for e in loc):
            return loc
    raise ValueError(f"{loc} is not a valid legend location.")


def validate_cycler(s):
    """Return a Cycler object from a string repr or the object itself."""
    if isinstance(s, str):
        try:
            s = _parse_cycler_string(s)
        except Exception as e:
            raise ValueError(f"{s!r} is not a valid cycler construction: {e}"
                             ) from e
    if isinstance(s, Cycler):
        cycler_inst = s
    else:
        raise ValueError(f"Object is not a string or Cycler instance: {s!r}")

    unknowns = cycler_inst.keys - (set(_prop_validators) | set(_prop_aliases))
    if unknowns:
        raise ValueError("Unknown artist properties: %s" % unknowns)

    # Not a full validation, but it'll at least normalize property names
    # A fuller validation would require v0.10 of cycler.
    checker = set()
    for prop in cycler_inst.keys:
        norm_prop = _prop_aliases.get(prop, prop)
        if norm_prop != prop and norm_prop in cycler_inst.keys:
            raise ValueError(f"Cannot specify both {norm_prop!r} and alias "
                             f"{prop!r} in the same prop_cycle")
        if norm_prop in checker:
            raise ValueError(f"Another property was already aliased to "
                             f"{norm_prop!r}. Collision normalizing {prop!r}.")
        checker.update([norm_prop])

    # This is just an extra-careful check, just in case there is some
    # edge-case I haven't thought of.
    assert len(checker) == len(cycler_inst.keys)

    # Now, it should be safe to mutate this cycler
    for prop in cycler_inst.keys:
        norm_prop = _prop_aliases.get(prop, prop)
        cycler_inst.change_key(prop, norm_prop)

    for key, vals in cycler_inst.by_key().items():
        _prop_validators[key](vals)

    return cycler_inst


def validate_hist_bins(s):
    valid_strs = ["auto", "sturges", "fd", "doane", "scott", "rice", "sqrt"]
    if isinstance(s, str) and s in valid_strs:
        return s
    try:
        return int(s)
    except (TypeError, ValueError):
        pass
    try:
        return validate_floatlist(s)
    except ValueError:
        pass
    raise ValueError(f"'hist.bins' must be one of {valid_strs}, an int or"
                     " a sequence of floats")


class _ignorecase(list):
    """A marker class indicating that a list-of-str is case-insensitive."""


def _convert_validator_spec(key, conv):
    if isinstance(conv, list):
        ignorecase = isinstance(conv, _ignorecase)
        return ValidateInStrings(key, conv, ignorecase=ignorecase)
    else:
        return conv


# Mapping of rcParams to validators.
# Converters given as lists or _ignorecase are converted to ValidateInStrings
# immediately below.
# The rcParams defaults are defined in lib/matplotlib/mpl-data/matplotlibrc, which
# gets copied to matplotlib/mpl-data/matplotlibrc by the setup script.
_validators = {
    "backend":           validate_backend,
    "backend_fallback":  validate_bool,
    "figure.hooks":      validate_stringlist,
    "toolbar":           _validate_toolbar,
    "interactive":       validate_bool,
    "timezone":          validate_string,

    "webagg.port":            validate_int,
    "webagg.address":         validate_string,
    "webagg.open_in_browser": validate_bool,
    "webagg.port_retries":    validate_int,

    # line props
    "lines.linewidth":       validate_float,  # line width in points
    "lines.linestyle":       _validate_linestyle,  # solid line
    "lines.color":           validate_color,  # first color in color cycle
    "lines.marker":          _validate_marker,  # marker name
    "lines.markerfacecolor": validate_color_or_auto,  # default color
    "lines.markeredgecolor": validate_color_or_auto,  # default color
    "lines.markeredgewidth": validate_float,
    "lines.markersize":      validate_float,  # markersize, in points
    "lines.antialiased":     validate_bool,  # antialiased (no jaggies)
    "lines.dash_joinstyle":  JoinStyle,
    "lines.solid_joinstyle": JoinStyle,
    "lines.dash_capstyle":   CapStyle,
    "lines.solid_capstyle":  CapStyle,
    "lines.dashed_pattern":  validate_floatlist,
    "lines.dashdot_pattern": validate_floatlist,
    "lines.dotted_pattern":  validate_floatlist,
    "lines.scale_dashes":    validate_bool,

    # marker props
    "markers.fillstyle": validate_fillstyle,

    ## pcolor(mesh) props:
    "pcolor.shading": ["auto", "flat", "nearest", "gouraud"],
    "pcolormesh.snap": validate_bool,

    ## patch props
    "patch.linewidth":       validate_float,  # line width in points
    "patch.edgecolor":       validate_color,
    "patch.force_edgecolor": validate_bool,
    "patch.facecolor":       validate_color,  # first color in cycle
    "patch.antialiased":     validate_bool,  # antialiased (no jaggies)

    ## hatch props
    "hatch.color":     _validate_color_or_edge,
    "hatch.linewidth": validate_float,

    ## Histogram properties
    "hist.bins": validate_hist_bins,

    ## Boxplot properties
    "boxplot.notch":       validate_bool,
    "boxplot.vertical":    validate_bool,
    "boxplot.whiskers":    validate_whiskers,
    "boxplot.bootstrap":   validate_int_or_None,
    "boxplot.patchartist": validate_bool,
    "boxplot.showmeans":   validate_bool,
    "boxplot.showcaps":    validate_bool,
    "boxplot.showbox":     validate_bool,
    "boxplot.showfliers":  validate_bool,
    "boxplot.meanline":    validate_bool,

    "boxplot.flierprops.color":           validate_color,
    "boxplot.flierprops.marker":          _validate_marker,
    "boxplot.flierprops.markerfacecolor": validate_color_or_auto,
    "boxplot.flierprops.markeredgecolor": validate_color,
    "boxplot.flierprops.markeredgewidth": validate_float,
    "boxplot.flierprops.markersize":      validate_float,
    "boxplot.flierprops.linestyle":       _validate_linestyle,
    "boxplot.flierprops.linewidth":       validate_float,

    "boxplot.boxprops.color":     validate_color,
    "boxplot.boxprops.linewidth": validate_float,
    "boxplot.boxprops.linestyle": _validate_linestyle,

    "boxplot.whiskerprops.color":     validate_color,
    "boxplot.whiskerprops.linewidth": validate_float,
    "boxplot.whiskerprops.linestyle": _validate_linestyle,

    "boxplot.capprops.color":     validate_color,
    "boxplot.capprops.linewidth": validate_float,
    "boxplot.capprops.linestyle": _validate_linestyle,

    "boxplot.medianprops.color":     validate_color,
    "boxplot.medianprops.linewidth": validate_float,
    "boxplot.medianprops.linestyle": _validate_linestyle,

    "boxplot.meanprops.color":           validate_color,
    "boxplot.meanprops.marker":          _validate_marker,
    "boxplot.meanprops.markerfacecolor": validate_color,
    "boxplot.meanprops.markeredgecolor": validate_color,
    "boxplot.meanprops.markersize":      validate_float,
    "boxplot.meanprops.linestyle":       _validate_linestyle,
    "boxplot.meanprops.linewidth":       validate_float,

    ## font props
    "font.enable_last_resort":     validate_bool,
    "font.family":     validate_stringlist,  # used by text object
    "font.style":      validate_string,
    "font.variant":    validate_string,
    "font.stretch":    validate_fontstretch,
    "font.weight":     validate_fontweight,
    "font.size":       validate_float,  # Base font size in points
    "font.serif":      validate_stringlist,
    "font.sans-serif": validate_stringlist,
    "font.cursive":    validate_stringlist,
    "font.fantasy":    validate_stringlist,
    "font.monospace":  validate_stringlist,

    # text props
    "text.color":          validate_color,
    "text.usetex":         validate_bool,
    "text.latex.engine":   ["latex", "latex+dvipng"],
    "text.latex.preamble": validate_string,
    "text.hinting":        ["default", "no_autohint", "force_autohint",
                            "no_hinting", "auto", "native", "either", "none"],
    "text.hinting_factor": validate_int_or_None,
    "text.kerning_factor": validate_int_or_None,
    "text.antialiased":    validate_bool,
    "text.parse_math":     validate_bool,
    "text.language":       validate_string_or_None,

    "mathtext.cal":            validate_font_properties,
    "mathtext.rm":             validate_font_properties,
    "mathtext.tt":             validate_font_properties,
    "mathtext.it":             validate_font_properties,
    "mathtext.bf":             validate_font_properties,
    "mathtext.bfit":           validate_font_properties,
    "mathtext.sf":             validate_font_properties,
    "mathtext.fontset":        ["dejavusans", "dejavuserif", "cm", "stix",
                                "stixsans", "custom"],
    "mathtext.default":        ["rm", "cal", "bfit", "it", "tt", "sf", "bf", "default",
                                "bb", "frak", "scr", "regular", "normal"],
    "mathtext.fallback":       _validate_mathtext_fallback,

    "image.aspect":              validate_aspect,  # equal, auto, a number
    "image.interpolation":       validate_string,
    "image.interpolation_stage": ["auto", "data", "rgba"],
    "image.cmap":                _validate_cmap,  # gray, jet, etc.
    "image.lut":                 validate_int,  # lookup table
    "image.origin":              ["upper", "lower"],
    "image.resample":            validate_bool,
    # Specify whether vector graphics backends will combine all images on a
    # set of Axes into a single composite image
    "image.composite_image": validate_bool,

    # contour props
    "contour.negative_linestyle": _validate_linestyle,
    "contour.corner_mask":        validate_bool,
    "contour.linewidth":          validate_float_or_None,
    "contour.algorithm":          ["mpl2005", "mpl2014", "serial", "threaded"],

    # errorbar props
    "errorbar.capsize": validate_float,
    "errorbar.capthick": validate_float_or_None,
    "errorbar.elinewidth": validate_float_or_None,

    # axis props
    # alignment of x/y axis title
    "xaxis.labellocation": ["left", "center", "right"],
    "yaxis.labellocation": ["bottom", "center", "top"],

    # Axes props
    "axes.axisbelow":        validate_axisbelow,
    "axes.facecolor":        validate_color,  # background color
    "axes.edgecolor":        validate_color,  # edge color
    "axes.linewidth":        validate_float,  # edge linewidth

    "axes.spines.left":      validate_bool,  # Set visibility of axes spines,
    "axes.spines.right":     validate_bool,  # i.e., the lines around the chart
    "axes.spines.bottom":    validate_bool,  # denoting data boundary.
    "axes.spines.top":       validate_bool,

    "axes.titlesize":     validate_fontsize,  # Axes title fontsize
    "axes.titlelocation": ["left", "center", "right"],  # Axes title alignment
    "axes.titleweight":   validate_fontweight,  # Axes title font weight
    "axes.titlecolor":    validate_color_or_auto,  # Axes title font color
    # title location, axes units, None means auto
    "axes.titley":        validate_float_or_None,
    # pad from Axes top decoration to title in points
    "axes.titlepad":      validate_float,
    "axes.grid":          validate_bool,  # display grid or not
    "axes.grid.which":    ["minor", "both", "major"],  # which grids are drawn
    "axes.grid.axis":     ["x", "y", "both"],  # grid type
    "axes.labelsize":     validate_fontsize,  # fontsize of x & y labels
    "axes.labelpad":      validate_float,  # space between label and axis
    "axes.labelweight":   validate_fontweight,  # fontsize of x & y labels
    "axes.labelcolor":    validate_color,  # color of axis label
    # use scientific notation if log10 of the axis range is smaller than the
    # first or larger than the second
    "axes.formatter.limits": validate_intlist,
    # use current locale to format ticks
    "axes.formatter.use_locale": validate_bool,
    "axes.formatter.use_mathtext": validate_bool,
    # minimum exponent to format in scientific notation
    "axes.formatter.min_exponent": validate_int,
    "axes.formatter.useoffset": validate_bool,
    "axes.formatter.offset_threshold": validate_int,
    "axes.unicode_minus": validate_bool,
    # This entry can be either a cycler object or a string repr of a
    # cycler-object, which is parsed safely via AST.
    "axes.prop_cycle": validate_cycler,
    # If "data", axes limits are set close to the data.
    # If "round_numbers" axes limits are set to the nearest round numbers.
    "axes.autolimit_mode": ["data", "round_numbers"],
    "axes.xmargin": _validate_greaterthan_minushalf,  # margin added to xaxis
    "axes.ymargin": _validate_greaterthan_minushalf,  # margin added to yaxis
    "axes.zmargin": _validate_greaterthan_minushalf,  # margin added to zaxis

    "polaraxes.grid":    validate_bool,  # display polar grid or not
    "axes3d.grid":       validate_bool,  # display 3d grid
    "axes3d.automargin": validate_bool,  # automatically add margin when
                                         # manually setting 3D axis limits

    "axes3d.xaxis.panecolor":    validate_color,  # 3d background pane
    "axes3d.yaxis.panecolor":    validate_color,  # 3d background pane
    "axes3d.zaxis.panecolor":    validate_color,  # 3d background pane

    "axes3d.depthshade": validate_bool,  # depth shade for 3D scatter plots
    "axes3d.depthshade_minalpha": validate_float,  # min alpha value for depth shading

    "axes3d.mouserotationstyle": ["azel", "trackball", "sphere", "arcball"],
    "axes3d.trackballsize": validate_float,
    "axes3d.trackballborder": validate_float,
    "axes3d.snap_rotation": validate_float,

    # scatter props
    "scatter.marker":     _validate_marker,
    "scatter.edgecolors": validate_string,

    "date.epoch": _validate_date,
    "date.autoformatter.year":        validate_string,
    "date.autoformatter.month":       validate_string,
    "date.autoformatter.day":         validate_string,
    "date.autoformatter.hour":        validate_string,
    "date.autoformatter.minute":      validate_string,
    "date.autoformatter.second":      validate_string,
    "date.autoformatter.microsecond": validate_string,

    'date.converter':          ['auto', 'concise'],
    # for auto date locator, choose interval_multiples
    'date.interval_multiples': validate_bool,

    # legend properties
    "legend.fancybox": validate_bool,
    "legend.loc": _validate_legend_loc,

    # the number of points in the legend line
    "legend.numpoints":      validate_int,
    # the number of points in the legend line for scatter
    "legend.scatterpoints":  validate_int,
    "legend.fontsize":       validate_fontsize,
    "legend.title_fontsize": validate_fontsize_None,
    # color of the legend
    "legend.labelcolor":     _validate_color_or_linecolor,
    # the relative size of legend markers vs. original
    "legend.markerscale":    validate_float,
    # using dict in rcParams not yet supported, so make sure it is bool
    "legend.shadow":         validate_bool,
    # whether or not to draw a frame around legend
    "legend.frameon":        validate_bool,
    # alpha value of the legend frame
    "legend.framealpha":     validate_float_or_None,
    # linewidth of legend frame
    "legend.linewidth": validate_float_or_None,

    ## the following dimensions are in fraction of the font size
    "legend.borderpad":      validate_float,  # units are fontsize
    # the vertical space between the legend entries
    "legend.labelspacing":   validate_float,
    # the length of the legend lines
    "legend.handlelength":   validate_float,
    # the length of the legend lines
    "legend.handleheight":   validate_float,
    # the space between the legend line and legend text
    "legend.handletextpad":  validate_float,
    # the border between the Axes and legend edge
    "legend.borderaxespad":  validate_float,
    # the border between the Axes and legend edge
    "legend.columnspacing":  validate_float,
    "legend.facecolor":      validate_color_or_inherit,
    "legend.edgecolor":      validate_color_or_inherit,

    # tick properties
    "xtick.top":           validate_bool,      # draw ticks on top side
    "xtick.bottom":        validate_bool,      # draw ticks on bottom side
    "xtick.labeltop":      validate_bool,      # draw label on top
    "xtick.labelbottom":   validate_bool,      # draw label on bottom
    "xtick.major.size":    validate_float,     # major xtick size in points
    "xtick.minor.size":    validate_float,     # minor xtick size in points
    "xtick.major.width":   validate_float,     # major xtick width in points
    "xtick.minor.width":   validate_float,     # minor xtick width in points
    "xtick.major.pad":     validate_float,     # distance to label in points
    "xtick.minor.pad":     validate_float,     # distance to label in points
    "xtick.color":         validate_color,     # color of xticks
    "xtick.labelcolor":    validate_color_or_inherit,  # color of xtick labels
    "xtick.minor.visible": validate_bool,      # visibility of minor xticks
    "xtick.minor.top":     validate_bool,      # draw top minor xticks
    "xtick.minor.bottom":  validate_bool,      # draw bottom minor xticks
    "xtick.major.top":     validate_bool,      # draw top major xticks
    "xtick.major.bottom":  validate_bool,      # draw bottom major xticks
    # number of minor xticks
    "xtick.minor.ndivs":   _validate_minor_tick_ndivs,
    "xtick.labelsize":     validate_fontsize,  # fontsize of xtick labels
    "xtick.direction":     ["out", "in", "inout"],  # direction of xticks
    "xtick.alignment":     ["center", "right", "left"],

    "ytick.left":          validate_bool,      # draw ticks on left side
    "ytick.right":         validate_bool,      # draw ticks on right side
    "ytick.labelleft":     validate_bool,      # draw tick labels on left side
    "ytick.labelright":    validate_bool,      # draw tick labels on right side
    "ytick.major.size":    validate_float,     # major ytick size in points
    "ytick.minor.size":    validate_float,     # minor ytick size in points
    "ytick.major.width":   validate_float,     # major ytick width in points
    "ytick.minor.width":   validate_float,     # minor ytick width in points
    "ytick.major.pad":     validate_float,     # distance to label in points
    "ytick.minor.pad":     validate_float,     # distance to label in points
    "ytick.color":         validate_color,     # color of yticks
    "ytick.labelcolor":    validate_color_or_inherit,  # color of ytick labels
    "ytick.minor.visible": validate_bool,      # visibility of minor yticks
    "ytick.minor.left":    validate_bool,      # draw left minor yticks
    "ytick.minor.right":   validate_bool,      # draw right minor yticks
    "ytick.major.left":    validate_bool,      # draw left major yticks
    "ytick.major.right":   validate_bool,      # draw right major yticks
    # number of minor yticks
    "ytick.minor.ndivs":   _validate_minor_tick_ndivs,
    "ytick.labelsize":     validate_fontsize,  # fontsize of ytick labels
    "ytick.direction":     ["out", "in", "inout"],  # direction of yticks
    "ytick.alignment":     [
        "center", "top", "bottom", "baseline", "center_baseline"],

    "grid.color":        validate_color,  # grid color
    "grid.linestyle":    _validate_linestyle,  # solid
    "grid.linewidth":    validate_float,     # in points
    "grid.alpha":        validate_float,

    "grid.major.color":        _validate_color_or_None,  # grid color
    "grid.major.linestyle":    _validate_linestyle_or_None,  # solid
    "grid.major.linewidth":    validate_float_or_None,     # in points
    "grid.major.alpha":        validate_float_or_None,

    "grid.minor.color":        _validate_color_or_None,  # grid color
    "grid.minor.linestyle":    _validate_linestyle_or_None,  # solid
    "grid.minor.linewidth":    validate_float_or_None,     # in points
    "grid.minor.alpha":        validate_float_or_None,

    ## figure props
    # figure title
    "figure.titlesize":   validate_fontsize,
    "figure.titleweight": validate_fontweight,

    # figure labels
    "figure.labelsize":   validate_fontsize,
    "figure.labelweight": validate_fontweight,

    # figure size in inches: width by height
    "figure.figsize":          _listify_validator(validate_float, n=2),
    "figure.dpi":              validate_float,
    "figure.facecolor":        validate_color,
    "figure.edgecolor":        validate_color,
    "figure.frameon":          validate_bool,
    "figure.autolayout":       validate_bool,
    "figure.max_open_warning": validate_int,
    "figure.raise_window":     validate_bool,
    "macosx.window_mode":      ["system", "tab", "window"],

    "figure.subplot.left":   validate_float,
    "figure.subplot.right":  validate_float,
    "figure.subplot.bottom": validate_float,
    "figure.subplot.top":    validate_float,
    "figure.subplot.wspace": validate_float,
    "figure.subplot.hspace": validate_float,

    "figure.constrained_layout.use": validate_bool,  # run constrained_layout?
    # wspace and hspace are fraction of adjacent subplots to use for space.
    # Much smaller than above because we don't need room for the text.
    "figure.constrained_layout.hspace": validate_float,
    "figure.constrained_layout.wspace": validate_float,
    # buffer around the Axes, in inches.
    "figure.constrained_layout.h_pad": validate_float,
    "figure.constrained_layout.w_pad": validate_float,

    ## Saving figure's properties
    'savefig.dpi':          validate_dpi,
    'savefig.facecolor':    validate_color_or_auto,
    'savefig.edgecolor':    validate_color_or_auto,
    'savefig.orientation':  ['landscape', 'portrait'],
    "savefig.format":       validate_string,
    "savefig.bbox":         validate_bbox,  # "tight", or "standard" (= None)
    "savefig.pad_inches":   validate_float,
    # default directory in savefig dialog box
    "savefig.directory":    _validate_pathlike,
    "savefig.transparent":  validate_bool,

    "tk.window_focus": validate_bool,  # Maintain shell focus for TkAgg

    # Set the papersize/type
    "ps.papersize":       _ignorecase(
                                ["figure", "letter", "legal", "ledger",
                                 *[f"{ab}{i}" for ab in "ab" for i in range(11)]]),
    "ps.useafm":          validate_bool,
    # use ghostscript or xpdf to distill ps output
    "ps.usedistiller":    validate_ps_distiller,
    "ps.distiller.res":   validate_int,  # dpi
    "ps.fonttype":        validate_fonttype,  # 3 (Type3) or 42 (Truetype)
    "pdf.compression":    validate_int,  # 0-9 compression level; 0 to disable
    "pdf.inheritcolor":   validate_bool,  # skip color setting commands
    # use only the 14 PDF core fonts embedded in every PDF viewing application
    "pdf.use14corefonts": validate_bool,
    "pdf.fonttype":       validate_fonttype,  # 3 (Type3) or 42 (Truetype)

    "pgf.texsystem": ["xelatex", "lualatex", "pdflatex"],  # latex variant used
    "pgf.rcfonts":   validate_bool,  # use mpl's rc settings for font config
    "pgf.preamble":  validate_string,  # custom LaTeX preamble

    # write raster image data into the svg file
    "svg.image_inline": validate_bool,
    "svg.fonttype": ["none", "path"],  # save text as text ("none") or "paths"
    "svg.hashsalt": validate_string_or_None,
    "svg.id": validate_string_or_None,

    # set this when you want to generate hardcopy docstring
    "docstring.hardcopy": validate_bool,

    "path.simplify":           validate_bool,
    "path.simplify_threshold": _validate_greaterequal0_lessequal1,
    "path.snap":               validate_bool,
    "path.sketch":             validate_sketch,
    "path.effects":            validate_anylist,
    "agg.path.chunksize":      validate_int,  # 0 to disable chunking

    # key-mappings (multi-character mappings should be a list/tuple)
    "keymap.fullscreen": validate_stringlist,
    "keymap.home":       validate_stringlist,
    "keymap.back":       validate_stringlist,
    "keymap.forward":    validate_stringlist,
    "keymap.pan":        validate_stringlist,
    "keymap.zoom":       validate_stringlist,
    "keymap.save":       validate_stringlist,
    "keymap.quit":       validate_stringlist,
    "keymap.quit_all":   validate_stringlist,  # e.g.: "W", "cmd+W", "Q"
    "keymap.grid":       validate_stringlist,
    "keymap.grid_minor": validate_stringlist,
    "keymap.yscale":     validate_stringlist,
    "keymap.xscale":     validate_stringlist,
    "keymap.help":       validate_stringlist,
    "keymap.copy":       validate_stringlist,

    # Animation settings
    "animation.html":         ["html5", "jshtml", "none"],
    # Limit, in MB, of size of base64 encoded animation in HTML
    # (i.e. IPython notebook)
    "animation.embed_limit":  validate_float,
    "animation.writer":       validate_string,
    "animation.codec":        validate_string,
    "animation.bitrate":      validate_int,
    # Controls image format when frames are written to disk
    "animation.frame_format": ["png", "jpeg", "tiff", "raw", "rgba", "ppm",
                               "sgi", "bmp", "pbm", "svg"],
    # Path to ffmpeg binary. If just binary name, subprocess uses $PATH.
    "animation.ffmpeg_path":  _validate_pathlike,
    # Additional arguments for ffmpeg movie writer (using pipes)
    "animation.ffmpeg_args":  validate_stringlist,
     # Path to convert binary. If just binary name, subprocess uses $PATH.
    "animation.convert_path": _validate_pathlike,
     # Additional arguments for convert movie writer (using pipes)
    "animation.convert_args": validate_stringlist,

    # Classic (pre 2.0) compatibility mode
    # This is used for things that are hard to make backward compatible
    # with a sane rcParam alone.  This does *not* turn on classic mode
    # altogether.  For that use `matplotlib.style.use("classic")`.
    "_internal.classic_mode": validate_bool
}
_hardcoded_defaults = {  # Defaults not inferred from
    # lib/matplotlib/mpl-data/matplotlibrc...
    # ... because they are private:
    "_internal.classic_mode": False,
    # ... because they are deprecated:
    # No current deprecations.
    # backend is handled separately when constructing rcParamsDefault.
}
_validators = {k: _convert_validator_spec(k, conv)
               for k, conv in _validators.items()}


@dataclass
class _Param:
    name: str
    default: Any
    validator: Callable[[Any], Any]
    description: str = None


@dataclass
class _Section:
    title: str
    description: str = None


@dataclass
class _Subsection:
    title: str
    description: str = None


# Definition of all rcParams. This is currently only used to generate the documentation.
#
# Actual runtime values do still come from the historic sources:
# - available parameters and defaults: lib/matplotlib/mpl-data/matplotlibrc
# - validators: _validators, see above
#
# The structure and format of this definition is not fixed and may change in the future.
# It's a work-in-progress state towards a consistent and more structured definition of
# rcParams that can be used both for documentation and runtime. The goal is to
# eventually eliminate the old sources of defaults and validators and have this be the
# single source of truth.
#
# In the transition phase, consistency is ensured via tests.
_DEFINITION = [
    _Section("Backends"),
    _Param(
        "webagg.port",
        default=8988,
        validator=validate_int,
        description="The port to use for the web server in the WebAgg backend."
    ),
    _Param(
        "webagg.address",
        default="127.0.0.1",
        validator=validate_string,
        description="The address on which the WebAgg web server should be reachable."
    ),
    _Param(
        "webagg.port_retries",
        default=50,
        validator=validate_int,
        description="If webagg.port is unavailable, a number of other random ports "
                    "will be tried until one that is available is found."
    ),
    _Param(
        "webagg.open_in_browser",
        default=True,
        validator=validate_bool,
        description="When True, open the web browser to the plot that is shown"
    ),
    _Param(
        "backend_fallback",
        default=True,
        validator=validate_bool,
        description="If you are running pyplot inside a GUI and your backend choice "
                    "conflicts, we will automatically try to find a compatible one for "
                    "you if backend_fallback is True"
    ),
    _Param(
        "interactive",
        default=False,
        validator=validate_bool
    ),
    _Param(
        "figure.hooks",
        default=[],
        validator=validate_stringlist,
        description="list of dotted.module.name:dotted.callable.name"
    ),
    _Param(
        "toolbar",
        default="toolbar2",
        validator=_validate_toolbar,
        description="{None, toolbar2, toolmanager}"
    ),
    _Param(
        "timezone",
        default="UTC",
        validator=validate_string,
        description="a pytz timezone string, e.g., US/Central or Europe/Paris"
    ),
    _Section(
        "Lines",
        description="Default properties for line objects, such as those returned by "
                    "plot()."
    ),
    _Param(
        "lines.linewidth",
        default=1.5,
        validator=validate_float,
        description="line width in points"
    ),
    _Param(
        "lines.linestyle",
        default="-",
        validator=_validate_linestyle,
        description="solid line"
    ),
    _Param(
        "lines.color",
        default="C0",
        validator=validate_color,
        description="has no affect on plot(); see axes.prop_cycle"
    ),
    _Param(
        "lines.marker",
        default="None",
        validator=_validate_marker,
        description="the default marker"
    ),
    _Param(
        "lines.markerfacecolor",
        default="auto",
        validator=validate_color_or_auto,
        description="the default marker face color"
    ),
    _Param(
        "lines.markeredgecolor",
        default="auto",
        validator=validate_color_or_auto,
        description="the default marker edge color"
    ),
    _Param(
        "lines.markeredgewidth",
        default=1.0,
        validator=validate_float,
        description="the line width around the marker symbol"
    ),
    _Param(
        "lines.markersize",
        default=6.0,
        validator=validate_float,
        description="marker size, in points"
    ),
    _Param(
        "lines.dash_joinstyle",
        default="round",
        validator=JoinStyle,
        description="{miter, round, bevel}"
    ),
    _Param(
        "lines.dash_capstyle",
        default="butt",
        validator=CapStyle,
        description="{butt, round, projecting}"
    ),
    _Param(
        "lines.solid_joinstyle",
        default="round",
        validator=JoinStyle,
        description="{miter, round, bevel}"
    ),
    _Param(
        "lines.solid_capstyle",
        default="projecting",
        validator=CapStyle,
        description="{butt, round, projecting}"
    ),
    _Param(
        "lines.antialiased",
        default=True,
        validator=validate_bool,
        description="render lines in antialiased (no jaggies)"
    ),
    _Param(
        "lines.dashed_pattern",
        default=[3.7, 1.6],
        validator=validate_floatlist,
        description="The dash pattern for linestyle 'dashed'"
    ),
    _Param(
        "lines.dashdot_pattern",
        default=[6.4, 1.6, 1.0, 1.6],
        validator=validate_floatlist,
        description="The dash pattern for linestyle 'dashdot'"
    ),
    _Param(
        "lines.dotted_pattern",
        default=[1.0, 1.65],
        validator=validate_floatlist,
        description="The dash pattern for linestyle 'dotted'"
    ),
    _Param(
        "lines.scale_dashes",
        default=True,
        validator=validate_bool
    ),
    _Param(
        "markers.fillstyle",
        default="full",
        validator=validate_fillstyle,
        description="{full, left, right, bottom, top, none}"
    ),
    _Param(
        "pcolor.shading",
        default="auto",
        validator=["auto", "flat", "nearest", "gouraud"]
    ),
    _Param(
        "pcolormesh.snap",
        default=True,
        validator=validate_bool,
        description="Whether to snap the mesh to pixel boundaries. This is provided "
                    "solely to allow old test images to remain unchanged. Set to False "
                    "to obtain the previous behavior."
    ),
    _Section("Patches"),
    _Param(
        "patch.linewidth",
        default=1.0,
        validator=validate_float,
        description="edge width in points."
    ),
    _Param(
        "patch.facecolor",
        default="C0",
        validator=validate_color
    ),
    _Param(
        "patch.edgecolor",
        default="black",
        validator=validate_color,
        description='By default, Patches and Collections do not draw edges. This value '
                    'is only used if facecolor is "none" (an Artist without facecolor '
                    'and edgecolor would be invisible)  or if patch.force_edgecolor '
                    'is True.'
    ),
    _Param(
        "patch.force_edgecolor",
        default=False,
        validator=validate_bool,
        description="By default, Patches and Collections do not draw edges. Set this "
                    "to True to draw edges with patch.edgedcolor as the default "
                    "edgecolor. This is mainly relevant for styles."
    ),
    _Param(
        "patch.antialiased",
        default=True,
        validator=validate_bool,
        description="render patches in antialiased (no jaggies)"
    ),
    _Section("Hatches"),
    _Param("hatch.color", "edge", _validate_color_or_edge),
    _Param("hatch.linewidth", 1.0, validate_float),
    _Section("Boxplot"),
    _Param("boxplot.notch", False, validate_bool),
    _Param("boxplot.vertical", True, validate_bool),
    _Param("boxplot.whiskers", 1.5, validate_whiskers),
    _Param("boxplot.bootstrap", None, validate_int_or_None),
    _Param("boxplot.patchartist", False, validate_bool),
    _Param("boxplot.showmeans", False, validate_bool),
    _Param("boxplot.showcaps", True, validate_bool),
    _Param("boxplot.showbox", True, validate_bool),
    _Param("boxplot.showfliers", True, validate_bool),
    _Param("boxplot.meanline", False, validate_bool),
    _Param("boxplot.flierprops.color", "black", validate_color),
    _Param("boxplot.flierprops.marker", "o", _validate_marker),
    _Param("boxplot.flierprops.markerfacecolor", "none", validate_color_or_auto),
    _Param("boxplot.flierprops.markeredgecolor", "black", validate_color),
    _Param("boxplot.flierprops.markeredgewidth", 1.0, validate_float),
    _Param("boxplot.flierprops.markersize", 6.0, validate_float),
    _Param("boxplot.flierprops.linestyle", "none", _validate_linestyle),
    _Param("boxplot.flierprops.linewidth", 1.0, validate_float),
    _Param("boxplot.boxprops.color", "black", validate_color),
    _Param("boxplot.boxprops.linewidth", 1.0, validate_float),
    _Param("boxplot.boxprops.linestyle", "-", _validate_linestyle),
    _Param("boxplot.whiskerprops.color", "black", validate_color),
    _Param("boxplot.whiskerprops.linewidth", 1.0, validate_float),
    _Param("boxplot.whiskerprops.linestyle", "-", _validate_linestyle),
    _Param("boxplot.capprops.color", "black", validate_color),
    _Param("boxplot.capprops.linewidth", 1.0, validate_float),
    _Param("boxplot.capprops.linestyle", "-", _validate_linestyle),
    _Param("boxplot.medianprops.color", "C1", validate_color),
    _Param("boxplot.medianprops.linewidth", 1.0, validate_float),
    _Param("boxplot.medianprops.linestyle", "-", _validate_linestyle),
    _Param("boxplot.meanprops.color", "C2", validate_color),
    _Param("boxplot.meanprops.marker", "^", _validate_marker),
    _Param("boxplot.meanprops.markerfacecolor", "C2", validate_color),
    _Param("boxplot.meanprops.markeredgecolor", "C2", validate_color),
    _Param("boxplot.meanprops.markersize", 6.0, validate_float),
    _Param("boxplot.meanprops.linestyle", "--", _validate_linestyle),
    _Param("boxplot.meanprops.linewidth", 1.0, validate_float),
    _Section(
        "Font",
        description="The font properties used by `.Text` "
                    "See https://matplotlib.org/stable/api/font_manager_api.html for "
                    "more information on font properties. The 6 font properties used "
                    "for font matching are given below with their default values."
    ),
    _Param("font.family", ["sans-serif"], validate_stringlist),
    _Param("font.style", "normal", validate_string),
    _Param("font.variant", "normal", validate_string),
    _Param("font.weight", "normal", validate_fontweight),
    _Param("font.stretch", "normal", validate_fontstretch),
    _Param("font.size", 10.0, validate_float),
    _Param(
        "font.serif",
        default=[
            "DejaVu Serif", "Bitstream Vera Serif", "Computer Modern Roman",
            "New Century Schoolbook", "Century Schoolbook L", "Utopia", "ITC Bookman",
            "Bookman", "Nimbus Roman No9 L", "Times New Roman", "Times", "Palatino",
            "Charter", "serif",
        ],
        validator=validate_stringlist
    ),
    _Param(
        "font.sans-serif",
        default=[
            "DejaVu Sans", "Bitstream Vera Sans", "Computer Modern Sans Serif",
            "Lucida Grande", "Verdana", "Geneva", "Lucid", "Arial", "Helvetica",
            "Avant Garde", "sans-serif",
        ],
        validator=validate_stringlist
    ),
    _Param(
        "font.cursive",
        default=[
            "Apple Chancery", "Textile", "Zapf Chancery", "Sand", "Script MT", "Felipa",
            "Comic Neue", "Comic Sans MS", "cursive",
        ],
        validator=validate_stringlist
    ),
    _Param(
        "font.fantasy",
        default=["Chicago", "Charcoal", "Impact", "Western", "xkcd script", "fantasy"],
        validator=validate_stringlist
    ),
    _Param(
        "font.monospace",
        default=[
            "DejaVu Sans Mono", "Bitstream Vera Sans Mono",
            "Computer Modern Typewriter", "Andale Mono", "Nimbus Mono L", "Courier New",
            "Courier", "Fixed", "Terminal", "monospace",
        ],
        validator=validate_stringlist
    ),
    _Param(
        "font.enable_last_resort",
        default=True,
        validator=validate_bool,
        description="If True, then Unicode Consortium's Last Resort font will be "
                    "appended to all font selections. This ensures that there will "
                    "always be a glyph displayed."
    ),
    _Section("Text properties"),
    _Param(
        "text.color",
        default="black",
        validator=validate_color
    ),
    _Param(
        "text.language",
        default=None,
        validator=validate_string_or_None,
        description="The language of the text in a format accepted by libraqm, namely "
                    "`a BCP47 language code "
                    "<https://www.w3.org/International/articles/language-tags/>`_. If "
                    "None, then no particular language will be implied, and default "
                    "font settings will be used."
    ),
    _Param(
        "text.hinting",
        default="default",
        validator=[
            "default", "no_autohint", "force_autohint", "no_hinting", "auto", "native",
            "either", "none",
        ],
        description="FreeType hinting flag (\"foo\" corresponds to FT_LOAD_FOO); may "
                    "be one of the following (Proprietary Matplotlib-specific synonyms "
                    "are given in parentheses, but their use is discouraged): "
                    "- default: Use the font's native hinter if possible, else "
                    "  FreeType's auto-hinter. (\"either\" is a synonym)."
                    "- no_autohint: Use the font's native hinter if possible, else "
                    "  don't hint. (\"native\" is a synonym.)"
                    "- force_autohint: Use FreeType's auto-hinter. (\"auto\" is a "
                    "  synonym.)"
                    "- no_hinting: Disable hinting. (\"none\" is a synonym.)"
    ),
    _Param(
        "text.hinting_factor",
        default=None,
        validator=validate_int_or_None,
        description="[DEPRECATED] This setting has no effect."
    ),
    _Param(
        "text.kerning_factor",
        default=None,
        validator=validate_int_or_None,
        description="[DEPRECATED] Specifies the scaling factor for kerning values. "
                    "This is provided solely to allow old test images to remain "
                    "unchanged. Set to 6 to obtain previous behavior. Values other "
                    "than 0 or 6 have no defined meaning."
    ),
    _Param(
        "text.antialiased",
        default=True,
        validator=validate_bool,
        description="If True (default), the text will be antialiased. This only "
                    "affects raster outputs."
    ),
    _Param(
        "text.parse_math",
        default=True,
        validator=validate_bool,
        description="Use mathtext if there is an even number of unescaped dollar signs."

    ),
    _Section("Mathtext and LaTeX"),
    _Param(
        "text.usetex",
        default=False,
        validator=validate_bool,
        description="use latex for all text handling. The following fonts are "
                    "supported through the usual rc parameter settings: "
                    "new century schoolbook, bookman, times, palatino, zapf chancery, "
                    "charter, serif, sans-serif, helvetica, avant garde, courier, "
                    "monospace, computer modern roman, computer modern sans serif, "
                    "computer modern typewriter"
    ),
    _Param(
        "text.latex.engine",
        default="latex",
        validator=["latex", "latex+dvipng"],
        description=(
            "The TeX engine/format to use.  The following values are supported:\n"
            "- 'latex': The classic TeX engine (the current default).  All backends "
            "render TeX's output by parsing the DVI output into glyphs and boxes and "
            "emitting those one by one.\n"
            "- 'latex+dvipng': The same as 'latex', with the exception that Agg-based "
            "backends rely on dvipng to rasterize TeX's output.  This value was the "
            "default up to Matplotlib 3.10."
        )
    ),
    _Param(
        "text.latex.preamble",
        default="",
        validator=validate_string,
        description='IMPROPER USE OF THIS FEATURE WILL LEAD TO LATEX FAILURES AND IS '
                    'THEREFORE UNSUPPORTED. PLEASE DO NOT ASK FOR HELP IF THIS FEATURE '
                    'DOES NOT DO WHAT YOU EXPECT IT TO. text.latex.preamble is a '
                    'single line of LaTeX code that will be passed on to the LaTeX '
                    'system. It may contain any code that is valid for the LaTeX '
                    '"preamble", i.e. between the "\\documentclass" and '
                    '"\\begin{document}" statements. Note that it has to be put on a '
                    'single line, which may become quite long. The following packages '
                    'are always loaded with usetex, so beware of package collisions: '
                    '   color, fix-cm, geometry, graphicx, textcomp. PostScript '
                    '(PSNFSS) font packages may also be loaded, depending on your font '
                    'settings.'
    ),
    _Param(
        "mathtext.fontset",
        default="dejavusans",
        validator=["dejavusans", "dejavuserif", "cm", "stix", "stixsans", "custom"],
        description="Should be 'dejavusans' (default), 'dejavuserif', "
                    "'cm' (Computer Modern), 'stix', 'stixsans' or 'custom'"
    ),
    _Param("mathtext.bf", "sans:bold", validate_font_properties),
    _Param("mathtext.bfit", "sans:italic:bold", validate_font_properties),
    _Param("mathtext.cal", "cursive", validate_font_properties),
    _Param("mathtext.it", "sans:italic", validate_font_properties),
    _Param("mathtext.rm", "sans", validate_font_properties),
    _Param("mathtext.sf", "sans", validate_font_properties),
    _Param("mathtext.tt", "monospace", validate_font_properties),
    _Param(
        "mathtext.fallback",
        default="cm",
        validator=_validate_mathtext_fallback,
        description="Select fallback font from ['cm' (Computer Modern), 'stix', "
                    "'stixsans'] when a symbol cannot be found in one of the custom "
                    "math fonts. Select 'None' to not perform fallback and replace the "
                    "missing character by a dummy symbol."
    ),
    _Param("mathtext.default", "normal",
           ["rm", "cal", "bfit", "it", "tt", "sf", "bf", "default", "bb", "frak", "scr",
            "regular", "normal"],
           description='The default font to use for math. Can be any of the LaTeX font '
                    'names, including the special name "regular" for the same font '
                    'used in regular text.',
           ),
    _Section("Axes"),
    _Param(
        "axes.facecolor",
        default="white",
        validator=validate_color,
        description="axes background color"
    ),
    _Param(
        "axes.edgecolor",
        default="black",
        validator=validate_color,
        description="axes edge color"
    ),
    _Param(
        "axes.linewidth",
        default=0.8,
        validator=validate_float,
        description="edge line width"
    ),
    _Param(
        "axes.grid",
        default=False,
        validator=validate_bool,
        description="display grid or not"
    ),
    _Param(
        "axes.grid.axis",
        default="both",
        validator=["x", "y", "both"],
        description="which axis the grid should apply to"
    ),
    _Param(
        "axes.grid.which",
        default="major",
        validator=["minor", "both", "major"],
        description="grid lines at {major, minor, both} ticks"
    ),
    _Param(
        "axes.titlelocation",
        default="center",
        validator=["left", "center", "right"],
        description="alignment of the title: {left, right, center}"
    ),
    _Param(
        "axes.titlesize",
        default="large",
        validator=validate_fontsize,
        description="font size of the axes title"
    ),
    _Param(
        "axes.titleweight",
        default="normal",
        validator=validate_fontweight,
        description="font weight of title"
    ),
    _Param(
        "axes.titlecolor",
        default="auto",
        validator=validate_color_or_auto,
        description="color of the axes title, auto falls back to text.color as default "
                    "value"
    ),
    _Param(
        "axes.titley",
        default=None,
        validator=validate_float_or_None,
        description="position title (axes relative units).  None implies auto"
    ),
    _Param(
        "axes.titlepad",
        default=6.0,
        validator=validate_float,
        description="pad between axes and title in points"
    ),
    _Param(
        "axes.labelsize",
        default="medium",
        validator=validate_fontsize,
        description="font size of the x and y labels"
    ),
    _Param(
        "axes.labelpad",
        default=4.0,
        validator=validate_float,
        description="space between label and axis"
    ),
    _Param(
        "axes.labelweight",
        default="normal",
        validator=validate_fontweight,
        description="weight of the x and y labels"
    ),
    _Param(
        "axes.labelcolor",
        default="black",
        validator=validate_color
    ),
    _Param(
        "axes.axisbelow",
        default="line",
        validator=validate_axisbelow,
        description="draw axis gridlines and ticks: "
                    "- below patches (True) "
                    "- above patches but below lines ('line') "
                    "- above all (False)"
    ),
    _Param(
        "axes.formatter.limits",
        default=[-5, 6],
        validator=validate_intlist,
        description="use scientific notation if log10 of the axis range is smaller "
                    "than the first or larger than the second"
    ),
    _Param(
        "axes.formatter.use_locale",
        default=False,
        validator=validate_bool,
        description="When True, format tick labels according to the user's locale. "
                    "For example, use ',' as a decimal separator in the fr_FR locale."
    ),
    _Param(
        "axes.formatter.use_mathtext",
        default=False,
        validator=validate_bool,
        description="When True, use mathtext for scientific notation."
    ),
    _Param(
        "axes.formatter.min_exponent",
        default=0,
        validator=validate_int,
        description="minimum exponent to format in scientific notation"
    ),
    _Param(
        "axes.formatter.useoffset",
        default=True,
        validator=validate_bool,
        description="If True, the tick label formatter will default to labeling ticks "
                    "relative to an offset when the data range is small compared to "
                    "the minimum absolute value of the data."
    ),
    _Param(
        "axes.formatter.offset_threshold",
        default=4,
        validator=validate_int,
        description="When useoffset is True, the offset will be used when it can "
                    "remove at least this number of significant digits from tick "
                    "labels."
    ),
    _Param(
        "axes.spines.left",
        default=True,
        validator=validate_bool,
        description="display axis spines"
    ),
    _Param("axes.spines.bottom", True, validate_bool),
    _Param("axes.spines.top", True, validate_bool),
    _Param(
        "axes.spines.right",
        default=True,
        validator=validate_bool
    ),
    _Param(
        "axes.unicode_minus",
        default=True,
        validator=validate_bool,
        description="use Unicode for the minus symbol rather than hyphen. See "
                    "https://en.wikipedia.org/wiki/Plus_and_minus_signs#Character_codes"

    ),
    _Param("axes.prop_cycle",
        default=cycler(
            "color",
            [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
             (1.0, 0.4980392156862745, 0.054901960784313725),
             (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
             (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
             (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
             (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
             (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
             (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
             (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
             (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
            ],
        ),
        validator=validate_cycler
    ),
    _Param(
        "axes.xmargin",
        default=0.05,
        validator=_validate_greaterthan_minushalf,
        description="x margin.  See `~.axes.Axes.margins`"
    ),
    _Param(
        "axes.ymargin",
        default=0.05,
        validator=_validate_greaterthan_minushalf,
        description="y margin.  See `~.axes.Axes.margins`"
    ),
    _Param(
        "axes.zmargin",
        default=0.05,
        validator=_validate_greaterthan_minushalf,
        description="z margin.  See `~.axes.Axes.margins`"
    ),
    _Param(
        "axes.autolimit_mode",
        default="data",
        validator=["data", "round_numbers"],
        description='If "data", use axes.xmargin and axes.ymargin as is. If '
                    '"round_numbers", after application of margins, axis limits are '
                    'further expanded to the nearest "round" number.',
    ),
    _Subsection("Polar Axes"),
    _Param(
        "polaraxes.grid",
        default=True,
        validator=validate_bool,
        description="display grid on polar axes"
    ),
    _Subsection("3D Axes"),
    _Param(
        "axes3d.grid",
        default=True,
        validator=validate_bool, description="display grid on 3D axes"
    ),
    _Param(
        "axes3d.automargin",
        default=False,
        validator=validate_bool,
        description="automatically add margin when manually setting 3D axis limits"
    ),
    _Param(
        "axes3d.xaxis.panecolor",
        default=(0.95, 0.95, 0.95, 0.5),
        validator=validate_color,
        description="background pane on 3D axes"
    ),
    _Param(
        "axes3d.yaxis.panecolor",
        default=(0.9, 0.9, 0.9, 0.5),
        validator=validate_color,
        description="background pane on 3D axes"
    ),
    _Param(
        "axes3d.zaxis.panecolor",
        default=(0.925, 0.925, 0.925, 0.5),
        validator=validate_color,
        description="background pane on 3D axes"
    ),
    _Param(
        "axes3d.depthshade",
        default=True,
        validator=validate_bool,
        description="depth shade for 3D scatter plots"
    ),
    _Param(
        "axes3d.depthshade_minalpha",
        default=0.3,
        validator=validate_float,
        description="minimum alpha value for depth shading"
    ),
    _Param(
        "axes3d.mouserotationstyle",
        default="arcball",
        validator=["azel", "trackball", "sphere", "arcball"],
        description="{azel, trackball, sphere, arcball} See also "
                    "https://matplotlib.org/stable/api/toolkits/mplot3d/view_angles.html#rotation-with-mouse"),  # noqa
    _Param(
        "axes3d.trackballsize",
        default=0.667,
        validator=validate_float,
        description="trackball diameter, in units of the Axes bbox"
    ),
    _Param(
        "axes3d.trackballborder",
        default=0.2,
        validator=validate_float,
        description="trackball border width, in units of the Axes bbox (only for "
                    "'sphere' and 'arcball' style)"
    ),
    _Section("Axis"),
    _Param(
        "axes3d.snap_rotation",
        default=5.0,
        validator=validate_float,
        description="Snap angle (in degrees) for 3D rotation when holding Control."
   ),
    _Param(
        "xaxis.labellocation",
        default="center",
        validator=["left", "center", "right"],
        description="alignment of the xaxis label: {left, right, center}"
    ),
    _Param(
        "yaxis.labellocation",
        default="center",
        validator=["bottom", "center", "top"],
        description="alignment of the yaxis label: {bottom, top, center}"
    ),
    _Section(
        "Dates",
        description="Default properties for date tick labels. These are used by the "
                    "`.AutoDateFormatter` when the appropriate time unit is detected."
                    "See "
                    "https://matplotlib.org/stable/api/dates_api.html#date-formatters "
                    "for more information."
    ),
    _Param("date.autoformatter.year", "%Y", validate_string),
    _Param("date.autoformatter.month", "%Y-%m", validate_string),
    _Param("date.autoformatter.day", "%Y-%m-%d", validate_string),
    _Param("date.autoformatter.hour", "%m-%d %H", validate_string),
    _Param("date.autoformatter.minute", "%d %H:%M", validate_string),
    _Param("date.autoformatter.second", "%H:%M:%S", validate_string),
    _Param("date.autoformatter.microsecond", "%M:%S.%f", validate_string),
    _Param(
        "date.epoch",
        default="1970-01-01T00:00:00",
        validator=_validate_date,
        description="The reference date for Matplotlib's internal date representation. "
                    "See https://matplotlib.org/stable/gallery/ticks/date_precision_and_epochs.html"),  #noqa
    _Param(
        "date.converter",
        default="auto",
        validator=["auto", "concise"],
        description="'auto', 'concise'"
    ),
    _Param(
        "date.interval_multiples",
        default=True,
        validator=validate_bool,
        description="For auto converter whether to use interval_multiples"
    ),
    _Section("Ticks"),
    _Param(
        "xtick.top",
        default=False,
        validator=validate_bool,
        description="draw ticks on the top side"
    ),
    _Param(
        "xtick.bottom",
        default=True,
        validator=validate_bool,
        description="draw ticks on the bottom side"
    ),
    _Param(
        "xtick.labeltop",
        default=False,
        validator=validate_bool,
        description="draw label on the top"
    ),
    _Param(
        "xtick.labelbottom",
        default=True,
        validator=validate_bool,
        description="draw label on the bottom"
    ),
    _Param(
        "xtick.major.size",
        default=3.5,
        validator=validate_float,
        description="major tick size in points"
    ),
    _Param(
        "xtick.minor.size",
        default=2.0,
        validator=validate_float,
        description="minor tick size in points"
    ),
    _Param(
        "xtick.major.width",
        default=0.8,
        validator=validate_float,
        description="major tick width in points"
    ),
    _Param(
        "xtick.minor.width",
        default=0.6,
        validator=validate_float,
        description="minor tick width in points"
    ),
    _Param(
        "xtick.major.pad",
        default=3.5,
        validator=validate_float,
        description="distance to major tick label in points"
    ),
    _Param(
        "xtick.minor.pad",
        default=3.4,
        validator=validate_float,
        description="distance to the minor tick label in points"
    ),
    _Param(
        "xtick.color",
        default="black",
        validator=validate_color,
        description="color of the ticks"
    ),
    _Param(
        "xtick.labelcolor",
        default="inherit",
        validator=validate_color_or_inherit,
        description="color of the tick labels or inherit from xtick.color"
    ),
    _Param(
        "xtick.labelsize",
        default="medium",
        validator=validate_fontsize,
        description="font size of the tick labels"
    ),
    _Param(
        "xtick.direction",
        default="out",
        validator=["out", "in", "inout"],
        description="direction: {in, out, inout}"
    ),
    _Param(
        "xtick.minor.visible",
        default=False,
        validator=validate_bool,
        description="visibility of minor ticks on x-axis"
    ),
    _Param(
        "xtick.major.top",
        default=True,
        validator=validate_bool,
        description="draw x axis top major ticks"
    ),
    _Param(
        "xtick.major.bottom",
        default=True,
        validator=validate_bool,
        description="draw x axis bottom major ticks"
    ),
    _Param(
        "xtick.minor.top",
        default=True,
        validator=validate_bool,
        description="draw x axis top minor ticks"
    ),
    _Param(
        "xtick.minor.bottom",
        default=True,
        validator=validate_bool,
        description="draw x axis bottom minor ticks"
    ),
    _Param(
        "xtick.minor.ndivs",
        default="auto",
        validator=_validate_minor_tick_ndivs,
        description="number of minor ticks between the major ticks on x-axis"
    ),
    _Param(
        "xtick.alignment",
        default="center",
        validator=["center", "right", "left"],
        description="alignment of xticks"
    ),
    _Param(
        "ytick.left",
        default=True,
        validator=validate_bool,
        description="draw ticks on the left side"
    ),
    _Param(
        "ytick.right",
        default=False,
        validator=validate_bool,
        description="draw ticks on the right side"
    ),
    _Param(
        "ytick.labelleft",
        default=True,
        validator=validate_bool,
        description="draw tick labels on the left side"
    ),
    _Param(
        "ytick.labelright",
        default=False,
        validator=validate_bool,
        description="draw tick labels on the right side"
    ),
    _Param(
        "ytick.major.size",
        default=3.5,
        validator=validate_float,
        description="major tick size in points"
    ),
    _Param(
        "ytick.minor.size",
        default=2.0,
        validator=validate_float,
        description="minor tick size in points"
    ),
    _Param(
        "ytick.major.width",
        default=0.8,
        validator=validate_float,
        description="major tick width in points"
    ),
    _Param(
        "ytick.minor.width",
        default=0.6,
        validator=validate_float,
        description="minor tick width in points"
    ),
    _Param(
        "ytick.major.pad",
        default=3.5,
        validator=validate_float,
        description="distance to major tick label in points"
    ),
    _Param(
        "ytick.minor.pad",
        default=3.4,
        validator=validate_float,
        description="distance to the minor tick label in points"
    ),
    _Param(
        "ytick.color",
        default="black",
        validator=validate_color,
        description="color of the ticks"
    ),
    _Param(
        "ytick.labelcolor",
        default="inherit",
        validator=validate_color_or_inherit,
        description="color of the tick labels or inherit from ytick.color"
    ),
    _Param(
        "ytick.labelsize",
        default="medium",
        validator=validate_fontsize,
        description="font size of the tick labels"
    ),
    _Param(
        "ytick.direction",
        default="out",
        validator=["out", "in", "inout"],
        description="direction: {in, out, inout}"
    ),
    _Param(
        "ytick.minor.visible",
        default=False,
        validator=validate_bool,
        description="visibility of minor ticks on y-axis"
    ),
    _Param(
        "ytick.major.left",
        default=True,
        validator=validate_bool,
        description="draw y axis left major ticks"
    ),
    _Param(
        "ytick.major.right",
        default=True,
        validator=validate_bool,
        description="draw y axis right major ticks"
    ),
    _Param(
        "ytick.minor.left",
        default=True,
        validator=validate_bool,
        description="draw y axis left minor ticks"
    ),
    _Param(
        "ytick.minor.right",
        default=True,
        validator=validate_bool,
        description="draw y axis right minor ticks"
    ),
    _Param(
        "ytick.minor.ndivs",
        default="auto",
        validator=_validate_minor_tick_ndivs,
        description="number of minor ticks between the major ticks on y-axis"
    ),
    _Param("ytick.alignment", "center_baseline",
           ["center", "top", "bottom", "baseline", "center_baseline"],
           description="alignment of yticks"
           ),
    _Section("Grid"),
    _Param(
        "grid.color",
        default="#b0b0b0",
        validator=validate_color,
        description='b0b0b0"  # grid color'
    ),
    _Param(
        "grid.linestyle",
        default="-",
        validator=_validate_linestyle,
        description="solid"
    ),
    _Param(
        "grid.linewidth",
        default=0.8,
        validator=validate_float,
        description="in points"
    ),
    _Param(
        "grid.alpha",
        default=1.0,
        validator=validate_float,
        description="transparency, between 0.0 and 1.0"
    ),
    _Param(
        "grid.major.color",
        default=None,
        validator=_validate_color_or_None,
        description="If None defaults to grid.color"
    ),
    _Param(
        "grid.major.linestyle",
        default=None,
        validator=_validate_linestyle_or_None,
        description="If None defaults to grid.linestyle"
    ),
    _Param(
        "grid.major.linewidth",
        default=None,
        validator=validate_float_or_None,
        description="If None defaults to grid.linewidth"
    ),
    _Param(
        "grid.major.alpha",
        default=None,
        validator=validate_float_or_None,
        description="If None defaults to grid.alpha"
    ),
    _Param(
        "grid.minor.color",
        default=None,
        validator=_validate_color_or_None,
        description="If None defaults to grid.color"
    ),
    _Param(
        "grid.minor.linestyle",
        default=None,
        validator=_validate_linestyle_or_None,
        description="If None defaults to grid.linestyle"
    ),
    _Param(
        "grid.minor.linewidth",
        default=None,
        validator=validate_float_or_None,
        description="If None defaults to grid.linewidth"
    ),
    _Param(
        "grid.minor.alpha",
        default=None,
        validator=validate_float_or_None,
        description="If None defaults to grid.alpha"
    ),
    _Section("Legend"),
    _Param(
        "legend.loc",
        default="best",
        validator=_validate_legend_loc
    ),
    _Param(
        "legend.frameon",
        default=True,
        validator=validate_bool,
        description="if True, draw the legend on a background patch"
    ),
    _Param(
        "legend.framealpha",
        default=0.8,
        validator=validate_float_or_None,
        description="legend patch transparency"
    ),
    _Param(
        "legend.facecolor",
        default="inherit",
        validator=validate_color_or_inherit,
        description="inherit from axes.facecolor; or color spec"
    ),
    _Param(
        "legend.edgecolor",
        default="0.8",
        validator=validate_color_or_inherit,
        description="background patch boundary color"
    ),
    _Param(
        "legend.linewidth",
        default=None,
        validator=validate_float_or_None,
        description="line width of the legend frame, None means inherit from "
                    "patch.linewidth"
    ),
    _Param(
        "legend.fancybox",
        default=True,
        validator=validate_bool,
        description="if True, use a rounded box for the legend background, else a "
                    "rectangle"
    ),
    _Param(
        "legend.shadow",
        default=False,
        validator=validate_bool,
        description="if True, give background a shadow effect"
    ),
    _Param(
        "legend.numpoints",
        default=1,
        validator=validate_int,
        description="the number of marker points in the legend line"
    ),
    _Param(
        "legend.scatterpoints",
        default=1,
        validator=validate_int,
        description="number of scatter points"
    ),
    _Param(
        "legend.markerscale",
        default=1.0,
        validator=validate_float,
        description="the relative size of legend markers vs. original"
    ),
    _Param(
        "legend.fontsize",
        default="medium",
        validator=validate_fontsize
    ),
    _Param(
        "legend.labelcolor",
        default="None",
        validator=_validate_color_or_linecolor
    ),
    _Param(
        "legend.title_fontsize",
        default=None,
        validator=validate_fontsize_None,
        description="None sets to the same as the default axes."
    ),
    _Param(
        "legend.borderpad",
        default=0.4,
        validator=validate_float,
        description="border whitespace"
    ),
    _Param(
        "legend.labelspacing",
        default=0.5,
        validator=validate_float,
        description="the vertical space between the legend entries"
    ),
    _Param(
        "legend.handlelength",
        default=2.0,
        validator=validate_float,
        description="the length of the legend lines"
    ),
    _Param(
        "legend.handleheight",
        default=0.7,
        validator=validate_float,
        description="the height of the legend handle"
    ),
    _Param(
        "legend.handletextpad",
        default=0.8,
        validator=validate_float,
        description="the space between the legend line and legend text"
    ),
    _Param(
        "legend.borderaxespad",
        default=0.5,
        validator=validate_float,
        description="the border between the axes and legend edge"
    ),
    _Param(
        "legend.columnspacing",
        default=2.0,
        validator=validate_float, description="column separation"
    ),
    _Section("Figure"),
    _Param(
        "figure.titlesize",
        default="large",
        validator=validate_fontsize,
        description="size of the figure title (``Figure.suptitle()``)"
    ),
    _Param(
        "figure.titleweight",
        default="normal",
        validator=validate_fontweight,
        description="weight of the figure title"
    ),
    _Param(
        "figure.labelsize",
        default="large",
        validator=validate_fontsize,
        description="size of the figure label (``Figure.sup[x|y]label()``)"
    ),
    _Param(
        "figure.labelweight",
        default="normal",
        validator=validate_fontweight,
        description="weight of the figure label"
    ),
    _Param(
        "figure.figsize",
        default=[6.4, 4.8],
        validator=_listify_validator(validate_float, n=2),
        description="figure size in inches"
    ),
    _Param(
        "figure.dpi",
        default=100.0,
        validator=validate_float, description="figure dots per inch"
    ),
    _Param(
        "figure.facecolor",
        default="white",
        validator=validate_color, description="figure face color"
    ),
    _Param(
        "figure.edgecolor",
        default="white",
        validator=validate_color, description="figure edge color"
    ),
    _Param(
        "figure.frameon",
        default=True,
        validator=validate_bool, description="enable figure frame"
    ),
    _Param(
        "figure.max_open_warning",
        default=20,
        validator=validate_int,
        description="The maximum number of figures to open through the pyplot "
                    "interface before emitting a warning. If less than one this "
                    "feature is disabled."
    ),
    _Param(
        "figure.raise_window",
        default=True,
        validator=validate_bool,
        description="Raise the GUI window to front when show() is called. If set to "
                    "False, we currently do not take any further actions and whether "
                    "the window appears on the front may depend on the GUI framework "
                    "and window manager."
    ),
    _Param(
        "figure.subplot.left",
        default=0.125,
        validator=validate_float,
        description="the left side of the subplots of the figure"
    ),
    _Param(
        "figure.subplot.right",
        default=0.9,
        validator=validate_float,
        description="the right side of the subplots of the figure"
    ),
    _Param(
        "figure.subplot.bottom",
        default=0.11,
        validator=validate_float,
        description="the bottom of the subplots of the figure"
    ),
    _Param(
        "figure.subplot.top",
        default=0.88,
        validator=validate_float,
        description="the top of the subplots of the figure"
    ),
    _Param(
        "figure.subplot.wspace",
        default=0.2,
        validator=validate_float,
        description="the amount of width reserved for space between subplots, "
                    "expressed as a fraction of the average axis width"
    ),
    _Param(
        "figure.subplot.hspace",
        default=0.2,
        validator=validate_float,
        description="the amount of height reserved for space between subplots, "
                    "expressed as a fraction of the average axis height"
    ),
    _Param(
        "figure.autolayout",
        default=False,
        validator=validate_bool,
        description="When True, automatically adjust subplot parameters to make the "
                    "plot fit the figure using `~.Figure.tight_layout`"
    ),
    _Param(
        "figure.constrained_layout.use",
        default=False,
        validator=validate_bool,
        description="When True, automatically make plot elements fit on the figure. "
                    '(Not compatible with "figure.autolayout", above).'
    ),
    _Param(
        "figure.constrained_layout.h_pad",
        default=0.04167,
        validator=validate_float,
        description="Padding (in inches) around axes; defaults to 3/72 inches, "
                    "i.e. 3 points"
    ),
    _Param(
        "figure.constrained_layout.w_pad",
        default=0.04167,
        validator=validate_float,
        description="Padding (in inches) around axes; defaults to 3/72 inches, "
                    "i.e. 3 points"
    ),
    _Param(
        "figure.constrained_layout.hspace",
        default=0.02,
        validator=validate_float,
        description="Spacing between subplots, relative to the subplot sizes.  Much "
                    "smaller than for tight_layout (figure.subplot.hspace, "
                    "figure.subplot.wspace) as constrained_layout already takes "
                    "surrounding texts (titles, labels, # ticklabels) into account."
    ),
    _Param(
        "figure.constrained_layout.wspace",
        default=0.02,
        validator=validate_float,
        description="Spacing between subplots, relative to the subplot sizes.  Much "
                    "smaller than for tight_layout (figure.subplot.hspace, "
                    "figure.subplot.wspace) as constrained_layout already takes "
                    "surrounding texts (titles, labels, # ticklabels) into account."
    ),
    _Section("Images"),
    _Param(
        "image.aspect",
        default="equal",
        validator=validate_aspect,
        description="{equal, auto} or a number"
    ),
    _Param(
        "image.interpolation",
        default="auto",
        validator=validate_string,
        description="see help(imshow) for options"
    ),
    _Param(
        "image.interpolation_stage",
        default="auto",
        validator=["auto", "data", "rgba"],
        description="see help(imshow) for options"
    ),
    _Param(
        "image.cmap",
        default="viridis",
        validator=_validate_cmap,
        description="A colormap name (plasma, magma, etc.)"
    ),
    _Param(
        "image.lut",
        default=256,
        validator=validate_int,
        description="the size of the colormap lookup table"
    ),
    _Param(
        "image.origin",
        default="upper",
        validator=["upper", "lower"], description="{lower, upper}"
    ),
    _Param(
        "image.resample",
        default=True,
        validator=validate_bool
    ),
    _Param(
        "image.composite_image",
        default=True,
        validator=validate_bool,
        description="When True, all the images on a set of axes are combined into a "
                    "single composite image before saving a figure as a vector "
                    "graphics file, such as a PDF."
    ),
    _Section("Contour plots"),
    _Param(
        "contour.negative_linestyle",
        default="dashed",
        validator=_validate_linestyle,
        description="string or on-off ink sequence"
    ),
    _Param(
        "contour.corner_mask",
        default=True,
        validator=validate_bool, description="{True, False}"
    ),
    _Param(
        "contour.linewidth",
        default=None,
        validator=validate_float_or_None,
        description="{float, None} Size of the contour line widths. If set to None, it "
                    'falls back to "line.linewidth".'
    ),
    _Param(
        "contour.algorithm",
        default="mpl2014",
        validator=["mpl2005", "mpl2014", "serial", "threaded"],
        description="{mpl2005, mpl2014, serial, threaded}"
    ),
    _Section("Errorbar plots"),
    _Param(
        "errorbar.capsize",
        default=0.0,
        validator=validate_float,
        description="length of end cap on error bars in pixels"
    ),
    _Param(
        "errorbar.capthick",
        default=None,
        validator=validate_float_or_None,
        description="thickness of end cap on error bars in points."),
    _Param(
        "errorbar.elinewidth",
        default=None,
        validator=validate_float_or_None,
        description="line width of the error bar lines in points."
    ),
    _Section("Histogram plots"),
    _Param(
        "hist.bins",
        default=10,
        validator=validate_hist_bins,
        description="The default number of histogram bins or 'auto'."
    ),
    _Section("Scatter plots"),
    _Param(
        "scatter.marker",
        default="o",
        validator=_validate_marker,
        description="The default marker type for scatter plots."
    ),
    _Param(
        "scatter.edgecolors",
        default="face",
        validator=validate_string,
        description="The default edge colors for scatter plots."
    ),
    _Section("AGG rendering"),
    _Param(
        "agg.path.chunksize",
        default=0,
        validator=validate_int,
        description="0 to disable; values in the range 10000 to 100000 can improve "
                    "speed slightly and prevent an Agg rendering failure when plotting "
                    "very large data sets, especially if they are very gappy. It may "
                    "cause minor artifacts, though. A value of 20000 is probably a "
                    "good starting point."
    ),
    _Section("Paths"),
    _Param(
        "path.simplify",
        default=True,
        validator=validate_bool,
        description='When True, simplify paths by removing "invisible" points to '
                    'reduce file size and increase rendering speed',
    ),
    _Param(
        "path.simplify_threshold",
        default=0.111111111111,
        validator=_validate_greaterequal0_lessequal1,
        description="The threshold of similarity below which vertices will be removed "
                    "in the simplification process."
    ),
    _Param(
        "path.snap",
        default=True,
        validator=validate_bool,
        description="When True, rectilinear axis-aligned paths will be snapped to the "
                    "nearest pixel when certain criteria are met. When False, paths "
                    "will never be snapped."
    ),
    _Param(
        "path.sketch",
        default=None,
        validator=validate_sketch,
        description="May be None, or a tuple of the form:"
                    "path.sketch: (scale, length, randomness)"
                    "- *scale* is the amplitude of the wiggle perpendicular to the line"
                    "  (in pixels)."
                    "- *length* is the length of the wiggle along the line (in pixels)."
                    "- *randomness* is the factor by which the length is  randomly "
                    "  scaled."
    ),
    _Param(
        "path.effects",
        default=[],
        validator=validate_anylist
    ),
    _Section("Saving figures"),
    _Param(
        "savefig.dpi",
        default="figure",
        validator=validate_dpi,
        description="figure dots per inch or 'figure'"
    ),
    _Param(
        "savefig.facecolor",
        default="auto",
        validator=validate_color_or_auto,
        description="figure face color when saving"
    ),
    _Param(
        "savefig.edgecolor",
        default="auto",
        validator=validate_color_or_auto,
        description="figure edge color when saving"
    ),
    _Param(
        "savefig.format",
        default="png",
        validator=validate_string, description="{png, ps, pdf, svg}"
    ),
    _Param(
        "savefig.bbox",
        default=None,
        validator=validate_bbox,
        description="{tight, standard} 'tight' is incompatible with generating frames "
                    "for animation"
    ),
    _Param(
        "savefig.pad_inches",
        default=0.1,
        validator=validate_float,
        description="padding to be used, when bbox is set to 'tight'"
    ),
    _Param(
        "savefig.directory",
        default="~",
        validator=_validate_pathlike,
        description="default directory in savefig dialog, gets updated after "
                    "interactive saves, unless set to the empty string (i.e. the "
                    "current directory); use '.' to start at the current directory but "
                    "update after interactive saves"
    ),
    _Param(
        "savefig.transparent",
        default=False,
        validator=validate_bool,
        description="whether figures are saved with a transparent background by default"

    ),
    _Param(
        "savefig.orientation",
        default="portrait",
        validator=["landscape", "portrait"],
        description="orientation of saved figure, for PostScript output only"
    ),
    _Subsection("Mac OSX backend parameters"),
    _Param(
        "macosx.window_mode",
        default="system",
        validator=["system", "tab", "window"],
        description="How to open new figures (system, tab, window) system uses "
                    "the MacOS system preferences"
    ),
    _Subsection("Tk backend parameters"),
    _Param(
        "tk.window_focus",
        default=False,
        validator=validate_bool,
        description="Maintain shell focus for TkAgg"
    ),
    _Subsection("PS backend parameters"),
    _Param(
        "ps.papersize",
        default="letter",
        validator=_ignorecase(
            ["figure", "letter", "legal", "ledger",
             *[f"{ab}{i}" for ab in "ab" for i in range(11)],
             ],
        ),
        description="{figure, letter, legal, ledger, A0-A10, B0-B10}"
    ),
    _Param(
        "ps.useafm",
        default=False,
        validator=validate_bool,
        description="use AFM fonts, results in small files"
    ),
    _Param(
        "ps.usedistiller",
        default=None,
        validator=validate_ps_distiller,
        description="{ghostscript, xpdf, None} Experimental: may produce smaller "
                    "files. xpdf intended for production of publication quality files, "
                    "but requires ghostscript, xpdf and ps2eps"
    ),
    _Param(
        "ps.distiller.res",
        default=6000,
        validator=validate_int, description="dpi"
    ),
    _Param(
        "ps.fonttype",
        default=3,
        validator=validate_fonttype,
        description="Output Type 3 (Type3) or Type 42 (TrueType)"
    ),
    _Subsection("PDF backend parameters"),
    _Param(
        "pdf.compression",
        default=6,
        validator=validate_int,
        description="integer from 0 to 9 0 disables compression (good for debugging)"
    ),
    _Param(
        "pdf.fonttype",
        default=3,
        validator=validate_fonttype,
        description="Output Type 3 (Type3) or Type 42 (TrueType)"
    ),
    _Param(
        "pdf.use14corefonts",
        default=False,
        validator=validate_bool
    ),
    _Param(
        "pdf.inheritcolor",
        default=False,
        validator=validate_bool
    ),
    _Subsection("SVG backend parameters"),
    _Param(
        "svg.image_inline",
        default=True,
        validator=validate_bool,
        description="Write raster image data directly into the SVG file"
    ),
    _Param(
        "svg.fonttype",
        default="path",
        validator=["none", "path"],
        description="How to handle SVG fonts: "
                    "path: Embed characters as paths -- supported by most SVG "
                    "      renderers"
                    "none: Assume fonts are installed on the machine where the SVG "
                    "will be viewed."
    ),
    _Param(
        "svg.hashsalt",
        default=None,
        validator=validate_string_or_None,
        description="If not None, use this string as hash salt instead of uuid4"
    ),
    _Param(
        "svg.id",
        default=None,
        validator=validate_string_or_None,
        description="If not None, use this string as the value for the `id` attribute "
                    "in the top <svg> tag"
    ),
    _Subsection("PGF parameters"),
    _Param(
        "pgf.rcfonts",
        default=True,
        validator=validate_bool
    ),
    _Param(
        "pgf.preamble",
        default="",
        validator=validate_string,
        description="See text.latex.preamble for documentation"
    ),
    _Param(
        "pgf.texsystem",
        default="xelatex",
        validator=["xelatex", "lualatex", "pdflatex"]
    ),
    _Subsection("Docstring parameters"),
    _Param(
        "docstring.hardcopy",
        default=False,
        validator=validate_bool,
        description="set this when you want to generate hardcopy docstring"
    ),
    _Section(
        "Interactive keymaps",
        description="Default key mappings for interactive navigation. See "
                    ":ref:`key-event-handling`."
    ),
    _Param(
        "keymap.fullscreen",
        default=["f", "ctrl+f"],
        validator=validate_stringlist,
        description="toggling"
    ),
    _Param(
        "keymap.home",
        default=["h", "r", "home"],
        validator=validate_stringlist,
        description="home or reset mnemonic"
    ),
    _Param(
        "keymap.back",
        default=["left", "c", "backspace", "MouseButton.BACK"],
        validator=validate_stringlist, description="forward / backward keys"
    ),
    _Param(
        "keymap.forward",
        default=["right", "v", "MouseButton.FORWARD"],
        validator=validate_stringlist,
        description="for quick navigation"
    ),
    _Param(
        "keymap.pan",
        default=["p"],
        validator=validate_stringlist, description="pan mnemonic"
    ),
    _Param(
        "keymap.zoom",
        default=["o"],
        validator=validate_stringlist, description="zoom mnemonic"
    ),
    _Param(
        "keymap.save",
        default=["s", "ctrl+s"],
        validator=validate_stringlist,
        description="saving current figure"
    ),
    _Param(
        "keymap.help",
        default=["f1"],
        validator=validate_stringlist,
        description="display help about active tools"
    ),
    _Param(
        "keymap.quit",
        default=["ctrl+w", "cmd+w", "q"],
        validator=validate_stringlist,
        description="close the current figure"
    ),
    _Param(
        "keymap.quit_all",
        default=[],
        validator=validate_stringlist, description="close all figures"
    ),
    _Param(
        "keymap.grid",
        default=["g"],
        validator=validate_stringlist,
        description="switching on/off major grids in current axes"
    ),
    _Param(
        "keymap.grid_minor",
        default=["G"],
        validator=validate_stringlist,
        description="switching on/off minor grids in current axes"
    ),
    _Param(
        "keymap.yscale",
        default=["l"],
        validator=validate_stringlist,
        description="toggle scaling of y-axes ('log'/'linear')"
    ),
    _Param(
        "keymap.xscale",
        default=["k", "L"],
        validator=validate_stringlist,
        description="toggle scaling of x-axes ('log'/'linear')"
    ),
    _Param(
        "keymap.copy",
        default=["ctrl+c", "cmd+c"],
        validator=validate_stringlist,
        description="copy figure to clipboard"
    ),
    _Section("Animation"),
    _Param(
        "animation.html",
        default="none",
        validator=["html5", "jshtml", "none"],
        description="How to display the animation as HTML in the IPython notebook: "
                    "- 'html5' uses HTML5 video tag "
                    "- 'jshtml' creates a JavaScript animation"
    ),
    _Param(
        "animation.writer",
        default="ffmpeg",
        validator=validate_string,
        description="MovieWriter 'backend' to use"
    ),
    _Param(
        "animation.codec",
        default="h264",
        validator=validate_string,
        description="Codec to use for writing movie"
    ),
    _Param(
        "animation.bitrate",
        default=-1,
        validator=validate_int,
        description="Controls size/quality trade-off for movie. -1 implies let "
                    "utility auto-determine"
    ),
    _Param("animation.frame_format", "png",
           ["png", "jpeg", "tiff", "raw", "rgba", "ppm", "sgi", "bmp", "pbm", "svg"],
           description="Controls frame format used by temp files"
           ),
    _Param(
        "animation.ffmpeg_path",
        default="ffmpeg",
        validator=_validate_pathlike,
        description="Path to ffmpeg binary.  Unqualified paths are resolved by "
                    "subprocess.Popen."
    ),
    _Param(
        "animation.ffmpeg_args",
        default=[],
        validator=validate_stringlist,
        description="Additional arguments to pass to ffmpeg"
    ),
    _Param(
        "animation.convert_path",
        default="convert",
        validator=_validate_pathlike,
        description="Path to ImageMagick's convert binary.  Unqualified paths are "
                    "resolved by subprocess.Popen, except that on Windows, we look up "
                    "an install of ImageMagick in the registry (as convert is also the "
                    "name of a system tool)."
    ),
    _Param(
        "animation.convert_args",
        default=["-layers", "OptimizePlus"],
        validator=validate_stringlist,
        description="Additional arguments to pass to convert"
    ),
    _Param(
        "animation.embed_limit",
        default=20.0,
        validator=validate_float,
        description="Limit, in MB, of size of base64 encoded animation in HTML (i.e. "
                    "IPython notebook)"
    ),
    _Param(
        "_internal.classic_mode",
        default=False,
        validator=validate_bool
    ),
    _Param("backend", None, validate_backend),
]


def _params_list():
    return [elem for elem in _DEFINITION if isinstance(elem, _Param)]
