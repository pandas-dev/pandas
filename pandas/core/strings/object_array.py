import re
import textwrap
from typing import Pattern, Union
import unicodedata
import warnings

import numpy as np

import pandas._libs.lib as lib
import pandas._libs.missing as libmissing
import pandas._libs.ops as libops
from pandas._typing import Scalar

from pandas.core.dtypes.common import is_re, is_scalar
from pandas.core.dtypes.missing import isna

from pandas.core.accessor import CachedAccessor
from pandas.core.arrays.numpy_ import PandasArray
from pandas.core.strings.base import BaseStringArrayMethods


class ObjectArrayMethods(BaseStringArrayMethods):
    def _map(self, f, na_value=None, dtype=None):
        arr = self._array  # object-dtype ndarray.
        if dtype is None:
            dtype = np.dtype("object")
        if na_value is None:
            na_value = np.nan

        if not len(arr):
            return np.ndarray(0, dtype=dtype)
        if na_value is None:
            na_value = np.nan

        if not isinstance(arr, np.ndarray):
            arr = np.asarray(arr, dtype=object)
        mask = isna(arr)
        convert = not np.all(mask)
        try:
            result = lib.map_infer_mask(arr, f, mask.view(np.uint8), convert)
        except (TypeError, AttributeError) as e:
            # Reraise the exception if callable `f` got wrong number of args.
            # The user may want to be warned by this, instead of getting NaN
            p_err = (
                r"((takes)|(missing)) (?(2)from \d+ to )?\d+ "
                r"(?(3)required )positional arguments?"
            )

            if len(e.args) >= 1 and re.search(p_err, e.args[0]):
                # FIXME: this should be totally avoidable
                raise e

            def g(x):
                # This type of fallback behavior can be removed once
                # we remove object-dtype .str accessor.
                try:
                    return f(x)
                except (TypeError, AttributeError):
                    return na_value

            return self._map(g, na_value=na_value, dtype=dtype)
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.slice(start=key.start, stop=key.stop, step=key.step)
        else:
            return self.get(key)

    def count(self, pat, flags=0):
        regex = re.compile(pat, flags=flags)
        f = lambda x: len(regex.findall(x))
        return self._map(f, dtype="int64")

    def pad(self, width, side="left", fillchar=" "):
        if side == "left":
            f = lambda x: x.rjust(width, fillchar)
        elif side == "right":
            f = lambda x: x.ljust(width, fillchar)
        elif side == "both":
            f = lambda x: x.center(width, fillchar)
        else:  # pragma: no cover
            raise ValueError("Invalid side")
        return self._map(f)

    def contains(self, pat, case=True, flags=0, na=np.nan, regex=True):
        if regex:
            if not case:
                flags |= re.IGNORECASE

            regex = re.compile(pat, flags=flags)

            if regex.groups > 0:
                warnings.warn(
                    "This pattern has match groups. To actually get the "
                    "groups, use str.extract.",
                    UserWarning,
                    stacklevel=3,
                )

            f = lambda x: regex.search(x) is not None
        else:
            if case:
                f = lambda x: pat in x
            else:
                upper_pat = pat.upper()
                f = lambda x: upper_pat in x.upper()
        return self._map(f, na, dtype=np.dtype("bool"))

    def startswith(self, pat, na=np.nan):
        f = lambda x: x.startswith(pat)
        return self._map(f, na_value=na, dtype=np.dtype(bool))

    def endswith(self, pat, na=np.nan):
        f = lambda x: x.endswith(pat)
        return self._map(f, na_value=na, dtype=np.dtype(bool))

    def replace(self, pat, repl, n=-1, case=None, flags=0, regex=True):
        # Check whether repl is valid (GH 13438, GH 15055)
        if not (isinstance(repl, str) or callable(repl)):
            raise TypeError("repl must be a string or callable")

        is_compiled_re = is_re(pat)
        if regex:
            if is_compiled_re:
                if (case is not None) or (flags != 0):
                    raise ValueError(
                        "case and flags cannot be set when pat is a compiled regex"
                    )
            else:
                # not a compiled regex
                # set default case
                if case is None:
                    case = True

                # add case flag, if provided
                if case is False:
                    flags |= re.IGNORECASE
            if is_compiled_re or len(pat) > 1 or flags or callable(repl):
                n = n if n >= 0 else 0
                compiled = re.compile(pat, flags=flags)
                f = lambda x: compiled.sub(repl=repl, string=x, count=n)
            else:
                f = lambda x: x.replace(pat, repl, n)
        else:
            if is_compiled_re:
                raise ValueError(
                    "Cannot use a compiled regex as replacement pattern with "
                    "regex=False"
                )
            if callable(repl):
                raise ValueError("Cannot use a callable replacement when regex=False")
            f = lambda x: x.replace(pat, repl, n)

        return self._map(f, dtype=str)

    def repeat(self, repeats):
        if is_scalar(repeats):

            def scalar_rep(x):
                try:
                    return bytes.__mul__(x, repeats)
                except TypeError:
                    return str.__mul__(x, repeats)

            return self._map(scalar_rep, dtype=str)
        else:
            from pandas.core.arrays.string_ import StringArray

            def rep(x, r):
                if x is libmissing.NA:
                    return x
                try:
                    return bytes.__mul__(x, r)
                except TypeError:
                    return str.__mul__(x, r)

            repeats = np.asarray(repeats, dtype=object)
            result = libops.vec_binop(np.asarray(self._array), repeats, rep)
            if isinstance(self._array, StringArray):
                # Not going through map, so we have to do this here.
                result = StringArray._from_sequence(result)
            return result

    def match(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = np.nan,
    ):
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        f = lambda x: regex.match(x) is not None
        return self._map(f, na_value=na, dtype=np.dtype(bool))

    def fullmatch(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = np.nan,
    ):
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        f = lambda x: regex.fullmatch(x) is not None
        return self._map(f, na_value=na, dtype=np.dtype(bool))

    def encode(self, encoding, errors="strict"):
        f = lambda x: x.encode(encoding, errors=errors)
        return self._map(f, dtype=object)

    def find(self, sub, start=0, end=None):
        return self._find(sub, start, end, side="left")

    def rfind(self, sub, start=0, end=None):
        return self._find(sub, start, end, side="right")

    def _find(self, sub, start, end, side):
        if side == "left":
            method = "find"
        elif side == "right":
            method = "rfind"
        else:  # pragma: no cover
            raise ValueError("Invalid side")

        if end is None:
            f = lambda x: getattr(x, method)(sub, start)
        else:
            f = lambda x: getattr(x, method)(sub, start, end)
        return self._map(f, dtype="int64")

    def findall(self, pat, flags=0):
        regex = re.compile(pat, flags=flags)
        return self._map(regex.findall, dtype="object")

    def get(self, i):
        def f(x):
            if isinstance(x, dict):
                return x.get(i)
            elif len(x) > i >= -len(x):
                return x[i]
            return np.nan

        return self._map(f)

    def index(self, sub, start=0, end=None):
        if end:
            f = lambda x: x.index(sub, start, end)
        else:
            f = lambda x: x.index(sub, start, end)
        return self._map(f, dtype="int64")

    def rindex(self, sub, start=0, end=None):
        if end:
            f = lambda x: x.rindex(sub, start, end)
        else:
            f = lambda x: x.rindex(sub, start, end)
        return self._map(f, dtype="int64")

    def join(self, sep):
        return self._map(sep.join)

    def partition(self, sep, expand):
        result = self._map(lambda x: x.partition(sep), dtype="object")
        return result

    def rpartition(self, sep, expand):
        return self._map(lambda x: x.rpartition(sep), dtype="object")

    def len(self):
        return self._map(len, dtype="int64")

    def slice(self, start=None, stop=None, step=None):
        obj = slice(start, stop, step)
        return self._map(lambda x: x[obj])

    def slice_replace(self, start=None, stop=None, repl=None):
        if repl is None:
            repl = ""

        def f(x):
            if x[start:stop] == "":
                local_stop = start
            else:
                local_stop = stop
            y = ""
            if start is not None:
                y += x[:start]
            y += repl
            if stop is not None:
                y += x[local_stop:]
            return y

        return self._map(f)

    def split(self, pat=None, n=-1, expand=False):
        if pat is None:
            if n is None or n == 0:
                n = -1
            f = lambda x: x.split(pat, n)
        else:
            if len(pat) == 1:
                if n is None or n == 0:
                    n = -1
                f = lambda x: x.split(pat, n)
            else:
                if n is None or n == -1:
                    n = 0
                regex = re.compile(pat)
                f = lambda x: regex.split(x, maxsplit=n)
        return self._map(f, dtype=object)

    def rsplit(self, pat=None, n=-1):
        if n is None or n == 0:
            n = -1
        f = lambda x: x.rsplit(pat, n)
        return self._map(f, dtype="object")

    def translate(self, table):
        return self._map(lambda x: x.translate(table))

    def wrap(self, width, **kwargs):
        kwargs["width"] = width
        tw = textwrap.TextWrapper(**kwargs)
        return self._map(lambda s: "\n".join(tw.wrap(s)))

    def get_dummies(self, sep="|"):
        from pandas import Series

        arr = Series(self._array).fillna("")
        try:
            arr = sep + arr + sep
        except TypeError:
            arr = sep + arr.astype(str) + sep

        tags = set()
        for ts in Series(arr).str.split(sep):
            tags.update(ts)
        tags = sorted(tags - {""})

        dummies = np.empty((len(arr), len(tags)), dtype=np.int64)

        for i, t in enumerate(tags):
            pat = sep + t + sep
            dummies[:, i] = lib.map_infer(arr.to_numpy(), lambda x: pat in x)
        return dummies, tags

    def upper(self):
        return self._map(lambda x: x.upper())

    def isalnum(self):
        return self._map(str.isalnum, dtype="bool")

    def isalpha(self):
        return self._map(str.isalpha, dtype="bool")

    def isdecimal(self):
        return self._map(str.isdecimal, dtype="bool")

    def isdigit(self):
        return self._map(str.isdigit, dtype="bool")

    def islower(self):
        return self._map(str.islower, dtype="bool")

    def isnumeric(self):
        return self._map(str.isnumeric, dtype="bool")

    def isspace(self):
        return self._map(str.isspace, dtype="bool")

    def istitle(self):
        return self._map(str.istitle, dtype="bool")

    def isupper(self):
        return self._map(str.isupper, dtype="bool")

    def capitalize(self):
        return self._map(str.capitalize)

    def casefold(self):
        return self._map(str.casefold)

    def title(self):
        return self._map(str.title)

    def swapcase(self):
        return self._map(str.swapcase)

    def lower(self):
        return self._map(str.lower)

    def normalize(self, form):
        f = lambda x: unicodedata.normalize(form, x)
        return self._map(f)

    def strip(self, to_strip=None):
        return self._map(lambda x: x.strip(to_strip))

    def lstrip(self, to_strip=None):
        return self._map(lambda x: x.lstrip(to_strip))

    def rstrip(self, to_strip=None):
        return self._map(lambda x: x.rstrip(to_strip))


class ObjectProxy(PandasArray):
    _str = CachedAccessor("str", ObjectArrayMethods)
