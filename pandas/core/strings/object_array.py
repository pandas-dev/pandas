import re
import textwrap
from typing import Pattern, Set, Union, cast
import unicodedata
import warnings

import numpy as np

import pandas._libs.lib as lib
import pandas._libs.missing as libmissing
import pandas._libs.ops as libops
from pandas._typing import Scalar

from pandas.core.dtypes.common import is_re, is_scalar
from pandas.core.dtypes.missing import isna

from pandas.core.strings.base import BaseStringArrayMethods


class ObjectStringArrayMixin(BaseStringArrayMethods):
    """
    String Methods operating on object-dtype ndarrays.
    """

    _str_na_value = np.nan

    def __len__(self):
        # For typing, _str_map relies on the object being sized.
        raise NotImplementedError

    def _str_map(self, f, na_value=None, dtype=None):
        """
        Map a callable over valid element of the array.

        Parameters
        ----------
        f : Callable
            A function to call on each non-NA element.
        na_value : Scalar, optional
            The value to set for NA values. Might also be used for the
            fill value if the callable `f` raises an exception.
            This defaults to ``self._str_na_value`` which is ``np.nan``
            for object-dtype and Categorical and ``pd.NA`` for StringArray.
        dtype : Dtype, optional
            The dtype of the result array.
        """
        arr = self
        if dtype is None:
            dtype = np.dtype("object")
        if na_value is None:
            na_value = self._str_na_value

        if not len(arr):
            return np.ndarray(0, dtype=dtype)

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

            return self._str_map(g, na_value=na_value, dtype=dtype)
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result

    def _str_count(self, pat, flags=0):
        regex = re.compile(pat, flags=flags)
        f = lambda x: len(regex.findall(x))
        return self._str_map(f, dtype="int64")

    def _str_pad(self, width, side="left", fillchar=" "):
        if side == "left":
            f = lambda x: x.rjust(width, fillchar)
        elif side == "right":
            f = lambda x: x.ljust(width, fillchar)
        elif side == "both":
            f = lambda x: x.center(width, fillchar)
        else:  # pragma: no cover
            raise ValueError("Invalid side")
        return self._str_map(f)

    def _str_contains(self, pat, case=True, flags=0, na=np.nan, regex=True):
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
        return self._str_map(f, na, dtype=np.dtype("bool"))

    def _str_startswith(self, pat, na=None):
        f = lambda x: x.startswith(pat)
        return self._str_map(f, na_value=na, dtype=np.dtype(bool))

    def _str_endswith(self, pat, na=None):
        f = lambda x: x.endswith(pat)
        return self._str_map(f, na_value=na, dtype=np.dtype(bool))

    def _str_replace(self, pat, repl, n=-1, case=None, flags=0, regex=True):
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

        return self._str_map(f, dtype=str)

    def _str_repeat(self, repeats):
        if is_scalar(repeats):

            def scalar_rep(x):
                try:
                    return bytes.__mul__(x, repeats)
                except TypeError:
                    return str.__mul__(x, repeats)

            return self._str_map(scalar_rep, dtype=str)
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
            result = libops.vec_binop(np.asarray(self), repeats, rep)
            if isinstance(self, StringArray):
                # Not going through map, so we have to do this here.
                result = StringArray._from_sequence(result)
            return result

    def _str_match(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = None,
    ):
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        f = lambda x: regex.match(x) is not None
        return self._str_map(f, na_value=na, dtype=np.dtype(bool))

    def _str_fullmatch(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = None,
    ):
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        f = lambda x: regex.fullmatch(x) is not None
        return self._str_map(f, na_value=na, dtype=np.dtype(bool))

    def _str_encode(self, encoding, errors="strict"):
        f = lambda x: x.encode(encoding, errors=errors)
        return self._str_map(f, dtype=object)

    def _str_find(self, sub, start=0, end=None):
        return self._str_find_(sub, start, end, side="left")

    def _str_rfind(self, sub, start=0, end=None):
        return self._str_find_(sub, start, end, side="right")

    def _str_find_(self, sub, start, end, side):
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
        return self._str_map(f, dtype="int64")

    def _str_findall(self, pat, flags=0):
        regex = re.compile(pat, flags=flags)
        return self._str_map(regex.findall, dtype="object")

    def _str_get(self, i):
        def f(x):
            if isinstance(x, dict):
                return x.get(i)
            elif len(x) > i >= -len(x):
                return x[i]
            return self._str_na_value

        return self._str_map(f)

    def _str_index(self, sub, start=0, end=None):
        if end:
            f = lambda x: x.index(sub, start, end)
        else:
            f = lambda x: x.index(sub, start, end)
        return self._str_map(f, dtype="int64")

    def _str_rindex(self, sub, start=0, end=None):
        if end:
            f = lambda x: x.rindex(sub, start, end)
        else:
            f = lambda x: x.rindex(sub, start, end)
        return self._str_map(f, dtype="int64")

    def _str_join(self, sep):
        return self._str_map(sep.join)

    def _str_partition(self, sep, expand):
        result = self._str_map(lambda x: x.partition(sep), dtype="object")
        return result

    def _str_rpartition(self, sep, expand):
        return self._str_map(lambda x: x.rpartition(sep), dtype="object")

    def _str_len(self):
        return self._str_map(len, dtype="int64")

    def _str_slice(self, start=None, stop=None, step=None):
        obj = slice(start, stop, step)
        return self._str_map(lambda x: x[obj])

    def _str_slice_replace(self, start=None, stop=None, repl=None):
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

        return self._str_map(f)

    def _str_split(self, pat=None, n=-1, expand=False):
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
        return self._str_map(f, dtype=object)

    def _str_rsplit(self, pat=None, n=-1):
        if n is None or n == 0:
            n = -1
        f = lambda x: x.rsplit(pat, n)
        return self._str_map(f, dtype="object")

    def _str_translate(self, table):
        return self._str_map(lambda x: x.translate(table))

    def _str_wrap(self, width, **kwargs):
        kwargs["width"] = width
        tw = textwrap.TextWrapper(**kwargs)
        return self._str_map(lambda s: "\n".join(tw.wrap(s)))

    def _str_get_dummies(self, sep="|"):
        from pandas import Series

        arr = Series(self).fillna("")
        try:
            arr = sep + arr + sep
        except TypeError:
            arr = cast(Series, arr)
            arr = sep + arr.astype(str) + sep
        arr = cast(Series, arr)

        tags: Set[str] = set()
        for ts in Series(arr).str.split(sep):
            tags.update(ts)
        tags2 = sorted(tags - {""})

        dummies = np.empty((len(arr), len(tags2)), dtype=np.int64)

        for i, t in enumerate(tags2):
            pat = sep + t + sep
            dummies[:, i] = lib.map_infer(arr.to_numpy(), lambda x: pat in x)
        return dummies, tags2

    def _str_upper(self):
        return self._str_map(lambda x: x.upper())

    def _str_isalnum(self):
        return self._str_map(str.isalnum, dtype="bool")

    def _str_isalpha(self):
        return self._str_map(str.isalpha, dtype="bool")

    def _str_isdecimal(self):
        return self._str_map(str.isdecimal, dtype="bool")

    def _str_isdigit(self):
        return self._str_map(str.isdigit, dtype="bool")

    def _str_islower(self):
        return self._str_map(str.islower, dtype="bool")

    def _str_isnumeric(self):
        return self._str_map(str.isnumeric, dtype="bool")

    def _str_isspace(self):
        return self._str_map(str.isspace, dtype="bool")

    def _str_istitle(self):
        return self._str_map(str.istitle, dtype="bool")

    def _str_isupper(self):
        return self._str_map(str.isupper, dtype="bool")

    def _str_capitalize(self):
        return self._str_map(str.capitalize)

    def _str_casefold(self):
        return self._str_map(str.casefold)

    def _str_title(self):
        return self._str_map(str.title)

    def _str_swapcase(self):
        return self._str_map(str.swapcase)

    def _str_lower(self):
        return self._str_map(str.lower)

    def _str_normalize(self, form):
        f = lambda x: unicodedata.normalize(form, x)
        return self._str_map(f)

    def _str_strip(self, to_strip=None):
        return self._str_map(lambda x: x.strip(to_strip))

    def _str_lstrip(self, to_strip=None):
        return self._str_map(lambda x: x.lstrip(to_strip))

    def _str_rstrip(self, to_strip=None):
        return self._str_map(lambda x: x.rstrip(to_strip))
