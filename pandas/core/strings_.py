import codecs
from functools import wraps
import re
import textwrap
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Pattern, Type, Union
import unicodedata
import warnings

import numpy as np

import pandas._libs.lib as lib
import pandas._libs.missing as libmissing
import pandas._libs.ops as libops
from pandas._typing import ArrayLike, Dtype, Scalar
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    ensure_object,
    is_bool_dtype,
    is_categorical_dtype,
    is_extension_array_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_re,
    is_scalar,
    is_string_dtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndexClass,
    ABCMultiIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

from pandas.core.accessor import CachedAccessor
from pandas.core.algorithms import take_1d
from pandas.core.arrays.numpy_ import PandasArray

if TYPE_CHECKING:
    from pandas.arrays import StringArray

_cpython_optimized_encoders = (
    "utf-8",
    "utf8",
    "latin-1",
    "latin1",
    "iso-8859-1",
    "mbcs",
    "ascii",
)
_cpython_optimized_decoders = _cpython_optimized_encoders + ("utf-16", "utf-32")

_shared_docs: Dict[str, str] = dict()


def cat_core(list_of_columns: List, sep: str):
    """
    Auxiliary function for :meth:`str.cat`

    Parameters
    ----------
    list_of_columns : list of numpy arrays
        List of arrays to be concatenated with sep;
        these arrays may not contain NaNs!
    sep : string
        The separator string for concatenating the columns.

    Returns
    -------
    nd.array
        The concatenation of list_of_columns with sep.
    """
    if sep == "":
        # no need to interleave sep if it is empty
        arr_of_cols = np.asarray(list_of_columns, dtype=object)
        return np.sum(arr_of_cols, axis=0)
    list_with_sep = [sep] * (2 * len(list_of_columns) - 1)
    list_with_sep[::2] = list_of_columns
    arr_with_sep = np.asarray(list_with_sep, dtype=object)
    return np.sum(arr_with_sep, axis=0)


def cat_safe(list_of_columns: List, sep: str):
    """
    Auxiliary function for :meth:`str.cat`.

    Same signature as cat_core, but handles TypeErrors in concatenation, which
    happen if the arrays in list_of columns have the wrong dtypes or content.

    Parameters
    ----------
    list_of_columns : list of numpy arrays
        List of arrays to be concatenated with sep;
        these arrays may not contain NaNs!
    sep : string
        The separator string for concatenating the columns.

    Returns
    -------
    nd.array
        The concatenation of list_of_columns with sep.
    """
    try:
        result = cat_core(list_of_columns, sep)
    except TypeError:
        # if there are any non-string values (wrong dtype or hidden behind
        # object dtype), np.sum will fail; catch and return with better message
        for column in list_of_columns:
            dtype = lib.infer_dtype(column, skipna=True)
            if dtype not in ["string", "empty"]:
                raise TypeError(
                    "Concatenation requires list-likes containing only "
                    "strings (or missing values). Offending values found in "
                    f"column {dtype}"
                ) from None
    return result


def str_count(arr, pat, flags=0):
    regex = re.compile(pat, flags=flags)
    f = lambda x: len(regex.findall(x))
    return _na_map(f, arr, dtype="int64")


def str_contains(arr, pat, case=True, flags=0, na=np.nan, regex=True):
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
            f = lambda x: upper_pat in x
            uppered = _na_map(lambda x: x.upper(), arr)
            return _na_map(f, uppered, na, dtype=np.dtype(bool))
    return _na_map(f, arr, na, dtype=np.dtype(bool))


def str_startswith(arr, pat, na=np.nan):
    f = lambda x: x.startswith(pat)
    return _na_map(f, arr, na, dtype=np.dtype(bool))


def str_endswith(arr, pat, na=np.nan):
   f = lambda x: x.endswith(pat)
    return _na_map(f, arr, na, dtype=np.dtype(bool))


def str_replace(arr, pat, repl, n=-1, case=None, flags=0, regex=True):
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
                "Cannot use a compiled regex as replacement pattern with regex=False"
            )
        if callable(repl):
            raise ValueError("Cannot use a callable replacement when regex=False")
        f = lambda x: x.replace(pat, repl, n)

    return _na_map(f, arr, dtype=str)


def str_repeat(arr, repeats):
    if is_scalar(repeats):

        def scalar_rep(x):
            try:
                return bytes.__mul__(x, repeats)
            except TypeError:
                return str.__mul__(x, repeats)

        return _na_map(scalar_rep, arr, dtype=str)
    else:

        def rep(x, r):
            if x is libmissing.NA:
                return x
            try:
                return bytes.__mul__(x, r)
            except TypeError:
                return str.__mul__(x, r)

        repeats = np.asarray(repeats, dtype=object)
        result = libops.vec_binop(np.asarray(arr), repeats, rep)
        return result


def str_match(
    arr: ArrayLike,
    pat: Union[str, Pattern],
    case: bool = True,
    flags: int = 0,
    na: Scalar = np.nan,
):
    if not case:
        flags |= re.IGNORECASE

    regex = re.compile(pat, flags=flags)

    f = lambda x: regex.match(x) is not None

    return _na_map(f, arr, na, dtype=np.dtype(bool))


def str_fullmatch(
    arr: ArrayLike,
    pat: Union[str, Pattern],
    case: bool = True,
    flags: int = 0,
    na: Scalar = np.nan,
):
   if not case:
        flags |= re.IGNORECASE

    regex = re.compile(pat, flags=flags)

    f = lambda x: regex.fullmatch(x) is not None

    return _na_map(f, arr, na, dtype=np.dtype(bool))


def _get_single_group_name(rx):
    try:
        return list(rx.groupindex.keys()).pop()
    except IndexError:
        return None


def _groups_or_na_fun(regex):
    """Used in both extract_noexpand and extract_frame"""
    if regex.groups == 0:
        raise ValueError("pattern contains no capture groups")
    empty_row = [np.nan] * regex.groups

    def f(x):
        if not isinstance(x, str):
            return empty_row
        m = regex.search(x)
        if m:
            return [np.nan if item is None else item for item in m.groups()]
        else:
            return empty_row

    return f


def _result_dtype(arr):
    # workaround #27953
    # ideally we just pass `dtype=arr.dtype` unconditionally, but this fails
    # when the list of values is empty.
    from pandas.core.arrays.string_ import StringDtype

    if isinstance(arr.dtype, StringDtype):
        return arr.dtype.name
    else:
        return object


def _str_extract_noexpand(arr, pat, flags=0):
    """
    Find groups in each string in the Series using passed regular
    expression. This function is called from
    str_extract(expand=False), and can return Series, DataFrame, or
    Index.

    """
    from pandas import DataFrame

    regex = re.compile(pat, flags=flags)
    groups_or_na = _groups_or_na_fun(regex)

    if regex.groups == 1:
        result = np.array([groups_or_na(val)[0] for val in arr], dtype=object)
        name = _get_single_group_name(regex)
    else:
        if isinstance(arr, ABCIndexClass):
            raise ValueError("only one regex group is supported with Index")
        name = None
        names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
        columns = [names.get(1 + i, i) for i in range(regex.groups)]
        if arr.empty:
            result = DataFrame(columns=columns, dtype=object)
        else:
            dtype = _result_dtype(arr)
            result = DataFrame(
                [groups_or_na(val) for val in arr],
                columns=columns,
                index=arr.index,
                dtype=dtype,
            )
    return result, name


def _str_extract_frame(arr, pat, flags=0):
    """
    For each subject string in the Series, extract groups from the
    first match of regular expression pat. This function is called from
    str_extract(expand=True), and always returns a DataFrame.

    """
    from pandas import DataFrame

    regex = re.compile(pat, flags=flags)
    groups_or_na = _groups_or_na_fun(regex)
    names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
    columns = [names.get(1 + i, i) for i in range(regex.groups)]

    if len(arr) == 0:
        return DataFrame(columns=columns, dtype=object)
    try:
        result_index = arr.index
    except AttributeError:
        result_index = None
    dtype = _result_dtype(arr)
    return DataFrame(
        [groups_or_na(val) for val in arr],
        columns=columns,
        index=result_index,
        dtype=dtype,
    )


def str_extract(arr, pat, flags=0, expand=True):
   if not isinstance(expand, bool):
        raise ValueError("expand must be True or False")
    if expand:
        return _str_extract_frame(arr._orig, pat, flags=flags)
    else:
        result, name = _str_extract_noexpand(arr._parent, pat, flags=flags)
        return arr._wrap_result(result, name=name, expand=expand)


def str_extractall(arr, pat, flags=0):
   regex = re.compile(pat, flags=flags)
    # the regex must contain capture groups.
    if regex.groups == 0:
        raise ValueError("pattern contains no capture groups")

    if isinstance(arr, ABCIndexClass):
        arr = arr.to_series().reset_index(drop=True)

    names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
    columns = [names.get(1 + i, i) for i in range(regex.groups)]
    match_list = []
    index_list = []
    is_mi = arr.index.nlevels > 1

    for subject_key, subject in arr.items():
        if isinstance(subject, str):

            if not is_mi:
                subject_key = (subject_key,)

            for match_i, match_tuple in enumerate(regex.findall(subject)):
                if isinstance(match_tuple, str):
                    match_tuple = (match_tuple,)
                na_tuple = [np.NaN if group == "" else group for group in match_tuple]
                match_list.append(na_tuple)
                result_key = tuple(subject_key + (match_i,))
                index_list.append(result_key)

    from pandas import MultiIndex

    index = MultiIndex.from_tuples(index_list, names=arr.index.names + ["match"])
    dtype = _result_dtype(arr)

    result = arr._constructor_expanddim(
        match_list, index=index, columns=columns, dtype=dtype
    )
    return result


def str_get_dummies(arr, sep="|"):
    arr = arr.fillna("")
    try:
        arr = sep + arr + sep
    except TypeError:
        arr = sep + arr.astype(str) + sep

    tags = set()
    for ts in arr.str.split(sep):
        tags.update(ts)
    tags = sorted(tags - {""})

    dummies = np.empty((len(arr), len(tags)), dtype=np.int64)

    for i, t in enumerate(tags):
        pat = sep + t + sep
        dummies[:, i] = lib.map_infer(arr.to_numpy(), lambda x: pat in x)
    return dummies, tags


def str_join(arr, sep):
    return _na_map(sep.join, arr, dtype=str)


def str_findall(arr, pat, flags=0):
    regex = re.compile(pat, flags=flags)
    return _na_map(regex.findall, arr)


def str_find(arr, sub, start=0, end=None, side="left"):
   if not isinstance(sub, str):
        msg = f"expected a string object, not {type(sub).__name__}"
        raise TypeError(msg)

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

    return _na_map(f, arr, dtype=np.dtype("int64"))


def str_index(arr, sub, start=0, end=None, side="left"):
    if not isinstance(sub, str):
        msg = f"expected a string object, not {type(sub).__name__}"
        raise TypeError(msg)

    if side == "left":
        method = "index"
    elif side == "right":
        method = "rindex"
    else:  # pragma: no cover
        raise ValueError("Invalid side")

    if end is None:
        f = lambda x: getattr(x, method)(sub, start)
    else:
        f = lambda x: getattr(x, method)(sub, start, end)

    return _na_map(f, arr, dtype=np.dtype("int64"))


def str_pad(arr, width, side="left", fillchar=" "):
    if not isinstance(fillchar, str):
        msg = f"fillchar must be a character, not {type(fillchar).__name__}"
        raise TypeError(msg)

    if len(fillchar) != 1:
        raise TypeError("fillchar must be a character, not str")

    if not is_integer(width):
        msg = f"width must be of integer type, not {type(width).__name__}"
        raise TypeError(msg)

    if side == "left":
        f = lambda x: x.rjust(width, fillchar)
    elif side == "right":
        f = lambda x: x.ljust(width, fillchar)
    elif side == "both":
        f = lambda x: x.center(width, fillchar)
    else:  # pragma: no cover
        raise ValueError("Invalid side")

    return _na_map(f, arr, dtype=str)


def str_split(arr, pat=None, n=None):

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
    res = _na_map(f, arr)
    return res


def str_rsplit(arr, pat=None, n=None):

    if n is None or n == 0:
        n = -1
    f = lambda x: x.rsplit(pat, n)
    res = _na_map(f, arr)
    return res


def str_slice(arr, start=None, stop=None, step=None):
    obj = slice(start, stop, step)
    f = lambda x: x[obj]
    return _na_map(f, arr, dtype=str)


def str_slice_replace(arr, start=None, stop=None, repl=None):
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

    return _na_map(f, arr, dtype=str)


def str_strip(arr, to_strip=None, side="both"):
    """
    Strip whitespace (including newlines) from each string in the
    Series/Index.

    Parameters
    ----------
    to_strip : str or unicode
    side : {'left', 'right', 'both'}, default 'both'

    Returns
    -------
    Series or Index
    """
    if side == "both":
        f = lambda x: x.strip(to_strip)
    elif side == "left":
        f = lambda x: x.lstrip(to_strip)
    elif side == "right":
        f = lambda x: x.rstrip(to_strip)
    else:  # pragma: no cover
        raise ValueError("Invalid side")
    return _na_map(f, arr, dtype=str)


def str_wrap(arr, width, **kwargs):
    kwargs["width"] = width

    tw = textwrap.TextWrapper(**kwargs)

    return _na_map(lambda s: "\n".join(tw.wrap(s)), arr, dtype=str)


def str_translate(arr, table):
     return _na_map(lambda x: x.translate(table), arr, dtype=str)


def str_get(arr, i):
    def f(x):
        if isinstance(x, dict):
            return x.get(i)
        elif len(x) > i >= -len(x):
            return x[i]
        return np.nan

    return _na_map(f, arr)


def str_decode(arr, encoding, errors="strict"):
    """
    Decode character string in the Series/Index using indicated encoding.

    Equivalent to :meth:`str.decode` in python2 and :meth:`bytes.decode` in
    python3.

    Parameters
    ----------
    encoding : str
    errors : str, optional

    Returns
    -------
    Series or Index
    """
    if encoding in _cpython_optimized_decoders:
        # CPython optimized implementation
        f = lambda x: x.decode(encoding, errors)
    else:
        decoder = codecs.getdecoder(encoding)
        f = lambda x: decoder(x, errors)[0]
    return _na_map(f, arr)


def str_encode(arr, encoding, errors="strict"):
    if encoding in _cpython_optimized_encoders:
        # CPython optimized implementation
        f = lambda x: x.encode(encoding, errors)
    else:
        encoder = codecs.getencoder(encoding)
        f = lambda x: encoder(x, errors)[0]
    return _na_map(f, arr)


def forbid_nonstring_types(forbidden, name=None):
    """
    Decorator to forbid specific types for a method of StringMethods.

    For calling `.str.{method}` on a Series or Index, it is necessary to first
    initialize the :class:`StringMethods` object, and then call the method.
    However, different methods allow different input types, and so this can not
    be checked during :meth:`StringMethods.__init__`, but must be done on a
    per-method basis. This decorator exists to facilitate this process, and
    make it explicit which (inferred) types are disallowed by the method.

    :meth:`StringMethods.__init__` allows the *union* of types its different
    methods allow (after skipping NaNs; see :meth:`StringMethods._validate`),
    namely: ['string', 'empty', 'bytes', 'mixed', 'mixed-integer'].

    The default string types ['string', 'empty'] are allowed for all methods.
    For the additional types ['bytes', 'mixed', 'mixed-integer'], each method
    then needs to forbid the types it is not intended for.

    Parameters
    ----------
    forbidden : list-of-str or None
        List of forbidden non-string types, may be one or more of
        `['bytes', 'mixed', 'mixed-integer']`.
    name : str, default None
        Name of the method to use in the error message. By default, this is
        None, in which case the name from the method being wrapped will be
        copied. However, for working with further wrappers (like _pat_wrapper
        and _noarg_wrapper), it is necessary to specify the name.

    Returns
    -------
    func : wrapper
        The method to which the decorator is applied, with an added check that
        enforces the inferred type to not be in the list of forbidden types.

    Raises
    ------
    TypeError
        If the inferred type of the underlying data is in `forbidden`.
    """
    # deal with None
    forbidden = [] if forbidden is None else forbidden

    allowed_types = {"string", "empty", "bytes", "mixed", "mixed-integer"} - set(
        forbidden
    )

    def _forbid_nonstring_types(func):
        func_name = func.__name__ if name is None else name

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            if self._inferred_dtype not in allowed_types:
                msg = (
                    f"Cannot use .str.{func_name} with values of "
                    f"inferred dtype '{self._inferred_dtype}'."
                )
                raise TypeError(msg)
            return func(self, *args, **kwargs)

        wrapper.__name__ = func_name
        return wrapper

    return _forbid_nonstring_types


def _noarg_wrapper(
    f,
    name=None,
    docstring=None,
    forbidden_types=["bytes"],
    returns_string=True,
    **kwargs,
):
    @forbid_nonstring_types(forbidden_types, name=name)
    def wrapper(self):
        result = _na_map(f, self._parent, **kwargs)
        return self._wrap_result(result, returns_string=returns_string)

    wrapper.__name__ = f.__name__ if name is None else name
    if docstring is not None:
        wrapper.__doc__ = docstring
    else:
        raise ValueError("Provide docstring")

    return wrapper


def _pat_wrapper(
    f,
    flags=False,
    na=False,
    name=None,
    forbidden_types=["bytes"],
    returns_string=True,
    **kwargs,
):
    @forbid_nonstring_types(forbidden_types, name=name)
    def wrapper1(self, pat):
        result = f(self._parent, pat)
        return self._wrap_result(result, returns_string=returns_string)

    @forbid_nonstring_types(forbidden_types, name=name)
    def wrapper2(self, pat, flags=0, **kwargs):
        result = f(self._parent, pat, flags=flags, **kwargs)
        return self._wrap_result(result, returns_string=returns_string)

    @forbid_nonstring_types(forbidden_types, name=name)
    def wrapper3(self, pat, na=np.nan):
        result = f(self._parent, pat, na=na)
        return self._wrap_result(result, returns_string=returns_string)

    wrapper = wrapper3 if na else wrapper2 if flags else wrapper1

    wrapper.__name__ = f.__name__ if name is None else name
    if f.__doc__:
        wrapper.__doc__ = f.__doc__

    return wrapper


def copy(source):
    """Copy a docstring from another source function (if present)"""

    def do_copy(target):
        if source.__doc__:
            target.__doc__ = source.__doc__
        return target

    return do_copy
