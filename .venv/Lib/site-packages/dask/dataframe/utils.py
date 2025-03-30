from __future__ import annotations

import math
import re
import sys
import textwrap
import traceback
from collections.abc import Callable, Iterable, Mapping, Sequence
from contextlib import contextmanager
from numbers import Number
from typing import TypeVar, overload

import numpy as np
import pandas as pd
from pandas.api.types import is_dtype_equal

import dask
from dask.base import is_dask_collection
from dask.core import get_deps
from dask.dataframe._compat import tm  # noqa: F401
from dask.dataframe.dispatch import (  # noqa : F401
    is_categorical_dtype_dispatch,
    make_meta,
    make_meta_obj,
    meta_nonempty,
)
from dask.dataframe.extensions import make_scalar
from dask.typing import NoDefault, no_default
from dask.utils import (
    asciitable,
    is_dataframe_like,
    is_index_like,
    is_series_like,
    typename,
)

meta_object_types: tuple[type, ...] = (pd.Series, pd.DataFrame, pd.Index, pd.MultiIndex)
try:
    import scipy.sparse as sp

    meta_object_types += (sp.spmatrix,)
except ImportError:
    pass


def is_scalar(x):
    # np.isscalar does not work for some pandas scalars, for example pd.NA
    if isinstance(x, (Sequence, Iterable)) and not isinstance(x, str):
        return False
    elif hasattr(x, "dtype"):
        return isinstance(x, np.ScalarType)
    if isinstance(x, dict):
        return False
    if isinstance(x, (str, int)) or x is None:
        return True

    from dask.dataframe.dask_expr._expr import Expr

    return not isinstance(x, Expr)


def is_integer_na_dtype(t):
    dtype = getattr(t, "dtype", t)
    types = (
        pd.Int8Dtype,
        pd.Int16Dtype,
        pd.Int32Dtype,
        pd.Int64Dtype,
        pd.UInt8Dtype,
        pd.UInt16Dtype,
        pd.UInt32Dtype,
        pd.UInt64Dtype,
    )
    return isinstance(dtype, types)


def is_float_na_dtype(t):
    dtype = getattr(t, "dtype", t)
    types = (
        pd.Float32Dtype,
        pd.Float64Dtype,
    )
    return isinstance(dtype, types)


_META_TYPES = "meta : pd.DataFrame, pd.Series, dict, iterable, tuple, optional"
_META_DESCRIPTION = """\
An empty ``pd.DataFrame`` or ``pd.Series`` that matches the dtypes and
column names of the output. This metadata is necessary for many algorithms
in dask dataframe to work.  For ease of use, some alternative inputs are
also available. Instead of a ``DataFrame``, a ``dict`` of ``{name: dtype}``
or iterable of ``(name, dtype)`` can be provided (note that the order of
the names should match the order of the columns). Instead of a series, a
tuple of ``(name, dtype)`` can be used. If not provided, dask will try to
infer the metadata. This may lead to unexpected results, so providing
``meta`` is recommended. For more information, see
``dask.dataframe.utils.make_meta``.
"""

T = TypeVar("T", bound=Callable)


@overload
def insert_meta_param_description(func: T) -> T: ...


@overload
def insert_meta_param_description(pad: int) -> Callable[[T], T]: ...


def insert_meta_param_description(*args, **kwargs):
    """Replace `$META` in docstring with param description.

    If pad keyword is provided, will pad description by that number of
    spaces (default is 8)."""
    if not args:
        return lambda f: insert_meta_param_description(f, **kwargs)
    f = args[0]
    indent = " " * kwargs.get("pad", 8)
    body = textwrap.wrap(
        _META_DESCRIPTION, initial_indent=indent, subsequent_indent=indent, width=78
    )
    descr = "{}\n{}".format(_META_TYPES, "\n".join(body))
    if f.__doc__:
        if "$META" in f.__doc__:
            f.__doc__ = f.__doc__.replace("$META", descr)
        else:
            # Put it at the end of the parameters section
            parameter_header = "Parameters\n%s----------" % indent[4:]
            first, last = re.split("Parameters\\n[ ]*----------", f.__doc__)
            parameters, rest = last.split("\n\n", 1)
            f.__doc__ = "{}{}{}\n{}{}\n\n{}".format(
                first, parameter_header, parameters, indent[4:], descr, rest
            )
    return f


@contextmanager
def raise_on_meta_error(funcname=None, udf=False):
    """Reraise errors in this block to show metadata inference failure.

    Parameters
    ----------
    funcname : str, optional
        If provided, will be added to the error message to indicate the
        name of the method that failed.
    """
    try:
        yield
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = "".join(traceback.format_tb(exc_traceback))
        msg = "Metadata inference failed{0}.\n\n"
        if udf:
            msg += (
                "You have supplied a custom function and Dask is unable to \n"
                "determine the type of output that that function returns. \n\n"
                "To resolve this please provide a meta= keyword.\n"
                "The docstring of the Dask function you ran should have more information.\n\n"
            )
        msg += (
            "Original error is below:\n"
            "------------------------\n"
            "{1}\n\n"
            "Traceback:\n"
            "---------\n"
            "{2}"
        )
        msg = msg.format(f" in `{funcname}`" if funcname else "", repr(e), tb)
        raise ValueError(msg) from e


UNKNOWN_CATEGORIES = "__UNKNOWN_CATEGORIES__"


def has_known_categories(x):
    """Returns whether the categories in `x` are known.

    Parameters
    ----------
    x : Series or CategoricalIndex
    """
    x = getattr(x, "_meta", x)
    if is_series_like(x):
        return UNKNOWN_CATEGORIES not in x.cat.categories
    elif is_index_like(x) and hasattr(x, "categories"):
        return UNKNOWN_CATEGORIES not in x.categories
    raise TypeError("Expected Series or CategoricalIndex")


def strip_unknown_categories(x, just_drop_unknown=False):
    """Replace any unknown categoricals with empty categoricals.

    Useful for preventing ``UNKNOWN_CATEGORIES`` from leaking into results.
    """
    if isinstance(x, (pd.Series, pd.DataFrame)):
        x = x.copy()
        if isinstance(x, pd.DataFrame):
            cat_mask = x.dtypes == "category"
            if cat_mask.any():
                cats = cat_mask[cat_mask].index
                for c in cats:
                    if not has_known_categories(x[c]):
                        if just_drop_unknown:
                            x[c].cat.remove_categories(UNKNOWN_CATEGORIES, inplace=True)
                        else:
                            x[c] = x[c].cat.set_categories([])
        elif isinstance(x, pd.Series):
            if isinstance(x.dtype, pd.CategoricalDtype) and not has_known_categories(x):
                x = x.cat.set_categories([])
        if isinstance(x.index, pd.CategoricalIndex) and not has_known_categories(
            x.index
        ):
            x.index = x.index.set_categories([])
    elif isinstance(x, pd.CategoricalIndex) and not has_known_categories(x):
        x = x.set_categories([])
    return x


def clear_known_categories(x, cols=None, index=True, dtype_backend=None):
    """Set categories to be unknown.

    Parameters
    ----------
    x : DataFrame, Series, Index
    cols : iterable, optional
        If x is a DataFrame, set only categoricals in these columns to unknown.
        By default, all categorical columns are set to unknown categoricals
    index : bool, optional
        If True and x is a Series or DataFrame, set the clear known categories
        in the index as well.
    dtype_backend : string, optional
        If set to PyArrow, the categorical dtype is implemented as a PyArrow
        dictionary
    """
    if dtype_backend == "pyarrow":
        # Right now Categorical with PyArrow is implemented as dictionary and
        # categorical accessor is not yet available
        return x

    if not is_index_like(x):
        x = x.copy()
        if is_dataframe_like(x):
            mask = x.dtypes == "category"
            if cols is None:
                cols = mask[mask].index
            elif not mask.loc[cols].all():
                raise ValueError("Not all columns are categoricals")
            for c in cols:
                x[c] = x[c].cat.set_categories([UNKNOWN_CATEGORIES])
        elif is_series_like(x):
            if is_categorical_dtype_dispatch(x.dtype):
                x = x.cat.set_categories([UNKNOWN_CATEGORIES])
        if index and is_categorical_dtype_dispatch(x.index.dtype):
            x.index = x.index.set_categories([UNKNOWN_CATEGORIES])
    elif is_categorical_dtype_dispatch(x.dtype):
        x = x.set_categories([UNKNOWN_CATEGORIES])
    return x


def _empty_series(name, dtype, index=None):
    if isinstance(dtype, str) and dtype == "category":
        s = pd.Series(pd.Categorical([UNKNOWN_CATEGORIES]), name=name).iloc[:0]
        if index is not None:
            s.index = make_meta(index)
        return s
    return pd.Series([], dtype=dtype, name=name, index=index)


_simple_fake_mapping = {
    "b": np.bool_(True),
    "V": np.void(b" "),
    "M": np.datetime64("1970-01-01"),
    "m": np.timedelta64(1),
    "S": np.str_("foo"),
    "a": np.str_("foo"),
    "U": np.str_("foo"),
    "O": "foo",
}


def _scalar_from_dtype(dtype):
    if dtype.kind in ("i", "f", "u"):
        return dtype.type(1)
    elif dtype.kind == "c":
        return dtype.type(complex(1, 0))
    elif dtype.kind in _simple_fake_mapping:
        o = _simple_fake_mapping[dtype.kind]
        return o.astype(dtype) if dtype.kind in ("m", "M") else o
    else:
        raise TypeError(f"Can't handle dtype: {dtype}")


def _nonempty_scalar(x):
    if type(x) in make_scalar._lookup:
        return make_scalar(x)

    if np.isscalar(x):
        dtype = x.dtype if hasattr(x, "dtype") else np.dtype(type(x))
        return make_scalar(dtype)

    if x is pd.NA:
        return pd.NA

    raise TypeError(f"Can't handle meta of type '{typename(type(x))}'")


def check_meta(x, meta, funcname=None, numeric_equal=True):
    """Check that the dask metadata matches the result.

    If metadata matches, ``x`` is passed through unchanged. A nice error is
    raised if metadata doesn't match.

    Parameters
    ----------
    x : DataFrame, Series, or Index
    meta : DataFrame, Series, or Index
        The expected metadata that ``x`` should match
    funcname : str, optional
        The name of the function in which the metadata was specified. If
        provided, the function name will be included in the error message to be
        more helpful to users.
    numeric_equal : bool, optionl
        If True, integer and floating dtypes compare equal. This is useful due
        to panda's implicit conversion of integer to floating upon encountering
        missingness, which is hard to infer statically.
    """
    eq_types = {"i", "f", "u"} if numeric_equal else set()

    def equal_dtypes(a, b):
        if isinstance(a, pd.CategoricalDtype) != isinstance(b, pd.CategoricalDtype):
            return False
        if isinstance(a, str) and a == "-" or isinstance(b, str) and b == "-":
            return False
        if isinstance(a, pd.CategoricalDtype) and isinstance(b, pd.CategoricalDtype):
            if UNKNOWN_CATEGORIES in a.categories or UNKNOWN_CATEGORIES in b.categories:
                return True
            return a == b
        return (a.kind in eq_types and b.kind in eq_types) or is_dtype_equal(a, b)

    if not (
        is_dataframe_like(meta) or is_series_like(meta) or is_index_like(meta)
    ) or is_dask_collection(meta):
        raise TypeError(
            "Expected partition to be DataFrame, Series, or "
            "Index, got `%s`" % typename(type(meta))
        )

    # Notice, we use .__class__ as opposed to type() in order to support
    # object proxies see <https://github.com/dask/dask/pull/6981>
    if x.__class__ != meta.__class__:
        errmsg = "Expected partition of type `{}` but got `{}`".format(
            typename(type(meta)),
            typename(type(x)),
        )
    elif is_dataframe_like(meta):
        dtypes = pd.concat([x.dtypes, meta.dtypes], axis=1, sort=True)
        bad_dtypes = [
            (repr(col), a, b)
            for col, a, b in dtypes.fillna("-").itertuples()
            if not equal_dtypes(a, b)
        ]
        if bad_dtypes:
            errmsg = "Partition type: `{}`\n{}".format(
                typename(type(meta)),
                asciitable(["Column", "Found", "Expected"], bad_dtypes),
            )
        else:
            check_matching_columns(meta, x)
            return x
    else:
        if equal_dtypes(x.dtype, meta.dtype):
            return x
        errmsg = "Partition type: `{}`\n{}".format(
            typename(type(meta)),
            asciitable(["", "dtype"], [("Found", x.dtype), ("Expected", meta.dtype)]),
        )

    raise ValueError(
        "Metadata mismatch found%s.\n\n"
        "%s" % ((" in `%s`" % funcname if funcname else ""), errmsg)
    )


def check_matching_columns(meta, actual):
    import dask.dataframe.methods as methods

    # Need nan_to_num otherwise nan comparison gives False
    if not np.array_equal(np.nan_to_num(meta.columns), np.nan_to_num(actual.columns)):
        extra = methods.tolist(actual.columns.difference(meta.columns))
        missing = methods.tolist(meta.columns.difference(actual.columns))
        if extra or missing:
            extra_info = f"  Extra:   {extra}\n  Missing: {missing}"
        else:
            extra_info = (
                f"Order of columns does not match."
                f"\nActual:   {actual.columns.tolist()}"
                f"\nExpected: {meta.columns.tolist()}"
            )
        raise ValueError(
            "The columns in the computed data do not match"
            " the columns in the provided metadata.\n"
            f"{extra_info}"
        )


def index_summary(idx, name=None):
    """Summarized representation of an Index."""
    n = len(idx)
    if name is None:
        name = idx.__class__.__name__
    if n:
        head = idx[0]
        tail = idx[-1]
        summary = f", {head} to {tail}"
    else:
        summary = ""

    return f"{name}: {n} entries{summary}"


###############################################################
# Testing
###############################################################


def _check_dask(dsk, check_names=True, check_dtypes=True, result=None, scheduler=None):
    import dask.dataframe as dd

    if hasattr(dsk, "__dask_graph__"):
        graph = dsk.__dask_graph__()
        if hasattr(graph, "validate"):
            graph.validate()
        if result is None:
            result = dsk.compute(scheduler=scheduler)
        if isinstance(dsk, dd.Index) or is_index_like(dsk._meta):
            assert "Index" in type(result).__name__, type(result)
            # assert type(dsk._meta) == type(result), type(dsk._meta)
            if check_names:
                assert dsk.name == result.name
                assert dsk._meta.name == result.name
                if isinstance(result, pd.MultiIndex):
                    assert result.names == dsk._meta.names
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
        elif isinstance(dsk, dd.Series) or is_series_like(dsk._meta):
            assert "Series" in type(result).__name__, type(result)
            assert type(dsk._meta) == type(result), type(dsk._meta)
            if check_names:
                assert dsk.name == result.name, (dsk.name, result.name)
                assert dsk._meta.name == result.name
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
            _check_dask(
                dsk.index,
                check_names=check_names,
                check_dtypes=check_dtypes,
                result=result.index,
            )
        elif isinstance(dsk, dd.DataFrame) or is_dataframe_like(dsk._meta):
            assert "DataFrame" in type(result).__name__, type(result)
            assert isinstance(dsk.columns, pd.Index), type(dsk.columns)
            assert type(dsk._meta) == type(result), type(dsk._meta)
            if check_names:
                tm.assert_index_equal(dsk.columns, result.columns)
                tm.assert_index_equal(dsk._meta.columns, result.columns)
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
            _check_dask(
                dsk.index,
                check_names=check_names,
                check_dtypes=check_dtypes,
                result=result.index,
            )
        else:
            if not np.isscalar(result) and not isinstance(
                result, (pd.Timestamp, pd.Timedelta)
            ):
                raise TypeError(
                    "Expected object of type dataframe, series, index, or scalar.\n"
                    "    Got: " + str(type(result))
                )
            if check_dtypes:
                assert_dask_dtypes(dsk, result)
        return result
    return dsk


def _maybe_sort(a, check_index: bool):
    import dask.dataframe.methods as methods

    # sort by value, then index
    try:
        if is_dataframe_like(a):
            if set(a.index.names) & set(a.columns):
                a.index.names = [
                    "-overlapped-index-name-%d" % i for i in range(len(a.index.names))
                ]
            a = a.sort_values(by=methods.tolist(a.columns))
        else:
            a = a.sort_values()
    except (TypeError, IndexError, ValueError):
        pass
    return a.sort_index() if check_index else a


def _maybe_convert_string(a, b):
    if pyarrow_strings_enabled():
        from dask.dataframe._pyarrow import to_pyarrow_string

        if isinstance(a, (pd.DataFrame, pd.Series, pd.Index)):
            a = to_pyarrow_string(a)

        if isinstance(b, (pd.DataFrame, pd.Series, pd.Index)):
            b = to_pyarrow_string(b)

    return a, b


def assert_eq_dtypes(a, b):
    a, b = _maybe_convert_string(a, b)
    tm.assert_series_equal(a.dtypes.value_counts(), b.dtypes.value_counts())


def assert_eq(
    a,
    b,
    check_names=True,
    check_dtype=True,
    check_divisions=True,
    check_index=True,
    sort_results=True,
    scheduler="sync",
    **kwargs,
):
    if check_divisions:
        assert_divisions(a, scheduler=scheduler)
        assert_divisions(b, scheduler=scheduler)
        if hasattr(a, "divisions") and hasattr(b, "divisions"):
            at = type(np.asarray(a.divisions).tolist()[0])  # numpy to python
            bt = type(np.asarray(b.divisions).tolist()[0])  # scalar conversion
            assert at == bt, (at, bt)
    assert_sane_keynames(a)
    assert_sane_keynames(b)
    a = _check_dask(
        a, check_names=check_names, check_dtypes=check_dtype, scheduler=scheduler
    )
    b = _check_dask(
        b, check_names=check_names, check_dtypes=check_dtype, scheduler=scheduler
    )
    if hasattr(a, "to_pandas"):
        a = a.to_pandas()
    if hasattr(b, "to_pandas"):
        b = b.to_pandas()

    a, b = _maybe_convert_string(a, b)

    if isinstance(a, (pd.DataFrame, pd.Series)) and sort_results:
        a = _maybe_sort(a, check_index)
        b = _maybe_sort(b, check_index)
    if not check_index:
        a = a.reset_index(drop=True)
        b = b.reset_index(drop=True)
    if isinstance(a, pd.DataFrame):
        tm.assert_frame_equal(
            a, b, check_names=check_names, check_dtype=check_dtype, **kwargs
        )
    elif isinstance(a, pd.Series):
        tm.assert_series_equal(
            a, b, check_names=check_names, check_dtype=check_dtype, **kwargs
        )
    elif isinstance(a, pd.Index):
        tm.assert_index_equal(a, b, exact=check_dtype, **kwargs)
    else:
        if a == b:
            return True
        else:
            if np.isnan(a):
                assert np.isnan(b)
            else:
                assert np.allclose(a, b)
    return True


def assert_dask_graph(dask, label):
    if hasattr(dask, "dask"):
        dask = dask.dask
    assert isinstance(dask, Mapping)
    for k in dask:
        if isinstance(k, tuple):
            k = k[0]
        if k.startswith(label):
            return True
    raise AssertionError(f"given dask graph doesn't contain label: {label}")


def assert_divisions(ddf, scheduler=None):
    if not hasattr(ddf, "divisions"):
        return

    assert isinstance(ddf.divisions, tuple)

    if not getattr(ddf, "known_divisions", False):
        return

    ddf.enforce_runtime_divisions().compute(scheduler=scheduler)


def assert_sane_keynames(ddf):
    if not hasattr(ddf, "dask"):
        return
    for k in ddf.dask.keys():
        while isinstance(k, tuple):
            k = k[0]
        assert isinstance(k, (str, bytes))
        assert len(k) < 100
        assert " " not in k
        assert k.split("-")[0].isidentifier(), k


def assert_dask_dtypes(ddf, res, numeric_equal=True):
    """Check that the dask metadata matches the result.

    If `numeric_equal`, integer and floating dtypes compare equal. This is
    useful due to the implicit conversion of integer to floating upon
    encountering missingness, which is hard to infer statically."""

    eq_type_sets = [{"O", "S", "U", "a"}]  # treat object and strings alike
    if numeric_equal:
        eq_type_sets.append({"i", "f", "u"})

    def eq_dtypes(a, b):
        return any(
            a.kind in eq_types and b.kind in eq_types for eq_types in eq_type_sets
        ) or (a == b)

    if not is_dask_collection(res) and is_dataframe_like(res):
        for a, b in pd.concat([ddf._meta.dtypes, res.dtypes], axis=1).itertuples(
            index=False
        ):
            assert eq_dtypes(a, b)
    elif not is_dask_collection(res) and (is_index_like(res) or is_series_like(res)):
        a = ddf._meta.dtype
        b = res.dtype
        assert eq_dtypes(a, b)
    else:
        if hasattr(ddf._meta, "dtype"):
            a = ddf._meta.dtype
            if not hasattr(res, "dtype"):
                assert np.isscalar(res)
                b = np.dtype(type(res))
            else:
                b = res.dtype
            assert eq_dtypes(a, b)
        else:
            assert type(ddf._meta) == type(res)


def assert_max_deps(x, n, eq=True):
    dependencies, dependents = get_deps(x.dask)
    if eq:
        assert max(map(len, dependencies.values())) == n
    else:
        assert max(map(len, dependencies.values())) <= n


def valid_divisions(divisions):
    """Are the provided divisions valid?

    Examples
    --------
    >>> valid_divisions([1, 2, 3])
    True
    >>> valid_divisions([3, 2, 1])
    False
    >>> valid_divisions([1, 1, 1])
    False
    >>> valid_divisions([0, 1, 1])
    True
    >>> valid_divisions((1, 2, 3))
    True
    >>> valid_divisions(123)
    False
    >>> valid_divisions([0, float('nan'), 1])
    False
    """
    if not isinstance(divisions, (tuple, list)):
        return False

    # Cast tuples to lists as `pd.isnull` treats them differently
    # https://github.com/pandas-dev/pandas/issues/52283
    if isinstance(divisions, tuple):
        divisions = list(divisions)

    if pd.isnull(divisions).any():
        return False

    for i, x in enumerate(divisions[:-2]):
        if x >= divisions[i + 1]:
            return False
        if isinstance(x, Number) and math.isnan(x):
            return False

    for x in divisions[-2:]:
        if isinstance(x, Number) and math.isnan(x):
            return False

    if divisions[-2] > divisions[-1]:
        return False

    return True


def drop_by_shallow_copy(df, columns, errors="raise"):
    """Use shallow copy to drop columns in place"""
    df2 = df.copy(deep=False)
    if not pd.api.types.is_list_like(columns):
        columns = [columns]
    df2.drop(columns=columns, inplace=True, errors=errors)
    return df2


class AttributeNotImplementedError(NotImplementedError, AttributeError):
    """NotImplementedError and AttributeError"""


def meta_frame_constructor(like):
    """Return a serial DataFrame constructor

    Parameters
    ----------
    like :
        Any series-like, Index-like or dataframe-like object.
    """
    if is_dask_collection(like):
        try:
            like = like._meta
        except AttributeError:
            raise TypeError(f"{type(like)} not supported by meta_frame_constructor")
    if is_dataframe_like(like):
        return like._constructor
    elif is_series_like(like):
        return like._constructor_expanddim
    elif is_index_like(like):
        return like.to_frame()._constructor
    else:
        raise TypeError(f"{type(like)} not supported by meta_frame_constructor")


def meta_series_constructor(like):
    """Return a serial Series constructor

    Parameters
    ----------
    like :
        Any series-like, Index-like or dataframe-like object.
    """
    if is_dask_collection(like):
        try:
            like = like._meta
        except AttributeError:
            raise TypeError(f"{type(like)} not supported by meta_series_constructor")
    if is_dataframe_like(like):
        return like._constructor_sliced
    elif is_series_like(like):
        return like._constructor
    elif is_index_like(like):
        return like.to_frame()._constructor_sliced
    else:
        raise TypeError(f"{type(like)} not supported by meta_series_constructor")


def get_string_dtype():
    """Depending on config setting, we might convert objects to pyarrow strings"""
    return pd.StringDtype("pyarrow") if pyarrow_strings_enabled() else object


def pyarrow_strings_enabled() -> bool:
    """Config setting to convert objects to pyarrow strings"""
    convert_string = dask.config.get("dataframe.convert-string")
    if convert_string is None:
        convert_string = True
    return convert_string


def get_numeric_only_kwargs(numeric_only: bool | NoDefault) -> dict:
    return {} if numeric_only is no_default else {"numeric_only": numeric_only}
