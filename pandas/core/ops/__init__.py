"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
import datetime
import operator
from typing import Any, Callable, Tuple, Union

import numpy as np

from pandas._libs import Timedelta, Timestamp, lib, ops as libops
from pandas.errors import NullFrequencyError
from pandas.util._decorators import Appender

from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike
from pandas.core.dtypes.common import (
    ensure_object,
    is_bool_dtype,
    is_datetime64_dtype,
    is_extension_array_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeArray,
    ABCDatetimeIndex,
    ABCExtensionArray,
    ABCIndexClass,
    ABCSeries,
    ABCTimedeltaArray,
    ABCTimedeltaIndex,
)
from pandas.core.dtypes.missing import isna, notna

from pandas._typing import ArrayLike
from pandas.core.construction import array, extract_array
from pandas.core.ops.array_ops import (
    comp_method_OBJECT_ARRAY,
    define_na_arithmetic_op,
    na_arithmetic_op,
)
from pandas.core.ops.docstrings import (
    _arith_doc_FRAME,
    _flex_comp_doc_FRAME,
    _make_flex_doc,
    _op_descriptions,
)
from pandas.core.ops.invalid import invalid_comparison
from pandas.core.ops.methods import (  # noqa:F401
    add_flex_arithmetic_methods,
    add_special_arithmetic_methods,
)
from pandas.core.ops.roperator import (  # noqa:F401
    radd,
    rand_,
    rdiv,
    rdivmod,
    rfloordiv,
    rmod,
    rmul,
    ror_,
    rpow,
    rsub,
    rtruediv,
    rxor,
)

# -----------------------------------------------------------------------------
# Ops Wrapping Utilities


def get_op_result_name(left, right):
    """
    Find the appropriate name to pin to an operation result.  This result
    should always be either an Index or a Series.

    Parameters
    ----------
    left : {Series, Index}
    right : object

    Returns
    -------
    name : object
        Usually a string
    """
    # `left` is always a Series when called from within ops
    if isinstance(right, (ABCSeries, ABCIndexClass)):
        name = _maybe_match_name(left, right)
    else:
        name = left.name
    return name


def _maybe_match_name(a, b):
    """
    Try to find a name to attach to the result of an operation between
    a and b.  If only one of these has a `name` attribute, return that
    name.  Otherwise return a consensus name if they match of None if
    they have different names.

    Parameters
    ----------
    a : object
    b : object

    Returns
    -------
    name : str or None

    See Also
    --------
    pandas.core.common.consensus_name_attr
    """
    a_has = hasattr(a, "name")
    b_has = hasattr(b, "name")
    if a_has and b_has:
        if a.name == b.name:
            return a.name
        else:
            # TODO: what if they both have np.nan for their names?
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None


def maybe_upcast_for_op(obj, shape: Tuple[int, ...]):
    """
    Cast non-pandas objects to pandas types to unify behavior of arithmetic
    and comparison operations.

    Parameters
    ----------
    obj: object
    shape : tuple[int]

    Returns
    -------
    out : object

    Notes
    -----
    Be careful to call this *after* determining the `name` attribute to be
    attached to the result of the arithmetic operation.
    """
    from pandas.core.arrays import DatetimeArray, TimedeltaArray

    if type(obj) is datetime.timedelta:
        # GH#22390  cast up to Timedelta to rely on Timedelta
        # implementation; otherwise operation against numeric-dtype
        # raises TypeError
        return Timedelta(obj)
    elif isinstance(obj, np.datetime64):
        # GH#28080 numpy casts integer-dtype to datetime64 when doing
        #  array[int] + datetime64, which we do not allow
        if isna(obj):
            # Avoid possible ambiguities with pd.NaT
            obj = obj.astype("datetime64[ns]")
            right = np.broadcast_to(obj, shape)
            return DatetimeArray(right)

        return Timestamp(obj)

    elif isinstance(obj, np.timedelta64):
        if isna(obj):
            # wrapping timedelta64("NaT") in Timedelta returns NaT,
            #  which would incorrectly be treated as a datetime-NaT, so
            #  we broadcast and wrap in a TimedeltaArray
            obj = obj.astype("timedelta64[ns]")
            right = np.broadcast_to(obj, shape)
            return TimedeltaArray(right)

        # In particular non-nanosecond timedelta64 needs to be cast to
        #  nanoseconds, or else we get undesired behavior like
        #  np.timedelta64(3, 'D') / 2 == np.timedelta64(1, 'D')
        return Timedelta(obj)

    elif isinstance(obj, np.ndarray) and is_timedelta64_dtype(obj.dtype):
        # GH#22390 Unfortunately we need to special-case right-hand
        # timedelta64 dtypes because numpy casts integer dtypes to
        # timedelta64 when operating with timedelta64
        return TimedeltaArray._from_sequence(obj)
    return obj


# -----------------------------------------------------------------------------


def _gen_eval_kwargs(name):
    """
    Find the keyword arguments to pass to numexpr for the given operation.

    Parameters
    ----------
    name : str

    Returns
    -------
    eval_kwargs : dict

    Examples
    --------
    >>> _gen_eval_kwargs("__add__")
    {}

    >>> _gen_eval_kwargs("rtruediv")
    {'reversed': True, 'truediv': True}
    """
    kwargs = {}

    # Series appear to only pass __add__, __radd__, ...
    # but DataFrame gets both these dunder names _and_ non-dunder names
    # add, radd, ...
    name = name.replace("__", "")

    if name.startswith("r"):
        if name not in ["radd", "rand", "ror", "rxor"]:
            # Exclude commutative operations
            kwargs["reversed"] = True

    return kwargs


def _get_frame_op_default_axis(name):
    """
    Only DataFrame cares about default_axis, specifically:
    special methods have default_axis=None and flex methods
    have default_axis='columns'.

    Parameters
    ----------
    name : str

    Returns
    -------
    default_axis: str or None
    """
    if name.replace("__r", "__") in ["__and__", "__or__", "__xor__"]:
        # bool methods
        return "columns"
    elif name.startswith("__"):
        # __add__, __mul__, ...
        return None
    else:
        # add, mul, ...
        return "columns"


def _get_opstr(op):
    """
    Find the operation string, if any, to pass to numexpr for this
    operation.

    Parameters
    ----------
    op : binary operator

    Returns
    -------
    op_str : string or None
    """

    return {
        operator.add: "+",
        radd: "+",
        operator.mul: "*",
        rmul: "*",
        operator.sub: "-",
        rsub: "-",
        operator.truediv: "/",
        rtruediv: "/",
        operator.floordiv: "//",
        rfloordiv: "//",
        operator.mod: None,  # TODO: Why None for mod but '%' for rmod?
        rmod: "%",
        operator.pow: "**",
        rpow: "**",
        operator.eq: "==",
        operator.ne: "!=",
        operator.le: "<=",
        operator.lt: "<",
        operator.ge: ">=",
        operator.gt: ">",
        operator.and_: "&",
        rand_: "&",
        operator.or_: "|",
        ror_: "|",
        operator.xor: "^",
        rxor: "^",
        divmod: None,
        rdivmod: None,
    }[op]


def _get_op_name(op, special):
    """
    Find the name to attach to this method according to conventions
    for special and non-special methods.

    Parameters
    ----------
    op : binary operator
    special : bool

    Returns
    -------
    op_name : str
    """
    opname = op.__name__.strip("_")
    if special:
        opname = "__{opname}__".format(opname=opname)
    return opname


# -----------------------------------------------------------------------------
# Masking NA values and fallbacks for operations numpy does not support


def fill_binop(left, right, fill_value):
    """
    If a non-None fill_value is given, replace null entries in left and right
    with this value, but only in positions where _one_ of left/right is null,
    not both.

    Parameters
    ----------
    left : array-like
    right : array-like
    fill_value : object

    Returns
    -------
    left : array-like
    right : array-like

    Notes
    -----
    Makes copies if fill_value is not None
    """
    # TODO: can we make a no-copy implementation?
    if fill_value is not None:
        left_mask = isna(left)
        right_mask = isna(right)
        left = left.copy()
        right = right.copy()

        # one but not both
        mask = left_mask ^ right_mask
        left[left_mask & mask] = fill_value
        right[right_mask & mask] = fill_value
    return left, right


def mask_cmp_op(x, y, op):
    """
    Apply the function `op` to only non-null points in x and y.

    Parameters
    ----------
    x : array-like
    y : array-like
    op : binary operation

    Returns
    -------
    result : ndarray[bool]
    """
    xrav = x.ravel()
    result = np.empty(x.size, dtype=bool)
    if isinstance(y, (np.ndarray, ABCSeries)):
        yrav = y.ravel()
        mask = notna(xrav) & notna(yrav)
        result[mask] = op(np.array(list(xrav[mask])), np.array(list(yrav[mask])))
    else:
        mask = notna(xrav)
        result[mask] = op(np.array(list(xrav[mask])), y)

    if op == operator.ne:  # pragma: no cover
        np.putmask(result, ~mask, True)
    else:
        np.putmask(result, ~mask, False)
    result = result.reshape(x.shape)
    return result


# -----------------------------------------------------------------------------
# Dispatch logic


def should_extension_dispatch(left: ABCSeries, right: Any) -> bool:
    """
    Identify cases where Series operation should use dispatch_to_extension_op.

    Parameters
    ----------
    left : Series
    right : object

    Returns
    -------
    bool
    """
    if (
        is_extension_array_dtype(left.dtype)
        or is_datetime64_dtype(left.dtype)
        or is_timedelta64_dtype(left.dtype)
    ):
        return True

    if not is_scalar(right) and is_extension_array_dtype(right):
        # GH#22378 disallow scalar to exclude e.g. "category", "Int64"
        return True

    return False


def should_series_dispatch(left, right, op):
    """
    Identify cases where a DataFrame operation should dispatch to its
    Series counterpart.

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    op : binary operator

    Returns
    -------
    override : bool
    """
    if left._is_mixed_type or right._is_mixed_type:
        return True

    if not len(left.columns) or not len(right.columns):
        # ensure obj.dtypes[0] exists for each obj
        return False

    ldtype = left.dtypes.iloc[0]
    rdtype = right.dtypes.iloc[0]

    if (is_timedelta64_dtype(ldtype) and is_integer_dtype(rdtype)) or (
        is_timedelta64_dtype(rdtype) and is_integer_dtype(ldtype)
    ):
        # numpy integer dtypes as timedelta64 dtypes in this scenario
        return True

    if is_datetime64_dtype(ldtype) and is_object_dtype(rdtype):
        # in particular case where right is an array of DateOffsets
        return True

    return False


def dispatch_to_series(left, right, func, str_rep=None, axis=None):
    """
    Evaluate the frame operation func(left, right) by evaluating
    column-by-column, dispatching to the Series implementation.

    Parameters
    ----------
    left : DataFrame
    right : scalar or DataFrame
    func : arithmetic or comparison operator
    str_rep : str or None, default None
    axis : {None, 0, 1, "index", "columns"}

    Returns
    -------
    DataFrame
    """
    # Note: we use iloc to access columns for compat with cases
    #       with non-unique columns.
    import pandas.core.computation.expressions as expressions

    right = lib.item_from_zerodim(right)
    if lib.is_scalar(right) or np.ndim(right) == 0:

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b) for i in range(len(a.columns))}

    elif isinstance(right, ABCDataFrame):
        assert right._indexed_same(left)

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b.iloc[:, i]) for i in range(len(a.columns))}

    elif isinstance(right, ABCSeries) and axis == "columns":
        # We only get here if called via left._combine_match_columns,
        # in which case we specifically want to operate row-by-row
        assert right.index.equals(left.columns)

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b.iloc[i]) for i in range(len(a.columns))}

    elif isinstance(right, ABCSeries):
        assert right.index.equals(left.index)  # Handle other cases later

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b) for i in range(len(a.columns))}

    else:
        # Remaining cases have less-obvious dispatch rules
        raise NotImplementedError(right)

    new_data = expressions.evaluate(column_op, str_rep, left, right)
    return new_data


def dispatch_to_extension_op(
    op,
    left: Union[ABCExtensionArray, np.ndarray],
    right: Any,
    keep_null_freq: bool = False,
):
    """
    Assume that left or right is a Series backed by an ExtensionArray,
    apply the operator defined by op.

    Parameters
    ----------
    op : binary operator
    left : ExtensionArray or np.ndarray
    right : object
    keep_null_freq : bool, default False
        Whether to re-raise a NullFrequencyError unchanged, as opposed to
        catching and raising TypeError.

    Returns
    -------
    ExtensionArray or np.ndarray
        2-tuple of these if op is divmod or rdivmod
    """
    # NB: left and right should already be unboxed, so neither should be
    #  a Series or Index.

    if left.dtype.kind in "mM" and isinstance(left, np.ndarray):
        # We need to cast datetime64 and timedelta64 ndarrays to
        #  DatetimeArray/TimedeltaArray.  But we avoid wrapping others in
        #  PandasArray as that behaves poorly with e.g. IntegerArray.
        left = array(left)

    # The op calls will raise TypeError if the op is not defined
    # on the ExtensionArray

    try:
        res_values = op(left, right)
    except NullFrequencyError:
        # DatetimeIndex and TimedeltaIndex with freq == None raise ValueError
        # on add/sub of integers (or int-like).  We re-raise as a TypeError.
        if keep_null_freq:
            # TODO: remove keep_null_freq after Timestamp+int deprecation
            #  GH#22535 is enforced
            raise
        raise TypeError(
            "incompatible type for a datetime/timedelta "
            "operation [{name}]".format(name=op.__name__)
        )
    return res_values


# -----------------------------------------------------------------------------
# Series


def _align_method_SERIES(left, right, align_asobject=False):
    """ align lhs and rhs Series """

    # ToDo: Different from _align_method_FRAME, list, tuple and ndarray
    # are not coerced here
    # because Series has inconsistencies described in #13637

    if isinstance(right, ABCSeries):
        # avoid repeated alignment
        if not left.index.equals(right.index):

            if align_asobject:
                # to keep original value's dtype for bool ops
                left = left.astype(object)
                right = right.astype(object)

            left, right = left.align(right, copy=False)

    return left, right


def _construct_result(left, result, index, name, dtype=None):
    """
    If the raw op result has a non-None name (e.g. it is an Index object) and
    the name argument is None, then passing name to the constructor will
    not be enough; we still need to override the name attribute.
    """
    out = left._constructor(result, index=index, dtype=dtype)
    out = out.__finalize__(left)
    out.name = name
    return out


def _construct_divmod_result(left, result, index, name, dtype=None):
    """divmod returns a tuple of like indexed series instead of a single series.
    """
    return (
        _construct_result(left, result[0], index=index, name=name, dtype=dtype),
        _construct_result(left, result[1], index=index, name=name, dtype=dtype),
    )


def _arith_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)
    eval_kwargs = _gen_eval_kwargs(op_name)
    construct_result = (
        _construct_divmod_result if op in [divmod, rdivmod] else _construct_result
    )

    def wrapper(left, right):
        if isinstance(right, ABCDataFrame):
            return NotImplemented

        left, right = _align_method_SERIES(left, right)
        res_name = get_op_result_name(left, right)

        keep_null_freq = isinstance(
            right,
            (
                ABCDatetimeIndex,
                ABCDatetimeArray,
                ABCTimedeltaIndex,
                ABCTimedeltaArray,
                Timestamp,
            ),
        )

        lvalues = extract_array(left, extract_numpy=True)
        rvalues = extract_array(right, extract_numpy=True)

        rvalues = maybe_upcast_for_op(rvalues, lvalues.shape)

        if should_extension_dispatch(left, rvalues) or isinstance(
            rvalues, (ABCTimedeltaArray, ABCDatetimeArray, Timestamp)
        ):
            result = dispatch_to_extension_op(op, lvalues, rvalues, keep_null_freq)

        else:
            with np.errstate(all="ignore"):
                result = na_arithmetic_op(lvalues, rvalues, op, str_rep, eval_kwargs)

        # We do not pass dtype to ensure that the Series constructor
        #  does inference in the case where `result` has object-dtype.
        return construct_result(left, result, index=left.index, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _comp_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    def wrapper(self, other):

        res_name = get_op_result_name(self, other)

        # TODO: shouldn't we be applying finalize whenever
        #  not isinstance(other, ABCSeries)?
        finalizer = (
            lambda x: x.__finalize__(self)
            if isinstance(other, (np.ndarray, ABCIndexClass))
            else x
        )

        if isinstance(other, ABCDataFrame):  # pragma: no cover
            # Defer to DataFrame implementation; fail early
            return NotImplemented

        if isinstance(other, ABCSeries) and not self._indexed_same(other):
            raise ValueError("Can only compare identically-labeled Series objects")

        other = lib.item_from_zerodim(other)
        if isinstance(other, list):
            # TODO: same for tuples?
            other = np.asarray(other)

        if isinstance(other, (np.ndarray, ABCExtensionArray, ABCIndexClass)):
            # TODO: make this treatment consistent across ops and classes.
            #  We are not catching all listlikes here (e.g. frozenset, tuple)
            #  The ambiguous case is object-dtype.  See GH#27803
            if len(self) != len(other):
                raise ValueError("Lengths must match to compare")

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        if should_extension_dispatch(lvalues, rvalues):
            res_values = dispatch_to_extension_op(op, lvalues, rvalues)

        elif is_scalar(rvalues) and isna(rvalues):
            # numpy does not like comparisons vs None
            if op is operator.ne:
                res_values = np.ones(len(lvalues), dtype=bool)
            else:
                res_values = np.zeros(len(lvalues), dtype=bool)

        elif is_object_dtype(lvalues.dtype):
            res_values = comp_method_OBJECT_ARRAY(op, lvalues, rvalues)

        else:
            op_name = "__{op}__".format(op=op.__name__)
            method = getattr(lvalues, op_name)
            with np.errstate(all="ignore"):
                res_values = method(rvalues)

            if res_values is NotImplemented:
                res_values = invalid_comparison(lvalues, rvalues, op)
            if is_scalar(res_values):
                raise TypeError(
                    "Could not compare {typ} type with Series".format(typ=type(rvalues))
                )

        result = self._constructor(res_values, index=self.index)
        result = finalizer(result)

        # Set the result's name after finalizer is called because finalizer
        #  would set it back to self.name
        result.name = res_name
        return result

    wrapper.__name__ = op_name
    return wrapper


def _bool_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            assert not isinstance(y, (list, ABCSeries, ABCIndexClass))
            if isinstance(y, np.ndarray):
                # bool-bool dtype operations should be OK, should not get here
                assert not (is_bool_dtype(x.dtype) and is_bool_dtype(y.dtype))
                x = ensure_object(x)
                y = ensure_object(y)
                result = libops.vec_binop(x, y, op)
            else:
                # let null fall thru
                assert lib.is_scalar(y)
                if not isna(y):
                    y = bool(y)
                try:
                    result = libops.scalar_binop(x, y, op)
                except (
                    TypeError,
                    ValueError,
                    AttributeError,
                    OverflowError,
                    NotImplementedError,
                ):
                    raise TypeError(
                        "cannot compare a dtyped [{dtype}] array "
                        "with a scalar of type [{typ}]".format(
                            dtype=x.dtype, typ=type(y).__name__
                        )
                    )

        return result

    fill_int = lambda x: x

    def fill_bool(x, left=None):
        # if `left` is specifically not-boolean, we do not cast to bool
        if x.dtype.kind in ["c", "f", "O"]:
            # dtypes that can hold NA
            mask = isna(x)
            if mask.any():
                x = x.astype(object)
                x[mask] = False

        if left is None or is_bool_dtype(left.dtype):
            x = x.astype(bool)
        return x

    def wrapper(self, other):
        is_self_int_dtype = is_integer_dtype(self.dtype)

        self, other = _align_method_SERIES(self, other, align_asobject=True)
        res_name = get_op_result_name(self, other)

        # TODO: shouldn't we be applying finalize whenever
        #  not isinstance(other, ABCSeries)?
        finalizer = (
            lambda x: x.__finalize__(self)
            if not isinstance(other, (ABCSeries, ABCIndexClass))
            else x
        )

        if isinstance(other, ABCDataFrame):
            # Defer to DataFrame implementation; fail early
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if is_list_like(other) and not hasattr(other, "dtype"):
            # e.g. list, tuple
            other = construct_1d_object_array_from_listlike(other)

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        if should_extension_dispatch(self, rvalues):
            res_values = dispatch_to_extension_op(op, lvalues, rvalues)

        else:
            if isinstance(rvalues, (ABCSeries, ABCIndexClass, np.ndarray)):
                is_other_int_dtype = is_integer_dtype(rvalues.dtype)
                rvalues = rvalues if is_other_int_dtype else fill_bool(rvalues, lvalues)

            else:
                # i.e. scalar
                is_other_int_dtype = lib.is_integer(rvalues)

            # For int vs int `^`, `|`, `&` are bitwise operators and return
            #   integer dtypes.  Otherwise these are boolean ops
            filler = fill_int if is_self_int_dtype and is_other_int_dtype else fill_bool

            res_values = na_op(lvalues, rvalues)
            res_values = filler(res_values)

        result = self._constructor(res_values, index=self.index, name=res_name)
        return finalizer(result)

    wrapper.__name__ = op_name
    return wrapper


def _flex_method_SERIES(cls, op, special):
    name = _get_op_name(op, special)
    doc = _make_flex_doc(name, "series")

    @Appender(doc)
    def flex_wrapper(self, other, level=None, fill_value=None, axis=0):
        # validate axis
        if axis is not None:
            self._get_axis_number(axis)
        if isinstance(other, ABCSeries):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple)):
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
            other = self._constructor(other, self.index)
            return self._binop(other, op, level=level, fill_value=fill_value)
        else:
            if fill_value is not None:
                self = self.fillna(fill_value)

            return self._constructor(op(self, other), self.index).__finalize__(self)

    flex_wrapper.__name__ = name
    return flex_wrapper


# -----------------------------------------------------------------------------
# DataFrame


def _combine_series_frame(self, other, func, fill_value=None, axis=None, level=None):
    """
    Apply binary operator `func` to self, other using alignment and fill
    conventions determined by the fill_value, axis, and level kwargs.

    Parameters
    ----------
    self : DataFrame
    other : Series
    func : binary operator
    fill_value : object, default None
    axis : {0, 1, 'columns', 'index', None}, default None
    level : int or None, default None

    Returns
    -------
    result : DataFrame
    """
    if fill_value is not None:
        raise NotImplementedError(
            "fill_value {fill} not supported.".format(fill=fill_value)
        )

    if axis is not None:
        axis = self._get_axis_number(axis)
        if axis == 0:
            return self._combine_match_index(other, func, level=level)
        else:
            return self._combine_match_columns(other, func, level=level)
    else:
        if not len(other):
            return self * np.nan

        if not len(self):
            # Ambiguous case, use _series so works with DataFrame
            return self._constructor(
                data=self._series, index=self.index, columns=self.columns
            )

        # default axis is columns
        return self._combine_match_columns(other, func, level=level)


def _align_method_FRAME(left, right, axis):
    """ convert rhs to meet lhs dims if input is list, tuple or np.ndarray """

    def to_series(right):
        msg = "Unable to coerce to Series, length must be {req_len}: given {given_len}"
        if axis is not None and left._get_axis_name(axis) == "index":
            if len(left.index) != len(right):
                raise ValueError(
                    msg.format(req_len=len(left.index), given_len=len(right))
                )
            right = left._constructor_sliced(right, index=left.index)
        else:
            if len(left.columns) != len(right):
                raise ValueError(
                    msg.format(req_len=len(left.columns), given_len=len(right))
                )
            right = left._constructor_sliced(right, index=left.columns)
        return right

    if isinstance(right, np.ndarray):

        if right.ndim == 1:
            right = to_series(right)

        elif right.ndim == 2:
            if right.shape == left.shape:
                right = left._constructor(right, index=left.index, columns=left.columns)

            elif right.shape[0] == left.shape[0] and right.shape[1] == 1:
                # Broadcast across columns
                right = np.broadcast_to(right, left.shape)
                right = left._constructor(right, index=left.index, columns=left.columns)

            elif right.shape[1] == left.shape[1] and right.shape[0] == 1:
                # Broadcast along rows
                right = to_series(right[0, :])

            else:
                raise ValueError(
                    "Unable to coerce to DataFrame, shape "
                    "must be {req_shape}: given {given_shape}".format(
                        req_shape=left.shape, given_shape=right.shape
                    )
                )

        elif right.ndim > 2:
            raise ValueError(
                "Unable to coerce to Series/DataFrame, dim "
                "must be <= 2: {dim}".format(dim=right.shape)
            )

    elif is_list_like(right) and not isinstance(right, (ABCSeries, ABCDataFrame)):
        # GH17901
        right = to_series(right)

    return right


def _arith_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)
    eval_kwargs = _gen_eval_kwargs(op_name)
    default_axis = _get_frame_op_default_axis(op_name)

    na_op = define_na_arithmetic_op(op, str_rep, eval_kwargs)

    if op_name in _op_descriptions:
        # i.e. include "add" but not "__add__"
        doc = _make_flex_doc(op_name, "dataframe")
    else:
        doc = _arith_doc_FRAME % op_name

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None, fill_value=None):

        other = _align_method_FRAME(self, other, axis)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            pass_op = op if should_series_dispatch(self, other, op) else na_op
            return self._combine_frame(other, pass_op, fill_value, level)
        elif isinstance(other, ABCSeries):
            # For these values of `axis`, we end up dispatching to Series op,
            # so do not want the masked op.
            pass_op = op if axis in [0, "columns", None] else na_op
            return _combine_series_frame(
                self, other, pass_op, fill_value=fill_value, axis=axis, level=level
            )
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                self = self.fillna(fill_value)

            return self._combine_const(other, op)

    f.__name__ = op_name

    return f


def _flex_comp_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)

    def na_op(x, y):
        try:
            with np.errstate(invalid="ignore"):
                result = op(x, y)
        except TypeError:
            result = mask_cmp_op(x, y, op)
        return result

    doc = _flex_comp_doc_FRAME.format(
        op_name=op_name, desc=_op_descriptions[op_name]["desc"]
    )

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None):

        other = _align_method_FRAME(self, other, axis)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            if not self._indexed_same(other):
                self, other = self.align(other, "outer", level=level, copy=False)
            new_data = dispatch_to_series(self, other, na_op, str_rep)
            return self._construct_result(other, new_data, na_op)

        elif isinstance(other, ABCSeries):
            return _combine_series_frame(
                self, other, na_op, fill_value=None, axis=axis, level=level
            )
        else:
            # in this case we always have `np.ndim(other) == 0`
            return self._combine_const(other, na_op)

    f.__name__ = op_name

    return f


def _comp_method_FRAME(cls, func, special):
    str_rep = _get_opstr(func)
    op_name = _get_op_name(func, special)

    @Appender("Wrapper for comparison method {name}".format(name=op_name))
    def f(self, other):

        other = _align_method_FRAME(self, other, axis=None)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            if not self._indexed_same(other):
                raise ValueError(
                    "Can only compare identically-labeled DataFrame objects"
                )
            new_data = dispatch_to_series(self, other, func, str_rep)
            return self._construct_result(other, new_data, func)

        elif isinstance(other, ABCSeries):
            return _combine_series_frame(
                self, other, func, fill_value=None, axis=None, level=None
            )
        else:

            # straight boolean comparisons we want to allow all columns
            # (regardless of dtype to pass thru) See #4537 for discussion.
            res = self._combine_const(other, func)
            return res

    f.__name__ = op_name

    return f


# -----------------------------------------------------------------------------
# Sparse


def maybe_dispatch_ufunc_to_dunder_op(
    self: ArrayLike, ufunc: Callable, method: str, *inputs: ArrayLike, **kwargs: Any
):
    """
    Dispatch a ufunc to the equivalent dunder method.

    Parameters
    ----------
    self : ArrayLike
        The array whose dunder method we dispatch to
    ufunc : Callable
        A NumPy ufunc
    method : {'reduce', 'accumulate', 'reduceat', 'outer', 'at', '__call__'}
    inputs : ArrayLike
        The input arrays.
    kwargs : Any
        The additional keyword arguments, e.g. ``out``.

    Returns
    -------
    result : Any
        The result of applying the ufunc
    """
    # special has the ufuncs we dispatch to the dunder op on
    special = {
        "add",
        "sub",
        "mul",
        "pow",
        "mod",
        "floordiv",
        "truediv",
        "divmod",
        "eq",
        "ne",
        "lt",
        "gt",
        "le",
        "ge",
        "remainder",
        "matmul",
    }
    aliases = {
        "subtract": "sub",
        "multiply": "mul",
        "floor_divide": "floordiv",
        "true_divide": "truediv",
        "power": "pow",
        "remainder": "mod",
        "divide": "div",
        "equal": "eq",
        "not_equal": "ne",
        "less": "lt",
        "less_equal": "le",
        "greater": "gt",
        "greater_equal": "ge",
    }

    # For op(., Array) -> Array.__r{op}__
    flipped = {
        "lt": "__gt__",
        "le": "__ge__",
        "gt": "__lt__",
        "ge": "__le__",
        "eq": "__eq__",
        "ne": "__ne__",
    }

    op_name = ufunc.__name__
    op_name = aliases.get(op_name, op_name)

    def not_implemented(*args, **kwargs):
        return NotImplemented

    if method == "__call__" and op_name in special and kwargs.get("out") is None:
        if isinstance(inputs[0], type(self)):
            name = "__{}__".format(op_name)
            return getattr(self, name, not_implemented)(inputs[1])
        else:
            name = flipped.get(op_name, "__r{}__".format(op_name))
            return getattr(self, name, not_implemented)(inputs[0])
    else:
        return NotImplemented
