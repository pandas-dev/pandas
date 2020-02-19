"""
Helpers for making specific test assertions.
"""
from typing import Union, cast

import numpy as np

import pandas._libs.testing as _testing

from pandas.core.dtypes.common import (
    is_bool,
    is_categorical_dtype,
    is_datetime64tz_dtype,
    is_extension_array_dtype,
    is_interval_dtype,
    is_number,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import array_equivalent

import pandas as pd
from pandas import Categorical, DataFrame, Index, MultiIndex, Series
from pandas.core.algorithms import take_1d
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    IntervalArray,
    PeriodArray,
    TimedeltaArray,
)

from pandas.io.formats.printing import pprint_thing


def _check_isinstance(left, right, cls):
    """
    Helper method for our assert_* methods that ensures that
    the two objects being compared have the right type before
    proceeding with the comparison.

    Parameters
    ----------
    left : The first object being compared.
    right : The second object being compared.
    cls : The class type to check against.

    Raises
    ------
    AssertionError : Either `left` or `right` is not an instance of `cls`.
    """
    cls_name = cls.__name__

    if not isinstance(left, cls):
        raise AssertionError(
            f"{cls_name} Expected type {cls}, found {type(left)} instead"
        )
    if not isinstance(right, cls):
        raise AssertionError(
            f"{cls_name} Expected type {cls}, found {type(right)} instead"
        )


def assert_almost_equal(
    left,
    right,
    check_dtype: Union[bool, str] = "equiv",
    check_less_precise: Union[bool, int] = False,
    **kwargs,
):
    """
    Check that the left and right objects are approximately equal.

    By approximately equal, we refer to objects that are numbers or that
    contain numbers which may be equivalent to specific levels of precision.

    Parameters
    ----------
    left : object
    right : object
    check_dtype : bool or {'equiv'}, default 'equiv'
        Check dtype if both a and b are the same type. If 'equiv' is passed in,
        then `RangeIndex` and `Int64Index` are also considered equivalent
        when doing type checking.
    check_less_precise : bool or int, default False
        Specify comparison precision. 5 digits (False) or 3 digits (True)
        after decimal points are compared. If int, then specify the number
        of digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    """
    if isinstance(left, pd.Index):
        assert_index_equal(
            left,
            right,
            check_exact=False,
            exact=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    elif isinstance(left, pd.Series):
        assert_series_equal(
            left,
            right,
            check_exact=False,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    elif isinstance(left, pd.DataFrame):
        assert_frame_equal(
            left,
            right,
            check_exact=False,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    else:
        # Other sequences.
        if check_dtype:
            if is_number(left) and is_number(right):
                # Do not compare numeric classes, like np.float64 and float.
                pass
            elif is_bool(left) and is_bool(right):
                # Do not compare bool classes, like np.bool_ and bool.
                pass
            else:
                if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
                    obj = "numpy array"
                else:
                    obj = "Input"
                assert_class_equal(left, right, obj=obj)
        _testing.assert_almost_equal(
            left,
            right,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )


def assert_dict_equal(left, right, compare_keys: bool = True):

    _check_isinstance(left, right, dict)
    _testing.assert_dict_equal(left, right, compare_keys=compare_keys)


def assert_index_equal(
    left: Index,
    right: Index,
    exact: Union[bool, str] = "equiv",
    check_names: bool = True,
    check_less_precise: Union[bool, int] = False,
    check_exact: bool = True,
    check_categorical: bool = True,
    obj: str = "Index",
) -> None:
    """
    Check that left and right Index are equal.

    Parameters
    ----------
    left : Index
    right : Index
    exact : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    check_names : bool, default True
        Whether to check the names attribute.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_exact : bool, default True
        Whether to compare number exactly.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    obj : str, default 'Index'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    def _check_types(l, r, obj="Index"):
        if exact:
            assert_class_equal(l, r, exact=exact, obj=obj)

            # Skip exact dtype checking when `check_categorical` is False
            if check_categorical:
                assert_attr_equal("dtype", l, r, obj=obj)

            # allow string-like to have different inferred_types
            if l.inferred_type in ("string"):
                assert r.inferred_type in ("string")
            else:
                assert_attr_equal("inferred_type", l, r, obj=obj)

    def _get_ilevel_values(index, level):
        # accept level number only
        unique = index.levels[level]
        level_codes = index.codes[level]
        filled = take_1d(unique._values, level_codes, fill_value=unique._na_value)
        values = unique._shallow_copy(filled, name=index.names[level])
        return values

    # instance validation
    _check_isinstance(left, right, Index)

    # class / dtype comparison
    _check_types(left, right, obj=obj)

    # level comparison
    if left.nlevels != right.nlevels:
        msg1 = f"{obj} levels are different"
        msg2 = f"{left.nlevels}, {left}"
        msg3 = f"{right.nlevels}, {right}"
        raise_assert_detail(obj, msg1, msg2, msg3)

    # length comparison
    if len(left) != len(right):
        msg1 = f"{obj} length are different"
        msg2 = f"{len(left)}, {left}"
        msg3 = f"{len(right)}, {right}"
        raise_assert_detail(obj, msg1, msg2, msg3)

    # MultiIndex special comparison for little-friendly error messages
    if left.nlevels > 1:
        left = cast(MultiIndex, left)
        right = cast(MultiIndex, right)

        for level in range(left.nlevels):
            # cannot use get_level_values here because it can change dtype
            llevel = _get_ilevel_values(left, level)
            rlevel = _get_ilevel_values(right, level)

            lobj = f"MultiIndex level [{level}]"
            assert_index_equal(
                llevel,
                rlevel,
                exact=exact,
                check_names=check_names,
                check_less_precise=check_less_precise,
                check_exact=check_exact,
                obj=lobj,
            )
            # get_level_values may change dtype
            _check_types(left.levels[level], right.levels[level], obj=obj)

    # skip exact index checking when `check_categorical` is False
    if check_exact and check_categorical:
        if not left.equals(right):
            diff = np.sum((left.values != right.values).astype(int)) * 100.0 / len(left)
            msg = f"{obj} values are different ({np.round(diff, 5)} %)"
            raise_assert_detail(obj, msg, left, right)
    else:
        _testing.assert_almost_equal(
            left.values,
            right.values,
            check_less_precise=check_less_precise,
            check_dtype=exact,
            obj=obj,
            lobj=left,
            robj=right,
        )

    # metadata comparison
    if check_names:
        assert_attr_equal("names", left, right, obj=obj)
    if isinstance(left, pd.PeriodIndex) or isinstance(right, pd.PeriodIndex):
        assert_attr_equal("freq", left, right, obj=obj)
    if isinstance(left, pd.IntervalIndex) or isinstance(right, pd.IntervalIndex):
        assert_interval_array_equal(left.values, right.values)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(left.values, right.values, obj=f"{obj} category")


def assert_class_equal(left, right, exact: Union[bool, str] = True, obj="Input"):
    """
    Checks classes are equal.
    """
    __tracebackhide__ = True

    def repr_class(x):
        if isinstance(x, Index):
            # return Index as it is to include values in the error message
            return x

        try:
            return type(x).__name__
        except AttributeError:
            return repr(type(x))

    if exact == "equiv":
        if type(left) != type(right):
            # allow equivalence of Int64Index/RangeIndex
            types = {type(left).__name__, type(right).__name__}
            if len(types - {"Int64Index", "RangeIndex"}):
                msg = f"{obj} classes are not equivalent"
                raise_assert_detail(obj, msg, repr_class(left), repr_class(right))
    elif exact:
        if type(left) != type(right):
            msg = f"{obj} classes are different"
            raise_assert_detail(obj, msg, repr_class(left), repr_class(right))


def assert_attr_equal(attr, left, right, obj="Attributes"):
    """
    checks attributes are equal. Both objects must have attribute.

    Parameters
    ----------
    attr : str
        Attribute name being compared.
    left : object
    right : object
    obj : str, default 'Attributes'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    __tracebackhide__ = True

    left_attr = getattr(left, attr)
    right_attr = getattr(right, attr)

    if left_attr is right_attr:
        return True
    elif (
        is_number(left_attr)
        and np.isnan(left_attr)
        and is_number(right_attr)
        and np.isnan(right_attr)
    ):
        # np.nan
        return True

    try:
        result = left_attr == right_attr
    except TypeError:
        # datetimetz on rhs may raise TypeError
        result = False
    if not isinstance(result, bool):
        result = result.all()

    if result:
        return True
    else:
        msg = f'Attribute "{attr}" are different'
        raise_assert_detail(obj, msg, left_attr, right_attr)


def assert_categorical_equal(
    left, right, check_dtype=True, check_category_order=True, obj="Categorical"
):
    """
    Test that Categoricals are equivalent.

    Parameters
    ----------
    left : Categorical
    right : Categorical
    check_dtype : bool, default True
        Check that integer dtype of the codes are the same
    check_category_order : bool, default True
        Whether the order of the categories should be compared, which
        implies identical integer codes.  If False, only the resulting
        values are compared.  The ordered attribute is
        checked regardless.
    obj : str, default 'Categorical'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    _check_isinstance(left, right, Categorical)

    if check_category_order:
        assert_index_equal(left.categories, right.categories, obj=f"{obj}.categories")
        assert_numpy_array_equal(
            left.codes, right.codes, check_dtype=check_dtype, obj=f"{obj}.codes",
        )
    else:
        assert_index_equal(
            left.categories.sort_values(),
            right.categories.sort_values(),
            obj=f"{obj}.categories",
        )
        assert_index_equal(
            left.categories.take(left.codes),
            right.categories.take(right.codes),
            obj=f"{obj}.values",
        )

    assert_attr_equal("ordered", left, right, obj=obj)


def assert_interval_array_equal(left, right, exact="equiv", obj="IntervalArray"):
    """
    Test that two IntervalArrays are equivalent.

    Parameters
    ----------
    left, right : IntervalArray
        The IntervalArrays to compare.
    exact : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    obj : str, default 'IntervalArray'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    _check_isinstance(left, right, IntervalArray)

    assert_index_equal(left.left, right.left, exact=exact, obj=f"{obj}.left")
    assert_index_equal(left.right, right.right, exact=exact, obj=f"{obj}.left")
    assert_attr_equal("closed", left, right, obj=obj)


def assert_period_array_equal(left, right, obj="PeriodArray"):
    _check_isinstance(left, right, PeriodArray)

    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}.values")
    assert_attr_equal("freq", left, right, obj=obj)


def assert_datetime_array_equal(left, right, obj="DatetimeArray"):
    __tracebackhide__ = True
    _check_isinstance(left, right, DatetimeArray)

    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}._data")
    assert_attr_equal("freq", left, right, obj=obj)
    assert_attr_equal("tz", left, right, obj=obj)


def assert_timedelta_array_equal(left, right, obj="TimedeltaArray"):
    __tracebackhide__ = True
    _check_isinstance(left, right, TimedeltaArray)
    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}._data")
    assert_attr_equal("freq", left, right, obj=obj)


def assert_sp_array_equal(
    left,
    right,
    check_dtype=True,
    check_kind=True,
    check_fill_value=True,
    consolidate_block_indices=False,
):
    """
    Check that the left and right SparseArray are equal.

    Parameters
    ----------
    left : SparseArray
    right : SparseArray
    check_dtype : bool, default True
        Whether to check the data dtype is identical.
    check_kind : bool, default True
        Whether to just the kind of the sparse index for each column.
    check_fill_value : bool, default True
        Whether to check that left.fill_value matches right.fill_value
    consolidate_block_indices : bool, default False
        Whether to consolidate contiguous blocks for sparse arrays with
        a BlockIndex. Some operations, e.g. concat, will end up with
        block indices that could be consolidated. Setting this to true will
        create a new BlockIndex for that array, with consolidated
        block indices.
    """
    _check_isinstance(left, right, pd.arrays.SparseArray)

    assert_numpy_array_equal(left.sp_values, right.sp_values, check_dtype=check_dtype)

    # SparseIndex comparison
    assert isinstance(left.sp_index, pd._libs.sparse.SparseIndex)
    assert isinstance(right.sp_index, pd._libs.sparse.SparseIndex)

    if not check_kind:
        left_index = left.sp_index.to_block_index()
        right_index = right.sp_index.to_block_index()
    else:
        left_index = left.sp_index
        right_index = right.sp_index

    if consolidate_block_indices and left.kind == "block":
        # we'll probably remove this hack...
        left_index = left_index.to_int_index().to_block_index()
        right_index = right_index.to_int_index().to_block_index()

    if not left_index.equals(right_index):
        raise_assert_detail(
            "SparseArray.index", "index are not equal", left_index, right_index
        )
    else:
        # Just ensure a
        pass

    if check_fill_value:
        assert_attr_equal("fill_value", left, right)
    if check_dtype:
        assert_attr_equal("dtype", left, right)
    assert_numpy_array_equal(left.to_dense(), right.to_dense(), check_dtype=check_dtype)


def raise_assert_detail(obj, message, left, right, diff=None):
    __tracebackhide__ = True

    if isinstance(left, np.ndarray):
        left = pprint_thing(left)
    elif is_categorical_dtype(left):
        left = repr(left)

    if isinstance(right, np.ndarray):
        right = pprint_thing(right)
    elif is_categorical_dtype(right):
        right = repr(right)

    msg = f"""{obj} are different

{message}
[left]:  {left}
[right]: {right}"""

    if diff is not None:
        msg += f"\n[diff]: {diff}"

    raise AssertionError(msg)


def assert_numpy_array_equal(
    left,
    right,
    strict_nan=False,
    check_dtype=True,
    err_msg=None,
    check_same=None,
    obj="numpy array",
):
    """
    Check that 'np.ndarray' is equivalent.

    Parameters
    ----------
    left, right : numpy.ndarray or iterable
        The two arrays to be compared.
    strict_nan : bool, default False
        If True, consider NaN and None to be different.
    check_dtype : bool, default True
        Check dtype if both a and b are np.ndarray.
    err_msg : str, default None
        If provided, used as assertion message.
    check_same : None|'copy'|'same', default None
        Ensure left and right refer/do not refer to the same memory area.
    obj : str, default 'numpy array'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    # instance validation
    # Show a detailed error message when classes are different
    assert_class_equal(left, right, obj=obj)
    # both classes must be an np.ndarray
    _check_isinstance(left, right, np.ndarray)

    def _get_base(obj):
        return obj.base if getattr(obj, "base", None) is not None else obj

    left_base = _get_base(left)
    right_base = _get_base(right)

    if check_same == "same":
        if left_base is not right_base:
            raise AssertionError(f"{repr(left_base)} is not {repr(right_base)}")
    elif check_same == "copy":
        if left_base is right_base:
            raise AssertionError(f"{repr(left_base)} is {repr(right_base)}")

    def _raise(left, right, err_msg):
        if err_msg is None:
            if left.shape != right.shape:
                raise_assert_detail(
                    obj, f"{obj} shapes are different", left.shape, right.shape,
                )

            diff = 0
            for l, r in zip(left, right):
                # count up differences
                if not array_equivalent(l, r, strict_nan=strict_nan):
                    diff += 1

            diff = diff * 100.0 / left.size
            msg = f"{obj} values are different ({np.round(diff, 5)} %)"
            raise_assert_detail(obj, msg, left, right)

        raise AssertionError(err_msg)

    # compare shape and values
    if not array_equivalent(left, right, strict_nan=strict_nan):
        _raise(left, right, err_msg)

    if check_dtype:
        if isinstance(left, np.ndarray) and isinstance(right, np.ndarray):
            assert_attr_equal("dtype", left, right, obj=obj)


def assert_extension_array_equal(
    left, right, check_dtype=True, check_less_precise=False, check_exact=False
):
    """
    Check that left and right ExtensionArrays are equal.

    Parameters
    ----------
    left, right : ExtensionArray
        The two arrays to compare.
    check_dtype : bool, default True
        Whether to check if the ExtensionArray dtypes are identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_exact : bool, default False
        Whether to compare number exactly.

    Notes
    -----
    Missing values are checked separately from valid values.
    A mask of missing values is computed for each and checked to match.
    The remaining all-valid values are cast to object dtype and checked.
    """
    assert isinstance(left, ExtensionArray), "left is not an ExtensionArray"
    assert isinstance(right, ExtensionArray), "right is not an ExtensionArray"
    if check_dtype:
        assert_attr_equal("dtype", left, right, obj="ExtensionArray")

    if hasattr(left, "asi8") and type(right) == type(left):
        # Avoid slow object-dtype comparisons
        assert_numpy_array_equal(left.asi8, right.asi8)  # type: ignore
        return

    left_na = np.asarray(left.isna())
    right_na = np.asarray(right.isna())
    assert_numpy_array_equal(left_na, right_na, obj="ExtensionArray NA mask")

    left_valid = np.asarray(left[~left_na].astype(object))
    right_valid = np.asarray(right[~right_na].astype(object))
    if check_exact:
        assert_numpy_array_equal(left_valid, right_valid, obj="ExtensionArray")
    else:
        _testing.assert_almost_equal(
            left_valid,
            right_valid,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            obj="ExtensionArray",
        )


# This could be refactored to use the NDFrame.equals method
def assert_series_equal(
    left,
    right,
    check_dtype=True,
    check_index_type="equiv",
    check_series_type=True,
    check_less_precise=False,
    check_names=True,
    check_exact=False,
    check_datetimelike_compat=False,
    check_categorical=True,
    check_category_order=True,
    obj="Series",
):
    """
    Check that left and right Series are equal.

    Parameters
    ----------
    left : Series
    right : Series
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    check_index_type : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_series_type : bool, default True
        Whether to check the Series class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    check_names : bool, default True
        Whether to check the Series and Index names attribute.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_datetimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    check_category_order : bool, default True
        Whether to compare category order of internal Categoricals

        .. versionadded:: 1.0.2
    obj : str, default 'Series'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    # instance validation
    _check_isinstance(left, right, Series)

    if check_series_type:
        # ToDo: There are some tests using rhs is sparse
        # lhs is dense. Should use assert_class_equal in future
        assert isinstance(left, type(right))
        # assert_class_equal(left, right, obj=obj)

    # length comparison
    if len(left) != len(right):
        msg1 = f"{len(left)}, {left.index}"
        msg2 = f"{len(right)}, {right.index}"
        raise_assert_detail(obj, "Series length are different", msg1, msg2)

    # index comparison
    assert_index_equal(
        left.index,
        right.index,
        exact=check_index_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.index",
    )

    if check_dtype:
        # We want to skip exact dtype checking when `check_categorical`
        # is False. We'll still raise if only one is a `Categorical`,
        # regardless of `check_categorical`
        if (
            is_categorical_dtype(left)
            and is_categorical_dtype(right)
            and not check_categorical
        ):
            pass
        else:
            assert_attr_equal("dtype", left, right, obj=f"Attributes of {obj}")

    if check_exact:
        assert_numpy_array_equal(
            left._internal_get_values(),
            right._internal_get_values(),
            check_dtype=check_dtype,
            obj=str(obj),
        )
    elif check_datetimelike_compat:
        # we want to check only if we have compat dtypes
        # e.g. integer and M|m are NOT compat, but we can simply check
        # the values in that case
        if needs_i8_conversion(left) or needs_i8_conversion(right):

            # datetimelike may have different objects (e.g. datetime.datetime
            # vs Timestamp) but will compare equal
            if not Index(left.values).equals(Index(right.values)):
                msg = (
                    f"[datetimelike_compat=True] {left.values} "
                    f"is not equal to {right.values}."
                )
                raise AssertionError(msg)
        else:
            assert_numpy_array_equal(
                left._internal_get_values(),
                right._internal_get_values(),
                check_dtype=check_dtype,
            )
    elif is_interval_dtype(left) or is_interval_dtype(right):
        assert_interval_array_equal(left.array, right.array)
    elif is_extension_array_dtype(left.dtype) and is_datetime64tz_dtype(left.dtype):
        # .values is an ndarray, but ._values is the ExtensionArray.
        # TODO: Use .array
        assert is_extension_array_dtype(right.dtype)
        assert_extension_array_equal(left._values, right._values)
    elif (
        is_extension_array_dtype(left)
        and not is_categorical_dtype(left)
        and is_extension_array_dtype(right)
        and not is_categorical_dtype(right)
    ):
        assert_extension_array_equal(left.array, right.array)
    else:
        _testing.assert_almost_equal(
            left._internal_get_values(),
            right._internal_get_values(),
            check_less_precise=check_less_precise,
            check_dtype=check_dtype,
            obj=str(obj),
        )

    # metadata comparison
    if check_names:
        assert_attr_equal("name", left, right, obj=obj)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(
                left.values,
                right.values,
                obj=f"{obj} category",
                check_category_order=check_category_order,
            )


# This could be refactored to use the NDFrame.equals method
def assert_frame_equal(
    left,
    right,
    check_dtype=True,
    check_index_type="equiv",
    check_column_type="equiv",
    check_frame_type=True,
    check_less_precise=False,
    check_names=True,
    by_blocks=False,
    check_exact=False,
    check_datetimelike_compat=False,
    check_categorical=True,
    check_like=False,
    obj="DataFrame",
):
    """
    Check that left and right DataFrame are equal.

    This function is intended to compare two DataFrames and output any
    differences. Is is mostly intended for use in unit tests.
    Additional parameters allow varying the strictness of the
    equality checks performed.

    Parameters
    ----------
    left : DataFrame
        First DataFrame to compare.
    right : DataFrame
        Second DataFrame to compare.
    check_dtype : bool, default True
        Whether to check the DataFrame dtype is identical.
    check_index_type : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_column_type : bool or {'equiv'}, default 'equiv'
        Whether to check the columns class, dtype and inferred_type
        are identical. Is passed as the ``exact`` argument of
        :func:`assert_index_equal`.
    check_frame_type : bool, default True
        Whether to check the DataFrame class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    check_names : bool, default True
        Whether to check that the `names` attribute for both the `index`
        and `column` attributes of the DataFrame is identical.
    by_blocks : bool, default False
        Specify how to compare internal data. If False, compare by columns.
        If True, compare by blocks.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_datetimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    check_like : bool, default False
        If True, ignore the order of index & columns.
        Note: index labels must match their respective rows
        (same as in columns) - same labels must be with the same data.
    obj : str, default 'DataFrame'
        Specify object name being compared, internally used to show appropriate
        assertion message.

    See Also
    --------
    assert_series_equal : Equivalent method for asserting Series equality.
    DataFrame.equals : Check DataFrame equality.

    Examples
    --------
    This example shows comparing two DataFrames that are equal
    but with columns of differing dtypes.

    >>> from pandas._testing import assert_frame_equal
    >>> df1 = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})
    >>> df2 = pd.DataFrame({'a': [1, 2], 'b': [3.0, 4.0]})

    df1 equals itself.

    >>> assert_frame_equal(df1, df1)

    df1 differs from df2 as column 'b' is of a different type.

    >>> assert_frame_equal(df1, df2)
    Traceback (most recent call last):
    ...
    AssertionError: Attributes of DataFrame.iloc[:, 1] (column name="b") are different

    Attribute "dtype" are different
    [left]:  int64
    [right]: float64

    Ignore differing dtypes in columns with check_dtype.

    >>> assert_frame_equal(df1, df2, check_dtype=False)
    """
    __tracebackhide__ = True

    # instance validation
    _check_isinstance(left, right, DataFrame)

    if check_frame_type:
        assert isinstance(left, type(right))
        # assert_class_equal(left, right, obj=obj)

    # shape comparison
    if left.shape != right.shape:
        raise_assert_detail(
            obj, f"{obj} shape mismatch", f"{repr(left.shape)}", f"{repr(right.shape)}",
        )

    if check_like:
        left, right = left.reindex_like(right), right

    # index comparison
    assert_index_equal(
        left.index,
        right.index,
        exact=check_index_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.index",
    )

    # column comparison
    assert_index_equal(
        left.columns,
        right.columns,
        exact=check_column_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.columns",
    )

    # compare by blocks
    if by_blocks:
        rblocks = right._to_dict_of_blocks()
        lblocks = left._to_dict_of_blocks()
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            assert_frame_equal(
                lblocks[dtype], rblocks[dtype], check_dtype=check_dtype, obj=obj
            )

    # compare by columns
    else:
        for i, col in enumerate(left.columns):
            assert col in right
            lcol = left.iloc[:, i]
            rcol = right.iloc[:, i]
            assert_series_equal(
                lcol,
                rcol,
                check_dtype=check_dtype,
                check_index_type=check_index_type,
                check_less_precise=check_less_precise,
                check_exact=check_exact,
                check_names=check_names,
                check_datetimelike_compat=check_datetimelike_compat,
                check_categorical=check_categorical,
                obj=f'{obj}.iloc[:, {i}] (column name="{col}")',
            )


def assert_equal(left, right, **kwargs):
    """
    Wrapper for tm.assert_*_equal to dispatch to the appropriate test function.

    Parameters
    ----------
    left, right : Index, Series, DataFrame, ExtensionArray, or np.ndarray
        The two items to be compared.
    **kwargs
        All keyword arguments are passed through to the underlying assert method.
    """
    __tracebackhide__ = True

    if isinstance(left, pd.Index):
        assert_index_equal(left, right, **kwargs)
    elif isinstance(left, pd.Series):
        assert_series_equal(left, right, **kwargs)
    elif isinstance(left, pd.DataFrame):
        assert_frame_equal(left, right, **kwargs)
    elif isinstance(left, IntervalArray):
        assert_interval_array_equal(left, right, **kwargs)
    elif isinstance(left, PeriodArray):
        assert_period_array_equal(left, right, **kwargs)
    elif isinstance(left, DatetimeArray):
        assert_datetime_array_equal(left, right, **kwargs)
    elif isinstance(left, TimedeltaArray):
        assert_timedelta_array_equal(left, right, **kwargs)
    elif isinstance(left, ExtensionArray):
        assert_extension_array_equal(left, right, **kwargs)
    elif isinstance(left, np.ndarray):
        assert_numpy_array_equal(left, right, **kwargs)
    elif isinstance(left, str):
        assert kwargs == {}
        assert left == right
    else:
        raise NotImplementedError(type(left))
