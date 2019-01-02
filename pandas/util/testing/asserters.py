from contextlib import contextmanager
import re
import warnings

import numpy as np

from pandas._libs import testing as _testing
import pandas.compat as compat
from pandas.compat import PY2, raise_with_traceback, range, string_types, zip

from pandas.core.dtypes.common import (
    is_bool, is_categorical_dtype, is_datetime64tz_dtype,
    is_datetimelike_v_numeric, is_datetimelike_v_object,
    is_extension_array_dtype, is_interval_dtype, is_list_like, is_number,
    needs_i8_conversion)
from pandas.core.dtypes.missing import array_equivalent

import pandas as pd
from pandas import Categorical, DataFrame, Index, Series
from pandas.core.algorithms import take_1d
from pandas.core.arrays import (
    DatetimeArrayMixin as DatetimeArray, ExtensionArray, IntervalArray,
    PeriodArray, TimedeltaArrayMixin as TimedeltaArray)

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

    err_msg = "{name} Expected type {exp_type}, found {act_type} instead"
    cls_name = cls.__name__

    if not isinstance(left, cls):
        raise AssertionError(err_msg.format(name=cls_name, exp_type=cls,
                                            act_type=type(left)))
    if not isinstance(right, cls):
        raise AssertionError(err_msg.format(name=cls_name, exp_type=cls,
                                            act_type=type(right)))


def raise_assert_detail(obj, message, left, right, diff=None):
    __tracebackhide__ = True

    if isinstance(left, np.ndarray):
        left = pprint_thing(left)
    elif is_categorical_dtype(left):
        left = repr(left)

    if PY2 and isinstance(left, string_types):
        # left needs to be printable in native text type in python2
        left = left.encode('utf-8')

    if isinstance(right, np.ndarray):
        right = pprint_thing(right)
    elif is_categorical_dtype(right):
        right = repr(right)

    if PY2 and isinstance(right, string_types):
        # right needs to be printable in native text type in python2
        right = right.encode('utf-8')

    msg = """{obj} are different

{message}
[left]:  {left}
[right]: {right}""".format(obj=obj, message=message, left=left, right=right)

    if diff is not None:
        msg += "\n[diff]: {diff}".format(diff=diff)

    raise AssertionError(msg)


def assert_almost_equal(left, right, check_dtype="equiv",
                        check_less_precise=False, **kwargs):
    """
    Check that the left and right objects are approximately equal.

    By approximately equal, we refer to objects that are numbers or that
    contain numbers which may be equivalent to specific levels of precision.

    Parameters
    ----------
    left : object
    right : object
    check_dtype : bool / string {'equiv'}, default 'equiv'
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
        return assert_index_equal(left, right,
                                  check_exact=False,
                                  exact=check_dtype,
                                  check_less_precise=check_less_precise,
                                  **kwargs)

    elif isinstance(left, pd.Series):
        return assert_series_equal(left, right,
                                   check_exact=False,
                                   check_dtype=check_dtype,
                                   check_less_precise=check_less_precise,
                                   **kwargs)

    elif isinstance(left, pd.DataFrame):
        return assert_frame_equal(left, right,
                                  check_exact=False,
                                  check_dtype=check_dtype,
                                  check_less_precise=check_less_precise,
                                  **kwargs)

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
                if (isinstance(left, np.ndarray) or
                        isinstance(right, np.ndarray)):
                    obj = "numpy array"
                else:
                    obj = "Input"
                assert_class_equal(left, right, obj=obj)
        return _testing.assert_almost_equal(
            left, right,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs)


def assert_dict_equal(left, right, compare_keys=True):
    _check_isinstance(left, right, dict)
    return _testing.assert_dict_equal(left, right, compare_keys=compare_keys)


def assert_index_equal(left, right, exact='equiv', check_names=True,
                       check_less_precise=False, check_exact=True,
                       check_categorical=True, obj='Index'):
    """Check that left and right Index are equal.

    Parameters
    ----------
    left : Index
    right : Index
    exact : bool / string {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    check_names : bool, default True
        Whether to check the names attribute.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare
    check_exact : bool, default True
        Whether to compare number exactly.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    obj : str, default 'Index'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    __tracebackhide__ = True

    def _check_types(l, r, obj='Index'):
        if exact:
            assert_class_equal(l, r, exact=exact, obj=obj)

            # Skip exact dtype checking when `check_categorical` is False
            if check_categorical:
                assert_attr_equal('dtype', l, r, obj=obj)

            # allow string-like to have different inferred_types
            if l.inferred_type in ('string', 'unicode'):
                assert r.inferred_type in ('string', 'unicode')
            else:
                assert_attr_equal('inferred_type', l, r, obj=obj)

    def _get_ilevel_values(index, level):
        # accept level number only
        unique = index.levels[level]
        labels = index.codes[level]
        filled = take_1d(unique.values, labels, fill_value=unique._na_value)
        values = unique._shallow_copy(filled, name=index.names[level])
        return values

    # instance validation
    _check_isinstance(left, right, Index)

    # class / dtype comparison
    _check_types(left, right, obj=obj)

    # level comparison
    if left.nlevels != right.nlevels:
        msg1 = '{obj} levels are different'.format(obj=obj)
        msg2 = '{nlevels}, {left}'.format(nlevels=left.nlevels, left=left)
        msg3 = '{nlevels}, {right}'.format(nlevels=right.nlevels, right=right)
        raise_assert_detail(obj, msg1, msg2, msg3)

    # length comparison
    if len(left) != len(right):
        msg1 = '{obj} length are different'.format(obj=obj)
        msg2 = '{length}, {left}'.format(length=len(left), left=left)
        msg3 = '{length}, {right}'.format(length=len(right), right=right)
        raise_assert_detail(obj, msg1, msg2, msg3)

    # MultiIndex special comparison for little-friendly error messages
    if left.nlevels > 1:
        for level in range(left.nlevels):
            # cannot use get_level_values here because it can change dtype
            llevel = _get_ilevel_values(left, level)
            rlevel = _get_ilevel_values(right, level)

            lobj = 'MultiIndex level [{level}]'.format(level=level)
            assert_index_equal(llevel, rlevel,
                               exact=exact, check_names=check_names,
                               check_less_precise=check_less_precise,
                               check_exact=check_exact, obj=lobj)
            # get_level_values may change dtype
            _check_types(left.levels[level], right.levels[level], obj=obj)

    # skip exact index checking when `check_categorical` is False
    if check_exact and check_categorical:
        if not left.equals(right):
            diff = np.sum((left.values != right.values)
                          .astype(int)) * 100.0 / len(left)
            msg = '{obj} values are different ({pct} %)'.format(
                obj=obj, pct=np.round(diff, 5))
            raise_assert_detail(obj, msg, left, right)
    else:
        _testing.assert_almost_equal(left.values, right.values,
                                     check_less_precise=check_less_precise,
                                     check_dtype=exact,
                                     obj=obj, lobj=left, robj=right)

    # metadata comparison
    if check_names:
        assert_attr_equal('names', left, right, obj=obj)
    if isinstance(left, pd.PeriodIndex) or isinstance(right, pd.PeriodIndex):
        assert_attr_equal('freq', left, right, obj=obj)
    if (isinstance(left, pd.IntervalIndex) or
            isinstance(right, pd.IntervalIndex)):
        assert_interval_array_equal(left.values, right.values)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(left.values, right.values,
                                     obj='{obj} category'.format(obj=obj))


def assert_class_equal(left, right, exact=True, obj='Input'):
    """checks classes are equal."""
    __tracebackhide__ = True

    def repr_class(x):
        if isinstance(x, Index):
            # return Index as it is to include values in the error message
            return x

        try:
            return x.__class__.__name__
        except AttributeError:
            return repr(type(x))

    if exact == 'equiv':
        if type(left) != type(right):
            # allow equivalence of Int64Index/RangeIndex
            types = {type(left).__name__, type(right).__name__}
            if len(types - {'Int64Index', 'RangeIndex'}):
                msg = '{obj} classes are not equivalent'.format(obj=obj)
                raise_assert_detail(obj, msg, repr_class(left),
                                    repr_class(right))
    elif exact:
        if type(left) != type(right):
            msg = '{obj} classes are different'.format(obj=obj)
            raise_assert_detail(obj, msg, repr_class(left),
                                repr_class(right))


def assert_attr_equal(attr, left, right, obj='Attributes'):
    """checks attributes are equal. Both objects must have attribute.

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
    elif (is_number(left_attr) and np.isnan(left_attr) and
          is_number(right_attr) and np.isnan(right_attr)):
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
        msg = 'Attribute "{attr}" are different'.format(attr=attr)
        raise_assert_detail(obj, msg, left_attr, right_attr)


def assert_is_valid_plot_return_object(objs):
    import matplotlib.pyplot as plt
    if isinstance(objs, (pd.Series, np.ndarray)):
        for el in objs.ravel():
            msg = ("one of 'objs' is not a matplotlib Axes instance, type "
                   "encountered {name!r}").format(name=el.__class__.__name__)
            assert isinstance(el, (plt.Axes, dict)), msg
    else:
        assert isinstance(objs, (plt.Artist, tuple, dict)), (
            'objs is neither an ndarray of Artist instances nor a '
            'single Artist instance, tuple, or dict, "objs" is a {name!r}'
            .format(name=objs.__class__.__name__))


def assert_categorical_equal(left, right, check_dtype=True,
                             check_category_order=True, obj='Categorical'):
    """Test that Categoricals are equivalent.

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
        assert_index_equal(left.categories, right.categories,
                           obj='{obj}.categories'.format(obj=obj))
        assert_numpy_array_equal(left.codes, right.codes,
                                 check_dtype=check_dtype,
                                 obj='{obj}.codes'.format(obj=obj))
    else:
        assert_index_equal(left.categories.sort_values(),
                           right.categories.sort_values(),
                           obj='{obj}.categories'.format(obj=obj))
        assert_index_equal(left.categories.take(left.codes),
                           right.categories.take(right.codes),
                           obj='{obj}.values'.format(obj=obj))

    assert_attr_equal('ordered', left, right, obj=obj)


def assert_interval_array_equal(left, right, exact='equiv',
                                obj='IntervalArray'):
    """Test that two IntervalArrays are equivalent.

    Parameters
    ----------
    left, right : IntervalArray
        The IntervalArrays to compare.
    exact : bool / string {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    obj : str, default 'IntervalArray'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    _check_isinstance(left, right, IntervalArray)

    assert_index_equal(left.left, right.left, exact=exact,
                       obj='{obj}.left'.format(obj=obj))
    assert_index_equal(left.right, right.right, exact=exact,
                       obj='{obj}.left'.format(obj=obj))
    assert_attr_equal('closed', left, right, obj=obj)


def assert_period_array_equal(left, right, obj='PeriodArray'):
    _check_isinstance(left, right, PeriodArray)

    assert_numpy_array_equal(left._data, right._data,
                             obj='{obj}.values'.format(obj=obj))
    assert_attr_equal('freq', left, right, obj=obj)


def assert_datetime_array_equal(left, right, obj='DatetimeArray'):
    __tracebackhide__ = True
    _check_isinstance(left, right, DatetimeArray)

    assert_numpy_array_equal(left._data, right._data,
                             obj='{obj}._data'.format(obj=obj))
    assert_attr_equal('freq', left, right, obj=obj)
    assert_attr_equal('tz', left, right, obj=obj)


def assert_timedelta_array_equal(left, right, obj='TimedeltaArray'):
    __tracebackhide__ = True
    _check_isinstance(left, right, TimedeltaArray)
    assert_numpy_array_equal(left._data, right._data,
                             obj='{obj}._data'.format(obj=obj))
    assert_attr_equal('freq', left, right, obj=obj)


def assert_numpy_array_equal(left, right, strict_nan=False,
                             check_dtype=True, err_msg=None,
                             check_same=None, obj='numpy array'):
    """ Checks that 'np.ndarray' is equivalent

    Parameters
    ----------
    left : np.ndarray or iterable
    right : np.ndarray or iterable
    strict_nan : bool, default False
        If True, consider NaN and None to be different.
    check_dtype: bool, default True
        check dtype if both a and b are np.ndarray
    err_msg : str, default None
        If provided, used as assertion message
    check_same : None|'copy'|'same', default None
        Ensure left and right refer/do not refer to the same memory area
    obj : str, default 'numpy array'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    __tracebackhide__ = True

    # instance validation
    # Show a detailed error message when classes are different
    assert_class_equal(left, right, obj=obj)
    # both classes must be an np.ndarray
    _check_isinstance(left, right, np.ndarray)

    def _get_base(obj):
        return obj.base if getattr(obj, 'base', None) is not None else obj

    left_base = _get_base(left)
    right_base = _get_base(right)

    if check_same == 'same':
        if left_base is not right_base:
            msg = "{left!r} is not {right!r}".format(
                left=left_base, right=right_base)
            raise AssertionError(msg)
    elif check_same == 'copy':
        if left_base is right_base:
            msg = "{left!r} is {right!r}".format(
                left=left_base, right=right_base)
            raise AssertionError(msg)

    def _raise(left, right, err_msg):
        if err_msg is None:
            if left.shape != right.shape:
                raise_assert_detail(obj, '{obj} shapes are different'
                                    .format(obj=obj), left.shape, right.shape)

            diff = 0
            for l, r in zip(left, right):
                # count up differences
                if not array_equivalent(l, r, strict_nan=strict_nan):
                    diff += 1

            diff = diff * 100.0 / left.size
            msg = '{obj} values are different ({pct} %)'.format(
                obj=obj, pct=np.round(diff, 5))
            raise_assert_detail(obj, msg, left, right)

        raise AssertionError(err_msg)

    # compare shape and values
    if not array_equivalent(left, right, strict_nan=strict_nan):
        _raise(left, right, err_msg)

    if check_dtype:
        if isinstance(left, np.ndarray) and isinstance(right, np.ndarray):
            assert_attr_equal('dtype', left, right, obj=obj)

    return True


def assert_extension_array_equal(left, right, check_dtype=True,
                                 check_less_precise=False,
                                 check_exact=False):
    """Check that left and right ExtensionArrays are equal.

    Parameters
    ----------
    left, right : ExtensionArray
        The two arrays to compare
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
    assert isinstance(left, ExtensionArray), 'left is not an ExtensionArray'
    assert isinstance(right, ExtensionArray), 'right is not an ExtensionArray'
    if check_dtype:
        assert_attr_equal('dtype', left, right, obj='ExtensionArray')

    left_na = np.asarray(left.isna())
    right_na = np.asarray(right.isna())
    assert_numpy_array_equal(left_na, right_na, obj='ExtensionArray NA mask')

    left_valid = np.asarray(left[~left_na].astype(object))
    right_valid = np.asarray(right[~right_na].astype(object))
    if check_exact:
        assert_numpy_array_equal(left_valid, right_valid, obj='ExtensionArray')
    else:
        _testing.assert_almost_equal(left_valid, right_valid,
                                     check_dtype=check_dtype,
                                     check_less_precise=check_less_precise,
                                     obj='ExtensionArray')


# This could be refactored to use the NDFrame.equals method
def assert_series_equal(left, right, check_dtype=True,
                        check_index_type='equiv',
                        check_series_type=True,
                        check_less_precise=False,
                        check_names=True,
                        check_exact=False,
                        check_datetimelike_compat=False,
                        check_categorical=True,
                        obj='Series'):
    """Check that left and right Series are equal.

    Parameters
    ----------
    left : Series
    right : Series
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    check_index_type : bool / string {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_series_type : bool, default True
        Whether to check the Series class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_names : bool, default True
        Whether to check the Series and Index names attribute.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_datetimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
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
        msg1 = '{len}, {left}'.format(len=len(left), left=left.index)
        msg2 = '{len}, {right}'.format(len=len(right), right=right.index)
        raise_assert_detail(obj, 'Series length are different', msg1, msg2)

    # index comparison
    assert_index_equal(left.index, right.index, exact=check_index_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise,
                       check_exact=check_exact,
                       check_categorical=check_categorical,
                       obj='{obj}.index'.format(obj=obj))

    if check_dtype:
        # We want to skip exact dtype checking when `check_categorical`
        # is False. We'll still raise if only one is a `Categorical`,
        # regardless of `check_categorical`
        if (is_categorical_dtype(left) and is_categorical_dtype(right) and
                not check_categorical):
            pass
        else:
            assert_attr_equal('dtype', left, right)

    if check_exact:
        assert_numpy_array_equal(left.get_values(), right.get_values(),
                                 check_dtype=check_dtype,
                                 obj='{obj}'.format(obj=obj),)
    elif check_datetimelike_compat:
        # we want to check only if we have compat dtypes
        # e.g. integer and M|m are NOT compat, but we can simply check
        # the values in that case
        if (is_datetimelike_v_numeric(left, right) or
            is_datetimelike_v_object(left, right) or
            needs_i8_conversion(left) or
                needs_i8_conversion(right)):

            # datetimelike may have different objects (e.g. datetime.datetime
            # vs Timestamp) but will compare equal
            if not Index(left.values).equals(Index(right.values)):
                msg = ('[datetimelike_compat=True] {left} is not equal to '
                       '{right}.').format(left=left.values, right=right.values)
                raise AssertionError(msg)
        else:
            assert_numpy_array_equal(left.get_values(), right.get_values(),
                                     check_dtype=check_dtype)
    elif is_interval_dtype(left) or is_interval_dtype(right):
        assert_interval_array_equal(left.array, right.array)

    elif (is_extension_array_dtype(left.dtype) and
          is_datetime64tz_dtype(left.dtype)):
        # .values is an ndarray, but ._values is the ExtensionArray.
        # TODO: Use .array
        assert is_extension_array_dtype(right.dtype)
        return assert_extension_array_equal(left._values, right._values)

    elif (is_extension_array_dtype(left) and not is_categorical_dtype(left) and
          is_extension_array_dtype(right) and not is_categorical_dtype(right)):
        return assert_extension_array_equal(left.array, right.array)

    else:
        _testing.assert_almost_equal(left.get_values(), right.get_values(),
                                     check_less_precise=check_less_precise,
                                     check_dtype=check_dtype,
                                     obj='{obj}'.format(obj=obj))

    # metadata comparison
    if check_names:
        assert_attr_equal('name', left, right, obj=obj)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(left.values, right.values,
                                     obj='{obj} category'.format(obj=obj))


# This could be refactored to use the NDFrame.equals method
def assert_frame_equal(left, right, check_dtype=True,
                       check_index_type='equiv',
                       check_column_type='equiv',
                       check_frame_type=True,
                       check_less_precise=False,
                       check_names=True,
                       by_blocks=False,
                       check_exact=False,
                       check_datetimelike_compat=False,
                       check_categorical=True,
                       check_like=False,
                       obj='DataFrame'):
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
    check_index_type : bool / string {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_column_type : bool / string {'equiv'}, default 'equiv'
        Whether to check the columns class, dtype and inferred_type
        are identical. Is passed as the ``exact`` argument of
        :func:`assert_index_equal`.
    check_frame_type : bool, default True
        Whether to check the DataFrame class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_names : bool, default True
        Whether to check that the `names` attribute for both the `index`
        and `column` attributes of the DataFrame is identical, i.e.

        * left.index.names == right.index.names
        * left.columns.names == right.columns.names
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

    >>> from pandas.util.testing import assert_frame_equal
    >>> df1 = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})
    >>> df2 = pd.DataFrame({'a': [1, 2], 'b': [3.0, 4.0]})

    df1 equals itself.
    >>> assert_frame_equal(df1, df1)

    df1 differs from df2 as column 'b' is of a different type.
    >>> assert_frame_equal(df1, df2)
    Traceback (most recent call last):
    AssertionError: Attributes are different

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
        # ToDo: There are some tests using rhs is SparseDataFrame
        # lhs is DataFrame. Should use assert_class_equal in future
        assert isinstance(left, type(right))
        # assert_class_equal(left, right, obj=obj)

    # shape comparison
    if left.shape != right.shape:
        raise_assert_detail(obj,
                            'DataFrame shape mismatch',
                            '{shape!r}'.format(shape=left.shape),
                            '{shape!r}'.format(shape=right.shape))

    if check_like:
        left, right = left.reindex_like(right), right

    # index comparison
    assert_index_equal(left.index, right.index, exact=check_index_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise,
                       check_exact=check_exact,
                       check_categorical=check_categorical,
                       obj='{obj}.index'.format(obj=obj))

    # column comparison
    assert_index_equal(left.columns, right.columns, exact=check_column_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise,
                       check_exact=check_exact,
                       check_categorical=check_categorical,
                       obj='{obj}.columns'.format(obj=obj))

    # compare by blocks
    if by_blocks:
        rblocks = right._to_dict_of_blocks()
        lblocks = left._to_dict_of_blocks()
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            assert_frame_equal(lblocks[dtype], rblocks[dtype],
                               check_dtype=check_dtype, obj='DataFrame.blocks')

    # compare by columns
    else:
        for i, col in enumerate(left.columns):
            assert col in right
            lcol = left.iloc[:, i]
            rcol = right.iloc[:, i]
            assert_series_equal(
                lcol, rcol, check_dtype=check_dtype,
                check_index_type=check_index_type,
                check_less_precise=check_less_precise,
                check_exact=check_exact, check_names=check_names,
                check_datetimelike_compat=check_datetimelike_compat,
                check_categorical=check_categorical,
                obj='DataFrame.iloc[:, {idx}]'.format(idx=i))


def assert_panel_equal(left, right,
                       check_dtype=True,
                       check_panel_type=False,
                       check_less_precise=False,
                       check_names=False,
                       by_blocks=False,
                       obj='Panel'):
    """Check that left and right Panels are equal.

    Parameters
    ----------
    left : Panel (or nd)
    right : Panel (or nd)
    check_dtype : bool, default True
        Whether to check the Panel dtype is identical.
    check_panel_type : bool, default False
        Whether to check the Panel class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare
    check_names : bool, default True
        Whether to check the Index names attribute.
    by_blocks : bool, default False
        Specify how to compare internal data. If False, compare by columns.
        If True, compare by blocks.
    obj : str, default 'Panel'
        Specify the object name being compared, internally used to show
        the appropriate assertion message.
    """

    if check_panel_type:
        assert_class_equal(left, right, obj=obj)

    for axis in left._AXIS_ORDERS:
        left_ind = getattr(left, axis)
        right_ind = getattr(right, axis)
        assert_index_equal(left_ind, right_ind, check_names=check_names)

    if by_blocks:
        rblocks = right._to_dict_of_blocks()
        lblocks = left._to_dict_of_blocks()
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            array_equivalent(lblocks[dtype].values, rblocks[dtype].values)
    else:

        # can potentially be slow
        for i, item in enumerate(left._get_axis(0)):
            msg = "non-matching item (right) '{item}'".format(item=item)
            assert item in right, msg
            litem = left.iloc[i]
            ritem = right.iloc[i]
            assert_frame_equal(litem, ritem,
                               check_less_precise=check_less_precise,
                               check_names=check_names)

        for i, item in enumerate(right._get_axis(0)):
            msg = "non-matching item (left) '{item}'".format(item=item)
            assert item in left, msg


def assert_equal(left, right, **kwargs):
    """
    Wrapper for tm.assert_*_equal to dispatch to the appropriate test function.

    Parameters
    ----------
    left : Index, Series, DataFrame, ExtensionArray, or np.ndarray
    right : Index, Series, DataFrame, ExtensionArray, or np.ndarray
    **kwargs
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
    else:
        raise NotImplementedError(type(left))


# -----------------------------------------------------------------------------
# Sparse


def assert_sp_array_equal(left, right, check_dtype=True, check_kind=True,
                          check_fill_value=True,
                          consolidate_block_indices=False):
    """Check that the left and right SparseArray are equal.

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

    _check_isinstance(left, right, pd.SparseArray)

    assert_numpy_array_equal(left.sp_values, right.sp_values,
                             check_dtype=check_dtype)

    # SparseIndex comparison
    assert isinstance(left.sp_index, pd._libs.sparse.SparseIndex)
    assert isinstance(right.sp_index, pd._libs.sparse.SparseIndex)

    if not check_kind:
        left_index = left.sp_index.to_block_index()
        right_index = right.sp_index.to_block_index()
    else:
        left_index = left.sp_index
        right_index = right.sp_index

    if consolidate_block_indices and left.kind == 'block':
        # we'll probably remove this hack...
        left_index = left_index.to_int_index().to_block_index()
        right_index = right_index.to_int_index().to_block_index()

    if not left_index.equals(right_index):
        raise_assert_detail('SparseArray.index', 'index are not equal',
                            left_index, right_index)
    else:
        # Just ensure a
        pass

    if check_fill_value:
        assert_attr_equal('fill_value', left, right)
    if check_dtype:
        assert_attr_equal('dtype', left, right)
    assert_numpy_array_equal(left.values, right.values,
                             check_dtype=check_dtype)


def assert_sp_series_equal(left, right, check_dtype=True, exact_indices=True,
                           check_series_type=True, check_names=True,
                           check_kind=True,
                           check_fill_value=True,
                           consolidate_block_indices=False,
                           obj='SparseSeries'):
    """Check that the left and right SparseSeries are equal.

    Parameters
    ----------
    left : SparseSeries
    right : SparseSeries
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    exact_indices : bool, default True
    check_series_type : bool, default True
        Whether to check the SparseSeries class is identical.
    check_names : bool, default True
        Whether to check the SparseSeries name attribute.
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
    obj : str, default 'SparseSeries'
        Specify the object name being compared, internally used to show
        the appropriate assertion message.
    """
    _check_isinstance(left, right, pd.SparseSeries)

    if check_series_type:
        assert_class_equal(left, right, obj=obj)

    assert_index_equal(left.index, right.index,
                       obj='{obj}.index'.format(obj=obj))

    assert_sp_array_equal(left.values, right.values,
                          check_kind=check_kind,
                          check_fill_value=check_fill_value,
                          consolidate_block_indices=consolidate_block_indices)

    if check_names:
        assert_attr_equal('name', left, right)
    if check_dtype:
        assert_attr_equal('dtype', left, right)

    assert_numpy_array_equal(np.asarray(left.values),
                             np.asarray(right.values))


def assert_sp_frame_equal(left, right, check_dtype=True, exact_indices=True,
                          check_frame_type=True, check_kind=True,
                          check_fill_value=True,
                          consolidate_block_indices=False,
                          obj='SparseDataFrame'):
    """Check that the left and right SparseDataFrame are equal.

    Parameters
    ----------
    left : SparseDataFrame
    right : SparseDataFrame
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    exact_indices : bool, default True
        SparseSeries SparseIndex objects must be exactly the same,
        otherwise just compare dense representations.
    check_frame_type : bool, default True
        Whether to check the SparseDataFrame class is identical.
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
    obj : str, default 'SparseDataFrame'
        Specify the object name being compared, internally used to show
        the appropriate assertion message.
    """
    _check_isinstance(left, right, pd.SparseDataFrame)

    if check_frame_type:
        assert_class_equal(left, right, obj=obj)

    assert_index_equal(left.index, right.index,
                       obj='{obj}.index'.format(obj=obj))
    assert_index_equal(left.columns, right.columns,
                       obj='{obj}.columns'.format(obj=obj))

    if check_fill_value:
        assert_attr_equal('default_fill_value', left, right, obj=obj)

    for col, series in compat.iteritems(left):
        assert (col in right)
        # trade-off?

        if exact_indices:
            assert_sp_series_equal(
                series, right[col],
                check_dtype=check_dtype,
                check_kind=check_kind,
                check_fill_value=check_fill_value,
                consolidate_block_indices=consolidate_block_indices
            )
        else:
            assert_series_equal(series.to_dense(), right[col].to_dense(),
                                check_dtype=check_dtype)

    # do I care?
    # assert(left.default_kind == right.default_kind)

    for col in right:
        assert col in left


def assert_raises_regex(_exception, _regexp, _callable=None,
                        *args, **kwargs):
    r"""
    Check that the specified Exception is raised and that the error message
    matches a given regular expression pattern. This may be a regular
    expression object or a string containing a regular expression suitable
    for use by `re.search()`. This is a port of the `assertRaisesRegexp`
    function from unittest in Python 2.7.

    .. deprecated:: 0.24.0
        Use `pytest.raises` instead.

    Examples
    --------
    >>> assert_raises_regex(ValueError, 'invalid literal for.*XYZ', int, 'XYZ')
    >>> import re
    >>> assert_raises_regex(ValueError, re.compile('literal'), int, 'XYZ')

    If an exception of a different type is raised, it bubbles up.

    >>> assert_raises_regex(TypeError, 'literal', int, 'XYZ')
    Traceback (most recent call last):
        ...
    ValueError: invalid literal for int() with base 10: 'XYZ'
    >>> dct = dict()
    >>> assert_raises_regex(KeyError, 'pear', dct.__getitem__, 'apple')
    Traceback (most recent call last):
        ...
    AssertionError: "pear" does not match "'apple'"

    You can also use this in a with statement.

    >>> with assert_raises_regex(TypeError, r'unsupported operand type\(s\)'):
    ...     1 + {}
    >>> with assert_raises_regex(TypeError, 'banana'):
    ...     'apple'[0] = 'b'
    Traceback (most recent call last):
        ...
    AssertionError: "banana" does not match "'str' object does not support \
item assignment"
    """
    warnings.warn(("assert_raises_regex has been deprecated and will "
                   "be removed in the next release. Please use "
                   "`pytest.raises` instead."), FutureWarning, stacklevel=2)

    manager = _AssertRaisesContextmanager(exception=_exception, regexp=_regexp)
    if _callable is not None:
        with manager:
            _callable(*args, **kwargs)
    else:
        return manager


class _AssertRaisesContextmanager(object):
    """
    Context manager behind `assert_raises_regex`.
    """

    def __init__(self, exception, regexp=None):
        """
        Initialize an _AssertRaisesContextManager instance.

        Parameters
        ----------
        exception : class
            The expected Exception class.
        regexp : str, default None
            The regex to compare against the Exception message.
        """

        self.exception = exception

        if regexp is not None and not hasattr(regexp, "search"):
            regexp = re.compile(regexp, re.DOTALL)

        self.regexp = regexp

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, trace_back):
        expected = self.exception

        if not exc_type:
            exp_name = getattr(expected, "__name__", str(expected))
            raise AssertionError("{name} not raised.".format(name=exp_name))

        return self.exception_matches(exc_type, exc_value, trace_back)

    def exception_matches(self, exc_type, exc_value, trace_back):
        """
        Check that the Exception raised matches the expected Exception
        and expected error message regular expression.

        Parameters
        ----------
        exc_type : class
            The type of Exception raised.
        exc_value : Exception
            The instance of `exc_type` raised.
        trace_back : stack trace object
            The traceback object associated with `exc_value`.

        Returns
        -------
        is_matched : bool
            Whether or not the Exception raised matches the expected
            Exception class and expected error message regular expression.

        Raises
        ------
        AssertionError : The error message provided does not match
                         the expected error message regular expression.
        """

        if issubclass(exc_type, self.exception):
            if self.regexp is not None:
                val = str(exc_value)

                if not self.regexp.search(val):
                    msg = '"{pat}" does not match "{val}"'.format(
                        pat=self.regexp.pattern, val=val)
                    e = AssertionError(msg)
                    raise_with_traceback(e, trace_back)

            return True
        else:
            # Failed, so allow Exception to bubble up.
            return False


@contextmanager
def assert_produces_warning(expected_warning=Warning, filter_level="always",
                            clear=None, check_stacklevel=True):
    """
    Context manager for running code expected to either raise a specific
    warning, or not raise any warnings. Verifies that the code raises the
    expected warning, and that it does not raise any other unexpected
    warnings. It is basically a wrapper around ``warnings.catch_warnings``.

    Parameters
    ----------
    expected_warning : {Warning, False, None}, default Warning
        The type of Exception raised. ``exception.Warning`` is the base
        class for all warnings. To check that no warning is returned,
        specify ``False`` or ``None``.
    filter_level : str, default "always"
        Specifies whether warnings are ignored, displayed, or turned
        into errors.
        Valid values are:

        * "error" - turns matching warnings into exceptions
        * "ignore" - discard the warning
        * "always" - always emit a warning
        * "default" - print the warning the first time it is generated
          from each location
        * "module" - print the warning the first time it is generated
          from each module
        * "once" - print the warning the first time it is generated

    clear : str, default None
        If not ``None`` then remove any previously raised warnings from
        the ``__warningsregistry__`` to ensure that no warning messages are
        suppressed by this context manager. If ``None`` is specified,
        the ``__warningsregistry__`` keeps track of which warnings have been
        shown, and does not show them again.
    check_stacklevel : bool, default True
        If True, displays the line that called the function containing
        the warning to show were the function is called. Otherwise, the
        line that implements the function is displayed.

    Examples
    --------
    >>> import warnings
    >>> with assert_produces_warning():
    ...     warnings.warn(UserWarning())
    ...
    >>> with assert_produces_warning(False):
    ...     warnings.warn(RuntimeWarning())
    ...
    Traceback (most recent call last):
        ...
    AssertionError: Caused unexpected warning(s): ['RuntimeWarning'].
    >>> with assert_produces_warning(UserWarning):
    ...     warnings.warn(RuntimeWarning())
    Traceback (most recent call last):
        ...
    AssertionError: Did not see expected warning of class 'UserWarning'.

    ..warn:: This is *not* thread-safe.
    """
    __tracebackhide__ = True

    with warnings.catch_warnings(record=True) as w:

        if clear is not None:
            # make sure that we are clearing these warnings
            # if they have happened before
            # to guarantee that we will catch them
            if not is_list_like(clear):
                clear = [clear]
            for m in clear:
                try:
                    m.__warningregistry__.clear()
                except Exception:
                    pass

        saw_warning = False
        warnings.simplefilter(filter_level)
        yield w
        extra_warnings = []

        for actual_warning in w:
            if (expected_warning and issubclass(actual_warning.category,
                                                expected_warning)):
                saw_warning = True

                if check_stacklevel and issubclass(actual_warning.category,
                                                   (FutureWarning,
                                                    DeprecationWarning)):
                    from inspect import getframeinfo, stack
                    caller = getframeinfo(stack()[2][0])
                    msg = ("Warning not set with correct stacklevel. "
                           "File where warning is raised: {actual} != "
                           "{caller}. Warning message: {message}"
                           ).format(actual=actual_warning.filename,
                                    caller=caller.filename,
                                    message=actual_warning.message)
                    assert actual_warning.filename == caller.filename, msg
            else:
                extra_warnings.append((actual_warning.category.__name__,
                                       actual_warning.message,
                                       actual_warning.filename,
                                       actual_warning.lineno))
        if expected_warning:
            msg = "Did not see expected warning of class {name!r}.".format(
                name=expected_warning.__name__)
            assert saw_warning, msg
        assert not extra_warnings, ("Caused unexpected warning(s): {extra!r}."
                                    ).format(extra=extra_warnings)
