import warnings

from pandas._testing import *  # noqa

warnings.warn(
    (
        "pandas.util.testing is deprecated. Use the functions in the "
        "public API at pandas.testing instead."
    ),
    FutureWarning,
    stacklevel=2,
)


# def assert_almost_equal(
#     left, right, check_dtype="equiv", check_less_precise=False, **kwargs
# ):
#     """
#     Check that the left and right objects are approximately equal.

#     By approximately equal, we refer to objects that are numbers or that
#     contain numbers which may be equivalent to specific levels of precision.

#     Parameters
#     ----------
#     left : object
#     right : object
#     check_dtype : bool or {'equiv'}, default 'equiv'
#         Check dtype if both a and b are the same type. If 'equiv' is passed in,
#         then `RangeIndex` and `Int64Index` are also considered equivalent
#         when doing type checking.
#     check_less_precise : bool or int, default False
#         Specify comparison precision. 5 digits (False) or 3 digits (True)
#         after decimal points are compared. If int, then specify the number
#         of digits to compare.

#         When comparing two numbers, if the first number has magnitude less
#         than 1e-5, we compare the two numbers directly and check whether
#         they are equivalent within the specified precision. Otherwise, we
#         compare the **ratio** of the second number to the first number and
#         check whether it is equivalent to 1 within the specified precision.
#     """

#     if isinstance(left, pd.Index):
#         assert_index_equal(
#             left,
#             right,
#             check_exact=False,
#             exact=check_dtype,
#             check_less_precise=check_less_precise,
#             **kwargs,
#         )

#     elif isinstance(left, pd.Series):
#         assert_series_equal(
#             left,
#             right,
#             check_exact=False,
#             check_dtype=check_dtype,
#             check_less_precise=check_less_precise,
#             **kwargs,
#         )

#     elif isinstance(left, pd.DataFrame):
#         assert_frame_equal(
#             left,
#             right,
#             check_exact=False,
#             check_dtype=check_dtype,
#             check_less_precise=check_less_precise,
#             **kwargs,
#         )

#     else:
#         # Other sequences.
#         if check_dtype:
#             if is_number(left) and is_number(right):
#                 # Do not compare numeric classes, like np.float64 and float.
#                 pass
#             elif is_bool(left) and is_bool(right):
#                 # Do not compare bool classes, like np.bool_ and bool.
#                 pass
#             else:
#                 if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
#                     obj = "numpy array"
#                 else:
#                     obj = "Input"
#                 assert_class_equal(left, right, obj=obj)
#         _testing.assert_almost_equal(
#             left,
#             right,
#             check_dtype=check_dtype,
#             check_less_precise=check_less_precise,
#             **kwargs,
#         )


# def assert_class_equal(left, right, exact=True, obj="Input"):
#     """checks classes are equal."""
#     __tracebackhide__ = True

#     def repr_class(x):
#         if isinstance(x, Index):
#             # return Index as it is to include values in the error message
#             return x

#         try:
#             return x.__class__.__name__
#         except AttributeError:
#             return repr(type(x))

#     if exact == "equiv":
#         if type(left) != type(right):
#             # allow equivalence of Int64Index/RangeIndex
#             types = {type(left).__name__, type(right).__name__}
#             if len(types - {"Int64Index", "RangeIndex"}):
#                 msg = "{obj} classes are not equivalent".format(obj=obj)
#                 raise_assert_detail(obj, msg, repr_class(left), repr_class(right))
#     elif exact:
#         if type(left) != type(right):
#             msg = "{obj} classes are different".format(obj=obj)
#             raise_assert_detail(obj, msg, repr_class(left), repr_class(right))


# def assert_attr_equal(attr, left, right, obj="Attributes"):
#     """checks attributes are equal. Both objects must have attribute.

#     Parameters
#     ----------
#     attr : str
#         Attribute name being compared.
#     left : object
#     right : object
#     obj : str, default 'Attributes'
#         Specify object name being compared, internally used to show appropriate
#         assertion message
#     """
#     __tracebackhide__ = True

#     left_attr = getattr(left, attr)
#     right_attr = getattr(right, attr)

#     if left_attr is right_attr:
#         return True
#     elif (
#         is_number(left_attr)
#         and np.isnan(left_attr)
#         and is_number(right_attr)
#         and np.isnan(right_attr)
#     ):
#         # np.nan
#         return True

#     try:
#         result = left_attr == right_attr
#     except TypeError:
#         # datetimetz on rhs may raise TypeError
#         result = False
#     if not isinstance(result, bool):
#         result = result.all()

#     if result:
#         return True
#     else:
#         msg = 'Attribute "{attr}" are different'.format(attr=attr)
#         raise_assert_detail(obj, msg, left_attr, right_attr)


# def assert_categorical_equal(
#     left: Categorical,
#     right: Categorical,
#     check_dtype: bool = True,
#     check_category_order: bool = True,
#     obj: str = "Categorical",
# ) -> None:
#     """Test that Categoricals are equivalent.

#     Parameters
#     ----------
#     left : Categorical
#     right : Categorical
#     check_dtype : bool, default True
#         Check that integer dtype of the codes are the same
#     check_category_order : bool, default True
#         Whether the order of the categories should be compared, which
#         implies identical integer codes.  If False, only the resulting
#         values are compared.  The ordered attribute is
#         checked regardless.
#     obj : str, default 'Categorical'
#         Specify object name being compared, internally used to show appropriate
#         assertion message
#     """
#     _check_isinstance(left, right, Categorical)

#     if check_category_order:
#         assert_index_equal(
#             left.categories, right.categories, obj="{obj}.categories".format(obj=obj)
#         )
#         assert_numpy_array_equal(
#             left.codes,
#             right.codes,
#             check_dtype=check_dtype,
#             obj="{obj}.codes".format(obj=obj),
#         )
#     else:
#         assert_index_equal(
#             left.categories.sort_values(),
#             right.categories.sort_values(),
#             obj="{obj}.categories".format(obj=obj),
#         )
#         assert_index_equal(
#             left.categories.take(left.codes),
#             right.categories.take(right.codes),
#             obj="{obj}.values".format(obj=obj),
#         )

#     assert_attr_equal("ordered", left, right, obj=obj)


# def assert_interval_array_equal(
#     left: IntervalArray,
#     right: IntervalArray,
#     exact: str = "equiv",
#     obj: str = "IntervalArray",
# ) -> None:
#     """Test that two IntervalArrays are equivalent.

#     Parameters
#     ----------
#     left, right : IntervalArray
#         The IntervalArrays to compare.
#     exact : bool or {'equiv'}, default 'equiv'
#         Whether to check the Index class, dtype and inferred_type
#         are identical. If 'equiv', then RangeIndex can be substituted for
#         Int64Index as well.
#     obj : str, default 'IntervalArray'
#         Specify object name being compared, internally used to show appropriate
#         assertion message
#     """
#     assert_index_equal(
#         left.left, right.left, exact=exact, obj="{obj}.left".format(obj=obj)
#     )
#     assert_index_equal(
#         left.right, right.right, exact=exact, obj="{obj}.left".format(obj=obj)
#     )
#     assert_attr_equal("closed", left, right, obj=obj)


# def assert_period_array_equal(
#     left: PeriodArray, right: PeriodArray, obj: str = "PeriodArray"
# ) -> None:
#     _check_isinstance(left, right, PeriodArray)

#     assert_numpy_array_equal(
#         left._data, right._data, obj="{obj}.values".format(obj=obj)
#     )
#     assert_attr_equal("freq", left, right, obj=obj)


# def assert_datetime_array_equal(
#     left: DatetimeArray, right: DatetimeArray, obj: str = "DatetimeArray"
# ) -> None:
#     __tracebackhide__ = True
#     _check_isinstance(left, right, DatetimeArray)

#     assert_numpy_array_equal(left._data, right._data, obj="{obj}._data".format(obj=obj))
#     assert_attr_equal("freq", left, right, obj=obj)
#     assert_attr_equal("tz", left, right, obj=obj)


# def assert_timedelta_array_equal(
#     left: TimedeltaArray, right: TimedeltaArray, obj: str = "TimedeltaArray"
# ) -> None:
#     __tracebackhide__ = True
#     _check_isinstance(left, right, TimedeltaArray)
#     assert_numpy_array_equal(left._data, right._data, obj="{obj}._data".format(obj=obj))
#     assert_attr_equal("freq", left, right, obj=obj)


# def assert_numpy_array_equal(
#     left: np.ndarray,
#     right: np.ndarray,
#     strict_nan: bool = False,
#     check_dtype: bool = True,
#     err_msg: Optional[str] = None,
#     check_same: Optional[str] = None,
#     obj: str = "numpy array",
# ) -> None:
#     """ Checks that 'np.ndarray' is equivalent

#     Parameters
#     ----------
#     left : np.ndarray or iterable
#     right : np.ndarray or iterable
#     strict_nan : bool, default False
#         If True, consider NaN and None to be different.
#     check_dtype: bool, default True
#         check dtype if both a and b are np.ndarray
#     err_msg : str, default None
#         If provided, used as assertion message
#     check_same : None|'copy'|'same', default None
#         Ensure left and right refer/do not refer to the same memory area
#     obj : str, default 'numpy array'
#         Specify object name being compared, internally used to show appropriate
#         assertion message
#     """
#     __tracebackhide__ = True

#     # instance validation
#     # Show a detailed error message when classes are different
#     assert_class_equal(left, right, obj=obj)
#     # both classes must be an np.ndarray
#     _check_isinstance(left, right, np.ndarray)

#     def _get_base(obj):
#         return obj.base if getattr(obj, "base", None) is not None else obj

#     left_base = _get_base(left)
#     right_base = _get_base(right)

#     if check_same == "same":
#         if left_base is not right_base:
#             msg = "{left!r} is not {right!r}".format(left=left_base, right=right_base)
#             raise AssertionError(msg)
#     elif check_same == "copy":
#         if left_base is right_base:
#             msg = "{left!r} is {right!r}".format(left=left_base, right=right_base)
#             raise AssertionError(msg)

#     def _raise(left, right, err_msg):
#         if err_msg is None:
#             if left.shape != right.shape:
#                 raise_assert_detail(
#                     obj,
#                     "{obj} shapes are different".format(obj=obj),
#                     left.shape,
#                     right.shape,
#                 )

#             diff = 0
#             for l, r in zip(left, right):
#                 # count up differences
#                 if not array_equivalent(l, r, strict_nan=strict_nan):
#                     diff += 1

#             diff = diff * 100.0 / left.size
#             msg = "{obj} values are different ({pct} %)".format(
#                 obj=obj, pct=np.round(diff, 5)
#             )
#             raise_assert_detail(obj, msg, left, right)

#         raise AssertionError(err_msg)

#     # compare shape and values
#     if not array_equivalent(left, right, strict_nan=strict_nan):
#         _raise(left, right, err_msg)

#     if check_dtype:
#         if isinstance(left, np.ndarray) and isinstance(right, np.ndarray):
#             assert_attr_equal("dtype", left, right, obj=obj)


# def assert_extension_array_equal(
#     left, right, check_dtype=True, check_less_precise=False, check_exact=False
# ):
#     """Check that left and right ExtensionArrays are equal.

#     Parameters
#     ----------
#     left, right : ExtensionArray
#         The two arrays to compare
#     check_dtype : bool, default True
#         Whether to check if the ExtensionArray dtypes are identical.
#     check_less_precise : bool or int, default False
#         Specify comparison precision. Only used when check_exact is False.
#         5 digits (False) or 3 digits (True) after decimal points are compared.
#         If int, then specify the digits to compare.
#     check_exact : bool, default False
#         Whether to compare number exactly.

#     Notes
#     -----
#     Missing values are checked separately from valid values.
#     A mask of missing values is computed for each and checked to match.
#     The remaining all-valid values are cast to object dtype and checked.
#     """
#     assert isinstance(left, ExtensionArray), "left is not an ExtensionArray"
#     assert isinstance(right, ExtensionArray), "right is not an ExtensionArray"
#     if check_dtype:
#         assert_attr_equal("dtype", left, right, obj="ExtensionArray")

#     if hasattr(left, "asi8") and type(right) == type(left):
#         # Avoid slow object-dtype comparisons
#         assert_numpy_array_equal(left.asi8, right.asi8)
#         return

#     left_na = np.asarray(left.isna())
#     right_na = np.asarray(right.isna())
#     assert_numpy_array_equal(left_na, right_na, obj="ExtensionArray NA mask")

#     left_valid = np.asarray(left[~left_na].astype(object))
#     right_valid = np.asarray(right[~right_na].astype(object))
#     if check_exact:
#         assert_numpy_array_equal(left_valid, right_valid, obj="ExtensionArray")
#     else:
#         _testing.assert_almost_equal(
#             left_valid,
#             right_valid,
#             check_dtype=check_dtype,
#             check_less_precise=check_less_precise,
#             obj="ExtensionArray",
#         )


# # This could be refactored to use the NDFrame.equals method
# def assert_series_equal(
#     left: Series,
#     right: Series,
#     check_dtype: bool = True,
#     check_index_type: str = "equiv",
#     check_series_type: bool = True,
#     check_less_precise: bool = False,
#     check_names: bool = True,
#     check_exact: bool = False,
#     check_datetimelike_compat: bool = False,
#     check_categorical: bool = True,
#     obj: str = "Series",
# ) -> None:
#     """
#     Check that left and right Series are equal.

#     Parameters
#     ----------
#     left : Series
#     right : Series
#     check_dtype : bool, default True
#         Whether to check the Series dtype is identical.
#     check_index_type : bool or {'equiv'}, default 'equiv'
#         Whether to check the Index class, dtype and inferred_type
#         are identical.
#     check_series_type : bool, default True
#         Whether to check the Series class is identical.
#     check_less_precise : bool or int, default False
#         Specify comparison precision. Only used when check_exact is False.
#         5 digits (False) or 3 digits (True) after decimal points are compared.
#         If int, then specify the digits to compare.

#         When comparing two numbers, if the first number has magnitude less
#         than 1e-5, we compare the two numbers directly and check whether
#         they are equivalent within the specified precision. Otherwise, we
#         compare the **ratio** of the second number to the first number and
#         check whether it is equivalent to 1 within the specified precision.
#     check_names : bool, default True
#         Whether to check the Series and Index names attribute.
#     check_exact : bool, default False
#         Whether to compare number exactly.
#     check_datetimelike_compat : bool, default False
#         Compare datetime-like which is comparable ignoring dtype.
#     check_categorical : bool, default True
#         Whether to compare internal Categorical exactly.
#     obj : str, default 'Series'
#         Specify object name being compared, internally used to show appropriate
#         assertion message.
#     """
#     __tracebackhide__ = True

#     # instance validation
#     _check_isinstance(left, right, Series)

#     if check_series_type:
#         # ToDo: There are some tests using rhs is sparse
#         # lhs is dense. Should use assert_class_equal in future
#         assert isinstance(left, type(right))
#         # assert_class_equal(left, right, obj=obj)

#     # length comparison
#     if len(left) != len(right):
#         msg1 = "{len}, {left}".format(len=len(left), left=left.index)
#         msg2 = "{len}, {right}".format(len=len(right), right=right.index)
#         raise_assert_detail(obj, "Series length are different", msg1, msg2)

#     # index comparison
#     assert_index_equal(
#         left.index,
#         right.index,
#         exact=check_index_type,
#         check_names=check_names,
#         check_less_precise=check_less_precise,
#         check_exact=check_exact,
#         check_categorical=check_categorical,
#         obj="{obj}.index".format(obj=obj),
#     )

#     if check_dtype:
#         # We want to skip exact dtype checking when `check_categorical`
#         # is False. We'll still raise if only one is a `Categorical`,
#         # regardless of `check_categorical`
#         if (
#             is_categorical_dtype(left)
#             and is_categorical_dtype(right)
#             and not check_categorical
#         ):
#             pass
#         else:
#             assert_attr_equal(
#                 "dtype", left, right, obj="Attributes of {obj}".format(obj=obj)
#             )

#     if check_exact:
#         assert_numpy_array_equal(
#             left._internal_get_values(),
#             right._internal_get_values(),
#             check_dtype=check_dtype,
#             obj="{obj}".format(obj=obj),
#         )
#     elif check_datetimelike_compat:
#         # we want to check only if we have compat dtypes
#         # e.g. integer and M|m are NOT compat, but we can simply check
#         # the values in that case
#         if needs_i8_conversion(left) or needs_i8_conversion(right):

#             # datetimelike may have different objects (e.g. datetime.datetime
#             # vs Timestamp) but will compare equal
#             if not Index(left.values).equals(Index(right.values)):
#                 msg = (
#                     "[datetimelike_compat=True] {left} is not equal to {right}."
#                 ).format(left=left.values, right=right.values)
#                 raise AssertionError(msg)
#         else:
#             assert_numpy_array_equal(
#                 left._internal_get_values(),
#                 right._internal_get_values(),
#                 check_dtype=check_dtype,
#             )
#     elif is_interval_dtype(left) or is_interval_dtype(right):
#         left_array = cast(IntervalArray, left.array)
#         right_array = cast(IntervalArray, right.array)
#         assert_interval_array_equal(left_array, right_array)
#     elif is_extension_array_dtype(left.dtype) and is_datetime64tz_dtype(left.dtype):
#         # .values is an ndarray, but ._values is the ExtensionArray.
#         # TODO: Use .array
#         assert is_extension_array_dtype(right.dtype)
#         assert_extension_array_equal(left._values, right._values)
#     elif (
#         is_extension_array_dtype(left)
#         and not is_categorical_dtype(left)
#         and is_extension_array_dtype(right)
#         and not is_categorical_dtype(right)
#     ):
#         assert_extension_array_equal(left.array, right.array)
#     else:
#         _testing.assert_almost_equal(
#             left._internal_get_values(),
#             right._internal_get_values(),
#             check_less_precise=check_less_precise,
#             check_dtype=check_dtype,
#             obj="{obj}".format(obj=obj),
#         )

#     # metadata comparison
#     if check_names:
#         assert_attr_equal("name", left, right, obj=obj)

#     if check_categorical:
#         if is_categorical_dtype(left) or is_categorical_dtype(right):
#             assert_categorical_equal(
#                 left.values, right.values, obj="{obj} category".format(obj=obj)
#             )


# # This could be refactored to use the NDFrame.equals method
# def assert_frame_equal(
#     left: Any,
#     right: Any,
#     check_dtype: bool = True,
#     check_index_type: str = "equiv",
#     check_column_type: str = "equiv",
#     check_frame_type: bool = True,
#     check_less_precise: bool = False,
#     check_names: bool = True,
#     by_blocks: bool = False,
#     check_exact: bool = False,
#     check_datetimelike_compat: bool = False,
#     check_categorical: bool = True,
#     check_like: bool = False,
#     obj: str = "DataFrame",
# ) -> None:
#     """
#     Check that left and right DataFrame are equal.

#     This function is intended to compare two DataFrames and output any
#     differences. Is is mostly intended for use in unit tests.
#     Additional parameters allow varying the strictness of the
#     equality checks performed.

#     Parameters
#     ----------
#     left : Any
#         First DataFrame to compare.
#     right : Any
#         Second DataFrame to compare.
#     check_dtype : bool, default True
#         Whether to check the DataFrame dtype is identical.
#     check_index_type : bool or {'equiv'}, default 'equiv'
#         Whether to check the Index class, dtype and inferred_type
#         are identical.
#     check_column_type : bool or {'equiv'}, default 'equiv'
#         Whether to check the columns class, dtype and inferred_type
#         are identical. Is passed as the ``exact`` argument of
#         :func:`assert_index_equal`.
#     check_frame_type : bool, default True
#         Whether to check the DataFrame class is identical.
#     check_less_precise : bool or int, default False
#         Specify comparison precision. Only used when check_exact is False.
#         5 digits (False) or 3 digits (True) after decimal points are compared.
#         If int, then specify the digits to compare.

#         When comparing two numbers, if the first number has magnitude less
#         than 1e-5, we compare the two numbers directly and check whether
#         they are equivalent within the specified precision. Otherwise, we
#         compare the **ratio** of the second number to the first number and
#         check whether it is equivalent to 1 within the specified precision.
#     check_names : bool, default True
#         Whether to check that the `names` attribute for both the `index`
#         and `column` attributes of the DataFrame is identical.
#     by_blocks : bool, default False
#         Specify how to compare internal data. If False, compare by columns.
#         If True, compare by blocks.
#     check_exact : bool, default False
#         Whether to compare number exactly.
#     check_datetimelike_compat : bool, default False
#         Compare datetime-like which is comparable ignoring dtype.
#     check_categorical : bool, default True
#         Whether to compare internal Categorical exactly.
#     check_like : bool, default False
#         If True, ignore the order of index & columns.
#         Note: index labels must match their respective rows
#         (same as in columns) - same labels must be with the same data.
#     obj : str, default 'DataFrame'
#         Specify object name being compared, internally used to show appropriate
#         assertion message.

#     See Also
#     --------
#     assert_series_equal : Equivalent method for asserting Series equality.
#     DataFrame.equals : Check DataFrame equality.

#     Examples
#     --------
#     This example shows comparing two DataFrames that are equal
#     but with columns of differing dtypes.

#     >>> from pandas.util.testing import assert_frame_equal
#     >>> df1 = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})
#     >>> df2 = pd.DataFrame({'a': [1, 2], 'b': [3.0, 4.0]})

#     df1 equals itself.

#     >>> assert_frame_equal(df1, df1)

#     df1 differs from df2 as column 'b' is of a different type.

#     >>> assert_frame_equal(df1, df2)
#     Traceback (most recent call last):
#     ...
#     AssertionError: Attributes of DataFrame.iloc[:, 1] are different

#     Attribute "dtype" are different
#     [left]:  int64
#     [right]: float64

#     Ignore differing dtypes in columns with check_dtype.

#     >>> assert_frame_equal(df1, df2, check_dtype=False)
#     """
#     __tracebackhide__ = True

#     # instance validation
#     _check_isinstance(left, right, DataFrame)

#     if check_frame_type:
#         assert isinstance(left, type(right))
#         # assert_class_equal(left, right, obj=obj)

#     # shape comparison
#     if left.shape != right.shape:
#         raise_assert_detail(
#             obj,
#             "{obj} shape mismatch".format(obj=obj),
#             "{shape!r}".format(shape=left.shape),
#             "{shape!r}".format(shape=right.shape),
#         )

#     if check_like:
#         left, right = left.reindex_like(right), right

#     # index comparison
#     assert_index_equal(
#         left.index,
#         right.index,
#         exact=check_index_type,
#         check_names=check_names,
#         check_less_precise=check_less_precise,
#         check_exact=check_exact,
#         check_categorical=check_categorical,
#         obj="{obj}.index".format(obj=obj),
#     )

#     # column comparison
#     assert_index_equal(
#         left.columns,
#         right.columns,
#         exact=check_column_type,
#         check_names=check_names,
#         check_less_precise=check_less_precise,
#         check_exact=check_exact,
#         check_categorical=check_categorical,
#         obj="{obj}.columns".format(obj=obj),
#     )

#     # compare by blocks
#     if by_blocks:
#         rblocks = right._to_dict_of_blocks()
#         lblocks = left._to_dict_of_blocks()
#         for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
#             assert dtype in lblocks
#             assert dtype in rblocks
#             assert_frame_equal(
#                 lblocks[dtype], rblocks[dtype], check_dtype=check_dtype, obj=obj
#             )

#     # compare by columns
#     else:
#         for i, col in enumerate(left.columns):
#             assert col in right
#             lcol = left.iloc[:, i]
#             rcol = right.iloc[:, i]
#             assert_series_equal(
#                 lcol,
#                 rcol,
#                 check_dtype=check_dtype,
#                 check_index_type=check_index_type,
#                 check_less_precise=check_less_precise,
#                 check_exact=check_exact,
#                 check_names=check_names,
#                 check_datetimelike_compat=check_datetimelike_compat,
#                 check_categorical=check_categorical,
#                 obj="{obj}.iloc[:, {idx}]".format(obj=obj, idx=i),
#             )


# def assert_equal(
#     left: Union[DataFrame, AnyArrayLike],
#     right: Union[DataFrame, AnyArrayLike],
#     **kwargs,
# ) -> None:
#     """
#     Wrapper for tm.assert_*_equal to dispatch to the appropriate test function.

#     Parameters
#     ----------
#     left : Index, Series, DataFrame, ExtensionArray, or np.ndarray
#     right : Index, Series, DataFrame, ExtensionArray, or np.ndarray
#     **kwargs
#     """
#     __tracebackhide__ = True

#     if isinstance(left, Index):
#         right = cast(Index, right)
#         assert_index_equal(left, right, **kwargs)
#     elif isinstance(left, Series):
#         right = cast(Series, right)
#         assert_series_equal(left, right, **kwargs)
#     elif isinstance(left, DataFrame):
#         right = cast(DataFrame, right)
#         assert_frame_equal(left, right, **kwargs)
#     elif isinstance(left, IntervalArray):
#         right = cast(IntervalArray, right)
#         assert_interval_array_equal(left, right, **kwargs)
#     elif isinstance(left, PeriodArray):
#         right = cast(PeriodArray, right)
#         assert_period_array_equal(left, right, **kwargs)
#     elif isinstance(left, DatetimeArray):
#         right = cast(DatetimeArray, right)
#         assert_datetime_array_equal(left, right, **kwargs)
#     elif isinstance(left, TimedeltaArray):
#         right = cast(TimedeltaArray, right)
#         assert_timedelta_array_equal(left, right, **kwargs)
#     elif isinstance(left, ExtensionArray):
#         right = cast(ExtensionArray, right)
#         assert_extension_array_equal(left, right, **kwargs)
#     elif isinstance(left, np.ndarray):
#         right = cast(np.ndarray, right)
#         assert_numpy_array_equal(left, right, **kwargs)
#     elif isinstance(left, str):
#         assert kwargs == {}
#         assert left == right
#     else:
#         raise NotImplementedError(type(left))


# def assert_sp_array_equal(
#     left: pd.SparseArray,
#     right: pd.SparseArray,
#     check_dtype: bool = True,
#     check_kind: bool = True,
#     check_fill_value: bool = True,
#     consolidate_block_indices: bool = False,
# ):
#     """Check that the left and right SparseArray are equal.

#     Parameters
#     ----------
#     left : SparseArray
#     right : SparseArray
#     check_dtype : bool, default True
#         Whether to check the data dtype is identical.
#     check_kind : bool, default True
#         Whether to just the kind of the sparse index for each column.
#     check_fill_value : bool, default True
#         Whether to check that left.fill_value matches right.fill_value
#     consolidate_block_indices : bool, default False
#         Whether to consolidate contiguous blocks for sparse arrays with
#         a BlockIndex. Some operations, e.g. concat, will end up with
#         block indices that could be consolidated. Setting this to true will
#         create a new BlockIndex for that array, with consolidated
#         block indices.
#     """

#     _check_isinstance(left, right, pd.SparseArray)

#     assert_numpy_array_equal(left.sp_values, right.sp_values, check_dtype=check_dtype)

#     # SparseIndex comparison
#     assert isinstance(left.sp_index, pd._libs.sparse.SparseIndex)
#     assert isinstance(right.sp_index, pd._libs.sparse.SparseIndex)

#     if not check_kind:
#         left_index = left.sp_index.to_block_index()
#         right_index = right.sp_index.to_block_index()
#     else:
#         left_index = left.sp_index
#         right_index = right.sp_index

#     if consolidate_block_indices and left.kind == "block":
#         # we'll probably remove this hack...
#         left_index = left_index.to_int_index().to_block_index()
#         right_index = right_index.to_int_index().to_block_index()

#     if not left_index.equals(right_index):
#         raise_assert_detail(
#             "SparseArray.index", "index are not equal", left_index, right_index
#         )
#     else:
#         # Just ensure a
#         pass

#     if check_fill_value:
#         assert_attr_equal("fill_value", left, right)
#     if check_dtype:
#         assert_attr_equal("dtype", left, right)
#     assert_numpy_array_equal(left.to_dense(), right.to_dense(), check_dtype=check_dtype)
