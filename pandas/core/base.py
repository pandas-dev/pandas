"""
Base and utility classes for pandas objects.
"""

import builtins
import textwrap
from typing import Dict, FrozenSet, List, Optional, Union

import numpy as np

import pandas._libs.lib as lib
from pandas._typing import T
from pandas.compat import PYPY
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, Substitution, cache_readonly, doc
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import is_nested_object
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_dict_like,
    is_extension_array_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
    needs_i8_conversion,
)
from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core import algorithms, common as com
from pandas.core.accessor import DirNamesMixin
from pandas.core.algorithms import duplicated, unique1d, value_counts
from pandas.core.arrays import ExtensionArray
from pandas.core.construction import create_series_with_explicit_dtype
import pandas.core.nanops as nanops

_shared_docs: Dict[str, str] = dict()
_indexops_doc_kwargs = dict(
    klass="IndexOpsMixin",
    inplace="",
    unique="IndexOpsMixin",
    duplicated="IndexOpsMixin",
)


class PandasObject(DirNamesMixin):
    """
    Baseclass for various pandas objects.
    """

    @property
    def _constructor(self):
        """
        Class constructor (for this class it's just `__class__`.
        """
        return type(self)

    def __repr__(self) -> str:
        """
        Return a string representation for a particular object.
        """
        # Should be overwritten by base classes
        return object.__repr__(self)

    def _reset_cache(self, key=None):
        """
        Reset cached properties. If ``key`` is passed, only clears that key.
        """
        if getattr(self, "_cache", None) is None:
            return
        if key is None:
            self._cache.clear()
        else:
            self._cache.pop(key, None)

    def __sizeof__(self):
        """
        Generates the total memory usage for an object that returns
        either a value or Series of values
        """
        if hasattr(self, "memory_usage"):
            mem = self.memory_usage(deep=True)
            return int(mem if is_scalar(mem) else mem.sum())

        # no memory_usage attribute, so fall back to object's 'sizeof'
        return super().__sizeof__()

    def _ensure_type(self: T, obj) -> T:
        """
        Ensure that an object has same type as self.

        Used by type checkers.
        """
        assert isinstance(obj, type(self)), type(obj)
        return obj


class NoNewAttributesMixin:
    """
    Mixin which prevents adding new attributes.

    Prevents additional attributes via xxx.attribute = "something" after a
    call to `self.__freeze()`. Mainly used to prevent the user from using
    wrong attributes on an accessor (`Series.cat/.str/.dt`).

    If you really want to add a new attribute at a later time, you need to use
    `object.__setattr__(self, key, value)`.
    """

    def _freeze(self):
        """
        Prevents setting additional attributes.
        """
        object.__setattr__(self, "__frozen", True)

    # prevent adding any attribute via s.xxx.new_attribute = ...
    def __setattr__(self, key: str, value):
        # _cache is used by a decorator
        # We need to check both 1.) cls.__dict__ and 2.) getattr(self, key)
        # because
        # 1.) getattr is false for attributes that raise errors
        # 2.) cls.__dict__ doesn't traverse into base classes
        if getattr(self, "__frozen", False) and not (
            key == "_cache"
            or key in type(self).__dict__
            or getattr(self, key, None) is not None
        ):
            raise AttributeError(f"You cannot add any new attribute '{key}'")
        object.__setattr__(self, key, value)


class GroupByError(Exception):
    pass


class DataError(GroupByError):
    pass


class SpecificationError(GroupByError):
    pass


class SelectionMixin:
    """
    mixin implementing the selection & aggregation interface on a group-like
    object sub-classes need to define: obj, exclusions
    """

    _selection = None
    _internal_names = ["_cache", "__setstate__"]
    _internal_names_set = set(_internal_names)

    _builtin_table = {builtins.sum: np.sum, builtins.max: np.max, builtins.min: np.min}

    _cython_table = {
        builtins.sum: "sum",
        builtins.max: "max",
        builtins.min: "min",
        np.all: "all",
        np.any: "any",
        np.sum: "sum",
        np.nansum: "sum",
        np.mean: "mean",
        np.nanmean: "mean",
        np.prod: "prod",
        np.nanprod: "prod",
        np.std: "std",
        np.nanstd: "std",
        np.var: "var",
        np.nanvar: "var",
        np.median: "median",
        np.nanmedian: "median",
        np.max: "max",
        np.nanmax: "max",
        np.min: "min",
        np.nanmin: "min",
        np.cumprod: "cumprod",
        np.nancumprod: "cumprod",
        np.cumsum: "cumsum",
        np.nancumsum: "cumsum",
    }

    @property
    def _selection_name(self):
        """
        Return a name for myself;

        This would ideally be called the 'name' property,
        but we cannot conflict with the Series.name property which can be set.
        """
        return self._selection

    @property
    def _selection_list(self):
        if not isinstance(
            self._selection, (list, tuple, ABCSeries, ABCIndexClass, np.ndarray)
        ):
            return [self._selection]
        return self._selection

    @cache_readonly
    def _selected_obj(self):
        if self._selection is None or isinstance(self.obj, ABCSeries):
            return self.obj
        else:
            return self.obj[self._selection]

    @cache_readonly
    def ndim(self) -> int:
        return self._selected_obj.ndim

    @cache_readonly
    def _obj_with_exclusions(self):
        if self._selection is not None and isinstance(self.obj, ABCDataFrame):
            return self.obj.reindex(columns=self._selection_list)

        if len(self.exclusions) > 0:
            return self.obj.drop(self.exclusions, axis=1)
        else:
            return self.obj

    def __getitem__(self, key):
        if self._selection is not None:
            raise IndexError(f"Column(s) {self._selection} already selected")

        if isinstance(key, (list, tuple, ABCSeries, ABCIndexClass, np.ndarray)):
            if len(self.obj.columns.intersection(key)) != len(key):
                bad_keys = list(set(key).difference(self.obj.columns))
                raise KeyError(f"Columns not found: {str(bad_keys)[1:-1]}")
            return self._gotitem(list(key), ndim=2)

        elif not getattr(self, "as_index", False):
            if key not in self.obj.columns:
                raise KeyError(f"Column not found: {key}")
            return self._gotitem(key, ndim=2)

        else:
            if key not in self.obj:
                raise KeyError(f"Column not found: {key}")
            return self._gotitem(key, ndim=1)

    def _gotitem(self, key, ndim: int, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : str / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """
        raise AbstractMethodError(self)

    def aggregate(self, func, *args, **kwargs):
        raise AbstractMethodError(self)

    agg = aggregate

    def _try_aggregate_string_function(self, arg: str, *args, **kwargs):
        """
        if arg is a string, then try to operate on it:
        - try to find a function (or attribute) on ourselves
        - try to find a numpy function
        - raise
        """
        assert isinstance(arg, str)

        f = getattr(self, arg, None)
        if f is not None:
            if callable(f):
                return f(*args, **kwargs)

            # people may try to aggregate on a non-callable attribute
            # but don't let them think they can pass args to it
            assert len(args) == 0
            assert len([kwarg for kwarg in kwargs if kwarg not in ["axis"]]) == 0
            return f

        f = getattr(np, arg, None)
        if f is not None:
            if hasattr(self, "__array__"):
                # in particular exclude Window
                return f(self, *args, **kwargs)

        raise AttributeError(
            f"'{arg}' is not a valid function for '{type(self).__name__}' object"
        )

    def _aggregate(self, arg, *args, **kwargs):
        """
        provide an implementation for the aggregators

        Parameters
        ----------
        arg : string, dict, function
        *args : args to pass on to the function
        **kwargs : kwargs to pass on to the function

        Returns
        -------
        tuple of result, how

        Notes
        -----
        how can be a string describe the required post-processing, or
        None if not required
        """
        is_aggregator = lambda x: isinstance(x, (list, tuple, dict))

        _axis = kwargs.pop("_axis", None)
        if _axis is None:
            _axis = getattr(self, "axis", 0)

        if isinstance(arg, str):
            return self._try_aggregate_string_function(arg, *args, **kwargs), None

        if isinstance(arg, dict):
            # aggregate based on the passed dict
            if _axis != 0:  # pragma: no cover
                raise ValueError("Can only pass dict with axis=0")

            obj = self._selected_obj

            # if we have a dict of any non-scalars
            # eg. {'A' : ['mean']}, normalize all to
            # be list-likes
            if any(is_aggregator(x) for x in arg.values()):
                new_arg = {}
                for k, v in arg.items():
                    if not isinstance(v, (tuple, list, dict)):
                        new_arg[k] = [v]
                    else:
                        new_arg[k] = v

                    # the keys must be in the columns
                    # for ndim=2, or renamers for ndim=1

                    # ok for now, but deprecated
                    # {'A': { 'ra': 'mean' }}
                    # {'A': { 'ra': ['mean'] }}
                    # {'ra': ['mean']}

                    # not ok
                    # {'ra' : { 'A' : 'mean' }}
                    if isinstance(v, dict):
                        raise SpecificationError("nested renamer is not supported")
                    elif isinstance(obj, ABCSeries):
                        raise SpecificationError("nested renamer is not supported")
                    elif isinstance(obj, ABCDataFrame) and k not in obj.columns:
                        raise KeyError(f"Column '{k}' does not exist!")

                arg = new_arg

            else:
                # deprecation of renaming keys
                # GH 15931
                keys = list(arg.keys())
                if isinstance(obj, ABCDataFrame) and len(
                    obj.columns.intersection(keys)
                ) != len(keys):
                    raise SpecificationError("nested renamer is not supported")

            from pandas.core.reshape.concat import concat

            def _agg_1dim(name, how, subset=None):
                """
                aggregate a 1-dim with how
                """
                colg = self._gotitem(name, ndim=1, subset=subset)
                if colg.ndim != 1:
                    raise SpecificationError(
                        "nested dictionary is ambiguous in aggregation"
                    )
                return colg.aggregate(how)

            def _agg_2dim(name, how):
                """
                aggregate a 2-dim with how
                """
                colg = self._gotitem(self._selection, ndim=2, subset=obj)
                return colg.aggregate(how)

            def _agg(arg, func):
                """
                run the aggregations over the arg with func
                return a dict
                """
                result = {}
                for fname, agg_how in arg.items():
                    result[fname] = func(fname, agg_how)
                return result

            # set the final keys
            keys = list(arg.keys())
            result = {}

            if self._selection is not None:

                sl = set(self._selection_list)

                # we are a Series like object,
                # but may have multiple aggregations
                if len(sl) == 1:

                    result = _agg(
                        arg, lambda fname, agg_how: _agg_1dim(self._selection, agg_how)
                    )

                # we are selecting the same set as we are aggregating
                elif not len(sl - set(keys)):

                    result = _agg(arg, _agg_1dim)

                # we are a DataFrame, with possibly multiple aggregations
                else:

                    result = _agg(arg, _agg_2dim)

            # no selection
            else:

                try:
                    result = _agg(arg, _agg_1dim)
                except SpecificationError:

                    # we are aggregating expecting all 1d-returns
                    # but we have 2d
                    result = _agg(arg, _agg_2dim)

            # combine results

            def is_any_series() -> bool:
                # return a boolean if we have *any* nested series
                return any(isinstance(r, ABCSeries) for r in result.values())

            def is_any_frame() -> bool:
                # return a boolean if we have *any* nested series
                return any(isinstance(r, ABCDataFrame) for r in result.values())

            if isinstance(result, list):
                return concat(result, keys=keys, axis=1, sort=True), True

            elif is_any_frame():
                # we have a dict of DataFrames
                # return a MI DataFrame

                return concat([result[k] for k in keys], keys=keys, axis=1), True

            elif isinstance(self, ABCSeries) and is_any_series():

                # we have a dict of Series
                # return a MI Series
                try:
                    result = concat(result)
                except TypeError:
                    # we want to give a nice error here if
                    # we have non-same sized objects, so
                    # we don't automatically broadcast

                    raise ValueError(
                        "cannot perform both aggregation "
                        "and transformation operations "
                        "simultaneously"
                    )

                return result, True

            # fall thru
            from pandas import DataFrame, Series

            try:
                result = DataFrame(result)
            except ValueError:

                # we have a dict of scalars
                result = Series(result, name=getattr(self, "name", None))

            return result, True
        elif is_list_like(arg):
            # we require a list, but not an 'str'
            return self._aggregate_multiple_funcs(arg, _axis=_axis), None
        else:
            result = None

        f = self._get_cython_func(arg)
        if f and not args and not kwargs:
            return getattr(self, f)(), None

        # caller can react
        return result, True

    def _aggregate_multiple_funcs(self, arg, _axis):
        from pandas.core.reshape.concat import concat

        if _axis != 0:
            raise NotImplementedError("axis other than 0 is not supported")

        if self._selected_obj.ndim == 1:
            obj = self._selected_obj
        else:
            obj = self._obj_with_exclusions

        results = []
        keys = []

        # degenerate case
        if obj.ndim == 1:
            for a in arg:
                colg = self._gotitem(obj.name, ndim=1, subset=obj)
                try:
                    new_res = colg.aggregate(a)

                except TypeError:
                    pass
                else:
                    results.append(new_res)

                    # make sure we find a good name
                    name = com.get_callable_name(a) or a
                    keys.append(name)

        # multiples
        else:
            for index, col in enumerate(obj):
                colg = self._gotitem(col, ndim=1, subset=obj.iloc[:, index])
                try:
                    new_res = colg.aggregate(arg)
                except (TypeError, DataError):
                    pass
                except ValueError as err:
                    # cannot aggregate
                    if "Must produce aggregated value" in str(err):
                        # raised directly in _aggregate_named
                        pass
                    elif "no results" in str(err):
                        # raised direcly in _aggregate_multiple_funcs
                        pass
                    else:
                        raise
                else:
                    results.append(new_res)
                    keys.append(col)

        # if we are empty
        if not len(results):
            raise ValueError("no results")

        try:
            return concat(results, keys=keys, axis=1, sort=False)
        except TypeError:

            # we are concatting non-NDFrame objects,
            # e.g. a list of scalars

            from pandas import Series

            result = Series(results, index=keys, name=self.name)
            if is_nested_object(result):
                raise ValueError("cannot combine transform and aggregation operations")
            return result

    def _get_cython_func(self, arg: str) -> Optional[str]:
        """
        if we define an internal function for this argument, return it
        """
        return self._cython_table.get(arg)

    def _is_builtin_func(self, arg):
        """
        if we define an builtin function for this argument, return it,
        otherwise return the arg
        """
        return self._builtin_table.get(arg, arg)


class ShallowMixin:
    _attributes: List[str] = []

    def _shallow_copy(self, obj, **kwargs):
        """
        return a new object with the replacement attributes
        """
        if isinstance(obj, self._constructor):
            obj = obj.obj
        for attr in self._attributes:
            if attr not in kwargs:
                kwargs[attr] = getattr(self, attr)
        return self._constructor(obj, **kwargs)


class IndexOpsMixin:
    """
    Common ops mixin to support a unified interface / docs for Series / Index
    """

    # ndarray compatibility
    __array_priority__ = 1000
    _deprecations: FrozenSet[str] = frozenset(
        ["tolist"]  # tolist is not deprecated, just suppressed in the __dir__
    )

    @property
    def _values(self) -> Union[ExtensionArray, np.ndarray]:
        # must be defined here as a property for mypy
        raise AbstractMethodError(self)

    def transpose(self, *args, **kwargs):
        """
        Return the transpose, which is by definition self.

        Returns
        -------
        %(klass)s
        """
        nv.validate_transpose(args, kwargs)
        return self

    T = property(
        transpose,
        doc="""
        Return the transpose, which is by definition self.
        """,
    )

    @property
    def shape(self):
        """
        Return a tuple of the shape of the underlying data.
        """
        return self._values.shape

    def __len__(self) -> int:
        # We need this defined here for mypy
        raise AbstractMethodError(self)

    @property
    def ndim(self) -> int:
        """
        Number of dimensions of the underlying data, by definition 1.
        """
        return 1

    def item(self):
        """
        Return the first element of the underlying data as a python scalar.

        Returns
        -------
        scalar
            The first element of %(klass)s.

        Raises
        ------
        ValueError
            If the data is not length-1.
        """
        if not (
            is_extension_array_dtype(self.dtype) or needs_i8_conversion(self.dtype)
        ):
            # numpy returns ints instead of datetime64/timedelta64 objects,
            #  which we need to wrap in Timestamp/Timedelta/Period regardless.
            return self.values.item()

        if len(self) == 1:
            return next(iter(self))
        raise ValueError("can only convert an array of size 1 to a Python scalar")

    @property
    def nbytes(self) -> int:
        """
        Return the number of bytes in the underlying data.
        """
        return self._values.nbytes

    @property
    def size(self) -> int:
        """
        Return the number of elements in the underlying data.
        """
        return len(self._values)

    @property
    def array(self) -> ExtensionArray:
        """
        The ExtensionArray of the data backing this Series or Index.

        .. versionadded:: 0.24.0

        Returns
        -------
        ExtensionArray
            An ExtensionArray of the values stored within. For extension
            types, this is the actual array. For NumPy native types, this
            is a thin (no copy) wrapper around :class:`numpy.ndarray`.

            ``.array`` differs ``.values`` which may require converting the
            data to a different form.

        See Also
        --------
        Index.to_numpy : Similar method that always returns a NumPy array.
        Series.to_numpy : Similar method that always returns a NumPy array.

        Notes
        -----
        This table lays out the different array types for each extension
        dtype within pandas.

        ================== =============================
        dtype              array type
        ================== =============================
        category           Categorical
        period             PeriodArray
        interval           IntervalArray
        IntegerNA          IntegerArray
        string             StringArray
        boolean            BooleanArray
        datetime64[ns, tz] DatetimeArray
        ================== =============================

        For any 3rd-party extension types, the array type will be an
        ExtensionArray.

        For all remaining dtypes ``.array`` will be a
        :class:`arrays.NumpyExtensionArray` wrapping the actual ndarray
        stored within. If you absolutely need a NumPy array (possibly with
        copying / coercing data), then use :meth:`Series.to_numpy` instead.

        Examples
        --------
        For regular NumPy types like int, and float, a PandasArray
        is returned.

        >>> pd.Series([1, 2, 3]).array
        <PandasArray>
        [1, 2, 3]
        Length: 3, dtype: int64

        For extension types, like Categorical, the actual ExtensionArray
        is returned

        >>> ser = pd.Series(pd.Categorical(['a', 'b', 'a']))
        >>> ser.array
        [a, b, a]
        Categories (2, object): [a, b]
        """
        raise AbstractMethodError(self)

    def to_numpy(self, dtype=None, copy=False, na_value=lib.no_default, **kwargs):
        """
        A NumPy ndarray representing the values in this Series or Index.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        dtype : str or numpy.dtype, optional
            The dtype to pass to :meth:`numpy.asarray`.
        copy : bool, default False
            Whether to ensure that the returned value is a not a view on
            another array. Note that ``copy=False`` does not *ensure* that
            ``to_numpy()`` is no-copy. Rather, ``copy=True`` ensure that
            a copy is made, even if not strictly necessary.
        na_value : Any, optional
            The value to use for missing values. The default value depends
            on `dtype` and the type of the array.

            .. versionadded:: 1.0.0

        **kwargs
            Additional keywords passed through to the ``to_numpy`` method
            of the underlying array (for extension arrays).

            .. versionadded:: 1.0.0

        Returns
        -------
        numpy.ndarray

        See Also
        --------
        Series.array : Get the actual data stored within.
        Index.array : Get the actual data stored within.
        DataFrame.to_numpy : Similar method for DataFrame.

        Notes
        -----
        The returned array will be the same up to equality (values equal
        in `self` will be equal in the returned array; likewise for values
        that are not equal). When `self` contains an ExtensionArray, the
        dtype may be different. For example, for a category-dtype Series,
        ``to_numpy()`` will return a NumPy array and the categorical dtype
        will be lost.

        For NumPy dtypes, this will be a reference to the actual data stored
        in this Series or Index (assuming ``copy=False``). Modifying the result
        in place will modify the data stored in the Series or Index (not that
        we recommend doing that).

        For extension types, ``to_numpy()`` *may* require copying data and
        coercing the result to a NumPy type (possibly object), which may be
        expensive. When you need a no-copy reference to the underlying data,
        :attr:`Series.array` should be used instead.

        This table lays out the different dtypes and default return types of
        ``to_numpy()`` for various dtypes within pandas.

        ================== ================================
        dtype              array type
        ================== ================================
        category[T]        ndarray[T] (same dtype as input)
        period             ndarray[object] (Periods)
        interval           ndarray[object] (Intervals)
        IntegerNA          ndarray[object]
        datetime64[ns]     datetime64[ns]
        datetime64[ns, tz] ndarray[object] (Timestamps)
        ================== ================================

        Examples
        --------
        >>> ser = pd.Series(pd.Categorical(['a', 'b', 'a']))
        >>> ser.to_numpy()
        array(['a', 'b', 'a'], dtype=object)

        Specify the `dtype` to control how datetime-aware data is represented.
        Use ``dtype=object`` to return an ndarray of pandas :class:`Timestamp`
        objects, each with the correct ``tz``.

        >>> ser = pd.Series(pd.date_range('2000', periods=2, tz="CET"))
        >>> ser.to_numpy(dtype=object)
        array([Timestamp('2000-01-01 00:00:00+0100', tz='CET', freq='D'),
               Timestamp('2000-01-02 00:00:00+0100', tz='CET', freq='D')],
              dtype=object)

        Or ``dtype='datetime64[ns]'`` to return an ndarray of native
        datetime64 values. The values are converted to UTC and the timezone
        info is dropped.

        >>> ser.to_numpy(dtype="datetime64[ns]")
        ... # doctest: +ELLIPSIS
        array(['1999-12-31T23:00:00.000000000', '2000-01-01T23:00:00...'],
              dtype='datetime64[ns]')
        """
        if is_extension_array_dtype(self.dtype):
            return self.array.to_numpy(dtype, copy=copy, na_value=na_value, **kwargs)
        elif kwargs:
            bad_keys = list(kwargs.keys())[0]
            raise TypeError(
                f"to_numpy() got an unexpected keyword argument '{bad_keys}'"
            )

        result = np.asarray(self._values, dtype=dtype)
        # TODO(GH-24345): Avoid potential double copy
        if copy or na_value is not lib.no_default:
            result = result.copy()
            if na_value is not lib.no_default:
                result[self.isna()] = na_value
        return result

    @property
    def _ndarray_values(self) -> np.ndarray:
        """
        The data as an ndarray, possibly losing information.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.

        - categorical -> codes
        """
        if is_extension_array_dtype(self):
            return self.array._ndarray_values
        # As a mixin, we depend on the mixing class having values.
        # Special mixin syntax may be developed in the future:
        # https://github.com/python/typing/issues/246
        return self.values  # type: ignore

    @property
    def empty(self):
        return not self.size

    def max(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the maximum value of the Index.

        Parameters
        ----------
        axis : int, optional
            For compatibility with NumPy. Only 0 or None are allowed.
        skipna : bool, default True

        Returns
        -------
        scalar
            Maximum value.

        See Also
        --------
        Index.min : Return the minimum value in an Index.
        Series.max : Return the maximum value in a Series.
        DataFrame.max : Return the maximum values in a DataFrame.

        Examples
        --------
        >>> idx = pd.Index([3, 2, 1])
        >>> idx.max()
        3

        >>> idx = pd.Index(['c', 'b', 'a'])
        >>> idx.max()
        'c'

        For a MultiIndex, the maximum is determined lexicographically.

        >>> idx = pd.MultiIndex.from_product([('a', 'b'), (2, 1)])
        >>> idx.max()
        ('b', 2)
        """
        nv.validate_minmax_axis(axis)
        nv.validate_max(args, kwargs)
        return nanops.nanmax(self._values, skipna=skipna)

    def argmax(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return int position of the largest value in the Series.

        If the maximum is achieved in multiple locations,
        the first row position is returned.

        Parameters
        ----------
        axis : {None}
            Dummy argument for consistency with Series.
        skipna : bool, default True
            Exclude NA/null values when showing the result.
        *args, **kwargs
            Additional arguments and keywords for compatibility with NumPy.

        Returns
        -------
        int
            Row position of the maximum values.

        See Also
        --------
        numpy.ndarray.argmax : Equivalent method for numpy arrays.
        Series.argmin : Similar method, but returning the minimum.
        Series.idxmax : Return index label of the maximum values.
        Series.idxmin : Return index label of the minimum values.

        Examples
        --------
        Consider dataset containing cereal calories

        >>> s = pd.Series({'Corn Flakes': 100.0, 'Almond Delight': 110.0,
        ...                'Cinnamon Toast Crunch': 120.0, 'Cocoa Puff': 110.0})
        >>> s
        Corn Flakes              100.0
        Almond Delight           110.0
        Cinnamon Toast Crunch    120.0
        Cocoa Puff               110.0
        dtype: float64

        >>> s.argmax()
        2

        The maximum cereal calories is in the third element,
        since series is zero-indexed.
        """
        nv.validate_minmax_axis(axis)
        nv.validate_argmax_with_skipna(skipna, args, kwargs)
        return nanops.nanargmax(self._values, skipna=skipna)

    def min(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the minimum value of the Index.

        Parameters
        ----------
        axis : {None}
            Dummy argument for consistency with Series.
        skipna : bool, default True

        Returns
        -------
        scalar
            Minimum value.

        See Also
        --------
        Index.max : Return the maximum value of the object.
        Series.min : Return the minimum value in a Series.
        DataFrame.min : Return the minimum values in a DataFrame.

        Examples
        --------
        >>> idx = pd.Index([3, 2, 1])
        >>> idx.min()
        1

        >>> idx = pd.Index(['c', 'b', 'a'])
        >>> idx.min()
        'a'

        For a MultiIndex, the minimum is determined lexicographically.

        >>> idx = pd.MultiIndex.from_product([('a', 'b'), (2, 1)])
        >>> idx.min()
        ('a', 1)
        """
        nv.validate_minmax_axis(axis)
        nv.validate_min(args, kwargs)
        return nanops.nanmin(self._values, skipna=skipna)

    def argmin(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return a ndarray of the minimum argument indexer.

        Parameters
        ----------
        axis : {None}
            Dummy argument for consistency with Series.
        skipna : bool, default True

        Returns
        -------
        numpy.ndarray

        See Also
        --------
        numpy.ndarray.argmin
        """
        nv.validate_minmax_axis(axis)
        nv.validate_argmax_with_skipna(skipna, args, kwargs)
        return nanops.nanargmin(self._values, skipna=skipna)

    def tolist(self):
        """
        Return a list of the values.

        These are each a scalar type, which is a Python scalar
        (for str, int, float) or a pandas scalar
        (for Timestamp/Timedelta/Interval/Period)

        Returns
        -------
        list

        See Also
        --------
        numpy.ndarray.tolist
        """
        if not isinstance(self._values, np.ndarray):
            # check for ndarray instead of dtype to catch DTA/TDA
            return list(self._values)
        return self._values.tolist()

    to_list = tolist

    def __iter__(self):
        """
        Return an iterator of the values.

        These are each a scalar type, which is a Python scalar
        (for str, int, float) or a pandas scalar
        (for Timestamp/Timedelta/Interval/Period)

        Returns
        -------
        iterator
        """
        # We are explicitly making element iterators.
        if not isinstance(self._values, np.ndarray):
            # Check type instead of dtype to catch DTA/TDA
            return iter(self._values)
        else:
            return map(self._values.item, range(self._values.size))

    @cache_readonly
    def hasnans(self):
        """
        Return if I have any nans; enables various perf speedups.
        """
        return bool(isna(self).any())

    def _reduce(
        self,
        op,
        name: str,
        axis=0,
        skipna=True,
        numeric_only=None,
        filter_type=None,
        **kwds,
    ):
        """
        Perform the reduction type operation if we can.
        """
        func = getattr(self, name, None)
        if func is None:
            raise TypeError(
                f"{type(self).__name__} cannot perform the operation {name}"
            )
        return func(skipna=skipna, **kwds)

    def _map_values(self, mapper, na_action=None):
        """
        An internal function that maps values using the input
        correspondence (which can be a dict, Series, or function).

        Parameters
        ----------
        mapper : function, dict, or Series
            The input correspondence object
        na_action : {None, 'ignore'}
            If 'ignore', propagate NA values, without passing them to the
            mapping function

        Returns
        -------
        Union[Index, MultiIndex], inferred
            The output of the mapping function applied to the index.
            If the function returns a tuple with more than one element
            a MultiIndex will be returned.
        """
        # we can fastpath dict/Series to an efficient map
        # as we know that we are not going to have to yield
        # python types
        if is_dict_like(mapper):
            if isinstance(mapper, dict) and hasattr(mapper, "__missing__"):
                # If a dictionary subclass defines a default value method,
                # convert mapper to a lookup function (GH #15999).
                dict_with_default = mapper
                mapper = lambda x: dict_with_default[x]
            else:
                # Dictionary does not have a default. Thus it's safe to
                # convert to an Series for efficiency.
                # we specify the keys here to handle the
                # possibility that they are tuples

                # The return value of mapping with an empty mapper is
                # expected to be pd.Series(np.nan, ...). As np.nan is
                # of dtype float64 the return value of this method should
                # be float64 as well
                mapper = create_series_with_explicit_dtype(
                    mapper, dtype_if_empty=np.float64
                )

        if isinstance(mapper, ABCSeries):
            # Since values were input this means we came from either
            # a dict or a series and mapper should be an index
            if is_categorical_dtype(self._values):
                # use the built in categorical series mapper which saves
                # time by mapping the categories instead of all values
                return self._values.map(mapper)
            if is_extension_array_dtype(self.dtype):
                values = self._values
            else:
                values = self.values

            indexer = mapper.index.get_indexer(values)
            new_values = algorithms.take_1d(mapper._values, indexer)

            return new_values

        # we must convert to python types
        if is_extension_array_dtype(self.dtype) and hasattr(self._values, "map"):
            # GH#23179 some EAs do not have `map`
            values = self._values
            if na_action is not None:
                raise NotImplementedError
            map_f = lambda values, f: values.map(f)
        else:
            values = self.astype(object)
            values = getattr(values, "values", values)
            if na_action == "ignore":

                def map_f(values, f):
                    return lib.map_infer_mask(values, f, isna(values).view(np.uint8))

            else:
                map_f = lib.map_infer

        # mapper is a function
        new_values = map_f(values, mapper)

        return new_values

    def value_counts(
        self, normalize=False, sort=True, ascending=False, bins=None, dropna=True
    ):
        """
        Return a Series containing counts of unique values.

        The resulting object will be in descending order so that the
        first element is the most frequently-occurring element.
        Excludes NA values by default.

        Parameters
        ----------
        normalize : bool, default False
            If True then the object returned will contain the relative
            frequencies of the unique values.
        sort : bool, default True
            Sort by frequencies.
        ascending : bool, default False
            Sort in ascending order.
        bins : int, optional
            Rather than count values, group them into half-open bins,
            a convenience for ``pd.cut``, only works with numeric data.
        dropna : bool, default True
            Don't include counts of NaN.

        Returns
        -------
        Series

        See Also
        --------
        Series.count: Number of non-NA elements in a Series.
        DataFrame.count: Number of non-NA elements in a DataFrame.
        DataFrame.value_counts: Equivalent method on DataFrames.

        Examples
        --------
        >>> index = pd.Index([3, 1, 2, 3, 4, np.nan])
        >>> index.value_counts()
        3.0    2
        4.0    1
        2.0    1
        1.0    1
        dtype: int64

        With `normalize` set to `True`, returns the relative frequency by
        dividing all values by the sum of values.

        >>> s = pd.Series([3, 1, 2, 3, 4, np.nan])
        >>> s.value_counts(normalize=True)
        3.0    0.4
        4.0    0.2
        2.0    0.2
        1.0    0.2
        dtype: float64

        **bins**

        Bins can be useful for going from a continuous variable to a
        categorical variable; instead of counting unique
        apparitions of values, divide the index in the specified
        number of half-open bins.

        >>> s.value_counts(bins=3)
        (2.0, 3.0]      2
        (0.996, 2.0]    2
        (3.0, 4.0]      1
        dtype: int64

        **dropna**

        With `dropna` set to `False` we can also see NaN index values.

        >>> s.value_counts(dropna=False)
        3.0    2
        NaN    1
        4.0    1
        2.0    1
        1.0    1
        dtype: int64
        """
        result = value_counts(
            self,
            sort=sort,
            ascending=ascending,
            normalize=normalize,
            bins=bins,
            dropna=dropna,
        )
        return result

    def unique(self):
        values = self._values

        if hasattr(values, "unique"):

            result = values.unique()
            if self.dtype.kind in ["m", "M"] and isinstance(self, ABCSeries):
                # GH#31182 Series._values returns EA, unpack for backward-compat
                if getattr(self.dtype, "tz", None) is None:
                    result = np.asarray(result)
        else:
            result = unique1d(values)

        return result

    def nunique(self, dropna: bool = True) -> int:
        """
        Return number of unique elements in the object.

        Excludes NA values by default.

        Parameters
        ----------
        dropna : bool, default True
            Don't include NaN in the count.

        Returns
        -------
        int

        See Also
        --------
        DataFrame.nunique: Method nunique for DataFrame.
        Series.count: Count non-NA/null observations in the Series.

        Examples
        --------
        >>> s = pd.Series([1, 3, 5, 7, 7])
        >>> s
        0    1
        1    3
        2    5
        3    7
        4    7
        dtype: int64

        >>> s.nunique()
        4
        """
        uniqs = self.unique()
        n = len(uniqs)
        if dropna and isna(uniqs).any():
            n -= 1
        return n

    @property
    def is_unique(self) -> bool:
        """
        Return boolean if values in the object are unique.

        Returns
        -------
        bool
        """
        return self.nunique(dropna=False) == len(self)

    @property
    def is_monotonic(self) -> bool:
        """
        Return boolean if values in the object are
        monotonic_increasing.

        Returns
        -------
        bool
        """
        from pandas import Index

        return Index(self).is_monotonic

    @property
    def is_monotonic_increasing(self) -> bool:
        """
        Alias for is_monotonic.
        """
        # mypy complains if we alias directly
        return self.is_monotonic

    @property
    def is_monotonic_decreasing(self) -> bool:
        """
        Return boolean if values in the object are
        monotonic_decreasing.

        Returns
        -------
        bool
        """
        from pandas import Index

        return Index(self).is_monotonic_decreasing

    def memory_usage(self, deep=False):
        """
        Memory usage of the values.

        Parameters
        ----------
        deep : bool
            Introspect the data deeply, interrogate
            `object` dtypes for system-level memory consumption.

        Returns
        -------
        bytes used

        See Also
        --------
        numpy.ndarray.nbytes

        Notes
        -----
        Memory usage does not include memory consumed by elements that
        are not components of the array if deep=False or if used on PyPy
        """
        if hasattr(self.array, "memory_usage"):
            return self.array.memory_usage(deep=deep)

        v = self.array.nbytes
        if deep and is_object_dtype(self) and not PYPY:
            v += lib.memory_usage_of_objects(self.array)
        return v

    @doc(
        algorithms.factorize,
        values="",
        order="",
        size_hint="",
        sort=textwrap.dedent(
            """\
            sort : bool, default False
                Sort `uniques` and shuffle `codes` to maintain the
                relationship.
            """
        ),
    )
    def factorize(self, sort=False, na_sentinel=-1):
        return algorithms.factorize(self, sort=sort, na_sentinel=na_sentinel)

    _shared_docs[
        "searchsorted"
    ] = """
        Find indices where elements should be inserted to maintain order.

        Find the indices into a sorted %(klass)s `self` such that, if the
        corresponding elements in `value` were inserted before the indices,
        the order of `self` would be preserved.

        .. note::

            The %(klass)s *must* be monotonically sorted, otherwise
            wrong locations will likely be returned. Pandas does *not*
            check this for you.

        Parameters
        ----------
        value : array_like
            Values to insert into `self`.
        side : {'left', 'right'}, optional
            If 'left', the index of the first suitable location found is given.
            If 'right', return the last such index.  If there is no suitable
            index, return either 0 or N (where N is the length of `self`).
        sorter : 1-D array_like, optional
            Optional array of integer indices that sort `self` into ascending
            order. They are typically the result of ``np.argsort``.

        Returns
        -------
        int or array of int
            A scalar or array of insertion points with the
            same shape as `value`.

            .. versionchanged:: 0.24.0
                If `value` is a scalar, an int is now always returned.
                Previously, scalar inputs returned an 1-item array for
                :class:`Series` and :class:`Categorical`.

        See Also
        --------
        sort_values
        numpy.searchsorted

        Notes
        -----
        Binary search is used to find the required insertion points.

        Examples
        --------
        >>> x = pd.Series([1, 2, 3])
        >>> x
        0    1
        1    2
        2    3
        dtype: int64

        >>> x.searchsorted(4)
        3

        >>> x.searchsorted([0, 4])
        array([0, 3])

        >>> x.searchsorted([1, 3], side='left')
        array([0, 2])

        >>> x.searchsorted([1, 3], side='right')
        array([1, 3])

        >>> x = pd.Categorical(['apple', 'bread', 'bread',
                                'cheese', 'milk'], ordered=True)
        [apple, bread, bread, cheese, milk]
        Categories (4, object): [apple < bread < cheese < milk]

        >>> x.searchsorted('bread')
        1

        >>> x.searchsorted(['bread'], side='right')
        array([3])

        If the values are not monotonically sorted, wrong locations
        may be returned:

        >>> x = pd.Series([2, 1, 3])
        >>> x.searchsorted(1)
        0  # wrong result, correct would be 1
        """

    @Substitution(klass="Index")
    @Appender(_shared_docs["searchsorted"])
    def searchsorted(self, value, side="left", sorter=None) -> np.ndarray:
        return algorithms.searchsorted(self._values, value, side=side, sorter=sorter)

    def drop_duplicates(self, keep="first", inplace=False):
        inplace = validate_bool_kwarg(inplace, "inplace")
        if isinstance(self, ABCIndexClass):
            if self.is_unique:
                return self._shallow_copy()

        duplicated = self.duplicated(keep=keep)
        result = self[np.logical_not(duplicated)]
        if inplace:
            return self._update_inplace(result)
        else:
            return result

    def duplicated(self, keep="first"):
        if isinstance(self, ABCIndexClass):
            if self.is_unique:
                return np.zeros(len(self), dtype=np.bool)
            return duplicated(self, keep=keep)
        else:
            return self._constructor(
                duplicated(self, keep=keep), index=self.index
            ).__finalize__(self)

    # ----------------------------------------------------------------------
    # abstracts

    def _update_inplace(self, result, verify_is_copy=True, **kwargs):
        raise AbstractMethodError(self)
