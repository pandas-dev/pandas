"""
Data structure for 1-dimensional cross-sectional and time series data
"""
from __future__ import division

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import types
import warnings

from numpy import nan, ndarray
import numpy as np
import numpy.ma as ma

from pandas.core.common import (isnull, notnull, is_bool_indexer,
                                _default_index, _maybe_upcast,
                                _asarray_tuplesafe, _infer_dtype_from_scalar,
                                is_list_like, _values_from_object,
                                is_categorical_dtype,
                                _possibly_cast_to_datetime,
                                _possibly_castable, _possibly_convert_platform,
                                _try_sort, is_extension_type, is_datetimetz,
                                _maybe_match_name, ABCSparseArray,
                                _coerce_to_dtype, SettingWithCopyError,
                                _maybe_box_datetimelike, ABCDataFrame,
                                _dict_compat, is_integer)
from pandas.core.index import (Index, MultiIndex, InvalidIndexError,
                               Float64Index, _ensure_index)
from pandas.core.indexing import check_bool_indexer, maybe_convert_indices
from pandas.core import generic, base
from pandas.core.internals import SingleBlockManager
from pandas.core.categorical import Categorical, CategoricalAccessor
import pandas.core.strings as strings
from pandas.tseries.common import (maybe_to_datetimelike,
                                   CombinedDatetimelikeProperties)
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.tdi import TimedeltaIndex
from pandas.tseries.period import PeriodIndex, Period
from pandas import compat
from pandas.util.terminal import get_terminal_size
from pandas.compat import zip, u, OrderedDict, StringIO
from pandas.compat.numpy import function as nv

import pandas.core.ops as ops
import pandas.core.algorithms as algos

import pandas.core.common as com
import pandas.core.datetools as datetools
import pandas.core.nanops as nanops
import pandas.formats.format as fmt
from pandas.util.decorators import Appender, deprecate_kwarg, Substitution

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.index as _index

from pandas.core.config import get_option

from pandas import _np_version_under1p9

__all__ = ['Series']

_shared_doc_kwargs = dict(
    axes='index', klass='Series', axes_single_arg="{0, 'index'}",
    inplace="""inplace : boolean, default False
        If True, performs operation inplace and returns None.""",
    duplicated='Series')


def _coerce_method(converter):
    """ install the scalar coercion methods """

    def wrapper(self):
        if len(self) == 1:
            return converter(self.iloc[0])
        raise TypeError("cannot convert the series to "
                        "{0}".format(str(converter)))

    return wrapper

# ----------------------------------------------------------------------
# Series class


class Series(base.IndexOpsMixin, strings.StringAccessorMixin,
             generic.NDFrame,):
    """
    One-dimensional ndarray with axis labels (including time series).

    Labels need not be unique but must be any hashable type. The object
    supports both integer- and label-based indexing and provides a host of
    methods for performing operations involving the index. Statistical
    methods from ndarray have been overridden to automatically exclude
    missing data (currently represented as NaN)

    Operations between Series (+, -, /, *, **) align values based on their
    associated index values-- they need not be the same length. The result
    index will be the sorted union of the two indexes.

    Parameters
    ----------
    data : array-like, dict, or scalar value
        Contains data stored in Series
    index : array-like or Index (1d)
        Values must be unique and hashable, same length as data. Index
        object (or other iterable of same length as data) Will default to
        RangeIndex(len(data)) if not provided. If both a dict and index
        sequence are used, the index will override the keys found in the
        dict.
    dtype : numpy.dtype or None
        If None, dtype will be inferred
    copy : boolean, default False
        Copy input data
    """
    _metadata = ['name']
    _accessors = frozenset(['dt', 'cat', 'str'])
    _allow_index_ops = True

    def __init__(self, data=None, index=None, dtype=None, name=None,
                 copy=False, fastpath=False):

        # we are called internally, so short-circuit
        if fastpath:

            # data is an ndarray, index is defined
            if not isinstance(data, SingleBlockManager):
                data = SingleBlockManager(data, index, fastpath=True)
            if copy:
                data = data.copy()
            if index is None:
                index = data.index

        else:

            if index is not None:
                index = _ensure_index(index)

            if data is None:
                data = {}
            if dtype is not None:
                dtype = self._validate_dtype(dtype)

            if isinstance(data, MultiIndex):
                raise NotImplementedError("initializing a Series from a "
                                          "MultiIndex is not supported")
            elif isinstance(data, Index):
                # need to copy to avoid aliasing issues
                if name is None:
                    name = data.name

                data = data._to_embed(keep_tz=True)
                copy = True
            elif isinstance(data, np.ndarray):
                pass
            elif isinstance(data, Series):
                if name is None:
                    name = data.name
                if index is None:
                    index = data.index
                else:
                    data = data.reindex(index, copy=copy)
                data = data._data
            elif isinstance(data, dict):
                if index is None:
                    if isinstance(data, OrderedDict):
                        index = Index(data)
                    else:
                        index = Index(_try_sort(data))
                try:
                    if isinstance(index, DatetimeIndex):
                        if len(data):
                            # coerce back to datetime objects for lookup
                            data = _dict_compat(data)
                            data = lib.fast_multiget(data, index.astype('O'),
                                                     default=np.nan)
                        else:
                            data = np.nan
                    # GH #12169
                    elif isinstance(index, (PeriodIndex, TimedeltaIndex)):
                        data = ([data.get(i, nan) for i in index]
                                if data else np.nan)
                    else:
                        data = lib.fast_multiget(data, index.values,
                                                 default=np.nan)
                except TypeError:
                    data = ([data.get(i, nan) for i in index]
                            if data else np.nan)

            elif isinstance(data, SingleBlockManager):
                if index is None:
                    index = data.index
                else:
                    data = data.reindex(index, copy=copy)
            elif isinstance(data, Categorical):
                # GH12574: Allow dtype=category only, otherwise error
                if ((dtype is not None) and
                        not is_categorical_dtype(dtype)):
                    raise ValueError("cannot specify a dtype with a "
                                     "Categorical unless "
                                     "dtype='category'")
            elif (isinstance(data, types.GeneratorType) or
                  (compat.PY3 and isinstance(data, map))):
                data = list(data)
            elif isinstance(data, (set, frozenset)):
                raise TypeError("{0!r} type is unordered"
                                "".format(data.__class__.__name__))
            else:

                # handle sparse passed here (and force conversion)
                if isinstance(data, ABCSparseArray):
                    data = data.to_dense()

            if index is None:
                if not is_list_like(data):
                    data = [data]
                index = _default_index(len(data))

            # create/copy the manager
            if isinstance(data, SingleBlockManager):
                if dtype is not None:
                    data = data.astype(dtype=dtype, raise_on_error=False)
                elif copy:
                    data = data.copy()
            else:
                data = _sanitize_array(data, index, dtype, copy,
                                       raise_cast_failure=True)

                data = SingleBlockManager(data, index, fastpath=True)

        generic.NDFrame.__init__(self, data, fastpath=True)

        self.name = name
        self._set_axis(0, index, fastpath=True)

    @classmethod
    def from_array(cls, arr, index=None, name=None, dtype=None, copy=False,
                   fastpath=False):
        # return a sparse series here
        if isinstance(arr, ABCSparseArray):
            from pandas.sparse.series import SparseSeries
            cls = SparseSeries

        return cls(arr, index=index, name=name, dtype=dtype, copy=copy,
                   fastpath=fastpath)

    @property
    def _constructor(self):
        return Series

    @property
    def _constructor_expanddim(self):
        from pandas.core.frame import DataFrame
        return DataFrame

    # types
    @property
    def _can_hold_na(self):
        return self._data._can_hold_na

    @property
    def is_time_series(self):
        warnings.warn("is_time_series is deprecated. Please use "
                      "Series.index.is_all_dates", FutureWarning, stacklevel=2)
        # return self._subtyp in ['time_series', 'sparse_time_series']
        return self.index.is_all_dates

    _index = None

    def _set_axis(self, axis, labels, fastpath=False):
        """ override generic, we want to set the _typ here """

        if not fastpath:
            labels = _ensure_index(labels)

        is_all_dates = labels.is_all_dates
        if is_all_dates:

            if not isinstance(labels,
                              (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
                labels = DatetimeIndex(labels)

                # need to set here becuase we changed the index
                if fastpath:
                    self._data.set_axis(axis, labels)
        self._set_subtyp(is_all_dates)

        object.__setattr__(self, '_index', labels)
        if not fastpath:
            self._data.set_axis(axis, labels)

    def _set_subtyp(self, is_all_dates):
        if is_all_dates:
            object.__setattr__(self, '_subtyp', 'time_series')
        else:
            object.__setattr__(self, '_subtyp', 'series')

    def _update_inplace(self, result, **kwargs):
        # we want to call the generic version and not the IndexOpsMixin
        return generic.NDFrame._update_inplace(self, result, **kwargs)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value is not None and not com.is_hashable(value):
            raise TypeError('Series.name must be a hashable type')
        object.__setattr__(self, '_name', value)

    # ndarray compatibility
    @property
    def dtype(self):
        """ return the dtype object of the underlying data """
        return self._data.dtype

    @property
    def dtypes(self):
        """ return the dtype object of the underlying data """
        return self._data.dtype

    @property
    def ftype(self):
        """ return if the data is sparse|dense """
        return self._data.ftype

    @property
    def ftypes(self):
        """ return if the data is sparse|dense """
        return self._data.ftype

    @property
    def values(self):
        """
        Return Series as ndarray or ndarray-like
        depending on the dtype

        Returns
        -------
        arr : numpy.ndarray or ndarray-like

        Examples
        --------
        >>> pd.Series([1, 2, 3]).values
        array([1, 2, 3])

        >>> pd.Series(list('aabc')).values
        array(['a', 'a', 'b', 'c'], dtype=object)

        >>> pd.Series(list('aabc')).astype('category').values
        [a, a, b, c]
        Categories (3, object): [a, b, c]

        Timezone aware datetime data is converted to UTC:

        >>> pd.Series(pd.date_range('20130101', periods=3,
                                    tz='US/Eastern')).values
        array(['2013-01-01T00:00:00.000000000-0500',
               '2013-01-02T00:00:00.000000000-0500',
               '2013-01-03T00:00:00.000000000-0500'], dtype='datetime64[ns]')

        """
        return self._data.external_values()

    @property
    def _values(self):
        """ return the internal repr of this data """
        return self._data.internal_values()

    def get_values(self):
        """ same as values (but handles sparseness conversions); is a view """
        return self._data.get_values()

    @property
    def asobject(self):
        """
        return object Series which contains boxed values

        *this is an internal non-public method*
        """
        return self._data.asobject

    # ops
    def ravel(self, order='C'):
        """
        Return the flattened underlying data as an ndarray

        See also
        --------
        numpy.ndarray.ravel
        """
        return self._values.ravel(order=order)

    def compress(self, condition, *args, **kwargs):
        """
        Return selected slices of an array along given axis as a Series

        See also
        --------
        numpy.ndarray.compress
        """
        nv.validate_compress(args, kwargs)
        return self[condition]

    def nonzero(self):
        """
        Return the indices of the elements that are non-zero

        This method is equivalent to calling `numpy.nonzero` on the
        series data. For compatability with NumPy, the return value is
        the same (a tuple with an array of indices for each dimension),
        but it will always be a one-item tuple because series only have
        one dimension.

        Examples
        --------
        >>> s = pd.Series([0, 3, 0, 4])
        >>> s.nonzero()
        (array([1, 3]),)
        >>> s.iloc[s.nonzero()[0]]
        1    3
        3    4
        dtype: int64

        See Also
        --------
        numpy.nonzero
        """
        return self._values.nonzero()

    def put(self, *args, **kwargs):
        """
        Applies the `put` method to its `values` attribute
        if it has one.

        See also
        --------
        numpy.ndarray.put
        """
        self._values.put(*args, **kwargs)

    def __len__(self):
        """
        return the length of the Series
        """
        return len(self._data)

    def view(self, dtype=None):
        return self._constructor(self._values.view(dtype),
                                 index=self.index).__finalize__(self)

    def __array__(self, result=None):
        """
        the array interface, return my values
        """
        return self.get_values()

    def __array_wrap__(self, result, context=None):
        """
        Gets called after a ufunc
        """
        return self._constructor(result, index=self.index,
                                 copy=False).__finalize__(self)

    def __array_prepare__(self, result, context=None):
        """
        Gets called prior to a ufunc
        """

        # nice error message for non-ufunc types
        if context is not None and not isinstance(self._values, np.ndarray):
            obj = context[1][0]
            raise TypeError("{obj} with dtype {dtype} cannot perform "
                            "the numpy op {op}".format(
                                obj=type(obj).__name__,
                                dtype=getattr(obj, 'dtype', None),
                                op=context[0].__name__))
        return result

    # complex
    @property
    def real(self):
        return self.values.real

    @real.setter
    def real(self, v):
        self.values.real = v

    @property
    def imag(self):
        return self.values.imag

    @imag.setter
    def imag(self, v):
        self.values.imag = v

    # coercion
    __float__ = _coerce_method(float)
    __long__ = _coerce_method(int)
    __int__ = _coerce_method(int)

    def _unpickle_series_compat(self, state):
        if isinstance(state, dict):
            self._data = state['_data']
            self.name = state['name']
            self.index = self._data.index

        elif isinstance(state, tuple):

            # < 0.12 series pickle

            nd_state, own_state = state

            # recreate the ndarray
            data = np.empty(nd_state[1], dtype=nd_state[2])
            np.ndarray.__setstate__(data, nd_state)

            # backwards compat
            index, name = own_state[0], None
            if len(own_state) > 1:
                name = own_state[1]

            # recreate
            self._data = SingleBlockManager(data, index, fastpath=True)
            self._index = index
            self.name = name

        else:
            raise Exception("cannot unpickle legacy formats -> [%s]" % state)

    # indexers
    @property
    def axes(self):
        """Return a list of the row axis labels"""
        return [self.index]

    def _ixs(self, i, axis=0):
        """
        Return the i-th value or values in the Series by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Returns
        -------
        value : scalar (int) or Series (slice, sequence)
        """
        try:

            # dispatch to the values if we need
            values = self._values
            if isinstance(values, np.ndarray):
                return _index.get_value_at(values, i)
            else:
                return values[i]
        except IndexError:
            raise
        except:
            if isinstance(i, slice):
                indexer = self.index._convert_slice_indexer(i, kind='iloc')
                return self._get_values(indexer)
            else:
                label = self.index[i]
                if isinstance(label, Index):
                    return self.take(i, axis=axis, convert=True)
                else:
                    return _index.get_value_at(self, i)

    @property
    def _is_mixed_type(self):
        return False

    def _slice(self, slobj, axis=0, kind=None):
        slobj = self.index._convert_slice_indexer(slobj,
                                                  kind=kind or 'getitem')
        return self._get_values(slobj)

    def __getitem__(self, key):
        key = com._apply_if_callable(key, self)
        try:
            result = self.index.get_value(self, key)

            if not lib.isscalar(result):
                if is_list_like(result) and not isinstance(result, Series):

                    # we need to box if we have a non-unique index here
                    # otherwise have inline ndarray/lists
                    if not self.index.is_unique:
                        result = self._constructor(
                            result, index=[key] * len(result),
                            dtype=self.dtype).__finalize__(self)

            return result
        except InvalidIndexError:
            pass
        except (KeyError, ValueError):
            if isinstance(key, tuple) and isinstance(self.index, MultiIndex):
                # kludge
                pass
            elif key is Ellipsis:
                return self
            elif is_bool_indexer(key):
                pass
            else:

                # we can try to coerce the indexer (or this will raise)
                new_key = self.index._convert_scalar_indexer(key,
                                                             kind='getitem')
                if type(new_key) != type(key):
                    return self.__getitem__(new_key)
                raise

        except Exception:
            raise

        if com.is_iterator(key):
            key = list(key)

        if is_bool_indexer(key):
            key = check_bool_indexer(self.index, key)

        return self._get_with(key)

    def _get_with(self, key):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind='getitem')
            return self._get_values(indexer)
        elif isinstance(key, ABCDataFrame):
            raise TypeError('Indexing a Series with DataFrame is not '
                            'supported, use the appropriate DataFrame column')
        else:
            if isinstance(key, tuple):
                try:
                    return self._get_values_tuple(key)
                except:
                    if len(key) == 1:
                        key = key[0]
                        if isinstance(key, slice):
                            return self._get_values(key)
                    raise

            # pragma: no cover
            if not isinstance(key, (list, np.ndarray, Series, Index)):
                key = list(key)

            if isinstance(key, Index):
                key_type = key.inferred_type
            else:
                key_type = lib.infer_dtype(key)

            if key_type == 'integer':
                if self.index.is_integer() or self.index.is_floating():
                    return self.reindex(key)
                else:
                    return self._get_values(key)
            elif key_type == 'boolean':
                return self._get_values(key)
            else:
                try:
                    # handle the dup indexing case (GH 4246)
                    if isinstance(key, (list, tuple)):
                        return self.ix[key]

                    return self.reindex(key)
                except Exception:
                    # [slice(0, 5, None)] will break if you convert to ndarray,
                    # e.g. as requested by np.median
                    # hack
                    if isinstance(key[0], slice):
                        return self._get_values(key)
                    raise

    def _get_values_tuple(self, key):
        # mpl hackaround
        if any(k is None for k in key):
            return self._get_values(key)

        if not isinstance(self.index, MultiIndex):
            raise ValueError('Can only tuple-index with a MultiIndex')

        # If key is contained, would have returned by now
        indexer, new_index = self.index.get_loc_level(key)
        return self._constructor(self._values[indexer],
                                 index=new_index).__finalize__(self)

    def _get_values(self, indexer):
        try:
            return self._constructor(self._data.get_slice(indexer),
                                     fastpath=True).__finalize__(self)
        except Exception:
            return self._values[indexer]

    def __setitem__(self, key, value):
        key = com._apply_if_callable(key, self)

        def setitem(key, value):
            try:
                self._set_with_engine(key, value)
                return
            except (SettingWithCopyError):
                raise
            except (KeyError, ValueError):
                values = self._values
                if (is_integer(key) and
                        not self.index.inferred_type == 'integer'):

                    values[key] = value
                    return
                elif key is Ellipsis:
                    self[:] = value
                    return
                elif is_bool_indexer(key):
                    pass
                elif com.is_timedelta64_dtype(self.dtype):
                    # reassign a null value to iNaT
                    if isnull(value):
                        value = tslib.iNaT

                        try:
                            self.index._engine.set_value(self._values, key,
                                                         value)
                            return
                        except TypeError:
                            pass

                self.loc[key] = value
                return

            except TypeError as e:
                if (isinstance(key, tuple) and
                        not isinstance(self.index, MultiIndex)):
                    raise ValueError("Can only tuple-index with a MultiIndex")

                # python 3 type errors should be raised
                if 'unorderable' in str(e):  # pragma: no cover
                    raise IndexError(key)

            if is_bool_indexer(key):
                key = check_bool_indexer(self.index, key)
                try:
                    self.where(~key, value, inplace=True)
                    return
                except InvalidIndexError:
                    pass

            self._set_with(key, value)

        # do the setitem
        cacher_needs_updating = self._check_is_chained_assignment_possible()
        setitem(key, value)
        if cacher_needs_updating:
            self._maybe_update_cacher()

    def _set_with_engine(self, key, value):
        values = self._values
        try:
            self.index._engine.set_value(values, key, value)
            return
        except KeyError:
            values[self.index.get_loc(key)] = value
            return

    def _set_with(self, key, value):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind='getitem')
            return self._set_values(indexer, value)
        else:
            if isinstance(key, tuple):
                try:
                    self._set_values(key, value)
                except Exception:
                    pass

            if not isinstance(key, (list, Series, np.ndarray, Series)):
                try:
                    key = list(key)
                except:
                    key = [key]

            if isinstance(key, Index):
                key_type = key.inferred_type
            else:
                key_type = lib.infer_dtype(key)

            if key_type == 'integer':
                if self.index.inferred_type == 'integer':
                    self._set_labels(key, value)
                else:
                    return self._set_values(key, value)
            elif key_type == 'boolean':
                self._set_values(key.astype(np.bool_), value)
            else:
                self._set_labels(key, value)

    def _set_labels(self, key, value):
        if isinstance(key, Index):
            key = key.values
        else:
            key = _asarray_tuplesafe(key)
        indexer = self.index.get_indexer(key)
        mask = indexer == -1
        if mask.any():
            raise ValueError('%s not contained in the index' % str(key[mask]))
        self._set_values(indexer, value)

    def _set_values(self, key, value):
        if isinstance(key, Series):
            key = key._values
        self._data = self._data.setitem(indexer=key, value=value)
        self._maybe_update_cacher()

    def repeat(self, reps, *args, **kwargs):
        """
        Repeat elements of an Series. Refer to `numpy.ndarray.repeat`
        for more information about the `reps` argument.

        See also
        --------
        numpy.ndarray.repeat
        """
        nv.validate_repeat(args, kwargs)
        new_index = self.index.repeat(reps)
        new_values = self._values.repeat(reps)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def reshape(self, *args, **kwargs):
        """
        Return the values attribute of `self` with shape `args`.
        However, if the specified shape matches exactly the current
        shape, `self` is returned for compatibility reasons.

        See also
        --------
        numpy.ndarray.reshape
        """
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            shape = args[0]
        else:
            shape = args

        if tuple(shape) == self.shape:
            # XXX ignoring the "order" keyword.
            nv.validate_reshape(tuple(), kwargs)
            return self

        return self._values.reshape(shape, **kwargs)

    def iget_value(self, i, axis=0):
        """
        DEPRECATED. Use ``.iloc[i]`` or ``.iat[i]`` instead
        """
        warnings.warn("iget_value(i) is deprecated. Please use .iloc[i] or "
                      ".iat[i]", FutureWarning, stacklevel=2)
        return self._ixs(i)

    def iget(self, i, axis=0):
        """
        DEPRECATED. Use ``.iloc[i]`` or ``.iat[i]`` instead
        """

        warnings.warn("iget(i) is deprecated. Please use .iloc[i] or .iat[i]",
                      FutureWarning, stacklevel=2)
        return self._ixs(i)

    def irow(self, i, axis=0):
        """
        DEPRECATED. Use ``.iloc[i]`` or ``.iat[i]`` instead
        """
        warnings.warn("irow(i) is deprecated. Please use .iloc[i] or .iat[i]",
                      FutureWarning, stacklevel=2)
        return self._ixs(i)

    def get_value(self, label, takeable=False):
        """
        Quickly retrieve single value at passed index label

        Parameters
        ----------
        index : label
        takeable : interpret the index as indexers, default False

        Returns
        -------
        value : scalar value
        """
        if takeable is True:
            return _maybe_box_datetimelike(self._values[label])
        return self.index.get_value(self._values, label)

    def set_value(self, label, value, takeable=False):
        """
        Quickly set single value at passed label. If label is not contained, a
        new object is created with the label placed at the end of the result
        index

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed
        value : object
            Scalar value
        takeable : interpret the index as indexers, default False

        Returns
        -------
        series : Series
            If label is contained, will be reference to calling Series,
            otherwise a new object
        """
        try:
            if takeable:
                self._values[label] = value
            else:
                self.index._engine.set_value(self._values, label, value)
            return self
        except KeyError:

            # set using a non-recursive method
            self.loc[label] = value
            return self

    def reset_index(self, level=None, drop=False, name=None, inplace=False):
        """
        Analogous to the :meth:`pandas.DataFrame.reset_index` function, see
        docstring there.

        Parameters
        ----------
        level : int, str, tuple, or list, default None
            Only remove the given levels from the index. Removes all levels by
            default
        drop : boolean, default False
            Do not try to insert index into dataframe columns
        name : object, default None
            The name of the column corresponding to the Series values
        inplace : boolean, default False
            Modify the Series in place (do not create a new object)

        Returns
        ----------
        resetted : DataFrame, or Series if drop == True
        """
        if drop:
            new_index = _default_index(len(self))
            if level is not None and isinstance(self.index, MultiIndex):
                if not isinstance(level, (tuple, list)):
                    level = [level]
                level = [self.index._get_level_number(lev) for lev in level]
                if len(level) < len(self.index.levels):
                    new_index = self.index.droplevel(level)

            if inplace:
                self.index = new_index
                # set name if it was passed, otherwise, keep the previous name
                self.name = name or self.name
            else:
                return self._constructor(self._values.copy(),
                                         index=new_index).__finalize__(self)
        elif inplace:
            raise TypeError('Cannot reset_index inplace on a Series '
                            'to create a DataFrame')
        else:
            df = self.to_frame(name)
            return df.reset_index(level=level, drop=drop)

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        buf = StringIO(u(""))
        width, height = get_terminal_size()
        max_rows = (height if get_option("display.max_rows") == 0 else
                    get_option("display.max_rows"))

        self.to_string(buf=buf, name=self.name, dtype=self.dtype,
                       max_rows=max_rows)
        result = buf.getvalue()

        return result

    def to_string(self, buf=None, na_rep='NaN', float_format=None, header=True,
                  index=True, length=False, dtype=False, name=False,
                  max_rows=None):
        """
        Render a string representation of the Series

        Parameters
        ----------
        buf : StringIO-like, optional
            buffer to write to
        na_rep : string, optional
            string representation of NAN to use, default 'NaN'
        float_format : one-parameter function, optional
            formatter function to apply to columns' elements if they are floats
            default None
        header: boolean, default True
            Add the Series header (index name)
        index : bool, optional
            Add index (row) labels, default True
        length : boolean, default False
            Add the Series length
        dtype : boolean, default False
            Add the Series dtype
        name : boolean, default False
            Add the Series name if not None
        max_rows : int, optional
            Maximum number of rows to show before truncating. If None, show
            all.

        Returns
        -------
        formatted : string (if not buffer passed)
        """

        the_repr = self._get_repr(float_format=float_format, na_rep=na_rep,
                                  header=header, index=index, length=length,
                                  dtype=dtype, name=name, max_rows=max_rows)

        # catch contract violations
        if not isinstance(the_repr, compat.text_type):
            raise AssertionError("result must be of type unicode, type"
                                 " of result is {0!r}"
                                 "".format(the_repr.__class__.__name__))

        if buf is None:
            return the_repr
        else:
            try:
                buf.write(the_repr)
            except AttributeError:
                with open(buf, 'w') as f:
                    f.write(the_repr)

    def _get_repr(self, name=False, header=True, index=True, length=True,
                  dtype=True, na_rep='NaN', float_format=None, max_rows=None):
        """

        Internal function, should always return unicode string
        """
        formatter = fmt.SeriesFormatter(self, name=name, length=length,
                                        header=header, index=index,
                                        dtype=dtype, na_rep=na_rep,
                                        float_format=float_format,
                                        max_rows=max_rows)
        result = formatter.to_string()

        # TODO: following check prob. not neces.
        if not isinstance(result, compat.text_type):
            raise AssertionError("result must be of type unicode, type"
                                 " of result is {0!r}"
                                 "".format(result.__class__.__name__))
        return result

    def __iter__(self):
        """ provide iteration over the values of the Series
        box values if necessary """
        if com.is_datetimelike(self):
            return (_maybe_box_datetimelike(x) for x in self._values)
        else:
            return iter(self._values)

    def iteritems(self):
        """
        Lazily iterate over (index, value) tuples
        """
        return zip(iter(self.index), iter(self))

    if compat.PY3:  # pragma: no cover
        items = iteritems

    # ----------------------------------------------------------------------
    # Misc public methods

    def keys(self):
        """Alias for index"""
        return self.index

    def tolist(self):
        """ Convert Series to a nested list """
        return list(self)

    def to_dict(self):
        """
        Convert Series to {label -> value} dict

        Returns
        -------
        value_dict : dict
        """
        return dict(compat.iteritems(self))

    def to_frame(self, name=None):
        """
        Convert Series to DataFrame

        Parameters
        ----------
        name : object, default None
            The passed name should substitute for the series name (if it has
            one).

        Returns
        -------
        data_frame : DataFrame
        """
        if name is None:
            df = self._constructor_expanddim(self)
        else:
            df = self._constructor_expanddim({name: self})

        return df

    def to_sparse(self, kind='block', fill_value=None):
        """
        Convert Series to SparseSeries

        Parameters
        ----------
        kind : {'block', 'integer'}
        fill_value : float, defaults to NaN (missing)

        Returns
        -------
        sp : SparseSeries
        """
        from pandas.core.sparse import SparseSeries
        return SparseSeries(self, kind=kind,
                            fill_value=fill_value).__finalize__(self)

    def _set_name(self, name, inplace=False):
        """
        Set the Series name.

        Parameters
        ----------
        name : str
        inplace : bool
            whether to modify `self` directly or return a copy
        """
        ser = self if inplace else self.copy()
        ser.name = name
        return ser

    # ----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck

    def count(self, level=None):
        """
        Return number of non-NA/null observations in the Series

        Parameters
        ----------
        level : int or level name, default None
            If the axis is a MultiIndex (hierarchical), count along a
            particular level, collapsing into a smaller Series

        Returns
        -------
        nobs : int or Series (if level specified)
        """
        from pandas.core.index import _get_na_value

        if level is None:
            return notnull(_values_from_object(self)).sum()

        if isinstance(level, compat.string_types):
            level = self.index._get_level_number(level)

        lev = self.index.levels[level]
        lab = np.array(self.index.labels[level], subok=False, copy=True)

        mask = lab == -1
        if mask.any():
            lab[mask] = cnt = len(lev)
            lev = lev.insert(cnt, _get_na_value(lev.dtype.type))

        obs = lab[notnull(self.values)]
        out = np.bincount(obs, minlength=len(lev) or None)
        return self._constructor(out, index=lev,
                                 dtype='int64').__finalize__(self)

    def mode(self):
        """Returns the mode(s) of the dataset.

        Empty if nothing occurs at least 2 times.  Always returns Series even
        if only one value.

        Parameters
        ----------
        sort : bool, default True
            If True, will lexicographically sort values, if False skips
            sorting. Result ordering when ``sort=False`` is not defined.

        Returns
        -------
        modes : Series (sorted)
        """
        # TODO: Add option for bins like value_counts()
        return algos.mode(self)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last',
                                                   False: 'first'})
    @Appender(base._shared_docs['drop_duplicates'] % _shared_doc_kwargs)
    def drop_duplicates(self, keep='first', inplace=False):
        return super(Series, self).drop_duplicates(keep=keep, inplace=inplace)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last',
                                                   False: 'first'})
    @Appender(base._shared_docs['duplicated'] % _shared_doc_kwargs)
    def duplicated(self, keep='first'):
        return super(Series, self).duplicated(keep=keep)

    def idxmin(self, axis=None, skipna=True, *args, **kwargs):
        """
        Index of first occurrence of minimum of values.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        idxmin : Index of minimum of values

        Notes
        -----
        This method is the Series version of ``ndarray.argmin``.

        See Also
        --------
        DataFrame.idxmin
        numpy.ndarray.argmin
        """
        skipna = nv.validate_argmin_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmin(_values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    def idxmax(self, axis=None, skipna=True, *args, **kwargs):
        """
        Index of first occurrence of maximum of values.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        idxmax : Index of maximum of values

        Notes
        -----
        This method is the Series version of ``ndarray.argmax``.

        See Also
        --------
        DataFrame.idxmax
        numpy.ndarray.argmax
        """
        skipna = nv.validate_argmax_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmax(_values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    # ndarray compat
    argmin = idxmin
    argmax = idxmax

    def round(self, decimals=0, *args, **kwargs):
        """
        Round each value in a Series to the given number of decimals.

        Parameters
        ----------
        decimals : int
            Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.

        Returns
        -------
        Series object

        See Also
        --------
        numpy.around
        DataFrame.round

        """
        nv.validate_round(args, kwargs)
        result = _values_from_object(self).round(decimals)
        result = self._constructor(result, index=self.index).__finalize__(self)

        return result

    def quantile(self, q=0.5, interpolation='linear'):
        """
        Return value at the given quantile, a la numpy.percentile.

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)
            0 <= q <= 1, the quantile(s) to compute
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            .. versionadded:: 0.18.0

            This optional parameter specifies the interpolation method to use,
            when the desired quantile lies between two data points `i` and `j`:

                * linear: `i + (j - i) * fraction`, where `fraction` is the
                  fractional part of the index surrounded by `i` and `j`.
                * lower: `i`.
                * higher: `j`.
                * nearest: `i` or `j` whichever is nearest.
                * midpoint: (`i` + `j`) / 2.

        Returns
        -------
        quantile : float or Series
            if ``q`` is an array, a Series will be returned where the
            index is ``q`` and the values are the quantiles.

        Examples
        --------
        >>> s = Series([1, 2, 3, 4])
        >>> s.quantile(.5)
        2.5
        >>> s.quantile([.25, .5, .75])
        0.25    1.75
        0.50    2.50
        0.75    3.25
        dtype: float64

        """

        self._check_percentile(q)

        if _np_version_under1p9:
            if interpolation != 'linear':
                raise ValueError("Interpolation methods other than linear "
                                 "are not supported in numpy < 1.9.")

        kwargs = dict()
        if not _np_version_under1p9:
            kwargs.update({'interpolation': interpolation})

        result = self._data.quantile(qs=q, **kwargs)

        if com.is_list_like(result):
            # explicitly use Float64Index to coerce empty result to float dtype
            index = Float64Index(q)
            return self._constructor(result, index=index, name=self.name)
        else:
            # scalar
            return result

    def corr(self, other, method='pearson', min_periods=None):
        """
        Compute correlation with `other` Series, excluding missing values

        Parameters
        ----------
        other : Series
        method : {'pearson', 'kendall', 'spearman'}
            * pearson : standard correlation coefficient
            * kendall : Kendall Tau correlation coefficient
            * spearman : Spearman rank correlation
        min_periods : int, optional
            Minimum number of observations needed to have a valid result


        Returns
        -------
        correlation : float
        """
        this, other = self.align(other, join='inner', copy=False)
        if len(this) == 0:
            return np.nan
        return nanops.nancorr(this.values, other.values, method=method,
                              min_periods=min_periods)

    def cov(self, other, min_periods=None):
        """
        Compute covariance with Series, excluding missing values

        Parameters
        ----------
        other : Series
        min_periods : int, optional
            Minimum number of observations needed to have a valid result

        Returns
        -------
        covariance : float

        Normalized by N-1 (unbiased estimator).
        """
        this, other = self.align(other, join='inner', copy=False)
        if len(this) == 0:
            return np.nan
        return nanops.nancov(this.values, other.values,
                             min_periods=min_periods)

    def diff(self, periods=1):
        """
        1st discrete difference of object

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for forming difference

        Returns
        -------
        diffed : Series
        """
        result = algos.diff(_values_from_object(self), periods)
        return self._constructor(result, index=self.index).__finalize__(self)

    def autocorr(self, lag=1):
        """
        Lag-N autocorrelation

        Parameters
        ----------
        lag : int, default 1
            Number of lags to apply before performing autocorrelation.

        Returns
        -------
        autocorr : float
        """
        return self.corr(self.shift(lag))

    def dot(self, other):
        """
        Matrix multiplication with DataFrame or inner-product with Series
        objects

        Parameters
        ----------
        other : Series or DataFrame

        Returns
        -------
        dot_product : scalar or Series
        """
        from pandas.core.frame import DataFrame
        if isinstance(other, (Series, DataFrame)):
            common = self.index.union(other.index)
            if (len(common) > len(self.index) or
                    len(common) > len(other.index)):
                raise ValueError('matrices are not aligned')

            left = self.reindex(index=common, copy=False)
            right = other.reindex(index=common, copy=False)
            lvals = left.values
            rvals = right.values
        else:
            left = self
            lvals = self.values
            rvals = np.asarray(other)
            if lvals.shape[0] != rvals.shape[0]:
                raise Exception('Dot product shape mismatch, %s vs %s' %
                                (lvals.shape, rvals.shape))

        if isinstance(other, DataFrame):
            return self._constructor(np.dot(lvals, rvals),
                                     index=other.columns).__finalize__(self)
        elif isinstance(other, Series):
            return np.dot(lvals, rvals)
        elif isinstance(rvals, np.ndarray):
            return np.dot(lvals, rvals)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    @Substitution(klass='Series', value='v')
    @Appender(base._shared_docs['searchsorted'])
    def searchsorted(self, v, side='left', sorter=None):
        if sorter is not None:
            sorter = com._ensure_platform_int(sorter)
        return self._values.searchsorted(Series(v)._values,
                                         side=side, sorter=sorter)

    # -------------------------------------------------------------------
    # Combination

    def append(self, to_append, verify_integrity=False):
        """
        Concatenate two or more Series.

        Parameters
        ----------
        to_append : Series or list/tuple of Series
        verify_integrity : boolean, default False
            If True, raise Exception on creating index with duplicates

        Returns
        -------
        appended : Series

        Examples
        --------
        >>> s1 = pd.Series([1, 2, 3])
        >>> s2 = pd.Series([4, 5, 6])
        >>> s3 = pd.Series([4, 5, 6], index=[3,4,5])
        >>> s1.append(s2)
        0    1
        1    2
        2    3
        0    4
        1    5
        2    6
        dtype: int64

        >>> s1.append(s3)
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        With `verify_integrity` set to True:

        >>> s1.append(s2, verify_integrity=True)
        ValueError: Indexes have overlapping values: [0, 1, 2]


        """
        from pandas.tools.merge import concat

        if isinstance(to_append, (list, tuple)):
            to_concat = [self] + to_append
        else:
            to_concat = [self, to_append]
        return concat(to_concat, ignore_index=False,
                      verify_integrity=verify_integrity)

    def _binop(self, other, func, level=None, fill_value=None):
        """
        Perform generic binary operation with optional fill value

        Parameters
        ----------
        other : Series
        func : binary operator
        fill_value : float or object
            Value to substitute for NA/null values. If both Series are NA in a
            location, the result will be NA regardless of the passed fill value
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level

        Returns
        -------
        combined : Series
        """
        if not isinstance(other, Series):
            raise AssertionError('Other operand must be Series')

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            this, other = self.align(other, level=level, join='outer',
                                     copy=False)
            new_index = this.index

        this_vals = this.values
        other_vals = other.values

        if fill_value is not None:
            this_mask = isnull(this_vals)
            other_mask = isnull(other_vals)
            this_vals = this_vals.copy()
            other_vals = other_vals.copy()

            # one but not both
            mask = this_mask ^ other_mask
            this_vals[this_mask & mask] = fill_value
            other_vals[other_mask & mask] = fill_value

        result = func(this_vals, other_vals)
        name = _maybe_match_name(self, other)
        result = self._constructor(result, index=new_index, name=name)
        result = result.__finalize__(self)
        if name is None:
            # When name is None, __finalize__ overwrites current name
            result.name = None
        return result

    def combine(self, other, func, fill_value=nan):
        """
        Perform elementwise binary operation on two Series using given function
        with optional fill value when an index is missing from one Series or
        the other

        Parameters
        ----------
        other : Series or scalar value
        func : function
        fill_value : scalar value

        Returns
        -------
        result : Series
        """
        if isinstance(other, Series):
            new_index = self.index.union(other.index)
            new_name = _maybe_match_name(self, other)
            new_values = np.empty(len(new_index), dtype=self.dtype)
            for i, idx in enumerate(new_index):
                lv = self.get(idx, fill_value)
                rv = other.get(idx, fill_value)
                new_values[i] = func(lv, rv)
        else:
            new_index = self.index
            new_values = func(self._values, other)
            new_name = self.name
        return self._constructor(new_values, index=new_index, name=new_name)

    def combine_first(self, other):
        """
        Combine Series values, choosing the calling Series's values
        first. Result index will be the union of the two indexes

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
        """
        new_index = self.index.union(other.index)
        this = self.reindex(new_index, copy=False)
        other = other.reindex(new_index, copy=False)
        # TODO: do we need name?
        name = _maybe_match_name(self, other)  # noqa
        rs_vals = com._where_compat(isnull(this), other._values, this._values)
        return self._constructor(rs_vals, index=new_index).__finalize__(self)

    def update(self, other):
        """
        Modify Series in place using non-NA values from passed
        Series. Aligns on index

        Parameters
        ----------
        other : Series
        """
        other = other.reindex_like(self)
        mask = notnull(other)

        self._data = self._data.putmask(mask=mask, new=other, inplace=True)
        self._maybe_update_cacher()

    # ----------------------------------------------------------------------
    # Reindexing, sorting

    @Appender(generic._shared_docs['sort_values'] % _shared_doc_kwargs)
    def sort_values(self, axis=0, ascending=True, inplace=False,
                    kind='quicksort', na_position='last'):

        axis = self._get_axis_number(axis)

        # GH 5856/5853
        if inplace and self._is_cached:
            raise ValueError("This Series is a view of some other array, to "
                             "sort in-place you must create a copy")

        def _try_kind_sort(arr):
            # easier to ask forgiveness than permission
            try:
                # if kind==mergesort, it can fail for object dtype
                return arr.argsort(kind=kind)
            except TypeError:
                # stable sort not available for object dtype
                # uses the argsort default quicksort
                return arr.argsort(kind='quicksort')

        arr = self._values
        sortedIdx = np.empty(len(self), dtype=np.int32)

        bad = isnull(arr)

        good = ~bad
        idx = _default_index(len(self))

        argsorted = _try_kind_sort(arr[good])

        if not ascending:
            argsorted = argsorted[::-1]

        if na_position == 'last':
            n = good.sum()
            sortedIdx[:n] = idx[good][argsorted]
            sortedIdx[n:] = idx[bad]
        elif na_position == 'first':
            n = bad.sum()
            sortedIdx[n:] = idx[good][argsorted]
            sortedIdx[:n] = idx[bad]
        else:
            raise ValueError('invalid na_position: {!r}'.format(na_position))

        result = self._constructor(arr[sortedIdx], index=self.index[sortedIdx])

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    @Appender(generic._shared_docs['sort_index'] % _shared_doc_kwargs)
    def sort_index(self, axis=0, level=None, ascending=True, inplace=False,
                   sort_remaining=True):

        axis = self._get_axis_number(axis)
        index = self.index
        if level is not None:
            new_index, indexer = index.sortlevel(level, ascending=ascending,
                                                 sort_remaining=sort_remaining)
        elif isinstance(index, MultiIndex):
            from pandas.core.groupby import _lexsort_indexer
            indexer = _lexsort_indexer(index.labels, orders=ascending)
            indexer = com._ensure_platform_int(indexer)
            new_index = index.take(indexer)
        else:
            new_index, indexer = index.sort_values(return_indexer=True,
                                                   ascending=ascending)

        new_values = self._values.take(indexer)
        result = self._constructor(new_values, index=new_index)

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    def sort(self, axis=0, ascending=True, kind='quicksort',
             na_position='last', inplace=True):
        """
        DEPRECATED: use :meth:`Series.sort_values(inplace=True)` for INPLACE
        sorting

        Sort values and index labels by value. This is an inplace sort by
        default.  Series.order is the equivalent but returns a new Series.

        Parameters
        ----------
        axis : int (can only be zero)
        ascending : boolean, default True
            Sort ascending. Passing False sorts descending
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm
        na_position : {'first', 'last'} (optional, default='last')
            'first' puts NaNs at the beginning
            'last' puts NaNs at the end
        inplace : boolean, default True
            Do operation in place.

        See Also
        --------
        Series.sort_values
        """
        warnings.warn("sort is deprecated, use sort_values(inplace=True) for "
                      "INPLACE sorting", FutureWarning, stacklevel=2)

        return self.sort_values(ascending=ascending, kind=kind,
                                na_position=na_position, inplace=inplace)

    def order(self, na_last=None, ascending=True, kind='quicksort',
              na_position='last', inplace=False):
        """
        DEPRECATED: use :meth:`Series.sort_values`

        Sorts Series object, by value, maintaining index-value link.
        This will return a new Series by default. Series.sort is the equivalent
        but as an inplace method.

        Parameters
        ----------
        na_last : boolean (optional, default=True)--DEPRECATED; use na_position
            Put NaN's at beginning or end
        ascending : boolean, default True
            Sort ascending. Passing False sorts descending
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm
        na_position : {'first', 'last'} (optional, default='last')
            'first' puts NaNs at the beginning
            'last' puts NaNs at the end
        inplace : boolean, default False
            Do operation in place.

        Returns
        -------
        y : Series

        See Also
        --------
        Series.sort_values
        """
        warnings.warn("order is deprecated, use sort_values(...)",
                      FutureWarning, stacklevel=2)

        return self.sort_values(ascending=ascending, kind=kind,
                                na_position=na_position, inplace=inplace)

    def argsort(self, axis=0, kind='quicksort', order=None):
        """
        Overrides ndarray.argsort. Argsorts the value, omitting NA/null values,
        and places the result in the same locations as the non-NA values

        Parameters
        ----------
        axis : int (can only be zero)
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm
        order : ignored

        Returns
        -------
        argsorted : Series, with -1 indicated where nan values are present

        See also
        --------
        numpy.ndarray.argsort
        """
        values = self._values
        mask = isnull(values)

        if mask.any():
            result = Series(-1, index=self.index, name=self.name,
                            dtype='int64')
            notmask = ~mask
            result[notmask] = np.argsort(values[notmask], kind=kind)
            return self._constructor(result,
                                     index=self.index).__finalize__(self)
        else:
            return self._constructor(
                np.argsort(values, kind=kind), index=self.index,
                dtype='int64').__finalize__(self)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last',
                                                   False: 'first'})
    def nlargest(self, n=5, keep='first'):
        """Return the largest `n` elements.

        Parameters
        ----------
        n : int
            Return this many descending sorted values
        keep : {'first', 'last', False}, default 'first'
            Where there are duplicate values:
            - ``first`` : take the first occurrence.
            - ``last`` : take the last occurrence.
        take_last : deprecated

        Returns
        -------
        top_n : Series
            The n largest values in the Series, in sorted order

        Notes
        -----
        Faster than ``.sort_values(ascending=False).head(n)`` for small `n`
        relative to the size of the ``Series`` object.

        See Also
        --------
        Series.nsmallest

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> s = pd.Series(np.random.randn(1e6))
        >>> s.nlargest(10)  # only sorts up to the N requested
        """
        return algos.select_n(self, n=n, keep=keep, method='nlargest')

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last',
                                                   False: 'first'})
    def nsmallest(self, n=5, keep='first'):
        """Return the smallest `n` elements.

        Parameters
        ----------
        n : int
            Return this many ascending sorted values
        keep : {'first', 'last', False}, default 'first'
            Where there are duplicate values:
            - ``first`` : take the first occurrence.
            - ``last`` : take the last occurrence.
        take_last : deprecated

        Returns
        -------
        bottom_n : Series
            The n smallest values in the Series, in sorted order

        Notes
        -----
        Faster than ``.sort_values().head(n)`` for small `n` relative to
        the size of the ``Series`` object.

        See Also
        --------
        Series.nlargest

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> s = pd.Series(np.random.randn(1e6))
        >>> s.nsmallest(10)  # only sorts up to the N requested
        """
        return algos.select_n(self, n=n, keep=keep, method='nsmallest')

    def sortlevel(self, level=0, ascending=True, sort_remaining=True):
        """
        Sort Series with MultiIndex by chosen level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

        Parameters
        ----------
        level : int or level name, default None
        ascending : bool, default True

        Returns
        -------
        sorted : Series

        See Also
        --------
        Series.sort_index(level=...)

        """
        return self.sort_index(level=level, ascending=ascending,
                               sort_remaining=sort_remaining)

    def swaplevel(self, i=-2, j=-1, copy=True):
        """
        Swap levels i and j in a MultiIndex

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : Series

        .. versionchanged:: 0.18.1

           The indexes ``i`` and ``j`` are now optional, and default to
           the two innermost levels of the index.

        """
        new_index = self.index.swaplevel(i, j)
        return self._constructor(self._values, index=new_index,
                                 copy=copy).__finalize__(self)

    def reorder_levels(self, order):
        """
        Rearrange index levels using input order. May not drop or duplicate
        levels

        Parameters
        ----------
        order: list of int representing new level order.
               (reference level by number or key)
        axis: where to reorder levels

        Returns
        -------
        type of caller (new object)
        """
        if not isinstance(self.index, MultiIndex):  # pragma: no cover
            raise Exception('Can only reorder levels on a hierarchical axis.')

        result = self.copy()
        result.index = result.index.reorder_levels(order)
        return result

    def unstack(self, level=-1, fill_value=None):
        """
        Unstack, a.k.a. pivot, Series with MultiIndex to produce DataFrame.
        The level involved will automatically get sorted.

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to unstack, can pass level name
        fill_value : replace NaN with this value if the unstack produces
            missing values

            .. versionadded: 0.18.0

        Examples
        --------
        >>> s
        one  a   1.
        one  b   2.
        two  a   3.
        two  b   4.

        >>> s.unstack(level=-1)
             a   b
        one  1.  2.
        two  3.  4.

        >>> s.unstack(level=0)
           one  two
        a  1.   2.
        b  3.   4.

        Returns
        -------
        unstacked : DataFrame
        """
        from pandas.core.reshape import unstack
        return unstack(self, level, fill_value)

    # ----------------------------------------------------------------------
    # function application

    def map(self, arg, na_action=None):
        """
        Map values of Series using input correspondence (which can be
        a dict, Series, or function)

        Parameters
        ----------
        arg : function, dict, or Series
        na_action : {None, 'ignore'}
            If 'ignore', propagate NA values

        Examples
        --------
        >>> x
        one   1
        two   2
        three 3

        >>> y
        1  foo
        2  bar
        3  baz

        >>> x.map(y)
        one   foo
        two   bar
        three baz

        Returns
        -------
        y : Series
            same index as caller
        """

        if is_extension_type(self.dtype):
            values = self._values
            if na_action is not None:
                raise NotImplementedError
            map_f = lambda values, f: values.map(f)
        else:
            values = self.asobject

            if na_action == 'ignore':
                def map_f(values, f):
                    return lib.map_infer_mask(values, f,
                                              isnull(values).view(np.uint8))
            else:
                map_f = lib.map_infer

        if isinstance(arg, (dict, Series)):
            if isinstance(arg, dict):
                arg = self._constructor(arg, index=arg.keys())

            indexer = arg.index.get_indexer(values)
            new_values = algos.take_1d(arg._values, indexer)
        else:
            new_values = map_f(values, arg)

        return self._constructor(new_values,
                                 index=self.index).__finalize__(self)

    def apply(self, func, convert_dtype=True, args=(), **kwds):
        """
        Invoke function on values of Series. Can be ufunc (a NumPy function
        that applies to the entire Series) or a Python function that only works
        on single values

        Parameters
        ----------
        func : function
        convert_dtype : boolean, default True
            Try to find better dtype for elementwise function results. If
            False, leave as dtype=object
        args : tuple
            Positional arguments to pass to function in addition to the value
        Additional keyword arguments will be passed as keywords to the function

        Returns
        -------
        y : Series or DataFrame if func returns a Series

        See also
        --------
        Series.map: For element-wise operations

        Examples
        --------

        Create a series with typical summer temperatures for each city.

        >>> import pandas as pd
        >>> import numpy as np
        >>> series = pd.Series([20, 21, 12], index=['London',
        ... 'New York','Helsinki'])
        London      20
        New York    21
        Helsinki    12
        dtype: int64

        Square the values by defining a function and passing it as an
        argument to ``apply()``.

        >>> def square(x):
        ...     return x**2
        >>> series.apply(square)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Square the values by passing an anonymous function as an
        argument to ``apply()``.

        >>> series.apply(lambda x: x**2)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Define a custom function that needs additional positional
        arguments and pass these additional arguments using the
        ``args`` keyword.

        >>> def subtract_custom_value(x, custom_value):
        ...     return x-custom_value

        >>> series.apply(subtract_custom_value, args=(5,))
        London      15
        New York    16
        Helsinki     7
        dtype: int64

        Define a custom function that takes keyword arguments
        and pass these arguments to ``apply``.

        >>> def add_custom_values(x, **kwargs):
        ...     for month in kwargs:
        ...         x+=kwargs[month]
        ...         return x

        >>> series.apply(add_custom_values, june=30, july=20, august=25)
        London      95
        New York    96
        Helsinki    87
        dtype: int64

        Use a function from the Numpy library.

        >>> series.apply(np.log)
        London      2.995732
        New York    3.044522
        Helsinki    2.484907
        dtype: float64


        """
        if len(self) == 0:
            return self._constructor(dtype=self.dtype,
                                     index=self.index).__finalize__(self)

        if kwds or args and not isinstance(func, np.ufunc):
            f = lambda x: func(x, *args, **kwds)
        else:
            f = func

        if isinstance(f, np.ufunc):
            return f(self)

        if is_extension_type(self.dtype):
            mapped = self._values.map(f)
        else:
            values = self.asobject
            mapped = lib.map_infer(values, f, convert=convert_dtype)

        if len(mapped) and isinstance(mapped[0], Series):
            from pandas.core.frame import DataFrame
            return DataFrame(mapped.tolist(), index=self.index)
        else:
            return self._constructor(mapped,
                                     index=self.index).__finalize__(self)

    def _reduce(self, op, name, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        """
        perform a reduction operation

        if we have an ndarray as a value, then simply perform the operation,
        otherwise delegate to the object

        """
        delegate = self._values
        if isinstance(delegate, np.ndarray):
            # Validate that 'axis' is consistent with Series's single axis.
            self._get_axis_number(axis)
            if numeric_only:
                raise NotImplementedError('Series.{0} does not implement '
                                          'numeric_only.'.format(name))
            return op(delegate, skipna=skipna, **kwds)

        return delegate._reduce(op=op, name=name, axis=axis, skipna=skipna,
                                numeric_only=numeric_only,
                                filter_type=filter_type, **kwds)

    def _reindex_indexer(self, new_index, indexer, copy):
        if indexer is None:
            if copy:
                return self.copy()
            return self

        # be subclass-friendly
        new_values = algos.take_1d(self.get_values(), indexer)
        return self._constructor(new_values, index=new_index)

    def _needs_reindex_multi(self, axes, method, level):
        """ check if we do need a multi reindex; this is for compat with
        higher dims
        """
        return False

    @Appender(generic._shared_docs['align'] % _shared_doc_kwargs)
    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0,
              broadcast_axis=None):
        return super(Series, self).align(other, join=join, axis=axis,
                                         level=level, copy=copy,
                                         fill_value=fill_value, method=method,
                                         limit=limit, fill_axis=fill_axis,
                                         broadcast_axis=broadcast_axis)

    @Appender(generic._shared_docs['rename'] % _shared_doc_kwargs)
    def rename(self, index=None, **kwargs):
        non_mapping = lib.isscalar(index) or (com.is_list_like(index) and
                                              not com.is_dict_like(index))
        if non_mapping:
            return self._set_name(index, inplace=kwargs.get('inplace'))
        return super(Series, self).rename(index=index, **kwargs)

    @Appender(generic._shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, index=None, **kwargs):
        return super(Series, self).reindex(index=index, **kwargs)

    @Appender(generic._shared_docs['fillna'] % _shared_doc_kwargs)
    def fillna(self, value=None, method=None, axis=None, inplace=False,
               limit=None, downcast=None, **kwargs):
        return super(Series, self).fillna(value=value, method=method,
                                          axis=axis, inplace=inplace,
                                          limit=limit, downcast=downcast,
                                          **kwargs)

    @Appender(generic._shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0):
        return super(Series, self).shift(periods=periods, freq=freq, axis=axis)

    def reindex_axis(self, labels, axis=0, **kwargs):
        """ for compatibility with higher dims """
        if axis != 0:
            raise ValueError("cannot reindex series on non-zero axis!")
        return self.reindex(index=labels, **kwargs)

    def memory_usage(self, index=True, deep=False):
        """Memory usage of the Series

        Parameters
        ----------
        index : bool
            Specifies whether to include memory usage of Series index
        deep : bool
            Introspect the data deeply, interrogate
            `object` dtypes for system-level memory consumption

        Returns
        -------
        scalar bytes of memory consumed

        Notes
        -----
        Memory usage does not include memory consumed by elements that
        are not components of the array if deep=False

        See Also
        --------
        numpy.ndarray.nbytes
        """
        v = super(Series, self).memory_usage(deep=deep)
        if index:
            v += self.index.memory_usage(deep=deep)
        return v

    def take(self, indices, axis=0, convert=True, is_copy=False, **kwargs):
        """
        return Series corresponding to requested indices

        Parameters
        ----------
        indices : list / array of ints
        convert : translate negative to positive indices (default)

        Returns
        -------
        taken : Series

        See also
        --------
        numpy.ndarray.take
        """
        nv.validate_take(tuple(), kwargs)

        # check/convert indicies here
        if convert:
            indices = maybe_convert_indices(indices, len(self._get_axis(axis)))

        indices = com._ensure_platform_int(indices)
        new_index = self.index.take(indices)
        new_values = self._values.take(indices)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def isin(self, values):
        """
        Return a boolean :class:`~pandas.Series` showing whether each element
        in the :class:`~pandas.Series` is exactly contained in the passed
        sequence of ``values``.

        Parameters
        ----------
        values : set or list-like
            The sequence of values to test. Passing in a single string will
            raise a ``TypeError``. Instead, turn a single string into a
            ``list`` of one element.

            .. versionadded:: 0.18.1

            Support for values as a set

        Returns
        -------
        isin : Series (bool dtype)

        Raises
        ------
        TypeError
          * If ``values`` is a string

        See Also
        --------
        pandas.DataFrame.isin

        Examples
        --------

        >>> s = pd.Series(list('abc'))
        >>> s.isin(['a', 'c', 'e'])
        0     True
        1    False
        2     True
        dtype: bool

        Passing a single string as ``s.isin('a')`` will raise an error. Use
        a list of one element instead:

        >>> s.isin(['a'])
        0     True
        1    False
        2    False
        dtype: bool

        """
        result = algos.isin(_values_from_object(self), values)
        return self._constructor(result, index=self.index).__finalize__(self)

    def between(self, left, right, inclusive=True):
        """
        Return boolean Series equivalent to left <= series <= right. NA values
        will be treated as False

        Parameters
        ----------
        left : scalar
            Left boundary
        right : scalar
            Right boundary

        Returns
        -------
        is_between : Series
        """
        if inclusive:
            lmask = self >= left
            rmask = self <= right
        else:
            lmask = self > left
            rmask = self < right

        return lmask & rmask

    @classmethod
    def from_csv(cls, path, sep=',', parse_dates=True, header=None,
                 index_col=0, encoding=None, infer_datetime_format=False):
        """
        Read CSV file (DISCOURAGED, please use :func:`pandas.read_csv`
        instead).

        It is preferable to use the more powerful :func:`pandas.read_csv`
        for most general purposes, but ``from_csv`` makes for an easy
        roundtrip to and from a file (the exact counterpart of
        ``to_csv``), especially with a time Series.

        This method only differs from :func:`pandas.read_csv` in some defaults:

        - `index_col` is ``0`` instead of ``None`` (take first column as index
          by default)
        - `header` is ``None`` instead of ``0`` (the first row is not used as
          the column names)
        - `parse_dates` is ``True`` instead of ``False`` (try parsing the index
          as datetime by default)

        With :func:`pandas.read_csv`, the option ``squeeze=True`` can be used
        to return a Series like ``from_csv``.

        Parameters
        ----------
        path : string file path or file handle / StringIO
        sep : string, default ','
            Field delimiter
        parse_dates : boolean, default True
            Parse dates. Different default from read_table
        header : int, default None
            Row to use as header (skip prior rows)
        index_col : int or sequence, default 0
            Column to use for index. If a sequence is given, a MultiIndex
            is used. Different default from read_table
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        infer_datetime_format: boolean, default False
            If True and `parse_dates` is True for a column, try to infer the
            datetime format based on the first datetime string. If the format
            can be inferred, there often will be a large parsing speed-up.

        See also
        --------
        pandas.read_csv

        Returns
        -------
        y : Series
        """
        from pandas.core.frame import DataFrame
        df = DataFrame.from_csv(path, header=header, index_col=index_col,
                                sep=sep, parse_dates=parse_dates,
                                encoding=encoding,
                                infer_datetime_format=infer_datetime_format)
        result = df.iloc[:, 0]
        if header is None:
            result.index.name = result.name = None

        return result

    def to_csv(self, path, index=True, sep=",", na_rep='', float_format=None,
               header=False, index_label=None, mode='w', nanRep=None,
               encoding=None, date_format=None, decimal='.'):
        """
        Write Series to a comma-separated values (csv) file

        Parameters
        ----------
        path : string file path or file handle / StringIO. If None is provided
            the result is returned as a string.
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        header : boolean, default False
            Write out series name
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
        mode : Python write mode, default 'w'
        sep : character, default ","
            Field delimiter for the output file.
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        date_format: string, default None
            Format string for datetime objects.
        decimal: string, default '.'
            Character recognized as decimal separator. E.g. use ',' for
            European data
        """
        from pandas.core.frame import DataFrame
        df = DataFrame(self)
        # result is only a string if no path provided, otherwise None
        result = df.to_csv(path, index=index, sep=sep, na_rep=na_rep,
                           float_format=float_format, header=header,
                           index_label=index_label, mode=mode, nanRep=nanRep,
                           encoding=encoding, date_format=date_format,
                           decimal=decimal)
        if path is None:
            return result

    def dropna(self, axis=0, inplace=False, **kwargs):
        """
        Return Series without null values

        Returns
        -------
        valid : Series
        inplace : boolean, default False
            Do operation in place.
        """
        kwargs.pop('how', None)
        if kwargs:
            raise TypeError('dropna() got an unexpected keyword '
                            'argument "{0}"'.format(list(kwargs.keys())[0]))

        axis = self._get_axis_number(axis or 0)

        if self._can_hold_na:
            result = remove_na(self)
            if inplace:
                self._update_inplace(result)
            else:
                return result
        else:
            if inplace:
                # do nothing
                pass
            else:
                return self.copy()

    valid = lambda self, inplace=False, **kwargs: self.dropna(inplace=inplace,
                                                              **kwargs)

    def first_valid_index(self):
        """
        Return label for first non-NA/null value
        """
        if len(self) == 0:
            return None

        mask = isnull(self._values)
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[i]

    def last_valid_index(self):
        """
        Return label for last non-NA/null value
        """
        if len(self) == 0:
            return None

        mask = isnull(self._values[::-1])
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[len(self) - i - 1]

    # ----------------------------------------------------------------------
    # Time series-oriented methods

    def asof(self, where):
        """
        Return last good (non-NaN) value in Series if value is NaN for
        requested date.

        If there is no good value, NaN is returned.

        Parameters
        ----------
        where : date or array of dates

        Notes
        -----
        Dates are assumed to be sorted

        Returns
        -------
        value or NaN
        """
        if isinstance(where, compat.string_types):
            where = datetools.to_datetime(where)

        values = self._values

        if not hasattr(where, '__iter__'):
            start = self.index[0]
            if isinstance(self.index, PeriodIndex):
                where = Period(where, freq=self.index.freq).ordinal
                start = start.ordinal

            if where < start:
                return np.nan
            loc = self.index.searchsorted(where, side='right')
            if loc > 0:
                loc -= 1
            while isnull(values[loc]) and loc > 0:
                loc -= 1
            return values[loc]

        if not isinstance(where, Index):
            where = Index(where)

        locs = self.index.asof_locs(where, notnull(values))
        new_values = algos.take_1d(values, locs)
        return self._constructor(new_values, index=where).__finalize__(self)

    def to_timestamp(self, freq=None, how='start', copy=True):
        """
        Cast to datetimeindex of timestamps, at *beginning* of period

        Parameters
        ----------
        freq : string, default frequency of PeriodIndex
            Desired frequency
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end

        Returns
        -------
        ts : Series with DatetimeIndex
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        new_index = self.index.to_timestamp(freq=freq, how=how)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def to_period(self, freq=None, copy=True):
        """
        Convert Series from DatetimeIndex to PeriodIndex with desired
        frequency (inferred from index if not passed)

        Parameters
        ----------
        freq : string, default

        Returns
        -------
        ts : Series with PeriodIndex
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        new_index = self.index.to_period(freq=freq)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    # -------------------------------------------------------------------------
    # Datetimelike delegation methods

    def _make_dt_accessor(self):
        try:
            return maybe_to_datetimelike(self)
        except Exception:
            raise AttributeError("Can only use .dt accessor with datetimelike "
                                 "values")

    dt = base.AccessorProperty(CombinedDatetimelikeProperties,
                               _make_dt_accessor)

    # -------------------------------------------------------------------------
    # Categorical methods

    def _make_cat_accessor(self):
        if not is_categorical_dtype(self.dtype):
            raise AttributeError("Can only use .cat accessor with a "
                                 "'category' dtype")
        return CategoricalAccessor(self.values, self.index)

    cat = base.AccessorProperty(CategoricalAccessor, _make_cat_accessor)

    def _dir_deletions(self):
        return self._accessors

    def _dir_additions(self):
        rv = set()
        for accessor in self._accessors:
            try:
                getattr(self, accessor)
                rv.add(accessor)
            except AttributeError:
                pass
        return rv


Series._setup_axes(['index'], info_axis=0, stat_axis=0, aliases={'rows': 0})
Series._add_numeric_operations()
Series._add_series_only_operations()
Series._add_series_or_dataframe_operations()
_INDEX_TYPES = ndarray, Index, list, tuple

# -----------------------------------------------------------------------------
# Supplementary functions


def remove_na(series):
    """
    Return series containing only true/non-NaN values, possibly empty.
    """
    return series[notnull(_values_from_object(series))]


def _sanitize_index(data, index, copy=False):
    """ sanitize an index type to return an ndarray of the underlying, pass
    thru a non-Index
    """

    if index is None:
        return data

    if len(data) != len(index):
        raise ValueError('Length of values does not match length of ' 'index')

    if isinstance(data, PeriodIndex):
        data = data.asobject
    elif isinstance(data, DatetimeIndex):
        data = data._to_embed(keep_tz=True)
        if copy:
            data = data.copy()
    elif isinstance(data, np.ndarray):

        # coerce datetimelike types
        if data.dtype.kind in ['M', 'm']:
            data = _sanitize_array(data, index, copy=copy)

    return data


def _sanitize_array(data, index, dtype=None, copy=False,
                    raise_cast_failure=False):
    """ sanitize input data to an ndarray, copy if specified, coerce to the
    dtype if specified
    """

    if dtype is not None:
        dtype = _coerce_to_dtype(dtype)

    if isinstance(data, ma.MaskedArray):
        mask = ma.getmaskarray(data)
        if mask.any():
            data, fill_value = _maybe_upcast(data, copy=True)
            data[mask] = fill_value
        else:
            data = data.copy()

    def _try_cast(arr, take_fast_path):

        # perf shortcut as this is the most common case
        if take_fast_path:
            if _possibly_castable(arr) and not copy and dtype is None:
                return arr

        try:
            subarr = _possibly_cast_to_datetime(arr, dtype)
            if not is_extension_type(subarr):
                subarr = np.array(subarr, dtype=dtype, copy=copy)
        except (ValueError, TypeError):
            if is_categorical_dtype(dtype):
                subarr = Categorical(arr)
            elif dtype is not None and raise_cast_failure:
                raise
            else:
                subarr = np.array(arr, dtype=object, copy=copy)
        return subarr

    # GH #846
    if isinstance(data, (np.ndarray, Index, Series)):

        if dtype is not None:
            subarr = np.array(data, copy=False)

            # possibility of nan -> garbage
            if com.is_float_dtype(data.dtype) and com.is_integer_dtype(dtype):
                if not isnull(data).any():
                    subarr = _try_cast(data, True)
                elif copy:
                    subarr = data.copy()
            else:
                subarr = _try_cast(data, True)
        elif isinstance(data, Index):
            # don't coerce Index types
            # e.g. indexes can have different conversions (so don't fast path
            # them)
            # GH 6140
            subarr = _sanitize_index(data, index, copy=True)
        else:
            subarr = _try_cast(data, True)

        if copy:
            subarr = data.copy()

    elif isinstance(data, Categorical):
        subarr = data

        if copy:
            subarr = data.copy()
        return subarr

    elif isinstance(data, list) and len(data) > 0:
        if dtype is not None:
            try:
                subarr = _try_cast(data, False)
            except Exception:
                if raise_cast_failure:  # pragma: no cover
                    raise
                subarr = np.array(data, dtype=object, copy=copy)
                subarr = lib.maybe_convert_objects(subarr)

        else:
            subarr = _possibly_convert_platform(data)

        subarr = _possibly_cast_to_datetime(subarr, dtype)

    else:
        subarr = _try_cast(data, False)

    def create_from_value(value, index, dtype):
        # return a new empty value suitable for the dtype

        if is_datetimetz(dtype):
            subarr = DatetimeIndex([value] * len(index), dtype=dtype)
        elif is_categorical_dtype(dtype):
            subarr = Categorical([value] * len(index))
        else:
            if not isinstance(dtype, (np.dtype, type(np.dtype))):
                dtype = dtype.dtype
            subarr = np.empty(len(index), dtype=dtype)
            subarr.fill(value)

        return subarr

    # scalar like
    if subarr.ndim == 0:
        if isinstance(data, list):  # pragma: no cover
            subarr = np.array(data, dtype=object)
        elif index is not None:
            value = data

            # figure out the dtype from the value (upcast if necessary)
            if dtype is None:
                dtype, value = _infer_dtype_from_scalar(value)
            else:
                # need to possibly convert the value here
                value = _possibly_cast_to_datetime(value, dtype)

            subarr = create_from_value(value, index, dtype)

        else:
            return subarr.item()

    # the result that we want
    elif subarr.ndim == 1:
        if index is not None:

            # a 1-element ndarray
            if len(subarr) != len(index) and len(subarr) == 1:
                subarr = create_from_value(subarr[0], index,
                                           subarr.dtype)

    elif subarr.ndim > 1:
        if isinstance(data, np.ndarray):
            raise Exception('Data must be 1-dimensional')
        else:
            subarr = _asarray_tuplesafe(data, dtype=dtype)

    # This is to prevent mixed-type Series getting all casted to
    # NumPy string type, e.g. NaN --> '-1#IND'.
    if issubclass(subarr.dtype.type, compat.string_types):
        subarr = np.array(data, dtype=object, copy=copy)

    return subarr


# backwards compatiblity
class TimeSeries(Series):
    def __init__(self, *args, **kwargs):
        # deprecation TimeSeries, #10890
        warnings.warn("TimeSeries is deprecated. Please use Series",
                      FutureWarning, stacklevel=2)

        super(TimeSeries, self).__init__(*args, **kwargs)

# ----------------------------------------------------------------------
# Add plotting methods to Series

import pandas.tools.plotting as _gfx  # noqa

Series.plot = base.AccessorProperty(_gfx.SeriesPlotMethods,
                                    _gfx.SeriesPlotMethods)
Series.hist = _gfx.hist_series

# Add arithmetic!
ops.add_flex_arithmetic_methods(Series, **ops.series_flex_funcs)
ops.add_special_arithmetic_methods(Series, **ops.series_special_funcs)
