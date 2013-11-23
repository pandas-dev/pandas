"""
Data structure for 1-dimensional cross-sectional and time series data
"""
from __future__ import division

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import operator
import types
import warnings

from numpy import nan, ndarray
import numpy as np
import numpy.ma as ma

from pandas.core.common import (isnull, notnull, _is_bool_indexer,
                                _default_index, _maybe_promote, _maybe_upcast,
                                _asarray_tuplesafe, is_integer_dtype,
                                _NS_DTYPE, _TD_DTYPE,
                                _infer_dtype_from_scalar, is_list_like,
                                _values_from_object,
                                _possibly_cast_to_datetime, _possibly_castable,
                                _possibly_convert_platform,
                                ABCSparseArray, _maybe_match_name,
                                _ensure_object, SettingWithCopyError)

from pandas.core.index import (Index, MultiIndex, InvalidIndexError,
                               _ensure_index, _handle_legacy_indexes)
from pandas.core.indexing import (
    _check_bool_indexer, _check_slice_bounds,
    _is_index_slice, _maybe_convert_indices)
from pandas.core import generic
from pandas.core.internals import SingleBlockManager
from pandas.core.categorical import Categorical
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex, Period
from pandas import compat
from pandas.util.terminal import get_terminal_size
from pandas.compat import zip, lzip, u, OrderedDict

import pandas.core.array as pa
import pandas.core.ops as ops

import pandas.core.common as com
import pandas.core.datetools as datetools
import pandas.core.format as fmt
import pandas.core.nanops as nanops
from pandas.util.decorators import Appender, Substitution, cache_readonly

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.index as _index

from pandas.compat.scipy import scoreatpercentile as _quantile
from pandas.core.config import get_option

__all__ = ['Series']

_shared_doc_kwargs = dict(
    axes='index',
    klass='Series',
    axes_single_arg="{0,'index'}"
)


def _coerce_method(converter):
    """ install the scalar coercion methods """

    def wrapper(self):
        if len(self) == 1:
            return converter(self.iloc[0])
        raise TypeError(
            "cannot convert the series to {0}".format(str(converter)))
    return wrapper


def _unbox(func):
    @Appender(func.__doc__)
    def f(self, *args, **kwargs):
        result = func(self.values, *args, **kwargs)
        if isinstance(result, (pa.Array, Series)) and result.ndim == 0:
            # return NumPy type
            return result.dtype.type(result.item())
        else:  # pragma: no cover
            return result
    f.__name__ = func.__name__
    return f

#----------------------------------------------------------------------
# Series class


class Series(generic.NDFrame):

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
        np.arange(len(data)) if not provided. If both a dict and index
        sequence are used, the index will override the keys found in the
        dict.
    dtype : numpy.dtype or None
        If None, dtype will be inferred
    copy : boolean, default False
        Copy input data
    """
    _metadata = ['name']

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
                raise NotImplementedError
            elif isinstance(data, Index):
                # need to copy to avoid aliasing issues
                if name is None:
                    name = data.name
                data = data.values
                copy = True
            elif isinstance(data, pa.Array):
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
                        index = Index(sorted(data))
                try:
                    if isinstance(index, DatetimeIndex):
                        # coerce back to datetime objects for lookup
                        data = lib.fast_multiget(data, index.astype('O'),
                                                 default=pa.NA)
                    elif isinstance(index, PeriodIndex):
                        data = [data.get(i, nan) for i in index]
                    else:
                        data = lib.fast_multiget(data, index.values,
                                                 default=pa.NA)
                except TypeError:
                    data = [data.get(i, nan) for i in index]

            elif isinstance(data, SingleBlockManager):
                if index is None:
                    index = data.index
                else:
                    data = data.reindex(index, copy=copy)
            elif isinstance(data, Categorical):
                if name is None:
                    name = data.name
                data = np.asarray(data)
            elif isinstance(data, types.GeneratorType):
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
                    data = data.astype(dtype, raise_on_error=False)
                elif copy:
                    data = data.copy()
            else:
                data = _sanitize_array(data, index, dtype, copy,
                                       raise_cast_failure=True)

                data = SingleBlockManager(data, index, fastpath=True)

        generic.NDFrame.__init__(self, data, fastpath=True)

        object.__setattr__(self, 'name', name)
        self._set_axis(0, index, fastpath=True)

    @classmethod
    def from_array(cls, arr, index=None, name=None, copy=False,
                   fastpath=False):
        # return a sparse series here
        if isinstance(arr, ABCSparseArray):
            from pandas.sparse.series import SparseSeries
            cls = SparseSeries

        return cls(arr, index=index, name=name, copy=copy, fastpath=fastpath)

    @property
    def _constructor(self):
        return Series

    # types
    @property
    def _can_hold_na(self):
        return self._data._can_hold_na

    @property
    def is_time_series(self):
        return self._subtyp in ['time_series', 'sparse_time_series']

    _index = None

    def _set_axis(self, axis, labels, fastpath=False):
        """ override generic, we want to set the _typ here """

        if not fastpath:
            labels = _ensure_index(labels)

        is_all_dates = labels.is_all_dates
        if is_all_dates:
            from pandas.tseries.index import DatetimeIndex
            from pandas.tseries.period import PeriodIndex
            if not isinstance(labels, (DatetimeIndex, PeriodIndex)):
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

    # ndarray compatibility
    def item(self):
        return self.values.item()

    @property
    def data(self):
        return self.values.data

    @property
    def strides(self):
        return self.values.strides

    @property
    def size(self):
        return self.values.size

    @property
    def flags(self):
        return self.values.flags

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def dtypes(self):
        """ for compat """
        return self._data.dtype

    @property
    def ftype(self):
        return self._data.ftype

    @property
    def ftypes(self):
        """ for compat """
        return self._data.ftype

    @property
    def shape(self):
        return self._data.shape

    @property
    def ndim(self):
        return 1

    @property
    def base(self):
        return self.values.base

    def ravel(self, order='C'):
        return self.values.ravel(order=order)

    def transpose(self):
        """ support for compatiblity """
        return self

    T = property(transpose)

    def nonzero(self):
        """ numpy like, returns same as nonzero """
        return self.values.nonzero()

    def put(self, *args, **kwargs):
        self.values.put(*args, **kwargs)

    def __len__(self):
        return len(self._data)

    def view(self, dtype=None):
        return self._constructor(self.values.view(dtype),
                                 index=self.index).__finalize__(self)

    def __array__(self, result=None):
        """ the array interface, return my values """
        return self.values

    def __array_wrap__(self, result):
        """
        Gets called prior to a ufunc (and after)
        """
        return self._constructor(result, index=self.index,
                                 copy=False).__finalize__(self)

    def __contains__(self, key):
        return key in self.index

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

    # we are preserving name here
    def __getstate__(self):
        return dict(_data=self._data, name=self.name)

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
            index = _handle_legacy_indexes([index])[0]

            # recreate
            self._data = SingleBlockManager(data, index, fastpath=True)
            self.index = index
            self.name = name

        else:
            raise Exception("cannot unpickle legacy formats -> [%s]" % state)

    # indexers
    @property
    def axes(self):
        return [self.index]

    def _maybe_box(self, values):
        """ genericically box the values """

        if isinstance(values, self.__class__):
            return values
        elif not hasattr(values, '__iter__'):
            v = lib.infer_dtype([values])
            if v == 'datetime':
                return lib.Timestamp(v)
            return values

        v = lib.infer_dtype(values)
        if v == 'datetime':
            return lib.map_infer(values, lib.Timestamp)

        if isinstance(values, np.ndarray):
            return self.__class__(values)

        return values

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
            return _index.get_value_at(self.values, i)
        except IndexError:
            raise
        except:
            if isinstance(i, slice):
                indexer = self.index._convert_slice_indexer(i, typ='iloc')
                return self._get_values(indexer)
            else:
                label = self.index[i]
                if isinstance(label, Index):
                    i = _maybe_convert_indices(i, len(self))
                    return self.reindex(i, takeable=True)
                else:
                    return _index.get_value_at(self, i)

    @property
    def _is_mixed_type(self):
        return False

    def _slice(self, slobj, axis=0, raise_on_error=False, typ=None):
        if raise_on_error:
            _check_slice_bounds(slobj, self.values)
        slobj = self.index._convert_slice_indexer(slobj, typ=typ or 'getitem')
        return self._constructor(self.values[slobj],
                                 index=self.index[slobj]).__finalize__(self)

    def __getitem__(self, key):
        try:
            result = self.index.get_value(self, key)
            if isinstance(result, np.ndarray):
                return self._constructor(result,index=[key]*len(result)).__finalize__(self)
            return result
        except InvalidIndexError:
            pass
        except (KeyError, ValueError):
            if isinstance(key, tuple) and isinstance(self.index, MultiIndex):
                # kludge
                pass
            elif key is Ellipsis:
                return self
            elif _is_bool_indexer(key):
                pass
            else:

                # we can try to coerce the indexer (or this will raise)
                new_key = self.index._convert_scalar_indexer(key)
                if type(new_key) != type(key):
                    return self.__getitem__(new_key)
                raise

        except Exception:
            raise

        if com.is_iterator(key):
            key = list(key)

        if _is_bool_indexer(key):
            key = _check_bool_indexer(self.index, key)

        return self._get_with(key)

    def _get_with(self, key):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, typ='getitem')
            return self._get_values(indexer)
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
            if not isinstance(key, (list, pa.Array, Series)):
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
        return self._constructor(self.values[indexer],
                                 index=new_index).__finalize__(self)

    def _get_values(self, indexer):
        try:
            return self._constructor(self._data.get_slice(indexer),
                                     fastpath=True).__finalize__(self)
        except Exception:
            return self.values[indexer]

    def __setitem__(self, key, value):
        try:
            self._set_with_engine(key, value)
            return
        except (SettingWithCopyError):
            raise
        except (KeyError, ValueError):
            values = self.values
            if (com.is_integer(key)
                    and not self.index.inferred_type == 'integer'):

                values[key] = value
                return
            elif key is Ellipsis:
                self[:] = value
                return
            elif _is_bool_indexer(key):
                pass
            elif com.is_timedelta64_dtype(self.dtype):
                # reassign a null value to iNaT
                if isnull(value):
                    value = tslib.iNaT

                    try:
                        self.index._engine.set_value(self.values, key, value)
                        return
                    except (TypeError):
                        pass

            self.loc[key] = value
            return

        except TypeError as e:
            if isinstance(key, tuple) and not isinstance(self.index,
                                                         MultiIndex):
                raise ValueError("Can only tuple-index with a MultiIndex")

            # python 3 type errors should be raised
            if 'unorderable' in str(e):  # pragma: no cover
                raise IndexError(key)

        if _is_bool_indexer(key):
            key = _check_bool_indexer(self.index, key)
            try:
                self.where(~key, value, inplace=True)
                return
            except (InvalidIndexError):
                pass

        self._set_with(key, value)

    def _set_with_engine(self, key, value):
        values = self.values
        try:
            self.index._engine.set_value(values, key, value)
            self._check_setitem_copy()
            return
        except KeyError:
            values[self.index.get_loc(key)] = value
            return

    def _set_with(self, key, value):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, typ='getitem')
            return self._set_values(indexer, value)
        else:
            if isinstance(key, tuple):
                try:
                    self._set_values(key, value)
                except Exception:
                    pass

            if not isinstance(key, (list, Series, pa.Array, Series)):
                key = list(key)

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
            raise ValueError('%s not contained in the index'
                             % str(key[mask]))
        self._set_values(indexer, value)

    def _set_values(self, key, value):
        if isinstance(key, Series):
            key = key.values
        self._data = self._data.setitem(key, value)

    # help out SparseSeries
    _get_val_at = ndarray.__getitem__

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        if j < 0:
            j = 0
        slobj = slice(i, j)
        return self._slice(slobj)

    def __setslice__(self, i, j, value):
        """Set slice equal to given value(s)"""
        if i < 0:
            i = 0
        if j < 0:
            j = 0
        slobj = slice(i, j)
        return self.__setitem__(slobj, value)

    def repeat(self, reps):
        """
        See ndarray.repeat
        """
        new_index = self.index.repeat(reps)
        new_values = self.values.repeat(reps)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def reshape(self, *args, **kwargs):
        """
        See numpy.ndarray.reshape
        """
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            shape = args[0]
        else:
            shape = args

        if tuple(shape) == self.shape:
            # XXX ignoring the "order" keyword.
            return self

        return self.values.reshape(shape, **kwargs)

    def get(self, label, default=None):
        """
        Returns value occupying requested label, default to specified
        missing value if not present. Analogous to dict.get

        Parameters
        ----------
        label : object
            Label value looking for
        default : object, optional
            Value to return if label not in index

        Returns
        -------
        y : scalar
        """
        try:
            return self.get_value(label)
        except KeyError:
            return default

    iget_value = _ixs
    iget = _ixs
    irow = _ixs

    def get_value(self, label):
        """
        Quickly retrieve single value at passed index label

        Parameters
        ----------
        index : label

        Returns
        -------
        value : scalar value
        """
        return self.index.get_value(self.values, label)

    def set_value(self, label, value):
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

        Returns
        -------
        series : Series
            If label is contained, will be reference to calling Series,
            otherwise a new object
        """
        try:
            self.index._engine.set_value(self.values, label, value)
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
            new_index = pa.arange(len(self))
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
                return self._constructor(self.values.copy(),
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
        width, height = get_terminal_size()
        max_rows = (height if get_option("display.max_rows") == 0
                    else get_option("display.max_rows"))
        if len(self.index) > (max_rows or 1000):
            result = self._tidy_repr(min(30, max_rows - 4))
        elif len(self.index) > 0:
            result = self._get_repr(print_header=True,
                                    length=len(self) > 50,
                                    name=True,
                                    dtype=True)
        elif self.name is None:
            result = u('Series([], dtype: %s)') % (self.dtype)
        else:
            result = u('Series([], name: %s, dtype: %s)') % (self.name,
                                                             self.dtype)
        return result

    def _tidy_repr(self, max_vals=20):
        """

        Internal function, should always return unicode string
        """
        num = max_vals // 2
        head = self.iloc[:num]._get_repr(print_header=True, length=False,
                                         dtype=False, name=False)
        tail = self.iloc[-(max_vals - num):]._get_repr(print_header=False,
                                                       length=False,
                                                       name=False,
                                                       dtype=False)
        result = head + '\n...\n' + tail
        result = '%s\n%s' % (result, self._repr_footer())

        return compat.text_type(result)

    def _repr_footer(self):

        # time series
        if self.is_time_series:
            if self.index.freq is not None:
                freqstr = u('Freq: %s, ') % self.index.freqstr
            else:
                freqstr = u('')

            namestr = u("Name: %s, ") % com.pprint_thing(
                self.name) if self.name is not None else ""
            return u('%s%sLength: %d') % (freqstr, namestr, len(self))

        # reg series
        namestr = u("Name: %s, ") % com.pprint_thing(
            self.name) if self.name is not None else ""
        return u('%sLength: %d, dtype: %s') % (namestr,
                                               len(self),
                                               str(self.dtype.name))

    def to_string(self, buf=None, na_rep='NaN', float_format=None,
                  nanRep=None, length=False, dtype=False, name=False):
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
        length : boolean, default False
            Add the Series length
        dtype : boolean, default False
            Add the Series dtype
        name : boolean, default False
            Add the Series name (which may be None)

        Returns
        -------
        formatted : string (if not buffer passed)
        """

        if nanRep is not None:  # pragma: no cover
            warnings.warn("nanRep is deprecated, use na_rep", FutureWarning)
            na_rep = nanRep

        the_repr = self._get_repr(float_format=float_format, na_rep=na_rep,
                                  length=length, dtype=dtype, name=name)

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

    def _get_repr(
        self, name=False, print_header=False, length=True, dtype=True,
            na_rep='NaN', float_format=None):
        """

        Internal function, should always return unicode string
        """

        formatter = fmt.SeriesFormatter(self, name=name, header=print_header,
                                        length=length, dtype=dtype,
                                        na_rep=na_rep,
                                        float_format=float_format)
        result = formatter.to_string()

        # TODO: following check prob. not neces.
        if not isinstance(result, compat.text_type):
            raise AssertionError("result must be of type unicode, type"
                                 " of result is {0!r}"
                                 "".format(result.__class__.__name__))
        return result

    def __iter__(self):
        if np.issubdtype(self.dtype, np.datetime64):
            return (lib.Timestamp(x) for x in self.values)
        else:
            return iter(self.values)

    def iteritems(self):
        """
        Lazily iterate over (index, value) tuples
        """
        return lzip(iter(self.index), iter(self))

    if compat.PY3:  # pragma: no cover
        items = iteritems

    # inversion
    def __neg__(self):
        arr = operator.neg(self.values)
        return self._constructor(arr, self.index).__finalize__(self)

    def __invert__(self):
        arr = operator.inv(self.values)
        return self._constructor(arr, self.index).__finalize__(self)

    #----------------------------------------------------------------------
    # unbox reductions

    all = _unbox(pa.Array.all)
    any = _unbox(pa.Array.any)

    #----------------------------------------------------------------------
    # Misc public methods

    def keys(self):
        "Alias for index"
        return self.index

    @property
    def values(self):
        """
        Return Series as ndarray

        Returns
        -------
        arr : numpy.ndarray
        """
        return self._data.values

    def get_values(self):
        """ same as values (but handles sparseness conversions); is a view """
        return self._data.values

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
        from pandas.core.frame import DataFrame
        if name is None:
            df = DataFrame(self)
        else:
            df = DataFrame({name: self})

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

    #----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck

    def count(self, level=None):
        """
        Return number of non-NA/null observations in the Series

        Parameters
        ----------
        level : int, default None
            If the axis is a MultiIndex (hierarchical), count along a
            particular level, collapsing into a smaller Series

        Returns
        -------
        nobs : int or Series (if level specified)
        """
        if level is not None:
            mask = notnull(self.values)

            if isinstance(level, compat.string_types):
                level = self.index._get_level_number(level)

            level_index = self.index.levels[level]

            if len(self) == 0:
                return self._constructor(0, index=level_index)\
                           .__finalize__(self)

            # call cython function
            max_bin = len(level_index)
            labels = com._ensure_int64(self.index.labels[level])
            counts = lib.count_level_1d(mask.view(pa.uint8),
                                        labels, max_bin)
            return self._constructor(counts,
                                     index=level_index).__finalize__(self)

        return notnull(_values_from_object(self)).sum()

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None):
        """
        Returns Series containing counts of unique values. The resulting Series
        will be in descending order so that the first element is the most
        frequently-occurring element. Excludes NA values

        Parameters
        ----------
        normalize : boolean, default False
            If True then the Series returned will contain the relative
            frequencies of the unique values.
        sort : boolean, default True
            Sort by values
        ascending : boolean, default False
            Sort in ascending order
        bins : integer, optional
            Rather than count values, group them into half-open bins,
            a convenience for pd.cut, only works with numeric data

        Returns
        -------
        counts : Series
        """
        from pandas.core.algorithms import value_counts
        return value_counts(self.values, sort=sort, ascending=ascending,
                            normalize=normalize, bins=bins)

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
        from pandas.core.algorithms import mode
        return mode(self)

    def unique(self):
        """
        Return array of unique values in the Series. Significantly faster than
        numpy.unique

        Returns
        -------
        uniques : ndarray
        """
        return nanops.unique1d(self.values)

    def nunique(self):
        """
        Return count of unique elements in the Series

        Returns
        -------
        nunique : int
        """
        return len(self.value_counts())

    def drop_duplicates(self, take_last=False, inplace=False):
        """
        Return Series with duplicate values removed

        Parameters
        ----------
        take_last : boolean, default False
            Take the last observed index in a group. Default first
        inplace : boolean, default False
            If True, performs operation inplace and returns None.

        Returns
        -------
        deduplicated : Series
        """
        duplicated = self.duplicated(take_last=take_last)
        result = self[-duplicated]
        if inplace:
            return self._update_inplace(result)
        else:
            return result

    def duplicated(self, take_last=False):
        """
        Return boolean Series denoting duplicate values

        Parameters
        ----------
        take_last : boolean, default False
            Take the last observed index in a group. Default first

        Returns
        -------
        duplicated : Series
        """
        keys = _ensure_object(self.values)
        duplicated = lib.duplicated(keys, take_last=take_last)
        return self._constructor(duplicated,
                                 index=self.index).__finalize__(self)

    def idxmin(self, axis=None, out=None, skipna=True):
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
        """
        i = nanops.nanargmin(_values_from_object(self), skipna=skipna)
        if i == -1:
            return pa.NA
        return self.index[i]

    def idxmax(self, axis=None, out=None, skipna=True):
        """
        Index of first occurrence of maximum of values.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        idxmax : Index of minimum of values

        Notes
        -----
        This method is the Series version of ``ndarray.argmax``.

        See Also
        --------
        DataFrame.idxmax
        """
        i = nanops.nanargmax(_values_from_object(self), skipna=skipna)
        if i == -1:
            return pa.NA
        return self.index[i]

    # ndarray compat
    argmin = idxmin
    argmax = idxmax

    @Appender(pa.Array.round.__doc__)
    def round(self, decimals=0, out=None):
        """

        """
        result = _values_from_object(self).round(decimals, out=out)
        if out is None:
            result = self._constructor(result,
                                       index=self.index).__finalize__(self)

        return result

    def quantile(self, q=0.5):
        """
        Return value at the given quantile, a la scoreatpercentile in
        scipy.stats

        Parameters
        ----------
        q : quantile
            0 <= q <= 1

        Returns
        -------
        quantile : float
        """
        valid_values = self.dropna().values
        if len(valid_values) == 0:
            return pa.NA
        result = _quantile(valid_values, q * 100)
        if result.dtype == _TD_DTYPE:
            from pandas.tseries.timedeltas import to_timedelta
            return to_timedelta(result)

        return result

    def ptp(self, axis=None, out=None):
        return _values_from_object(self).ptp(axis, out)

    def describe(self, percentile_width=50):
        """
        Generate various summary statistics of Series, excluding NaN
        values. These include: count, mean, std, min, max, and
        lower%/50%/upper% percentiles

        Parameters
        ----------
        percentile_width : float, optional
            width of the desired uncertainty interval, default is 50,
            which corresponds to lower=25, upper=75

        Returns
        -------
        desc : Series
        """
        from pandas.compat import Counter

        if self.dtype == object:
            names = ['count', 'unique']
            objcounts = Counter(self.dropna().values)
            data = [self.count(), len(objcounts)]
            if data[1] > 0:
                names += ['top', 'freq']
                top, freq = objcounts.most_common(1)[0]
                data += [top, freq]

        elif issubclass(self.dtype.type, np.datetime64):
            names = ['count', 'unique']
            asint = self.dropna().values.view('i8')
            objcounts = Counter(asint)
            data = [self.count(), len(objcounts)]
            if data[1] > 0:
                top, freq = objcounts.most_common(1)[0]
                names += ['first', 'last', 'top', 'freq']
                data += [lib.Timestamp(asint.min()),
                         lib.Timestamp(asint.max()),
                         lib.Timestamp(top), freq]
        else:

            lb = .5 * (1. - percentile_width / 100.)
            ub = 1. - lb

            def pretty_name(x):
                x *= 100
                if x == int(x):
                    return '%.0f%%' % x
                else:
                    return '%.1f%%' % x

            names = ['count']
            data = [self.count()]
            names += ['mean', 'std', 'min', pretty_name(lb), '50%',
                      pretty_name(ub), 'max']
            data += [self.mean(), self.std(), self.min(),
                     self.quantile(
                         lb), self.median(), self.quantile(ub),
                     self.max()]

        return self._constructor(data, index=names).__finalize__(self)

    def corr(self, other, method='pearson',
             min_periods=None):
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
            return pa.NA
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
        this, other = self.align(other, join='inner')
        if len(this) == 0:
            return pa.NA
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
        result = com.diff(_values_from_object(self), periods)
        return self._constructor(result, index=self.index).__finalize__(self)

    def autocorr(self):
        """
        Lag-1 autocorrelation

        Returns
        -------
        autocorr : float
        """
        return self.corr(self.shift(1))

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

#------------------------------------------------------------------------------
# Combination

    def append(self, to_append, verify_integrity=False):
        """
        Concatenate two or more Series. The indexes must not overlap

        Parameters
        ----------
        to_append : Series or list/tuple of Series
        verify_integrity : boolean, default False
            If True, raise Exception on creating index with duplicates

        Returns
        -------
        appended : Series
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
        level : int or name
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
            this, other = self.align(other, level=level, join='outer')
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
        return self._constructor(result, index=new_index).__finalize__(self)

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
            new_index = self.index + other.index
            new_name = _maybe_match_name(self, other)
            new_values = pa.empty(len(new_index), dtype=self.dtype)
            for i, idx in enumerate(new_index):
                lv = self.get(idx, fill_value)
                rv = other.get(idx, fill_value)
                new_values[i] = func(lv, rv)
        else:
            new_index = self.index
            new_values = func(self.values, other)
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
        new_index = self.index + other.index
        this = self.reindex(new_index, copy=False)
        other = other.reindex(new_index, copy=False)
        name = _maybe_match_name(self, other)
        rs_vals = com._where_compat(isnull(this), other.values, this.values)
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

        self._data = self._data.putmask(mask, other, inplace=True)
        self._maybe_update_cacher()

    #----------------------------------------------------------------------
    # Reindexing, sorting

    def sort(self, axis=0, kind='quicksort', order=None, ascending=True):
        """
        Sort values and index labels by value, in place. For compatibility with
        ndarray API. No return value

        Parameters
        ----------
        axis : int (can only be zero)
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm
        order : ignored
        ascending : boolean, default True
            Sort ascending. Passing False sorts descending

        See Also
        --------
        Series.order
        """

        # GH 5856/5863
        if self._is_cached:
            raise ValueError("This Series is a view of some other array, to "
                             "sort in-place you must create a copy")

        result = self.order(na_last=True, kind=kind,
                            ascending=ascending)

        self._update_inplace(result)

    def sort_index(self, ascending=True):
        """
        Sort object by labels (along an axis)

        Parameters
        ----------
        ascending : boolean or list, default True
            Sort ascending vs. descending. Specify list for multiple sort
            orders

        Examples
        --------
        >>> result1 = s.sort_index(ascending=False)
        >>> result2 = s.sort_index(ascending=[1, 0])

        Returns
        -------
        sorted_obj : Series
        """
        index = self.index
        if isinstance(index, MultiIndex):
            from pandas.core.groupby import _lexsort_indexer
            indexer = _lexsort_indexer(index.labels, orders=ascending)
            indexer = com._ensure_platform_int(indexer)
            new_labels = index.take(indexer)
        else:
            new_labels, indexer = index.order(return_indexer=True,
                                              ascending=ascending)

        new_values = self.values.take(indexer)
        return self._constructor(new_values,
                                 index=new_labels).__finalize__(self)

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

        """
        values = self.values
        mask = isnull(values)

        if mask.any():
            result = Series(
                -1, index=self.index, name=self.name, dtype='int64')
            notmask = -mask
            result[notmask] = np.argsort(values[notmask], kind=kind)
            return self._constructor(result,
                                     index=self.index).__finalize__(self)
        else:
            return self._constructor(
                np.argsort(values, kind=kind), index=self.index,
                dtype='int64').__finalize__(self)

    def rank(self, method='average', na_option='keep', ascending=True):
        """
        Compute data ranks (1 through n). Equal values are assigned a rank that
        is the average of the ranks of those values

        Parameters
        ----------
        method : {'average', 'min', 'max', 'first'}
            * average: average rank of group
            * min: lowest rank in group
            * max: highest rank in group
            * first: ranks assigned in order they appear in the array
        na_option : {'keep'}
            keep: leave NA values where they are
        ascending : boolean, default True
            False for ranks by high (1) to low (N)

        Returns
        -------
        ranks : Series
        """
        from pandas.core.algorithms import rank
        ranks = rank(self.values, method=method, na_option=na_option,
                     ascending=ascending)
        return self._constructor(ranks, index=self.index).__finalize__(self)

    def order(self, na_last=True, ascending=True, kind='mergesort'):
        """
        Sorts Series object, by value, maintaining index-value link

        Parameters
        ----------
        na_last : boolean (optional, default=True)
            Put NaN's at beginning or end
        ascending : boolean, default True
            Sort ascending. Passing False sorts descending
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'mergesort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm

        Returns
        -------
        y : Series
        """
        def _try_kind_sort(arr):
            # easier to ask forgiveness than permission
            try:
                # if kind==mergesort, it can fail for object dtype
                return arr.argsort(kind=kind)
            except TypeError:
                # stable sort not available for object dtype
                # uses the argsort default quicksort
                return arr.argsort(kind='quicksort')

        arr = self.values
        sortedIdx = pa.empty(len(self), dtype=np.int32)

        bad = isnull(arr)

        good = -bad
        idx = pa.arange(len(self))

        argsorted = _try_kind_sort(arr[good])

        if not ascending:
            argsorted = argsorted[::-1]

        if na_last:
            n = good.sum()
            sortedIdx[:n] = idx[good][argsorted]
            sortedIdx[n:] = idx[bad]
        else:
            n = bad.sum()
            sortedIdx[n:] = idx[good][argsorted]
            sortedIdx[:n] = idx[bad]

        return self._constructor(arr[sortedIdx], index=self.index[sortedIdx])\
                   .__finalize__(self)

    def sortlevel(self, level=0, ascending=True):
        """
        Sort Series with MultiIndex by chosen level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

        Parameters
        ----------
        level : int
        ascending : bool, default True

        Returns
        -------
        sorted : Series
        """
        if not isinstance(self.index, MultiIndex):
            raise TypeError('can only sort by level with a hierarchical index')

        new_index, indexer = self.index.sortlevel(level, ascending=ascending)
        new_values = self.values.take(indexer)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def swaplevel(self, i, j, copy=True):
        """
        Swap levels i and j in a MultiIndex

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : Series
        """
        new_index = self.index.swaplevel(i, j)
        return self._constructor(self.values, index=new_index,
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

    def unstack(self, level=-1):
        """
        Unstack, a.k.a. pivot, Series with MultiIndex to produce DataFrame

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to unstack, can pass level name

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
        return unstack(self, level)

    #----------------------------------------------------------------------
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
        values = self.values
        if com.is_datetime64_dtype(values.dtype):
            values = lib.map_infer(values, lib.Timestamp)

        if na_action == 'ignore':
            mask = isnull(values)

            def map_f(values, f):
                return lib.map_infer_mask(values, f, mask.view(pa.uint8))
        else:
            map_f = lib.map_infer

        if isinstance(arg, (dict, Series)):
            if isinstance(arg, dict):
                arg = self._constructor(arg)

            indexer = arg.index.get_indexer(values)
            new_values = com.take_1d(arg.values, indexer)
            return self._constructor(new_values,
                                     index=self.index).__finalize__(self)
        else:
            mapped = map_f(values, arg)
            return self._constructor(mapped,
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

        See also
        --------
        Series.map: For element-wise operations

        Returns
        -------
        y : Series or DataFrame if func returns a Series
        """
        if len(self) == 0:
            return Series()

        if kwds or args and not isinstance(func, np.ufunc):
            f = lambda x: func(x, *args, **kwds)
        else:
            f = func

        if isinstance(f, np.ufunc):
            return f(self)

        values = _values_from_object(self)
        if com.is_datetime64_dtype(values.dtype):
            values = lib.map_infer(values, lib.Timestamp)

        mapped = lib.map_infer(values, f, convert=convert_dtype)
        if len(mapped) and isinstance(mapped[0], Series):
            from pandas.core.frame import DataFrame
            return DataFrame(mapped.tolist(), index=self.index)
        else:
            return self._constructor(mapped,
                                     index=self.index).__finalize__(self)

    def _reduce(self, op, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        """ perform a reduction operation """
        return op(_values_from_object(self), skipna=skipna, **kwds)

    def _reindex_indexer(self, new_index, indexer, copy):
        if indexer is None:
            if copy:
                return self.copy()
            return self

        # be subclass-friendly
        new_values = com.take_1d(self.get_values(), indexer)
        return self._constructor(new_values, index=new_index)

    def _needs_reindex_multi(self, axes, method, level):
        """ check if we do need a multi reindex; this is for compat with
        higher dims
        """
        return False

    @Appender(generic._shared_docs['rename'] % _shared_doc_kwargs)
    def rename(self, index=None, **kwargs):
        return super(Series, self).rename(index=index, **kwargs)

    @Appender(generic._shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, index=None, **kwargs):
        return super(Series, self).reindex(index=index, **kwargs)

    def reindex_axis(self, labels, axis=0, **kwargs):
        """ for compatibility with higher dims """
        if axis != 0:
            raise ValueError("cannot reindex series on non-zero axis!")
        return self.reindex(index=labels, **kwargs)

    def take(self, indices, axis=0, convert=True):
        """
        Analogous to ndarray.take, return Series corresponding to requested
        indices

        Parameters
        ----------
        indices : list / array of ints
        convert : translate negative to positive indices (default)

        Returns
        -------
        taken : Series
        """
        # check/convert indicies here
        if convert:
            indices = _maybe_convert_indices(
                indices, len(self._get_axis(axis)))

        indices = com._ensure_platform_int(indices)
        new_index = self.index.take(indices)
        new_values = self.values.take(indices)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def isin(self, values):
        """
        Return a boolean :class:`~pandas.Series` showing whether each element
        in the :class:`~pandas.Series` is exactly contained in the passed
        sequence of ``values``.

        Parameters
        ----------
        values : list-like
            The sequence of values to test. Passing in a single string will
            raise a ``TypeError``. Instead, turn a single string into a
            ``list`` of one element.

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
        if not com.is_list_like(values):
            raise TypeError("only list-like objects are allowed to be passed"
                            " to Series.isin(), you passed a "
                            "{0!r}".format(type(values).__name__))

        # may need i8 conversion for proper membership testing
        comps = _values_from_object(self)
        if com.is_datetime64_dtype(self):
            from pandas.tseries.tools import to_datetime
            values = Series(to_datetime(values)).values.view('i8')
            comps = comps.view('i8')
        elif com.is_timedelta64_dtype(self):
            from pandas.tseries.timedeltas import to_timedelta
            values = Series(to_timedelta(values)).values.view('i8')
            comps = comps.view('i8')

        value_set = set(values)
        result = lib.ismember(comps, value_set)
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
        Read delimited file into Series

        Parameters
        ----------
        path : string file path or file handle / StringIO
        sep : string, default ','
            Field delimiter
        parse_dates : boolean, default True
            Parse dates. Different default from read_table
        header : int, default 0
            Row to use at header (skip prior rows)
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

        Returns
        -------
        y : Series
        """
        from pandas.core.frame import DataFrame
        df = DataFrame.from_csv(path, header=header, index_col=index_col,
                                sep=sep, parse_dates=parse_dates,
                                encoding=encoding,
                                infer_datetime_format=infer_datetime_format)
        result = df.icol(0)
        result.index.name = result.name = None
        return result

    def to_csv(self, path, index=True, sep=",", na_rep='',
               float_format=None, header=False,
               index_label=None, mode='w', nanRep=None, encoding=None,
               date_format=None):
        """
        Write Series to a comma-separated values (csv) file

        Parameters
        ----------
        path : string file path or file handle / StringIO
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
        """
        from pandas.core.frame import DataFrame
        df = DataFrame(self)
        df.to_csv(path, index=index, sep=sep, na_rep=na_rep,
                  float_format=float_format, header=header,
                  index_label=index_label, mode=mode, nanRep=nanRep,
                  encoding=encoding, date_format=date_format)

    def dropna(self, axis=0, inplace=False, **kwargs):
        """
        Return Series without null values

        Returns
        -------
        valid : Series
        inplace : boolean, default False
            Do operation in place.
        """
        axis = self._get_axis_number(axis or 0)
        result = remove_na(self)
        if inplace:
            self._update_inplace(result)
        else:
            return result

    valid = lambda self, inplace=False, **kwargs: self.dropna(inplace=inplace,
                                                              **kwargs)

    def first_valid_index(self):
        """
        Return label for first non-NA/null value
        """
        if len(self) == 0:
            return None

        mask = isnull(self.values)
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

        mask = isnull(self.values[::-1])
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[len(self) - i - 1]

    #----------------------------------------------------------------------
    # Time series-oriented methods

    def asof(self, where):
        """
        Return last good (non-NaN) value in TimeSeries if value is NaN for
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

        values = self.values

        if not hasattr(where, '__iter__'):
            start = self.index[0]
            if isinstance(self.index, PeriodIndex):
                where = Period(where, freq=self.index.freq).ordinal
                start = start.ordinal

            if where < start:
                return pa.NA
            loc = self.index.searchsorted(where, side='right')
            if loc > 0:
                loc -= 1
            while isnull(values[loc]) and loc > 0:
                loc -= 1
            return values[loc]

        if not isinstance(where, Index):
            where = Index(where)

        locs = self.index.asof_locs(where, notnull(values))
        new_values = com.take_1d(values, locs)
        return self._constructor(new_values, index=where).__finalize__(self)

    @property
    def weekday(self):
        return self._constructor([d.weekday() for d in self.index],
                                 index=self.index).__finalize__(self)

    def tz_convert(self, tz, copy=True):
        """
        Convert TimeSeries to target time zone

        Parameters
        ----------
        tz : string or pytz.timezone object
        copy : boolean, default True
            Also make a copy of the underlying data

        Returns
        -------
        converted : TimeSeries
        """
        new_index = self.index.tz_convert(tz)

        new_values = self.values
        if copy:
            new_values = new_values.copy()

        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def tz_localize(self, tz, copy=True, infer_dst=False):
        """
        Localize tz-naive TimeSeries to target time zone
        Entries will retain their "naive" value but will be annotated as
        being relative to the specified tz.

        After localizing the TimeSeries, you may use tz_convert() to
        get the Datetime values recomputed to a different tz.

        Parameters
        ----------
        tz : string or pytz.timezone object
        copy : boolean, default True
            Also make a copy of the underlying data
        infer_dst : boolean, default False
            Attempt to infer fall dst-transition hours based on order

        Returns
        -------
        localized : TimeSeries
        """
        from pandas.tseries.index import DatetimeIndex

        if not isinstance(self.index, DatetimeIndex):
            if len(self.index) > 0:
                raise Exception('Cannot tz-localize non-time series')

            new_index = DatetimeIndex([], tz=tz)
        else:
            new_index = self.index.tz_localize(tz, infer_dst=infer_dst)

        new_values = self.values
        if copy:
            new_values = new_values.copy()

        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    @cache_readonly
    def str(self):
        from pandas.core.strings import StringMethods
        return StringMethods(self)

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
        ts : TimeSeries with DatetimeIndex
        """
        new_values = self.values
        if copy:
            new_values = new_values.copy()

        new_index = self.index.to_timestamp(freq=freq, how=how)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def to_period(self, freq=None, copy=True):
        """
        Convert TimeSeries from DatetimeIndex to PeriodIndex with desired
        frequency (inferred from index if not passed)

        Parameters
        ----------
        freq : string, default

        Returns
        -------
        ts : TimeSeries with PeriodIndex
        """
        new_values = self.values
        if copy:
            new_values = new_values.copy()

        if freq is None:
            freq = self.index.freqstr or self.index.inferred_freq
        new_index = self.index.to_period(freq=freq)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

Series._setup_axes(['index'], info_axis=0, stat_axis=0,
                   aliases={'rows': 0})
Series._add_numeric_operations()
_INDEX_TYPES = ndarray, Index, list, tuple

#------------------------------------------------------------------------------
# Supplementary functions


def remove_na(series):
    """
    Return series containing only true/non-NaN values, possibly empty.
    """
    return series[notnull(_values_from_object(series))]


def _sanitize_array(data, index, dtype=None, copy=False,
                    raise_cast_failure=False):
    if dtype is not None:
        dtype = np.dtype(dtype)

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
            arr = _possibly_cast_to_datetime(arr, dtype)
            subarr = pa.array(arr, dtype=dtype, copy=copy)
        except (ValueError, TypeError):
            if dtype is not None and raise_cast_failure:
                raise
            else:  # pragma: no cover
                subarr = pa.array(arr, dtype=object, copy=copy)
        return subarr

    # GH #846
    if isinstance(data, (pa.Array, Series)):
        subarr = np.array(data, copy=False)
        if dtype is not None:

            # possibility of nan -> garbage
            if com.is_float_dtype(data.dtype) and com.is_integer_dtype(dtype):
                if not isnull(data).any():
                    subarr = _try_cast(data, True)
                elif copy:
                    subarr = data.copy()
            else:
                if (com.is_datetime64_dtype(data.dtype) and
                        not com.is_datetime64_dtype(dtype)):
                    if dtype == object:
                        ints = np.asarray(data).view('i8')
                        subarr = tslib.ints_to_pydatetime(ints)
                    elif raise_cast_failure:
                        raise TypeError('Cannot cast datetime64 to %s' % dtype)
                else:
                    subarr = _try_cast(data, True)
        else:
            # don't coerce Index types
            # e.g. indexes can have different conversions (so don't fast path them)
            # GH 6140
            subarr = _try_cast(data, not isinstance(data, Index))

        if copy:
            subarr = data.copy()

    elif isinstance(data, list) and len(data) > 0:
        if dtype is not None:
            try:
                subarr = _try_cast(data, False)
            except Exception:
                if raise_cast_failure:  # pragma: no cover
                    raise
                subarr = pa.array(data, dtype=object, copy=copy)
                subarr = lib.maybe_convert_objects(subarr)

        else:
            subarr = _possibly_convert_platform(data)

        subarr = _possibly_cast_to_datetime(subarr, dtype)

    else:
        subarr = _try_cast(data, False)

    # scalar like
    if subarr.ndim == 0:
        if isinstance(data, list):  # pragma: no cover
            subarr = pa.array(data, dtype=object)
        elif index is not None:
            value = data

            # figure out the dtype from the value (upcast if necessary)
            if dtype is None:
                dtype, value = _infer_dtype_from_scalar(value)
            else:
                # need to possibly convert the value here
                value = _possibly_cast_to_datetime(value, dtype)

            subarr = pa.empty(len(index), dtype=dtype)
            subarr.fill(value)

        else:
            return subarr.item()

    # the result that we want
    elif subarr.ndim == 1:
        if index is not None:

            # a 1-element ndarray
            if len(subarr) != len(index) and len(subarr) == 1:
                value = subarr[0]
                subarr = pa.empty(len(index), dtype=subarr.dtype)
                subarr.fill(value)

    elif subarr.ndim > 1:
        if isinstance(data, pa.Array):
            raise Exception('Data must be 1-dimensional')
        else:
            subarr = _asarray_tuplesafe(data, dtype=dtype)

    # This is to prevent mixed-type Series getting all casted to
    # NumPy string type, e.g. NaN --> '-1#IND'.
    if issubclass(subarr.dtype.type, compat.string_types):
        subarr = pa.array(data, dtype=object, copy=copy)

    return subarr

# backwards compatiblity
TimeSeries = Series

#----------------------------------------------------------------------
# Add plotting methods to Series

import pandas.tools.plotting as _gfx

Series.plot = _gfx.plot_series
Series.hist = _gfx.hist_series

# Add arithmetic!
ops.add_flex_arithmetic_methods(Series, **ops.series_flex_funcs)
ops.add_special_arithmetic_methods(Series, **ops.series_special_funcs)
