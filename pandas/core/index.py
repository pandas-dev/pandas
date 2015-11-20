# pylint: disable=E1101,E1103,W0232
import datetime
import warnings
import operator
from functools import partial
from sys import getsizeof

import numpy as np
import pandas.tslib as tslib
import pandas.lib as lib
import pandas.algos as _algos
import pandas.index as _index
from pandas.lib import Timestamp, Timedelta, is_datetime_array

from pandas.compat import range, zip, lrange, lzip, u, map
from pandas import compat
from pandas.core import algorithms
from pandas.core.base import PandasObject, FrozenList, FrozenNDArray, IndexOpsMixin, PandasDelegate
import pandas.core.base as base
from pandas.util.decorators import (Appender, Substitution, cache_readonly,
                                    deprecate, deprecate_kwarg)
import pandas.core.common as com
from pandas.core.missing import _clean_reindex_fill_method
from pandas.core.common import (isnull, array_equivalent, is_dtype_equal, is_object_dtype,
                                is_datetimetz, ABCSeries, ABCCategorical, ABCPeriodIndex,
                                _values_from_object, is_float, is_integer, is_iterator, is_categorical_dtype,
                                _ensure_object, _ensure_int64, is_bool_indexer,
                                is_list_like, is_bool_dtype, is_null_slice, is_integer_dtype)
from pandas.core.strings import StringAccessorMixin
from pandas.core.config import get_option
from pandas.io.common import PerformanceWarning


# simplify
default_pprint = lambda x, max_seq_items=None: com.pprint_thing(x,
                                                                escape_chars=('\t', '\r', '\n'),
                                                                quote_strings=True,
                                                                max_seq_items=max_seq_items)


__all__ = ['Index']


_unsortable_types = frozenset(('mixed', 'mixed-integer'))

_index_doc_kwargs = dict(klass='Index', inplace='',
                         duplicated='np.array')
_index_shared_docs = dict()


def _try_get_item(x):
    try:
        return x.item()
    except AttributeError:
        return x

class InvalidIndexError(Exception):
    pass

_o_dtype = np.dtype(object)
_Identity = object

def _new_Index(cls, d):
    """ This is called upon unpickling, rather than the default which doesn't have arguments
        and breaks __new__ """
    return cls.__new__(cls, **d)

class Index(IndexOpsMixin, StringAccessorMixin, PandasObject):

    """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: object)
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index
    tupleize_cols : bool (default: True)
        When True, attempt to create a MultiIndex if possible

    Notes
    -----
    An Index instance can **only** contain hashable objects
    """
    # To hand over control to subclasses
    _join_precedence = 1

    # Cython methods
    _groupby = _algos.groupby_object
    _arrmap = _algos.arrmap_object
    _left_indexer_unique = _algos.left_join_indexer_unique_object
    _left_indexer = _algos.left_join_indexer_object
    _inner_indexer = _algos.inner_join_indexer_object
    _outer_indexer = _algos.outer_join_indexer_object
    _box_scalars = False

    _typ = 'index'
    _data = None
    _id = None
    name = None
    asi8 = None
    _comparables = ['name']
    _attributes = ['name']
    _allow_index_ops = True
    _allow_datetime_index_ops = False
    _allow_period_index_ops = False
    _is_numeric_dtype = False
    _can_hold_na = True

    _engine_type = _index.ObjectEngine

    def __new__(cls, data=None, dtype=None, copy=False, name=None, fastpath=False,
                tupleize_cols=True, **kwargs):

        # no class inference!
        if fastpath:
            return cls._simple_new(data, name)

        if is_categorical_dtype(data) or is_categorical_dtype(dtype):
            return CategoricalIndex(data, copy=copy, name=name, **kwargs)

        if isinstance(data, (np.ndarray, Index, ABCSeries)):
            if issubclass(data.dtype.type, np.datetime64) or is_datetimetz(data):
                from pandas.tseries.index import DatetimeIndex
                result = DatetimeIndex(data, copy=copy, name=name, **kwargs)
                if dtype is not None and _o_dtype == dtype:
                    return Index(result.to_pydatetime(), dtype=_o_dtype)
                else:
                    return result
            elif issubclass(data.dtype.type, np.timedelta64):
                from pandas.tseries.tdi import TimedeltaIndex
                result = TimedeltaIndex(data, copy=copy, name=name, **kwargs)
                if dtype is not None and _o_dtype == dtype:
                    return Index(result.to_pytimedelta(), dtype=_o_dtype)
                else:
                    return result

            if dtype is not None:
                try:
                    data = np.array(data, dtype=dtype, copy=copy)
                except (TypeError, ValueError):
                    pass

            # maybe coerce to a sub-class
            from pandas.tseries.period import PeriodIndex
            if isinstance(data, PeriodIndex):
                return PeriodIndex(data, copy=copy, name=name, **kwargs)
            if issubclass(data.dtype.type, np.integer):
                return Int64Index(data, copy=copy, dtype=dtype, name=name)
            elif issubclass(data.dtype.type, np.floating):
                return Float64Index(data, copy=copy, dtype=dtype, name=name)
            elif issubclass(data.dtype.type, np.bool) or is_bool_dtype(data):
                subarr = data.astype('object')
            else:
                subarr = com._asarray_tuplesafe(data, dtype=object)

            # _asarray_tuplesafe does not always copy underlying data,
            # so need to make sure that this happens
            if copy:
                subarr = subarr.copy()

            if dtype is None:
                inferred = lib.infer_dtype(subarr)
                if inferred == 'integer':
                    return Int64Index(subarr.astype('i8'), copy=copy, name=name)
                elif inferred in ['floating', 'mixed-integer-float']:
                    return Float64Index(subarr, copy=copy, name=name)
                elif inferred == 'boolean':
                    # don't support boolean explicity ATM
                    pass
                elif inferred != 'string':
                    if (inferred.startswith('datetime') or
                        tslib.is_timestamp_array(subarr)):
                        from pandas.tseries.index import DatetimeIndex
                        return DatetimeIndex(subarr, copy=copy, name=name, **kwargs)
                    elif (inferred.startswith('timedelta') or
                          lib.is_timedelta_array(subarr)):
                        from pandas.tseries.tdi import TimedeltaIndex
                        return TimedeltaIndex(subarr, copy=copy, name=name, **kwargs)
                    elif inferred == 'period':
                        return PeriodIndex(subarr, name=name, **kwargs)
            return cls._simple_new(subarr, name)

        elif hasattr(data, '__array__'):
            return Index(np.asarray(data), dtype=dtype, copy=copy, name=name,
                         **kwargs)
        elif data is None or np.isscalar(data):
            cls._scalar_data_error(data)
        else:
            if tupleize_cols and isinstance(data, list) and data and isinstance(data[0], tuple):

                # we must be all tuples, otherwise don't construct
                # 10697
                if all( isinstance(e, tuple) for e in data ):
                    try:
                        # must be orderable in py3
                        if compat.PY3:
                            sorted(data)
                        return MultiIndex.from_tuples(
                            data, names=name or kwargs.get('names'))
                    except (TypeError, KeyError):
                        # python2 - MultiIndex fails on mixed types
                        pass
            # other iterable of some kind
            subarr = com._asarray_tuplesafe(data, dtype=object)
            return Index(subarr, dtype=dtype, copy=copy, name=name, **kwargs)

    @classmethod
    def _simple_new(cls, values, name=None, dtype=None, **kwargs):
        """
        we require the we have a dtype compat for the values
        if we are passed a non-dtype compat, then coerce using the constructor

        Must be careful not to recurse.
        """
        if not hasattr(values, 'dtype'):
            if values is None and dtype is not None:
                values = np.empty(0, dtype=dtype)
            else:
                values = np.array(values,copy=False)
                if is_object_dtype(values):
                    values = cls(values, name=name, dtype=dtype, **kwargs)._values

        result = object.__new__(cls)
        result._data = values
        result.name = name
        for k, v in compat.iteritems(kwargs):
            setattr(result,k,v)
        result._reset_identity()
        return result

    def _update_inplace(self, result, **kwargs):
        # guard when called from IndexOpsMixin
        raise TypeError("Index can't be updated inplace")

    def is_(self, other):
        """
        More flexible, faster check like ``is`` but that works through views

        Note: this is *not* the same as ``Index.identical()``, which checks
        that metadata is also the same.

        Parameters
        ----------
        other : object
            other object to compare against.

        Returns
        -------
        True if both have same underlying data, False otherwise : bool
        """
        # use something other than None to be clearer
        return self._id is getattr(other, '_id', Ellipsis)

    def _reset_identity(self):
        """Initializes or resets ``_id`` attribute with new object"""
        self._id = _Identity()

    # ndarray compat
    def __len__(self):
        """
        return the length of the Index
        """
        return len(self._data)

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return self._data.view(np.ndarray)

    def __array_wrap__(self, result, context=None):
        """
        Gets called after a ufunc
        """
        if is_bool_dtype(result):
            return result

        attrs = self._get_attributes_dict()
        attrs = self._maybe_update_attributes(attrs)
        return Index(result, **attrs)

    @cache_readonly
    def dtype(self):
        """ return the dtype object of the underlying data """
        return self._data.dtype

    @cache_readonly
    def dtype_str(self):
        """ return the dtype str of the underlying data """
        return str(self.dtype)

    @property
    def values(self):
        """ return the underlying data as an ndarray """
        return self._data.view(np.ndarray)

    def get_values(self):
        """ return the underlying data as an ndarray """
        return self.values

    # ops compat
    def tolist(self):
        """
        return a list of the Index values
        """
        return list(self.values)

    def repeat(self, n):
        """
        return a new Index of the values repeated n times

        See also
        --------
        numpy.ndarray.repeat
        """
        return self._shallow_copy(self._values.repeat(n))

    def ravel(self, order='C'):
        """
        return an ndarray of the flattened values of the underlying data

        See also
        --------
        numpy.ndarray.ravel
        """
        return self._values.ravel(order=order)

    # construction helpers
    @classmethod
    def _scalar_data_error(cls, data):
        raise TypeError(
            '{0}(...) must be called with a collection of some kind, {1} was '
            'passed'.format(cls.__name__, repr(data))
        )

    @classmethod
    def _string_data_error(cls, data):
        raise TypeError('String dtype not supported, you may need '
                        'to explicitly cast to a numeric type')

    @classmethod
    def _coerce_to_ndarray(cls, data):
        """coerces data to ndarray, raises on scalar data. Converts other
        iterables to list first and then to array. Does not touch ndarrays."""

        if not isinstance(data, (np.ndarray, Index)):
            if data is None or np.isscalar(data):
                cls._scalar_data_error(data)

            # other iterable of some kind
            if not isinstance(data, (ABCSeries, list, tuple)):
                data = list(data)
            data = np.asarray(data)
        return data

    def _get_attributes_dict(self):
        """ return an attributes dict for my class """
        return dict([ (k,getattr(self,k,None)) for k in self._attributes])

    def view(self, cls=None):

        # we need to see if we are subclassing an
        # index type here
        if cls is not None and not hasattr(cls,'_typ'):
            result = self._data.view(cls)
        else:
            result = self._shallow_copy()
        if isinstance(result, Index):
            result._id = self._id
        return result

    def _shallow_copy(self, values=None, infer=False, **kwargs):
        """
        create a new Index, don't copy the data, use the same object attributes
        with passed in attributes taking precedence

        *this is an internal non-public method*

        Parameters
        ----------
        values : the values to create the new Index, optional
        infer : boolean, default False
            if True, infer the new type of the passed values
        kwargs : updates the default attributes for this Index
        """
        if values is None:
            values = self.values
        attributes = self._get_attributes_dict()
        attributes.update(kwargs)

        if infer:
            attributes['copy'] = False
            return Index(values, **attributes)

        return self.__class__._simple_new(values,**attributes)

    def _coerce_scalar_to_index(self, item):
        """
        we need to coerce a scalar to a compat for our index type

        Parameters
        ----------
        item : scalar item to coerce
        """
        return Index([item], dtype=self.dtype, **self._get_attributes_dict())

    def copy(self, names=None, name=None, dtype=None, deep=False):
        """
        Make a copy of this object.  Name and dtype sets those attributes on
        the new object.

        Parameters
        ----------
        name : string, optional
        dtype : numpy dtype or pandas type

        Returns
        -------
        copy : Index

        Notes
        -----
        In most cases, there should be no functional difference from using
        ``deep``, but if ``deep`` is passed it will attempt to deepcopy.
        """
        if names is not None and name is not None:
            raise TypeError("Can only provide one of `names` and `name`")
        if deep:
            from copy import deepcopy
            new_index = self._shallow_copy(self._data.copy())
            name = name or deepcopy(self.name)
        else:
            new_index = self._shallow_copy()
            name = self.name
        if name is not None:
            names = [name]
        if names:
            new_index = new_index.set_names(names)
        if dtype:
            new_index = new_index.astype(dtype)
        return new_index

    __copy__ = copy

    def __unicode__(self):
        """
        Return a string representation for this object.

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        klass = self.__class__.__name__
        data = self._format_data()
        attrs = self._format_attrs()
        space = self._format_space()

        prepr = (u(",%s") % space).join([u("%s=%s") % (k, v)
                                          for k, v in attrs])

        # no data provided, just attributes
        if data is None:
            data = ''

        res = u("%s(%s%s)") % (klass,
                               data,
                               prepr)

        return res

    def _format_space(self):

        # using space here controls if the attributes
        # are line separated or not (the default)

        #max_seq_items = get_option('display.max_seq_items')
        #if len(self) > max_seq_items:
        #    space = "\n%s" % (' ' * (len(klass) + 1))
        return " "

    @property
    def _formatter_func(self):
        """
        Return the formatted data as a unicode string
        """
        return default_pprint

    def _format_data(self):
        """
        Return the formatted data as a unicode string
        """
        from pandas.core.format import get_console_size, _get_adjustment
        display_width, _ = get_console_size()
        if display_width is None:
            display_width = get_option('display.width') or 80

        space1 = "\n%s" % (' ' * (len(self.__class__.__name__) + 1))
        space2 = "\n%s" % (' ' * (len(self.__class__.__name__) + 2))

        n = len(self)
        sep = ','
        max_seq_items = get_option('display.max_seq_items') or n
        formatter = self._formatter_func

        # do we want to justify (only do so for non-objects)
        is_justify = not (self.inferred_type in ('string', 'unicode') or
                          (self.inferred_type == 'categorical' and
                           is_object_dtype(self.categories)))

        # are we a truncated display
        is_truncated = n > max_seq_items

        # adj can optionaly handle unicode eastern asian width
        adj = _get_adjustment()

        def _extend_line(s, line, value, display_width, next_line_prefix):

            if adj.len(line.rstrip()) + adj.len(value.rstrip()) >= display_width:
                s += line.rstrip()
                line = next_line_prefix
            line += value
            return s, line

        def best_len(values):
            if values:
                return max([adj.len(x) for x in values])
            else:
                return 0

        if n == 0:
            summary = '[], '
        elif n == 1:
            first = formatter(self[0])
            summary = '[%s], ' % first
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = '[%s, %s], ' % (first, last)
        else:

            if n > max_seq_items:
                n = min(max_seq_items//2,10)
                head = [ formatter(x) for x in self[:n] ]
                tail = [ formatter(x) for x in self[-n:] ]
            else:
                head = []
                tail = [ formatter(x) for x in self ]

            # adjust all values to max length if needed
            if is_justify:

                # however, if we are not truncated and we are only a single line, then don't justify
                if is_truncated or not (len(', '.join(head)) < display_width and len(', '.join(tail)) < display_width):
                    max_len = max(best_len(head), best_len(tail))
                    head = [x.rjust(max_len) for x in head]
                    tail = [x.rjust(max_len) for x in tail]

            summary = ""
            line = space2

            for i in range(len(head)):
                word = head[i] + sep + ' '
                summary, line = _extend_line(summary, line, word,
                                             display_width, space2)

            if is_truncated:
                # remove trailing space of last line
                summary += line.rstrip() + space2 + '...'
                line = space2

            for i in range(len(tail)-1):
                word = tail[i] + sep + ' '
                summary, line = _extend_line(summary, line, word,
                                             display_width, space2)

            # last value: no sep added + 1 space of width used for trailing ','
            summary, line = _extend_line(summary, line, tail[-1],
                                         display_width - 2, space2)
            summary += line
            summary += '],'

            if len(summary) > (display_width):
                summary += space1
            else:  # one row
                summary += ' '

            # remove initial space
            summary = '[' + summary[len(space2):]

        return summary

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        attrs = []
        attrs.append(('dtype',"'%s'" % self.dtype))
        if self.name is not None:
            attrs.append(('name',default_pprint(self.name)))
        max_seq_items = get_option('display.max_seq_items') or len(self)
        if len(self) > max_seq_items:
            attrs.append(('length',len(self)))
        return attrs

    def to_series(self, **kwargs):
        """
        Create a Series with both index and values equal to the index keys
        useful with map for returning an indexer based on an index

        Returns
        -------
        Series : dtype will be based on the type of the Index values.
        """

        from pandas import Series
        return Series(self._to_embed(), index=self, name=self.name)

    def _to_embed(self, keep_tz=False):
        """
        *this is an internal non-public method*

        return an array repr of this object, potentially casting to object

        """
        return self.values.copy()

    def astype(self, dtype):
        return Index(self.values.astype(dtype), name=self.name,
                     dtype=dtype)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self

    def to_datetime(self, dayfirst=False):
        """
        For an Index containing strings or datetime.datetime objects, attempt
        conversion to DatetimeIndex
        """
        from pandas.tseries.index import DatetimeIndex
        if self.inferred_type == 'string':
            from dateutil.parser import parse
            parser = lambda x: parse(x, dayfirst=dayfirst)
            parsed = lib.try_parse_dates(self.values, parser=parser)
            return DatetimeIndex(parsed)
        else:
            return DatetimeIndex(self.values)

    def _assert_can_do_setop(self, other):
        if not com.is_list_like(other):
            raise TypeError('Input must be Index or array-like')
        return True

    def _convert_can_do_setop(self, other):
        if not isinstance(other, Index):
            other = Index(other, name=self.name)
            result_name = self.name
        else:
            result_name = self.name if self.name == other.name else None
        return other, result_name

    @property
    def nlevels(self):
        return 1

    def _get_names(self):
        return FrozenList((self.name,))

    def _set_names(self, values, level=None):
        if len(values) != 1:
            raise ValueError('Length of new names must be 1, got %d'
                             % len(values))
        self.name = values[0]

    names = property(fset=_set_names, fget=_get_names)

    def set_names(self, names, level=None, inplace=False):
        """
        Set new names on index. Defaults to returning new index.

        Parameters
        ----------
        names : str or sequence
            name(s) to set
        level : int or level name, or sequence of int / level names (default None)
            If the index is a MultiIndex (hierarchical), level(s) to set (None for all levels)
            Otherwise level must be None
        inplace : bool
            if True, mutates in place

        Returns
        -------
        new index (of same type and class...etc) [if inplace, returns None]

        Examples
        --------
        >>> Index([1, 2, 3, 4]).set_names('foo')
        Int64Index([1, 2, 3, 4], dtype='int64')
        >>> Index([1, 2, 3, 4]).set_names(['foo'])
        Int64Index([1, 2, 3, 4], dtype='int64')
        >>> idx = MultiIndex.from_tuples([(1, u'one'), (1, u'two'),
                                          (2, u'one'), (2, u'two')],
                                          names=['foo', 'bar'])
        >>> idx.set_names(['baz', 'quz'])
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'baz', u'quz'])
        >>> idx.set_names('baz', level=0)
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'baz', u'bar'])
        """
        if level is not None and self.nlevels == 1:
            raise ValueError('Level must be None for non-MultiIndex')

        if level is not None and not is_list_like(level) and is_list_like(names):
            raise TypeError("Names must be a string")

        if not is_list_like(names) and level is None and self.nlevels > 1:
            raise TypeError("Must pass list-like as `names`.")

        if not is_list_like(names):
            names = [names]
        if level is not None and not is_list_like(level):
            level = [level]

        if inplace:
            idx = self
        else:
            idx = self._shallow_copy()
        idx._set_names(names, level=level)
        if not inplace:
            return idx

    def rename(self, name, inplace=False):
        """
        Set new names on index. Defaults to returning new index.

        Parameters
        ----------
        name : str or list
            name to set
        inplace : bool
            if True, mutates in place

        Returns
        -------
        new index (of same type and class...etc) [if inplace, returns None]
        """
        return self.set_names([name], inplace=inplace)

    @property
    def _has_complex_internals(self):
        # to disable groupby tricks in MultiIndex
        return False

    def summary(self, name=None):
        if len(self) > 0:
            head = self[0]
            if hasattr(head, 'format') and\
               not isinstance(head, compat.string_types):
                head = head.format()
            tail = self[-1]
            if hasattr(tail, 'format') and\
               not isinstance(tail, compat.string_types):
                tail = tail.format()
            index_summary = ', %s to %s' % (com.pprint_thing(head),
                                            com.pprint_thing(tail))
        else:
            index_summary = ''

        if name is None:
            name = type(self).__name__
        return '%s: %s entries%s' % (name, len(self), index_summary)

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return self.values

    _na_value = np.nan
    """The expected NA value to use with this index."""

    @property
    def is_monotonic(self):
        """ alias for is_monotonic_increasing (deprecated) """
        return self._engine.is_monotonic_increasing

    @property
    def is_monotonic_increasing(self):
        """
        return if the index is monotonic increasing (only equal or
        increasing) values.
        """
        return self._engine.is_monotonic_increasing

    @property
    def is_monotonic_decreasing(self):
        """
        return if the index is monotonic decreasing (only equal or
        decreasing) values.
        """
        return self._engine.is_monotonic_decreasing

    def is_lexsorted_for_tuple(self, tup):
        return True

    @cache_readonly(allow_setting=True)
    def is_unique(self):
        """ return if the index has unique values """
        return self._engine.is_unique

    @property
    def has_duplicates(self):
        return not self.is_unique

    def is_boolean(self):
        return self.inferred_type in ['boolean']

    def is_integer(self):
        return self.inferred_type in ['integer']

    def is_floating(self):
        return self.inferred_type in ['floating', 'mixed-integer-float']

    def is_numeric(self):
        return self.inferred_type in ['integer', 'floating']

    def is_object(self):
        return is_object_dtype(self.dtype)

    def is_categorical(self):
        return self.inferred_type in ['categorical']

    def is_mixed(self):
        return 'mixed' in self.inferred_type

    def holds_integer(self):
        return self.inferred_type in ['integer', 'mixed-integer']

    def _convert_scalar_indexer(self, key, kind=None):
        """
        convert a scalar indexer

        Parameters
        ----------
        key : label of the slice bound
        kind : optional, type of the indexing operation (loc/ix/iloc/None)

        right now we are converting
        floats -> ints if the index supports it
        """

        def to_int():
            ikey = int(key)
            if ikey != key:
                return self._invalid_indexer('label', key)
            return ikey

        if kind == 'iloc':
            if is_integer(key):
                return key
            elif is_float(key):
                key = to_int()
                warnings.warn("scalar indexers for index type {0} should be integers and not floating point".format(
                    type(self).__name__), FutureWarning, stacklevel=5)
                return key
            return self._invalid_indexer('label', key)

        if is_float(key):
            if isnull(key):
                return self._invalid_indexer('label', key)
            warnings.warn("scalar indexers for index type {0} should be integers and not floating point".format(
                type(self).__name__), FutureWarning, stacklevel=3)
            return to_int()

        return key

    def _convert_slice_indexer_getitem(self, key, is_index_slice=False):
        """ called from the getitem slicers, determine how to treat the key
            whether positional or not """
        if self.is_integer() or is_index_slice:
            return key
        return self._convert_slice_indexer(key)

    def _convert_slice_indexer(self, key, kind=None):
        """
        convert a slice indexer. disallow floats in the start/stop/step

        Parameters
        ----------
        key : label of the slice bound
        kind : optional, type of the indexing operation (loc/ix/iloc/None)
        """

        # if we are not a slice, then we are done
        if not isinstance(key, slice):
            return key

        # validate iloc
        if kind == 'iloc':

            # need to coerce to_int if needed
            def f(c):
                v = getattr(key,c)
                if v is None or is_integer(v):
                    return v

                # warn if it's a convertible float
                if v == int(v):
                    warnings.warn("slice indexers when using iloc should be integers "
                                  "and not floating point", FutureWarning, stacklevel=7)
                    return int(v)

                self._invalid_indexer('slice {0} value'.format(c), v)

            return slice(*[ f(c) for c in ['start','stop','step']])

        # validate slicers
        def validate(v):
            if v is None or is_integer(v):
                return True

            # dissallow floats (except for .ix)
            elif is_float(v):
                if kind == 'ix':
                    return True

                return False

            return True
        for c in ['start','stop','step']:
            v = getattr(key,c)
            if not validate(v):
                self._invalid_indexer('slice {0} value'.format(c), v)

        # figure out if this is a positional indexer
        start, stop, step = key.start, key.stop, key.step

        def is_int(v):
            return v is None or is_integer(v)

        is_null_slicer = start is None and stop is None
        is_index_slice = is_int(start) and is_int(stop)
        is_positional = is_index_slice and not self.is_integer()

        if kind == 'getitem':
            return self._convert_slice_indexer_getitem(
                key, is_index_slice=is_index_slice)

        # convert the slice to an indexer here

        # if we are mixed and have integers
        try:
            if is_positional and self.is_mixed():
                if start is not None:
                    i = self.get_loc(start)
                if stop is not None:
                    j = self.get_loc(stop)
                is_positional = False
        except KeyError:
            if self.inferred_type == 'mixed-integer-float':
                raise

        if is_null_slicer:
            indexer = key
        elif is_positional:
            indexer = key
        else:
            try:
                indexer = self.slice_indexer(start, stop, step)
            except Exception:
                if is_index_slice:
                    if self.is_integer():
                        raise
                    else:
                        indexer = key
                else:
                    raise

        return indexer

    def _convert_list_indexer(self, keyarr, kind=None):
        """
        passed a key that is tuplesafe that is integer based
        and we have a mixed index (e.g. number/labels). figure out
        the indexer. return None if we can't help
        """
        if kind in [None, 'iloc', 'ix'] and is_integer_dtype(keyarr) \
           and not self.is_floating() and not isinstance(keyarr, ABCPeriodIndex):

            if self.inferred_type == 'mixed-integer':
                indexer = self.get_indexer(keyarr)
                if (indexer >= 0).all():
                    return indexer
                # missing values are flagged as -1 by get_indexer and negative indices are already
                # converted to positive indices in the above if-statement, so the negative flags are changed to
                # values outside the range of indices so as to trigger an IndexError in maybe_convert_indices
                indexer[indexer < 0] = len(self)
                from pandas.core.indexing import maybe_convert_indices
                return maybe_convert_indices(indexer, len(self))

            elif not self.inferred_type == 'integer':
                keyarr = np.where(keyarr < 0,
                                  len(self) + keyarr, keyarr)
                return keyarr

        return None

    def _invalid_indexer(self, form, key):
        """ consistent invalid indexer message """
        raise TypeError("cannot do {form} indexing on {klass} with these "
                        "indexers [{key}] of {kind}".format(form=form,
                                                           klass=type(self),
                                                           key=key,
                                                           kind=type(key)))

    def get_duplicates(self):
        from collections import defaultdict
        counter = defaultdict(lambda: 0)
        for k in self.values:
            counter[k] += 1
        return sorted(k for k, v in compat.iteritems(counter) if v > 1)

    _get_duplicates = get_duplicates

    def _cleanup(self):
        self._engine.clear_mapping()

    @cache_readonly
    def _engine(self):
        # property, for now, slow to look up
        return self._engine_type(lambda: self.values, len(self))

    def _validate_index_level(self, level):
        """
        Validate index level.

        For single-level Index getting level number is a no-op, but some
        verification must be done like in MultiIndex.

        """
        if isinstance(level, int):
            if level < 0 and level != -1:
                raise IndexError("Too many levels: Index has only 1 level,"
                                 " %d is not a valid level number" % (level,))
            elif level > 0:
                raise IndexError("Too many levels:"
                                 " Index has only 1 level, not %d" %
                                 (level + 1))
        elif level != self.name:
            raise KeyError('Level %s must be same as name (%s)'
                           % (level, self.name))

    def _get_level_number(self, level):
        self._validate_index_level(level)
        return 0

    @cache_readonly
    def inferred_type(self):
        """ return a string of the type inferred from the values """
        return lib.infer_dtype(self)

    def is_type_compatible(self, kind):
        return kind == self.inferred_type

    @cache_readonly
    def is_all_dates(self):
        if self._data is None:
            return False
        return is_datetime_array(_ensure_object(self.values))

    def __iter__(self):
        return iter(self.values)

    def __reduce__(self):
        d = dict(data=self._data)
        d.update(self._get_attributes_dict())
        return _new_Index, (self.__class__, d), None

    def __setstate__(self, state):
        """Necessary for making this object picklable"""

        if isinstance(state, dict):
            self._data = state.pop('data')
            for k, v in compat.iteritems(state):
                setattr(self, k, v)

        elif isinstance(state, tuple):

            if len(state) == 2:
                nd_state, own_state = state
                data = np.empty(nd_state[1], dtype=nd_state[2])
                np.ndarray.__setstate__(data, nd_state)
                self.name = own_state[0]

            else:  # pragma: no cover
                data = np.empty(state)
                np.ndarray.__setstate__(data, state)

            self._data = data
            self._reset_identity()
        else:
            raise Exception("invalid pickle state")
    _unpickle_compat = __setstate__

    def __deepcopy__(self, memo={}):
        return self.copy(deep=True)

    def __nonzero__(self):
        raise ValueError("The truth value of a {0} is ambiguous. "
                         "Use a.empty, a.bool(), a.item(), a.any() or a.all()."
                         .format(self.__class__.__name__))

    __bool__ = __nonzero__

    def __contains__(self, key):
        hash(key)
        # work around some kind of odd cython bug
        try:
            return key in self._engine
        except TypeError:
            return False

    def __hash__(self):
        raise TypeError("unhashable type: %r" % type(self).__name__)

    def __setitem__(self, key, value):
        raise TypeError("Index does not support mutable operations")

    def __getitem__(self, key):
        """
        Override numpy.ndarray's __getitem__ method to work as desired.

        This function adds lists and Series as valid boolean indexers
        (ndarrays only supports ndarray with dtype=bool).

        If resulting ndim != 1, plain ndarray is returned instead of
        corresponding `Index` subclass.

        """
        # There's no custom logic to be implemented in __getslice__, so it's
        # not overloaded intentionally.
        getitem = self._data.__getitem__
        promote = self._shallow_copy

        if np.isscalar(key):
            return getitem(key)

        if isinstance(key, slice):
            # This case is separated from the conditional above to avoid
            # pessimization of basic indexing.
            return promote(getitem(key))

        if is_bool_indexer(key):
            key = np.asarray(key)

        key = _values_from_object(key)
        result = getitem(key)
        if not np.isscalar(result):
            return promote(result)
        else:
            return result

    def _ensure_compat_append(self, other):
        """
        prepare the append

        Returns
        -------
        list of to_concat, name of result Index
        """
        name = self.name
        to_concat = [self]

        if isinstance(other, (list, tuple)):
            to_concat = to_concat + list(other)
        else:
            to_concat.append(other)

        for obj in to_concat:
            if (isinstance(obj, Index) and
                obj.name != name and
                obj.name is not None):
                name = None
                break

        to_concat = self._ensure_compat_concat(to_concat)
        to_concat = [x._values if isinstance(x, Index) else x
                     for x in to_concat]
        return to_concat, name

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        to_concat, name = self._ensure_compat_append(other)
        attribs = self._get_attributes_dict()
        attribs['name'] = name
        return self._shallow_copy(np.concatenate(to_concat), infer=True, **attribs)

    @staticmethod
    def _ensure_compat_concat(indexes):
        from pandas.tseries.api import DatetimeIndex, PeriodIndex, TimedeltaIndex
        klasses = DatetimeIndex, PeriodIndex, TimedeltaIndex

        is_ts = [isinstance(idx, klasses) for idx in indexes]

        if any(is_ts) and not all(is_ts):
            return [_maybe_box(idx) for idx in indexes]

        return indexes

    def take(self, indices, axis=0, allow_fill=True, fill_value=None):
        """
        return a new Index of the values selected by the indexer

        For internal compatibility with numpy arrays.

        # filling must always be None/nan here
        # but is passed thru internally

        See also
        --------
        numpy.ndarray.take
        """

        indices = com._ensure_platform_int(indices)
        taken = self.values.take(indices)
        return self._shallow_copy(taken)

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        if self._can_hold_na:
            return isnull(self)
        else:
            # shouldn't reach to this condition by checking hasnans beforehand
            values = np.empty(len(self), dtype=np.bool_)
            values.fill(False)
            return values

    @cache_readonly
    def _nan_idxs(self):
        if self._can_hold_na:
            w, = self._isnan.nonzero()
            return w
        else:
            return np.array([], dtype=np.int64)

    @cache_readonly
    def hasnans(self):
        """ return if I have any nans; enables various perf speedups """
        if self._can_hold_na:
            return self._isnan.any()
        else:
            return False

    def _convert_for_op(self, value):
        """ Convert value to be insertable to ndarray """
        return value

    def _assert_can_do_op(self, value):
        """ Check value is valid for scalar op """
        if not lib.isscalar(value):
            msg = "'value' must be a scalar, passed: {0}"
            raise TypeError(msg.format(type(value).__name__))

    def putmask(self, mask, value):
        """
        return a new Index of the values set with the mask

        See also
        --------
        numpy.ndarray.putmask
        """
        values = self.values.copy()
        try:
            np.putmask(values, mask, self._convert_for_op(value))
            return self._shallow_copy(values)
        except (ValueError, TypeError):
            # coerces to object
            return self.astype(object).putmask(mask, value)

    def format(self, name=False, formatter=None, **kwargs):
        """
        Render a string representation of the Index
        """
        header = []
        if name:
            header.append(com.pprint_thing(self.name,
                                           escape_chars=('\t', '\r', '\n'))
                          if self.name is not None else '')

        if formatter is not None:
            return header + list(self.map(formatter))

        return self._format_with_header(header, **kwargs)

    def _format_with_header(self, header, na_rep='NaN', **kwargs):
        values = self.values

        from pandas.core.format import format_array

        if is_categorical_dtype(values.dtype):
            values = np.array(values)
        elif is_object_dtype(values.dtype):
            values = lib.maybe_convert_objects(values, safe=1)

        if is_object_dtype(values.dtype):
            result = [com.pprint_thing(x, escape_chars=('\t', '\r', '\n'))
                      for x in values]

            # could have nans
            mask = isnull(values)
            if mask.any():
                result = np.array(result)
                result[mask] = na_rep
                result = result.tolist()

        else:
            result = _trim_front(format_array(values, None, justify='left'))
        return header + result

    def to_native_types(self, slicer=None, **kwargs):
        """ slice and dice then format """
        values = self
        if slicer is not None:
            values = values[slicer]
        return values._format_native_types(**kwargs)

    def _format_native_types(self, na_rep='', quoting=None, **kwargs):
        """ actually format my specific types """
        mask = isnull(self)
        if not self.is_object() and not quoting:
            values = np.asarray(self).astype(str)
        else:
            values = np.array(self, dtype=object, copy=True)

        values[mask] = na_rep
        return values

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, Index):
            return False

        return array_equivalent(_values_from_object(self), _values_from_object(other))

    def identical(self, other):
        """Similar to equals, but check that other comparable attributes are
        also equal
        """
        return (self.equals(other) and
                all((getattr(self, c, None) == getattr(other, c, None)
                     for c in self._comparables)) and
                type(self) == type(other))

    def asof(self, label):
        """
        For a sorted index, return the most recent label up to and including
        the passed label. Return NaN if not found.

        See also
        --------
        get_loc : asof is a thin wrapper around get_loc with method='pad'
        """
        try:
            loc = self.get_loc(label, method='pad')
        except KeyError:
            return _get_na_value(self.dtype)
        else:
            if isinstance(loc, slice):
                loc = loc.indices(len(self))[-1]
            return self[loc]

    def asof_locs(self, where, mask):
        """
        where : array of timestamps
        mask : array of booleans where data is not NA

        """
        locs = self.values[mask].searchsorted(where.values, side='right')

        locs = np.where(locs > 0, locs - 1, 0)
        result = np.arange(len(self))[mask].take(locs)

        first = mask.argmax()
        result[(locs == 0) & (where < self.values[first])] = -1

        return result

    def sort_values(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index
        """
        _as = self.argsort()
        if not ascending:
            _as = _as[::-1]

        sorted_index = self.take(_as)

        if return_indexer:
            return sorted_index, _as
        else:
            return sorted_index

    def order(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index

        DEPRECATED: use :meth:`Index.sort_values`
        """
        warnings.warn("order is deprecated, use sort_values(...)",
                      FutureWarning, stacklevel=2)
        return self.sort_values(return_indexer=return_indexer, ascending=ascending)

    def sort(self, *args, **kwargs):
        raise TypeError("cannot sort an Index object in-place, use sort_values instead")

    def sortlevel(self, level=None, ascending=True, sort_remaining=None):
        """

        For internal compatibility with with the Index API

        Sort the Index. This is for compat with MultiIndex

        Parameters
        ----------
        ascending : boolean, default True
            False to sort in descending order

        level, sort_remaining are compat paramaters

        Returns
        -------
        sorted_index : Index
        """
        return self.sort_values(return_indexer=True, ascending=ascending)

    def shift(self, periods=1, freq=None):
        """
        Shift Index containing datetime objects by input number of periods and
        DateOffset

        Returns
        -------
        shifted : Index
        """
        raise NotImplementedError("Not supported for type %s" % type(self).__name__)

    def argsort(self, *args, **kwargs):
        """
        return an ndarray indexer of the underlying data

        See also
        --------
        numpy.ndarray.argsort
        """
        result = self.asi8
        if result is None:
            result = np.array(self)
        return result.argsort(*args, **kwargs)

    def __add__(self, other):
        if com.is_list_like(other):
            warnings.warn("using '+' to provide set union with Indexes is deprecated, "
                          "use '|' or .union()", FutureWarning, stacklevel=2)
        if isinstance(other, Index):
            return self.union(other)
        return Index(np.array(self) + other)

    def __radd__(self, other):
        if is_list_like(other):
            warnings.warn("using '+' to provide set union with Indexes is deprecated, "
                          "use '|' or .union()", FutureWarning, stacklevel=2)
        return Index(other + np.array(self))

    __iadd__ = __add__

    def __sub__(self, other):
        warnings.warn("using '-' to provide set differences with Indexes is deprecated, "
                      "use .difference()",FutureWarning, stacklevel=2)
        return self.difference(other)

    def __and__(self, other):
        return self.intersection(other)

    def __or__(self, other):
        return self.union(other)

    def __xor__(self, other):
        return self.sym_diff(other)

    def union(self, other):
        """
        Form the union of two Index objects and sorts if possible.

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        union : Index

        Examples
        --------

        >>> idx1 = pd.Index([1, 2, 3, 4])
        >>> idx2 = pd.Index([3, 4, 5, 6])
        >>> idx1.union(idx2)
        Int64Index([1, 2, 3, 4, 5, 6], dtype='int64')

        """
        self._assert_can_do_setop(other)
        other = _ensure_index(other)

        if len(other) == 0 or self.equals(other):
            return self

        if len(self) == 0:
            return other

        if not is_dtype_equal(self.dtype,other.dtype):
            this = self.astype('O')
            other = other.astype('O')
            return this.union(other)

        if self.is_monotonic and other.is_monotonic:
            try:
                result = self._outer_indexer(self.values, other._values)[0]
            except TypeError:
                # incomparable objects
                result = list(self.values)

                # worth making this faster? a very unusual case
                value_set = set(self.values)
                result.extend([x for x in other._values if x not in value_set])
        else:
            indexer = self.get_indexer(other)
            indexer, = (indexer == -1).nonzero()

            if len(indexer) > 0:
                other_diff = com.take_nd(other._values, indexer,
                                         allow_fill=False)
                result = com._concat_compat((self.values, other_diff))

                try:
                    self.values[0] < other_diff[0]
                except TypeError as e:
                    warnings.warn("%s, sort order is undefined for "
                                  "incomparable objects" % e,
                                  RuntimeWarning,
                                  stacklevel=3)
                else:
                    types = frozenset((self.inferred_type,
                                       other.inferred_type))
                    if not types & _unsortable_types:
                        result.sort()

            else:
                result = self.values

                try:
                    result = np.sort(result)
                except TypeError as e:
                    warnings.warn("%s, sort order is undefined for "
                                  "incomparable objects" % e,
                                  RuntimeWarning,
                                  stacklevel=3)

        # for subclasses
        return self._wrap_union_result(other, result)

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        return self.__class__(data=result, name=name)

    def intersection(self, other):
        """
        Form the intersection of two Index objects.

        This returns a new Index with elements common to the index and `other`.
        Sortedness of the result is not guaranteed.

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        intersection : Index

        Examples
        --------

        >>> idx1 = pd.Index([1, 2, 3, 4])
        >>> idx2 = pd.Index([3, 4, 5, 6])
        >>> idx1.intersection(idx2)
        Int64Index([3, 4], dtype='int64')

        """
        self._assert_can_do_setop(other)
        other = _ensure_index(other)

        if self.equals(other):
            return self

        if not is_dtype_equal(self.dtype,other.dtype):
            this = self.astype('O')
            other = other.astype('O')
            return this.intersection(other)

        if self.is_monotonic and other.is_monotonic:
            try:
                result = self._inner_indexer(self.values, other._values)[0]
                return self._wrap_union_result(other, result)
            except TypeError:
                pass

        try:
            indexer = Index(self.values).get_indexer(other._values)
            indexer = indexer.take((indexer != -1).nonzero()[0])
        except:
            # duplicates
            indexer = Index(self.values).get_indexer_non_unique(other._values)[0].unique()
            indexer = indexer[indexer != -1]

        taken = self.take(indexer)
        if self.name != other.name:
            taken.name = None
        return taken

    def difference(self, other):
        """
        Return a new Index with elements from the index that are not in `other`.

        This is the sorted set difference of two Index objects.

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        difference : Index

        Examples
        --------

        >>> idx1 = pd.Index([1, 2, 3, 4])
        >>> idx2 = pd.Index([3, 4, 5, 6])
        >>> idx1.difference(idx2)
        Int64Index([1, 2], dtype='int64')

        """
        self._assert_can_do_setop(other)

        if self.equals(other):
            return Index([], name=self.name)

        other, result_name = self._convert_can_do_setop(other)

        theDiff = sorted(set(self) - set(other))
        return Index(theDiff, name=result_name)

    diff = deprecate('diff', difference)

    def sym_diff(self, other, result_name=None):
        """
        Compute the sorted symmetric difference of two Index objects.

        Parameters
        ----------
        other : Index or array-like
        result_name : str

        Returns
        -------
        sym_diff : Index

        Notes
        -----
        ``sym_diff`` contains elements that appear in either ``idx1`` or
        ``idx2`` but not both. Equivalent to the Index created by
        ``(idx1 - idx2) + (idx2 - idx1)`` with duplicates dropped.

        The sorting of a result containing ``NaN`` values is not guaranteed
        across Python versions. See GitHub issue #6444.

        Examples
        --------
        >>> idx1 = Index([1, 2, 3, 4])
        >>> idx2 = Index([2, 3, 4, 5])
        >>> idx1.sym_diff(idx2)
        Int64Index([1, 5], dtype='int64')

        You can also use the ``^`` operator:

        >>> idx1 ^ idx2
        Int64Index([1, 5], dtype='int64')
        """
        self._assert_can_do_setop(other)
        other, result_name_update = self._convert_can_do_setop(other)
        if result_name is None:
            result_name = result_name_update

        the_diff = sorted(set((self.difference(other)).union(other.difference(self))))
        attribs = self._get_attributes_dict()
        attribs['name'] = result_name
        if 'freq' in attribs:
            attribs['freq'] = None
        return self._shallow_copy(the_diff, infer=True, **attribs)

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Parameters
        ----------
        key : label
        method : {None, 'pad'/'ffill', 'backfill'/'bfill', 'nearest'}, optional
            * default: exact matches only.
            * pad / ffill: find the PREVIOUS index value if no exact match.
            * backfill / bfill: use NEXT index value if no exact match
            * nearest: use the NEAREST index value if no exact match. Tied
              distances are broken by preferring the larger index value.
        tolerance : optional
            Maximum distance from index value for inexact matches. The value of
            the index at the matching location most satisfy the equation
            ``abs(index[loc] - key) <= tolerance``.

            .. versionadded:: 0.17.0

        Returns
        -------
        loc : int if unique index, possibly slice or mask if not
        """
        if method is None:
            if tolerance is not None:
                raise ValueError('tolerance argument only valid if using pad, '
                                 'backfill or nearest lookups')
            key = _values_from_object(key)
            return self._engine.get_loc(key)

        indexer = self.get_indexer([key], method=method,
                                   tolerance=tolerance)
        if indexer.ndim > 1 or indexer.size > 1:
            raise TypeError('get_loc requires scalar valued input')
        loc = indexer.item()
        if loc == -1:
            raise KeyError(key)
        return loc

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """

        # if we have something that is Index-like, then
        # use this, e.g. DatetimeIndex
        s = getattr(series,'_values',None)
        if isinstance(s, Index) and lib.isscalar(key):
            return s[key]

        s = _values_from_object(series)
        k = _values_from_object(key)

        # prevent integer truncation bug in indexing
        if is_float(k) and not self.is_floating():
            raise KeyError

        try:
            return self._engine.get_value(s, k)
        except KeyError as e1:
            if len(self) > 0 and self.inferred_type in ['integer','boolean']:
                raise

            try:
                return tslib.get_value_box(s, key)
            except IndexError:
                raise
            except TypeError:
                # generator/iterator-like
                if is_iterator(key):
                    raise InvalidIndexError(key)
                else:
                    raise e1
            except Exception:  # pragma: no cover
                raise e1
        except TypeError:
            # python 3
            if np.isscalar(key):  # pragma: no cover
                raise IndexError(key)
            raise InvalidIndexError(key)

    def set_value(self, arr, key, value):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        self._engine.set_value(
            _values_from_object(arr), _values_from_object(key), value)

    def get_level_values(self, level):
        """
        Return vector of label values for requested level, equal to the length
        of the index

        Parameters
        ----------
        level : int

        Returns
        -------
        values : ndarray
        """
        # checks that level number is actually just 1
        self._validate_index_level(level)
        return self

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        """
        Compute indexer and mask for new index given the current index. The
        indexer should be then used as an input to ndarray.take to align the
        current data to the new index.

        Parameters
        ----------
        target : Index
        method : {None, 'pad'/'ffill', 'backfill'/'bfill', 'nearest'}, optional
            * default: exact matches only.
            * pad / ffill: find the PREVIOUS index value if no exact match.
            * backfill / bfill: use NEXT index value if no exact match
            * nearest: use the NEAREST index value if no exact match. Tied
              distances are broken by preferring the larger index value.
        limit : int, optional
            Maximum number of consecutive labels in ``target`` to match for
            inexact matches.
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            .. versionadded:: 0.17.0

        Examples
        --------
        >>> indexer = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)

        Returns
        -------
        indexer : ndarray of int
            Integers from 0 to n - 1 indicating that the index at these
            positions matches the corresponding target values. Missing values
            in the target are marked by -1.
        """
        method = _clean_reindex_fill_method(method)
        target = _ensure_index(target)
        if tolerance is not None:
            tolerance = self._convert_tolerance(tolerance)

        pself, ptarget = self._possibly_promote(target)
        if pself is not self or ptarget is not target:
            return pself.get_indexer(ptarget, method=method, limit=limit,
                                     tolerance=tolerance)

        if not is_dtype_equal(self.dtype, target.dtype):
            this = self.astype(object)
            target = target.astype(object)
            return this.get_indexer(target, method=method, limit=limit,
                                    tolerance=tolerance)

        if not self.is_unique:
            raise InvalidIndexError('Reindexing only valid with uniquely'
                                    ' valued Index objects')

        if method == 'pad' or method == 'backfill':
            indexer = self._get_fill_indexer(target, method, limit, tolerance)
        elif method == 'nearest':
            indexer = self._get_nearest_indexer(target, limit, tolerance)
        else:
            if tolerance is not None:
                raise ValueError('tolerance argument only valid if doing pad, '
                                 'backfill or nearest reindexing')
            if limit is not None:
                raise ValueError('limit argument only valid if doing pad, '
                                 'backfill or nearest reindexing')

            indexer = self._engine.get_indexer(target._values)

        return com._ensure_platform_int(indexer)

    def _convert_tolerance(self, tolerance):
        # override this method on subclasses
        return tolerance

    def _get_fill_indexer(self, target, method, limit=None, tolerance=None):
        if self.is_monotonic_increasing and target.is_monotonic_increasing:
            method = (self._engine.get_pad_indexer if method == 'pad'
                      else self._engine.get_backfill_indexer)
            indexer = method(target._values, limit)
        else:
            indexer = self._get_fill_indexer_searchsorted(target, method, limit)
        if tolerance is not None:
            indexer = self._filter_indexer_tolerance(
                target._values, indexer, tolerance)
        return indexer

    def _get_fill_indexer_searchsorted(self, target, method, limit=None):
        """
        Fallback pad/backfill get_indexer that works for monotonic decreasing
        indexes and non-monotonic targets
        """
        if limit is not None:
            raise ValueError('limit argument for %r method only well-defined '
                             'if index and target are monotonic' % method)

        side = 'left' if method == 'pad' else 'right'
        target = np.asarray(target)

        # find exact matches first (this simplifies the algorithm)
        indexer = self.get_indexer(target)
        nonexact = (indexer == -1)
        indexer[nonexact] = self._searchsorted_monotonic(target[nonexact], side)
        if side == 'left':
            # searchsorted returns "indices into a sorted array such that,
            # if the corresponding elements in v were inserted before the
            # indices, the order of a would be preserved".
            # Thus, we need to subtract 1 to find values to the left.
            indexer[nonexact] -= 1
            # This also mapped not found values (values of 0 from
            # np.searchsorted) to -1, which conveniently is also our
            # sentinel for missing values
        else:
            # Mark indices to the right of the largest value as not found
            indexer[indexer == len(self)] = -1
        return indexer

    def _get_nearest_indexer(self, target, limit, tolerance):
        """
        Get the indexer for the nearest index labels; requires an index with
        values that can be subtracted from each other (e.g., not strings or
        tuples).
        """
        left_indexer = self.get_indexer(target, 'pad', limit=limit)
        right_indexer = self.get_indexer(target, 'backfill', limit=limit)

        target = np.asarray(target)
        left_distances = abs(self.values[left_indexer] - target)
        right_distances = abs(self.values[right_indexer] - target)

        op = operator.lt if self.is_monotonic_increasing else operator.le
        indexer = np.where(op(left_distances, right_distances)
                           | (right_indexer == -1),
                           left_indexer, right_indexer)
        if tolerance is not None:
            indexer = self._filter_indexer_tolerance(
                target, indexer, tolerance)
        return indexer

    def _filter_indexer_tolerance(self, target, indexer, tolerance):
        distance = abs(self.values[indexer] - target)
        indexer = np.where(distance <= tolerance, indexer, -1)
        return indexer

    def get_indexer_non_unique(self, target):
        """ return an indexer suitable for taking from a non unique index
            return the labels in the same order as the target, and
            return a missing indexer into the target (missing are marked as -1
            in the indexer); target must be an iterable """
        target = _ensure_index(target)
        pself, ptarget = self._possibly_promote(target)
        if pself is not self or ptarget is not target:
            return pself.get_indexer_non_unique(ptarget)

        if self.is_all_dates:
            self = Index(self.asi8)
            tgt_values = target.asi8
        else:
            tgt_values = target._values

        indexer, missing = self._engine.get_indexer_non_unique(tgt_values)
        return Index(indexer), missing

    def get_indexer_for(self, target, **kwargs):
        """ guaranteed return of an indexer even when non-unique """
        if self.is_unique:
            return self.get_indexer(target, **kwargs)
        indexer, _ = self.get_indexer_non_unique(target, **kwargs)
        return indexer

    def _possibly_promote(self, other):
        # A hack, but it works
        from pandas.tseries.index import DatetimeIndex
        if self.inferred_type == 'date' and isinstance(other, DatetimeIndex):
            return DatetimeIndex(self), other
        elif self.inferred_type == 'boolean':
            if not is_object_dtype(self.dtype):
                return self.astype('object'), other.astype('object')
        return self, other

    def groupby(self, to_groupby):
        """
        Group the index labels by a given array of values.

        Parameters
        ----------
        to_groupby : array
            Values used to determine the groups.

        Returns
        -------
        groups : dict
            {group name -> group labels}

        """
        return self._groupby(self.values, _values_from_object(to_groupby))

    def map(self, mapper):
        return self._arrmap(self.values, mapper)

    def isin(self, values, level=None):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values.

        Parameters
        ----------
        values : set or sequence of values
            Sought values.
        level : str or int, optional
            Name or position of the index level to use (if the index is a
            MultiIndex).

        Notes
        -----
        If `level` is specified:

        - if it is the name of one *and only one* index level, use that level;
        - otherwise it should be a number indicating level position.

        Returns
        -------
        is_contained : ndarray (boolean dtype)

        """
        if level is not None:
            self._validate_index_level(level)
        return algorithms.isin(np.array(self), values)

    def _can_reindex(self, indexer):
        """
        *this is an internal non-public method*

        Check if we are allowing reindexing with this particular indexer

        Parameters
        ----------
        indexer : an integer indexer

        Raises
        ------
        ValueError if its a duplicate axis
        """

        # trying to reindex on an axis with duplicates
        if not self.is_unique and len(indexer):
            raise ValueError("cannot reindex from a duplicate axis")

    def reindex(self, target, method=None, level=None, limit=None,
                tolerance=None):
        """
        Create index with target's values (move/add/delete values as necessary)

        Parameters
        ----------
        target : an iterable

        Returns
        -------
        new_index : pd.Index
            Resulting index
        indexer : np.ndarray or None
            Indices of output values in original index

        """
        # GH6552: preserve names when reindexing to non-named target
        # (i.e. neither Index nor Series).
        preserve_names = not hasattr(target, 'name')

        # GH7774: preserve dtype/tz if target is empty and not an Index.
        target = _ensure_has_len(target)  # target may be an iterator
        if not isinstance(target, Index) and len(target) == 0:
            attrs = self._get_attributes_dict()
            attrs.pop('freq', None)  # don't preserve freq
            target = self._simple_new(None, dtype=self.dtype, **attrs)
        else:
            target = _ensure_index(target)

        if level is not None:
            if method is not None:
                raise TypeError('Fill method not supported if level passed')
            _, indexer, _ = self._join_level(target, level, how='right',
                                             return_indexers=True)
        else:
            if self.equals(target):
                indexer = None
            else:
                if self.is_unique:
                    indexer = self.get_indexer(target, method=method,
                                               limit=limit,
                                               tolerance=tolerance)
                else:
                    if method is not None or limit is not None:
                        raise ValueError("cannot reindex a non-unique index "
                                         "with a method or limit")
                    indexer, missing = self.get_indexer_non_unique(target)

        if preserve_names and target.nlevels == 1 and target.name != self.name:
            target = target.copy()
            target.name = self.name

        return target, indexer

    def _reindex_non_unique(self, target):
        """
        *this is an internal non-public method*

        Create a new index with target's values (move/add/delete values as necessary)
        use with non-unique Index and a possibly non-unique target

        Parameters
        ----------
        target : an iterable

        Returns
        -------
        new_index : pd.Index
            Resulting index
        indexer : np.ndarray or None
            Indices of output values in original index

        """

        target = _ensure_index(target)
        indexer, missing = self.get_indexer_non_unique(target)
        check = indexer != -1
        new_labels = self.take(indexer[check])
        new_indexer = None

        if len(missing):
            l = np.arange(len(indexer))

            missing = com._ensure_platform_int(missing)
            missing_labels = target.take(missing)
            missing_indexer = com._ensure_int64(l[~check])
            cur_labels = self.take(indexer[check])._values
            cur_indexer = com._ensure_int64(l[check])

            new_labels = np.empty(tuple([len(indexer)]), dtype=object)
            new_labels[cur_indexer] = cur_labels
            new_labels[missing_indexer] = missing_labels

            # a unique indexer
            if target.is_unique:

                # see GH5553, make sure we use the right indexer
                new_indexer = np.arange(len(indexer))
                new_indexer[cur_indexer] = np.arange(len(cur_labels))
                new_indexer[missing_indexer] = -1

            # we have a non_unique selector, need to use the original
            # indexer here
            else:

                # need to retake to have the same size as the indexer
                indexer = indexer._values
                indexer[~check] = 0

                # reset the new indexer to account for the new size
                new_indexer = np.arange(len(self.take(indexer)))
                new_indexer[~check] = -1

        return self._shallow_copy(new_labels), indexer, new_indexer

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        *this is an internal non-public method*

        Compute join_index and indexers to conform data
        structures to the new index.

        Parameters
        ----------
        other : Index
        how : {'left', 'right', 'inner', 'outer'}
        level : int or level name, default None
        return_indexers : boolean, default False

        Returns
        -------
        join_index, (left_indexer, right_indexer)
        """
        self_is_mi = isinstance(self, MultiIndex)
        other_is_mi = isinstance(other, MultiIndex)

        # try to figure out the join level
        # GH3662
        if (level is None and (self_is_mi or other_is_mi)):

            # have the same levels/names so a simple join
            if self.names == other.names:
                pass
            else:
                return self._join_multi(other, how=how, return_indexers=return_indexers)

        # join on the level
        if (level is not None and (self_is_mi or other_is_mi)):
            return self._join_level(other, level, how=how,
                                    return_indexers=return_indexers)

        other = _ensure_index(other)

        if len(other) == 0 and how in ('left', 'outer'):
            join_index = self._shallow_copy()
            if return_indexers:
                rindexer = np.repeat(-1, len(join_index))
                return join_index, None, rindexer
            else:
                return join_index

        if len(self) == 0 and how in ('right', 'outer'):
            join_index = other._shallow_copy()
            if return_indexers:
                lindexer = np.repeat(-1, len(join_index))
                return join_index, lindexer, None
            else:
                return join_index

        if self._join_precedence < other._join_precedence:
            how = {'right': 'left', 'left': 'right'}.get(how, how)
            result = other.join(self, how=how, level=level,
                                return_indexers=return_indexers)
            if return_indexers:
                x, y, z = result
                result = x, z, y
            return result

        if not is_dtype_equal(self.dtype,other.dtype):
            this = self.astype('O')
            other = other.astype('O')
            return this.join(other, how=how,
                             return_indexers=return_indexers)

        _validate_join_method(how)

        if not self.is_unique and not other.is_unique:
            return self._join_non_unique(other, how=how,
                                         return_indexers=return_indexers)
        elif not self.is_unique or not other.is_unique:
            if self.is_monotonic and other.is_monotonic:
                return self._join_monotonic(other, how=how,
                                            return_indexers=return_indexers)
            else:
                return self._join_non_unique(other, how=how,
                                             return_indexers=return_indexers)
        elif self.is_monotonic and other.is_monotonic:
            try:
                return self._join_monotonic(other, how=how,
                                            return_indexers=return_indexers)
            except TypeError:
                pass

        if how == 'left':
            join_index = self
        elif how == 'right':
            join_index = other
        elif how == 'inner':
            join_index = self.intersection(other)
        elif how == 'outer':
            join_index = self.union(other)

        if return_indexers:
            if join_index is self:
                lindexer = None
            else:
                lindexer = self.get_indexer(join_index)
            if join_index is other:
                rindexer = None
            else:
                rindexer = other.get_indexer(join_index)
            return join_index, lindexer, rindexer
        else:
            return join_index

    def _join_multi(self, other, how, return_indexers=True):

        self_is_mi = isinstance(self, MultiIndex)
        other_is_mi = isinstance(other, MultiIndex)

        # figure out join names
        self_names = [ n for n in self.names if n is not None ]
        other_names = [ n for n in other.names if n is not None ]
        overlap = list(set(self_names) & set(other_names))

        # need at least 1 in common, but not more than 1
        if not len(overlap):
            raise ValueError("cannot join with no level specified and no overlapping names")
        if len(overlap) > 1:
            raise NotImplementedError("merging with more than one level overlap on a multi-index is not implemented")
        jl = overlap[0]

        # make the indices into mi's that match
        if not (self_is_mi and other_is_mi):

            flip_order = False
            if self_is_mi:
                self, other = other, self
                flip_order = True
                # flip if join method is right or left
                how = {'right': 'left', 'left': 'right'}.get(how, how)

            level = other.names.index(jl)
            result = self._join_level(other, level, how=how,
                                      return_indexers=return_indexers)

            if flip_order:
                if isinstance(result, tuple):
                    return result[0], result[2], result[1]
            return result

        # 2 multi-indexes
        raise NotImplementedError("merging with both multi-indexes is not implemented")

    def _join_non_unique(self, other, how='left', return_indexers=False):
        from pandas.tools.merge import _get_join_indexers

        left_idx, right_idx = _get_join_indexers([self.values], [other._values],
                                                 how=how, sort=True)

        left_idx = com._ensure_platform_int(left_idx)
        right_idx = com._ensure_platform_int(right_idx)

        join_index = self.values.take(left_idx)
        mask = left_idx == -1
        np.putmask(join_index, mask, other._values.take(right_idx))

        join_index = self._wrap_joined_index(join_index, other)

        if return_indexers:
            return join_index, left_idx, right_idx
        else:
            return join_index

    def _join_level(self, other, level, how='left',
                    return_indexers=False,
                    keep_order=True):
        """
        The join method *only* affects the level of the resulting
        MultiIndex. Otherwise it just exactly aligns the Index data to the
        labels of the level in the MultiIndex. If `keep_order` == True, the
        order of the data indexed by the MultiIndex will not be changed;
        otherwise, it will tie out with `other`.
        """
        from pandas.algos import groupsort_indexer

        def _get_leaf_sorter(labels):
            '''
            returns sorter for the inner most level while preserving the
            order of higher levels
            '''
            if labels[0].size == 0:
                return np.empty(0, dtype='int64')

            if len(labels) == 1:
                lab = com._ensure_int64(labels[0])
                sorter, _ = groupsort_indexer(lab, 1 + lab.max())
                return sorter

            # find indexers of begining of each set of
            # same-key labels w.r.t all but last level
            tic = labels[0][:-1] != labels[0][1:]
            for lab in labels[1:-1]:
                tic |= lab[:-1] != lab[1:]

            starts = np.hstack(([True], tic, [True])).nonzero()[0]
            lab = com._ensure_int64(labels[-1])
            return lib.get_level_sorter(lab, com._ensure_int64(starts))

        if isinstance(self, MultiIndex) and isinstance(other, MultiIndex):
            raise TypeError('Join on level between two MultiIndex objects '
                            'is ambiguous')

        left, right = self, other

        flip_order = not isinstance(self, MultiIndex)
        if flip_order:
            left, right = right, left
            how = {'right': 'left', 'left': 'right'}.get(how, how)

        level = left._get_level_number(level)
        old_level = left.levels[level]

        if not right.is_unique:
            raise NotImplementedError('Index._join_level on non-unique index '
                                      'is not implemented')

        new_level, left_lev_indexer, right_lev_indexer = \
            old_level.join(right, how=how, return_indexers=True)

        if left_lev_indexer is None:
            if keep_order or len(left) == 0:
                left_indexer = None
                join_index = left
            else:  # sort the leaves
                left_indexer = _get_leaf_sorter(left.labels[:level + 1])
                join_index = left[left_indexer]

        else:
            left_lev_indexer = com._ensure_int64(left_lev_indexer)
            rev_indexer = lib.get_reverse_indexer(left_lev_indexer,
                                                  len(old_level))

            new_lev_labels = com.take_nd(rev_indexer, left.labels[level],
                                         allow_fill=False)

            new_labels = list(left.labels)
            new_labels[level] = new_lev_labels

            new_levels = list(left.levels)
            new_levels[level] = new_level

            if keep_order:  # just drop missing values. o.w. keep order
                left_indexer = np.arange(len(left))
                mask = new_lev_labels != -1
                if not mask.all():
                    new_labels = [lab[mask] for lab in new_labels]
                    left_indexer = left_indexer[mask]

            else:  # tie out the order with other
                if level == 0:  # outer most level, take the fast route
                    ngroups = 1 + new_lev_labels.max()
                    left_indexer, counts = groupsort_indexer(new_lev_labels,
                                                             ngroups)
                    # missing values are placed first; drop them!
                    left_indexer = left_indexer[counts[0]:]
                    new_labels = [lab[left_indexer] for lab in new_labels]

                else:  # sort the leaves
                    mask = new_lev_labels != -1
                    mask_all = mask.all()
                    if not mask_all:
                        new_labels = [lab[mask] for lab in new_labels]

                    left_indexer = _get_leaf_sorter(new_labels[:level + 1])
                    new_labels = [lab[left_indexer] for lab in new_labels]

                    # left_indexers are w.r.t masked frame.
                    # reverse to original frame!
                    if not mask_all:
                        left_indexer = mask.nonzero()[0][left_indexer]

            join_index = MultiIndex(levels=new_levels,
                                    labels=new_labels,
                                    names=left.names,
                                    verify_integrity=False)

        if right_lev_indexer is not None:
            right_indexer = com.take_nd(right_lev_indexer,
                                        join_index.labels[level],
                                        allow_fill=False)
        else:
            right_indexer = join_index.labels[level]

        if flip_order:
            left_indexer, right_indexer = right_indexer, left_indexer

        if return_indexers:
            return join_index, left_indexer, right_indexer
        else:
            return join_index

    def _join_monotonic(self, other, how='left', return_indexers=False):
        if self.equals(other):
            ret_index = other if how == 'right' else self
            if return_indexers:
                return ret_index, None, None
            else:
                return ret_index

        sv = self.values
        ov = other._values

        if self.is_unique and other.is_unique:
            # We can perform much better than the general case
            if how == 'left':
                join_index = self
                lidx = None
                ridx = self._left_indexer_unique(sv, ov)
            elif how == 'right':
                join_index = other
                lidx = self._left_indexer_unique(ov, sv)
                ridx = None
            elif how == 'inner':
                join_index, lidx, ridx = self._inner_indexer(sv, ov)
                join_index = self._wrap_joined_index(join_index, other)
            elif how == 'outer':
                join_index, lidx, ridx = self._outer_indexer(sv, ov)
                join_index = self._wrap_joined_index(join_index, other)
        else:
            if how == 'left':
                join_index, lidx, ridx = self._left_indexer(sv, ov)
            elif how == 'right':
                join_index, ridx, lidx = self._left_indexer(ov, sv)
            elif how == 'inner':
                join_index, lidx, ridx = self._inner_indexer(sv, ov)
            elif how == 'outer':
                join_index, lidx, ridx = self._outer_indexer(sv, ov)
            join_index = self._wrap_joined_index(join_index, other)

        if return_indexers:
            return join_index, lidx, ridx
        else:
            return join_index

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        return Index(joined, name=name)

    def slice_indexer(self, start=None, end=None, step=None, kind=None):
        """
        For an ordered Index, compute the slice indexer for input labels and
        step

        Parameters
        ----------
        start : label, default None
            If None, defaults to the beginning
        end : label, default None
            If None, defaults to the end
        step : int, default None
        kind : string, default None

        Returns
        -------
        indexer : ndarray or slice

        Notes
        -----
        This function assumes that the data is sorted, so use at your own peril
        """
        start_slice, end_slice = self.slice_locs(start, end, step=step, kind=kind)

        # return a slice
        if not lib.isscalar(start_slice):
            raise AssertionError("Start slice bound is non-scalar")
        if not lib.isscalar(end_slice):
            raise AssertionError("End slice bound is non-scalar")

        return slice(start_slice, end_slice, step)

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        This function should be overloaded in subclasses that allow non-trivial
        casting on label-slice bounds, e.g. datetime-like indices allowing
        strings containing formatted datetimes.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : string / None

        Returns
        -------
        label :  object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """

        # We are a plain index here (sub-class override this method if they
        # wish to have special treatment for floats/ints, e.g. Float64Index and
        # datetimelike Indexes
        # reject them
        if is_float(label):
            self._invalid_indexer('slice',label)

        # we are trying to find integer bounds on a non-integer based index
        # this is rejected (generally .loc gets you here)
        elif is_integer(label):
            self._invalid_indexer('slice',label)

        return label

    def _searchsorted_monotonic(self, label, side='left'):
        if self.is_monotonic_increasing:
            return self.searchsorted(label, side=side)
        elif self.is_monotonic_decreasing:
            # np.searchsorted expects ascending sort order, have to reverse
            # everything for it to work (element ordering, search side and
            # resulting value).
            pos = self[::-1].searchsorted(
                label, side='right' if side == 'left' else 'right')
            return len(self) - pos

        raise ValueError('index must be monotonic increasing or decreasing')

    def get_slice_bound(self, label, side, kind):
        """
        Calculate slice bound that corresponds to given label.

        Returns leftmost (one-past-the-rightmost if ``side=='right'``) position
        of given label.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : string / None, the type of indexer

        """
        if side not in ('left', 'right'):
            raise ValueError(
                "Invalid value for side kwarg,"
                " must be either 'left' or 'right': %s" % (side,))

        original_label = label

        # For datetime indices label may be a string that has to be converted
        # to datetime boundary according to its resolution.
        label = self._maybe_cast_slice_bound(label, side, kind)

        # we need to look up the label
        try:
            slc = self.get_loc(label)
        except KeyError as err:
            try:
                return self._searchsorted_monotonic(label, side)
            except ValueError:
                # raise the original KeyError
                raise err

        if isinstance(slc, np.ndarray):
            # get_loc may return a boolean array or an array of indices, which
            # is OK as long as they are representable by a slice.
            if is_bool_dtype(slc):
                slc = lib.maybe_booleans_to_slice(slc.view('u1'))
            else:
                slc = lib.maybe_indices_to_slice(slc.astype('i8'), len(self))
            if isinstance(slc, np.ndarray):
                raise KeyError(
                    "Cannot get %s slice bound for non-unique label:"
                    " %r" % (side, original_label))

        if isinstance(slc, slice):
            if side == 'left':
                return slc.start
            else:
                return slc.stop
        else:
            if side == 'right':
                return slc + 1
            else:
                return slc

    def slice_locs(self, start=None, end=None, step=None, kind=None):
        """
        Compute slice locations for input labels.

        Parameters
        ----------
        start : label, default None
            If None, defaults to the beginning
        end : label, default None
            If None, defaults to the end
        step : int, defaults None
            If None, defaults to 1
        kind : string, defaults None

        Returns
        -------
        start, end : int

        """
        inc = (step is None or step >= 0)

        if not inc:
            # If it's a reverse slice, temporarily swap bounds.
            start, end = end, start

        start_slice = None
        if start is not None:
            start_slice = self.get_slice_bound(start, 'left', kind)
        if start_slice is None:
            start_slice = 0

        end_slice = None
        if end is not None:
            end_slice = self.get_slice_bound(end, 'right', kind)
        if end_slice is None:
            end_slice = len(self)

        if not inc:
            # Bounds at this moment are swapped, swap them back and shift by 1.
            #
            # slice_locs('B', 'A', step=-1): s='B', e='A'
            #
            #              s='A'                 e='B'
            # AFTER SWAP:    |                     |
            #                v ------------------> V
            #           -----------------------------------
            #           | | |A|A|A|A| | | | | |B|B| | | | |
            #           -----------------------------------
            #              ^ <------------------ ^
            # SHOULD BE:   |                     |
            #           end=s-1              start=e-1
            #
            end_slice, start_slice = start_slice - 1, end_slice - 1

            # i == -1 triggers ``len(self) + i`` selection that points to the
            # last element, not before-the-first one, subtracting len(self)
            # compensates that.
            if end_slice == -1:
                end_slice -= len(self)
            if start_slice == -1:
                start_slice -= len(self)

        return start_slice, end_slice

    def delete(self, loc):
        """
        Make new Index with passed location(-s) deleted

        Returns
        -------
        new_index : Index
        """
        attribs = self._get_attributes_dict()
        return self._shallow_copy(np.delete(self._data, loc), **attribs)

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location. Follows
        Python list.append semantics for negative values

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : Index
        """
        _self = np.asarray(self)
        item = self._coerce_scalar_to_index(item)._values

        idx = np.concatenate(
            (_self[:loc], item, _self[loc:]))
        attribs = self._get_attributes_dict()
        return self._shallow_copy(idx, infer=True, **attribs)

    def drop(self, labels, errors='raise'):
        """
        Make new Index with passed list of labels deleted

        Parameters
        ----------
        labels : array-like
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and existing labels are dropped.

        Returns
        -------
        dropped : Index
        """
        labels = com._index_labels_to_array(labels)
        indexer = self.get_indexer(labels)
        mask = indexer == -1
        if mask.any():
            if errors != 'ignore':
                raise ValueError('labels %s not contained in axis' % labels[mask])
            indexer = indexer[~mask]
        return self.delete(indexer)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(base._shared_docs['drop_duplicates'] % _index_doc_kwargs)
    def drop_duplicates(self, keep='first'):
        return super(Index, self).drop_duplicates(keep=keep)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(base._shared_docs['duplicated'] % _index_doc_kwargs)
    def duplicated(self, keep='first'):
        return super(Index, self).duplicated(keep=keep)

    _index_shared_docs['fillna'] = """
        Fill NA/NaN values with the specified value

        Parameters
        ----------
        value : scalar
            Scalar value to use to fill holes (e.g. 0).
            This value cannot be a list-likes.
        downcast : dict, default is None
            a dict of item->dtype of what to downcast if possible,
            or the string 'infer' which will try to downcast to an appropriate
            equal type (e.g. float64 to int64 if possible)

        Returns
        -------
        filled : Index
        """

    @Appender(_index_shared_docs['fillna'])
    def fillna(self, value=None, downcast=None):
        self._assert_can_do_op(value)
        if self.hasnans:
            result = self.putmask(self._isnan, value)
            if downcast is None:
                # no need to care metadata other than name
                # because it can't have freq if
                return Index(result, name=self.name)

        return self._shallow_copy()

    def _evaluate_with_timedelta_like(self, other, op, opstr):
        raise TypeError("can only perform ops with timedelta like values")

    def _evaluate_with_datetime_like(self, other, op, opstr):
        raise TypeError("can only perform ops with datetime like values")

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """

        def _make_compare(op):

            def _evaluate_compare(self, other):
                if isinstance(other, (np.ndarray, Index, ABCSeries)):
                    if other.ndim > 0 and len(self) != len(other):
                        raise ValueError('Lengths must match to compare')
                func = getattr(self.values, op)
                result = func(np.asarray(other))

                # technically we could support bool dtyped Index
                # for now just return the indexing array directly
                if is_bool_dtype(result):
                    return result
                try:
                    return Index(result)
                except TypeError:
                    return result

            return _evaluate_compare

        cls.__eq__ = _make_compare('__eq__')
        cls.__ne__ = _make_compare('__ne__')
        cls.__lt__ = _make_compare('__lt__')
        cls.__gt__ = _make_compare('__gt__')
        cls.__le__ = _make_compare('__le__')
        cls.__ge__ = _make_compare('__ge__')

    @classmethod
    def _add_numericlike_set_methods_disabled(cls):
        """ add in the numeric set-like methods to disable """

        def _make_invalid_op(name):

            def invalid_op(self, other=None):
                raise TypeError("cannot perform {name} with this index type: {typ}".format(name=name,
                                                                                           typ=type(self)))
            invalid_op.__name__ = name
            return invalid_op

        cls.__add__ = cls.__radd__ = __iadd__ = _make_invalid_op('__add__')
        cls.__sub__ = __isub__ = _make_invalid_op('__sub__')

    @classmethod
    def _add_numeric_methods_disabled(cls):
        """ add in numeric methods to disable """

        def _make_invalid_op(name):

            def invalid_op(self, other=None):
                raise TypeError("cannot perform {name} with this index type: {typ}".format(name=name,
                                                                                           typ=type(self)))
            invalid_op.__name__ = name
            return invalid_op

        cls.__mul__ = cls.__rmul__ = _make_invalid_op('__mul__')
        cls.__floordiv__ = cls.__rfloordiv__ = _make_invalid_op('__floordiv__')
        cls.__truediv__ = cls.__rtruediv__ = _make_invalid_op('__truediv__')
        if not compat.PY3:
            cls.__div__ = cls.__rdiv__ = _make_invalid_op('__div__')
        cls.__neg__ = _make_invalid_op('__neg__')
        cls.__pos__ = _make_invalid_op('__pos__')
        cls.__abs__ = _make_invalid_op('__abs__')
        cls.__inv__ = _make_invalid_op('__inv__')

    def _maybe_update_attributes(self, attrs):
        """ Update Index attributes (e.g. freq) depending on op """
        return attrs

    @classmethod
    def _add_numeric_methods(cls):
        """ add in numeric methods """

        def _make_evaluate_binop(op, opstr, reversed=False):

            def _evaluate_numeric_binop(self, other):
                import pandas.tseries.offsets as offsets

                # if we are an inheritor of numeric, but not actually numeric (e.g. DatetimeIndex/PeriodInde)
                if not self._is_numeric_dtype:
                    raise TypeError("cannot evaluate a numeric op {opstr} for type: {typ}".format(opstr=opstr,
                                                                                                  typ=type(self)))

                if isinstance(other, Index):
                    if not other._is_numeric_dtype:
                        raise TypeError("cannot evaluate a numeric op {opstr} with type: {typ}".format(opstr=type(self),
                                                                                                       typ=type(other)))
                elif isinstance(other, np.ndarray) and not other.ndim:
                    other = other.item()

                if isinstance(other, (Index, ABCSeries, np.ndarray)):
                    if len(self) != len(other):
                        raise ValueError("cannot evaluate a numeric op with unequal lengths")
                    other = _values_from_object(other)
                    if other.dtype.kind not in ['f','i']:
                        raise TypeError("cannot evaluate a numeric op with a non-numeric dtype")
                elif isinstance(other, (offsets.DateOffset, np.timedelta64, Timedelta, datetime.timedelta)):
                    return self._evaluate_with_timedelta_like(other, op, opstr)
                elif isinstance(other, (Timestamp, np.datetime64)):
                    return self._evaluate_with_datetime_like(other, op, opstr)
                else:
                    if not (is_float(other) or is_integer(other)):
                        raise TypeError("can only perform ops with scalar values")

                # if we are a reversed non-communative op
                values = self.values
                if reversed:
                    values, other = other, values

                attrs = self._get_attributes_dict()
                attrs = self._maybe_update_attributes(attrs)
                return Index(op(values, other), **attrs)

            return _evaluate_numeric_binop

        def _make_evaluate_unary(op, opstr):

            def _evaluate_numeric_unary(self):

                # if we are an inheritor of numeric, but not actually numeric (e.g. DatetimeIndex/PeriodInde)
                if not self._is_numeric_dtype:
                    raise TypeError("cannot evaluate a numeric op {opstr} for type: {typ}".format(opstr=opstr,
                                                                                                  typ=type(self)))
                attrs = self._get_attributes_dict()
                attrs = self._maybe_update_attributes(attrs)
                return Index(op(self.values), **attrs)

            return _evaluate_numeric_unary

        cls.__add__ = cls.__radd__ = _make_evaluate_binop(operator.add,'__add__')
        cls.__sub__ = _make_evaluate_binop(operator.sub,'__sub__')
        cls.__rsub__ = _make_evaluate_binop(operator.sub,'__sub__',reversed=True)
        cls.__mul__ = cls.__rmul__ = _make_evaluate_binop(operator.mul,'__mul__')
        cls.__floordiv__ = _make_evaluate_binop(operator.floordiv,'__floordiv__')
        cls.__rfloordiv__ = _make_evaluate_binop(operator.floordiv,'__floordiv__',reversed=True)
        cls.__truediv__ = _make_evaluate_binop(operator.truediv,'__truediv__')
        cls.__rtruediv__ = _make_evaluate_binop(operator.truediv,'__truediv__',reversed=True)
        if not compat.PY3:
            cls.__div__ = _make_evaluate_binop(operator.div,'__div__')
            cls.__rdiv__ = _make_evaluate_binop(operator.div,'__div__',reversed=True)
        cls.__neg__ = _make_evaluate_unary(lambda x: -x,'__neg__')
        cls.__pos__ = _make_evaluate_unary(lambda x: x,'__pos__')
        cls.__abs__ = _make_evaluate_unary(lambda x: np.abs(x),'__abs__')
        cls.__inv__ = _make_evaluate_unary(lambda x: -x,'__inv__')

    @classmethod
    def _add_logical_methods(cls):
        """ add in logical methods """

        _doc = """

        %(desc)s

        Parameters
        ----------
        All arguments to numpy.%(outname)s are accepted.

        Returns
        -------
        %(outname)s : bool or array_like (if axis is specified)
            A single element array_like may be converted to bool."""

        def _make_logical_function(name, desc, f):

            @Substitution(outname=name, desc=desc)
            @Appender(_doc)
            def logical_func(self, *args, **kwargs):
                result = f(self.values)
                if isinstance(result, (np.ndarray, ABCSeries, Index)) \
                   and result.ndim == 0:
                    # return NumPy type
                    return result.dtype.type(result.item())
                else:  # pragma: no cover
                    return result
            logical_func.__name__ = name
            return logical_func

        cls.all = _make_logical_function(
            'all', 'Return whether all elements are True', np.all)
        cls.any = _make_logical_function(
            'any', 'Return whether any element is True', np.any)

    @classmethod
    def _add_logical_methods_disabled(cls):
        """ add in logical methods to disable """

        def _make_invalid_op(name):

            def invalid_op(self, other=None):
                raise TypeError("cannot perform {name} with this index type: {typ}".format(name=name,
                                                                                           typ=type(self)))
            invalid_op.__name__ = name
            return invalid_op

        cls.all = _make_invalid_op('all')
        cls.any = _make_invalid_op('any')


Index._add_numeric_methods_disabled()
Index._add_logical_methods()
Index._add_comparison_methods()

class CategoricalIndex(Index, PandasDelegate):
    """

    Immutable Index implementing an ordered, sliceable set. CategoricalIndex
    represents a sparsely populated Index with an underlying Categorical.

    .. versionadded:: 0.16.1

    Parameters
    ----------
    data : array-like or Categorical, (1-dimensional)
    categories : optional, array-like
        categories for the CategoricalIndex
    ordered : boolean,
        designating if the categories are ordered
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index

    """

    _typ = 'categoricalindex'
    _engine_type = _index.Int64Engine
    _attributes = ['name']

    def __new__(cls, data=None, categories=None, ordered=None, dtype=None, copy=False, name=None, fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name=name)

        if isinstance(data, ABCCategorical):
            data = cls._create_categorical(cls, data, categories, ordered)
        elif isinstance(data, CategoricalIndex):
            data = data._data
            data = cls._create_categorical(cls, data, categories, ordered)
        else:

            # don't allow scalars
            # if data is None, then categories must be provided
            if lib.isscalar(data):
                if data is not None or categories is None:
                    cls._scalar_data_error(data)
                data = []
            data = cls._create_categorical(cls, data, categories, ordered)

        if copy:
            data = data.copy()

        return cls._simple_new(data, name=name)

    def _create_from_codes(self, codes, categories=None, ordered=None, name=None):
        """
        *this is an internal non-public method*

        create the correct categorical from codes

        Parameters
        ----------
        codes : new codes
        categories : optional categories, defaults to existing
        ordered : optional ordered attribute, defaults to existing
        name : optional name attribute, defaults to existing

        Returns
        -------
        CategoricalIndex
        """

        from pandas.core.categorical import Categorical
        if categories is None:
            categories = self.categories
        if ordered is None:
            ordered = self.ordered
        if name is None:
            name = self.name
        cat = Categorical.from_codes(codes, categories=categories, ordered=self.ordered)
        return CategoricalIndex(cat, name=name)

    @staticmethod
    def _create_categorical(self, data, categories=None, ordered=None):
        """
        *this is an internal non-public method*

        create the correct categorical from data and the properties

        Parameters
        ----------
        data : data for new Categorical
        categories : optional categories, defaults to existing
        ordered : optional ordered attribute, defaults to existing

        Returns
        -------
        Categorical
        """

        if not isinstance(data, ABCCategorical):
            from pandas.core.categorical import Categorical
            data = Categorical(data, categories=categories, ordered=ordered)
        else:
            if categories is not None:
                data = data.set_categories(categories)
            if ordered is not None:
                data = data.set_ordered(ordered)
        return data

    @classmethod
    def _simple_new(cls, values, name=None, categories=None, ordered=None, **kwargs):
        result = object.__new__(cls)

        values = cls._create_categorical(cls, values, categories, ordered)
        result._data = values
        result.name = name
        for k, v in compat.iteritems(kwargs):
            setattr(result,k,v)

        result._reset_identity()
        return result

    def _is_dtype_compat(self, other):
        """
        *this is an internal non-public method*

        provide a comparison between the dtype of self and other (coercing if needed)

        Raises
        ------
        TypeError if the dtypes are not compatible
        """

        if is_categorical_dtype(other):
            if isinstance(other, CategoricalIndex):
                other = other._values
            if not other.is_dtype_equal(self):
                raise TypeError("categories must match existing categories when appending")
        else:
            values = other
            if not is_list_like(values):
                values = [ values ]
            other = CategoricalIndex(self._create_categorical(self, other, categories=self.categories, ordered=self.ordered))
            if not other.isin(values).all():
                raise TypeError("cannot append a non-category item to a CategoricalIndex")

        return other

    def equals(self, other):
        """
        Determines if two CategorialIndex objects contain the same elements.
        """
        if self.is_(other):
            return True

        try:
            other = self._is_dtype_compat(other)
            return array_equivalent(self._data, other)
        except (TypeError, ValueError):
            pass

        return False

    @property
    def _formatter_func(self):
        return self.categories._formatter_func

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        max_categories = (10 if get_option("display.max_categories") == 0
                    else get_option("display.max_categories"))
        attrs = [('categories', default_pprint(self.categories, max_seq_items=max_categories)),
                 ('ordered',self.ordered)]
        if self.name is not None:
            attrs.append(('name',default_pprint(self.name)))
        attrs.append(('dtype',"'%s'" % self.dtype))
        max_seq_items = get_option('display.max_seq_items') or len(self)
        if len(self) > max_seq_items:
            attrs.append(('length',len(self)))
        return attrs

    @property
    def inferred_type(self):
        return 'categorical'

    @property
    def values(self):
        """ return the underlying data, which is a Categorical """
        return self._data

    def get_values(self):
        """ return the underlying data as an ndarray """
        return self._data.get_values()

    @property
    def codes(self):
        return self._data.codes

    @property
    def categories(self):
        return self._data.categories

    @property
    def ordered(self):
        return self._data.ordered

    def __contains__(self, key):
        hash(key)
        return key in self.values

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return np.array(self._data, dtype=dtype)

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        return self._data.codes == -1

    @Appender(_index_shared_docs['fillna'])
    def fillna(self, value, downcast=None):
        self._assert_can_do_op(value)
        return CategoricalIndex(self._data.fillna(value), name=self.name)

    def argsort(self, *args, **kwargs):
        return self.values.argsort(*args, **kwargs)

    @cache_readonly
    def _engine(self):

        # we are going to look things up with the codes themselves
        return self._engine_type(lambda: self.codes.astype('i8'), len(self))

    @cache_readonly
    def is_unique(self):
        return not self.duplicated().any()

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(base._shared_docs['duplicated'] % _index_doc_kwargs)
    def duplicated(self, keep='first'):
        from pandas.hashtable import duplicated_int64
        return duplicated_int64(self.codes.astype('i8'), keep)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.astype('object')

    def get_loc(self, key, method=None):
        """
        Get integer location for requested label

        Parameters
        ----------
        key : label
        method : {None}
            * default: exact matches only.

        Returns
        -------
        loc : int if unique index, possibly slice or mask if not
        """
        codes = self.categories.get_loc(key)
        if (codes == -1):
            raise KeyError(key)
        indexer, _ = self._engine.get_indexer_non_unique(np.array([codes]))
        if (indexer==-1).any():
            raise KeyError(key)

        return indexer

    def _can_reindex(self, indexer):
        """ always allow reindexing """
        pass

    def reindex(self, target, method=None, level=None, limit=None,
                tolerance=None):
        """
        Create index with target's values (move/add/delete values as necessary)

        Returns
        -------
        new_index : pd.Index
            Resulting index
        indexer : np.ndarray or None
            Indices of output values in original index

        """

        if method is not None:
            raise NotImplementedError("argument method is not implemented for CategoricalIndex.reindex")
        if level is not None:
            raise NotImplementedError("argument level is not implemented for CategoricalIndex.reindex")
        if limit is not None:
            raise NotImplementedError("argument limit is not implemented for CategoricalIndex.reindex")

        target = _ensure_index(target)

        if not is_categorical_dtype(target) and not target.is_unique:
            raise ValueError("cannot reindex with a non-unique indexer")

        indexer, missing = self.get_indexer_non_unique(np.array(target))
        new_target = self.take(indexer)


        # filling in missing if needed
        if len(missing):
            cats = self.categories.get_indexer(target)
            if (cats==-1).any():

                # coerce to a regular index here!
                result = Index(np.array(self),name=self.name)
                new_target, indexer, _ = result._reindex_non_unique(np.array(target))

            else:

                codes = new_target.codes.copy()
                codes[indexer==-1] = cats[missing]
                new_target = self._create_from_codes(codes)

        # we always want to return an Index type here
        # to be consistent with .reindex for other index types (e.g. they don't coerce
        # based on the actual values, only on the dtype)
        # unless we had an inital Categorical to begin with
        # in which case we are going to conform to the passed Categorical
        new_target = np.asarray(new_target)
        if is_categorical_dtype(target):
            new_target = target._shallow_copy(new_target, name=self.name)
        else:
            new_target = Index(new_target, name=self.name)

        return new_target, indexer

    def _reindex_non_unique(self, target):
        """ reindex from a non-unique; which CategoricalIndex's are almost always """
        new_target, indexer = self.reindex(target)
        new_indexer = None

        check = indexer==-1
        if check.any():
            new_indexer = np.arange(len(self.take(indexer)))
            new_indexer[check] = -1

        return new_target, indexer, new_indexer

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        """
        Compute indexer and mask for new index given the current index. The
        indexer should be then used as an input to ndarray.take to align the
        current data to the new index. The mask determines whether labels are
        found or not in the current index

        Parameters
        ----------
        target : MultiIndex or Index (of tuples)
        method : {'pad', 'ffill', 'backfill', 'bfill'}
            pad / ffill: propagate LAST valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Notes
        -----
        This is a low-level method and probably should be used at your own risk

        Examples
        --------
        >>> indexer, mask = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)
        >>> new_values[-mask] = np.nan

        Returns
        -------
        (indexer, mask) : (ndarray, ndarray)
        """
        method = _clean_reindex_fill_method(method)
        target = _ensure_index(target)

        if isinstance(target, CategoricalIndex):
            target = target.categories

        if method == 'pad' or method == 'backfill':
            raise NotImplementedError("method='pad' and method='backfill' not implemented yet "
                                      'for CategoricalIndex')
        elif method == 'nearest':
            raise NotImplementedError("method='nearest' not implemented yet "
                                      'for CategoricalIndex')
        else:

            codes = self.categories.get_indexer(target)
            indexer, _ = self._engine.get_indexer_non_unique(codes)

        return com._ensure_platform_int(indexer)

    def get_indexer_non_unique(self, target):
        """ this is the same for a CategoricalIndex for get_indexer; the API returns the missing values as well """
        target = _ensure_index(target)

        if isinstance(target, CategoricalIndex):
            target = target.categories

        codes = self.categories.get_indexer(target)
        return self._engine.get_indexer_non_unique(codes)

    def _convert_list_indexer(self, keyarr, kind=None):
        """
        we are passed a list indexer.
        Return our indexer or raise if all of the values are not included in the categories
        """
        codes = self.categories.get_indexer(keyarr)
        if (codes==-1).any():
            raise KeyError("a list-indexer must only include values that are in the categories")

        return None

    def take(self, indexer, axis=0, allow_fill=True, fill_value=None):
        """
        For internal compatibility with numpy arrays.

        # filling must always be None/nan here
        # but is passed thru internally
        assert isnull(fill_value)

        See also
        --------
        numpy.ndarray.take
        """

        indexer = com._ensure_platform_int(indexer)
        taken = self.codes.take(indexer)
        return self._create_from_codes(taken)

    def delete(self, loc):
        """
        Make new Index with passed location(-s) deleted

        Returns
        -------
        new_index : Index
        """
        return self._create_from_codes(np.delete(self.codes, loc))

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location. Follows
        Python list.append semantics for negative values

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : Index

        Raises
        ------
        ValueError if the item is not in the categories

        """
        code = self.categories.get_indexer([item])
        if (code == -1):
            raise TypeError("cannot insert an item into a CategoricalIndex that is not already an existing category")

        codes = self.codes
        codes = np.concatenate(
            (codes[:loc], code, codes[loc:]))
        return self._create_from_codes(codes)

    def append(self, other):
        """
        Append a collection of CategoricalIndex options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index

        Raises
        ------
        ValueError if other is not in the categories
        """
        to_concat, name = self._ensure_compat_append(other)
        to_concat = [ self._is_dtype_compat(c) for c in to_concat ]
        codes = np.concatenate([ c.codes for c in to_concat ])
        return self._create_from_codes(codes, name=name)

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """

        def _make_compare(op):

            def _evaluate_compare(self, other):

                # if we have a Categorical type, then must have the same categories
                if isinstance(other, CategoricalIndex):
                    other = other._values
                elif isinstance(other, Index):
                    other = self._create_categorical(self, other._values, categories=self.categories, ordered=self.ordered)

                if isinstance(other, (ABCCategorical, np.ndarray, ABCSeries)):
                    if len(self.values) != len(other):
                        raise ValueError("Lengths must match to compare")

                if isinstance(other, ABCCategorical):
                    if not self.values.is_dtype_equal(other):
                        raise TypeError("categorical index comparisions must have the same categories and ordered attributes")

                return getattr(self.values, op)(other)

            return _evaluate_compare

        cls.__eq__ = _make_compare('__eq__')
        cls.__ne__ = _make_compare('__ne__')
        cls.__lt__ = _make_compare('__lt__')
        cls.__gt__ = _make_compare('__gt__')
        cls.__le__ = _make_compare('__le__')
        cls.__ge__ = _make_compare('__ge__')


    def _delegate_method(self, name, *args, **kwargs):
        """ method delegation to the ._values """
        method = getattr(self._values, name)
        if 'inplace' in kwargs:
            raise ValueError("cannot use inplace with CategoricalIndex")
        res = method(*args, **kwargs)
        if lib.isscalar(res):
            return res
        return CategoricalIndex(res, name=self.name)

    @classmethod
    def _add_accessors(cls):
        """ add in Categorical accessor methods """

        from pandas.core.categorical import Categorical
        CategoricalIndex._add_delegate_accessors(delegate=Categorical,
                                                 accessors=["rename_categories",
                                                            "reorder_categories",
                                                            "add_categories",
                                                            "remove_categories",
                                                            "remove_unused_categories",
                                                            "set_categories",
                                                            "as_ordered",
                                                            "as_unordered",
                                                            "min",
                                                            "max"],
                                                 typ='method',
                                                 overwrite=True)


CategoricalIndex._add_numericlike_set_methods_disabled()
CategoricalIndex._add_numeric_methods_disabled()
CategoricalIndex._add_logical_methods_disabled()
CategoricalIndex._add_comparison_methods()
CategoricalIndex._add_accessors()


class NumericIndex(Index):
    """
    Provide numeric type operations

    This is an abstract class

    """
    _is_numeric_dtype = True

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        This function should be overloaded in subclasses that allow non-trivial
        casting on label-slice bounds, e.g. datetime-like indices allowing
        strings containing formatted datetimes.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : string / None

        Returns
        -------
        label :  object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """

        # we are a numeric index, so we accept
        # integer/floats directly
        if not (is_integer(label) or is_float(label)):
            self._invalid_indexer('slice',label)

        return label

    def _convert_tolerance(self, tolerance):
        try:
            return float(tolerance)
        except ValueError:
            raise ValueError('tolerance argument for %s must be numeric: %r'
                             % (type(self).__name__, tolerance))


class Int64Index(NumericIndex):

    """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects. Int64Index is a special case
    of `Index` with purely integer labels. This is the default index type used
    by the DataFrame and Series ctors when no explicit index is provided by the
    user.

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: int64)
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index

    Notes
    -----
    An Index instance can **only** contain hashable objects
    """

    _typ = 'int64index'
    _groupby = _algos.groupby_int64
    _arrmap = _algos.arrmap_int64
    _left_indexer_unique = _algos.left_join_indexer_unique_int64
    _left_indexer = _algos.left_join_indexer_int64
    _inner_indexer = _algos.inner_join_indexer_int64
    _outer_indexer = _algos.outer_join_indexer_int64

    _can_hold_na = False

    _engine_type = _index.Int64Engine

    def __new__(cls, data=None, dtype=None, copy=False, name=None, fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name=name)

        # isscalar, generators handled in coerce_to_ndarray
        data = cls._coerce_to_ndarray(data)

        if issubclass(data.dtype.type, compat.string_types):
            cls._string_data_error(data)

        elif issubclass(data.dtype.type, np.integer):
            # don't force the upcast as we may be dealing
            # with a platform int
            if dtype is None or not issubclass(np.dtype(dtype).type,
                                               np.integer):
                dtype = np.int64

            subarr = np.array(data, dtype=dtype, copy=copy)
        else:
            subarr = np.array(data, dtype=np.int64, copy=copy)
            if len(data) > 0:
                if (subarr != data).any():
                    raise TypeError('Unsafe NumPy casting to integer, you must'
                                    ' explicitly cast')

        return cls._simple_new(subarr, name=name)

    @property
    def inferred_type(self):
        return 'integer'

    @property
    def asi8(self):
        # do not cache or you'll create a memory leak
        return self.values.view('i8')

    @property
    def is_all_dates(self):
        """
        Checks that all the labels are datetime objects
        """
        return False

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        # if not isinstance(other, Int64Index):
        #     return False

        try:
            return array_equivalent(_values_from_object(self), _values_from_object(other))
        except TypeError:
            # e.g. fails in numpy 1.6 with DatetimeIndex #1681
            return False

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        return Int64Index(joined, name=name)


Int64Index._add_numeric_methods()
Int64Index._add_logical_methods()


class Float64Index(NumericIndex):

    """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects. Float64Index is a special case
    of `Index` with purely floating point labels.

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: object)
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index

    Notes
    -----
    An Float64Index instance can **only** contain hashable objects
    """

    _typ = 'float64index'
    _engine_type = _index.Float64Engine
    _groupby = _algos.groupby_float64
    _arrmap = _algos.arrmap_float64
    _left_indexer_unique = _algos.left_join_indexer_unique_float64
    _left_indexer = _algos.left_join_indexer_float64
    _inner_indexer = _algos.inner_join_indexer_float64
    _outer_indexer = _algos.outer_join_indexer_float64

    def __new__(cls, data=None, dtype=None, copy=False, name=None, fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name)

        data = cls._coerce_to_ndarray(data)

        if issubclass(data.dtype.type, compat.string_types):
            cls._string_data_error(data)

        if dtype is None:
            dtype = np.float64

        try:
            subarr = np.array(data, dtype=dtype, copy=copy)
        except:
            raise TypeError('Unsafe NumPy casting, you must '
                            'explicitly cast')

        # coerce to float64 for storage
        if subarr.dtype != np.float64:
            subarr = subarr.astype(np.float64)

        return cls._simple_new(subarr, name)

    @property
    def inferred_type(self):
        return 'floating'

    def astype(self, dtype):
        if np.dtype(dtype) not in (np.object, np.float64):
            raise TypeError('Setting %s dtype to anything other than '
                            'float64 or object is not supported' %
                            self.__class__)
        return Index(self._values, name=self.name, dtype=dtype)

    def _convert_scalar_indexer(self, key, kind=None):
        """
        convert a scalar indexer

        Parameters
        ----------
        key : label of the slice bound
        kind : optional, type of the indexing operation (loc/ix/iloc/None)

        right now we are converting
        floats -> ints if the index supports it
        """

        if kind == 'iloc':
            if is_integer(key):
                return key
            return super(Float64Index, self)._convert_scalar_indexer(key, kind=kind)

        return key

    def _convert_slice_indexer(self, key, kind=None):
        """
        convert a slice indexer, by definition these are labels
        unless we are iloc

        Parameters
        ----------
        key : label of the slice bound
        kind : optional, type of the indexing operation (loc/ix/iloc/None)
        """

        # if we are not a slice, then we are done
        if not isinstance(key, slice):
            return key

        if kind == 'iloc':
            return super(Float64Index, self)._convert_slice_indexer(key,
                                                                    kind=kind)

        # translate to locations
        return self.slice_indexer(key.start, key.stop, key.step)

    def get_value(self, series, key):
        """ we always want to get an index value, never a value """
        if not np.isscalar(key):
            raise InvalidIndexError

        from pandas.core.indexing import maybe_droplevels
        from pandas.core.series import Series

        k = _values_from_object(key)
        loc = self.get_loc(k)
        new_values = _values_from_object(series)[loc]

        if np.isscalar(new_values) or new_values is None:
            return new_values

        new_index = self[loc]
        new_index = maybe_droplevels(new_index, k)
        return Series(new_values, index=new_index, name=series.name)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self is other:
            return True

        # need to compare nans locations and make sure that they are the same
        # since nans don't compare equal this is a bit tricky
        try:
            if not isinstance(other, Float64Index):
                other = self._constructor(other)
            if not is_dtype_equal(self.dtype,other.dtype) or self.shape != other.shape:
                return False
            left, right = self._values, other._values
            return ((left == right) | (self._isnan & other._isnan)).all()
        except TypeError:
            # e.g. fails in numpy 1.6 with DatetimeIndex #1681
            return False

    def __contains__(self, other):
        if super(Float64Index, self).__contains__(other):
            return True

        try:
            # if other is a sequence this throws a ValueError
            return np.isnan(other) and self.hasnans
        except ValueError:
            try:
                return len(other) <= 1 and _try_get_item(other) in self
            except TypeError:
                return False
        except:
            return False

    def get_loc(self, key, method=None, tolerance=None):
        try:
            if np.all(np.isnan(key)):
                nan_idxs = self._nan_idxs
                try:
                    return nan_idxs.item()
                except (ValueError, IndexError):
                    # should only need to catch ValueError here but on numpy
                    # 1.7 .item() can raise IndexError when NaNs are present
                    return nan_idxs
        except (TypeError, NotImplementedError):
            pass
        return super(Float64Index, self).get_loc(key, method=method,
                                                 tolerance=tolerance)

    @property
    def is_all_dates(self):
        """
        Checks that all the labels are datetime objects
        """
        return False

    @cache_readonly
    def is_unique(self):
        return super(Float64Index, self).is_unique and self._nan_idxs.size < 2

    @Appender(Index.isin.__doc__)
    def isin(self, values, level=None):
        value_set = set(values)
        if level is not None:
            self._validate_index_level(level)
        return lib.ismember_nans(np.array(self), value_set,
                                 isnull(list(value_set)).any())


Float64Index._add_numeric_methods()
Float64Index._add_logical_methods_disabled()


class MultiIndex(Index):

    """
    Implements multi-level, a.k.a. hierarchical, index object for pandas
    objects

    Parameters
    ----------
    levels : sequence of arrays
        The unique labels for each level
    labels : sequence of arrays
        Integers for each level designating which label at each location
    sortorder : optional int
        Level of sortedness (must be lexicographically sorted by that
        level)
    names : optional sequence of objects
        Names for each of the index levels. (name is accepted for compat)
    copy : boolean, default False
        Copy the meta-data
    verify_integrity : boolean, default True
        Check that the levels/labels are consistent and valid
    """

    # initialize to zero-length tuples to make everything work
    _typ = 'multiindex'
    _names = FrozenList()
    _levels = FrozenList()
    _labels = FrozenList()
    _comparables = ['names']
    rename = Index.set_names

    def __new__(cls, levels=None, labels=None, sortorder=None, names=None,
                copy=False, verify_integrity=True, _set_identity=True, name=None, **kwargs):

        # compat with Index
        if name is not None:
            names = name
        if levels is None or labels is None:
            raise TypeError("Must pass both levels and labels")
        if len(levels) != len(labels):
            raise ValueError('Length of levels and labels must be the same.')
        if len(levels) == 0:
            raise ValueError('Must pass non-zero number of levels/labels')
        if len(levels) == 1:
            if names:
                name = names[0]
            else:
                name = None
            return Index(levels[0], name=name, copy=True).take(labels[0])

        result = object.__new__(MultiIndex)

        # we've already validated levels and labels, so shortcut here
        result._set_levels(levels, copy=copy, validate=False)
        result._set_labels(labels, copy=copy, validate=False)

        if names is not None:
            # handles name validation
            result._set_names(names)

        if sortorder is not None:
            result.sortorder = int(sortorder)
        else:
            result.sortorder = sortorder

        if verify_integrity:
            result._verify_integrity()
        if _set_identity:
            result._reset_identity()

        return result

    def _verify_integrity(self):
        """Raises ValueError if length of levels and labels don't match or any
        label would exceed level bounds"""
        # NOTE: Currently does not check, among other things, that cached
        # nlevels matches nor that sortorder matches actually sortorder.
        labels, levels = self.labels, self.levels
        if len(levels) != len(labels):
            raise ValueError("Length of levels and labels must match. NOTE:"
                             " this index is in an inconsistent state.")
        label_length = len(self.labels[0])
        for i, (level, label) in enumerate(zip(levels, labels)):
            if len(label) != label_length:
                raise ValueError("Unequal label lengths: %s" % (
                                 [len(lab) for lab in labels]))
            if len(label) and label.max() >= len(level):
                raise ValueError("On level %d, label max (%d) >= length of"
                                 " level  (%d). NOTE: this index is in an"
                                 " inconsistent state" % (i, label.max(),
                                                          len(level)))

    def _get_levels(self):
        return self._levels

    def _set_levels(self, levels, level=None, copy=False, validate=True,
                    verify_integrity=False):
        # This is NOT part of the levels property because it should be
        # externally not allowed to set levels. User beware if you change
        # _levels directly
        if validate and len(levels) == 0:
            raise ValueError('Must set non-zero number of levels.')
        if validate and level is None and len(levels) != self.nlevels:
            raise ValueError('Length of levels must match number of levels.')
        if validate and level is not None and len(levels) != len(level):
            raise ValueError('Length of levels must match length of level.')

        if level is None:
            new_levels = FrozenList(_ensure_index(lev, copy=copy)._shallow_copy()
                                    for lev in levels)
        else:
            level = [self._get_level_number(l) for l in level]
            new_levels = list(self._levels)
            for l, v in zip(level, levels):
                new_levels[l] = _ensure_index(v, copy=copy)._shallow_copy()
            new_levels = FrozenList(new_levels)

        names = self.names
        self._levels = new_levels
        if any(names):
            self._set_names(names)

        self._tuples = None
        self._reset_cache()

        if verify_integrity:
            self._verify_integrity()

    def set_levels(self, levels, level=None, inplace=False, verify_integrity=True):
        """
        Set new levels on MultiIndex. Defaults to returning
        new index.

        Parameters
        ----------
        levels : sequence or list of sequence
            new level(s) to apply
        level : int or level name, or sequence of int / level names (default None)
            level(s) to set (None for all levels)
        inplace : bool
            if True, mutates in place
        verify_integrity : bool (default True)
            if True, checks that levels and labels are compatible

        Returns
        -------
        new index (of same type and class...etc)


        Examples
        --------
        >>> idx = MultiIndex.from_tuples([(1, u'one'), (1, u'two'),
                                          (2, u'one'), (2, u'two')],
                                          names=['foo', 'bar'])
        >>> idx.set_levels([['a','b'], [1,2]])
        MultiIndex(levels=[[u'a', u'b'], [1, 2]],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_levels(['a','b'], level=0)
        MultiIndex(levels=[[u'a', u'b'], [u'one', u'two']],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_levels(['a','b'], level='bar')
        MultiIndex(levels=[[1, 2], [u'a', u'b']],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_levels([['a','b'], [1,2]], level=[0,1])
        MultiIndex(levels=[[u'a', u'b'], [1, 2]],
                   labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                   names=[u'foo', u'bar'])
        """
        if level is not None and not is_list_like(level):
            if not is_list_like(levels):
                raise TypeError("Levels must be list-like")
            if is_list_like(levels[0]):
                raise TypeError("Levels must be list-like")
            level = [level]
            levels = [levels]
        elif level is None or is_list_like(level):
            if not is_list_like(levels) or not is_list_like(levels[0]):
                raise TypeError("Levels must be list of lists-like")

        if inplace:
            idx = self
        else:
            idx = self._shallow_copy()
        idx._reset_identity()
        idx._set_levels(levels, level=level, validate=True,
                        verify_integrity=verify_integrity)
        if not inplace:
            return idx

    # remove me in 0.14 and change to read only property
    __set_levels = deprecate("setting `levels` directly",
                             partial(set_levels, inplace=True,
                                     verify_integrity=True),
                             alt_name="set_levels")
    levels = property(fget=_get_levels, fset=__set_levels)

    def _get_labels(self):
        return self._labels

    def _set_labels(self, labels, level=None, copy=False, validate=True,
                    verify_integrity=False):

        if validate and level is None and len(labels) != self.nlevels:
            raise ValueError("Length of labels must match number of levels")
        if validate and level is not None and len(labels) != len(level):
            raise ValueError('Length of labels must match length of levels.')

        if level is None:
            new_labels = FrozenList(_ensure_frozen(lab, lev, copy=copy)._shallow_copy()
                                    for lev, lab in zip(self.levels, labels))
        else:
            level = [self._get_level_number(l) for l in level]
            new_labels = list(self._labels)
            for l, lev, lab in zip(level, self.levels, labels):
                new_labels[l] = _ensure_frozen(lab, lev, copy=copy)._shallow_copy()
            new_labels = FrozenList(new_labels)

        self._labels = new_labels
        self._tuples = None
        self._reset_cache()

        if verify_integrity:
            self._verify_integrity()

    def set_labels(self, labels, level=None, inplace=False, verify_integrity=True):
        """
        Set new labels on MultiIndex. Defaults to returning
        new index.

        Parameters
        ----------
        labels : sequence or list of sequence
            new labels to apply
        level : int or level name, or sequence of int / level names (default None)
            level(s) to set (None for all levels)
        inplace : bool
            if True, mutates in place
        verify_integrity : bool (default True)
            if True, checks that levels and labels are compatible

        Returns
        -------
        new index (of same type and class...etc)

        Examples
        --------
        >>> idx = MultiIndex.from_tuples([(1, u'one'), (1, u'two'),
                                          (2, u'one'), (2, u'two')],
                                          names=['foo', 'bar'])
        >>> idx.set_labels([[1,0,1,0], [0,0,1,1]])
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[1, 0, 1, 0], [0, 0, 1, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_labels([1,0,1,0], level=0)
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[1, 0, 1, 0], [0, 1, 0, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_labels([0,0,1,1], level='bar')
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[0, 0, 1, 1], [0, 0, 1, 1]],
                   names=[u'foo', u'bar'])
        >>> idx.set_labels([[1,0,1,0], [0,0,1,1]], level=[0,1])
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[1, 0, 1, 0], [0, 0, 1, 1]],
                   names=[u'foo', u'bar'])
        """
        if level is not None and not is_list_like(level):
            if not is_list_like(labels):
                raise TypeError("Labels must be list-like")
            if is_list_like(labels[0]):
                raise TypeError("Labels must be list-like")
            level = [level]
            labels = [labels]
        elif level is None or is_list_like(level):
            if not is_list_like(labels) or not is_list_like(labels[0]):
                raise TypeError("Labels must be list of lists-like")

        if inplace:
            idx = self
        else:
            idx = self._shallow_copy()
        idx._reset_identity()
        idx._set_labels(labels, level=level, verify_integrity=verify_integrity)
        if not inplace:
            return idx

    # remove me in 0.14 and change to readonly property
    __set_labels = deprecate("setting labels directly",
                             partial(set_labels, inplace=True,
                                     verify_integrity=True),
                             alt_name="set_labels")
    labels = property(fget=_get_labels, fset=__set_labels)

    def copy(self, names=None, dtype=None, levels=None, labels=None,
             deep=False, _set_identity=False):
        """
        Make a copy of this object. Names, dtype, levels and labels can be
        passed and will be set on new copy.

        Parameters
        ----------
        names : sequence, optional
        dtype : numpy dtype or pandas type, optional
        levels : sequence, optional
        labels : sequence, optional

        Returns
        -------
        copy : MultiIndex

        Notes
        -----
        In most cases, there should be no functional difference from using
        ``deep``, but if ``deep`` is passed it will attempt to deepcopy.
        This could be potentially expensive on large MultiIndex objects.
        """
        if deep:
            from copy import deepcopy
            levels = levels if levels is not None else deepcopy(self.levels)
            labels = labels if labels is not None else deepcopy(self.labels)
            names = names if names is not None else deepcopy(self.names)
        else:
            levels = self.levels
            labels = self.labels
            names = self.names
        return MultiIndex(levels=levels,
                          labels=labels,
                          names=names,
                          sortorder=self.sortorder,
                          verify_integrity=False,
                          _set_identity=_set_identity)

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return self.values

    def view(self, cls=None):
        """ this is defined as a copy with the same identity """
        result = self.copy()
        result._id = self._id
        return result

    def _shallow_copy(self, values=None, infer=False, **kwargs):
        if values is not None:
            if 'name' in kwargs:
                kwargs['names'] = kwargs.pop('name',None)
            return MultiIndex.from_tuples(values, **kwargs)
        return self.view()

    @cache_readonly
    def dtype(self):
        return np.dtype('O')

    @cache_readonly
    def nbytes(self):
        """ return the number of bytes in the underlying data """
        level_nbytes = sum(( i.nbytes for i in self.levels ))
        label_nbytes = sum(( i.nbytes for i in self.labels ))
        names_nbytes = sum(( getsizeof(i) for i in self.names ))
        return level_nbytes + label_nbytes + names_nbytes

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        attrs = [('levels', default_pprint(self._levels, max_seq_items=False)),
                 ('labels', default_pprint(self._labels, max_seq_items=False))]
        if not all(name is None for name in self.names):
            attrs.append(('names', default_pprint(self.names)))
        if self.sortorder is not None:
            attrs.append(('sortorder', default_pprint(self.sortorder)))
        return attrs

    def _format_space(self):
        return "\n%s" % (' ' * (len(self.__class__.__name__) + 1))

    def _format_data(self):
        # we are formatting thru the attributes
        return None

    def __len__(self):
        return len(self.labels[0])

    def _get_names(self):
        return FrozenList(level.name for level in self.levels)

    def _set_names(self, names, level=None, validate=True):
        """
        sets names on levels. WARNING: mutates!

        Note that you generally want to set this *after* changing levels, so
        that it only acts on copies
        """

        names = list(names)

        if validate and level is not None and len(names) != len(level):
            raise ValueError('Length of names must match length of level.')
        if validate and level is None and len(names) != self.nlevels:
            raise ValueError(
                'Length of names must match number of levels in MultiIndex.')

        if level is None:
            level = range(self.nlevels)
        else:
            level = [self._get_level_number(l) for l in level]

        # set the name
        for l, name in zip(level, names):
            self.levels[l].rename(name, inplace=True)

    names = property(
        fset=_set_names, fget=_get_names, doc="Names of levels in MultiIndex")

    def _reference_duplicate_name(self, name):
        """
        Returns True if the name refered to in self.names is duplicated.
        """
        # count the times name equals an element in self.names.
        return sum(name == n for n in self.names) > 1

    def _format_native_types(self, **kwargs):
        # we go through the levels and format them
        levels = [level._format_native_types(**kwargs)
                  for level in self.levels]
        mi = MultiIndex(levels=levels, labels=self.labels, names=self.names,
                        sortorder=self.sortorder, verify_integrity=False)
        return mi.values

    @property
    def _constructor(self):
        return MultiIndex.from_tuples

    @cache_readonly
    def inferred_type(self):
        return 'mixed'

    @staticmethod
    def _from_elements(values, labels=None, levels=None, names=None,
                       sortorder=None):
        return MultiIndex(levels, labels, names, sortorder=sortorder)

    def _get_level_number(self, level):
        try:
            count = self.names.count(level)
            if count > 1:
                raise ValueError('The name %s occurs multiple times, use a '
                                 'level number' % level)
            level = self.names.index(level)
        except ValueError:
            if not isinstance(level, int):
                raise KeyError('Level %s not found' % str(level))
            elif level < 0:
                level += self.nlevels
                if level < 0:
                    orig_level = level - self.nlevels
                    raise IndexError(
                        'Too many levels: Index has only %d levels, '
                        '%d is not a valid level number' % (self.nlevels, orig_level)
                    )
            # Note: levels are zero-based
            elif level >= self.nlevels:
                raise IndexError('Too many levels: Index has only %d levels, '
                                 'not %d' % (self.nlevels, level + 1))
        return level

    _tuples = None

    @property
    def values(self):
        if self._tuples is not None:
            return self._tuples

        values = []
        for lev, lab in zip(self.levels, self.labels):
            # Need to box timestamps, etc.
            box = hasattr(lev, '_box_values')
            # Try to minimize boxing.
            if box and len(lev) > len(lab):
                taken = lev._box_values(com.take_1d(lev._values, lab))
            elif box:
                taken = com.take_1d(lev._box_values(lev._values), lab,
                                    fill_value=_get_na_value(lev.dtype.type))
            else:
                taken = com.take_1d(np.asarray(lev._values), lab)
            values.append(taken)

        self._tuples = lib.fast_zip(values)
        return self._tuples

    # fml
    @property
    def _is_v1(self):
        return False

    @property
    def _is_v2(self):
        return False

    @property
    def _has_complex_internals(self):
        # to disable groupby tricks
        return True

    @cache_readonly
    def is_unique(self):
        return not self.duplicated().any()

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(base._shared_docs['duplicated'] % _index_doc_kwargs)
    def duplicated(self, keep='first'):
        from pandas.core.groupby import get_group_index
        from pandas.hashtable import duplicated_int64

        shape = map(len, self.levels)
        ids = get_group_index(self.labels, shape, sort=False, xnull=False)

        return duplicated_int64(ids, keep)

    @Appender(_index_shared_docs['fillna'])
    def fillna(self, value=None, downcast=None):
        # isnull is not implemented for MultiIndex
        raise NotImplementedError('isnull is not defined for MultiIndex')

    def get_value(self, series, key):
        # somewhat broken encapsulation
        from pandas.core.indexing import maybe_droplevels
        from pandas.core.series import Series

        # Label-based
        s = _values_from_object(series)
        k = _values_from_object(key)

        def _try_mi(k):
            # TODO: what if a level contains tuples??
            loc = self.get_loc(k)
            new_values = series._values[loc]
            new_index = self[loc]
            new_index = maybe_droplevels(new_index, k)
            return Series(new_values, index=new_index, name=series.name)

        try:
            return self._engine.get_value(s, k)
        except KeyError as e1:
            try:
                return _try_mi(key)
            except KeyError:
                pass

            try:
                return _index.get_value_at(s, k)
            except IndexError:
                raise
            except TypeError:
                # generator/iterator-like
                if is_iterator(key):
                    raise InvalidIndexError(key)
                else:
                    raise e1
            except Exception:  # pragma: no cover
                raise e1
        except TypeError:

            # a Timestamp will raise a TypeError in a multi-index
            # rather than a KeyError, try it here
            # note that a string that 'looks' like a Timestamp will raise
            # a KeyError! (GH5725)
            if isinstance(key, (datetime.datetime, np.datetime64)) or (
                    compat.PY3 and isinstance(key, compat.string_types)):
                try:
                    return _try_mi(key)
                except (KeyError):
                    raise
                except:
                    pass

                try:
                    return _try_mi(Timestamp(key))
                except:
                    pass

            raise InvalidIndexError(key)

    def get_level_values(self, level):
        """
        Return vector of label values for requested level, equal to the length
        of the index

        Parameters
        ----------
        level : int or level name

        Returns
        -------
        values : ndarray
        """
        num = self._get_level_number(level)
        unique = self.levels[num]  # .values
        labels = self.labels[num]
        filled = com.take_1d(unique._values, labels, fill_value=unique._na_value)
        values = unique._simple_new(filled, self.names[num],
                                    freq=getattr(unique, 'freq', None),
                                    tz=getattr(unique, 'tz', None))
        return values

    def format(self, space=2, sparsify=None, adjoin=True, names=False,
               na_rep=None, formatter=None):
        if len(self) == 0:
            return []

        stringified_levels = []
        for lev, lab in zip(self.levels, self.labels):
            na = na_rep if na_rep is not None else _get_na_rep(lev.dtype.type)

            if len(lev) > 0:

                formatted = lev.take(lab).format(formatter=formatter)

                # we have some NA
                mask = lab == -1
                if mask.any():
                    formatted = np.array(formatted, dtype=object)
                    formatted[mask] = na
                    formatted = formatted.tolist()

            else:
                # weird all NA case
                formatted = [com.pprint_thing(na if isnull(x) else x,
                                              escape_chars=('\t', '\r', '\n'))
                             for x in com.take_1d(lev._values, lab)]
            stringified_levels.append(formatted)

        result_levels = []
        for lev, name in zip(stringified_levels, self.names):
            level = []

            if names:
                level.append(com.pprint_thing(name,
                                              escape_chars=('\t', '\r', '\n'))
                             if name is not None else '')

            level.extend(np.array(lev, dtype=object))
            result_levels.append(level)

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        if sparsify:
            sentinel = ''
            # GH3547
            # use value of sparsify as sentinel,  unless it's an obvious
            # "Truthey" value
            if sparsify not in [True, 1]:
                sentinel = sparsify
            # little bit of a kludge job for #1217
            result_levels = _sparsify(result_levels,
                                      start=int(names),
                                      sentinel=sentinel)


        if adjoin:
            from pandas.core.format import  _get_adjustment
            adj = _get_adjustment()
            return adj.adjoin(space, *result_levels).split('\n')
        else:
            return result_levels

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.set_levels([ i._to_safe_for_reshape() for i in self.levels ])

    def to_hierarchical(self, n_repeat, n_shuffle=1):
        """
        Return a MultiIndex reshaped to conform to the
        shapes given by n_repeat and n_shuffle.

        Useful to replicate and rearrange a MultiIndex for combination
        with another Index with n_repeat items.

        Parameters
        ----------
        n_repeat : int
            Number of times to repeat the labels on self
        n_shuffle : int
            Controls the reordering of the labels. If the result is going
            to be an inner level in a MultiIndex, n_shuffle will need to be
            greater than one. The size of each label must divisible by
            n_shuffle.

        Returns
        -------
        MultiIndex

        Examples
        --------
        >>> idx = MultiIndex.from_tuples([(1, u'one'), (1, u'two'),
                                          (2, u'one'), (2, u'two')])
        >>> idx.to_hierarchical(3)
        MultiIndex(levels=[[1, 2], [u'one', u'two']],
                   labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                           [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
        """
        levels = self.levels
        labels = [np.repeat(x, n_repeat) for x in self.labels]
        # Assumes that each label is divisible by n_shuffle
        labels = [x.reshape(n_shuffle, -1).ravel(1) for x in labels]
        names = self.names
        return MultiIndex(levels=levels, labels=labels, names=names)

    @property
    def is_all_dates(self):
        return False

    def is_lexsorted(self):
        """
        Return True if the labels are lexicographically sorted
        """
        return self.lexsort_depth == self.nlevels

    def is_lexsorted_for_tuple(self, tup):
        """
        Return True if we are correctly lexsorted given the passed tuple
        """
        return len(tup) <= self.lexsort_depth

    @cache_readonly
    def lexsort_depth(self):
        if self.sortorder is not None:
            if self.sortorder == 0:
                return self.nlevels
            else:
                return 0

        int64_labels = [com._ensure_int64(lab) for lab in self.labels]
        for k in range(self.nlevels, 0, -1):
            if lib.is_lexsorted(int64_labels[:k]):
                return k

        return 0

    @classmethod
    def from_arrays(cls, arrays, sortorder=None, names=None):
        """
        Convert arrays to MultiIndex

        Parameters
        ----------
        arrays : list / sequence of array-likes
            Each array-like gives one level's value for each data point.
            len(arrays) is the number of levels.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level)

        Returns
        -------
        index : MultiIndex

        Examples
        --------
        >>> arrays = [[1, 1, 2, 2], ['red', 'blue', 'red', 'blue']]
        >>> MultiIndex.from_arrays(arrays, names=('number', 'color'))

        See Also
        --------
        MultiIndex.from_tuples : Convert list of tuples to MultiIndex
        MultiIndex.from_product : Make a MultiIndex from cartesian product
                                  of iterables
        """
        from pandas.core.categorical import Categorical

        if len(arrays) == 1:
            name = None if names is None else names[0]
            return Index(arrays[0], name=name)

        cats = [Categorical.from_array(arr, ordered=True) for arr in arrays]
        levels = [c.categories for c in cats]
        labels = [c.codes for c in cats]
        if names is None:
            names = [getattr(arr, "name", None) for arr in arrays]

        return MultiIndex(levels=levels, labels=labels,
                          sortorder=sortorder, names=names,
                          verify_integrity=False)

    @classmethod
    def from_tuples(cls, tuples, sortorder=None, names=None):
        """
        Convert list of tuples to MultiIndex

        Parameters
        ----------
        tuples : list / sequence of tuple-likes
            Each tuple is the index of one row/column.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level)

        Returns
        -------
        index : MultiIndex

        Examples
        --------
        >>> tuples = [(1, u'red'), (1, u'blue'),
                      (2, u'red'), (2, u'blue')]
        >>> MultiIndex.from_tuples(tuples, names=('number', 'color'))

        See Also
        --------
        MultiIndex.from_arrays : Convert list of arrays to MultiIndex
        MultiIndex.from_product : Make a MultiIndex from cartesian product
                                  of iterables
        """
        if len(tuples) == 0:
            # I think this is right? Not quite sure...
            raise TypeError('Cannot infer number of levels from empty list')

        if isinstance(tuples, (np.ndarray, Index)):
            if isinstance(tuples, Index):
                tuples = tuples._values

            arrays = list(lib.tuples_to_object_array(tuples).T)
        elif isinstance(tuples, list):
            arrays = list(lib.to_object_array_tuples(tuples).T)
        else:
            arrays = lzip(*tuples)

        return MultiIndex.from_arrays(arrays, sortorder=sortorder,
                                      names=names)

    @classmethod
    def from_product(cls, iterables, sortorder=None, names=None):
        """
        Make a MultiIndex from the cartesian product of multiple iterables

        Parameters
        ----------
        iterables : list / sequence of iterables
            Each iterable has unique labels for each level of the index.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level).
        names : list / sequence of strings or None
            Names for the levels in the index.

        Returns
        -------
        index : MultiIndex

        Examples
        --------
        >>> numbers = [0, 1, 2]
        >>> colors = [u'green', u'purple']
        >>> MultiIndex.from_product([numbers, colors],
                                     names=['number', 'color'])
        MultiIndex(levels=[[0, 1, 2], [u'green', u'purple']],
                   labels=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
                   names=[u'number', u'color'])

        See Also
        --------
        MultiIndex.from_arrays : Convert list of arrays to MultiIndex
        MultiIndex.from_tuples : Convert list of tuples to MultiIndex
        """
        from pandas.core.categorical import Categorical
        from pandas.tools.util import cartesian_product

        categoricals = [Categorical.from_array(it, ordered=True) for it in iterables]
        labels = cartesian_product([c.codes for c in categoricals])

        return MultiIndex(levels=[c.categories for c in categoricals],
                          labels=labels, sortorder=sortorder, names=names)

    @property
    def nlevels(self):
        return len(self.levels)

    @property
    def levshape(self):
        return tuple(len(x) for x in self.levels)

    def __contains__(self, key):
        hash(key)
        # work around some kind of odd cython bug
        try:
            self.get_loc(key)
            return True
        except LookupError:
            return False

    def __reduce__(self):
        """Necessary for making this object picklable"""
        d = dict(levels = [lev for lev in self.levels],
                 labels = [label for label in self.labels],
                 sortorder = self.sortorder,
                 names = list(self.names))
        return _new_Index, (self.__class__, d), None

    def __setstate__(self, state):
        """Necessary for making this object picklable"""

        if isinstance(state, dict):
            levels = state.get('levels')
            labels = state.get('labels')
            sortorder = state.get('sortorder')
            names = state.get('names')

        elif isinstance(state, tuple):

            nd_state, own_state = state
            levels, labels, sortorder, names = own_state

        self._set_levels([Index(x) for x in levels], validate=False)
        self._set_labels(labels)
        self._set_names(names)
        self.sortorder = sortorder
        self._verify_integrity()
        self._reset_identity()

    def __getitem__(self, key):
        if np.isscalar(key):
            retval = []
            for lev, lab in zip(self.levels, self.labels):
                if lab[key] == -1:
                    retval.append(np.nan)
                else:
                    retval.append(lev[lab[key]])

            return tuple(retval)
        else:
            if is_bool_indexer(key):
                key = np.asarray(key)
                sortorder = self.sortorder
            else:
                # cannot be sure whether the result will be sorted
                sortorder = None

            new_labels = [lab[key] for lab in self.labels]

            return MultiIndex(levels=self.levels,
                              labels=new_labels,
                              names=self.names,
                              sortorder=sortorder,
                              verify_integrity=False)

    def take(self, indexer, axis=None):
        indexer = com._ensure_platform_int(indexer)
        new_labels = [lab.take(indexer) for lab in self.labels]
        return MultiIndex(levels=self.levels, labels=new_labels,
                          names=self.names, verify_integrity=False)

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        if not isinstance(other, (list, tuple)):
            other = [other]

        if all((isinstance(o, MultiIndex) and o.nlevels >= self.nlevels) for o in other):
            arrays = []
            for i in range(self.nlevels):
                label = self.get_level_values(i)
                appended = [o.get_level_values(i) for o in other]
                arrays.append(label.append(appended))
            return MultiIndex.from_arrays(arrays, names=self.names)

        to_concat = (self.values,) + tuple(k._values for k in other)
        new_tuples = np.concatenate(to_concat)

        # if all(isinstance(x, MultiIndex) for x in other):
        try:
            return MultiIndex.from_tuples(new_tuples, names=self.names)
        except:
            return Index(new_tuples)

    def argsort(self, *args, **kwargs):
        return self.values.argsort(*args, **kwargs)

    def repeat(self, n):
        return MultiIndex(levels=self.levels,
                          labels=[label.view(np.ndarray).repeat(n) for label in self.labels],
                          names=self.names,
                          sortorder=self.sortorder,
                          verify_integrity=False)

    def drop(self, labels, level=None, errors='raise'):
        """
        Make new MultiIndex with passed list of labels deleted

        Parameters
        ----------
        labels : array-like
            Must be a list of tuples
        level : int or level name, default None

        Returns
        -------
        dropped : MultiIndex
        """
        if level is not None:
            return self._drop_from_level(labels, level)

        try:
            if not isinstance(labels, (np.ndarray, Index)):
                labels = com._index_labels_to_array(labels)
            indexer = self.get_indexer(labels)
            mask = indexer == -1
            if mask.any():
                if errors != 'ignore':
                    raise ValueError('labels %s not contained in axis'
                                     % labels[mask])
                indexer = indexer[~mask]
        except Exception:
            pass

        inds = []
        for label in labels:
            try:
                loc = self.get_loc(label)
                if isinstance(loc, int):
                    inds.append(loc)
                else:
                    inds.extend(lrange(loc.start, loc.stop))
            except KeyError:
                if errors != 'ignore':
                    raise

        return self.delete(inds)

    def _drop_from_level(self, labels, level):
        labels = com._index_labels_to_array(labels)
        i = self._get_level_number(level)
        index = self.levels[i]
        values = index.get_indexer(labels)

        mask = ~lib.ismember(self.labels[i], set(values))

        return self[mask]

    def droplevel(self, level=0):
        """
        Return Index with requested level removed. If MultiIndex has only 2
        levels, the result will be of Index type not MultiIndex.

        Parameters
        ----------
        level : int/level name or list thereof

        Notes
        -----
        Does not check if result index is unique or not

        Returns
        -------
        index : Index or MultiIndex
        """
        levels = level
        if not isinstance(levels, (tuple, list)):
            levels = [level]

        new_levels = list(self.levels)
        new_labels = list(self.labels)
        new_names = list(self.names)

        levnums = sorted(self._get_level_number(lev) for lev in levels)[::-1]

        for i in levnums:
            new_levels.pop(i)
            new_labels.pop(i)
            new_names.pop(i)

        if len(new_levels) == 1:

            # set nan if needed
            mask = new_labels[0] == -1
            result = new_levels[0].take(new_labels[0])
            if mask.any():
                result = result.putmask(mask, np.nan)

            result.name = new_names[0]
            return result
        else:
            return MultiIndex(levels=new_levels, labels=new_labels,
                              names=new_names, verify_integrity=False)

    def swaplevel(self, i, j):
        """
        Swap level i with level j. Do not change the ordering of anything

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : MultiIndex
        """
        new_levels = list(self.levels)
        new_labels = list(self.labels)
        new_names = list(self.names)

        i = self._get_level_number(i)
        j = self._get_level_number(j)

        new_levels[i], new_levels[j] = new_levels[j], new_levels[i]
        new_labels[i], new_labels[j] = new_labels[j], new_labels[i]
        new_names[i], new_names[j] = new_names[j], new_names[i]

        return MultiIndex(levels=new_levels, labels=new_labels,
                          names=new_names, verify_integrity=False)

    def reorder_levels(self, order):
        """
        Rearrange levels using input order. May not drop or duplicate levels

        Parameters
        ----------
        """
        order = [self._get_level_number(i) for i in order]
        if len(order) != self.nlevels:
            raise AssertionError(('Length of order must be same as '
                                  'number of levels (%d), got %d')
                                 % (self.nlevels, len(order)))
        new_levels = [self.levels[i] for i in order]
        new_labels = [self.labels[i] for i in order]
        new_names = [self.names[i] for i in order]

        return MultiIndex(levels=new_levels, labels=new_labels,
                          names=new_names, verify_integrity=False)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def sortlevel(self, level=0, ascending=True, sort_remaining=True):
        """
        Sort MultiIndex at the requested level. The result will respect the
        original ordering of the associated factor at that level.

        Parameters
        ----------
        level : list-like, int or str, default 0
            If a string is given, must be a name of the level
            If list-like must be names or ints of levels.
        ascending : boolean, default True
            False to sort in descending order
            Can also be a list to specify a directed ordering
        sort_remaining : sort by the remaining levels after level.

        Returns
        -------
        sorted_index : MultiIndex
        """
        from pandas.core.groupby import _indexer_from_factorized

        if isinstance(level, (compat.string_types, int)):
            level = [level]
        level = [self._get_level_number(lev) for lev in level]
        sortorder = None

        # we have a directed ordering via ascending
        if isinstance(ascending, list):
            if not len(level) == len(ascending):
                raise ValueError("level must have same length as ascending")

            from pandas.core.groupby import _lexsort_indexer
            indexer = _lexsort_indexer(self.labels, orders=ascending)

        # level ordering
        else:

            labels = list(self.labels)
            shape = list(self.levshape)

            # partition labels and shape
            primary = tuple(labels.pop(lev - i) for i, lev in enumerate(level))
            primshp = tuple(shape.pop(lev - i) for i, lev in enumerate(level))

            if sort_remaining:
                primary += primary + tuple(labels)
                primshp += primshp + tuple(shape)
            else:
                sortorder = level[0]

            indexer = _indexer_from_factorized(primary,
                                               primshp,
                                               compress=False)

            if not ascending:
                indexer = indexer[::-1]

        indexer = com._ensure_platform_int(indexer)
        new_labels = [lab.take(indexer) for lab in self.labels]

        new_index = MultiIndex(labels=new_labels, levels=self.levels,
                               names=self.names, sortorder=sortorder,
                               verify_integrity=False)

        return new_index, indexer

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        """
        Compute indexer and mask for new index given the current index. The
        indexer should be then used as an input to ndarray.take to align the
        current data to the new index. The mask determines whether labels are
        found or not in the current index

        Parameters
        ----------
        target : MultiIndex or Index (of tuples)
        method : {'pad', 'ffill', 'backfill', 'bfill'}
            pad / ffill: propagate LAST valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Notes
        -----
        This is a low-level method and probably should be used at your own risk

        Examples
        --------
        >>> indexer, mask = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)
        >>> new_values[-mask] = np.nan

        Returns
        -------
        (indexer, mask) : (ndarray, ndarray)
        """
        method = _clean_reindex_fill_method(method)

        target = _ensure_index(target)

        target_index = target
        if isinstance(target, MultiIndex):
            target_index = target._tuple_index

        if not is_object_dtype(target_index.dtype):
            return np.ones(len(target_index)) * -1

        if not self.is_unique:
            raise Exception('Reindexing only valid with uniquely valued Index '
                            'objects')

        self_index = self._tuple_index

        if method == 'pad' or method == 'backfill':
            if tolerance is not None:
                raise NotImplementedError("tolerance not implemented yet "
                                          'for MultiIndex')
            indexer = self_index._get_fill_indexer(target, method, limit)
        elif method == 'nearest':
            raise NotImplementedError("method='nearest' not implemented yet "
                                      'for MultiIndex; see GitHub issue 9365')
        else:
            indexer = self_index._engine.get_indexer(target._values)

        return com._ensure_platform_int(indexer)

    def reindex(self, target, method=None, level=None, limit=None,
                tolerance=None):
        """
        Create index with target's values (move/add/delete values as necessary)

        Returns
        -------
        new_index : pd.MultiIndex
            Resulting index
        indexer : np.ndarray or None
            Indices of output values in original index

        """
        # GH6552: preserve names when reindexing to non-named target
        # (i.e. neither Index nor Series).
        preserve_names = not hasattr(target, 'names')

        if level is not None:
            if method is not None:
                raise TypeError('Fill method not supported if level passed')

            # GH7774: preserve dtype/tz if target is empty and not an Index.
            target = _ensure_has_len(target)  # target may be an iterator
            if len(target) == 0 and not isinstance(target, Index):
                idx = self.levels[level]
                attrs = idx._get_attributes_dict()
                attrs.pop('freq', None)  # don't preserve freq
                target = type(idx)._simple_new(np.empty(0, dtype=idx.dtype),
                                               **attrs)
            else:
                target = _ensure_index(target)
            target, indexer, _ = self._join_level(target, level, how='right',
                                                  return_indexers=True,
                                                  keep_order=False)
        else:
            if self.equals(target):
                indexer = None
            else:
                if self.is_unique:
                    indexer = self.get_indexer(target, method=method,
                                               limit=limit,
                                               tolerance=tolerance)
                else:
                    raise Exception(
                        "cannot handle a non-unique multi-index!")

        if not isinstance(target, MultiIndex):
            if indexer is None:
                target = self
            elif (indexer >= 0).all():
                target = self.take(indexer)
            else:
                # hopefully?
                target = MultiIndex.from_tuples(target)

        if (preserve_names and target.nlevels == self.nlevels and
            target.names != self.names):
            target = target.copy(deep=False)
            target.names = self.names

        return target, indexer

    @cache_readonly
    def _tuple_index(self):
        """
        Convert MultiIndex to an Index of tuples

        Returns
        -------
        index : Index
        """
        return Index(self._values)

    def get_slice_bound(self, label, side, kind):
        if not isinstance(label, tuple):
            label = label,
        return self._partial_tup_index(label, side=side)

    def slice_locs(self, start=None, end=None, step=None, kind=None):
        """
        For an ordered MultiIndex, compute the slice locations for input
        labels. They can be tuples representing partial levels, e.g. for a
        MultiIndex with 3 levels, you can pass a single value (corresponding to
        the first level), or a 1-, 2-, or 3-tuple.

        Parameters
        ----------
        start : label or tuple, default None
            If None, defaults to the beginning
        end : label or tuple
            If None, defaults to the end
        step : int or None
            Slice step
        kind : string, optional, defaults None

        Returns
        -------
        (start, end) : (int, int)

        Notes
        -----
        This function assumes that the data is sorted by the first level
        """
        # This function adds nothing to its parent implementation (the magic
        # happens in get_slice_bound method), but it adds meaningful doc.
        return super(MultiIndex, self).slice_locs(start, end, step, kind=kind)

    def _partial_tup_index(self, tup, side='left'):
        if len(tup) > self.lexsort_depth:
            raise KeyError('Key length (%d) was greater than MultiIndex'
                           ' lexsort depth (%d)' %
                           (len(tup), self.lexsort_depth))

        n = len(tup)
        start, end = 0, len(self)
        zipped = zip(tup, self.levels, self.labels)
        for k, (lab, lev, labs) in enumerate(zipped):
            section = labs[start:end]

            if lab not in lev:
                if not lev.is_type_compatible(lib.infer_dtype([lab])):
                    raise TypeError('Level type mismatch: %s' % lab)

                # short circuit
                loc = lev.searchsorted(lab, side=side)
                if side == 'right' and loc >= 0:
                    loc -= 1
                return start + section.searchsorted(loc, side=side)

            idx = lev.get_loc(lab)
            if k < n - 1:
                end = start + section.searchsorted(idx, side='right')
                start = start + section.searchsorted(idx, side='left')
            else:
                return start + section.searchsorted(idx, side=side)

    def get_loc(self, key, method=None):
        """
        Get integer location, slice or boolean mask for requested label or tuple
        If the key is past the lexsort depth, the return may be a boolean mask
        array, otherwise it is always a slice or int.

        Parameters
        ----------
        key : label or tuple
        method : None

        Returns
        -------
        loc : int, slice object or boolean mask
        """
        if method is not None:
            raise NotImplementedError('only the default get_loc method is '
                                      'currently supported for MultiIndex')

        def _maybe_to_slice(loc):
            '''convert integer indexer to boolean mask or slice if possible'''
            if not isinstance(loc, np.ndarray) or loc.dtype != 'int64':
                return loc

            loc = lib.maybe_indices_to_slice(loc, len(self))
            if isinstance(loc, slice):
                return loc

            mask = np.empty(len(self), dtype='bool')
            mask.fill(False)
            mask[loc] = True
            return mask

        if not isinstance(key, tuple):
            loc = self._get_level_indexer(key, level=0)
            return _maybe_to_slice(loc)

        keylen = len(key)
        if self.nlevels < keylen:
            raise KeyError('Key length ({0}) exceeds index depth ({1})'
                    ''.format(keylen, self.nlevels))

        if keylen == self.nlevels and self.is_unique:
            def _maybe_str_to_time_stamp(key, lev):
                if lev.is_all_dates and not isinstance(key, Timestamp):
                    try:
                        return Timestamp(key, tz=getattr(lev, 'tz', None))
                    except Exception:
                        pass
                return key
            key = _values_from_object(key)
            key = tuple(map(_maybe_str_to_time_stamp, key, self.levels))
            return self._engine.get_loc(key)

        # -- partial selection or non-unique index
        # break the key into 2 parts based on the lexsort_depth of the index;
        # the first part returns a continuous slice of the index; the 2nd part
        # needs linear search within the slice
        i = self.lexsort_depth
        lead_key, follow_key = key[:i], key[i:]
        start, stop = self.slice_locs(lead_key, lead_key) \
                if lead_key else (0, len(self))

        if start == stop:
            raise KeyError(key)

        if not follow_key:
            return slice(start, stop)

        warnings.warn('indexing past lexsort depth may impact performance.',
                      PerformanceWarning, stacklevel=10)

        loc = np.arange(start, stop, dtype='int64')

        for i, k in enumerate(follow_key, len(lead_key)):
            mask = self.labels[i][loc] == self.levels[i].get_loc(k)
            if not mask.all():
                loc = loc[mask]
            if not len(loc):
                raise KeyError(key)

        return _maybe_to_slice(loc) \
                if len(loc) != stop - start \
                else slice(start, stop)

    def get_loc_level(self, key, level=0, drop_level=True):
        """
        Get integer location slice for requested label or tuple

        Parameters
        ----------
        key : label or tuple
        level : int/level name or list thereof

        Returns
        -------
        loc : int or slice object
        """
        def maybe_droplevels(indexer, levels, drop_level):
            if not drop_level:
                return self[indexer]
            # kludgearound
            orig_index = new_index = self[indexer]
            levels = [self._get_level_number(i) for i in levels]
            for i in sorted(levels, reverse=True):
                try:
                    new_index = new_index.droplevel(i)
                except:

                    # no dropping here
                    return orig_index
            return new_index

        if isinstance(level, (tuple, list)):
            if len(key) != len(level):
                raise AssertionError('Key for location must have same '
                                     'length as number of levels')
            result = None
            for lev, k in zip(level, key):
                loc, new_index = self.get_loc_level(k, level=lev)
                if isinstance(loc, slice):
                    mask = np.zeros(len(self), dtype=bool)
                    mask[loc] = True
                    loc = mask

                result = loc if result is None else result & loc

            return result, maybe_droplevels(result, level, drop_level)

        level = self._get_level_number(level)

        # kludge for #1796
        if isinstance(key, list):
            key = tuple(key)

        if isinstance(key, tuple) and level == 0:

            try:
                if key in self.levels[0]:
                    indexer = self._get_level_indexer(key, level=level)
                    new_index = maybe_droplevels(indexer, [0], drop_level)
                    return indexer, new_index
            except TypeError:
                pass

            if not any(isinstance(k, slice) for k in key):

                # partial selection
                # optionally get indexer to avoid re-calculation
                def partial_selection(key, indexer=None):
                    if indexer is None:
                        indexer = self.get_loc(key)
                    ilevels = [i for i in range(len(key))
                               if key[i] != slice(None, None)]
                    return indexer, maybe_droplevels(indexer, ilevels,
                                                     drop_level)

                if len(key) == self.nlevels:

                    if self.is_unique:

                        # here we have a completely specified key, but are
                        # using some partial string matching here
                        # GH4758
                        can_index_exactly = any([
                            (l.is_all_dates and
                             not isinstance(k, compat.string_types))
                            for k, l in zip(key, self.levels)
                        ])
                        if any([
                            l.is_all_dates for k, l in zip(key, self.levels)
                        ]) and not can_index_exactly:
                            indexer = self.get_loc(key)

                            # we have a multiple selection here
                            if not isinstance(indexer, slice) \
                                    or indexer.stop - indexer.start != 1:
                                return partial_selection(key, indexer)

                            key = tuple(self[indexer].tolist()[0])

                        return (self._engine.get_loc(_values_from_object(key)),
                                None)
                    else:
                        return partial_selection(key)
                else:
                    return partial_selection(key)
            else:
                indexer = None
                for i, k in enumerate(key):
                    if not isinstance(k, slice):
                        k = self._get_level_indexer(k, level=i)
                        if isinstance(k, slice):
                            # everything
                            if k.start == 0 and k.stop == len(self):
                                k = slice(None, None)
                        else:
                            k_index = k

                    if isinstance(k, slice):
                        if k == slice(None, None):
                            continue
                        else:
                            raise TypeError(key)

                    if indexer is None:
                        indexer = k_index
                    else:  # pragma: no cover
                        indexer &= k_index
                if indexer is None:
                    indexer = slice(None, None)
                ilevels = [i for i in range(len(key))
                           if key[i] != slice(None, None)]
                return indexer, maybe_droplevels(indexer, ilevels,
                                                 drop_level)
        else:
            indexer = self._get_level_indexer(key, level=level)
            return indexer, maybe_droplevels(indexer, [level], drop_level)

    def _get_level_indexer(self, key, level=0, indexer=None):
        # return an indexer, boolean array or a slice showing where the key is
        # in the totality of values
        # if the indexer is provided, then use this

        level_index = self.levels[level]
        labels = self.labels[level]

        def convert_indexer(start, stop, step, indexer=indexer, labels=labels):
            # given the inputs and the labels/indexer, compute an indexer set
            # if we have a provided indexer, then this need not consider
            # the entire labels set

            r = np.arange(start,stop,step)
            if indexer is not None and len(indexer) != len(labels):

                # we have an indexer which maps the locations in the labels that we
                # have already selected (and is not an indexer for the entire set)
                # otherwise this is wasteful
                # so we only need to examine locations that are in this set
                # the only magic here is that the result are the mappings to the
                # set that we have selected
                from pandas import Series
                mapper = Series(indexer)
                indexer = labels.take(com._ensure_platform_int(indexer))
                result = Series(Index(indexer).isin(r).nonzero()[0])
                m = result.map(mapper)._values

            else:
                m = np.zeros(len(labels),dtype=bool)
                m[np.in1d(labels,r,assume_unique=True)] = True

            return m

        if isinstance(key, slice):
            # handle a slice, returnig a slice if we can
            # otherwise a boolean indexer

            try:
                if key.start is not None:
                    start = level_index.get_loc(key.start)
                else:
                    start = 0
                if key.stop is not None:
                    stop  = level_index.get_loc(key.stop)
                else:
                    stop = len(level_index)-1
                step = key.step
            except (KeyError):

                # we have a partial slice (like looking up a partial date string)
                start = stop = level_index.slice_indexer(key.start, key.stop, key.step)
                step = start.step

            if isinstance(start,slice) or isinstance(stop,slice):
                # we have a slice for start and/or stop
                # a partial date slicer on a DatetimeIndex generates a slice
                # note that the stop ALREADY includes the stopped point (if
                # it was a string sliced)
                return convert_indexer(start.start,stop.stop,step)

            elif level > 0 or self.lexsort_depth == 0 or step is not None:
                # need to have like semantics here to right
                # searching as when we are using a slice
                # so include the stop+1 (so we include stop)
                return convert_indexer(start,stop+1,step)
            else:
                # sorted, so can return slice object -> view
                i = labels.searchsorted(start, side='left')
                j = labels.searchsorted(stop, side='right')
                return slice(i, j, step)

        else:

            loc = level_index.get_loc(key)
            if level > 0 or self.lexsort_depth == 0:
                return np.array(labels == loc,dtype=bool)
            else:
                # sorted, so can return slice object -> view
                i = labels.searchsorted(loc, side='left')
                j = labels.searchsorted(loc, side='right')
                return slice(i, j)

    def get_locs(self, tup):
        """
        Given a tuple of slices/lists/labels/boolean indexer to a level-wise spec
        produce an indexer to extract those locations

        Parameters
        ----------
        key : tuple of (slices/list/labels)

        Returns
        -------
        locs : integer list of locations or boolean indexer suitable
               for passing to iloc
        """

        # must be lexsorted to at least as many levels
        if not self.is_lexsorted_for_tuple(tup):
            raise KeyError('MultiIndex Slicing requires the index to be fully lexsorted'
                           ' tuple len ({0}), lexsort depth ({1})'.format(len(tup), self.lexsort_depth))

        # indexer
        # this is the list of all values that we want to select
        n = len(self)
        indexer = None

        def _convert_to_indexer(r):
            # return an indexer
            if isinstance(r, slice):
                m = np.zeros(n,dtype=bool)
                m[r] = True
                r = m.nonzero()[0]
            elif is_bool_indexer(r):
                if len(r) != n:
                    raise ValueError("cannot index with a boolean indexer that is"
                                     " not the same length as the index")
                r = r.nonzero()[0]
            return Int64Index(r)

        def _update_indexer(idxr, indexer=indexer):
            if indexer is None:
                indexer = Index(np.arange(n))
            if idxr is None:
                return indexer
            return indexer & idxr

        for i,k in enumerate(tup):

            if is_bool_indexer(k):
                # a boolean indexer, must be the same length!
                k = np.asarray(k)
                indexer = _update_indexer(_convert_to_indexer(k), indexer=indexer)

            elif is_list_like(k):
                # a collection of labels to include from this level (these are or'd)
                indexers = None
                for x in k:
                    try:
                        idxrs = _convert_to_indexer(self._get_level_indexer(x, level=i, indexer=indexer))
                        indexers = idxrs if indexers is None else indexers | idxrs
                    except (KeyError):

                        # ignore not founds
                        continue

                if indexers is not None:
                    indexer = _update_indexer(indexers, indexer=indexer)
                else:

                    # no matches we are done
                    return Int64Index([])._values

            elif is_null_slice(k):
                # empty slice
                indexer = _update_indexer(None, indexer=indexer)

            elif isinstance(k,slice):

                # a slice, include BOTH of the labels
                indexer = _update_indexer(_convert_to_indexer(self._get_level_indexer(k,level=i,indexer=indexer)), indexer=indexer)
            else:
                # a single label
                indexer = _update_indexer(_convert_to_indexer(self.get_loc_level(k,level=i,drop_level=False)[0]), indexer=indexer)

        # empty indexer
        if indexer is None:
            return Int64Index([])._values
        return indexer._values

    def truncate(self, before=None, after=None):
        """
        Slice index between two labels / tuples, return new MultiIndex

        Parameters
        ----------
        before : label or tuple, can be partial. Default None
            None defaults to start
        after : label or tuple, can be partial. Default None
            None defaults to end

        Returns
        -------
        truncated : MultiIndex
        """
        if after and before and after < before:
            raise ValueError('after < before')

        i, j = self.levels[0].slice_locs(before, after)
        left, right = self.slice_locs(before, after)

        new_levels = list(self.levels)
        new_levels[0] = new_levels[0][i:j]

        new_labels = [lab[left:right] for lab in self.labels]
        new_labels[0] = new_labels[0] - i

        return MultiIndex(levels=new_levels, labels=new_labels,
                          verify_integrity=False)

    def equals(self, other):
        """
        Determines if two MultiIndex objects have the same labeling information
        (the levels themselves do not necessarily have to be the same)

        See also
        --------
        equal_levels
        """
        if self.is_(other):
            return True

        if not isinstance(other, MultiIndex):
            return array_equivalent(self._values,
                                    _values_from_object(_ensure_index(other)))

        if self.nlevels != other.nlevels:
            return False

        if len(self) != len(other):
            return False

        for i in range(self.nlevels):
            svalues = com.take_nd(np.asarray(self.levels[i]._values), self.labels[i],
                                  allow_fill=False)
            ovalues = com.take_nd(np.asarray(other.levels[i]._values), other.labels[i],
                                  allow_fill=False)
            if not array_equivalent(svalues, ovalues):
                return False

        return True

    def equal_levels(self, other):
        """
        Return True if the levels of both MultiIndex objects are the same

        """
        if self.nlevels != other.nlevels:
            return False

        for i in range(self.nlevels):
            if not self.levels[i].equals(other.levels[i]):
                return False
        return True

    def union(self, other):
        """
        Form the union of two MultiIndex objects, sorting if possible

        Parameters
        ----------
        other : MultiIndex or array / Index of tuples

        Returns
        -------
        Index

        >>> index.union(index2)
        """
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if len(other) == 0 or self.equals(other):
            return self

        uniq_tuples = lib.fast_unique_multiple([self._values, other._values])
        return MultiIndex.from_arrays(lzip(*uniq_tuples), sortorder=0,
                                      names=result_names)

    def intersection(self, other):
        """
        Form the intersection of two MultiIndex objects, sorting if possible

        Parameters
        ----------
        other : MultiIndex or array / Index of tuples

        Returns
        -------
        Index
        """
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if self.equals(other):
            return self

        self_tuples = self._values
        other_tuples = other._values
        uniq_tuples = sorted(set(self_tuples) & set(other_tuples))
        if len(uniq_tuples) == 0:
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names, verify_integrity=False)
        else:
            return MultiIndex.from_arrays(lzip(*uniq_tuples), sortorder=0,
                                          names=result_names)

    def difference(self, other):
        """
        Compute sorted set difference of two MultiIndex objects

        Returns
        -------
        diff : MultiIndex
        """
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if len(other) == 0:
                return self

        if self.equals(other):
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names, verify_integrity=False)

        difference = sorted(set(self._values) - set(other._values))

        if len(difference) == 0:
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names, verify_integrity=False)
        else:
            return MultiIndex.from_tuples(difference, sortorder=0,
                                          names=result_names)

    def astype(self, dtype):
        if not is_object_dtype(np.dtype(dtype)):
            raise TypeError('Setting %s dtype to anything other than object '
                            'is not supported' % self.__class__)
        return self._shallow_copy()

    def _convert_can_do_setop(self, other):
        result_names = self.names

        if not hasattr(other, 'names'):
            if len(other) == 0:
                other = MultiIndex(levels=[[]] * self.nlevels,
                                   labels=[[]] * self.nlevels,
                                   verify_integrity=False)
            else:
                msg = 'other must be a MultiIndex or a list of tuples'
                try:
                    other = MultiIndex.from_tuples(other)
                except:
                    raise TypeError(msg)
        else:
            result_names = self.names if self.names == other.names else None
        return other, result_names

    def insert(self, loc, item):
        """
        Make new MultiIndex inserting new item at location

        Parameters
        ----------
        loc : int
        item : tuple
            Must be same length as number of levels in the MultiIndex

        Returns
        -------
        new_index : Index
        """
        # Pad the key with empty strings if lower levels of the key
        # aren't specified:
        if not isinstance(item, tuple):
            item = (item,) + ('',) * (self.nlevels - 1)
        elif len(item) != self.nlevels:
            raise ValueError(
                'Item must have length equal to number of levels.')

        new_levels = []
        new_labels = []
        for k, level, labels in zip(item, self.levels, self.labels):
            if k not in level:
                # have to insert into level
                # must insert at end otherwise you have to recompute all the
                # other labels
                lev_loc = len(level)
                level = level.insert(lev_loc, k)
            else:
                lev_loc = level.get_loc(k)

            new_levels.append(level)
            new_labels.append(np.insert(_ensure_int64(labels), loc, lev_loc))

        return MultiIndex(levels=new_levels, labels=new_labels,
                          names=self.names, verify_integrity=False)

    def delete(self, loc):
        """
        Make new index with passed location deleted

        Returns
        -------
        new_index : MultiIndex
        """
        new_labels = [np.delete(lab, loc) for lab in self.labels]
        return MultiIndex(levels=self.levels, labels=new_labels,
                          names=self.names, verify_integrity=False)

    get_major_bounds = slice_locs

    __bounds = None

    @property
    def _bounds(self):
        """
        Return or compute and return slice points for level 0, assuming
        sortedness
        """
        if self.__bounds is None:
            inds = np.arange(len(self.levels[0]))
            self.__bounds = self.labels[0].searchsorted(inds)

        return self.__bounds

    def _wrap_joined_index(self, joined, other):
        names = self.names if self.names == other.names else None
        return MultiIndex.from_tuples(joined, names=names)

    @Appender(Index.isin.__doc__)
    def isin(self, values, level=None):
        if level is None:
            return lib.ismember(np.array(self), set(values))
        else:
            num = self._get_level_number(level)
            levs = self.levels[num]
            labs = self.labels[num]

            sought_labels = levs.isin(values).nonzero()[0]
            if levs.size == 0:
                return np.zeros(len(labs), dtype=np.bool_)
            else:
                return np.lib.arraysetops.in1d(labs, sought_labels)


MultiIndex._add_numeric_methods_disabled()
MultiIndex._add_logical_methods_disabled()


# For utility purposes

def _sparsify(label_list, start=0, sentinel=''):
    pivoted = lzip(*label_list)
    k = len(label_list)

    result = pivoted[:start + 1]
    prev = pivoted[start]

    for cur in pivoted[start + 1:]:
        sparse_cur = []

        for i, (p, t) in enumerate(zip(prev, cur)):
            if i == k - 1:
                sparse_cur.append(t)
                result.append(sparse_cur)
                break

            if p == t:
                sparse_cur.append(sentinel)
            else:
                sparse_cur.extend(cur[i:])
                result.append(sparse_cur)
                break

        prev = cur

    return lzip(*result)


def _ensure_index(index_like, copy=False):
    if isinstance(index_like, Index):
        if copy:
            index_like = index_like.copy()
        return index_like
    if hasattr(index_like, 'name'):
        return Index(index_like, name=index_like.name, copy=copy)

    # must check for exactly list here because of strict type
    # check in clean_index_list
    if isinstance(index_like, list):
        if type(index_like) != list:
            index_like = list(index_like)
        # 2200 ?
        converted, all_arrays = lib.clean_index_list(index_like)

        if len(converted) > 0 and all_arrays:
            return MultiIndex.from_arrays(converted)
        else:
            index_like = converted
    else:
       # clean_index_list does the equivalent of copying
       # so only need to do this if not list instance
        if copy:
            from copy import copy
            index_like = copy(index_like)

    return Index(index_like)


def _ensure_frozen(array_like, categories, copy=False):
    array_like = com._coerce_indexer_dtype(array_like, categories)
    array_like = array_like.view(FrozenNDArray)
    if copy:
        array_like = array_like.copy()
    return array_like


def _validate_join_method(method):
    if method not in ['left', 'right', 'inner', 'outer']:
        raise ValueError('do not recognize join method %s' % method)


# TODO: handle index names!
def _get_combined_index(indexes, intersect=False):
    indexes = _get_distinct_indexes(indexes)
    if len(indexes) == 0:
        return Index([])
    if len(indexes) == 1:
        return indexes[0]
    if intersect:
        index = indexes[0]
        for other in indexes[1:]:
            index = index.intersection(other)
        return index
    union = _union_indexes(indexes)
    return _ensure_index(union)


def _get_distinct_indexes(indexes):
    return list(dict((id(x), x) for x in indexes).values())


def _union_indexes(indexes):
    if len(indexes) == 0:
        raise AssertionError('Must have at least 1 Index to union')
    if len(indexes) == 1:
        result = indexes[0]
        if isinstance(result, list):
            result = Index(sorted(result))
        return result

    indexes, kind = _sanitize_and_check(indexes)
    def _unique_indices(inds):
        def conv(i):
            if isinstance(i, Index):
                i = i.tolist()
            return i
        return Index(lib.fast_unique_multiple_list([ conv(i) for i in inds ]))

    if kind == 'special':
        result = indexes[0]

        if hasattr(result, 'union_many'):
            return result.union_many(indexes[1:])
        else:
            for other in indexes[1:]:
                result = result.union(other)
            return result
    elif kind == 'array':
        index = indexes[0]
        for other in indexes[1:]:
            if not index.equals(other):
                return _unique_indices(indexes)

        return index
    else:
        return _unique_indices(indexes)


def _trim_front(strings):
    """
    Trims zeros and decimal points
    """
    trimmed = strings
    while len(strings) > 0 and all([x[0] == ' ' for x in trimmed]):
        trimmed = [x[1:] for x in trimmed]
    return trimmed


def _sanitize_and_check(indexes):
    kinds = list(set([type(index) for index in indexes]))

    if list in kinds:
        if len(kinds) > 1:
            indexes = [Index(com._try_sort(x))
                       if not isinstance(x, Index) else x
                       for x in indexes]
            kinds.remove(list)
        else:
            return indexes, 'list'

    if len(kinds) > 1 or Index not in kinds:
        return indexes, 'special'
    else:
        return indexes, 'array'


def _get_consensus_names(indexes):

    # find the non-none names, need to tupleify to make
    # the set hashable, then reverse on return
    consensus_names = set([
        tuple(i.names) for i in indexes if all(n is not None for n in i.names)
    ])
    if len(consensus_names) == 1:
        return list(list(consensus_names)[0])
    return [None] * indexes[0].nlevels


def _maybe_box(idx):
    from pandas.tseries.api import DatetimeIndex, PeriodIndex, TimedeltaIndex
    klasses = DatetimeIndex, PeriodIndex, TimedeltaIndex

    if isinstance(idx, klasses):
        return idx.asobject
    return idx


def _all_indexes_same(indexes):
    first = indexes[0]
    for index in indexes[1:]:
        if not first.equals(index):
            return False
    return True


def _get_na_rep(dtype):
    return {np.datetime64: 'NaT', np.timedelta64: 'NaT'}.get(dtype, 'NaN')


def _get_na_value(dtype):
    return {np.datetime64: tslib.NaT, np.timedelta64: tslib.NaT}.get(dtype,
                                                                     np.nan)


def _ensure_has_len(seq):
    """If seq is an iterator, put its values into a list."""
    try:
        len(seq)
    except TypeError:
        return list(seq)
    else:
        return seq
