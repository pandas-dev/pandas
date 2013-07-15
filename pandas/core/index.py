# pylint: disable=E1101,E1103,W0232

from itertools import izip

import numpy as np

import pandas.tslib as tslib
import pandas.lib as lib
import pandas.algos as _algos
import pandas.index as _index
from pandas.lib import Timestamp
from pandas.core.base import PandasObject

from pandas.util.decorators import cache_readonly
from pandas.core.common import isnull
import pandas.core.common as com
from pandas.core.config import get_option


__all__ = ['Index']


def _indexOp(opname):
    """
    Wrapper function for index comparison operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        func = getattr(self.view(np.ndarray), opname)
        result = func(other)
        try:
            return result.view(np.ndarray)
        except:  # pragma: no cover
            return result
    return wrapper


class InvalidIndexError(Exception):
    pass

_o_dtype = np.dtype(object)


def _shouldbe_timestamp(obj):
    return (tslib.is_datetime_array(obj)
            or tslib.is_datetime64_array(obj)
            or tslib.is_timestamp_array(obj))


class Index(PandasObject, np.ndarray):
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

    Note
    ----
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

    name = None
    asi8 = None

    _engine_type = _index.ObjectEngine

    def __new__(cls, data, dtype=None, copy=False, name=None, **kwargs):
        from pandas.tseries.period import PeriodIndex
        if isinstance(data, np.ndarray):
            if issubclass(data.dtype.type, np.datetime64):
                from pandas.tseries.index import DatetimeIndex
                result = DatetimeIndex(data, copy=copy, name=name, **kwargs)
                if dtype is not None and _o_dtype == dtype:
                    return Index(result.to_pydatetime(), dtype=_o_dtype)
                else:
                    return result
            elif issubclass(data.dtype.type, np.timedelta64):
                return Int64Index(data, copy=copy, name=name)

            if dtype is not None:
                try:
                    data = np.array(data, dtype=dtype, copy=copy)
                except TypeError:
                    pass
            elif isinstance(data, PeriodIndex):
                return PeriodIndex(data, copy=copy, name=name, **kwargs)

            if issubclass(data.dtype.type, np.integer):
                return Int64Index(data, copy=copy, dtype=dtype, name=name)

            subarr = com._asarray_tuplesafe(data, dtype=object)
        elif np.isscalar(data):
            raise ValueError('Index(...) must be called with a collection '
                             'of some kind, %s was passed' % repr(data))
        else:
            # other iterable of some kind
            subarr = com._asarray_tuplesafe(data, dtype=object)

        if dtype is None:
            inferred = lib.infer_dtype(subarr)
            if inferred == 'integer':
                return Int64Index(subarr.astype('i8'), name=name)
            elif inferred != 'string':
                if (inferred.startswith('datetime') or
                        tslib.is_timestamp_array(subarr)):
                    from pandas.tseries.index import DatetimeIndex
                    return DatetimeIndex(subarr, copy=copy, name=name, **kwargs)

                elif inferred == 'period':
                    return PeriodIndex(subarr, name=name, **kwargs)

        subarr = subarr.view(cls)
        subarr.name = name
        return subarr

    def __array_finalize__(self, obj):
        if not isinstance(obj, type(self)):
            # Only relevant if array being created from an Index instance
            return

        self.name = getattr(obj, 'name', None)

    def _shallow_copy(self):
        return self.view()

    def __unicode__(self):
        """
        Return a string representation for a particular Index

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both py2/py3.
        """
        prepr = com.pprint_thing(self, escape_chars=('\t', '\r', '\n'),quote_strings=True)
        return '%s(%s, dtype=%s)' % (type(self).__name__, prepr, self.dtype)

    def to_series(self):
        """
        return a series with both index and values equal to the index keys
        useful with map for returning an indexer based on an index
        """
        import pandas as pd
        return pd.Series(self.values,index=self,name=self.name)

    def astype(self, dtype):
        return Index(self.values.astype(dtype), name=self.name,
                     dtype=dtype)

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
        return True

    def tolist(self):
        """
        Overridden version of ndarray.tolist
        """
        return list(self.values)

    @cache_readonly
    def dtype(self):
        return self.values.dtype

    @property
    def nlevels(self):
        return 1

    # for compat with multindex code

    def _get_names(self):
        return [self.name]

    def _set_names(self, values):
        if len(values) != 1:
            raise AssertionError('Length of new names must be 1, got %d'
                                 % len(values))
        self.name = values[0]

    names = property(fset=_set_names, fget=_get_names)

    @property
    def _has_complex_internals(self):
        # to disable groupby tricks in MultiIndex
        return False

    def summary(self, name=None):
        if len(self) > 0:
            head = self[0]
            if hasattr(head,'format'):
                head = head.format()
            tail = self[-1]
            if hasattr(tail,'format'):
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

    @property
    def values(self):
        return np.asarray(self)

    @property
    def is_monotonic(self):
        return self._engine.is_monotonic

    def is_lexsorted_for_tuple(self, tup):
        return True

    @cache_readonly(allow_setting=True)
    def is_unique(self):
        return self._engine.is_unique

    def is_numeric(self):
        return self.inferred_type in ['integer', 'floating']

    def holds_integer(self):
        return self.inferred_type in ['integer', 'mixed-integer']

    def get_duplicates(self):
        from collections import defaultdict
        counter = defaultdict(lambda: 0)
        for k in self.values:
            counter[k] += 1
        return sorted(k for k, v in counter.iteritems() if v > 1)

    _get_duplicates = get_duplicates

    def _cleanup(self):
        self._engine.clear_mapping()

    @cache_readonly
    def _engine(self):
        # property, for now, slow to look up
        return self._engine_type(lambda: self.values, len(self))

    def _get_level_number(self, level):
        if not isinstance(level, int):
            if level != self.name:
                raise AssertionError('Level %s must be same as name (%s)'
                                     % (level, self.name))
            level = 0
        return level

    @cache_readonly
    def inferred_type(self):
        return lib.infer_dtype(self)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type

    @cache_readonly
    def is_all_dates(self):
        return self.inferred_type == 'datetime'

    def __iter__(self):
        return iter(self.values)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = self.name,
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if len(state) == 2:
            nd_state, own_state = state
            np.ndarray.__setstate__(self, nd_state)
            self.name = own_state[0]
        else:  # pragma: no cover
            np.ndarray.__setstate__(self, state)

    def __deepcopy__(self, memo={}):
        """
        Index is not mutable, so disabling deepcopy
        """
        return self

    def __contains__(self, key):
        hash(key)
        # work around some kind of odd cython bug
        try:
            return key in self._engine
        except TypeError:
            return False

    def __hash__(self):
        return hash(self.view(np.ndarray))

    def __setitem__(self, key, value):
        raise Exception(str(self.__class__) + ' object is immutable')

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            return arr_idx[key]
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)

            result = arr_idx[key]
            if result.ndim > 1:
                return result

            return Index(result, name=self.name)

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
        name = self.name
        to_concat = [self]

        if isinstance(other, (list, tuple)):
            to_concat = to_concat + list(other)
        else:
            to_concat.append(other)

        for obj in to_concat:
            if isinstance(obj, Index) and obj.name != name:
                name = None
                break

        to_concat = self._ensure_compat_concat(to_concat)
        to_concat = [x.values if isinstance(x, Index) else x
                     for x in to_concat]

        return Index(np.concatenate(to_concat), name=name)

    @staticmethod
    def _ensure_compat_concat(indexes):
        from pandas.tseries.api import DatetimeIndex, PeriodIndex
        klasses = DatetimeIndex, PeriodIndex

        is_ts = [isinstance(idx, klasses) for idx in indexes]

        if any(is_ts) and not all(is_ts):
            return [_maybe_box(idx) for idx in indexes]

        return indexes

    def take(self, indexer, axis=0):
        """
        Analogous to ndarray.take
        """
        indexer = com._ensure_platform_int(indexer)
        taken = self.view(np.ndarray).take(indexer)
        return self._constructor(taken, name=self.name)

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

        if values.dtype == np.object_:
            values = lib.maybe_convert_objects(values, safe=1)

        if values.dtype == np.object_:
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

    def _format_native_types(self, na_rep='', **kwargs):
        """ actually format my specific types """
        mask = isnull(self)
        values = np.array(self,dtype=object,copy=True)
        values[mask] = na_rep
        return values.tolist()

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self is other:
            return True

        if not isinstance(other, Index):
            return False

        if type(other) != Index:
            return other.equals(self)

        return np.array_equal(self, other)

    def asof(self, label):
        """
        For a sorted index, return the most recent label up to and including
        the passed label. Return NaN if not found
        """
        if isinstance(label, (Index, np.ndarray)):
            raise TypeError('%s' % type(label))

        if label not in self:
            loc = self.searchsorted(label, side='left')
            if loc > 0:
                return self[loc - 1]
            else:
                return np.nan

        if not isinstance(label, Timestamp):
            label = Timestamp(label)
        return label

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

    def order(self, return_indexer=False, ascending=True):
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

    def sort(self, *args, **kwargs):
        raise Exception('Cannot sort an Index object')

    def shift(self, periods=1, freq=None):
        """
        Shift Index containing datetime objects by input number of periods and
        DateOffset

        Returns
        -------
        shifted : Index
        """
        if periods == 0:
            # OK because immutable
            return self

        offset = periods * freq
        return Index([idx + offset for idx in self], name=self.name)

    def argsort(self, *args, **kwargs):
        """
        See docstring for ndarray.argsort
        """
        return self.view(np.ndarray).argsort(*args, **kwargs)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        else:
            return Index(self.view(np.ndarray) + other)

    __eq__ = _indexOp('__eq__')
    __ne__ = _indexOp('__ne__')
    __lt__ = _indexOp('__lt__')
    __gt__ = _indexOp('__gt__')
    __le__ = _indexOp('__le__')
    __ge__ = _indexOp('__ge__')

    def __sub__(self, other):
        return self.diff(other)

    def __and__(self, other):
        return self.intersection(other)

    def __or__(self, other):
        return self.union(other)

    def union(self, other):
        """
        Form the union of two Index objects and sorts if possible

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        union : Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if len(other) == 0 or self.equals(other):
            return self

        if len(self) == 0:
            return _ensure_index(other)

        self._assert_can_do_setop(other)

        if self.dtype != other.dtype:
            this = self.astype('O')
            other = other.astype('O')
            return this.union(other)

        if self.is_monotonic and other.is_monotonic:
            try:
                result = self._outer_indexer(self, other.values)[0]
            except TypeError:
                # incomparable objects
                result = list(self.values)

                # worth making this faster? a very unusual case
                value_set = set(self.values)
                result.extend([x for x in other.values if x not in value_set])
        else:
            indexer = self.get_indexer(other)
            indexer = (indexer == -1).nonzero()[0]

            if len(indexer) > 0:
                other_diff = com.take_nd(other.values, indexer,
                                         allow_fill=False)
                result = com._concat_compat((self.values, other_diff))
                try:
                    result.sort()
                except Exception:
                    pass
            else:
                # contained in
                try:
                    result = np.sort(self.values)
                except TypeError:  # pragma: no cover
                    result = self.values

        # for subclasses
        return self._wrap_union_result(other, result)

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        return type(self)(data=result, name=name)

    def intersection(self, other):
        """
        Form the intersection of two Index objects. Sortedness of the result is
        not guaranteed

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        intersection : Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        self._assert_can_do_setop(other)

        other = _ensure_index(other)

        if self.equals(other):
            return self

        if self.dtype != other.dtype:
            this = self.astype('O')
            other = other.astype('O')
            return this.intersection(other)

        if self.is_monotonic and other.is_monotonic:
            try:
                result = self._inner_indexer(self, other.values)[0]
                return self._wrap_union_result(other, result)
            except TypeError:
                pass

        indexer = self.get_indexer(other.values)
        indexer = indexer.take((indexer != -1).nonzero()[0])
        return self.take(indexer)

    def diff(self, other):
        """
        Compute sorted set difference of two Index objects

        Notes
        -----
        One can do either of these and achieve the same result

        >>> index - index2
        >>> index.diff(index2)

        Returns
        -------
        diff : Index
        """

        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if self.equals(other):
            return Index([], name=self.name)

        if not isinstance(other, Index):
            other = np.asarray(other)
            result_name = self.name
        else:
            result_name = self.name if self.name == other.name else None

        theDiff = sorted(set(self) - set(other))
        return Index(theDiff, name=result_name)

    def unique(self):
        """
        Return array of unique values in the Index. Significantly faster than
        numpy.unique

        Returns
        -------
        uniques : ndarray
        """
        from pandas.core.nanops import unique1d
        return unique1d(self.values)

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int if unique index, possibly slice or mask if not
        """
        return self._engine.get_loc(key)

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        try:
            return self._engine.get_value(series, key)
        except KeyError, e1:
            if len(self) > 0 and self.inferred_type == 'integer':
                raise

            try:
                return tslib.get_value_box(series, key)
            except IndexError:
                raise
            except TypeError:
                # generator/iterator-like
                if com.is_iterator(key):
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
        self._engine.set_value(arr, key, value)

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
        num = self._get_level_number(level)
        return self

    def get_indexer(self, target, method=None, limit=None):
        """
        Compute indexer and mask for new index given the current index. The
        indexer should be then used as an input to ndarray.take to align the
        current data to the new index. The mask determines whether labels are
        found or not in the current index

        Parameters
        ----------
        target : Index
        method : {'pad', 'ffill', 'backfill', 'bfill'}
            pad / ffill: propagate LAST valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Notes
        -----
        This is a low-level method and probably should be used at your own risk

        Examples
        --------
        >>> indexer = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)

        Returns
        -------
        indexer : ndarray
        """
        method = self._get_method(method)
        target = _ensure_index(target)

        pself, ptarget = self._possibly_promote(target)
        if pself is not self or ptarget is not target:
            return pself.get_indexer(ptarget, method=method, limit=limit)

        if self.dtype != target.dtype:
            this = self.astype(object)
            target = target.astype(object)
            return this.get_indexer(target, method=method, limit=limit)

        if not self.is_unique:
            raise Exception('Reindexing only valid with uniquely valued Index '
                            'objects')

        if method == 'pad':
            if not self.is_monotonic:
                raise AssertionError('Must be monotonic for forward fill')
            indexer = self._engine.get_pad_indexer(target.values, limit)
        elif method == 'backfill':
            if not self.is_monotonic:
                raise AssertionError('Must be monotonic for backward fill')
            indexer = self._engine.get_backfill_indexer(target.values, limit)
        elif method is None:
            indexer = self._engine.get_indexer(target.values)
        else:
            raise ValueError('unrecognized method: %s' % method)

        return com._ensure_platform_int(indexer)

    def get_indexer_non_unique(self, target, **kwargs):
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
            tgt_values = target.values

        indexer, missing = self._engine.get_indexer_non_unique(tgt_values)
        return Index(indexer), missing

    def _possibly_promote(self, other):
        # A hack, but it works
        from pandas.tseries.index import DatetimeIndex
        if self.inferred_type == 'date' and isinstance(other, DatetimeIndex):
            return DatetimeIndex(self), other
        return self, other

    def groupby(self, to_groupby):
        return self._groupby(self.values, to_groupby)

    def map(self, mapper):
        return self._arrmap(self.values, mapper)

    def isin(self, values):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values

        Parameters
        ----------
        values : set or sequence of values

        Returns
        -------
        is_contained : ndarray (boolean dtype)
        """
        value_set = set(values)
        return lib.ismember(self._array_values(), value_set)

    def _array_values(self):
        return self

    def _get_method(self, method):
        if method:
            method = method.lower()

        aliases = {
            'ffill': 'pad',
            'bfill': 'backfill'
        }
        return aliases.get(method, method)

    def reindex(self, target, method=None, level=None, limit=None,
                copy_if_needed=False, takeable=False):
        """
        For Index, simply returns the new index and the results of
        get_indexer. Provided here to enable an interface that is amenable for
        subclasses of Index whose internals are different (like MultiIndex)

        Returns
        -------
        (new_index, indexer, mask) : tuple
        """
        target = _ensure_index(target)
        if level is not None:
            if method is not None:
                raise ValueError('Fill method not supported if level passed')
            _, indexer, _ = self._join_level(target, level, how='right',
                                             return_indexers=True)
        else:

            if self.equals(target):
                indexer = None

                # to avoid aliasing an existing index
                if copy_if_needed and target.name != self.name and self.name is not None:
                    if target.name is None:
                        target = self.copy()

            else:

                if takeable:
                    if method is not None or limit is not None:
                        raise ValueError("cannot do a takeable reindex with "
                                         "with a method or limit")
                    return self[target], target

                if self.is_unique:
                    indexer = self.get_indexer(target, method=method,
                                               limit=limit)
                else:
                    if method is not None or limit is not None:
                        raise ValueError("cannot reindex a non-unique index "
                                         "with a method or limit")
                    indexer, missing = self.get_indexer_non_unique(target)

        return target, indexer

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        Internal API method. Compute join_index and indexers to conform data
        structures to the new index.

        Parameters
        ----------
        other : Index
        how : {'left', 'right', 'inner', 'outer'}
        level :
        return_indexers : boolean, default False

        Returns
        -------
        join_index, (left_indexer, right_indexer)
        """
        if (level is not None and (isinstance(self, MultiIndex) or
                                   isinstance(other, MultiIndex))):
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

        if self.dtype != other.dtype:
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

    def _join_non_unique(self, other, how='left', return_indexers=False):
        from pandas.tools.merge import _get_join_indexers

        left_idx, right_idx = _get_join_indexers([self.values], [other.values],
                                                 how=how, sort=True)

        left_idx = com._ensure_platform_int(left_idx)
        right_idx = com._ensure_platform_int(right_idx)

        join_index = self.values.take(left_idx)
        mask = left_idx == -1
        np.putmask(join_index, mask, other.values.take(right_idx))

        join_index = self._wrap_joined_index(join_index, other)

        if return_indexers:
            return join_index, left_idx, right_idx
        else:
            return join_index

    def _join_level(self, other, level, how='left', return_indexers=False):
        """
        The join method *only* affects the level of the resulting
        MultiIndex. Otherwise it just exactly aligns the Index data to the
        labels of the level in the MultiIndex. The order of the data indexed by
        the MultiIndex will not be changed (currently)
        """
        if isinstance(self, MultiIndex) and isinstance(other, MultiIndex):
            raise Exception('Join on level between two MultiIndex objects '
                            'is ambiguous')

        left, right = self, other

        flip_order = not isinstance(self, MultiIndex)
        if flip_order:
            left, right = right, left
            how = {'right': 'left', 'left': 'right'}.get(how, how)

        level = left._get_level_number(level)
        old_level = left.levels[level]

        new_level, left_lev_indexer, right_lev_indexer = \
            old_level.join(right, how=how, return_indexers=True)

        if left_lev_indexer is not None:
            left_lev_indexer = com._ensure_int64(left_lev_indexer)
            rev_indexer = lib.get_reverse_indexer(left_lev_indexer,
                                                  len(old_level))

            new_lev_labels = com.take_nd(rev_indexer, left.labels[level],
                                         allow_fill=False)
            omit_mask = new_lev_labels != -1

            new_labels = list(left.labels)
            new_labels[level] = new_lev_labels

            if not omit_mask.all():
                new_labels = [lab[omit_mask] for lab in new_labels]

            new_levels = list(left.levels)
            new_levels[level] = new_level

            join_index = MultiIndex(levels=new_levels, labels=new_labels,
                                    names=left.names)
            left_indexer = np.arange(len(left))[new_lev_labels != -1]
        else:
            join_index = left
            left_indexer = None

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
        ov = other.values

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
                join_index, ridx, lidx = self._left_indexer(other, self)
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

    def slice_indexer(self, start=None, end=None, step=None):
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

        Returns
        -------
        indexer : ndarray or slice

        Notes
        -----
        This function assumes that the data is sorted, so use at your own peril
        """
        start_slice, end_slice = self.slice_locs(start, end)

        # return a slice
        if np.isscalar(start_slice) and np.isscalar(end_slice):
            return slice(start_slice, end_slice, step)

        # loc indexers
        return Index(start_slice) & Index(end_slice)

    def slice_locs(self, start=None, end=None):
        """
        For an ordered Index, compute the slice locations for input labels

        Parameters
        ----------
        start : label, default None
            If None, defaults to the beginning
        end : label, default None
            If None, defaults to the end

        Returns
        -------
        (start, end) : (int, int)

        Notes
        -----
        This function assumes that the data is sorted, so use at your own peril
        """

        is_unique = self.is_unique
        if start is None:
            start_slice = 0
        else:
            try:
                start_slice = self.get_loc(start)

                if not is_unique:

                    # get_loc will return a boolean array for non_uniques
                    # if we are not monotonic
                    if isinstance(start_slice,np.ndarray):
                        raise KeyError("cannot peform a slice operation "
                                       "on a non-unique non-monotonic index")

                if isinstance(start_slice, slice):
                    start_slice = start_slice.start

            except KeyError:
                if self.is_monotonic:
                    start_slice = self.searchsorted(start, side='left')
                else:
                    raise

        if end is None:
            end_slice = len(self)
        else:
            try:
                end_slice = self.get_loc(end)

                if not is_unique:

                    # get_loc will return a boolean array for non_uniques
                    if isinstance(end_slice,np.ndarray):
                        raise KeyError("cannot perform a slice operation "
                                       "on a non-unique non-monotonic index")

                if isinstance(end_slice, slice):
                    end_slice = end_slice.stop
                else:
                    end_slice += 1

            except KeyError:
                if self.is_monotonic:
                    end_slice = self.searchsorted(end, side='right')
                else:
                    raise

        return start_slice, end_slice

    def delete(self, loc):
        """
        Make new Index with passed location deleted

        Returns
        -------
        new_index : Index
        """
        arr = np.delete(self.values, loc)
        return Index(arr)

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : Index
        """
        index = np.asarray(self)
        # because numpy is fussy with tuples
        item_idx = Index([item], dtype=index.dtype)
        new_index = np.concatenate((index[:loc], item_idx, index[loc:]))
        return Index(new_index, name=self.name)

    def drop(self, labels):
        """
        Make new Index with passed list of labels deleted

        Parameters
        ----------
        labels : array-like

        Returns
        -------
        dropped : Index
        """
        labels = com._index_labels_to_array(labels)
        indexer = self.get_indexer(labels)
        mask = indexer == -1
        if mask.any():
            raise ValueError('labels %s not contained in axis' % labels[mask])
        return self.delete(indexer)


class Int64Index(Index):
    """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects. Int64Index is a special case of `Index`
    with purely integer labels. This is the default index type used by the DataFrame
    and Series ctors when no explicit index is provided by the user.

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: object)
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index

    Note
    ----
    An Index instance can **only** contain hashable objects
    """

    _groupby = _algos.groupby_int64
    _arrmap = _algos.arrmap_int64
    _left_indexer_unique = _algos.left_join_indexer_unique_int64
    _left_indexer = _algos.left_join_indexer_int64
    _inner_indexer = _algos.inner_join_indexer_int64
    _outer_indexer = _algos.outer_join_indexer_int64

    _engine_type = _index.Int64Engine

    def __new__(cls, data, dtype=None, copy=False, name=None):
        if not isinstance(data, np.ndarray):
            if np.isscalar(data):
                raise ValueError('Index(...) must be called with a collection '
                                 'of some kind, %s was passed' % repr(data))

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)
            data = np.asarray(data)

        if issubclass(data.dtype.type, basestring):
            raise TypeError('String dtype not supported, you may need '
                            'to explicitly cast to int')
        elif issubclass(data.dtype.type, np.integer):
            # don't force the upcast as we may be dealing
            # with a platform int
            if dtype is None or not issubclass(np.dtype(dtype).type, np.integer):
                dtype = np.int64

            subarr = np.array(data, dtype=dtype, copy=copy)
        else:
            subarr = np.array(data, dtype=np.int64, copy=copy)
            if len(data) > 0:
                if (subarr != data).any():
                    raise TypeError('Unsafe NumPy casting, you must '
                                    'explicitly cast')

        subarr = subarr.view(cls)
        subarr.name = name
        return subarr

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
        if self is other:
            return True

        # if not isinstance(other, Int64Index):
        #     return False

        try:
            return np.array_equal(self, other)
        except TypeError:
            # e.g. fails in numpy 1.6 with DatetimeIndex #1681
            return False

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        return Int64Index(joined, name=name)


class MultiIndex(Index):
    """
    Implements multi-level, a.k.a. hierarchical, index object for pandas
    objects

    Parameters
    ----------
    levels : list or tuple of arrays
        The unique labels for each level
    labels : list or tuple of arrays
        Integers for each level designating which label at each location
    sortorder : optional int
        Level of sortedness (must be lexicographically sorted by that
        level)
    names : optional sequence of objects
        Names for each of the index levels.
    """
    # shadow property
    names = None

    def __new__(cls, levels=None, labels=None, sortorder=None, names=None):
        if len(levels) != len(labels):
            raise AssertionError(
                'Length of levels and labels must be the same')
        if len(levels) == 0:
            raise Exception('Must pass non-zero number of levels/labels')

        if len(levels) == 1:
            if names:
                name = names[0]
            else:
                name = None

            return Index(levels[0], name=name).take(labels[0])

        levels = [_ensure_index(lev) for lev in levels]
        labels = [np.asarray(labs, dtype=np.int_) for labs in labels]

        # v3, 0.8.0
        subarr = np.empty(0, dtype=object).view(cls)
        subarr.levels = levels
        subarr.labels = labels

        if names is None:
            subarr.names = [None] * subarr.nlevels
        else:
            if len(names) != subarr.nlevels:
                raise AssertionError(('Length of names (%d) must be same as level '
                                      '(%d)') % (len(names),subarr.nlevels))

            subarr.names = list(names)

        # set the name
        for i, name in enumerate(subarr.names):
            subarr.levels[i].name = name

        if sortorder is not None:
            subarr.sortorder = int(sortorder)
        else:
            subarr.sortorder = sortorder

        return subarr

    def __array_finalize__(self, obj):
        """
        Update custom MultiIndex attributes when a new array is created by
        numpy, e.g. when calling ndarray.view()
        """
        if not isinstance(obj, type(self)):
            # Only relevant if this array is being created from an Index
            # instance.
            return

        self.levels = list(getattr(obj, 'levels', []))
        self.labels = list(getattr(obj, 'labels', []))
        self.names = list(getattr(obj, 'names', []))
        self.sortorder = getattr(obj, 'sortorder', None)

    def _array_values(self):
        # hack for various methods
        return self.values

    @cache_readonly
    def dtype(self):
        return np.dtype('O')

    def __unicode__(self):
        """
        Return a string representation for a particular Index

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both py2/py3.
        """
        output = 'MultiIndex\n%s'

        summary = com.pprint_thing(self, escape_chars=('\t', '\r', '\n'),
                                   quote_strings=True)

        return output % summary

    def __len__(self):
        return len(self.labels[0])

    def _format_native_types(self, **kwargs):
        return self.tolist()

    @property
    def _constructor(self):
        return MultiIndex.from_tuples

    @cache_readonly
    def inferred_type(self):
        return 'mixed'

    @staticmethod
    def _from_elements(values, labels=None, levels=None, names=None,
                       sortorder=None):
        index = values.view(MultiIndex)
        index.levels = levels
        index.labels = labels
        index.names = names
        index.sortorder = sortorder
        return index

    def _get_level_number(self, level):
        try:
            count = self.names.count(level)
            if count > 1:
                raise Exception('The name %s occurs multiple times, use a '
                                'level number' % level)
            level = self.names.index(level)
        except ValueError:
            if not isinstance(level, int):
                raise Exception('Level %s not found' % str(level))
            elif level < 0:
                level += self.nlevels
            # Note: levels are zero-based
            elif level >= self.nlevels:
                raise ValueError('Index has only %d levels, not %d'
                                 % (self.nlevels, level + 1))
        return level

    _tuples = None

    @property
    def values(self):
        if self._is_v2:
            return self.view(np.ndarray)
        else:
            if self._tuples is not None:
                return self._tuples

            values = []
            for lev, lab in zip(self.levels, self.labels):
                taken = com.take_1d(lev.values, lab)
                # Need to box timestamps, etc.
                if hasattr(lev, '_box_values'):
                    taken = lev._box_values(taken)
                values.append(taken)

            self._tuples = lib.fast_zip(values)
            return self._tuples

    # fml
    @property
    def _is_v1(self):
        contents = self.view(np.ndarray)
        return len(contents) > 0 and not isinstance(contents[0], tuple)

    @property
    def _is_v2(self):
        contents = self.view(np.ndarray)
        return len(contents) > 0 and isinstance(contents[0], tuple)

    @property
    def _has_complex_internals(self):
        # to disable groupby tricks
        return True

    @property
    def has_duplicates(self):
        """
        Return True if there are no unique groups
        """
        # has duplicates
        shape = [len(lev) for lev in self.levels]
        group_index = np.zeros(len(self), dtype='i8')
        for i in xrange(len(shape)):
            stride = np.prod([x for x in shape[i + 1:]], dtype='i8')
            group_index += self.labels[i] * stride

        if len(np.unique(group_index)) < len(group_index):
            return True

        return False

    def get_value(self, series, key):
        # somewhat broken encapsulation
        from pandas.core.indexing import _maybe_droplevels
        from pandas.core.series import Series

        # Label-based
        try:
            return self._engine.get_value(series, key)
        except KeyError, e1:
            try:
                # TODO: what if a level contains tuples??
                loc = self.get_loc(key)
                new_values = series.values[loc]
                new_index = self[loc]
                new_index = _maybe_droplevels(new_index, key)
                return Series(new_values, index=new_index, name=series.name)
            except KeyError:
                pass

            try:
                return _index.get_value_at(series, key)
            except IndexError:
                raise
            except TypeError:
                # generator/iterator-like
                if com.is_iterator(key):
                    raise InvalidIndexError(key)
                else:
                    raise e1
            except Exception:  # pragma: no cover
                raise e1
        except TypeError:
            raise InvalidIndexError(key)

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
        num = self._get_level_number(level)
        unique_vals = self.levels[num]  # .values
        labels = self.labels[num]
        values = unique_vals.take(labels)
        values.name = self.names[num]
        return values

    def format(self, space=2, sparsify=None, adjoin=True, names=False,
               na_rep='NaN', formatter=None):
        if len(self) == 0:
            return []

        stringified_levels = []
        for lev, lab in zip(self.levels, self.labels):
            if len(lev) > 0:

                formatted = lev.take(lab).format(formatter=formatter)

                # we have some NA
                mask = lab==-1
                if mask.any():
                    formatted = np.array(formatted,dtype=object)
                    formatted[mask] = na_rep
                    formatted = formatted.tolist()

            else:
                # weird all NA case
                formatted = [com.pprint_thing(na_rep if isnull(x) else x, escape_chars=('\t', '\r', '\n'))
                             for x in com.take_1d(lev.values, lab)]
            stringified_levels.append(formatted)

        result_levels = []
        for lev, name in zip(stringified_levels, self.names):
            level = []

            if names:
                level.append(com.pprint_thing(name, escape_chars=('\t', '\r', '\n'))
                             if name is not None else '')


            level.extend(np.array(lev, dtype=object))
            result_levels.append(level)

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        if sparsify:
            sentinal = ''
            # GH3547
            # use value of sparsify as sentinal,  unless it's an obvious "Truthey" value
            if sparsify not in [True,1]:
                sentinal = sparsify
            # little bit of a kludge job for #1217
            result_levels = _sparsify(result_levels,
                                      start=int(names),
                                      sentinal=sentinal)

        if adjoin:
            return com.adjoin(space, *result_levels).split('\n')
        else:
            return result_levels

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
        arrays : list / sequence
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level)

        Returns
        -------
        index : MultiIndex
        """
        from pandas.core.categorical import Categorical

        if len(arrays) == 1:
            name = None if names is None else names[0]
            return Index(arrays[0], name=name)

        cats = [Categorical.from_array(arr) for arr in arrays]
        levels = [c.levels for c in cats]
        labels = [c.labels for c in cats]
        if names is None:
            names = [c.name for c in cats]

        return MultiIndex(levels=levels, labels=labels,
                          sortorder=sortorder, names=names)

    @classmethod
    def from_tuples(cls, tuples, sortorder=None, names=None):
        """
        Convert list of tuples to MultiIndex

        Parameters
        ----------
        tuples : array-like
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level)

        Returns
        -------
        index : MultiIndex
        """
        if len(tuples) == 0:
            raise Exception('Cannot infer number of levels from empty list')

        if isinstance(tuples, np.ndarray):
            if isinstance(tuples, Index):
                tuples = tuples.values

            arrays = list(lib.tuples_to_object_array(tuples).T)
        elif isinstance(tuples, list):
            arrays = list(lib.to_object_array_tuples(tuples).T)
        else:
            arrays = zip(*tuples)

        return MultiIndex.from_arrays(arrays, sortorder=sortorder,
                                      names=names)

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
        except KeyError:
            return False

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = (self.levels, self.labels, self.sortorder, self.names)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        np.ndarray.__setstate__(self, nd_state)
        levels, labels, sortorder, names = own_state

        self.levels = [Index(x) for x in levels]
        self.labels = labels
        self.names = names
        self.sortorder = sortorder

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
            if com._is_bool_indexer(key):
                key = np.asarray(key)
                sortorder = self.sortorder
            else:
                # cannot be sure whether the result will be sorted
                sortorder = None

            result = np.empty(0, dtype=object).view(type(self))
            new_labels = [lab[key] for lab in self.labels]

            # an optimization
            result.levels = list(self.levels)
            result.labels = new_labels
            result.sortorder = sortorder
            result.names = self.names

            return result

    def take(self, indexer, axis=None):
        """
        Analogous to ndarray.take
        """
        indexer = com._ensure_platform_int(indexer)
        new_labels = [lab.take(indexer) for lab in self.labels]
        return MultiIndex(levels=self.levels, labels=new_labels,
                          names=self.names)

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

        to_concat = (self.values,) + tuple(k.values for k in other)
        new_tuples = np.concatenate(to_concat)

        # if all(isinstance(x, MultiIndex) for x in other):
        try:
            return MultiIndex.from_tuples(new_tuples, names=self.names)
        except:
            return Index(new_tuples)

    def argsort(self, *args, **kwargs):
        return self.values.argsort()

    def drop(self, labels, level=None):
        """
        Make new MultiIndex with passed list of labels deleted

        Parameters
        ----------
        labels : array-like
            Must be a list of tuples
        level : int or name, default None

        Returns
        -------
        dropped : MultiIndex
        """
        if level is not None:
            return self._drop_from_level(labels, level)

        try:
            if not isinstance(labels, np.ndarray):
                labels = com._index_labels_to_array(labels)
            indexer = self.get_indexer(labels)
            mask = indexer == -1
            if mask.any():
                raise ValueError('labels %s not contained in axis'
                                 % labels[mask])
            return self.delete(indexer)
        except Exception:
            pass

        inds = []
        for label in labels:
            loc = self.get_loc(label)
            if isinstance(loc, int):
                inds.append(loc)
            else:
                inds.extend(range(loc.start, loc.stop))

        return self.delete(inds)

    def _drop_from_level(self, labels, level):
        labels = com._index_labels_to_array(labels)
        i = self._get_level_number(level)
        index = self.levels[i]
        values = index.get_indexer(labels)

        mask = -lib.ismember(self.labels[i], set(values))

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
            result = new_levels[0].take(new_labels[0])
            result.name = new_names[0]
            return result
        else:
            return MultiIndex(levels=new_levels, labels=new_labels,
                              names=new_names)

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
                          names=new_names)

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
                          names=new_names)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def sortlevel(self, level=0, ascending=True):
        """
        Sort MultiIndex at the requested level. The result will respect the
        original ordering of the associated factor at that level.

        Parameters
        ----------
        level : int or str, default 0
            If a string is given, must be a name of the level
        ascending : boolean, default True
            False to sort in descending order

        Returns
        -------
        sorted_index : MultiIndex
        """
        from pandas.core.groupby import _indexer_from_factorized

        labels = list(self.labels)

        level = self._get_level_number(level)
        primary = labels.pop(level)

        shape = list(self.levshape)
        primshp = shape.pop(level)

        indexer = _indexer_from_factorized((primary,) + tuple(labels),
                                           (primshp,) + tuple(shape),
                                           compress=False)
        if not ascending:
            indexer = indexer[::-1]

        indexer = com._ensure_platform_int(indexer)
        new_labels = [lab.take(indexer) for lab in self.labels]

        new_index = MultiIndex(labels=new_labels, levels=self.levels,
                               names=self.names, sortorder=level)

        return new_index, indexer

    def get_indexer(self, target, method=None, limit=None):
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
        method = self._get_method(method)

        target = _ensure_index(target)

        target_index = target
        if isinstance(target, MultiIndex):
            target_index = target._tuple_index

        if target_index.dtype != object:
            return np.ones(len(target_index)) * -1

        if not self.is_unique:
            raise Exception('Reindexing only valid with uniquely valued Index '
                            'objects')

        self_index = self._tuple_index

        if method == 'pad':
            if not self.is_unique or not self.is_monotonic:
                raise AssertionError(('Must be unique and monotonic to '
                                      'use forward fill getting the indexer'))
            indexer = self_index._engine.get_pad_indexer(target_index,
                                                         limit=limit)
        elif method == 'backfill':
            if not self.is_unique or not self.is_monotonic:
                raise AssertionError(('Must be unique and monotonic to '
                                      'use backward fill getting the indexer'))
            indexer = self_index._engine.get_backfill_indexer(target_index,
                                                              limit=limit)
        else:
            indexer = self_index._engine.get_indexer(target_index)

        return com._ensure_platform_int(indexer)

    def reindex(self, target, method=None, level=None, limit=None,
                copy_if_needed=False, takeable=False):
        """
        Performs any necessary conversion on the input index and calls
        get_indexer. This method is here so MultiIndex and an Index of
        like-labeled tuples can play nice together

        Returns
        -------
        (new_index, indexer, mask) : (MultiIndex, ndarray, ndarray)
        """
        if level is not None:
            if method is not None:
                raise ValueError('Fill method not supported if level passed')
            target, indexer, _ = self._join_level(target, level, how='right',
                                                  return_indexers=True)
        else:
            if self.equals(target):
                indexer = None
            else:
                indexer = self.get_indexer(target, method=method,
                                           limit=limit)

        if not isinstance(target, MultiIndex):
            if indexer is None:
                target = self
            elif (indexer >= 0).all():
                target = self.take(indexer)
            else:
                # hopefully?
                target = MultiIndex.from_tuples(target)

        return target, indexer

    @cache_readonly
    def _tuple_index(self):
        """
        Convert MultiIndex to an Index of tuples

        Returns
        -------
        index : Index
        """
        return Index(self.values)

    def slice_locs(self, start=None, end=None, strict=False):
        """
        For an ordered MultiIndex, compute the slice locations for input
        labels. They can tuples representing partial levels, e.g. for a
        MultiIndex with 3 levels, you can pass a single value (corresponding to
        the first level), or a 1-, 2-, or 3-tuple.

        Parameters
        ----------
        start : label or tuple, default None
            If None, defaults to the beginning
        end : label or tuple
            If None, defaults to the end
        strict : boolean,

        Returns
        -------
        (start, end) : (int, int)

        Notes
        -----
        This function assumes that the data is sorted by the first level
        """
        if start is None:
            start_slice = 0
        else:
            if not isinstance(start, tuple):
                start = start,
            start_slice = self._partial_tup_index(start, side='left')

        if end is None:
            end_slice = len(self)
        else:
            if not isinstance(end, tuple):
                end = end,
            end_slice = self._partial_tup_index(end, side='right')

        return start_slice, end_slice

    def _partial_tup_index(self, tup, side='left'):
        if len(tup) > self.lexsort_depth:
            raise KeyError('MultiIndex lexsort depth %d, key was length %d' %
                           (self.lexsort_depth, len(tup)))

        n = len(tup)
        start, end = 0, len(self)
        zipped = izip(tup, self.levels, self.labels)
        for k, (lab, lev, labs) in enumerate(zipped):
            section = labs[start:end]

            if lab not in lev:
                if not lev.is_type_compatible(lib.infer_dtype([lab])):
                    raise Exception('Level type mismatch: %s' % lab)

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

    def get_loc(self, key):
        """
        Get integer location slice for requested label or tuple

        Parameters
        ----------
        key : label or tuple

        Returns
        -------
        loc : int or slice object
        """
        if isinstance(key, tuple):
            if len(key) == self.nlevels:
                if self.is_unique:
                    return self._engine.get_loc(key)
                else:
                    return slice(*self.slice_locs(key, key))
            else:
                # partial selection
                result = slice(*self.slice_locs(key, key))
                if result.start == result.stop:
                    raise KeyError(key)
                return result
        else:
            return self._get_level_indexer(key, level=0)

    def get_loc_level(self, key, level=0):
        """
        Get integer location slice for requested label or tuple

        Parameters
        ----------
        key : label or tuple

        Returns
        -------
        loc : int or slice object
        """
        def _drop_levels(indexer, levels):
            # kludgearound
            new_index = self[indexer]
            levels = [self._get_level_number(i) for i in levels]
            for i in sorted(levels, reverse=True):
                new_index = new_index.droplevel(i)
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
            return result, _drop_levels(result, level)

        level = self._get_level_number(level)

        # kludge for #1796
        if isinstance(key, list):
            key = tuple(key)

        if isinstance(key, tuple) and level == 0:

            try:
                if key in self.levels[0]:
                    indexer = self._get_level_indexer(key, level=level)
                    new_index = _drop_levels(indexer, [0])
                    return indexer, new_index
            except TypeError:
                pass

            if not any(isinstance(k, slice) for k in key):
                if len(key) == self.nlevels:
                    if self.is_unique:
                        return self._engine.get_loc(key), None
                    else:
                        indexer = slice(*self.slice_locs(key, key))
                        return indexer, self[indexer]
                else:
                    # partial selection
                    indexer = slice(*self.slice_locs(key, key))
                    if indexer.start == indexer.stop:
                        raise KeyError(key)
                    ilevels = [i for i in range(len(key))
                               if key[i] != slice(None, None)]
                    return indexer, _drop_levels(indexer, ilevels)
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
                return indexer, _drop_levels(indexer, ilevels)
        else:
            indexer = self._get_level_indexer(key, level=level)
            new_index = _drop_levels(indexer, [level])
            return indexer, new_index

    def _get_level_indexer(self, key, level=0):
        level_index = self.levels[level]
        loc = level_index.get_loc(key)
        labels = self.labels[level]

        if level > 0 or self.lexsort_depth == 0:
            return labels == loc
        else:
            # sorted, so can return slice object -> view
            i = labels.searchsorted(loc, side='left')
            j = labels.searchsorted(loc, side='right')
            return slice(i, j)

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

        return MultiIndex(levels=new_levels, labels=new_labels)

    def equals(self, other):
        """
        Determines if two MultiIndex objects have the same labeling information
        (the levels themselves do not necessarily have to be the same)

        See also
        --------
        equal_levels
        """
        if self is other:
            return True

        if not isinstance(other, MultiIndex):
            return np.array_equal(self.values, _ensure_index(other))

        if self.nlevels != other.nlevels:
            return False

        if len(self) != len(other):
            return False

        for i in xrange(self.nlevels):
            svalues = com.take_nd(self.levels[i].values, self.labels[i],
                                  allow_fill=False)
            ovalues = com.take_nd(other.levels[i].values, other.labels[i],
                                  allow_fill=False)
            if not np.array_equal(svalues, ovalues):
                return False

        return True

    def equal_levels(self, other):
        """
        Return True if the levels of both MultiIndex objects are the same

        """
        if self.nlevels != other.nlevels:
            return False

        for i in xrange(self.nlevels):
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
        """
        self._assert_can_do_setop(other)

        if len(other) == 0 or self.equals(other):
            return self

        result_names = self.names if self.names == other.names else None

        uniq_tuples = lib.fast_unique_multiple([self.values, other.values])
        return MultiIndex.from_arrays(zip(*uniq_tuples), sortorder=0,
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

        if self.equals(other):
            return self

        result_names = self.names if self.names == other.names else None

        self_tuples = self.values
        other_tuples = other.values
        uniq_tuples = sorted(set(self_tuples) & set(other_tuples))
        if len(uniq_tuples) == 0:
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names)
        else:
            return MultiIndex.from_arrays(zip(*uniq_tuples), sortorder=0,
                                          names=result_names)

    def diff(self, other):
        """
        Compute sorted set difference of two MultiIndex objects

        Returns
        -------
        diff : MultiIndex
        """
        self._assert_can_do_setop(other)

        if not isinstance(other, MultiIndex):
            if len(other) == 0:
                return self
            try:
                other = MultiIndex.from_tuples(other)
            except:
                raise TypeError("other should be a MultiIndex or a list of tuples")
            result_names = self.names
        else:
            result_names = self.names if self.names == other.names else None

        if self.equals(other):
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names)

        difference = sorted(set(self.values) - set(other.values))

        if len(difference) == 0:
            return MultiIndex(levels=[[]] * self.nlevels,
                              labels=[[]] * self.nlevels,
                              names=result_names)
        else:
            return MultiIndex.from_tuples(difference, sortorder=0,
                                          names=result_names)

    def _assert_can_do_setop(self, other):
        pass

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
            raise ValueError('Passed item incompatible tuple length')

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
            new_labels.append(np.insert(labels, loc, lev_loc))

        return MultiIndex(levels=new_levels, labels=new_labels,
                          names=self.names)

    def delete(self, loc):
        """
        Make new index with passed location deleted

        Returns
        -------
        new_index : MultiIndex
        """
        new_labels = [np.delete(lab, loc) for lab in self.labels]
        return MultiIndex(levels=self.levels, labels=new_labels,
                          names=self.names)

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


# For utility purposes

def _sparsify(label_list, start=0,sentinal=''):
    pivoted = zip(*label_list)
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
                sparse_cur.append(sentinal)
            else:
                sparse_cur.extend(cur[i:])
                result.append(sparse_cur)
                break

        prev = cur

    return zip(*result)


def _ensure_index(index_like):
    if isinstance(index_like, Index):
        return index_like
    if hasattr(index_like, 'name'):
        return Index(index_like, name=index_like.name)

    if isinstance(index_like, list):
        # #2200 ?
        converted, all_arrays = lib.clean_index_list(index_like)

        if len(converted) > 0 and all_arrays:
            return MultiIndex.from_arrays(converted)
        else:
            index_like = converted

    return Index(index_like)


def _validate_join_method(method):
    if method not in ['left', 'right', 'inner', 'outer']:
        raise Exception('do not recognize join method %s' % method)


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
    return dict((id(x), x) for x in indexes).values()


def _union_indexes(indexes):
    if len(indexes) == 0:
        raise AssertionError('Must have at least 1 Index to union')
    if len(indexes) == 1:
        result = indexes[0]
        if isinstance(result, list):
            result = Index(sorted(result))
        return result

    indexes, kind = _sanitize_and_check(indexes)

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
                return Index(lib.fast_unique_multiple(indexes))

        return index
    else:
        return Index(lib.fast_unique_multiple_list(indexes))


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


def _handle_legacy_indexes(indexes):
    from pandas.core.daterange import DateRange
    from pandas.tseries.index import DatetimeIndex

    converted = []
    for index in indexes:
        if isinstance(index, DateRange):
            if len(index) == 0:
                kwds = dict(data=[], freq=index.offset, tz=index.tzinfo)
            else:
                kwds = dict(start=index[0], end=index[-1],
                            freq=index.offset, tz=index.tzinfo)

            index = DatetimeIndex(**kwds)

        converted.append(index)

    return converted


def _get_consensus_names(indexes):

    # find the non-none names, need to tupleify to make
    # the set hashable, then reverse on return
    consensus_names = set([ tuple(i.names) for i in indexes if all(n is not None for n in i.names) ])
    if len(consensus_names) == 1:
        return list(list(consensus_names)[0])
    return [None] * indexes[0].nlevels


def _maybe_box(idx):
    from pandas.tseries.api import DatetimeIndex, PeriodIndex
    klasses = DatetimeIndex, PeriodIndex

    if isinstance(idx, klasses):
        return idx.asobject
    return idx


def _all_indexes_same(indexes):
    first = indexes[0]
    for index in indexes[1:]:
        if not first.equals(index):
            return False
    return True
