# pylint: disable=E1101,E1103,W0232

from datetime import time
from itertools import izip

import numpy as np

from pandas.util.decorators import cache_readonly
import pandas.core.common as com
import pandas._tseries as lib
import pandas._engines as _gin

__all__ = ['Index']


def _indexOp(opname):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        func = getattr(self.view(np.ndarray), opname)
        return func(other)
    return wrapper


class InvalidIndexError(Exception):
    pass


class Index(np.ndarray):
    """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: object)
    copy : bool
        Make a copy of input ndarray

    Note
    ----
    An Index instance can **only** contain hashable objects
    """
    # Cython methods
    _map_indices = lib.map_indices_object
    _is_monotonic = lib.is_monotonic_object
    _groupby = lib.groupby_object
    _arrmap = lib.arrmap_object
    _left_indexer = lib.left_join_indexer_object
    _inner_indexer = lib.inner_join_indexer_object
    _outer_indexer = lib.outer_join_indexer_object
    _merge_indexer = lib.merge_indexer_object
    _pad = lib.pad_object
    _backfill = lib.backfill_object

    name = None

    def __new__(cls, data, dtype=None, copy=False, name=None):
        if isinstance(data, np.ndarray):
            if dtype is None and issubclass(data.dtype.type, np.integer):
                return Int64Index(data, copy=copy, name=name)

            subarr = np.array(data, dtype=object, copy=copy)
        elif np.isscalar(data):
            raise ValueError('Index(...) must be called with a collection '
                             'of some kind, %s was passed' % repr(data))
        else:
            # other iterable of some kind
            subarr = com._asarray_tuplesafe(data, dtype=object)

        if lib.is_integer_array(subarr) and dtype is None:
            return Int64Index(subarr.astype('i8'), name=name)

        subarr = subarr.view(cls)
        subarr.name = name
        return subarr

    def __array_finalize__(self, obj):
        self.name = getattr(obj, 'name', None)

    def astype(self, dtype):
        return Index(self.values.astype(dtype), name=self.name,
                     dtype=dtype)

    @property
    def dtype(self):
        return self.values.dtype

    @property
    def nlevels(self):
        return 1

    # for compat with multindex code

    def _get_names(self):
        return [self.name]

    def _set_names(self, values):
        assert(len(values) == 1)
        self.name = values[0]

    names = property(fset=_set_names, fget=_get_names)

    @property
    def _constructor(self):
        return Index

    @property
    def _has_complex_internals(self):
        # to disable groupby tricks in MultiIndex
        return False

    def summary(self):
        if len(self) > 0:
            index_summary = ', %s to %s' % (str(self[0]), str(self[-1]))
        else:
            index_summary = ''

        name = type(self).__name__
        return '%s: %s entries%s' % (name, len(self), index_summary)

    def __str__(self):
        try:
            return np.array_repr(self.values)
        except UnicodeError:
            converted = u','.join(unicode(x) for x in self.values)
            return u'%s([%s], dtype=''%s'')' % (type(self).__name__, converted,
                                              str(self.values.dtype))

    @property
    def values(self):
        return np.asarray(self)

    @cache_readonly
    def is_monotonic(self):
        try:
            # wrong buffer type raises ValueError
            return self._is_monotonic(self.values)
        except TypeError:
            return False

    def get_duplicates(self):
        from collections import defaultdict
        counter = defaultdict(lambda: 0)
        for k in self.values:
            counter[k] += 1
        return sorted(k for k, v in counter.iteritems() if v > 1)

    _get_duplicates = get_duplicates

    @property
    def indexMap(self):
        "{label -> location}"
        return self._engine.get_mapping(1)

    def _cleanup(self):
        self._engine.clear_mapping()

    @cache_readonly
    def _engine(self):
        import weakref
        # property, for now, slow to look up
        return _gin.DictIndexEngine(weakref.ref(self),
                                    self._map_indices)

    def _get_level_number(self, level):
        if not isinstance(level, int):
            assert(level == self.name)
            level = 0
        return level

    def _verify_integrity(self):
        return self._engine.has_integrity

    @cache_readonly
    def inferred_type(self):
        return lib.infer_dtype(self)

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
        return key in self._engine

    def __hash__(self):
        return hash(self.view(np.ndarray))

    def __setitem__(self, key, value):
        """Disable the setting of values."""
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
        if isinstance(other, (list, tuple)):
            to_concat = (self.values,) + tuple(other)
        else:
            to_concat = self.values, other.values
        return Index(np.concatenate(to_concat))

    def take(self, *args, **kwargs):
        """
        Analogous to ndarray.take
        """
        taken = self.view(np.ndarray).take(*args, **kwargs)
        return self._constructor(taken, name=self.name)

    def format(self, name=False):
        """
        Render a string representation of the Index
        """
        from pandas.core.format import format_array

        header = []
        if name:
            header.append(str(self.name) if self.name is not None else '')

        if self.is_all_dates:
            zero_time = time(0, 0)
            result = []
            for dt in self:
                if dt.time() != zero_time or dt.tzinfo is not None:
                    return header + ['%s' % x for x in self]
                result.append(dt.strftime("%Y-%m-%d"))
            return header + result

        values = self.values
        if values.dtype == np.object_:
            values = lib.maybe_convert_objects(values)

        if values.dtype == np.object_:
            result = [com._stringify(x) for x in values]
        else:
            result = _trim_front(format_array(values, None, justify='left'))
        return header + result

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
        if label not in self:
            loc = self.searchsorted(label, side='left')
            if loc > 0:
                return self[loc - 1]
            else:
                return np.nan

        return label

    def sort(self, *args, **kwargs):
        raise Exception('Cannot sort an Index object')

    def shift(self, periods, offset):
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

        offset = periods * offset
        return Index([idx + offset for idx in self])

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
                other_diff = other.values.take(indexer)
                result = list(self) + list(other_diff)
                # timsort wins
                try:
                    result.sort()
                except Exception:
                    pass
            else:
                # contained in
                result = sorted(self)

        # for subclasses
        return self._wrap_union_result(other, result)

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        if type(self) == type(other):
            return type(self)(result, name=name)
        else:
            return Index(result, name=name)

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

        other = _ensure_index(other)

        if self.equals(other):
            return self

        if self.dtype != other.dtype:
            this = self.astype('O')
            other = other.astype('O')
            return this.intersection(other)

        if self.is_monotonic and other.is_monotonic:
            result = self._inner_indexer(self, other.values)[0]
            return self._wrap_union_result(other, result)
        else:
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
            return Index([])

        otherArr = np.asarray(other)
        theDiff = sorted(set(self) - set(otherArr))
        return Index(theDiff)

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
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
            if self.inferred_type == 'integer':
                raise

            try:
                return _gin.get_value_at(series, key)
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

    def set_value(self, arr, key, value):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        self._engine.set_value(arr, key, value)

    def get_indexer(self, target, method=None):
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
        >>> indexer, mask = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)
        >>> new_values[-mask] = np.nan

        Returns
        -------
        (indexer, mask) : (ndarray, ndarray)
        """
        method = self._get_method(method)

        target = _ensure_index(target)

        if self.dtype != target.dtype:
            this = Index(self, dtype=object)
            target = Index(target, dtype=object)
            return this.get_indexer(target, method=method)
        # if self.dtype != target.dtype:
        #     target = Index(target, dtype=object)

        if method == 'pad':
            indexer = self._pad(self, target, self.indexMap, target.indexMap)
        elif method == 'backfill':
            indexer = self._backfill(self, target, self.indexMap,
                                     target.indexMap)
        elif method is None:
            indexer = self._get_indexer_standard(target)
        else:
            raise ValueError('unrecognized method: %s' % method)
        return indexer

    def _get_indexer_standard(self, other):
        if (self.dtype != np.object_ and
            self.is_monotonic and other.is_monotonic):
            return self._left_indexer(other, self)
        else:
            return self._merge_indexer(other, self.indexMap)

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
        return lib.ismember(self, value_set)

    def _get_method(self, method):
        if method:
            method = method.lower()

        aliases = {
            'ffill': 'pad',
            'bfill': 'backfill'
        }
        return aliases.get(method, method)

    def reindex(self, target, method=None, level=None):
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
            _, indexer, _ = self._join_level(target, level, how='left',
                                                  return_indexers=True)
        else:
            if self.equals(target):
                indexer = None
            else:
                indexer = self.get_indexer(target, method=method)
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

        if self.dtype != other.dtype:
            this = self.astype('O')
            other = other.astype('O')
            return this.join(other, how=how,
                             return_indexers=return_indexers)

        _validate_join_method(how)

        if self.is_monotonic and other.is_monotonic:
            return self._join_monotonic(other, how=how,
                                        return_indexers=return_indexers)

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
            rev_indexer = lib.get_reverse_indexer(left_lev_indexer,
                                                  len(old_level))

            new_lev_labels = rev_indexer.take(left.labels[level])
            omit_mask = new_lev_labels != -1

            new_labels = list(left.labels)
            new_labels[level] = new_lev_labels

            if not omit_mask.all():
                new_labels = [lab[omit_mask] for lab in new_labels]

            new_levels = list(left.levels)
            new_levels[level] = new_level

            join_index = MultiIndex(levels=new_levels, labels=new_labels,
                                    names=left.names)
        else:
            join_index = left

        left_indexer = None

        if right_lev_indexer is not None:
            right_indexer = right_lev_indexer.take(join_index.labels[level])
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

        if how == 'left':
            join_index = self
            lidx = None
            ridx = self._left_indexer(self, other)
        elif how == 'right':
            join_index = other
            lidx = self._left_indexer(other, self)
            ridx = None
        elif how == 'inner':
            join_index, lidx, ridx = self._inner_indexer(self.values,
                                                         other.values)
            join_index = self._wrap_joined_index(join_index, other)
        elif how == 'outer':
            join_index, lidx, ridx = self._outer_indexer(self.values,
                                                         other.values)
            join_index = self._wrap_joined_index(join_index, other)

        if return_indexers:
            return join_index, lidx, ridx
        else:
            return join_index

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        return Index(joined, name=name)

    def slice_locs(self, start=None, end=None):
        """
        For an ordered Index, compute the slice locations for input labels

        Parameters
        ----------
        start : label, default None
            If None, defaults to the beginning
        end : label
            If None, defaults to the end

        Returns
        -------
        (begin, end) : (int, int)

        Notes
        -----
        This function assumes that the data is sorted, so use at your own peril
        """
        if start is None:
            beg_slice = 0
        else:
            try:
                beg_slice = self.get_loc(start)
            except KeyError:
                if self.is_monotonic:
                    beg_slice = self.searchsorted(start, side='left')
                else:
                    raise

        if end is None:
            end_slice = len(self)
        else:
            try:
                end_slice = self.get_loc(end) + 1
            except KeyError:
                if self.is_monotonic:
                    end_slice = self.searchsorted(end, side='right')
                else:
                    raise

        return beg_slice, end_slice

    def delete(self, loc):
        """
        Make new Index with passed location deleted

        Returns
        -------
        new_index : Index
        """
        arr = np.delete(np.asarray(self), loc)
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
        item_idx = Index([item])
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
        labels = np.asarray(list(labels), dtype=object)
        indexer = self.get_indexer(labels)
        mask = indexer == -1
        if mask.any():
            raise ValueError('labels %s not contained in axis' % labels[mask])
        return self.delete(indexer)

    def copy(self, order='C'):
        """
        Overridden ndarray.copy to copy over attributes

        Returns
        -------
        cp : Index
            Returns view on same base ndarray
        """
        cp = self.view(np.ndarray).view(type(self))
        cp.__dict__.update(self.__dict__)
        return cp


class Int64Index(Index):

    _map_indices = lib.map_indices_int64
    _is_monotonic = lib.is_monotonic_int64
    _groupby = lib.groupby_int64
    _arrmap = lib.arrmap_int64
    _left_indexer = lib.left_join_indexer_int64
    _inner_indexer = lib.inner_join_indexer_int64
    _outer_indexer = lib.outer_join_indexer_int64
    _merge_indexer = lib.merge_indexer_int64
    _pad = lib.pad_int64
    _backfill = lib.backfill_int64

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
            subarr = np.array(data, dtype=np.int64, copy=copy)
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
    def _constructor(self):
        return Int64Index

    @property
    def dtype(self):
        return np.dtype('int64')

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

        return np.array_equal(self, other)

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        return Int64Index(joined, name=name)


class DateIndex(Index):
    pass


class Factor(np.ndarray):
    """
    Represents a categorical variable in classic R / S-plus fashion

    Parameters
    ----------
    data : array-like

    Returns
    -------
    **Attributes**
      * labels : ndarray
      * levels : ndarray
    """
    def __new__(cls, data):
        data = np.asarray(data, dtype=object)
        levels, factor = unique_with_labels(data)
        factor = factor.view(Factor)
        factor.levels = levels
        return factor

    levels = None

    def __array_finalize__(self, obj):
        self.levels = getattr(obj, 'levels', None)

    @property
    def labels(self):
        return self.view(np.ndarray)

    def asarray(self):
        return np.asarray(self.levels).take(self.labels)

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Factor:\n%s\nLevels (%d): %s'
        values = self.asarray()
        return temp % (repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            return self.levels[i]
        else:
            return np.ndarray.__getitem__(self, key)


def unique_with_labels(values):
    rizer = lib.Factorizer(len(values))
    labels, _ = rizer.factorize(values, sort=False)
    uniques = Index(rizer.uniques)

    try:
        sorter = uniques.argsort()
        reverse_indexer = np.empty(len(sorter), dtype='i4')
        reverse_indexer.put(sorter, np.arange(len(sorter)))
        labels = reverse_indexer.take(labels)
        uniques = uniques.take(sorter)
    except TypeError:
        pass

    return uniques, labels


def unique_int64(values):
    if values.dtype != np.int64:
        values = values.astype('i8')

    table = lib.Int64HashTable(len(values))
    uniques = table.unique(values)
    return uniques


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
    """
    # shadow property
    names = None

    def __new__(cls, levels=None, labels=None, sortorder=None, names=None):
        assert(len(levels) == len(labels))
        if len(levels) == 0:
            raise Exception('Must pass non-zero number of levels/labels')

        if len(levels) == 1:
            if names:
                name = names[0]
            else:
                name = None

            return Index(levels[0], name=name).take(labels[0])

        levels = [_ensure_index(lev) for lev in levels]
        labels = [np.asarray(labs, dtype=np.int32) for labs in labels]

        values = [np.asarray(lev).take(lab)
                  for lev, lab in zip(levels, labels)]
        subarr = lib.fast_zip(values).view(cls)

        subarr.levels = levels
        subarr.labels = labels

        if names is None:
            subarr.names = [None] * subarr.nlevels
        else:
            assert(len(names) == subarr.nlevels)
            subarr.names = list(names)

        # set the name
        for i, name in enumerate(subarr.names):
            subarr.levels[i].name = name

        if sortorder is not None:
            subarr.sortorder = int(sortorder)
        else:
            subarr.sortorder = sortorder

        return subarr

    def copy(self, order='C'):
        """
        Overridden ndarray.copy to copy over attributes

        Returns
        -------
        cp : Index
            Returns view on same base ndarray
        """
        cp = self.view(np.ndarray).view(type(self))
        cp.levels = list(self.levels)
        cp.labels = list(self.labels)
        cp.names = list(self.names)
        cp.sortorder = self.sortorder
        return cp

    @property
    def dtype(self):
        return np.dtype('O')

    @property
    def _constructor(self):
        return MultiIndex.from_tuples

    @staticmethod
    def _from_elements(values, labels=None, levels=None, names=None,
                       sortorder=None):
        index = values.view(MultiIndex)
        index.levels = levels
        index.labels = labels
        index.names  = names
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
            elif level >= self.nlevels:
                raise ValueError('Index has only %d levels, not %d'
                                 % (self.nlevels, level))
        return level

    @property
    def values(self):
        if self._is_legacy_format:
            # for legacy MultiIndex
            values = [np.asarray(lev).take(lab)
                      for lev, lab in zip(self.levels, self.labels)]
            return lib.fast_zip(values)
        else:
            return self.view(np.ndarray)

    @property
    def _is_legacy_format(self):
        contents = self.view(np.ndarray)
        return len(contents) > 0 and not isinstance(contents[0], tuple)

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
            stride = np.prod([x for x in shape[i+1:]], dtype='i8')
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
                return _gin.get_value_at(series, key)
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
        unique_vals = self.levels[num].values
        labels = self.labels[num]
        return unique_vals.take(labels)

    def format(self, space=2, sparsify=True, adjoin=True, names=False):
        if len(self) == 0:
            return []

        stringified_levels = [lev.format() for lev in self.levels]

        result_levels = []
        for lab, lev, name in zip(self.labels, stringified_levels, self.names):
            level = []

            if names:
                level.append(str(name) if name is not None else '')

            level.extend(np.array(lev, dtype=object).take(lab))
            result_levels.append(level)

        if sparsify:
            result_levels = _sparsify(result_levels)

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

    @cache_readonly
    def lexsort_depth(self):
        if self.sortorder is not None:
            if self.sortorder == 0:
                return self.nlevels
            else:
                return 0

        for k in range(self.nlevels, 0, -1):
            if lib.is_lexsorted(self.labels[:k]):
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
        levels = []
        labels = []
        for arr in arrays:
            factor = Factor(arr)
            levels.append(factor.levels)
            labels.append(factor.labels)

        return MultiIndex(levels=levels, labels=labels, sortorder=sortorder,
                          names=names)

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
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            return tuple(lev[lab[key]]
                         for lev, lab in zip(self.levels, self.labels))
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)
                sortorder = self.sortorder
            else:
                # cannot be sure whether the result will be sorted
                sortorder = None

            new_tuples = arr_idx[key]
            new_labels = [lab[key] for lab in self.labels]

            # an optimization
            result = new_tuples.view(MultiIndex)
            result.levels = list(self.levels)
            result.labels = new_labels
            result.sortorder = sortorder
            result.names = self.names

            return result

    def take(self, *args, **kwargs):
        """
        Analogous to ndarray.take
        """
        new_labels = [lab.take(*args, **kwargs) for lab in self.labels]
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
        if isinstance(other, (list, tuple)):
            to_concat = (self.values,) + tuple(k.values for k in other)
        else:
            to_concat = self.values, other.values
        new_tuples = np.concatenate(to_concat)
        return MultiIndex.from_tuples(new_tuples, names=self.names)

    def argsort(self, *args, **kwargs):
        return self.values.argsort()

    def drop(self, labels):
        """
        Make new MultiIndex with passed list of labels deleted

        Parameters
        ----------
        labels : array-like
            Must be a list of tuples

        Returns
        -------
        dropped : MultiIndex
        """
        try:
            if not isinstance(labels, np.ndarray):
                labels = com._asarray_tuplesafe(labels)
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
        assert(len(order) == self.nlevels)
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
        from pandas.core.frame import _indexer_from_factorized

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

        new_labels = [lab.take(indexer) for lab in self.labels]

        new_index = MultiIndex._from_elements(self.values.take(indexer),
                                              labels=new_labels,
                                              levels=self.levels,
                                              names=self.names,
                                              sortorder=level)

        return new_index, indexer

    def get_indexer(self, target, method=None):
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
        if isinstance(target, MultiIndex) and target._is_legacy_format:
            target_index = target.get_tuple_index()

        if target_index.dtype != object:
            return np.ones(len(target_index)) * -1

        self_index = self
        if self._is_legacy_format:
            self_index = self.get_tuple_index()

        if method == 'pad':
            indexer = self._pad(self_index, target_index, self_index.indexMap,
                                target.indexMap)
        elif method == 'backfill':
            indexer = self._backfill(self_index, target_index,
                                     self_index.indexMap, target.indexMap)
        else:
            indexer = self._merge_indexer(target_index, self_index.indexMap)

        return indexer

    def reindex(self, target, method=None, level=None):
        """
        Performs any necessary conversion on the input index and calls
        get_indexer. This method is here so MultiIndex and an Index of
        like-labeled tuples can play nice together

        Returns
        -------
        (new_index, indexer, mask) : (MultiIndex, ndarray, ndarray)
        """
        if level is not None:
            target, _, indexer = self._join_level(target, level, how='left',
                                                  return_indexers=True)
        else:
            if self.equals(target):
                indexer = None
            else:
                indexer = self.get_indexer(target, method=method)

        if not isinstance(target, MultiIndex):
            if indexer is None:
                target = self
            elif (indexer >= 0).all():
                target = self.take(indexer)
            else:
                # hopefully?
                target = MultiIndex.from_tuples(target)

        return target, indexer

    def get_tuple_index(self):
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
        (begin, end) : (int, int)

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
                if lib.infer_dtype([lab]) != lev.inferred_type:
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
                return self._engine.get_loc(key)
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
            assert(len(key) == len(level))
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
                    return self._engine.get_loc(key), None
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
            svalues = np.asarray(self.levels[i]).take(self.labels[i])
            ovalues = np.asarray(other.levels[i]).take(other.labels[i])
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
        if not isinstance(other, MultiIndex):
            raise TypeError('can only call with other hierarchical '
                            'index objects')

        assert(self.nlevels == other.nlevels)

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
        if not isinstance(item, tuple) or len(item) != self.nlevels:
            raise Exception("%s cannot be inserted in this MultIndex"
                            % str(item))

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
        return MultiIndex(levels=self.levels, labels=new_labels)

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

NULL_INDEX = Index([])


def _sparsify(label_list):
    pivoted = zip(*label_list)
    k = len(label_list)

    result = [pivoted[0]]
    prev = pivoted[0]

    for cur in pivoted[1:]:
        sparse_cur = []

        for i, (p, t) in enumerate(zip(prev, cur)):
            if i == k - 1:
                sparse_cur.append(t)
                result.append(sparse_cur)
                break

            if p == t:
                sparse_cur.append('')
            else:
                sparse_cur.append(t)

        prev = cur

    return zip(*result)


def _ensure_index(index_like):
    if isinstance(index_like, Index):
        return index_like
    return Index(index_like)


def _validate_join_method(method):
    if method not in ['left', 'right', 'inner', 'outer']:
        raise Exception('do not recognize join method %s' % method)


# TODO: handle index names!
def _get_combined_index(indexes, intersect=False):
    indexes = _get_distinct_indexes(indexes)
    if len(indexes) == 0:
        return NULL_INDEX
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
    assert(len(indexes) > 0)
    if len(indexes) == 1:
        result = indexes[0]
        if isinstance(result, list):
            result = Index(sorted(result))
        return result

    indexes, kind = _sanitize_and_check(indexes)

    if kind == 'special':
        result = indexes[0]
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
