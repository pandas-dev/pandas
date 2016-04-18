import numpy as np
import pandas.lib as lib
import pandas.algos as _algos
import pandas.index as _index

from pandas import compat
from pandas.indexes.base import Index, InvalidIndexError
from pandas.util.decorators import Appender, cache_readonly
import pandas.core.common as com
from pandas.core.common import (is_dtype_equal, isnull, pandas_dtype,
                                is_float_dtype, is_object_dtype,
                                is_integer_dtype)
import pandas.indexes.base as ibase


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
        kind : {'ix', 'loc', 'getitem'}

        Returns
        -------
        label :  object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """
        assert kind in ['ix', 'loc', 'getitem', None]

        # we will try to coerce to integers
        return self._maybe_cast_indexer(label)

    def _convert_tolerance(self, tolerance):
        try:
            return float(tolerance)
        except ValueError:
            raise ValueError('tolerance argument for %s must be numeric: %r' %
                             (type(self).__name__, tolerance))


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

    def __new__(cls, data=None, dtype=None, copy=False, name=None,
                fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name=name)

        # isscalar, generators handled in coerce_to_ndarray
        data = cls._coerce_to_ndarray(data)

        if issubclass(data.dtype.type, compat.string_types):
            cls._string_data_error(data)

        elif issubclass(data.dtype.type, np.integer):
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

    def _convert_scalar_indexer(self, key, kind=None):
        """
        convert a scalar indexer

        Parameters
        ----------
        key : label of the slice bound
        kind : {'ix', 'loc', 'getitem'} or None
        """

        assert kind in ['ix', 'loc', 'getitem', 'iloc', None]

        # don't coerce ilocs to integers
        if kind != 'iloc':
            key = self._maybe_cast_indexer(key)
        return (super(Int64Index, self)
                ._convert_scalar_indexer(key, kind=kind))

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        try:
            return com.array_equivalent(com._values_from_object(self),
                                        com._values_from_object(other))
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

    def __new__(cls, data=None, dtype=None, copy=False, name=None,
                fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name)

        data = cls._coerce_to_ndarray(data)

        if issubclass(data.dtype.type, compat.string_types):
            cls._string_data_error(data)

        if dtype is None:
            dtype = np.float64
        dtype = np.dtype(dtype)

        # allow integer / object dtypes to be passed, but coerce to float64
        if dtype.kind in ['i', 'O', 'f']:
            dtype = np.float64

        else:
            raise TypeError("cannot support {0} dtype in "
                            "Float64Index".format(dtype))

        try:
            subarr = np.array(data, dtype=dtype, copy=copy)
        except:
            raise TypeError('Unsafe NumPy casting, you must explicitly cast')

        # coerce to float64 for storage
        if subarr.dtype != np.float64:
            subarr = subarr.astype(np.float64)

        return cls._simple_new(subarr, name)

    @property
    def inferred_type(self):
        return 'floating'

    def astype(self, dtype):
        dtype = pandas_dtype(dtype)
        if is_float_dtype(dtype) or is_integer_dtype(dtype):
            values = self._values.astype(dtype)
        elif is_object_dtype(dtype):
            values = self._values
        else:
            raise TypeError('Setting %s dtype to anything other than '
                            'float64 or object is not supported' %
                            self.__class__)
        return Index(values, name=self.name, dtype=dtype)

    def _convert_scalar_indexer(self, key, kind=None):
        """
        convert a scalar indexer

        Parameters
        ----------
        key : label of the slice bound
        kind : {'ix', 'loc', 'getitem'} or None
        """

        assert kind in ['ix', 'loc', 'getitem', 'iloc', None]

        if kind == 'iloc':
            return self._validate_indexer('positional', key, kind)

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
        return self.slice_indexer(key.start, key.stop, key.step, kind=kind)

    def _format_native_types(self, na_rep='', float_format=None, decimal='.',
                             quoting=None, **kwargs):
        from pandas.formats.format import FloatArrayFormatter
        formatter = FloatArrayFormatter(self.values, na_rep=na_rep,
                                        float_format=float_format,
                                        decimal=decimal, quoting=quoting,
                                        fixed_width=False)
        return formatter.get_result_as_array()

    def get_value(self, series, key):
        """ we always want to get an index value, never a value """
        if not lib.isscalar(key):
            raise InvalidIndexError

        from pandas.core.indexing import maybe_droplevels
        from pandas.core.series import Series

        k = com._values_from_object(key)
        loc = self.get_loc(k)
        new_values = com._values_from_object(series)[loc]

        if lib.isscalar(new_values) or new_values is None:
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
            if (not is_dtype_equal(self.dtype, other.dtype) or
                    self.shape != other.shape):
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
                return len(other) <= 1 and ibase._try_get_item(other) in self
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
