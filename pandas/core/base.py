"""
Base and utility classes for pandas objects.
"""
from pandas import compat
import numpy as np
from pandas.core import common as com

class StringMixin(object):

    """implements string methods so long as object defines a `__unicode__`
    method.

    Handles Python2/3 compatibility transparently.
    """
    # side note - this could be made into a metaclass if more than one
    #             object needs

    #----------------------------------------------------------------------
    # Formatting

    def __unicode__(self):
        raise NotImplementedError

    def __str__(self):
        """
        Return a string representation for a particular Object

        Invoked by str(df) in both py2/py3.
        Yields Bytestring in Py2, Unicode String in py3.
        """

        if compat.PY3:
            return self.__unicode__()
        return self.__bytes__()

    def __bytes__(self):
        """
        Return a string representation for a particular object.

        Invoked by bytes(obj) in py3 only.
        Yields a bytestring in both py2/py3.
        """
        from pandas.core.config import get_option

        encoding = get_option("display.encoding")
        return self.__unicode__().encode(encoding, 'replace')

    def __repr__(self):
        """
        Return a string representation for a particular object.

        Yields Bytestring in Py2, Unicode String in py3.
        """
        return str(self)


class PandasObject(StringMixin):

    """baseclass for various pandas objects"""

    @property
    def _constructor(self):
        """class constructor (for this class it's just `__class__`"""
        return self.__class__

    def __unicode__(self):
        """
        Return a string representation for a particular object.

        Invoked by unicode(obj) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        # Should be overwritten by base classes
        return object.__repr__(self)

    def _local_dir(self):
        """ provide addtional __dir__ for this object """
        return []

    def __dir__(self):
        """
        Provide method name lookup and completion
        Only provide 'public' methods
        """
        return list(sorted(list(set(dir(type(self)) + self._local_dir()))))

    def _reset_cache(self, key=None):
        """
        Reset cached properties. If ``key`` is passed, only clears that key.
        """
        if getattr(self, '_cache', None) is None:
            return
        if key is None:
            self._cache.clear()
        else:
            self._cache.pop(key, None)


class FrozenList(PandasObject, list):

    """
    Container that doesn't allow setting item *but*
    because it's technically non-hashable, will be used
    for lookups, appropriately, etc.
    """
    # Sidenote: This has to be of type list, otherwise it messes up PyTables
    #           typechecks

    def __add__(self, other):
        if isinstance(other, tuple):
            other = list(other)
        return self.__class__(super(FrozenList, self).__add__(other))

    __iadd__ = __add__

    # Python 2 compat
    def __getslice__(self, i, j):
        return self.__class__(super(FrozenList, self).__getslice__(i, j))

    def __getitem__(self, n):
        # Python 3 compat
        if isinstance(n, slice):
            return self.__class__(super(FrozenList, self).__getitem__(n))
        return super(FrozenList, self).__getitem__(n)

    def __radd__(self, other):
        if isinstance(other, tuple):
            other = list(other)
        return self.__class__(other + list(self))

    def __eq__(self, other):
        if isinstance(other, (tuple, FrozenList)):
            other = list(other)
        return super(FrozenList, self).__eq__(other)

    __req__ = __eq__

    def __mul__(self, other):
        return self.__class__(super(FrozenList, self).__mul__(other))

    __imul__ = __mul__

    def __reduce__(self):
        return self.__class__, (list(self),)

    def __hash__(self):
        return hash(tuple(self))

    def _disabled(self, *args, **kwargs):
        """This method will not function because object is immutable."""
        raise TypeError("'%s' does not support mutable operations." %
                        self.__class__.__name__)

    def __unicode__(self):
        from pandas.core.common import pprint_thing
        return pprint_thing(self, quote_strings=True,
                            escape_chars=('\t', '\r', '\n'))

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__,
                           str(self))

    __setitem__ = __setslice__ = __delitem__ = __delslice__ = _disabled
    pop = append = extend = remove = sort = insert = _disabled


class FrozenNDArray(PandasObject, np.ndarray):

    # no __array_finalize__ for now because no metadata
    def __new__(cls, data, dtype=None, copy=False):
        if copy is None:
            copy = not isinstance(data, FrozenNDArray)
        res = np.array(data, dtype=dtype, copy=copy).view(cls)
        return res

    def _disabled(self, *args, **kwargs):
        """This method will not function because object is immutable."""
        raise TypeError("'%s' does not support mutable operations." %
                        self.__class__)

    __setitem__ = __setslice__ = __delitem__ = __delslice__ = _disabled
    put = itemset = fill = _disabled

    def _shallow_copy(self):
        return self.view()

    def values(self):
        """returns *copy* of underlying array"""
        arr = self.view(np.ndarray).copy()
        return arr

    def __unicode__(self):
        """
        Return a string representation for this object.

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        prepr = com.pprint_thing(self, escape_chars=('\t', '\r', '\n'),
                                 quote_strings=True)
        return "%s(%s, dtype='%s')" % (type(self).__name__, prepr, self.dtype)


class IndexOpsMixin(object):
    """ common ops mixin to support a unified inteface / docs for Series / Index """

    def _is_allowed_index_op(self, name):
        if not self._allow_index_ops:
            raise TypeError("cannot perform an {name} operations on this type {typ}".format(
                name=name,typ=type(self._get_access_object())))

    def _ops_compat(self, name, op_accessor):

        obj = self._get_access_object()
        try:
            return self._wrap_access_object(getattr(obj,op_accessor))
        except AttributeError:
            raise TypeError("cannot perform an {name} operations on this type {typ}".format(
                name=name,typ=type(obj)))

    def _get_access_object(self):
        if isinstance(self, com.ABCSeries):
            return self.index
        return self

    def _wrap_access_object(self, obj):
        # we may need to coerce the input as we don't want non int64 if
        # we have an integer result
        if hasattr(obj,'dtype') and com.is_integer_dtype(obj):
            obj = obj.astype(np.int64)

        if isinstance(self, com.ABCSeries):
            return self._constructor(obj,index=self.index).__finalize__(self)

        return obj

    def max(self):
        """ The maximum value of the object """
        import pandas.core.nanops
        return pandas.core.nanops.nanmax(self.values)

    def min(self):
        """ The minimum value of the object """
        import pandas.core.nanops
        return pandas.core.nanops.nanmin(self.values)

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None):
        """
        Returns object containing counts of unique values. The resulting object
        will be in descending order so that the first element is the most
        frequently-occurring element. Excludes NA values.

        Parameters
        ----------
        normalize : boolean, default False
            If True then the object returned will contain the relative
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

    def unique(self):
        """
        Return array of unique values in the object. Significantly faster than
        numpy.unique. Includes NA values.

        Returns
        -------
        uniques : ndarray
        """
        from pandas.core.nanops import unique1d
        return unique1d(self.values)

    def nunique(self):
        """
        Return count of unique elements in the object. Excludes NA values.

        Returns
        -------
        nunique : int
        """
        return len(self.value_counts())

    def factorize(self, sort=False, na_sentinel=-1):
        """
        Encode the object as an enumerated type or categorical variable

        Parameters
        ----------
        sort : boolean, default False
            Sort by values
        na_sentinel: int, default -1
            Value to mark "not found"

        Returns
        -------
        labels : the indexer to the original array
        uniques : the unique Index
        """
        from pandas.core.algorithms import factorize
        return factorize(self, sort=sort, na_sentinel=na_sentinel)

# facilitate the properties on the wrapped ops
def _field_accessor(name, docstring=None):
    op_accessor = '_{0}'.format(name)
    def f(self):
        return self._ops_compat(name,op_accessor)

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)

class DatetimeIndexOpsMixin(object):
    """ common ops mixin to support a unified inteface datetimelike Index """

    def _is_allowed_datetime_index_op(self, name):
        if not self._allow_datetime_index_ops:
            raise TypeError("cannot perform an {name} operations on this type {typ}".format(
                name=name,typ=type(self._get_access_object())))

    def _is_allowed_period_index_op(self, name):
        if not self._allow_period_index_ops:
            raise TypeError("cannot perform an {name} operations on this type {typ}".format(
                name=name,typ=type(self._get_access_object())))

    def _ops_compat(self, name, op_accessor):

        from pandas.tseries.index import DatetimeIndex
        from pandas.tseries.period import PeriodIndex
        obj = self._get_access_object()
        if isinstance(obj, DatetimeIndex):
            self._is_allowed_datetime_index_op(name)
        elif isinstance(obj, PeriodIndex):
            self._is_allowed_period_index_op(name)
        try:
            return self._wrap_access_object(getattr(obj,op_accessor))
        except AttributeError:
            raise TypeError("cannot perform an {name} operations on this type {typ}".format(
                name=name,typ=type(obj)))

    date = _field_accessor('date','Returns numpy array of datetime.date. The date part of the Timestamps')
    time = _field_accessor('time','Returns numpy array of datetime.time. The time part of the Timestamps')
    year = _field_accessor('year', "The year of the datetime")
    month = _field_accessor('month', "The month as January=1, December=12")
    day = _field_accessor('day', "The days of the datetime")
    hour = _field_accessor('hour', "The hours of the datetime")
    minute = _field_accessor('minute', "The minutes of the datetime")
    second = _field_accessor('second', "The seconds of the datetime")
    microsecond = _field_accessor('microsecond', "The microseconds of the datetime")
    nanosecond = _field_accessor('nanosecond', "The nanoseconds of the datetime")
    weekofyear = _field_accessor('weekofyear', "The week ordinal of the year")
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', "The day of the week with Monday=0, Sunday=6")
    weekday = dayofweek
    dayofyear = _field_accessor('dayofyear', "The ordinal day of the year")
    quarter = _field_accessor('quarter', "The quarter of the date")
    qyear = _field_accessor('qyear')
    is_month_start = _field_accessor('is_month_start', "Logical indicating if first day of month (defined by frequency)")
    is_month_end = _field_accessor('is_month_end', "Logical indicating if last day of month (defined by frequency)")
    is_quarter_start = _field_accessor('is_quarter_start', "Logical indicating if first day of quarter (defined by frequency)")
    is_quarter_end = _field_accessor('is_quarter_end', "Logical indicating if last day of quarter (defined by frequency)")
    is_year_start = _field_accessor('is_year_start', "Logical indicating if first day of year (defined by frequency)")
    is_year_end = _field_accessor('is_year_end', "Logical indicating if last day of year (defined by frequency)")
