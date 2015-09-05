"""
Base and utility classes for pandas objects.
"""
from pandas import compat
import numpy as np
from pandas.core import common as com
import pandas.core.nanops as nanops
import pandas.lib as lib
from pandas.util.decorators import Appender, cache_readonly, deprecate_kwarg
from pandas.core.strings import StringMethods
from pandas.core.common import AbstractMethodError

_shared_docs = dict()
_indexops_doc_kwargs = dict(klass='IndexOpsMixin', inplace='',
                            duplicated='IndexOpsMixin')


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
        raise AbstractMethodError(self)

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

    def _dir_additions(self):
        """ add addtional __dir__ for this object """
        return set()

    def _dir_deletions(self):
        """ delete unwanted __dir__ for this object """
        return set()

    def __dir__(self):
        """
        Provide method name lookup and completion
        Only provide 'public' methods
        """
        rv = set(dir(type(self)))
        rv = (rv - self._dir_deletions()) | self._dir_additions()
        return sorted(rv)

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

class PandasDelegate(PandasObject):
    """ an abstract base class for delegating methods/properties """

    def _delegate_property_get(self, name, *args, **kwargs):
        raise TypeError("You cannot access the property {name}".format(name=name))

    def _delegate_property_set(self, name, value, *args, **kwargs):
        raise TypeError("The property {name} cannot be set".format(name=name))

    def _delegate_method(self, name, *args, **kwargs):
        raise TypeError("You cannot call method {name}".format(name=name))

    @classmethod
    def _add_delegate_accessors(cls, delegate, accessors, typ, overwrite=False):
        """
        add accessors to cls from the delegate class

        Parameters
        ----------
        cls : the class to add the methods/properties to
        delegate : the class to get methods/properties & doc-strings
        acccessors : string list of accessors to add
        typ : 'property' or 'method'
        overwrite : boolean, default False
           overwrite the method/property in the target class if it exists

        """

        def _create_delegator_property(name):

            def _getter(self):
                return self._delegate_property_get(name)
            def _setter(self, new_values):
                return self._delegate_property_set(name, new_values)

            _getter.__name__ = name
            _setter.__name__ = name

            return property(fget=_getter, fset=_setter, doc=getattr(delegate,name).__doc__)

        def _create_delegator_method(name):

            def f(self, *args, **kwargs):
                return self._delegate_method(name, *args, **kwargs)

            f.__name__ = name
            f.__doc__ = getattr(delegate,name).__doc__

            return f

        for name in accessors:

            if typ == 'property':
                f = _create_delegator_property(name)
            else:
                f = _create_delegator_method(name)

            # don't overwrite existing methods/properties
            if overwrite or not hasattr(cls, name):
                setattr(cls,name,f)


class AccessorProperty(object):
    """Descriptor for implementing accessor properties like Series.str
    """
    def __init__(self, accessor_cls, construct_accessor):
        self.accessor_cls = accessor_cls
        self.construct_accessor = construct_accessor
        self.__doc__ = accessor_cls.__doc__

    def __get__(self, instance, owner=None):
        if instance is None:
            # this ensures that Series.str.<method> is well defined
            return self.accessor_cls
        return self.construct_accessor(instance)

    def __set__(self, instance, value):
        raise AttributeError("can't set attribute")

    def __delete__(self, instance):
        raise AttributeError("can't delete attribute")


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

    # ndarray compatibility
    __array_priority__ = 1000

    def transpose(self):
        """ return the transpose, which is by definition self """
        return self

    T = property(transpose, doc="return the transpose, which is by definition self")

    @property
    def shape(self):
        """ return a tuple of the shape of the underlying data """
        return self.values.shape

    @property
    def ndim(self):
        """ return the number of dimensions of the underlying data, by definition 1 """
        return 1

    def item(self):
        """ return the first element of the underlying data as a python scalar """
        try:
            return self.values.item()
        except IndexError:
            # copy numpy's message here because Py26 raises an IndexError
            raise ValueError('can only convert an array of size 1 to a '
                             'Python scalar')

    @property
    def data(self):
        """ return the data pointer of the underlying data """
        return self.values.data

    @property
    def itemsize(self):
        """ return the size of the dtype of the item of the underlying data """
        return self.values.itemsize

    @property
    def nbytes(self):
        """ return the number of bytes in the underlying data """
        return self.values.nbytes

    @property
    def strides(self):
        """ return the strides of the underlying data """
        return self.values.strides

    @property
    def size(self):
        """ return the number of elements in the underlying data """
        return self.values.size

    @property
    def flags(self):
        """ return the ndarray.flags for the underlying data """
        return self.values.flags

    @property
    def base(self):
        """ return the base object if the memory of the underlying data is shared """
        return self.values.base

    @property
    def _values(self):
        """ the internal implementation """
        return self.values

    def max(self):
        """ The maximum value of the object """
        return nanops.nanmax(self.values)

    def argmax(self, axis=None):
        """
        return a ndarray of the maximum argument indexer

        See also
        --------
        numpy.ndarray.argmax
        """
        return nanops.nanargmax(self.values)

    def min(self):
        """ The minimum value of the object """
        return nanops.nanmin(self.values)

    def argmin(self, axis=None):
        """
        return a ndarray of the minimum argument indexer

        See also
        --------
        numpy.ndarray.argmin
        """
        return nanops.nanargmin(self.values)

    @cache_readonly
    def hasnans(self):
        """ return if I have any nans; enables various perf speedups """
        return com.isnull(self).any()

    def _reduce(self, op, name, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        """ perform the reduction type operation if we can """
        func = getattr(self,name,None)
        if func is None:
            raise TypeError("{klass} cannot perform the operation {op}".format(klass=self.__class__.__name__,op=name))
        return func(**kwds)

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None, dropna=True):
        """
        Returns object containing counts of unique values.

        The resulting object will be in descending order so that the
        first element is the most frequently-occurring element.
        Excludes NA values by default.

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
        dropna : boolean, default True
            Don't include counts of NaN.

        Returns
        -------
        counts : Series
        """
        from pandas.core.algorithms import value_counts
        from pandas.tseries.api import DatetimeIndex, PeriodIndex
        result = value_counts(self, sort=sort, ascending=ascending,
                              normalize=normalize, bins=bins, dropna=dropna)

        if isinstance(self, PeriodIndex):
            # preserve freq
            result.index = self._simple_new(result.index.values,
                                            freq=self.freq)
        elif isinstance(self, DatetimeIndex):
            result.index = self._simple_new(result.index.values,
                                            tz=getattr(self, 'tz', None))
        return result

    def unique(self):
        """
        Return array of unique values in the object. Significantly faster than
        numpy.unique. Includes NA values.

        Returns
        -------
        uniques : ndarray
        """
        from pandas.core.nanops import unique1d
        values = self.values
        if hasattr(values,'unique'):
            return values.unique()

        return unique1d(values)

    def nunique(self, dropna=True):
        """
        Return number of unique elements in the object.

        Excludes NA values by default.

        Parameters
        ----------
        dropna : boolean, default True
            Don't include NaN in the count.

        Returns
        -------
        nunique : int
        """
        uniqs = self.unique()
        n = len(uniqs)
        if dropna and com.isnull(uniqs).any():
            n -= 1
        return n


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

    def searchsorted(self, key, side='left'):
        """ np.ndarray searchsorted compat """

        ### FIXME in GH7447
        #### needs coercion on the key (DatetimeIndex does alreay)
        #### needs tests/doc-string
        return self.values.searchsorted(key, side=side)

    # string methods
    def _make_str_accessor(self):
        from pandas.core.series import Series
        from pandas.core.index import Index
        if isinstance(self, Series) and not com.is_object_dtype(self.dtype):
            # this really should exclude all series with any non-string values,
            # but that isn't practical for performance reasons until we have a
            # str dtype (GH 9343)
            raise AttributeError("Can only use .str accessor with string "
                                 "values, which use np.object_ dtype in "
                                 "pandas")
        elif isinstance(self, Index):
            # see scc/inferrence.pyx which can contain string values
            allowed_types = ('string', 'unicode', 'mixed', 'mixed-integer')
            if self.inferred_type not in allowed_types:
                message = ("Can only use .str accessor with string values "
                           "(i.e. inferred_type is 'string', 'unicode' or 'mixed')")
                raise AttributeError(message)
            if self.nlevels > 1:
                message = "Can only use .str accessor with Index, not MultiIndex"
                raise AttributeError(message)
        return StringMethods(self)

    str = AccessorProperty(StringMethods, _make_str_accessor)

    def _dir_additions(self):
        return set()

    def _dir_deletions(self):
        try:
            getattr(self, 'str')
        except AttributeError:
            return set(['str'])
        return set()

    _shared_docs['drop_duplicates'] = (
        """Return %(klass)s with duplicate values removed

        Parameters
        ----------

        keep : {'first', 'last', False}, default 'first'
            - ``first`` : Drop duplicates except for the first occurrence.
            - ``last`` : Drop duplicates except for the last occurrence.
            - False : Drop all duplicates.
        take_last : deprecated
        %(inplace)s

        Returns
        -------
        deduplicated : %(klass)s
        """)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(_shared_docs['drop_duplicates'] % _indexops_doc_kwargs)
    def drop_duplicates(self, keep='first', inplace=False):
        duplicated = self.duplicated(keep=keep)
        result = self[np.logical_not(duplicated)]
        if inplace:
            return self._update_inplace(result)
        else:
            return result

    _shared_docs['duplicated'] = (
        """Return boolean %(duplicated)s denoting duplicate values

        Parameters
        ----------
        keep : {'first', 'last', False}, default 'first'
            - ``first`` : Mark duplicates as ``True`` except for the first occurrence.
            - ``last`` : Mark duplicates as ``True`` except for the last occurrence.
            - False : Mark all duplicates as ``True``.
        take_last : deprecated

        Returns
        -------
        duplicated : %(duplicated)s
        """)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(_shared_docs['duplicated'] % _indexops_doc_kwargs)
    def duplicated(self, keep='first'):
        keys = com._values_from_object(com._ensure_object(self.values))
        duplicated = lib.duplicated(keys, keep=keep)
        try:
            return self._constructor(duplicated,
                                     index=self.index).__finalize__(self)
        except AttributeError:
            return np.array(duplicated, dtype=bool)

    #----------------------------------------------------------------------
    # abstracts

    def _update_inplace(self, result, **kwargs):
        raise AbstractMethodError(self)
