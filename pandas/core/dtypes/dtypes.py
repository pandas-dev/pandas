""" define extension dtypes """

import re
import numpy as np
from pandas import compat
from pandas.core.dtypes.generic import ABCIndexClass, ABCCategoricalIndex


class ExtensionDtype(object):
    """
    A np.dtype duck-typed class, suitable for holding a custom dtype.

    THIS IS NOT A REAL NUMPY DTYPE
    """
    name = None
    names = None
    type = None
    subdtype = None
    kind = None
    str = None
    num = 100
    shape = tuple()
    itemsize = 8
    base = None
    isbuiltin = 0
    isnative = 0
    _metadata = []
    _cache = {}

    def __unicode__(self):
        return self.name

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

    def __hash__(self):
        raise NotImplementedError("sub-classes should implement an __hash__ "
                                  "method")

    def __eq__(self, other):
        raise NotImplementedError("sub-classes should implement an __eq__ "
                                  "method")

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getstate__(self):
        # pickle support; we don't want to pickle the cache
        return {k: getattr(self, k, None) for k in self._metadata}

    @classmethod
    def reset_cache(cls):
        """ clear the cache """
        cls._cache = {}

    @classmethod
    def is_dtype(cls, dtype):
        """ Return a boolean if the passed type is an actual dtype that
        we can match (via string or type)
        """
        if hasattr(dtype, 'dtype'):
            dtype = dtype.dtype
        if isinstance(dtype, np.dtype):
            return False
        elif dtype is None:
            return False
        elif isinstance(dtype, cls):
            return True
        try:
            return cls.construct_from_string(dtype) is not None
        except:
            return False


class CategoricalDtypeType(type):
    """
    the type of CategoricalDtype, this metaclass determines subclass ability
    """
    pass


class CategoricalDtype(ExtensionDtype):
    """
    Type for categorical data with the categories and orderedness

    .. versionchanged:: 0.21.0

    Parameters
    ----------
    categories : sequence, optional
        Must be unique, and must not contain any nulls.
    ordered : bool, default False

    Notes
    -----
    This class is useful for specifying the type of a ``Categorical``
    independent of the values. See :ref:`categorical.categoricaldtype`
    for more.

    Examples
    --------
    >>> t = CategoricalDtype(categories=['b', 'a'], ordered=True)
    >>> pd.Series(['a', 'b', 'a', 'c'], dtype=t)
    0      a
    1      b
    2      a
    3    NaN
    dtype: category
    Categories (2, object): [b < a]

    See Also
    --------
    pandas.Categorical
    """
    # TODO: Document public vs. private API
    name = 'category'
    type = CategoricalDtypeType
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')
    _metadata = ['categories', 'ordered']
    _cache = {}

    def __init__(self, categories=None, ordered=False):
        self._finalize(categories, ordered, fastpath=False)

    @classmethod
    def _from_fastpath(cls, categories=None, ordered=False):
        self = cls.__new__(cls)
        self._finalize(categories, ordered, fastpath=True)
        return self

    @classmethod
    def _from_categorical_dtype(cls, dtype, categories=None, ordered=None):
        if categories is ordered is None:
            return dtype
        if categories is None:
            categories = dtype.categories
        if ordered is None:
            ordered = dtype.ordered
        return cls(categories, ordered)

    def _finalize(self, categories, ordered, fastpath=False):

        if ordered is None:
            ordered = False
        else:
            self._validate_ordered(ordered)

        if categories is not None:
            categories = self._validate_categories(categories,
                                                   fastpath=fastpath)

        self._categories = categories
        self._ordered = ordered

    def __setstate__(self, state):
        self._categories = state.pop('categories', None)
        self._ordered = state.pop('ordered', False)

    def __hash__(self):
        # _hash_categories returns a uint64, so use the negative
        # space for when we have unknown categories to avoid a conflict
        if self.categories is None:
            if self.ordered:
                return -1
            else:
                return -2
        # We *do* want to include the real self.ordered here
        return int(self._hash_categories(self.categories, self.ordered))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name

        if not (hasattr(other, 'ordered') and hasattr(other, 'categories')):
            return False
        elif self.categories is None or other.categories is None:
            # We're forced into a suboptimal corner thanks to math and
            # backwards compatibility. We require that `CDT(...) == 'category'`
            # for all CDTs **including** `CDT(None, ...)`. Therefore, *all*
            # CDT(., .) = CDT(None, False) and *all*
            # CDT(., .) = CDT(None, True).
            return True
        elif self.ordered:
            return other.ordered and self.categories.equals(other.categories)
        elif other.ordered:
            return False
        else:
            # both unordered; this could probably be optimized / cached
            return hash(self) == hash(other)

    def __repr__(self):
        tpl = u'CategoricalDtype(categories={}ordered={})'
        if self.categories is None:
            data = u"None, "
        else:
            data = self.categories._format_data(name=self.__class__.__name__)
        return tpl.format(data, self.ordered)

    @staticmethod
    def _hash_categories(categories, ordered=True):
        from pandas.core.util.hashing import (
            hash_array, _combine_hash_arrays, hash_tuples
        )

        if len(categories) and isinstance(categories[0], tuple):
            # assumes if any individual category is a tuple, then all our. ATM
            # I don't really want to support just some of the categories being
            # tuples.
            categories = list(categories)  # breaks if a np.array of categories
            cat_array = hash_tuples(categories)
        else:
            if categories.dtype == 'O':
                types = [type(x) for x in categories]
                if not len(set(types)) == 1:
                    # TODO: hash_array doesn't handle mixed types. It casts
                    # everything to a str first, which means we treat
                    # {'1', '2'} the same as {'1', 2}
                    # find a better solution
                    cat_array = np.array([hash(x) for x in categories])
                    hashed = hash((tuple(categories), ordered))
                    return hashed
            cat_array = hash_array(np.asarray(categories), categorize=False)
        if ordered:
            cat_array = np.vstack([
                cat_array, np.arange(len(cat_array), dtype=cat_array.dtype)
            ])
        else:
            cat_array = [cat_array]
        hashed = _combine_hash_arrays(iter(cat_array),
                                      num_items=len(cat_array))
        if len(hashed) == 0:
            # bug in Numpy<1.12 for length 0 arrays. Just return the correct
            # value of 0
            return 0
        else:
            return np.bitwise_xor.reduce(hashed)

    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if
        it's not possible """
        try:
            if string == 'category':
                return cls()
        except:
            pass

        raise TypeError("cannot construct a CategoricalDtype")

    @staticmethod
    def _validate_ordered(ordered):
        """
        Validates that we have a valid ordered parameter. If
        it is not a boolean, a TypeError will be raised.

        Parameters
        ----------
        ordered : object
            The parameter to be verified.

        Raises
        ------
        TypeError
            If 'ordered' is not a boolean.
        """
        from pandas.core.dtypes.common import is_bool
        if not is_bool(ordered):
            raise TypeError("'ordered' must either be 'True' or 'False'")

    @staticmethod
    def _validate_categories(categories, fastpath=False):
        """
        Validates that we have good categories

        Parameters
        ----------
        categories : array-like
        fastpath : bool
            Whether to skip nan and uniqueness checks

        Returns
        -------
        categories : Index
        """
        from pandas import Index

        if not isinstance(categories, ABCIndexClass):
            categories = Index(categories, tupleize_cols=False)

        if not fastpath:

            if categories.hasnans:
                raise ValueError('Categorial categories cannot be null')

            if not categories.is_unique:
                raise ValueError('Categorical categories must be unique')

        if isinstance(categories, ABCCategoricalIndex):
            categories = categories.categories

        return categories

    @property
    def categories(self):
        """
        An ``Index`` containing the unique categories allowed.
        """
        return self._categories

    @property
    def ordered(self):
        """Whether the categories have an ordered relationship"""
        return self._ordered


class DatetimeTZDtypeType(type):
    """
    the type of DatetimeTZDtype, this metaclass determines subclass ability
    """
    pass


class DatetimeTZDtype(ExtensionDtype):

    """
    A np.dtype duck-typed class, suitable for holding a custom datetime with tz
    dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of
    np.datetime64[ns]
    """
    type = DatetimeTZDtypeType
    kind = 'M'
    str = '|M8[ns]'
    num = 101
    base = np.dtype('M8[ns]')
    _metadata = ['unit', 'tz']
    _match = re.compile("(datetime64|M8)\[(?P<unit>.+), (?P<tz>.+)\]")
    _cache = {}

    def __new__(cls, unit=None, tz=None):
        """ Create a new unit if needed, otherwise return from the cache

        Parameters
        ----------
        unit : string unit that this represents, currently must be 'ns'
        tz : string tz that this represents
        """

        if isinstance(unit, DatetimeTZDtype):
            unit, tz = unit.unit, unit.tz

        elif unit is None:
            # we are called as an empty constructor
            # generally for pickle compat
            return object.__new__(cls)

        elif tz is None:

            # we were passed a string that we can construct
            try:
                m = cls._match.search(unit)
                if m is not None:
                    unit = m.groupdict()['unit']
                    tz = m.groupdict()['tz']
            except:
                raise ValueError("could not construct DatetimeTZDtype")

        elif isinstance(unit, compat.string_types):

            if unit != 'ns':
                raise ValueError("DatetimeTZDtype only supports ns units")

            unit = unit
            tz = tz

        if tz is None:
            raise ValueError("DatetimeTZDtype constructor must have a tz "
                             "supplied")

        # hash with the actual tz if we can
        # some cannot be hashed, so stringfy
        try:
            key = (unit, tz)
            hash(key)
        except TypeError:
            key = (unit, str(tz))

        # set/retrieve from cache
        try:
            return cls._cache[key]
        except KeyError:
            u = object.__new__(cls)
            u.unit = unit
            u.tz = tz
            cls._cache[key] = u
            return u

    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if
        it's not possible
        """
        try:
            return cls(unit=string)
        except ValueError:
            raise TypeError("could not construct DatetimeTZDtype")

    def __unicode__(self):
        # format the tz
        return "datetime64[{unit}, {tz}]".format(unit=self.unit, tz=self.tz)

    @property
    def name(self):
        return str(self)

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name

        return (isinstance(other, DatetimeTZDtype) and
                self.unit == other.unit and
                str(self.tz) == str(other.tz))


class PeriodDtypeType(type):
    """
    the type of PeriodDtype, this metaclass determines subclass ability
    """
    pass


class PeriodDtype(ExtensionDtype):
    __metaclass__ = PeriodDtypeType
    """
    A Period duck-typed class, suitable for holding a period with freq dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of np.int64.
    """
    type = PeriodDtypeType
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')
    num = 102
    _metadata = ['freq']
    _match = re.compile("(P|p)eriod\[(?P<freq>.+)\]")
    _cache = {}

    def __new__(cls, freq=None):
        """
        Parameters
        ----------
        freq : frequency
        """

        if isinstance(freq, PeriodDtype):
            return freq

        elif freq is None:
            # empty constructor for pickle compat
            return object.__new__(cls)

        from pandas.tseries.offsets import DateOffset
        if not isinstance(freq, DateOffset):
            freq = cls._parse_dtype_strict(freq)

        try:
            return cls._cache[freq.freqstr]
        except KeyError:
            u = object.__new__(cls)
            u.freq = freq
            cls._cache[freq.freqstr] = u
            return u

    @classmethod
    def _parse_dtype_strict(cls, freq):
        if isinstance(freq, compat.string_types):
            if freq.startswith('period[') or freq.startswith('Period['):
                m = cls._match.search(freq)
                if m is not None:
                    freq = m.group('freq')
            from pandas.tseries.frequencies import to_offset
            freq = to_offset(freq)
            if freq is not None:
                return freq

        raise ValueError("could not construct PeriodDtype")

    @classmethod
    def construct_from_string(cls, string):
        """
        attempt to construct this type from a string, raise a TypeError
        if its not possible
        """
        from pandas.tseries.offsets import DateOffset
        if isinstance(string, (compat.string_types, DateOffset)):
            # avoid tuple to be regarded as freq
            try:
                return cls(freq=string)
            except ValueError:
                pass
        raise TypeError("could not construct PeriodDtype")

    def __unicode__(self):
        return "period[{freq}]".format(freq=self.freq.freqstr)

    @property
    def name(self):
        return str(self)

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name or other == self.name.title()

        return isinstance(other, PeriodDtype) and self.freq == other.freq

    @classmethod
    def is_dtype(cls, dtype):
        """
        Return a boolean if we if the passed type is an actual dtype that we
        can match (via string or type)
        """

        if isinstance(dtype, compat.string_types):
            # PeriodDtype can be instanciated from freq string like "U",
            # but dosn't regard freq str like "U" as dtype.
            if dtype.startswith('period[') or dtype.startswith('Period['):
                try:
                    if cls._parse_dtype_strict(dtype) is not None:
                        return True
                    else:
                        return False
                except ValueError:
                    return False
            else:
                return False
        return super(PeriodDtype, cls).is_dtype(dtype)


class IntervalDtypeType(type):
    """
    the type of IntervalDtype, this metaclass determines subclass ability
    """
    pass


class IntervalDtype(ExtensionDtype):
    __metaclass__ = IntervalDtypeType
    """
    A Interval duck-typed class, suitable for holding an interval

    THIS IS NOT A REAL NUMPY DTYPE
    """
    type = IntervalDtypeType
    kind = None
    str = '|O08'
    base = np.dtype('O')
    num = 103
    _metadata = ['subtype']
    _match = re.compile("(I|i)nterval\[(?P<subtype>.+)\]")
    _cache = {}

    def __new__(cls, subtype=None):
        """
        Parameters
        ----------
        subtype : the dtype of the Interval
        """

        if isinstance(subtype, IntervalDtype):
            return subtype
        elif subtype is None:
            # we are called as an empty constructor
            # generally for pickle compat
            u = object.__new__(cls)
            u.subtype = None
            return u
        elif (isinstance(subtype, compat.string_types) and
              subtype == 'interval'):
            subtype = ''
        else:
            if isinstance(subtype, compat.string_types):
                m = cls._match.search(subtype)
                if m is not None:
                    subtype = m.group('subtype')

            from pandas.core.dtypes.common import pandas_dtype
            try:
                subtype = pandas_dtype(subtype)
            except TypeError:
                raise ValueError("could not construct IntervalDtype")

        if subtype is None:
            u = object.__new__(cls)
            u.subtype = None
            return u

        try:
            return cls._cache[str(subtype)]
        except KeyError:
            u = object.__new__(cls)
            u.subtype = subtype
            cls._cache[str(subtype)] = u
            return u

    @classmethod
    def construct_from_string(cls, string):
        """
        attempt to construct this type from a string, raise a TypeError
        if its not possible
        """
        if isinstance(string, compat.string_types):
            try:
                return cls(string)
            except ValueError:
                pass
        raise TypeError("could not construct IntervalDtype")

    def __unicode__(self):
        if self.subtype is None:
            return "interval"
        return "interval[{subtype}]".format(subtype=self.subtype)

    @property
    def name(self):
        return str(self)

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name or other == self.name.title()

        return (isinstance(other, IntervalDtype) and
                self.subtype == other.subtype)

    @classmethod
    def is_dtype(cls, dtype):
        """
        Return a boolean if we if the passed type is an actual dtype that we
        can match (via string or type)
        """

        if isinstance(dtype, compat.string_types):
            if dtype.lower().startswith('interval'):
                try:
                    if cls.construct_from_string(dtype) is not None:
                        return True
                    else:
                        return False
                except ValueError:
                    return False
            else:
                return False
        return super(IntervalDtype, cls).is_dtype(dtype)
