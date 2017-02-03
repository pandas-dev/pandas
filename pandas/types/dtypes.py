""" define extension dtypes """

import re
import reprlib
import numpy as np
from pandas import compat


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

    @classmethod
    def is_dtype(cls, dtype):
        """ Return a boolean if we if the passed type is an actual dtype that
        we can match (via string or type)
        """
        if hasattr(dtype, 'dtype'):
            dtype = dtype.dtype
        if isinstance(dtype, cls):
            return True
        elif isinstance(dtype, np.dtype):
            return False
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
    Type for categorical data with the categories and orderedness,
    but not the values

    .. versionadded:: 0.20.0

    Parameters
    ----------
    categories : list or None
    ordered : bool, default False

    Examples
    --------
    >>> t = CategoricalDtype(categories=['b', 'a'], ordered=True)
    >>> s = Series(['a', 'a', 'b', 'b', 'a'])
    >>> s.astype(t)
    0    a
    1    a
    2    b
    3    b
    4    a
    dtype: category
    Categories (2, object): [b < a]

    Notes
    -----
    An instance of ``CategoricalDtype`` compares equal with any other
    instance of ``CategoricalDtype``, regardless of categories or ordered.
    In addition they compare equal to the string ``'category'``.
    To check whether two instances of a ``CategoricalDtype`` match,
    use the ``is`` operator.

    >>> t1 = CategoricalDtype(['a', 'b'], ordered=True)
    >>> t2 = CategoricalDtype(['a', 'c'], ordered=False)
    >>> t1 == t2
    True
    >>> t1 == 'category'
    True
    >>> t1 is t2
    False
    >>> t1 is CategoricalDtype(['a', 'b'], ordered=True)
    True

    A np.dtype duck-typed class, suitable for holding a custom categorical
    dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of np.object
    """
    # TODO: Document public vs. private API
    name = 'category'
    type = CategoricalDtypeType
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')
    _cache = {}

    def __new__(cls, categories=None, ordered=False):
        categories_ = categories if categories is None else tuple(categories)
        t = (categories_, ordered)

        try:
            return cls._cache[t]
        except KeyError:
            c = object.__new__(cls)
            c.categories = categories
            c.ordered = ordered

            cls._cache[t] = c
            return c

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name

        return isinstance(other, CategoricalDtype)

    # def __unicode__(self):
    #     tpl = 'CategoricalDtype({!r}, ordered={})'
    #     return tpl.format(reprlib.repr(self.categories), self.ordered)

    # def __repr__(self):
    #     """ return the base repr for the categories """
    #     tpl = 'CategoricalDtype({!r}, ordered={})'
    #     return tpl.format(reprlib.repr(self.categories), self.ordered)

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

        # set/retrieve from cache
        key = (unit, str(tz))
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
