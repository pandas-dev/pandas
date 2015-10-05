""" define extension dtypes """

import re
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
        raise NotImplementedError("sub-classes should implement an __hash__ method")

    def __eq__(self, other):
        raise NotImplementedError("sub-classes should implement an __eq__ method")

    @classmethod
    def is_dtype(cls, dtype):
        """ Return a boolean if we if the passed type is an actual dtype that we can match (via string or type) """
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
    A np.dtype duck-typed class, suitable for holding a custom categorical dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of np.object
    """
    name = 'category'
    type = CategoricalDtypeType
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name

        return isinstance(other, CategoricalDtype)

    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if its not possible """
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
    A np.dtype duck-typed class, suitable for holding a custom datetime with tz dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of np.datetime64[ns]
    """
    type = DatetimeTZDtypeType
    kind = 'M'
    str = '|M8[ns]'
    num = 101
    base = np.dtype('M8[ns]')
    _metadata = ['unit','tz']
    _match = re.compile("(datetime64|M8)\[(?P<unit>.+), (?P<tz>.+)\]")

    def __init__(self, unit, tz=None):
        """
        Parameters
        ----------
        unit : string unit that this represents, currently must be 'ns'
        tz : string tz that this represents
        """

        if isinstance(unit, DatetimeTZDtype):
            self.unit, self.tz = unit.unit, unit.tz
            return

        if tz is None:

            # we were passed a string that we can construct
            try:
                m = self._match.search(unit)
                if m is not None:
                    self.unit = m.groupdict()['unit']
                    self.tz = m.groupdict()['tz']
                    return
            except:
                raise ValueError("could not construct DatetimeTZDtype")

            raise ValueError("DatetimeTZDtype constructor must have a tz supplied")

        if unit != 'ns':
            raise ValueError("DatetimeTZDtype only supports ns units")
        self.unit = unit
        self.tz = tz

    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if its not possible """
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

        return isinstance(other, DatetimeTZDtype) and self.unit == other.unit and self.tz == other.tz
