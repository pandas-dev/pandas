import enum
import re
import numpy as np

from numba.core.types.abstract import Dummy, Hashable, Literal, Number, Type
from functools import total_ordering, cached_property
from numba.core import utils
from numba.core.typeconv import Conversion
from numba.np import npdatetime_helpers


class Boolean(Hashable):
    pass

def parse_integer_bitwidth(name):
    bitwidth = int(re.findall(r'\d+', name)[-1])
    return bitwidth


def parse_integer_signed(name):
    signed = name.startswith('int')
    return signed


class Integer(Number):
    pass


class IntegerLiteral(Literal, Integer):
    pass


class BooleanLiteral(Literal, Boolean):
    pass


class Float(Number):
    pass


class Complex(Number):
    pass


class _NPDatetimeBase(Type):
    """
    Common base class for np.datetime64 and np.timedelta64.
    """

    def __init__(self, unit, *args, **kws):
        name = '%s[%s]' % (self.type_name, unit)
        self.unit = unit
        self.unit_code = npdatetime_helpers.DATETIME_UNITS[self.unit]
        super(_NPDatetimeBase, self).__init__(name, *args, **kws)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        # A coarser-grained unit is "smaller", i.e. less precise values
        # can be represented (but the magnitude of representable values is
        # also greater...).
        return self.unit_code < other.unit_code

    def cast_python_value(self, value):
        cls = getattr(np, self.type_name)
        if self.unit:
            return cls(value, self.unit)
        else:
            return cls(value)


@total_ordering
class NPTimedelta(_NPDatetimeBase):
    type_name = 'timedelta64'

@total_ordering
class NPDatetime(_NPDatetimeBase):
    type_name = 'datetime64'


class EnumClass(Dummy):
    """
    Type class for Enum classes.
    """
    basename = "Enum class"

    def __init__(self, cls, dtype):
        assert isinstance(cls, type)
        assert isinstance(dtype, Type)
        self.instance_class = cls
        self.dtype = dtype
        name = "%s<%s>(%s)" % (self.basename, self.dtype, self.instance_class.__name__)
        super(EnumClass, self).__init__(name)

    @property
    def key(self):
        return self.instance_class, self.dtype

    @cached_property
    def member_type(self):
        """
        The type of this class' members.
        """
        return EnumMember(self.instance_class, self.dtype)


class IntEnumClass(EnumClass):
    """
    Type class for IntEnum classes.
    """
    basename = "IntEnum class"

    @cached_property
    def member_type(self):
        """
        The type of this class' members.
        """
        return IntEnumMember(self.instance_class, self.dtype)


class EnumMember(Type):
    """
    Type class for Enum members.
    """
    basename = "Enum"
    class_type_class = EnumClass

    def __init__(self, cls, dtype):
        assert isinstance(cls, type)
        assert isinstance(dtype, Type)
        self.instance_class = cls
        self.dtype = dtype
        name = "%s<%s>(%s)" % (self.basename, self.dtype, self.instance_class.__name__)
        super(EnumMember, self).__init__(name)

    @property
    def key(self):
        return self.instance_class, self.dtype

    @property
    def class_type(self):
        """
        The type of this member's class.
        """
        return self.class_type_class(self.instance_class, self.dtype)


class IntEnumMember(EnumMember):
    """
    Type class for IntEnum members.
    """
    basename = "IntEnum"
    class_type_class = IntEnumClass

    def can_convert_to(self, typingctx, other):
        """
        Convert IntEnum members to plain integers.
        """
        if issubclass(self.instance_class, enum.IntEnum):
            conv = typingctx.can_convert(self.dtype, other)
            return max(conv, Conversion.safe)
