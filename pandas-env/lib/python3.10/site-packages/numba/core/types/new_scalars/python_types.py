"""
    Type definitions for Python types.
"""

from numba.core.types.abstract import Literal
from numba.core.types.new_scalars.scalars \
    import (Integer, IntegerLiteral, Boolean,
            BooleanLiteral, Float, Complex,
            parse_integer_signed)
from functools import total_ordering
from numba.core.typeconv import Conversion


@total_ordering
class PythonInteger(Integer):
    def __init__(self, name, bitwidth=None, signed=None):
        super(PythonInteger, self).__init__(name)
        if bitwidth is None:
            bitwidth = 64
        if signed is None:
            signed = parse_integer_signed(name)
        self.bitwidth = bitwidth
        self.signed = signed

    def cast_python_value(self, value):
        return int(value)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        if self.signed != other.signed:
            return NotImplemented
        return self.bitwidth < other.bitwidth

    @property
    def maxval(self):
        """
        The maximum value representable by this type.
        """
        if self.signed:
            return (1 << (self.bitwidth - 1)) - 1
        else:
            return (1 << self.bitwidth) - 1

    @property
    def minval(self):
        """
        The minimal value representable by this type.
        """
        if self.signed:
            return -(1 << (self.bitwidth - 1))
        else:
            return 0


class PythonIntegerLiteral(IntegerLiteral, PythonInteger):
    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[int]({})'.format(value)
        basetype = self.literal_type
        PythonInteger.__init__(self,
                               name=name,
                               bitwidth=basetype.bitwidth,
                               signed=basetype.signed,)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


Literal.ctor_map[int] = PythonIntegerLiteral


class PythonBoolean(Boolean):
    def cast_python_value(self, value):
        return bool(value)


class PythonBooleanLiteral(BooleanLiteral, PythonBoolean):

    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[bool]({})'.format(value)
        PythonBoolean.__init__(self, name=name)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


Literal.ctor_map[bool] = PythonBooleanLiteral


@total_ordering
class PythonFloat(Float):
    def __init__(self, *args, **kws):
        super(PythonFloat, self).__init__(*args, **kws)
        # Determine bitwidth
        assert self.name.startswith('py_float')
        bitwidth = 64
        self.bitwidth = bitwidth

    def cast_python_value(self, value):
        return float(value)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth


@total_ordering
class PythonComplex(Complex):
    def __init__(self, name, underlying_float, **kwargs):
        super(PythonComplex, self).__init__(name, **kwargs)
        self.underlying_float = underlying_float
        # Determine bitwidth
        assert self.name.startswith('py_complex')
        bitwidth = 128
        self.bitwidth = bitwidth

    def cast_python_value(self, value):
        return complex(value)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth
