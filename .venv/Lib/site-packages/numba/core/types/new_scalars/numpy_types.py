"""
    Type definitions for NumPy types.
"""

import numpy as np

from numba.core.types.abstract import Literal
from numba.core.types.new_scalars.scalars \
    import (Integer, IntegerLiteral, Boolean,
            BooleanLiteral, Float, Complex,
            parse_integer_bitwidth, parse_integer_signed)
from functools import total_ordering
from numba.core.typeconv import Conversion


@total_ordering
class NumPyInteger(Integer):
    def __init__(self, name, bitwidth=None, signed=None):
        super(NumPyInteger, self).__init__(name)
        if bitwidth is None:
            bitwidth = parse_integer_bitwidth(name)
        if signed is None:
            signed = parse_integer_signed(name)
        self.bitwidth = bitwidth
        self.signed = signed

    @classmethod
    def from_bitwidth(cls, bitwidth, signed=True):
        name = ('np_int%d' if signed else 'np_uint%d') % bitwidth
        return cls(name)

    def cast_python_value(self, value):
        sign_char = "" if self.signed else "u"
        return getattr(
            np,
            sign_char + "int" + str(self.bitwidth)
        )(value)

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


class NumPyIntegerLiteral(IntegerLiteral):
    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[int]({})'.format(value)
        basetype = self.literal_type
        NumPyInteger.__init__(self,
                              name=name,
                              bitwidth=basetype.bitwidth,
                              signed=basetype.signed,)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


Literal.ctor_map[np.integer] = NumPyIntegerLiteral


class NumPyBoolean(Boolean):
    def cast_python_value(self, value):
        return np.bool_(value)


class NumPyBooleanLiteral(BooleanLiteral, NumPyBoolean):

    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[np.bool_]({})'.format(value)
        NumPyBoolean.__init__(self,
                              name=name)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


Literal.ctor_map[np.bool_] = NumPyBooleanLiteral


@total_ordering
class NumPyFloat(Float):
    def __init__(self, *args, **kws):
        super(NumPyFloat, self).__init__(*args, **kws)
        # Determine bitwidth
        assert self.name.startswith('np_float')
        bitwidth = int(self.name[8:])
        self.bitwidth = bitwidth

    def cast_python_value(self, value):
        return getattr(np, "float" + str(self.bitwidth))(value)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth


@total_ordering
class NumPyComplex(Complex):
    def __init__(self, name, underlying_float, **kwargs):
        super(NumPyComplex, self).__init__(name, **kwargs)
        self.underlying_float = underlying_float
        # Determine bitwidth
        assert self.name.startswith('np_complex')
        bitwidth = int(self.name[10:])
        self.bitwidth = bitwidth

    def cast_python_value(self, value):
        return getattr(np, "complex" + str(self.bitwidth))(value)

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth
