"""
    Type definitions for machine types.
"""

from numba.core.types.new_scalars.scalars \
    import (Integer, IntegerLiteral, Boolean,
            BooleanLiteral, Float, Complex,
            parse_integer_bitwidth, parse_integer_signed)
from functools import total_ordering
from numba.core.typeconv import Conversion


@total_ordering
class MachineInteger(Integer):
    def __init__(self, name, bitwidth=None, signed=None):
        super(MachineInteger, self).__init__(name)
        if bitwidth is None:
            bitwidth = parse_integer_bitwidth(name)
        if signed is None:
            signed = parse_integer_signed(name)
        self.bitwidth = bitwidth
        self.signed = signed

    @classmethod
    def from_bitwidth(cls, bitwidth, signed=True):
        name = ('int%d' if signed else 'uint%d') % bitwidth
        return cls(name)

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


class MachineIntegerLiteral(IntegerLiteral, MachineInteger):
    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[machine_int]({})'.format(value)
        basetype = self.literal_type
        MachineInteger.__init__(self,
                                name=name,
                                bitwidth=basetype.bitwidth,
                                signed=basetype.signed,)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


class MachineBoolean(Boolean):
    pass


class MachineBooleanLiteral(BooleanLiteral, MachineBoolean):

    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[machine_bool]({})'.format(value)
        MachineBoolean.__init__(self,
                                name=name)

    def can_convert_to(self, typingctx, other):
        conv = typingctx.can_convert(self.literal_type, other)
        if conv is not None:
            return max(conv, Conversion.promote)


@total_ordering
class MachineFloat(Float):
    def __init__(self, *args, **kws):
        super(MachineFloat, self).__init__(*args, **kws)
        # Determine bitwidth
        assert self.name.startswith('c_float')
        bitwidth = int(self.name[8:])
        self.bitwidth = bitwidth

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth


@total_ordering
class MachineComplex(Complex):
    def __init__(self, name, underlying_float, **kwargs):
        super(MachineComplex, self).__init__(name, **kwargs)
        self.underlying_float = underlying_float
        # Determine bitwidth
        assert self.name.startswith('c_complex')
        bitwidth = int(self.name[10:])
        self.bitwidth = bitwidth

    def __lt__(self, other):
        if self.__class__ is not other.__class__:
            return NotImplemented
        return self.bitwidth < other.bitwidth
