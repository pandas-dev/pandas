import numbers


def _create_binary_propagating_op(name):
    def method(self, other):
        if isinstance(other, numbers.Number) or other is NA or isinstance(other, str):
            return NA

        return NotImplemented

    method.__name__ = name
    return method


def _create_unary_propagating_op(name):
    def method(self):
        return NA

    method.__name__ = name
    return method


class NAType:

    _instance = None

    def __new__(cls, *args, **kwargs):
        if NAType._instance is None:
            NAType._instance = object.__new__(cls, *args, **kwargs)
        return NAType._instance

    def __repr__(self) -> str:
        return "NA"

    def __str__(self) -> str:
        return "NA"

    def __bool__(self):
        raise TypeError("boolean value of NA is ambiguous")

    def __hash__(self):
        # TODO what should we use here to hash?
        return 0

    # Binary arithmetic and comparison ops -> propagate

    __add__ = _create_binary_propagating_op("__add__")
    __radd__ = _create_binary_propagating_op("__radd__")
    __sub__ = _create_binary_propagating_op("__sub__")
    __rsub__ = _create_binary_propagating_op("__rsub__")
    __mul__ = _create_binary_propagating_op("__mul__")
    __rmul__ = _create_binary_propagating_op("__rmul__")
    __matmul__ = _create_binary_propagating_op("__matmul__")
    __rmatmul__ = _create_binary_propagating_op("__rmatmul__")
    __truediv__ = _create_binary_propagating_op("__truediv__")
    __rtruediv__ = _create_binary_propagating_op("__rtruediv__")
    __floordiv__ = _create_binary_propagating_op("__floordiv__")
    __rfloordiv__ = _create_binary_propagating_op("__rfloordiv__")
    __mod__ = _create_binary_propagating_op("__mod__")
    __rmod__ = _create_binary_propagating_op("__rmod__")
    __divmod__ = _create_binary_propagating_op("__divmod__")
    __rdivmod__ = _create_binary_propagating_op("__rdivmod__")
    __pow__ = _create_binary_propagating_op("__pow__")
    __rpow__ = _create_binary_propagating_op("__rpow__")

    # __lshift__ = _create_binary_propagating_op("__add__")
    # __rshift__ = _create_binary_propagating_op("__add__")

    __eq__ = _create_binary_propagating_op("__eq__")
    __ne__ = _create_binary_propagating_op("__ne__")
    __le__ = _create_binary_propagating_op("__le__")
    __lt__ = _create_binary_propagating_op("__lt__")
    __gt__ = _create_binary_propagating_op("__gt__")
    __ge__ = _create_binary_propagating_op("__ge__")

    # Unary ops

    __neg__ = _create_unary_propagating_op("__neg__")
    __pos__ = _create_unary_propagating_op("__pos__")
    __abs__ = _create_unary_propagating_op("__abs__")
    __invert__ = _create_unary_propagating_op("__invert__")

    # Logical ops using Kleene logic

    def __and__(self, other):
        if other is False:
            return False
        elif other is True or other is NA:
            return NA
        else:
            return NotImplemented

    __rand__ = __and__

    def __or__(self, other):
        if other is True:
            return True
        elif other is False or other is NA:
            return NA
        else:
            return NotImplemented

    __ror__ = __or__

    def __xor__(self, other):
        if other is False or other is True or other is NA:
            return NA
        return NotImplemented

    __rxor__ = __xor__


NA = NAType()
