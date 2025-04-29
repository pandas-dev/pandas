"""
This implements the typing template for `dict()`.
"""

from .. import types, errors
from .templates import (
    AbstractTemplate,
    Registry,
    signature,
)

registry = Registry()
infer = registry.register
infer_global = registry.register_global
infer_getattr = registry.register_attr


_message_dict_support = """
Unsupported use of `dict()` with positional or keyword argument(s). \
The only supported uses are `dict()` or `dict(iterable)`.
""".strip()


@infer_global(dict)
class DictBuiltin(AbstractTemplate):
    def generic(self, args, kws):
        if kws:
            raise errors.TypingError(_message_dict_support)
        if args:
            iterable, = args
            if isinstance(iterable, types.IterableType):
                dtype = iterable.iterator_type.yield_type
                if isinstance(dtype, types.UniTuple):
                    length = dtype.count
                    if length != 2:
                        msg = ("dictionary update sequence element has length "
                               f"{length}; 2 is required")
                        raise errors.TypingError(msg)
                    k = v = dtype.key[0]
                elif isinstance(dtype, types.Tuple):
                    k, v = dtype.key
                else:
                    raise errors.TypingError(_message_dict_support)

                # dict key must be hashable
                if not isinstance(k, types.Hashable):
                    msg = f"Unhashable type: {k}"
                    raise errors.TypingError(msg)

                return signature(types.DictType(k, v), iterable)
            else:
                msg = ("Non-iterable args used in dict(iterable) "
                       f"constructor. Got 'dict({args[0]})'")
                raise errors.TypingError(msg)
        return signature(types.DictType(types.undefined, types.undefined))
