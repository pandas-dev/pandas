import functools
from .manager import DataModelManager


def register(dmm, typecls):
    """Used as decorator to simplify datamodel registration.
    Returns the object being decorated so that chaining is possible.
    """
    def wraps(fn):
        dmm.register(typecls, fn)
        return fn

    return wraps


default_manager = DataModelManager()

register_default = functools.partial(register, default_manager)
