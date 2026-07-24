"""
Serialization support for compiled functions.
"""
import sys
import abc
import io
import copyreg


import pickle
from numba import cloudpickle
from llvmlite import ir


#
# Pickle support
#

def _rebuild_reduction(cls, *args):
    """
    Global hook to rebuild a given class from its __reduce__ arguments.
    """
    return cls._rebuild(*args)


# Keep unpickled object via `numba_unpickle` alive.
_unpickled_memo = {}


def _numba_unpickle(address, bytedata, hashed):
    """Used by `numba_unpickle` from _helperlib.c

    Parameters
    ----------
    address : int
    bytedata : bytes
    hashed : bytes

    Returns
    -------
    obj : object
        unpickled object
    """
    key = (address, hashed)
    try:
        obj = _unpickled_memo[key]
    except KeyError:
        _unpickled_memo[key] = obj = cloudpickle.loads(bytedata)
    return obj


def dumps(obj):
    """Similar to `pickle.dumps()`. Returns the serialized object in bytes.
    """
    pickler = NumbaPickler
    with io.BytesIO() as buf:
        p = pickler(buf, protocol=4)
        p.dump(obj)
        pickled = buf.getvalue()

    return pickled


def runtime_build_excinfo_struct(static_exc, exc_args):
    exc, static_args, locinfo = cloudpickle.loads(static_exc)
    real_args = []
    exc_args_iter = iter(exc_args)
    for arg in static_args:
        if isinstance(arg, ir.Value):
            real_args.append(next(exc_args_iter))
        else:
            real_args.append(arg)
    return (exc, tuple(real_args), locinfo)


# Alias to pickle.loads to allow `serialize.loads()`
loads = cloudpickle.loads


class _CustomPickled:
    """A wrapper for objects that must be pickled with `NumbaPickler`.

    Standard `pickle` will pick up the implementation registered via `copyreg`.
    This will spawn a `NumbaPickler` instance to serialize the data.

    `NumbaPickler` overrides the handling of this type so as not to spawn a
    new pickler for the object when it is already being pickled by a
    `NumbaPickler`.
    """

    __slots__ = 'ctor', 'states'

    def __init__(self, ctor, states):
        self.ctor = ctor
        self.states = states

    def _reduce(self):
        return _CustomPickled._rebuild, (self.ctor, self.states)

    @classmethod
    def _rebuild(cls, ctor, states):
        return cls(ctor, states)


def _unpickle__CustomPickled(serialized):
    """standard unpickling for `_CustomPickled`.

    Uses `NumbaPickler` to load.
    """
    ctor, states = loads(serialized)
    return _CustomPickled(ctor, states)


def _pickle__CustomPickled(cp):
    """standard pickling for `_CustomPickled`.

    Uses `NumbaPickler` to dump.
    """
    serialized = dumps((cp.ctor, cp.states))
    return _unpickle__CustomPickled, (serialized,)


# Register custom pickling for the standard pickler.
copyreg.pickle(_CustomPickled, _pickle__CustomPickled)


def custom_reduce(cls, states):
    """For customizing object serialization in `__reduce__`.

    Object states provided here are used as keyword arguments to the
    `._rebuild()` class method.

    Parameters
    ----------
    states : dict
        Dictionary of object states to be serialized.

    Returns
    -------
    result : tuple
        This tuple conforms to the return type requirement for `__reduce__`.
    """
    return custom_rebuild, (_CustomPickled(cls, states),)


def custom_rebuild(custom_pickled):
    """Customized object deserialization.

    This function is referenced internally by `custom_reduce()`.
    """
    cls, states = custom_pickled.ctor, custom_pickled.states
    return cls._rebuild(**states)


def is_serialiable(obj):
    """Check if *obj* can be serialized.

    Parameters
    ----------
    obj : object

    Returns
    --------
    can_serialize : bool
    """
    with io.BytesIO() as fout:
        pickler = NumbaPickler(fout)
        try:
            pickler.dump(obj)
        except pickle.PicklingError:
            return False
        else:
            return True


def _no_pickle(obj):
    raise pickle.PicklingError(f"Pickling of {type(obj)} is unsupported")


def disable_pickling(typ):
    """This is called on a type to disable pickling
    """
    NumbaPickler.disabled_types.add(typ)
    # Return `typ` to allow use as a decorator
    return typ


class NumbaPickler(cloudpickle.CloudPickler):
    disabled_types = set()
    """A set of types that pickling cannot is disabled.
    """

    def reducer_override(self, obj):
        # Overridden to disable pickling of certain types
        if type(obj) in self.disabled_types:
            _no_pickle(obj)  # noreturn
        return super().reducer_override(obj)


def _custom_reduce__custompickled(cp):
    return cp._reduce()


NumbaPickler.dispatch_table[_CustomPickled] = _custom_reduce__custompickled


class ReduceMixin(abc.ABC):
    """A mixin class for objects that should be reduced by the NumbaPickler
    instead of the standard pickler.
    """
    # Subclass MUST override the below methods

    @abc.abstractmethod
    def _reduce_states(self):
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def _rebuild(cls, **kwargs):
        raise NotImplementedError

    # Subclass can override the below methods

    def _reduce_class(self):
        return self.__class__

    # Private methods

    def __reduce__(self):
        return custom_reduce(self._reduce_class(), self._reduce_states())


class PickleCallableByPath:
    """Wrap a callable object to be pickled by path to workaround limitation
    in pickling due to non-pickleable objects in function non-locals.

    Note:
    - Do not use this as a decorator.
    - Wrapped object must be a global that exist in its parent module and it
      can be imported by `from the_module import the_object`.

    Usage:

    >>> def my_fn(x):
    >>>     ...
    >>> wrapped_fn = PickleCallableByPath(my_fn)
    >>> # refer to `wrapped_fn` instead of `my_fn`
    """
    def __init__(self, fn):
        self._fn = fn

    def __call__(self, *args, **kwargs):
        return self._fn(*args, **kwargs)

    def __reduce__(self):
        return type(self)._rebuild, (self._fn.__module__, self._fn.__name__,)

    @classmethod
    def _rebuild(cls, modname, fn_path):
        return cls(getattr(sys.modules[modname], fn_path))
