import weakref
import importlib

from numba import _dynfunc


class Environment(_dynfunc.Environment):
    """Stores globals and constant pyobjects for runtime.

    It is often needed to convert b/w nopython objects and pyobjects.
    """
    __slots__ = ('env_name', '__weakref__')
    # A weak-value dictionary to store live environment with env_name as the
    # key.
    _memo = weakref.WeakValueDictionary()

    @classmethod
    def from_fndesc(cls, fndesc):
        try:
            # Avoid creating new Env
            return cls._memo[fndesc.env_name]
        except KeyError:
            inst = cls(fndesc.lookup_globals())
            inst.env_name = fndesc.env_name
            cls._memo[fndesc.env_name] = inst
            return inst

    def can_cache(self):
        is_dyn = '__name__' not in self.globals
        return not is_dyn

    def __reduce__(self):
        return _rebuild_env, (
            self.globals.get('__name__'),
            self.consts,
            self.env_name,
        )

    def __del__(self):
        return

    def __repr__(self):
        return f"<Environment {self.env_name!r} >"


def _rebuild_env(modname, consts, env_name):
    env = lookup_environment(env_name)
    if env is not None:
        return env

    mod = importlib.import_module(modname)
    env = Environment(mod.__dict__)
    env.consts[:] = consts
    env.env_name = env_name
    # Cache loaded object
    Environment._memo[env_name] = env
    return env


def lookup_environment(env_name):
    """Returns the Environment object for the given name;
    or None if not found
    """
    return Environment._memo.get(env_name)
