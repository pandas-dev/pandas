import inspect
import os
from collections.abc import MutableMapping

NO_DEFAULT = object()


# must not inherit from AttributeError, so not to mess with python's attribute-lookup flow
class EnvironmentVariableError(KeyError):
    pass


class TypedEnv(MutableMapping):
    """
    This object can be used in 'exploratory' mode:

        nv = TypedEnv()
        print(nv.HOME)

    It can also be used as a parser and validator of environment variables:

    class MyEnv(TypedEnv):
        username = TypedEnv.Str("USER", default='anonymous')
        path = TypedEnv.CSV("PATH", separator=":")
        tmp = TypedEnv.Str("TMP TEMP".split())  # support 'fallback' var-names

    nv = MyEnv()

    print(nv.username)

    for p in nv.path:
        print(p)

    try:
        print(p.tmp)
    except EnvironmentVariableError:
        print("TMP/TEMP is not defined")
    else:
        assert False
    """

    __slots__ = ["_env", "_defined_keys"]

    class _BaseVar:
        def __init__(self, name, default=NO_DEFAULT):
            self.names = tuple(name) if isinstance(name, (tuple, list)) else (name,)
            self.name = self.names[0]
            self.default = default

        def convert(self, value):  # pylint:disable=no-self-use
            return value

        def __get__(self, instance, owner):
            if not instance:
                return self
            try:
                return self.convert(instance._raw_get(*self.names))
            except EnvironmentVariableError:
                if self.default is NO_DEFAULT:
                    raise
                return self.default

        def __set__(self, instance, value):
            instance[self.name] = value

    class Str(_BaseVar):
        pass

    class Bool(_BaseVar):
        """
        Converts 'yes|true|1|no|false|0' to the appropriate boolean value.
        Case-insensitive. Throws a ``ValueError`` for any other value.
        """

        def convert(self, value):
            value = value.lower()
            if value not in {"yes", "no", "true", "false", "1", "0"}:
                raise ValueError(f"Unrecognized boolean value: {value!r}")
            return value in {"yes", "true", "1"}

        def __set__(self, instance, value):
            instance[self.name] = "yes" if value else "no"

    class Int(_BaseVar):
        convert = staticmethod(int)

    class Float(_BaseVar):
        convert = staticmethod(float)

    class CSV(_BaseVar):
        """
        Comma-separated-strings get split using the ``separator`` (',' by default) into
        a list of objects of type ``type`` (``str`` by default).
        """

        def __init__(self, name, default=NO_DEFAULT, type=str, separator=","):  # pylint:disable=redefined-builtin
            super().__init__(name, default=default)
            self.type = type
            self.separator = separator

        def __set__(self, instance, value):
            instance[self.name] = self.separator.join(map(str, value))

        def convert(self, value):
            return [self.type(v.strip()) for v in value.split(self.separator)]

    # =========

    def __init__(self, env=None):
        if env is None:
            env = os.environ
        self._env = env
        self._defined_keys = {
            k
            for (k, v) in inspect.getmembers(self.__class__)
            if isinstance(v, self._BaseVar)
        }

    def __iter__(self):
        return iter(dir(self))

    def __len__(self):
        return len(self._env)

    def __delitem__(self, name):
        del self._env[name]

    def __setitem__(self, name, value):
        self._env[name] = str(value)

    def _raw_get(self, *key_names):
        for key in key_names:
            value = self._env.get(key, NO_DEFAULT)
            if value is not NO_DEFAULT:
                return value
        raise EnvironmentVariableError(key_names[0])

    def __contains__(self, key):
        try:
            self._raw_get(key)
        except EnvironmentVariableError:
            return False
        return True

    def __getattr__(self, name):
        # if we're here then there was no descriptor defined
        try:
            return self._raw_get(name)
        except EnvironmentVariableError:
            raise AttributeError(
                f"{self.__class__} has no attribute {name!r}"
            ) from None

    def __getitem__(self, key):
        return getattr(self, key)  # delegate through the descriptors

    def get(self, key, default=None):
        try:
            return self[key]
        except EnvironmentVariableError:
            return default

    def __dir__(self):
        if self._defined_keys:
            # return only defined
            return sorted(self._defined_keys)
        # return whatever is in the environment (for convenience)
        members = set(self._env.keys())
        members.update(dir(self.__class__))
        return sorted(members)
