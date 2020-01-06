import importlib as _importlib
from importlib.abc import Loader as _Loader, MetaPathFinder
from importlib.machinery import ModuleSpec as _ModuleSpec
import sys as _sys
import warnings as _warnings

from pandas.util._decorators import Appender, Substitution, cache_readonly

from pandas.core.util.hashing import hash_array, hash_pandas_object

# Custom import hook for the deprecated pandas.util.testing module.
# The custom Finder only runs when python's standard import fail,
# since we're adding to the end of sys.meta_path. In TestingLoader.create_module
# we simply warn and then return the real testing module, pandas._testing.
# This differs from _DeprecatedModule by using the import system, which
# allows it to work with both `import foo; foo.bar` and `from foo import bar`.
# But because Python caches imports, the warning appears only once.


class _TestingFinder(MetaPathFinder):
    @classmethod
    def find_spec(cls, fullname, path=None, target=None):

        name_parts = fullname.split(".")
        if name_parts[:3] != ["pandas", "util", "testing"] or len(name_parts) > 3:
            return None
        else:
            # return ModuleSpec(fullname, DataPackageImporter())
            return _ModuleSpec(fullname, _TestingLoader())


class _TestingLoader(_Loader):
    @classmethod
    def create_module(cls, spec):
        module = _importlib.import_module("pandas._testing")
        _warnings.warn(
            (
                "pandas.util.testing is deprecated. Use the "
                "public methods from pandas.testing instead."
            ),
            FutureWarning,
            stacklevel=2,
        )
        return module

    @classmethod
    def exec_module(cls, module):
        pass


_sys.meta_path.append(_TestingFinder())


__all__ = [
    "Appender",
    "Substitution",
    "cache_readonly",
    "hash_array",
    "hash_pandas_object",
]
