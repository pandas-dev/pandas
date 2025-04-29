from __future__ import annotations

import importlib
import sys
import types
import warnings
from typing import Literal

if sys.version_info >= (3, 12):
    import importlib.metadata as importlib_metadata
else:
    import importlib_metadata
from packaging.version import Version

PY_VERSION = Version(".".join(map(str, sys.version_info[:3])))

EMSCRIPTEN = sys.platform == "emscripten"

LINUX = sys.platform == "linux"
MACOS = sys.platform == "darwin"
WINDOWS = sys.platform == "win32"


def entry_points(group=None):
    warnings.warn(
        "`dask._compatibility.entry_points` has been replaced by `importlib_metadata.entry_points` and will be removed "
        "in a future version. Please use `importlib_metadata.entry_points` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return importlib_metadata.entry_points(group=group)


VERSIONS = {
    "numpy": "1.21.0",
    "pandas": "2.0.0",
    "bokeh": "3.1.0",
    "jinja2": "2.10.3",
    "pyarrow": "14.0.1",
    "lz4": "4.3.2",
}

# A mapping from import name to package name (on PyPI) for packages where
# these two names are different.

INSTALL_MAPPING = {
    "sqlalchemy": "SQLAlchemy",
    "tables": "pytables",
}


def get_version(module: types.ModuleType) -> str:
    version = getattr(module, "__version__", None)

    if version is None:
        raise ImportError(f"Can't determine version for {module.__name__}")
    if module.__name__ == "psycopg2":
        # psycopg2 appends " (dt dec pq3 ext lo64)" to it's version
        version = version.split()[0]
    return version


def import_optional_dependency(
    name: str,
    extra: str = "",
    min_version: str | None = None,
    *,
    errors: Literal["raise", "warn", "ignore"] = "raise",
) -> types.ModuleType | None:
    """
    Import an optional dependency.

    By default, if a dependency is missing an ImportError with a nice
    message will be raised. If a dependency is present, but too old,
    we raise.

    Parameters
    ----------
    name : str
        The module name.
    extra : str
        Additional text to include in the ImportError message.
    errors : str {'raise', 'warn', 'ignore'}
        What to do when a dependency is not found or its version is too old.

        * raise : Raise an ImportError
        * warn : Only applicable when a module's version is to old.
          Warns that the version is too old and returns None
        * ignore: If the module is not installed, return None, otherwise,
          return the module, even if the version is too old.
          It's expected that users validate the version locally when
          using ``errors="ignore"`` (see. ``io/html.py``)
    min_version : str, default None
        Specify a minimum version that is different from the global pandas
        minimum version required.
    Returns
    -------
    maybe_module : Optional[ModuleType]
        The imported module, when found and the version is correct.
        None is returned when the package is not found and `errors`
        is False, or when the package's version is too old and `errors`
        is ``'warn'`` or ``'ignore'``.
    """
    assert errors in {"warn", "raise", "ignore"}

    package_name = INSTALL_MAPPING.get(name)
    install_name = package_name if package_name is not None else name

    msg = (
        f"Missing optional dependency '{install_name}'. {extra} "
        f"Use pip or conda to install {install_name}."
    )
    try:
        module = importlib.import_module(name)
    except ImportError as err:
        if errors == "raise":
            raise ImportError(msg) from err
        return None

    # Handle submodules: if we have submodule, grab parent module from sys.modules
    parent = name.split(".")[0]
    if parent != name:
        install_name = parent
        module_to_get = sys.modules[install_name]
    else:
        module_to_get = module
    minimum_version = min_version if min_version is not None else VERSIONS.get(parent)
    if minimum_version:
        version = get_version(module_to_get)
        if version and Version(version) < Version(minimum_version):
            msg = (
                f"Dask requires version '{minimum_version}' or newer of '{parent}' "
                f"(version '{version}' currently installed)."
            )
            if errors == "warn":
                warnings.warn(msg, UserWarning)
                return None
            elif errors == "raise":
                raise ImportError(msg)
            else:
                return None

    return module
