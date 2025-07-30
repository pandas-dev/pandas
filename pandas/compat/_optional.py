from __future__ import annotations

import importlib
import sys
from typing import (
    TYPE_CHECKING,
    Literal,
    overload,
)
import warnings

from pandas.util._exceptions import find_stack_level

from pandas.util.version import Version

if TYPE_CHECKING:
    import types

# Update install.rst, actions-310-minimum_versions.yaml,
# deps_minimum.toml & pyproject.toml when updating versions!

VERSIONS = {
    "adbc-driver-postgresql": "0.10.0",
    "adbc-driver-sqlite": "0.8.0",
    "bs4": "4.12.3",
    "bottleneck": "1.3.6",
    "fastparquet": "2024.2.0",
    "fsspec": "2023.12.2",
    "html5lib": "1.1",
    "hypothesis": "6.84.0",
    "gcsfs": "2023.12.2",
    "jinja2": "3.1.3",
    "lxml.etree": "4.9.2",
    "matplotlib": "3.8.3",
    "numba": "0.59.0",
    "numexpr": "2.9.0",
    "odfpy": "1.4.1",
    "openpyxl": "3.1.2",
    "psycopg2": "2.9.6",  # (dt dec pq3 ext lo64)
    "pymysql": "1.1.0",
    "pyarrow": "12.0.1",
    "pyiceberg": "0.7.1",
    "pyreadstat": "1.2.6",
    "pytest": "7.3.2",
    "python-calamine": "0.1.7",
    "pytz": "2023.4",
    "pyxlsb": "1.0.10",
    "s3fs": "2023.12.2",
    "scipy": "1.12.0",
    "sqlalchemy": "2.0.0",
    "tables": "3.8.0",
    "tabulate": "0.9.0",
    "xarray": "2024.1.1",
    "xlrd": "2.0.1",
    "xlsxwriter": "3.2.0",
    "zstandard": "0.22.0",
    "qtpy": "2.3.0",
    "pyqt5": "5.15.9",
}

# A mapping from import name to package name (on PyPI) for packages where
# these two names are different.

INSTALL_MAPPING = {
    "bs4": "beautifulsoup4",
    "bottleneck": "Bottleneck",
    "jinja2": "Jinja2",
    "lxml.etree": "lxml",
    "odf": "odfpy",
    "python_calamine": "python-calamine",
    "sqlalchemy": "SQLAlchemy",
    "tables": "pytables",
}

# Mapping of operation contexts to alternative dependencies
OPERATION_CONTEXTS = {
    "excel": {
        "alternatives": ["openpyxl", "xlsxwriter", "calamine", "xlrd", "pyxlsb", "odfpy"],
        "description": "Excel file operations",
    },
    "plotting": {
        "alternatives": ["matplotlib"],
        "description": "plotting operations",
        "fallback": "Use df.describe() for text-based data summaries",
    },
    "html": {
        "alternatives": ["lxml", "html5lib", "beautifulsoup4"],
        "description": "HTML parsing",
    },
    "xml": {
        "alternatives": ["lxml"],
        "description": "XML parsing", 
    },
    "sql": {
        "alternatives": ["sqlalchemy", "psycopg2", "pymysql"],
        "description": "SQL database operations",
    },
    "performance": {
        "alternatives": ["numexpr", "bottleneck", "numba"],
        "description": "performance acceleration",
        "fallback": "Operations will use standard implementations",
    },
    "parquet": {
        "alternatives": ["pyarrow", "fastparquet"],
        "description": "Parquet file operations",
    },
    "feather": {
        "alternatives": ["pyarrow"],
        "description": "Feather file operations",
    },
    "orc": {
        "alternatives": ["pyarrow"],
        "description": "ORC file operations",
    },
    "hdf5": {
        "alternatives": ["tables"],
        "description": "HDF5 file operations",
    },
    "spss": {
        "alternatives": ["pyreadstat"],
        "description": "SPSS file operations",
    },
    "style": {
        "alternatives": ["jinja2"],
        "description": "DataFrame styling operations",
    },
    "compression": {
        "alternatives": ["zstandard"],
        "description": "data compression operations",
    },
    "clipboard": {
        "alternatives": ["pyqt5", "qtpy"],
        "description": "clipboard operations",
    },
}


def _build_context_message(
    name: str, operation_context: str | None, extra: str, install_name: str
) -> str:
    """
    Build an enhanced error message with context-aware alternatives.
    
    Parameters
    ----------
    name : str
        The module name that failed to import.
    operation_context : str or None
        The operation context (e.g., 'excel', 'plotting').
    extra : str
        Additional text to include in the ImportError message.
    install_name : str
        The package name to install.
        
    Returns
    -------
    str
        The enhanced error message.
    """
    base_msg = f"Missing optional dependency '{install_name}'."
    if extra:
        base_msg += f" {extra}"
    
    if operation_context and operation_context in OPERATION_CONTEXTS:
        context_info = OPERATION_CONTEXTS[operation_context]
        # Filter out the failed dependency from alternatives
        alternatives = [
            alt for alt in context_info["alternatives"] 
            if alt != name and alt != install_name
        ]
        
        if alternatives:
            if len(alternatives) == 1:
                alt_msg = f" For {context_info['description']}, try installing {alternatives[0]}."
            elif len(alternatives) == 2:
                alt_msg = f" For {context_info['description']}, try installing {alternatives[0]} or {alternatives[1]}."
            else:
                alt_list = ", ".join(alternatives[:-1]) + f", or {alternatives[-1]}"
                alt_msg = f" For {context_info['description']}, try installing {alt_list}."
            base_msg += alt_msg
            
        if "fallback" in context_info:
            base_msg += f" {context_info['fallback']}."
    
    base_msg += f" Use pip or conda to install {install_name}."
    return base_msg


def get_version(module: types.ModuleType) -> str:
    version = getattr(module, "__version__", None)

    if version is None:
        raise ImportError(f"Can't determine version for {module.__name__}")
    if module.__name__ == "psycopg2":
        # psycopg2 appends " (dt dec pq3 ext lo64)" to it's version
        version = version.split()[0]
    return version


@overload
def import_optional_dependency(
    name: str,
    extra: str = ...,
    min_version: str | None = ...,
    *,
    errors: Literal["raise"] = ...,
    operation_context: str | None = ...,
) -> types.ModuleType: ...


@overload
def import_optional_dependency(
    name: str,
    extra: str = ...,
    min_version: str | None = ...,
    *,
    errors: Literal["warn", "ignore"],
    operation_context: str | None = ...,
) -> types.ModuleType | None: ...


def import_optional_dependency(
    name: str,
    extra: str = "",
    min_version: str | None = None,
    *,
    errors: Literal["raise", "warn", "ignore"] = "raise",
    operation_context: str | None = None,
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
    operation_context : str, default None
        Provide context about the operation requiring this dependency to show
        relevant alternatives in error messages. Supported contexts: 'excel',
        'plotting', 'html', 'xml', 'sql', 'performance', 'parquet', 'feather', 
        'orc', 'hdf5', 'spss', 'style', 'compression', 'clipboard'.
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

    msg = _build_context_message(name, operation_context, extra, install_name)
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
                f"Pandas requires version '{minimum_version}' or newer of '{parent}' "
                f"(version '{version}' currently installed)."
            )
            if errors == "warn":
                warnings.warn(
                    msg,
                    UserWarning,
                    stacklevel=find_stack_level(),
                )
                return None
            elif errors == "raise":
                raise ImportError(msg)
            else:
                return None

    return module
