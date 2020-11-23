import distutils.version
import importlib
import types
import warnings

# Update install.rst when updating versions!

VERSIONS = {
    "bs4": "4.6.0",
    "bottleneck": "1.2.1",
    "fsspec": "0.7.4",
    "fastparquet": "0.3.2",
    "gcsfs": "0.6.0",
    "lxml.etree": "4.3.0",
    "matplotlib": "2.2.3",
    "numexpr": "2.6.8",
    "odfpy": "1.3.0",
    "openpyxl": "2.5.7",
    "pandas_gbq": "0.12.0",
    "pyarrow": "0.15.0",
    "pytest": "5.0.1",
    "pyxlsb": "1.0.6",
    "s3fs": "0.4.0",
    "scipy": "1.2.0",
    "sqlalchemy": "1.2.8",
    "tables": "3.5.1",
    "tabulate": "0.8.3",
    "xarray": "0.12.3",
    "xlrd": "1.2.0",
    "xlwt": "1.3.0",
    "xlsxwriter": "1.0.2",
    "numba": "0.46.0",
}

# A mapping from import name to package name (on PyPI) for packages where
# these two names are different.

INSTALL_MAPPING = {
    "bs4": "beautifulsoup4",
    "bottleneck": "Bottleneck",
    "lxml.etree": "lxml",
    "odf": "odfpy",
    "pandas_gbq": "pandas-gbq",
    "sqlalchemy": "SQLAlchemy",
    "jinja2": "Jinja2",
}


def _get_version(module: types.ModuleType) -> str:
    version = getattr(module, "__version__", None)
    if version is None:
        # xlrd uses a capitalized attribute name
        version = getattr(module, "__VERSION__", None)

    if version is None:
        raise ImportError(f"Can't determine version for {module.__name__}")
    return version


def import_optional_dependency(
    name: str, extra: str = "", raise_on_missing: bool = True, on_version: str = "raise"
):
    """
    Import an optional dependency.

    By default, if a dependency is missing an ImportError with a nice
    message will be raised. If a dependency is present, but too old,
    we raise.

    Parameters
    ----------
    name : str
        The module name. This should be top-level only, so that the
        version may be checked.
    extra : str
        Additional text to include in the ImportError message.
    raise_on_missing : bool, default True
        Whether to raise if the optional dependency is not found.
        When False and the module is not present, None is returned.
    on_version : str {'raise', 'warn'}
        What to do when a dependency's version is too old.

        * raise : Raise an ImportError
        * warn : Warn that the version is too old. Returns None
        * ignore: Return the module, even if the version is too old.
          It's expected that users validate the version locally when
          using ``on_version="ignore"`` (see. ``io/html.py``)

    Returns
    -------
    maybe_module : Optional[ModuleType]
        The imported module, when found and the version is correct.
        None is returned when the package is not found and `raise_on_missing`
        is False, or when the package's version is too old and `on_version`
        is ``'warn'``.
    """

    package_name = INSTALL_MAPPING.get(name)
    install_name = package_name if package_name is not None else name

    msg = (
        f"Missing optional dependency '{install_name}'. {extra} "
        f"Use pip or conda to install {install_name}."
    )
    try:
        module = importlib.import_module(name)
    except ImportError:
        if raise_on_missing:
            raise ImportError(msg) from None
        else:
            return None

    minimum_version = VERSIONS.get(name)
    if minimum_version:
        version = _get_version(module)
        if distutils.version.LooseVersion(version) < minimum_version:
            assert on_version in {"warn", "raise", "ignore"}
            msg = (
                f"Pandas requires version '{minimum_version}' or newer of '{name}' "
                f"(version '{version}' currently installed)."
            )
            if on_version == "warn":
                warnings.warn(msg, UserWarning)
                return None
            elif on_version == "raise":
                raise ImportError(msg)

    return module
