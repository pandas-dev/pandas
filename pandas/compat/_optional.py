import distutils.version
import importlib
import types
import warnings

# Update install.rst when updating versions!

VERSIONS = {
    "bs4": "4.6.0",
    "bottleneck": "1.2.1",
    "fastparquet": "0.2.1",
    "gcsfs": "0.2.2",
    "lxml.etree": "3.8.0",
    "matplotlib": "2.2.2",
    "numexpr": "2.6.2",
    "odfpy": "1.3.0",
    "openpyxl": "2.4.8",
    "pandas_gbq": "0.8.0",
    "pyarrow": "0.9.0",
    "pytables": "3.4.2",
    "s3fs": "0.0.8",
    "scipy": "0.19.0",
    "sqlalchemy": "1.1.4",
    "tables": "3.4.2",
    "xarray": "0.8.2",
    "xlrd": "1.1.0",
    "xlwt": "1.2.0",
    "xlsxwriter": "0.9.8",
}

message = (
    "Missing optional dependency '{name}'. {extra} "
    "Use pip or conda to install {name}."
)
version_message = (
    "Pandas requires version '{minimum_version}' or newer of '{name}' "
    "(version '{actual_version}' currently installed)."
)


def _get_version(module: types.ModuleType) -> str:
    version = getattr(module, "__version__", None)
    if version is None:
        # xlrd uses a capitalized attribute name
        version = getattr(module, "__VERSION__", None)

    if version is None:
        raise ImportError("Can't determine version for {}".format(module.__name__))
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
    try:
        module = importlib.import_module(name)
    except ImportError:
        if raise_on_missing:
            raise ImportError(message.format(name=name, extra=extra)) from None
        else:
            return None

    minimum_version = VERSIONS.get(name)
    if minimum_version:
        version = _get_version(module)
        if distutils.version.LooseVersion(version) < minimum_version:
            assert on_version in {"warn", "raise", "ignore"}
            msg = version_message.format(
                minimum_version=minimum_version, name=name, actual_version=version
            )
            if on_version == "warn":
                warnings.warn(msg, UserWarning)
                return None
            elif on_version == "raise":
                raise ImportError(msg)

    return module
