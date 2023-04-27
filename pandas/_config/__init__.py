"""
pandas._config is considered explicitly upstream of everything else in pandas,
should have no intra-pandas dependencies.

importing `dates` and `display` ensures that keys needed by _libs
are initialized.
"""
__all__ = [
    "config",
    "detect_console_encoding",
    "get_option",
    "set_option",
    "reset_option",
    "describe_option",
    "option_context",
    "options",
    "using_copy_on_write",
]
from pandas._config import config
from pandas._config import dates  # pyright: ignore # noqa:F401
from pandas._config.config import (
    _global_config,
    describe_option,
    get_option,
    option_context,
    options,
    reset_option,
    set_option,
)
from pandas._config.display import detect_console_encoding


def using_copy_on_write() -> bool:
    """
    Determine if copy-on-write mode is being used.

    This function examines the global pandas configuration settings.
    Copy-on-write is a memory-saving technique that avoids copying data
    until it is actually modified, which can be particularly useful when
    working with large datasets.

    Parameters
    ----------
        None

    Returns
    -------
    bool
        'True' if copy-on-write mode is enabled and the data manager
        is set to "block", otherwise 'False'.

    Example
    -------
    >>> pd.set_option('mode.chained_assignment', 'raise')
    >>> pd.set_option('mode.copy', True)
    >>> pd.set_option('mode.use_inf_as_na', True)
    >>> pd.set_option('compute.use_bottleneck', False)
    >>> pd.set_option('compute.use_numexpr', False)
    >>> using_copy_on_write()
    True
    """
    _mode_options = _global_config["mode"]
    return _mode_options["copy_on_write"] and _mode_options["data_manager"] == "block"


def using_nullable_dtypes() -> bool:
    """
    Determine if nullable data types are being used.

    This function reads the `mode` configuration option from the global
    configuration object and indicate whether nullable data types are allowed.

    Parameters
    ----------
        None

    Returns
    -------
    bool
        'True' if pandas is using nullable data types, otherwise 'False'.

    Example
    -------
    >>> pd.set_option('mode.use_inf_as_na', True)
    >>> pd.set_option('mode.use_nullable_dtypes', True)
    >>> using_nullable_dtypes()
    True
    >>> pd.set_option('mode.use_nullable_dtypes', False)
    >>> using_nullable_dtypes()
    False

    Notes
    -----
    This function assumes that the global pandas configuration settings have
    already been initialized. If you are not sure whether the settings
    have been initialized, you can call the
    pandas.api.types.is_extension_type() function to force initialization.
    """
    _mode_options = _global_config["mode"]
    return _mode_options["nullable_dtypes"]
