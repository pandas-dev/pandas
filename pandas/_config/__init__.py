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
from pandas._config.config import _global_config
from pandas._config.config import describe_option
from pandas._config.config import get_option
from pandas._config.config import option_context
from pandas._config.config import options
from pandas._config.config import reset_option
from pandas._config.config import set_option
from pandas._config.display import detect_console_encoding


def using_copy_on_write():
    _mode_options = _global_config["mode"]
    return _mode_options["copy_on_write"] and _mode_options["data_manager"] == "block"


def using_nullable_dtypes():
    _mode_options = _global_config["mode"]
    return _mode_options["nullable_dtypes"]
