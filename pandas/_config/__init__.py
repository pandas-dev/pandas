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
]
from . import config
from . import dates  # noqa:F401
from .config import (
    describe_option,
    get_option,
    option_context,
    options,
    reset_option,
    set_option,
)
from .display import detect_console_encoding
