"""
pandas._config is considered explicitly upstream of everything else in pandas,
should have no intra-pandas dependencies.
"""
__all__ = ["config", "get_option", "set_option", "reset_option",
           "describe_option", "option_context", "options"]
from pandas._config import config
from pandas._config.config import (
    describe_option, get_option, option_context, options, reset_option,
    set_option)
