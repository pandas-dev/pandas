from .api import (
    MessageFailure as MessageFailure,
    add_message as add_message,
    debug as debug,
    error as error,
    get_level as get_level,
    get_messages as get_messages,
    info as info,
    set_level as set_level,
    success as success,
    warning as warning,
)
from .constants import (
    DEBUG as DEBUG,
    DEFAULT_LEVELS as DEFAULT_LEVELS,
    DEFAULT_TAGS as DEFAULT_TAGS,
    ERROR as ERROR,
    INFO as INFO,
    SUCCESS as SUCCESS,
    WARNING as WARNING,
)

default_app_config: str = ...
