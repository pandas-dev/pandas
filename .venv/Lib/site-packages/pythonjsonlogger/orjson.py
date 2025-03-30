"""JSON Formatter using [orjson](https://github.com/ijl/orjson)"""

### IMPORTS
### ============================================================================
## Future
from __future__ import annotations

## Standard Library
from typing import Any

## Installed

## Application
from . import core
from . import defaults as d
from .utils import package_is_available

# We import msgspec after checking it is available
package_is_available("orjson", throw_error=True)
import orjson  # pylint: disable=wrong-import-position,wrong-import-order


### FUNCTIONS
### ============================================================================
def orjson_default(obj: Any) -> Any:
    """orjson default encoder function for non-standard types"""
    if d.use_exception_default(obj):
        return d.exception_default(obj)
    if d.use_traceback_default(obj):
        return d.traceback_default(obj)
    if d.use_bytes_default(obj):
        return d.bytes_default(obj)
    if d.use_enum_default(obj):
        return d.enum_default(obj)
    if d.use_type_default(obj):
        return d.type_default(obj)
    return d.unknown_default(obj)


### CLASSES
### ============================================================================
class OrjsonFormatter(core.BaseJsonFormatter):
    """JSON formatter using [orjson](https://github.com/ijl/orjson) for encoding."""

    def __init__(
        self,
        *args,
        json_default: core.OptionalCallableOrStr = orjson_default,
        json_indent: bool = False,
        **kwargs,
    ) -> None:
        """
        Args:
            args: see [BaseJsonFormatter][pythonjsonlogger.core.BaseJsonFormatter]
            json_default: a function for encoding non-standard objects
            json_indent: indent output with 2 spaces.
            kwargs: see [BaseJsonFormatter][pythonjsonlogger.core.BaseJsonFormatter]
        """
        super().__init__(*args, **kwargs)

        self.json_default = core.str_to_object(json_default)
        self.json_indent = json_indent
        return

    def jsonify_log_record(self, log_record: core.LogRecord) -> str:
        """Returns a json string of the log record."""
        opt = orjson.OPT_NON_STR_KEYS
        if self.json_indent:
            opt |= orjson.OPT_INDENT_2

        return orjson.dumps(log_record, default=self.json_default, option=opt).decode("utf8")
