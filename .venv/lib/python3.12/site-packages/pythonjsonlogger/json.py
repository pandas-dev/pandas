"""JSON formatter using the standard library's `json` for encoding.

Module contains the `JsonFormatter` and a custom `JsonEncoder` which supports a greater
variety of types.
"""

### IMPORTS
### ============================================================================
## Future
from __future__ import annotations

## Standard Library
import datetime
import json
from typing import Any, Callable, Optional, Union
import warnings

## Application
from . import core
from . import defaults as d


### CLASSES
### ============================================================================
class JsonEncoder(json.JSONEncoder):
    """A custom encoder extending [json.JSONEncoder](https://docs.python.org/3/library/json.html#json.JSONEncoder)"""

    def default(self, o: Any) -> Any:
        if d.use_datetime_any(o):
            return self.format_datetime_obj(o)

        if d.use_exception_default(o):
            return d.exception_default(o)

        if d.use_traceback_default(o):
            return d.traceback_default(o)

        if d.use_enum_default(o):
            return d.enum_default(o)

        if d.use_bytes_default(o):
            return d.bytes_default(o)

        if d.use_dataclass_default(o):
            return d.dataclass_default(o)

        if d.use_type_default(o):
            return d.type_default(o)

        try:
            return super().default(o)
        except TypeError:
            return d.unknown_default(o)

    def format_datetime_obj(self, o: datetime.time | datetime.date | datetime.datetime) -> str:
        """Format datetime objects found in `self.default`

        This allows subclasses to change the datetime format without understanding the
        internals of the default method.
        """
        return d.datetime_any(o)


class JsonFormatter(core.BaseJsonFormatter):
    """JSON formatter using the standard library's [`json`](https://docs.python.org/3/library/json.html) for encoding"""

    def __init__(
        self,
        *args,
        json_default: Optional[Callable] = None,
        json_encoder: Optional[Callable] = None,
        json_serializer: Callable = json.dumps,
        json_indent: Optional[Union[int, str]] = None,
        json_ensure_ascii: bool = True,
        **kwargs,
    ) -> None:
        """
        Args:
            args: see [BaseJsonFormatter][pythonjsonlogger.core.BaseJsonFormatter]
            json_default: a function for encoding non-standard objects
            json_encoder: custom JSON encoder
            json_serializer: a [`json.dumps`](https://docs.python.org/3/library/json.html#json.dumps)-compatible callable
                that will be used to serialize the log record.
            json_indent: indent parameter for the `json_serializer`
            json_ensure_ascii: `ensure_ascii` parameter for the `json_serializer`
            kwargs: see [BaseJsonFormatter][pythonjsonlogger.core.BaseJsonFormatter]
        """
        super().__init__(*args, **kwargs)

        self.json_default = json_default
        self.json_encoder = json_encoder
        self.json_serializer = json_serializer
        self.json_indent = json_indent
        self.json_ensure_ascii = json_ensure_ascii
        if not self.json_encoder and not self.json_default:
            self.json_encoder = JsonEncoder
        return

    def jsonify_log_record(self, log_data: core.LogData) -> str:
        """Returns a json string of the log data."""
        return self.json_serializer(
            log_data,
            default=self.json_default,
            cls=self.json_encoder,
            indent=self.json_indent,
            ensure_ascii=self.json_ensure_ascii,
        )


### DEPRECATED COMPATIBILITY
### ============================================================================
def __getattr__(name: str):
    if name == "RESERVED_ATTRS":
        warnings.warn(
            "RESERVED_ATTRS has been moved to pythonjsonlogger.core",
            DeprecationWarning,
        )
        return core.RESERVED_ATTRS
    raise AttributeError(f"module {__name__} has no attribute {name}")
