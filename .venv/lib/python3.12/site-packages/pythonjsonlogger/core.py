"""Core functionality shared by all JSON loggers"""

### IMPORTS
### ============================================================================
## Future
from __future__ import annotations

## Standard Library
from datetime import datetime, timezone
import logging
import re
import sys
from typing import Optional, Union, List, Dict, Container, Any, Sequence

if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    from typing_extensions import TypeAlias

## Installed

## Application


### CONSTANTS
### ============================================================================
RESERVED_ATTRS: List[str] = [
    "args",
    "asctime",
    "created",
    "exc_info",
    "exc_text",
    "filename",
    "funcName",
    "levelname",
    "levelno",
    "lineno",
    "module",
    "msecs",
    "message",
    "msg",
    "name",
    "pathname",
    "process",
    "processName",
    "relativeCreated",
    "stack_info",
    "thread",
    "threadName",
]
"""Default reserved attributes.

These come from the [default attributes of `LogRecord` objects](http://docs.python.org/library/logging.html#logrecord-attributes).

Note:
    Although considered a constant, this list is dependent on the Python version due to
    different `LogRecord` objects having different attributes in different Python versions.

*Changed in 3.0*: `RESERVED_ATTRS` is now `list[str]` instead of `tuple[str, ...]`.
"""

if sys.version_info >= (3, 12):
    # taskName added in python 3.12
    RESERVED_ATTRS.append("taskName")
    RESERVED_ATTRS.sort()


STYLE_STRING_TEMPLATE_REGEX = re.compile(r"\$\{(.+?)\}", re.IGNORECASE)  # $ style
STYLE_STRING_FORMAT_REGEX = re.compile(r"\{(.+?)\}", re.IGNORECASE)  # { style
STYLE_PERCENT_REGEX = re.compile(r"%\((.+?)\)", re.IGNORECASE)  # % style

## Type Aliases
## -----------------------------------------------------------------------------
LogData: TypeAlias = Dict[str, Any]
"""Type alias

*Changed in 4.0*: renamed from `LogRecord` to `LogData`
"""


### FUNCTIONS
### ============================================================================
def merge_record_extra(
    record: logging.LogRecord,
    target: Dict,
    reserved: Container[str],
    rename_fields: Optional[Dict[str, str]] = None,
) -> Dict:
    """
    Merges extra attributes from LogRecord object into target dictionary

    Args:
        record: logging.LogRecord
        target: dict to update
        reserved: dict or list with reserved keys to skip
        rename_fields: an optional dict, used to rename field names in the output.
            e.g. Rename `levelname` to `log.level`: `{'levelname': 'log.level'}`

    *Changed in 3.1*: `reserved` is now `Container[str]`.
    """
    if rename_fields is None:
        rename_fields = {}
    for key, value in record.__dict__.items():
        # this allows to have numeric keys
        if key not in reserved and not (hasattr(key, "startswith") and key.startswith("_")):
            target[rename_fields.get(key, key)] = value
    return target


### CLASSES
### ============================================================================
class BaseJsonFormatter(logging.Formatter):
    """Base class for all formatters

    Must not be used directly.

    *New in 3.1*

    *Changed in 3.2*: `defaults` argument is no longer ignored.

    *Added in 3.3*: `exc_info_as_array` and `stack_info_as_array` options are added.
    """

    _style: Union[logging.PercentStyle, str]  # type: ignore[assignment]

    ## Parent Methods
    ## -------------------------------------------------------------------------
    # pylint: disable=too-many-arguments,super-init-not-called
    def __init__(
        self,
        fmt: Optional[Union[str, Sequence[str]]] = None,
        datefmt: Optional[str] = None,
        style: str = "%",
        validate: bool = True,
        *,
        prefix: str = "",
        rename_fields: Optional[Dict[str, str]] = None,
        rename_fields_keep_missing: bool = False,
        static_fields: Optional[Dict[str, Any]] = None,
        reserved_attrs: Optional[Sequence[str]] = None,
        timestamp: Union[bool, str] = False,
        defaults: Optional[Dict[str, Any]] = None,
        exc_info_as_array: bool = False,
        stack_info_as_array: bool = False,
    ) -> None:
        """
        Args:
            fmt: String format or `Sequence` of field names of fields to log.
            datefmt: format to use when formatting `asctime` field
            style: how to extract log fields from `fmt`. Ignored if `fmt` is a `Sequence[str]`.
            validate: validate `fmt` against style, if implementing a custom `style` you
                must set this to `False`. Ignored if `fmt` is a `Sequence[str]`.
            defaults: a dictionary containing default fields that are added before all other fields and
                may be overridden. The supplied fields are still subject to `rename_fields`.
            prefix: an optional string prefix added at the beginning of
                the formatted string
            rename_fields: an optional dict, used to rename field names in the output.
                Rename `message` to `@message`: `{'message': '@message'}`
            rename_fields_keep_missing: When renaming fields, include missing fields in the output.
            static_fields: an optional dict, used to add fields with static values to all logs
            reserved_attrs: an optional list of fields that will be skipped when
                outputting json log record. Defaults to [all log record attributes][pythonjsonlogger.core.RESERVED_ATTRS].
            timestamp: an optional string/boolean field to add a timestamp when
                outputting the json log record. If string is passed, timestamp will be added
                to log record using string as key. If True boolean is passed, timestamp key
                will be "timestamp". Defaults to False/off.
            exc_info_as_array: break the exc_info into a list of lines based on line breaks.
            stack_info_as_array: break the stack_info into a list of lines based on line breaks.

        *Changed in 3.1*:

        - you can now use custom values for style by setting validate to `False`.
          The value is stored in `self._style` as a string. The `parse` method will need to be
          overridden in order to support the new style.
        - Renaming fields now preserves the order that fields were added in and avoids adding
          missing fields. The original behaviour, missing fields have a value of `None`, is still
          available by setting `rename_fields_keep_missing` to `True`.

        *Added in 4.0*:

        - `fmt` now supports comma seperated lists (`style=","`). Note that this style is specific
          to `python-json-logger` and thus care should be taken to not to pass this format to other
          logging Formatter implementations.
        - `fmt` now supports sequences of strings (e.g. lists and tuples) of field names.
        """
        ## logging.Formatter compatibility
        ## ---------------------------------------------------------------------
        # Note: validate added in python 3.8, defaults added in 3.10
        if fmt is None or isinstance(fmt, str):
            if style in logging._STYLES:
                _style = logging._STYLES[style][0](fmt)  # type: ignore[operator]
                if validate:
                    _style.validate()
                self._style = _style
                self._fmt = _style._fmt

            elif style == "," or not validate:
                self._style = style
                self._fmt = fmt
                # TODO: Validate comma format

            else:
                raise ValueError("Style must be one of: '%{$,'")

            self._required_fields = self.parse()

        # Note: we do this check second as string is still a Sequence[str]
        elif isinstance(fmt, Sequence):
            self._style = "__sequence__"
            self._fmt = str(fmt)
            self._required_fields = list(fmt)

        self.datefmt = datefmt

        ## JSON Logging specific
        ## ---------------------------------------------------------------------
        self.prefix = prefix

        # We recreate the dict in rename_fields and static_fields to support internal/external
        # references which require getting the item to do the conversion.
        # For more details see: https://github.com/nhairs/python-json-logger/pull/45
        self.rename_fields = (
            {key: rename_fields[key] for key in rename_fields} if rename_fields is not None else {}
        )
        self.static_fields = (
            {key: static_fields[key] for key in static_fields} if static_fields is not None else {}
        )

        self.rename_fields_keep_missing = rename_fields_keep_missing
        self.reserved_attrs = set(reserved_attrs if reserved_attrs is not None else RESERVED_ATTRS)
        self.timestamp = timestamp

        self._skip_fields = set(self._required_fields)
        self._skip_fields.update(self.reserved_attrs)
        self.defaults = defaults if defaults is not None else {}
        self.exc_info_as_array = exc_info_as_array
        self.stack_info_as_array = stack_info_as_array
        return

    def format(self, record: logging.LogRecord) -> str:
        """Formats a log record and serializes to json

        Args:
            record: the record to format
        """
        message_dict: Dict[str, Any] = {}
        # TODO: logging.LogRecord.msg and logging.LogRecord.message in typeshed
        #        are always type of str. We shouldn't need to override that.
        if isinstance(record.msg, dict):
            message_dict = record.msg
            record.message = ""
        else:
            record.message = record.getMessage()

        # only format time if needed
        if "asctime" in self._required_fields:
            record.asctime = self.formatTime(record, self.datefmt)

        # Display formatted exception, but allow overriding it in the
        # user-supplied dict.
        if record.exc_info and not message_dict.get("exc_info"):
            message_dict["exc_info"] = self.formatException(record.exc_info)
        if not message_dict.get("exc_info") and record.exc_text:
            message_dict["exc_info"] = record.exc_text

        # Display formatted record of stack frames
        # default format is a string returned from :func:`traceback.print_stack`
        if record.stack_info and not message_dict.get("stack_info"):
            message_dict["stack_info"] = self.formatStack(record.stack_info)

        log_data: LogData = {}
        self.add_fields(log_data, record, message_dict)
        log_data = self.process_log_record(log_data)

        return self.serialize_log_record(log_data)

    ## JSON Formatter Specific Methods
    ## -------------------------------------------------------------------------
    def parse(self) -> List[str]:
        """Parses format string looking for substitutions

        This method is responsible for returning a list of fields (as strings)
        to include in all log messages.

        You can support custom styles by overriding this method.

        Returns:
            list of fields to be extracted and serialized
        """
        if self._fmt is None:
            return []

        if isinstance(self._style, str):
            if self._style == "__sequence__":
                raise RuntimeError("Must not call parse when fmt is a sequence of strings")

            if self._style == ",":
                return [field.strip() for field in self._fmt.split(",") if field.strip()]

            raise ValueError(f"Style {self._style!r} is not supported")

        if isinstance(self._style, logging.StringTemplateStyle):
            formatter_style_pattern = STYLE_STRING_TEMPLATE_REGEX

        elif isinstance(self._style, logging.StrFormatStyle):
            formatter_style_pattern = STYLE_STRING_FORMAT_REGEX

        elif isinstance(self._style, logging.PercentStyle):
            # PercentStyle is parent class of StringTemplateStyle and StrFormatStyle
            # so it must be checked last.
            formatter_style_pattern = STYLE_PERCENT_REGEX

        else:
            raise ValueError(f"Style {self._style!r} is not supported")

        return formatter_style_pattern.findall(self._fmt)

    def serialize_log_record(self, log_data: LogData) -> str:
        """Returns the final representation of the data to be logged

        Args:
            log_data: the data

        *Changed in 4.0*: `log_record` renamed to `log_data`
        """
        return self.prefix + self.jsonify_log_record(log_data)

    def add_fields(
        self,
        log_data: Dict[str, Any],
        record: logging.LogRecord,
        message_dict: Dict[str, Any],
    ) -> None:
        """Extract fields from a LogRecord for logging

        This method can be overridden to implement custom logic for adding fields.

        Args:
            log_data: data that will be logged
            record: the record to extract data from
            message_dict: dictionary that was logged instead of a message. e.g
                `logger.info({"is_this_message_dict": True})`

        *Changed in 4.0*: `log_record` renamed to `log_data`
        """
        for field in self.defaults:
            log_data[self._get_rename(field)] = self.defaults[field]

        for field in self._required_fields:
            log_data[self._get_rename(field)] = record.__dict__.get(field)

        for data_dict in [self.static_fields, message_dict]:
            for key, value in data_dict.items():
                log_data[self._get_rename(key)] = value

        merge_record_extra(
            record,
            log_data,
            reserved=self._skip_fields,
            rename_fields=self.rename_fields,
        )

        if self.timestamp:
            key = self.timestamp if isinstance(self.timestamp, str) else "timestamp"
            log_data[self._get_rename(key)] = datetime.fromtimestamp(
                record.created, tz=timezone.utc
            )

        if self.rename_fields_keep_missing:
            for field in self.rename_fields.values():
                if field not in log_data:
                    log_data[field] = None
        return

    def _get_rename(self, key: str) -> str:
        return self.rename_fields.get(key, key)

    # Child Methods
    # ..........................................................................
    def jsonify_log_record(self, log_data: LogData) -> str:
        """Convert the log data into a JSON string.

        Child classes MUST override this method.

        Args:
            log_data: the data to serialize

        *Changed in 4.0*: `log_record` renamed to `log_data`
        """
        raise NotImplementedError()

    def process_log_record(self, log_data: LogData) -> LogData:
        """Custom processing of the data to be logged.

        Child classes can override this method to alter the log record before it
        is serialized.

        Args:
            log_data: incoming data

        *Changed in 4.0*: `log_record` renamed to `log_data`
        """
        return log_data

    def formatException(self, ei) -> Union[str, list[str]]:  # type: ignore
        """Format and return the specified exception information.

        If exc_info_as_array is set to True, This method returns an array of strings.
        """
        exception_info_str = super().formatException(ei)
        return exception_info_str.splitlines() if self.exc_info_as_array else exception_info_str

    def formatStack(self, stack_info) -> Union[str, list[str]]:  # type: ignore
        """Format and return the specified stack information.

        If stack_info_as_array is set to True, This method returns an array of strings.
        """
        stack_info_str = super().formatStack(stack_info)
        return stack_info_str.splitlines() if self.stack_info_as_array else stack_info_str
