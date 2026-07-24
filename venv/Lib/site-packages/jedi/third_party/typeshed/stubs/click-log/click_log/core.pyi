import logging

LOGGER_KEY: str
DEFAULT_LEVEL: int
PY2: bool
text_type: type

class ColorFormatter(logging.Formatter):
    colors: dict[str, dict[str, str]]
    def format(self, record: logging.LogRecord) -> str: ...

class ClickHandler(logging.Handler):
    def emit(self, record: logging.LogRecord) -> None: ...

def basic_config(logger: logging.Logger | str | None = None) -> None: ...
