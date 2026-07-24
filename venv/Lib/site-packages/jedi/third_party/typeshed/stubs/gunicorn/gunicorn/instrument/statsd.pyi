import logging
import socket
from collections.abc import Mapping
from datetime import timedelta
from typing import Final

from gunicorn.config import Config
from gunicorn.glogging import Logger
from gunicorn.http import Request
from gunicorn.http.wsgi import Response

from .._types import _EnvironType
from ..glogging import _LogLevelType

METRIC_VAR: Final = "metric"
VALUE_VAR: Final = "value"
MTYPE_VAR: Final = "mtype"
GAUGE_TYPE: Final = "gauge"
COUNTER_TYPE: Final = "counter"
HISTOGRAM_TYPE: Final = "histogram"
TIMER_TYPE: Final = "timer"

class Statsd(Logger):
    prefix: str
    sock: socket.socket | None
    dogstatsd_tags: str | None
    cfg: Config

    def __init__(self, cfg: Config) -> None: ...
    def critical(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def error(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def warning(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def info(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def debug(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def exception(
        self,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = True,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def log(
        self,
        lvl: _LogLevelType,
        msg: object,
        *args: object,
        exc_info: logging._ExcInfoType = None,
        stack_info: bool = False,
        stacklevel: int = 1,
        extra: Mapping[str, object] | None = None,
    ) -> None: ...
    def access(self, resp: Response, req: Request, environ: _EnvironType, request_time: timedelta) -> None: ...
    def gauge(self, name: str, value: float) -> None: ...
    def increment(self, name: str, value: int, sampling_rate: float = 1.0) -> None: ...
    def decrement(self, name: str, value: int, sampling_rate: float = 1.0) -> None: ...
    def timer(self, name: str, value: float) -> None: ...
    def histogram(self, name: str, value: float) -> None: ...
