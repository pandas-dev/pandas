"""
logging helpers, supports vendoring
"""

from __future__ import annotations

import contextlib
import logging
import os
import sys

from typing import IO
from typing import Iterator
from typing import Mapping

log = logging.getLogger(__name__.rsplit(".", 1)[0])
log.propagate = False


class AlwaysStdErrHandler(logging.StreamHandler):  # type: ignore[type-arg]
    def __init__(self) -> None:
        super().__init__(sys.stderr)

    @property
    def stream(self) -> IO[str]:
        return sys.stderr

    @stream.setter
    def stream(self, value: IO[str]) -> None:
        assert value is sys.stderr


def make_default_handler() -> logging.Handler:
    try:
        from rich.console import Console

        console = Console(stderr=True)
        from rich.logging import RichHandler

        return RichHandler(console=console)
    except ImportError:
        last_resort = logging.lastResort
        assert last_resort is not None
        return last_resort


_default_handler = make_default_handler()

log.addHandler(_default_handler)


def _default_log_level(_env: Mapping[str, str] = os.environ) -> int:
    val: str | None = _env.get("SETUPTOOLS_SCM_DEBUG")
    return logging.WARNING if val is None else logging.DEBUG


log.setLevel(_default_log_level())


@contextlib.contextmanager
def defer_to_pytest() -> Iterator[None]:
    log.propagate = True
    old_level = log.level
    log.setLevel(logging.NOTSET)
    log.removeHandler(_default_handler)
    try:
        yield
    finally:
        log.addHandler(_default_handler)
        log.propagate = False
        log.setLevel(old_level)


@contextlib.contextmanager
def enable_debug(handler: logging.Handler = _default_handler) -> Iterator[None]:
    log.addHandler(handler)
    old_level = log.level
    log.setLevel(logging.DEBUG)
    old_handler_level = handler.level
    handler.setLevel(logging.DEBUG)
    try:
        yield
    finally:
        log.setLevel(old_level)
        handler.setLevel(old_handler_level)
        if handler is not _default_handler:
            log.removeHandler(handler)
