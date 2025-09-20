"""Fixtures for use with jupyter events."""
from __future__ import annotations

import io
import json
import logging
from typing import Any, Callable

import pytest

from jupyter_events import EventLogger


@pytest.fixture
def jp_event_sink() -> io.StringIO:
    """A stream for capture events."""
    return io.StringIO()


@pytest.fixture
def jp_event_handler(jp_event_sink: io.StringIO) -> logging.Handler:
    """A logging handler that captures any events emitted by the event handler"""
    return logging.StreamHandler(jp_event_sink)


@pytest.fixture
def jp_read_emitted_events(
    jp_event_handler: logging.Handler, jp_event_sink: io.StringIO
) -> Callable[..., list[str] | None]:
    """Reads list of events since last time it was called."""

    def _read() -> list[str] | None:
        jp_event_handler.flush()
        event_buf = jp_event_sink.getvalue().strip()
        output = [json.loads(item) for item in event_buf.split("\n")] if event_buf else None
        # Clear the sink.
        jp_event_sink.truncate(0)
        jp_event_sink.seek(0)
        return output

    return _read


@pytest.fixture
def jp_event_schemas() -> list[Any]:
    """A list of schema references.

    Each item should be one of the following:
    - string of serialized JSON/YAML content representing a schema
    - a pathlib.Path object pointing to a schema file on disk
    - a dictionary with the schema data.
    """
    return []


@pytest.fixture
def jp_event_logger(jp_event_handler: logging.Handler, jp_event_schemas: list[Any]) -> EventLogger:
    """A pre-configured event logger for tests."""
    logger = EventLogger()
    for schema in jp_event_schemas:
        logger.register_event_schema(schema)
    logger.register_handler(handler=jp_event_handler)
    return logger
