"""Defines the different custom formats in which mypy can output."""

import json
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mypy.errors import MypyError


class ErrorFormatter(ABC):
    """Base class to define how errors are formatted before being printed."""

    @abstractmethod
    def report_error(self, error: "MypyError") -> str:
        raise NotImplementedError


class JSONFormatter(ErrorFormatter):
    """Formatter for basic JSON output format."""

    def report_error(self, error: "MypyError") -> str:
        """Prints out the errors as simple, static JSON lines."""
        return json.dumps(
            {
                "file": error.file_path,
                "line": error.line,
                "column": error.column,
                "message": error.message,
                "hint": None if len(error.hints) == 0 else "\n".join(error.hints),
                "code": None if error.errorcode is None else error.errorcode.code,
                "severity": error.severity,
            }
        )


OUTPUT_CHOICES = {"json": JSONFormatter()}
