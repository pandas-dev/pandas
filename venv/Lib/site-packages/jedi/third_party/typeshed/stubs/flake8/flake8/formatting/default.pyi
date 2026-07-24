from _typeshed import Incomplete

from ..violation import Violation
from .base import BaseFormatter

COLORS: dict[str, str]
COLORS_OFF: dict[str, str]

class SimpleFormatter(BaseFormatter):
    error_format: str
    def format(self, error: Violation) -> str | None: ...

class Default(SimpleFormatter):
    error_format: str
    def after_init(self) -> None: ...

class Pylint(SimpleFormatter):
    error_format: str

class FilenameOnly(SimpleFormatter):
    error_format: str
    filenames_already_printed: Incomplete
    def after_init(self) -> None: ...
    def show_source(self, error: Violation) -> str | None: ...
    def format(self, error: Violation) -> str | None: ...

class Nothing(BaseFormatter):
    def format(self, error: Violation) -> str | None: ...
    def show_source(self, error: Violation) -> str | None: ...
