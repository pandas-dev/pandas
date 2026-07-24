import argparse
import enum
from _typeshed import Incomplete
from collections.abc import Generator, Sequence

from .formatting.base import BaseFormatter
from .statistics import Statistics

__all__ = ("StyleGuide",)

class Selected(enum.Enum):
    Explicitly = "explicitly selected"
    Implicitly = "implicitly selected"

class Ignored(enum.Enum):
    Explicitly = "explicitly ignored"
    Implicitly = "implicitly ignored"

class Decision(enum.Enum):
    Ignored = "ignored error"
    Selected = "selected error"

class DecisionEngine:
    cache: Incomplete
    selected_explicitly: Incomplete
    ignored_explicitly: Incomplete
    selected: Incomplete
    ignored: Incomplete
    def __init__(self, options: argparse.Namespace) -> None: ...
    def was_selected(self, code: str) -> Selected | Ignored: ...
    def was_ignored(self, code: str) -> Selected | Ignored: ...
    def make_decision(self, code: str) -> Decision: ...
    def decision_for(self, code: str) -> Decision: ...

class StyleGuideManager:
    options: Incomplete
    formatter: Incomplete
    stats: Incomplete
    decider: Incomplete
    style_guides: Incomplete
    default_style_guide: Incomplete
    style_guide_for: Incomplete
    def __init__(self, options: argparse.Namespace, formatter: BaseFormatter, decider: DecisionEngine | None = None) -> None: ...
    def populate_style_guides_with(self, options: argparse.Namespace) -> Generator[StyleGuide]: ...
    def processing_file(self, filename: str) -> Generator[StyleGuide]: ...
    def handle_error(
        self, code: str, filename: str, line_number: int, column_number: int, text: str, physical_line: str | None = None
    ) -> int: ...

class StyleGuide:
    options: Incomplete
    formatter: Incomplete
    stats: Incomplete
    decider: Incomplete
    filename: Incomplete
    def __init__(
        self,
        options: argparse.Namespace,
        formatter: BaseFormatter,
        stats: Statistics,
        filename: str | None = None,
        decider: DecisionEngine | None = None,
    ) -> None: ...
    def copy(self, filename: str | None = None, extend_ignore_with: Sequence[str] | None = None) -> StyleGuide: ...
    def processing_file(self, filename: str) -> Generator[StyleGuide]: ...
    def applies_to(self, filename: str) -> bool: ...
    def should_report_error(self, code: str) -> Decision: ...
    def handle_error(
        self, code: str, filename: str, line_number: int, column_number: int, text: str, physical_line: str | None = None
    ) -> int: ...
