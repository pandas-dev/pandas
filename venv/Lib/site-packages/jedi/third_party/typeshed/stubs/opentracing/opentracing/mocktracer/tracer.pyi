from typing import Any

from ..scope_manager import ScopeManager
from ..span import Span
from ..tracer import Reference, Tracer
from .context import SpanContext
from .propagator import Propagator
from .span import MockSpan

class MockTracer(Tracer):
    def __init__(self, scope_manager: ScopeManager | None = None) -> None: ...
    @property
    def active_span(self) -> MockSpan | None: ...
    def register_propagator(self, format: str, propagator: Propagator) -> None: ...
    def finished_spans(self) -> list[MockSpan]: ...
    def reset(self) -> None: ...
    def start_span(  # type: ignore[override]
        self,
        operation_name: str | None = None,
        child_of: Span | SpanContext | None = None,
        references: list[Reference] | None = None,
        tags: dict[Any, Any] | None = None,
        start_time: float | None = None,
        ignore_active_span: bool = False,
    ) -> MockSpan: ...
    def extract(self, format: str, carrier: dict[Any, Any]) -> SpanContext: ...
