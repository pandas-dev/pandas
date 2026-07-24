from typing import Any, NamedTuple

from .scope import Scope
from .scope_manager import ScopeManager
from .span import Span, SpanContext

class Tracer:
    def __init__(self, scope_manager: ScopeManager | None = None) -> None: ...
    @property
    def scope_manager(self) -> ScopeManager: ...
    @property
    def active_span(self) -> Span | None: ...
    def start_active_span(
        self,
        operation_name: str,
        child_of: Span | SpanContext | None = None,
        references: list[Reference] | None = None,
        tags: dict[Any, Any] | None = None,
        start_time: float | None = None,
        ignore_active_span: bool = False,
        finish_on_close: bool = True,
    ) -> Scope: ...
    def start_span(
        self,
        operation_name: str | None = None,
        child_of: Span | SpanContext | None = None,
        references: list[Reference] | None = None,
        tags: dict[Any, Any] | None = None,
        start_time: float | None = None,
        ignore_active_span: bool = False,
    ) -> Span: ...
    def inject(self, span_context: SpanContext, format: str, carrier: dict[Any, Any]) -> None: ...
    def extract(self, format: str, carrier: dict[Any, Any]) -> SpanContext: ...

class ReferenceType:
    CHILD_OF: str
    FOLLOWS_FROM: str

class Reference(NamedTuple):
    type: str
    referenced_context: SpanContext | None

def child_of(referenced_context: SpanContext | None = None) -> Reference: ...
def follows_from(referenced_context: SpanContext | None = None) -> Reference: ...
def start_child_span(
    parent_span: Span, operation_name: str, tags: dict[Any, Any] | None = None, start_time: float | None = None
) -> Span: ...
