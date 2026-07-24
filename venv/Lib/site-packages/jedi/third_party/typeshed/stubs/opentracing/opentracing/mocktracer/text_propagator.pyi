from typing import Any

from .context import SpanContext
from .propagator import Propagator

prefix_tracer_state: str
prefix_baggage: str
field_name_trace_id: str
field_name_span_id: str
field_count: int

class TextPropagator(Propagator):
    def inject(self, span_context: SpanContext, carrier: dict[Any, Any]) -> None: ...
    def extract(self, carrier: dict[Any, Any]) -> SpanContext: ...
