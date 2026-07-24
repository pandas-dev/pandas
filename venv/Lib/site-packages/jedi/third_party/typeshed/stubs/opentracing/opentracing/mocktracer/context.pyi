from typing_extensions import Self

import opentracing

class SpanContext(opentracing.SpanContext):
    trace_id: int | None
    span_id: int | None
    def __init__(
        self, trace_id: int | None = None, span_id: int | None = None, baggage: dict[str, str] | None = None
    ) -> None: ...
    @property
    def baggage(self) -> dict[str, str]: ...
    def with_baggage_item(self, key: str, value: str) -> Self: ...
