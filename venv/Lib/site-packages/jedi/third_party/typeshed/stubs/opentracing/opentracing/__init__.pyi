from .propagation import (
    Format as Format,
    InvalidCarrierException as InvalidCarrierException,
    SpanContextCorruptedException as SpanContextCorruptedException,
    UnsupportedFormatException as UnsupportedFormatException,
)
from .scope import Scope as Scope
from .scope_manager import ScopeManager as ScopeManager
from .span import Span as Span, SpanContext as SpanContext
from .tracer import (
    Reference as Reference,
    ReferenceType as ReferenceType,
    Tracer as Tracer,
    child_of as child_of,
    follows_from as follows_from,
    start_child_span as start_child_span,
)

tracer: Tracer
is_tracer_registered: bool

def global_tracer() -> Tracer: ...
def set_global_tracer(value: Tracer) -> None: ...
def is_global_tracer_registered() -> bool: ...
