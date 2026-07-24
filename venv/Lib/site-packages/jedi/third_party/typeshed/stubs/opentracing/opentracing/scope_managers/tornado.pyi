from typing import Any

from ..scope import Scope
from ..scope_managers import ThreadLocalScopeManager
from ..span import Span

class TornadoScopeManager(ThreadLocalScopeManager):
    def activate(self, span: Span, finish_on_close: bool) -> Scope: ...
    @property
    def active(self) -> Scope: ...

class ThreadSafeStackContext:
    contexts: Any
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

def tracer_stack_context() -> ThreadSafeStackContext: ...
