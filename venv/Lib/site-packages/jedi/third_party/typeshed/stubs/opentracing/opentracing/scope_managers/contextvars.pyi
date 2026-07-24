from ..scope import Scope
from ..scope_manager import ScopeManager
from ..span import Span

class ContextVarsScopeManager(ScopeManager):
    def activate(self, span: Span, finish_on_close: bool) -> Scope: ...
    @property
    def active(self) -> Scope: ...

def no_parent_scope() -> None: ...
