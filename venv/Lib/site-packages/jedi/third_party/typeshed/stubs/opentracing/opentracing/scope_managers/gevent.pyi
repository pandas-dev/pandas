from ..scope import Scope
from ..scope_manager import ScopeManager
from ..span import Span

class GeventScopeManager(ScopeManager):
    def activate(self, span: Span, finish_on_close: bool) -> Scope: ...
    @property
    def active(self) -> Scope: ...
