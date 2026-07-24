from _typeshed import Incomplete

signals_available: bool

class Namespace:
    def signal(self, name: str, doc: str | None = None) -> _FakeSignal: ...

class _FakeSignal:
    name: str
    __doc__: str | None
    def __init__(self, name: str, doc: str | None = None) -> None: ...
    send: Incomplete
    connect: Incomplete
    disconnect: Incomplete
    has_receivers_for: Incomplete
    receivers_for: Incomplete
    temporarily_connected_to: Incomplete
    connected_to: Incomplete

scope_changed: _FakeSignal
