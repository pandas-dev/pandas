from typing_extensions import Self

__tracebackhide__: bool

class SnapshotMixin:
    def snapshot(self, id: str | None = None, path: str = "__snapshots") -> Self: ...
