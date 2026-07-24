from _typeshed import Incomplete

from .shapes import Path, UserNode

class SvgPath(Path, UserNode):
    fillColor: Incomplete
    def __init__(self, s, isClipPath: int = 0, autoclose=None, fillMode=0, **kw) -> None: ...
    def provideNode(self): ...

__all__ = ("SvgPath",)
