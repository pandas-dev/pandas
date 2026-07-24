from _typeshed import Incomplete

class Chunk: ...

class TagChunk(Chunk):
    __slots__ = ("tag", "label")
    tag: Incomplete
    label: Incomplete
    def __init__(self, tag: str, label: str | None = None) -> None: ...

class TextChunk(Chunk):
    __slots__ = "text"
    text: Incomplete
    def __init__(self, text: str) -> None: ...
