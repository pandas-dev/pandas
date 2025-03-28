class Mark:
    name: str
    index: int
    line: int
    column: int
    buffer: str | None
    pointer: int
    def __init__(self, name: str, index: int, line: int, column: int, buffer: str | None, pointer: int) -> None: ...
    def get_snippet(self, indent: int = 4, max_length: int = 75) -> str | None: ...

class YAMLError(Exception): ...

class MarkedYAMLError(YAMLError):
    context: str | None
    context_mark: Mark | None
    problem: str | None
    problem_mark: Mark | None
    note: str | None
    def __init__(
        self,
        context: str | None = None,
        context_mark: Mark | None = None,
        problem: str | None = None,
        problem_mark: Mark | None = None,
        note: str | None = None,
    ) -> None: ...
