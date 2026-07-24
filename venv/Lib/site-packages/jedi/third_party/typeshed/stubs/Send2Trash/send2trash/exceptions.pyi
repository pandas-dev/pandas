from typing import Any

class TrashPermissionError(PermissionError):
    # Typed the same as `filename` in `PermissionError`:
    def __init__(self, filename: Any) -> None: ...
