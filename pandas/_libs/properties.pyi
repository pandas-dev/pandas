# pyright: reportIncompleteStub = false
from typing import Any

cache_readonly = property

def __getattr__(name: str) -> Any: ...  # incomplete
