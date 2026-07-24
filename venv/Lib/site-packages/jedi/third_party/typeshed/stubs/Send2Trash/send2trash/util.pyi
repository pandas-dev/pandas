from _typeshed import StrOrBytesPath
from typing import Any

# Should be consistent with `__init__.py`
def preprocess_paths(paths: list[Any] | StrOrBytesPath) -> list[str | bytes]: ...
