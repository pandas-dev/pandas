import sys
from typing import Final, LiteralString

# mypy<=1.20 workaround, see https://github.com/python/mypy/pull/20392
if sys.version_info >= (3, 14):
    __conditional_annotations__: Final[set[int]] = ...

docdict: Final[dict[str, str]] = ...  # undocumented

def get(name: LiteralString) -> str: ...  # undocumented
def add_newdoc(name: LiteralString, doc: str) -> None: ...  # undocumented
