from collections.abc import Sequence
from typing import Any

def main() -> None: ...
def get_kwargs(args: Sequence[str]) -> dict[str, Any]: ...  # value depends on the key
