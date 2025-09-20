from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Final

from mypy.checker_shared import TypeCheckerSharedApi

# This is global mutable state. Don't add anything here unless there's a very
# good reason.


class TypeCheckerState:
    # Wrap this in a class since it's faster that using a module-level attribute.

    def __init__(self, type_checker: TypeCheckerSharedApi | None) -> None:
        # Value varies by file being processed
        self.type_checker = type_checker

    @contextmanager
    def set(self, value: TypeCheckerSharedApi) -> Iterator[None]:
        saved = self.type_checker
        self.type_checker = value
        try:
            yield
        finally:
            self.type_checker = saved


checker_state: Final = TypeCheckerState(type_checker=None)
