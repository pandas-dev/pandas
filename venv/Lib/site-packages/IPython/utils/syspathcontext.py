from __future__ import annotations

import sys
from types import TracebackType
from typing import Literal, Self

import warnings


class prepended_to_syspath:
    """A context for prepending a directory to sys.path for a second."""

    dir: str
    added: bool

    def __init__(self, dir: str) -> None:
        self.dir = dir
        self.added = False

    def __enter__(self) -> Self:
        if self.dir not in sys.path:
            sys.path.insert(0, self.dir)
            self.added = True
        else:
            self.added = False
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> Literal[False]:
        if self.added:
            try:
                sys.path.remove(self.dir)
            except ValueError:
                pass
        # Returning False causes any exceptions to be re-raised.
        return False
