from collections.abc import Callable
from typing import Final
from typing_extensions import LiteralString, Unpack

from .watch import Watch

__version__: Final[LiteralString]

watch: Watch
unwatch: Final[Callable[[Unpack[tuple[object, ...]]], None]]
