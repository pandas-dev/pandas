import mailbox
from _typeshed import StrPath
from collections.abc import Callable
from typing import IO, Any

DEFAULT_MAXMEM: int

class mbox_readonlydir(mailbox.mbox):
    maxmem: int
    def __init__(
        self,
        path: StrPath,
        factory: Callable[[IO[Any]], mailbox.mboxMessage] | None = None,
        create: bool = True,
        maxmem: int = 1048576,
    ) -> None: ...
    def flush(self) -> None: ...
