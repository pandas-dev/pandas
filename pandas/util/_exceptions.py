from __future__ import annotations

import contextlib
from typing import Iterator


@contextlib.contextmanager
def rewrite_exception(old_name: str, new_name: str) -> Iterator[None]:
    """
    Rewrite the message of an exception.
    """
    try:
        yield
    except Exception as err:
        if not err.args:
            raise
        msg = str(err.args[0])
        msg = msg.replace(old_name, new_name)
        args: tuple[str, ...] = (msg,)
        if len(err.args) > 1:
            args = args + err.args[1:]
        err.args = args
        raise
