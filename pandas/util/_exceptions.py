import contextlib
from typing import Tuple


@contextlib.contextmanager
def rewrite_exception(old_name: str, new_name: str):
    """Rewrite the message of an exception."""
    try:
        yield
    except Exception as e:
        msg = e.args[0]
        msg = msg.replace(old_name, new_name)
        args: Tuple[str, ...] = (msg,)
        if len(e.args) > 1:
            args = args + e.args[1:]
        e.args = args
        raise
