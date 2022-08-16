from __future__ import annotations

import contextlib
import functools
import inspect
import os
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


@functools.lru_cache
def find_stack_level(frame) -> int:
    """
    Find the first place in the stack that is not inside pandas
    (tests notwithstanding).

    ``frame`` should be passed as ``inspect.currentframe()`` by the
    calling function.

    https://stackoverflow.com/questions/17407119/python-inspect-stack-is-slow
    """

    import pandas as pd

    pkg_dir = os.path.dirname(pd.__file__)
    test_dir = os.path.join(pkg_dir, "tests")

    n = 1
    while frame:
        fname = inspect.getfile(frame)
        if fname.startswith(pkg_dir) and not fname.startswith(test_dir):
            frame = frame.f_back
            n += 1
        else:
            break
    return n
