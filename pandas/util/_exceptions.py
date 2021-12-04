from __future__ import annotations

import contextlib
import inspect
import os


@contextlib.contextmanager
def rewrite_exception(old_name: str, new_name: str):
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


def find_stack_level() -> int:
    """
    Find the first place in the stack that is not inside pandas
    (tests notwithstanding).
    """
    stack = inspect.stack()

    import pandas as pd

    pkg_dir = os.path.dirname(pd.__file__)
    test_dir = os.path.join(pkg_dir, "tests")

    import json

    filename = "/home/richard/pandas/find_stack_level.json"
    if os.path.exists(filename):
        with open(filename) as f:
            calls = json.load(f)
    else:
        calls = {}

    with open("/home/richard/pandas/pytest_hack.txt") as f:
        pytest_id = f.read()

    mystack = inspect.stack()
    test_files = [
        f"{e.filename[len(pkg_dir):]}:{e.lineno}:{pytest_id}"
        for e in mystack
        if e.filename.startswith(test_dir)
    ]
    if len(test_files) == 0:
        tester = "UNKNOWN"
    else:
        tester = test_files[0]
    # caller = f"{stack[1].filename[len(pkg_dir):]}:{stack[1].lineno}"
    entry = calls.get(tester, {})
    entry["calls"] = entry.get("calls", 0) + 1
    calls[tester] = entry

    with open(filename, "w") as f:
        json.dump(calls, f)

    for n in range(len(stack)):
        fname = stack[n].filename
        if fname.startswith(pkg_dir) and not fname.startswith(test_dir):
            continue
        else:
            break
    return n
