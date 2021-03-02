import contextlib
import inspect
from typing import Tuple


@contextlib.contextmanager
def rewrite_exception(old_name: str, new_name: str):
    """
    Rewrite the message of an exception.
    """
    try:
        yield
    except Exception as err:
        msg = str(err.args[0])
        msg = msg.replace(old_name, new_name)
        args: Tuple[str, ...] = (msg,)
        if len(err.args) > 1:
            args = args + err.args[1:]
        err.args = args
        raise


def find_stack_level() -> int:
    """
    Find the appropriate stacklevel with which to issue a warning for astype.
    """
    stack = inspect.stack()

    # find the lowest-level "astype" call that got us here
    for n in range(2, 6):
        if stack[n].function == "astype":
            break

    while stack[n].function in ["astype", "apply", "astype_array_safe", "astype_array"]:
        # e.g.
        #  bump up Block.astype -> BlockManager.astype -> NDFrame.astype
        #  bump up Datetime.Array.astype -> DatetimeIndex.astype
        n += 1

    if stack[n].function == "__init__":
        # Series.__init__
        n += 1

    return n
