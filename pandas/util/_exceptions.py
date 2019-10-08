import contextlib


@contextlib.contextmanager
def rewrite_exception(old_name, new_name):
    """Rewrite the message of an exception."""
    try:
        yield
    except Exception as err:
        msg = err.args[0]
        msg = msg.replace(old_name, new_name)
        args = (msg,)
        if len(err.args) > 1:
            args = args + err.args[1:]
        err.args = args
        raise
