import contextlib


@contextlib.contextmanager
def rewrite_exception(old_name, new_name):
    """Rewrite the message of an exception."""
    try:
        yield
    except Exception as e:
        msg = e.args[0]
        msg = msg.replace(old_name, new_name)
        args = (msg,)
        if len(e.args) > 1:
            args = args + e.args[1:]
        e.args = args
        raise
