import atexit


def _encode_string(s):
    encoded = s.encode('utf-8')
    return encoded


def _decode_string(b):
    return b.decode('utf-8')


_encode_string.__doc__ = """Encode a string for use by LLVM."""
_decode_string.__doc__ = """Decode a LLVM character (byte)string."""


_shutting_down = [False]


def _at_shutdown():
    _shutting_down[0] = True


atexit.register(_at_shutdown)


def _is_shutting_down(_shutting_down=_shutting_down):
    """
    Whether the interpreter is currently shutting down.
    For use in finalizers, __del__ methods, and similar; it is advised
    to early bind this function rather than look it up when calling it,
    since at shutdown module globals may be cleared.
    """
    return _shutting_down[0]
