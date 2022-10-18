"""
Utilities used in various `compat` components.
"""

import pickle


def flatten_buffer(b: bytes | bytearray | memoryview | pickle.PickleBuffer):
    """
    Return some 1-D `uint8` typed buffer.

    Coerces anything that does not match that description to one that does
    without copying if possible (otherwise will copy).
    """

    if isinstance(b, (bytes, bytearray)):
        return b

    if not isinstance(b, pickle.PickleBuffer):
        b = pickle.PickleBuffer(b)

    try:
        # coerce to 1-D `uint8` C-contiguous `memoryview` zero-copy
        return b.raw()
    except BufferError:
        # perform in-memory copy if buffer is not contiguous
        return memoryview(b).tobytes()
