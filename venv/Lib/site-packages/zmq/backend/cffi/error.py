"""zmq error functions"""

# Copyright (C) PyZMQ Developers
# Distributed under the terms of the Modified BSD License.

from ._cffi import ffi
from ._cffi import lib as C


def strerror(errno):
    return ffi.string(C.zmq_strerror(errno)).decode()


zmq_errno = C.zmq_errno

__all__ = ['strerror', 'zmq_errno']
