"""Pair of ZMQ inproc sockets used for communication between threads."""

from __future__ import annotations

from typing import Any

import zmq
from tornado.ioloop import IOLoop
from zmq.eventloop.zmqstream import ZMQStream


class SocketPair:
    """Pair of ZMQ inproc sockets for one-direction communication between 2 threads.

    One of the threads is always the shell_channel_thread, the other may be the control
    thread, main thread or a subshell thread.

    .. versionadded:: 7
    """

    from_socket: zmq.Socket[Any]
    to_socket: zmq.Socket[Any]
    to_stream: ZMQStream | None = None

    def __init__(self, context: zmq.Context[Any], name: str):
        """Initialize the inproc socker pair."""
        self.from_socket = context.socket(zmq.PAIR)
        self.to_socket = context.socket(zmq.PAIR)
        address = self._address(name)
        self.from_socket.bind(address)
        self.to_socket.connect(address)  # Or do I need to do this in another thread?

    def close(self):
        """Close the inproc socker pair."""
        self.from_socket.close()

        if self.to_stream is not None:
            self.to_stream.close()
        self.to_socket.close()

    def on_recv(self, io_loop: IOLoop, on_recv_callback, copy: bool = False):
        """Set the callback used when a message is received on the to stream."""
        # io_loop is that of the 'to' thread.
        if self.to_stream is None:
            self.to_stream = ZMQStream(self.to_socket, io_loop)
        self.to_stream.on_recv(on_recv_callback, copy=copy)

    def _address(self, name) -> str:
        """Return the address used for this inproc socket pair."""
        return f"inproc://subshell{name}"
