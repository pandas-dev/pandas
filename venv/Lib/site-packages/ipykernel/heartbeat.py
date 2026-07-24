"""The client and server for a basic ping-pong style heartbeat."""

# -----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

import errno
import socket
from pathlib import Path
from threading import Thread

import zmq
from jupyter_client.localinterfaces import localhost

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class Heartbeat(Thread):
    """A simple ping-pong style heartbeat that runs in a thread."""

    def __init__(self, context, addr=None, *, curve_publickey=None, curve_secretkey=None):
        """Initialize the heartbeat thread.

        Parameters
        ----------
        context : zmq.Context
        addr : tuple, optional
            (transport, ip, port)
        curve_publickey : bytes, optional
            CurveZMQ public key (Z85). When provided together with
            *curve_secretkey*, the heartbeat socket will operate as a
            CurveZMQ server so that only authenticated clients can connect.
        curve_secretkey : bytes, optional
            CurveZMQ secret key (Z85, paired with *curve_publickey*).
        """
        if addr is None:
            addr = ("tcp", localhost(), 0)
        Thread.__init__(self, name="Heartbeat")
        self.context = context
        self.transport, self.ip, self.port = addr
        self.original_port = self.port
        self._curve_publickey = curve_publickey
        self._curve_secretkey = curve_secretkey
        if self.original_port == 0:
            self.pick_port()
        self.addr = (self.ip, self.port)
        self.daemon = True
        self.pydev_do_not_trace = True
        self.is_pydev_daemon_thread = True
        self.name = "Heartbeat"

    def pick_port(self):
        """Pick a port for the heartbeat."""
        if self.transport == "tcp":
            s = socket.socket()
            # '*' means all interfaces to 0MQ, which is '' to socket.socket
            s.bind(("" if self.ip == "*" else self.ip, 0))
            self.port = s.getsockname()[1]
            s.close()
        elif self.transport == "ipc":
            self.port = 1
            while Path(f"{self.ip}-{self.port}").exists():
                self.port = self.port + 1
        else:
            raise ValueError("Unrecognized zmq transport: %s" % self.transport)
        return self.port

    def _try_bind_socket(self):
        c = ":" if self.transport == "tcp" else "-"
        return self.socket.bind(f"{self.transport}://{self.ip}" + c + str(self.port))

    def _bind_socket(self):
        try:
            win_in_use = errno.WSAEADDRINUSE  # type:ignore[attr-defined]
        except AttributeError:
            win_in_use = None

        # Try up to 100 times to bind a port when in conflict to avoid
        # infinite attempts in bad setups
        max_attempts = 1 if self.original_port else 100
        for attempt in range(max_attempts):
            try:
                self._try_bind_socket()
            except zmq.ZMQError as ze:
                if attempt == max_attempts - 1:
                    raise
                # Raise if we have any error not related to socket binding
                if ze.errno != errno.EADDRINUSE and ze.errno != win_in_use:
                    raise
                # Raise if we have any error not related to socket binding
                if self.original_port == 0:
                    self.pick_port()
                else:
                    raise
            else:
                return

    def run(self):
        """Run the heartbeat thread."""
        self.name = "Heartbeat"
        self.socket = self.context.socket(zmq.ROUTER)
        self.socket.linger = 1000
        if self._curve_secretkey is not None:
            self.socket.curve_secretkey = self._curve_secretkey
            self.socket.curve_publickey = self._curve_publickey
            self.socket.curve_server = True
        try:
            self._bind_socket()
        except Exception:
            self.socket.close()
            raise

        while True:
            try:
                zmq.device(zmq.QUEUE, self.socket, self.socket)  # type:ignore[attr-defined]
            except zmq.ZMQError as e:
                if e.errno == errno.EINTR:
                    # signal interrupt, resume heartbeat
                    continue
                if e.errno == zmq.ETERM:
                    # context terminated, close socket and exit
                    try:
                        self.socket.close()
                    except zmq.ZMQError:
                        # suppress further errors during cleanup
                        # this shouldn't happen, though
                        pass
                    break
                if e.errno == zmq.ENOTSOCK:
                    # socket closed elsewhere, exit
                    break
                raise
            else:
                break

    _curve_publickey: bytes | None
    _curve_secretkey: bytes | None
