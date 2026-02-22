# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root
# for license information.

import socket
import sys
import threading

from debugpy.common import log
from debugpy.common.util import hide_thread_from_debugger

def can_bind_ipv4_localhost():
    """Check if we can bind to IPv4 localhost."""
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        # Try to bind to IPv4 localhost on port 0 (any available port)
        sock.bind(("127.0.0.1", 0))
        sock.close()
        return True
    except (socket.error, OSError, AttributeError):
        return False

def can_bind_ipv6_localhost():
    """Check if we can bind to IPv6 localhost."""
    try:
        sock = socket.socket(socket.AF_INET6, socket.SOCK_STREAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        # Try to bind to IPv6 localhost on port 0 (any available port)
        sock.bind(("::1", 0))
        sock.close()
        return True
    except (socket.error, OSError, AttributeError):
        return False

def get_default_localhost():
    """Get the default localhost address.
    Defaults to IPv4 '127.0.0.1', but falls back to IPv6 '::1' if IPv4 is unavailable.
    """
    # First try IPv4 (preferred default)
    if can_bind_ipv4_localhost():
        return "127.0.0.1"

    # Fall back to IPv6 if IPv4 is not available
    if can_bind_ipv6_localhost():
        return "::1"
    
    # If neither works, still return IPv4 as a last resort
    # (this is a very unusual situation)
    return "127.0.0.1"

def get_address(sock):
    """Gets the socket address host and port."""
    try:
        host, port = sock.getsockname()[:2]
    except Exception as exc:
        log.swallow_exception("Failed to get socket address:")
        raise RuntimeError(f"Failed to get socket address: {exc}") from exc

    return host, port

def create_server(host, port=0, backlog=socket.SOMAXCONN, timeout=None):
    """Return a local server socket listening on the given port."""

    assert backlog > 0
    if host is None:
        host = get_default_localhost()
    if port is None:
        port = 0
    ipv6 = host.count(":") > 1

    try:
        server = _new_sock(ipv6)
        if port != 0:
            # If binding to a specific port, make sure that the user doesn't have
            # to wait until the OS times out the socket to be able to use that port
            # again.if the server or the adapter crash or are force-killed.
            if sys.platform == "win32":
                server.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
            else:
                try:
                    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                except (AttributeError, OSError):  # pragma: no cover
                    pass  # Not available everywhere
        server.bind((host, port))
        if timeout is not None:
            server.settimeout(timeout)
        server.listen(backlog)
    except Exception:  # pragma: no cover
        server.close()
        raise
    return server


def create_client(ipv6=False):
    """Return a client socket that may be connected to a remote address."""
    return _new_sock(ipv6)


def _new_sock(ipv6=False):
    address_family = socket.AF_INET6 if ipv6 else socket.AF_INET
    sock = socket.socket(address_family, socket.SOCK_STREAM, socket.IPPROTO_TCP)

    # Set TCP keepalive on an open socket.
    # It activates after 1 second (TCP_KEEPIDLE,) of idleness,
    # then sends a keepalive ping once every 3 seconds (TCP_KEEPINTVL),
    # and closes the connection after 5 failed ping (TCP_KEEPCNT), or 15 seconds
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    except (AttributeError, OSError):  # pragma: no cover
        pass  # May not be available everywhere.
    try:
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 1)
    except (AttributeError, OSError):  # pragma: no cover
        pass  # May not be available everywhere.
    try:
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 3)
    except (AttributeError, OSError):  # pragma: no cover
        pass  # May not be available everywhere.
    try:
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPCNT, 5)
    except (AttributeError, OSError):  # pragma: no cover
        pass  # May not be available everywhere.
    return sock


def shut_down(sock, how=socket.SHUT_RDWR):
    """Shut down the given socket."""
    sock.shutdown(how)


def close_socket(sock):
    """Shutdown and close the socket."""
    try:
        shut_down(sock)
    except Exception:  # pragma: no cover
        pass
    sock.close()


def serve(name, handler, host, port=0, backlog=socket.SOMAXCONN, timeout=None):
    """Accepts TCP connections on the specified host and port, and invokes the
    provided handler function for every new connection.

    Returns the created server socket.
    """

    assert backlog > 0

    try:
        listener = create_server(host, port, backlog, timeout)
    except Exception:  # pragma: no cover
        log.reraise_exception(
            "Error listening for incoming {0} connections on {1}:{2}:", name, host, port
        )
    host, port = get_address(listener)
    log.info("Listening for incoming {0} connections on {1}:{2}...", name, host, port)

    def accept_worker():
        while True:
            try:
                sock, address = listener.accept()
                other_host, other_port = address[:2]
            except (OSError, socket.error):
                # Listener socket has been closed.
                break

            log.info(
                "Accepted incoming {0} connection from {1}:{2}.",
                name,
                other_host,
                other_port,
            )
            handler(sock)

    thread = threading.Thread(target=accept_worker)
    thread.daemon = True
    hide_thread_from_debugger(thread)
    thread.start()

    return listener
