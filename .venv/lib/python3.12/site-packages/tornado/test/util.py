import contextlib
import os
import platform
import socket
import sys
import sysconfig
import textwrap
import typing
import unittest
import warnings

from tornado.testing import bind_unused_port

_TestCaseType = typing.TypeVar("_TestCaseType", bound=typing.Type[unittest.TestCase])

skipIfNonUnix = unittest.skipIf(
    os.name != "posix" or sys.platform == "cygwin", "non-unix platform"
)

# Set the environment variable NO_NETWORK=1 to disable any tests that
# depend on an external network.
skipIfNoNetwork = unittest.skipIf("NO_NETWORK" in os.environ, "network access disabled")

skipNotCPython = unittest.skipIf(
    # "CPython" here essentially refers to the traditional synchronous refcounting GC,
    # so we skip these tests in free-threading builds of cpython too.
    platform.python_implementation() != "CPython"
    or sysconfig.get_config_var("Py_GIL_DISABLED"),
    "Not CPython implementation",
)


def _detect_ipv6():
    if not socket.has_ipv6:
        # socket.has_ipv6 check reports whether ipv6 was present at compile
        # time. It's usually true even when ipv6 doesn't work for other reasons.
        return False
    sock = None
    try:
        sock = socket.socket(socket.AF_INET6)
        sock.bind(("::1", 0))
    except OSError:
        return False
    finally:
        if sock is not None:
            sock.close()
    return True


skipIfNoIPv6 = unittest.skipIf(not _detect_ipv6(), "ipv6 support not present")


def refusing_port():
    """Returns a local port number that will refuse all connections.

    Return value is (cleanup_func, port); the cleanup function
    must be called to free the port to be reused.
    """
    # On travis-ci port numbers are reassigned frequently. To avoid
    # collisions with other tests, we use an open client-side socket's
    # ephemeral port number to ensure that nothing can listen on that
    # port.
    server_socket, port = bind_unused_port()
    server_socket.setblocking(True)
    client_socket = socket.socket()
    client_socket.connect(("127.0.0.1", port))
    conn, client_addr = server_socket.accept()
    conn.close()
    server_socket.close()
    return (client_socket.close, client_addr[1])


def exec_test(caller_globals, caller_locals, s):
    """Execute ``s`` in a given context and return the result namespace.

    Used to define functions for tests in particular python
    versions that would be syntax errors in older versions.
    """
    # Flatten the real global and local namespace into our fake
    # globals: it's all global from the perspective of code defined
    # in s.
    global_namespace = dict(caller_globals, **caller_locals)  # type: ignore
    local_namespace = {}  # type: typing.Dict[str, typing.Any]
    exec(textwrap.dedent(s), global_namespace, local_namespace)
    return local_namespace


@contextlib.contextmanager
def ignore_deprecation():
    """Context manager to ignore deprecation warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        yield


ABT_SKIP_MESSAGE = "abstract base class"


def abstract_base_test(cls: _TestCaseType) -> _TestCaseType:
    """Decorator to mark a test class as an "abstract" base class.

    This is different from a regular abstract base class because
    we do not limit instantiation of the class. (If we did, it would
    interfere with test discovery). Instead, we prevent the tests from
    being run.

    Subclasses of an abstract base test are run as normal. There is
    no support for the ``@abstractmethod`` decorator so there is no runtime
    check that all such methods are implemented.

    Note that while it is semantically cleaner to modify the test loader
    to exclude abstract base tests, this is more complicated and would
    interfere with third-party test runners. This approach degrades
    gracefully to other tools such as editor-integrated testing.
    """

    # Type-checking fails due to https://github.com/python/mypy/issues/14458
    # @functools.wraps(cls)
    class AbstractBaseWrapper(cls):  # type: ignore
        @classmethod
        def setUpClass(cls):
            if cls is AbstractBaseWrapper:
                raise unittest.SkipTest(ABT_SKIP_MESSAGE)
            super(AbstractBaseWrapper, cls).setUpClass()

    return AbstractBaseWrapper  # type: ignore
