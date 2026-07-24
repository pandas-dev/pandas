from socket import *

from gevent._hub_primitives import (
    wait_on_watcher,
    wait_read as wait_read,
    wait_readwrite as wait_readwrite,
    wait_write as wait_write,
)
from gevent._types import _Watcher

# This matches the stdlib socket module almost exactly, but contains a couple of extensions
# as a result we just pretend we import everything from socket, which is not entirely correct
# but it gets us most of the way there without having to write a really long list of imports
# with the same platform and version checks, just so we can properly distinguish this module's
# socket class from the native socket class (which could cause issues anyways, since functions
# that accept a socket should still accept the gevent implementation...)
# we can put in the work and do it properly once we have a use-case for it.
# the majority of the gevent implementation can be found in _socket3 and _socketcommon
# which also just imports a lot of symbols from the stdlib socket/_socket module

wait = wait_on_watcher

def cancel_wait(watcher: _Watcher, error: type[BaseException] | BaseException) -> None: ...
