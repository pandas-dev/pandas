import sys
from signal import _HANDLER, _SIGNUM

# technically the implementations will always be around, but since they always
# throw an exception on windows, due to the missing SIGCHLD, we might as well
# pretent they don't exist, but what is different, is that the parameters are
# named even pre 3.10, so we don't just import the symbol from stdlib signal
if sys.platform != "win32":
    def getsignal(signalnum: _SIGNUM) -> _HANDLER: ...
    def signal(signalnum: _SIGNUM, handler: _HANDLER) -> _HANDLER: ...
    def set_wakeup_fd(fd: int, /, *, warn_on_full_buffer: bool = True) -> int: ...

    __all__ = ["signal", "getsignal", "set_wakeup_fd"]
