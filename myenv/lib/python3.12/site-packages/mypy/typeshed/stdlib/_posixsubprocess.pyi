import sys
from _typeshed import StrOrBytesPath
from collections.abc import Callable, Sequence
from typing import SupportsIndex

if sys.platform != "win32":
    def cloexec_pipe() -> tuple[int, int]: ...
    def fork_exec(
        __args: Sequence[StrOrBytesPath] | None,
        __executable_list: Sequence[bytes],
        __close_fds: bool,
        __pass_fds: tuple[int, ...],
        __cwd: str,
        __env: Sequence[bytes] | None,
        __p2cread: int,
        __p2cwrite: int,
        __c2pread: int,
        __c2pwrite: int,
        __errread: int,
        __errwrite: int,
        __errpipe_read: int,
        __errpipe_write: int,
        __restore_signals: int,
        __call_setsid: int,
        __pgid_to_set: int,
        __gid: SupportsIndex | None,
        __extra_groups: list[int] | None,
        __uid: SupportsIndex | None,
        __child_umask: int,
        __preexec_fn: Callable[[], None],
        __allow_vfork: bool,
    ) -> int: ...
