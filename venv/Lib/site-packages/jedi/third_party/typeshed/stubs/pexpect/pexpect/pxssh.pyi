from _typeshed import FileDescriptorOrPath
from collections.abc import Mapping
from typing import AnyStr, Literal

from .exceptions import ExceptionPexpect
from .pty_spawn import spawn
from .spawnbase import _Logfile

__all__ = ["ExceptionPxssh", "pxssh"]

class ExceptionPxssh(ExceptionPexpect): ...

class pxssh(spawn[AnyStr]):
    name: str
    UNIQUE_PROMPT: str
    PROMPT: str
    PROMPT_SET_SH: str
    PROMPT_SET_CSH: str
    PROMPT_SET_ZSH: str
    SSH_OPTS: str
    force_password: bool
    debug_command_string: bool
    options: dict[str, str]
    def __init__(
        self,
        timeout: float | None = 30,
        maxread: int = 2000,
        searchwindowsize: int | None = None,
        logfile: _Logfile | None = None,
        cwd: FileDescriptorOrPath | None = None,
        env: Mapping[str, str] | None = None,
        ignore_sighup: bool = True,
        echo: bool = True,
        options: dict[str, str] = {},
        encoding: str | None = None,
        codec_errors: str = "strict",
        debug_command_string: bool = False,
        use_poll: bool = False,
    ) -> None: ...
    def levenshtein_distance(self, a, b): ...
    def try_read_prompt(self, timeout_multiplier): ...
    def sync_original_prompt(self, sync_multiplier: float = 1.0): ...
    def login(
        self,
        server,
        username: str | None = None,
        password: str = "",
        terminal_type: str = "ansi",
        original_prompt: str = "[#$]",
        login_timeout: float | None = 10,
        port: int | None = None,
        auto_prompt_reset: bool = True,
        ssh_key: FileDescriptorOrPath | Literal[True] | None = None,
        quiet: bool = True,
        sync_multiplier: int = 1,
        check_local_ip: bool = True,
        password_regex: str = "(?i)(?:password:)|(?:passphrase for key)",
        ssh_tunnels: dict[str, list[str | int]] = {},
        spawn_local_ssh: bool = True,
        sync_original_prompt: bool = True,
        ssh_config: FileDescriptorOrPath | None = None,
        cmd: str = "ssh",
    ): ...
    def logout(self) -> None: ...
    def prompt(self, timeout: float | None = -1): ...
    def set_unique_prompt(self): ...
