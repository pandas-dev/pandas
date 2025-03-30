# -*- coding: utf-8 -*-

"""Stub typing declarations for the native PTY object."""

# Standard library imports
from typing import Optional

# Local imports
from .enums import Backend, Encoding, MouseMode, AgentConfig


class PTY:
    def __init__(self, cols: int, rows: int,
                 backend: Optional[int] = None,
                 encoding: Optional[str] = Encoding.UTF8,
                 mouse_mode: int = MouseMode.WINPTY_MOUSE_MODE_NONE,
                 timeout: int = 30000,
                 agent_config: int = AgentConfig.WINPTY_FLAG_COLOR_ESCAPES):
        ...

    def spawn(self,
              appname: bytes,
              cmdline: Optional[bytes] = None,
              cwd: Optional[bytes] = None,
              env: Optional[bytes] = None) -> bool:
        ...

    def set_size(self, cols: int, rows: int): ...

    def read(self,
             length: Optional[int] = 1000,
             blocking: bool = False) -> bytes:
        ...

    def read_stderr(self,
             length: Optional[int] = 1000,
             blocking: bool = False) -> bytes:
        ...

    def write(self, to_write: bytes) -> int: ...

    def isalive(self) -> bool: ...

    def get_exitstatus(self) -> Optional[int]: ...

    def iseof(self) -> bool: ...

    @property
    def pid(self) -> Optional[int]: ...

    @property
    def fd(self) -> Optional[int]: ...
