# -*- coding: utf-8 -*-

"""Stub typing declarations for the native PTY object."""

# Standard library imports
from typing import Optional

# Local imports
from .enums import MouseMode, AgentConfig

__version__: str

class WinptyError(Exception): ...

class PTY:
    def __init__(
        self,
        cols: int,
        rows: int,
        backend: Optional[int] = None,
        mouse_mode: int = MouseMode.WINPTY_MOUSE_MODE_NONE,
        timeout: int = 30000,
        agent_config: int = AgentConfig.WINPTY_FLAG_COLOR_ESCAPES,
    ): ...
    def spawn(
        self,
        appname: str,
        cmdline: Optional[str] = None,
        cwd: Optional[str] = None,
        env: Optional[str] = None,
    ) -> bool: ...
    def set_size(self, cols: int, rows: int): ...
    def read(self, blocking: bool = False) -> str: ...
    def write(self, to_write: str) -> int: ...
    def isalive(self) -> bool: ...
    def get_exitstatus(self) -> Optional[int]: ...
    def iseof(self) -> bool: ...
    def cancel_io(self) -> bool: ...
    @property
    def pid(self) -> Optional[int]: ...
    @property
    def fd(self) -> Optional[int]: ...
