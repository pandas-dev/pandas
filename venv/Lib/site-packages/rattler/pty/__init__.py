"""
PTY (pseudoterminal) support for interactive shell sessions.

This module provides Python bindings to rattler_pty, allowing you to:
- Create interactive shell sessions in a PTY
- Send commands programmatically
- Hand over control to the user with interact()
- Manage process lifecycle

Note: PTY functionality is only available on Unix platforms (Linux, macOS).
On Windows or when the pty feature is not compiled, importing from this
module will raise an ImportError.
"""

from __future__ import annotations

try:
    from rattler.pty.pty_process import PtyProcess, PtyProcessOptions
    from rattler.pty.pty_session import PtySession

    _PTY_AVAILABLE = True
    __all__ = ["PtyProcess", "PtyProcessOptions", "PtySession"]
except ImportError:
    _PTY_AVAILABLE = False
    __all__ = []
