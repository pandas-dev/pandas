from __future__ import annotations

from typing import Optional

# Try to import PTY classes - they may not be available on Windows or without pty feature
try:
    from rattler.rattler import PyPtySession

    _PTY_AVAILABLE = True
except (ImportError, AttributeError):
    _PTY_AVAILABLE = False


class PtySession:
    """
    A pseudoterminal (PTY) session for interactive shell use.

    This is the higher-level PTY API built on top of PtyProcess.
    It provides convenient methods for interactive shell sessions with:
    - Easy command sending via send_line()
    - Interactive mode with wait_until pattern matching
    - Automatic buffering and flushing

    Use this for interactive shell sessions where you want to send commands
    and optionally hand over control to the user.
    """

    _inner: PyPtySession

    def __init__(self, command: list[str]) -> None:
        """
        Create a new PTY session with the given command.

        The PTY session is created with echo enabled by default, which is
        appropriate for interactive shell use.

        Arguments:
            command: A list of strings representing the command and its arguments.
                     The first element is the executable, subsequent elements are arguments.

        Raises:
            RuntimeError: If the PTY session could not be created.

        Examples
        --------
        ```python
        >>> from rattler.pty import PtySession
        >>> # Start an interactive bash session
        >>> session = PtySession(["bash"])
        >>> session.send_line("export MY_VAR=hello")
        >>> # Hand over control to the user (blocks until shell exits)
        >>> exit_code = session.interact()
        >>>
        ```

        If you want to send commands without using interact() (which blocks),
        use exit() to clean up when done:

        ```python
        >>> from rattler.pty import PtySession
        >>> session = PtySession(["bash"])
        >>> session.send_line("echo hello")
        >>> status = session.exit()  # Clean up the session
        >>>
        ```
        """
        if not _PTY_AVAILABLE:
            raise ImportError("PTY functionality is not available on this platform")
        self._inner = PyPtySession(command)

    @classmethod
    def _from_py_pty_session(cls, py_pty_session: PyPtySession) -> PtySession:
        """Construct from FFI PyPtySession object."""
        session = cls.__new__(cls)
        session._inner = py_pty_session
        return session

    def send_line(self, line: str) -> int:
        """
        Send a string followed by a newline to the PTY session.

        This is like typing a command and pressing Enter. The command is flushed
        immediately, so the shell will receive it right away.

        Arguments:
            line: The string to send (newline will be added automatically).

        Returns:
            The number of bytes written.

        Raises:
            RuntimeError: If the write operation fails.

        Examples
        --------
        ```python
        >>> from rattler.pty import PtySession
        >>> session = PtySession(["bash"])
        >>> session.send_line("export MY_VAR=hello")
        22
        >>> session.send_line("echo $MY_VAR")
        13
        >>> session.exit()  # Clean up when done
        'Signaled(SIGTERM)'
        >>>
        ```
        """
        return self._inner.send_line(line)

    def send(self, data: str) -> int:
        """
        Send a string to the PTY session without adding a newline.

        Note: You'll need to call flush() to ensure the data is sent.

        Arguments:
            data: The string to send.

        Returns:
            The number of bytes written.

        Raises:
            RuntimeError: If the write operation fails.

        Examples
        --------
        ```python
        >>> from rattler.pty import PtySession
        >>> session = PtySession(["bash"])
        >>> session.send("echo")
        4
        >>> session.send(" hello\\n")
        7
        >>> session.flush()
        >>> session.exit()  # Clean up when done
        'Signaled(SIGTERM)'
        >>>
        ```
        """
        return self._inner.send(data)

    def flush(self) -> None:
        """
        Flush any pending output to the PTY.

        This is automatically called by send_line(), but can be called manually
        if you use send().

        Raises:
            RuntimeError: If the flush operation fails.

        Examples
        --------
        ```python
        >>> from rattler.pty import PtySession
        >>> session = PtySession(["bash"])
        >>> session.send("echo hello\\n")
        11
        >>> session.flush()  # Make sure the command is sent
        >>> session.exit()  # Clean up when done
        'Signaled(SIGTERM)'
        >>>
        ```
        """
        self._inner.flush()

    def exit(self) -> str:
        """
        Exit the process gracefully by sending SIGTERM.

        This method blocks until the process has exited. Useful for cleaning up
        PTY sessions when you're done sending commands but don't want to use interact().

        Returns:
            A string describing the exit status.

        Raises:
            RuntimeError: If the process could not be terminated.

        Examples
        --------
        ```python
        >>> from rattler.pty import PtySession
        >>> session = PtySession(["bash"])
        >>> session.send_line("echo hello")
        11
        >>> status = session.exit()
        >>> print(status)
        Signaled(SIGTERM)
        >>>
        ```
        """
        return self._inner.exit()

    def interact(self, wait_until: Optional[str] = None) -> Optional[int]:
        """
        Start an interactive session, optionally waiting for a pattern first.

        This method:
        1. Sets the terminal to raw mode
        2. If wait_until is provided, buffers output until that pattern appears
        3. Then forwards all I/O between the user's terminal and the PTY
        4. Returns when the shell process exits

        Arguments:
            wait_until: Optional pattern to wait for before showing output.
                        Useful for waiting for shell initialization to complete.
                        If not found within 1 second, a warning is shown and
                        interaction begins anyway.

        Returns:
            The exit code of the shell process, or None if terminated by signal.

        Raises:
            RuntimeError: If interaction fails.

        Examples
        --------
        ```python
        >>> session = PtySession(["bash"])
        >>> # Send some setup commands
        >>> session.send_line("export MY_VAR=hello")
        >>> session.send_line("echo 'READY'")
        >>> # Wait for "READY" before handing control to user
        >>> exit_code = session.interact(wait_until="READY")
        >>> print(f"Session exited with code: {exit_code}")
        Session exited with code: 0
        >>>
        ```

        ```python
        >>> session = PtySession(["bash"])
        >>> # Interact immediately without waiting
        >>> exit_code = session.interact()
        >>>
        ```
        """
        return self._inner.interact(wait_until)

    def __repr__(self) -> str:
        return "PtySession()"
