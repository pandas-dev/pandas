"""Asyncio support for Plumbum commands.

This module provides async versions of Plumbum commands that can be used with
Python's asyncio framework. Commands can be awaited directly or used with
async context managers.

For async machines (AsyncLocalMachine, AsyncSshMachine), see plumbum.machines.local
and plumbum.machines.ssh_machine respectively

Design Philosophy
-----------------
This implementation uses **delegation over inheritance** to wrap existing sync
commands rather than reimplementing their functionality. This approach:

- Maximizes code reuse (~100 lines of logic delegated to sync commands)
- Ensures consistency between sync and async APIs
- Reduces maintenance burden (changes to sync code automatically apply)
- Enables automatic support for features like `with_env()` and `with_cwd()`

Why Delegation Instead of Inheritance?
---------------------------------------
Sync and async methods are fundamentally incompatible in Python:

- Sync methods return values directly: `def run() -> tuple[int, str, str]`
- Async methods return coroutines: `async def run() -> AsyncResult`

You cannot override a sync method with an async one in the same class hierarchy.
Therefore, we use separate classes that wrap and delegate to sync commands,
reusing all their formulation, binding, and pipeline logic.

Example Usage
-------------
::

    from plumbum import async_local

    async def main():
        # Simple command execution
        result = await async_local["ls"]("-la")
        print(result)

        # With explicit run method
        ls = async_local["ls"]
        result = await ls.run(["-la"])
        print(result.stdout)

        # Pipeline support
        result = await (async_local["ls"] | async_local["grep"]["py"])()
        print(result)

.. versionadded:: 2.0
"""

from __future__ import annotations

__lazy_modules__ = {
    "asyncio",
    "contextlib",
    "plumbum.commands.processes",
    "plumbum.machines.local",
    "typing_extensions",
}

import asyncio
import codecs
import contextlib
import os
import sys
from typing import TYPE_CHECKING, Any, TextIO

from plumbum.commands.base import BoundCommand, BoundEnvCommand, Pipeline
from plumbum.commands.processes import ProcessExecutionError, ProcessTimedOut
from plumbum.machines.local import local

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

if TYPE_CHECKING:
    from collections.abc import Container, Coroutine, Sequence

    from plumbum.commands.base import BaseCommand


class AsyncResult:
    """Result of an async command execution.

    Attributes:
        returncode: The exit code of the process
        stdout: Standard output as a string
        stderr: Standard error as a string
    """

    __slots__ = ("returncode", "stderr", "stdout")

    def __init__(self, returncode: int, stdout: str, stderr: str):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr

    def __str__(self) -> str:
        return self.stdout

    def __repr__(self) -> str:
        return f"AsyncResult(returncode={self.returncode}, stdout={self.stdout!r}, stderr={self.stderr!r})"


class AsyncPipelineProcess:
    """Proxy around the downstream process of an async pipeline.

    Output attributes (``stdout``, ``stderr``, ``pid``, ...) are delegated to the
    downstream :class:`asyncio.subprocess.Process`, while ``stdin`` is taken from
    the *upstream* process so writes feed the head of the pipeline (mirroring
    :meth:`plumbum.commands.base.Pipeline.popen`, which sets
    ``dstproc.stdin = srcproc.stdin``). :meth:`wait` and :meth:`communicate` reap
    both stages and report a combined return code -- the downstream code, or the
    upstream code if the downstream stage succeeded. :meth:`kill`,
    :meth:`terminate` and :meth:`send_signal` are propagated to every stage
    (recursing through a nested upstream pipeline) rather than only the last one.

    A separate proxy is needed (rather than patching ``wait`` on the downstream
    process) because :attr:`asyncio.subprocess.Process.returncode` is a read-only
    property and cannot be reassigned to the combined value.

    Instances are returned by :meth:`AsyncCommandMixin.popen` for pipelines; you
    don't normally construct them directly.

    .. versionadded:: 1.11
    """

    def __init__(
        self,
        dstproc: asyncio.subprocess.Process | AsyncPipelineProcess,
        srcproc: asyncio.subprocess.Process | AsyncPipelineProcess,
    ) -> None:
        # ``dstproc`` is always a real Process (a pipeline's last stage is a
        # single command); ``srcproc`` may be another proxy for nested pipelines.
        self._dstproc = dstproc
        self.srcproc = srcproc
        self._returncode: int | None = None

    def __getattr__(self, name: str) -> Any:
        # Only called when `name` isn't a real attribute on the proxy, so this
        # delegates the rest of the Process interface to the downstream process.
        return getattr(self._dstproc, name)

    @property
    def stdin(self) -> asyncio.StreamWriter | None:
        # The downstream stage reads from the OS pipe, so its own stdin is None;
        # the pipeline's stdin is the head (upstream) process's stdin.
        return self.srcproc.stdin

    @property
    def returncode(self) -> int | None:
        if self._returncode is not None:
            return self._returncode
        return self._dstproc.returncode

    def _combine(self, rc_dst: int | None, rc_src: int | None) -> int:
        self._returncode = (rc_dst or rc_src) or 0
        return self._returncode

    def _signal_both(self, method: str) -> None:
        # Signal both ends; for a nested upstream pipeline ``srcproc`` is itself
        # an ``AsyncPipelineProcess``, so this recurses through every stage.
        # Each stage may already have exited, so ignore "no such process".
        for proc in (self._dstproc, self.srcproc):
            with contextlib.suppress(ProcessLookupError):
                getattr(proc, method)()

    def kill(self) -> None:
        self._signal_both("kill")

    def terminate(self) -> None:
        self._signal_both("terminate")

    def send_signal(self, signal: int) -> None:
        for proc in (self._dstproc, self.srcproc):
            with contextlib.suppress(ProcessLookupError):
                proc.send_signal(signal)

    async def wait(self) -> int:
        rc_dst = await self._dstproc.wait()
        rc_src = await self.srcproc.wait()
        return self._combine(rc_dst, rc_src)

    async def _reap(self) -> tuple[bytes | None, bytes | None]:
        # Reap the whole pipeline, draining every stage via communicate() so a
        # full stderr pipe can't deadlock the wait. The last stage's stdout is
        # its captured output; an upstream stage's stdout is an OS pipe (no
        # StreamReader), so communicate() there just drains its stderr. Does not
        # feed stdin -- the head's stdin is fed by communicate() below.
        src = self.srcproc
        reap_src = (
            src._reap() if isinstance(src, AsyncPipelineProcess) else src.communicate()
        )
        (stdout, stderr), _ = await asyncio.gather(
            self._dstproc.communicate(), reap_src
        )
        self._combine(self._dstproc.returncode, src.returncode)
        return stdout, stderr

    async def communicate(
        self, input: bytes | None = None
    ) -> tuple[bytes | None, bytes | None]:
        # Feed the pipeline's stdin (the head stage) concurrently with reaping
        # every stage -- draining each one's stderr, and the last stage's
        # stdout/stderr, so a full pipe on any stage can't deadlock the wait.
        async def feed() -> None:
            writer = self.stdin
            if writer is None:
                return
            # As in asyncio's own Process.communicate(), a stage that exits
            # early closes the pipe, so ignore the resulting write errors.
            with contextlib.suppress(BrokenPipeError, ConnectionResetError):
                try:
                    if input:
                        writer.write(input)
                        await writer.drain()
                finally:
                    writer.close()
                    await writer.wait_closed()

        _, (stdout, stderr) = await asyncio.gather(feed(), self._reap())
        return stdout, stderr


class AsyncCommandMixin:
    """Mixin that adds async execution capabilities to BaseCommand.

    This mixin wraps a sync BaseCommand and provides async execution methods
    while reusing all the existing formulation, binding, and pipeline logic.

    The delegation pattern allows us to:
    - Reuse BaseCommand.formulate() for command-to-argv conversion
    - Reuse BaseCommand.__getitem__() for argument binding
    - Reuse BaseCommand.__or__() for pipeline creation
    - Reuse BoundEnvCommand for environment and cwd handling
    - Maintain consistency with the sync API
    """

    __slots__ = ("_base_cmd",)

    _base_cmd: BaseCommand

    def __init__(self, base_cmd: BaseCommand):
        """Initialize with a sync BaseCommand to wrap.

        Args:
            base_cmd: The sync command to wrap and delegate to
        """
        self._base_cmd = base_cmd

    @property
    def _concrete(self) -> BaseCommand:
        """The innermost command, unwrapping any ``Bound``/``BoundEnv`` layers.

        Binding arguments or an environment wraps the base command, so the
        concrete command (the one carrying ``executable``/``remote``) lives
        underneath. This is used by the subclass ``executable``/``remote``
        properties so they keep working after ``[...]``/``with_env``/``with_cwd``.
        """
        cmd = self._base_cmd
        while isinstance(cmd, (BoundCommand, BoundEnvCommand)):
            cmd = cmd.cmd
        return cmd

    def __getitem__(self, args: Any) -> Self:
        """Bind arguments using the base command's logic.

        This delegates to the sync command's __getitem__ method, which handles
        all the argument binding logic, then re-wraps the result in the same
        async type as ``self`` (preserving e.g. ``AsyncLocalCommand``).
        """
        bound = self._base_cmd[args]
        return self.__class__(bound)

    def with_env(self, **env: str) -> Self:
        """Return a new async command with the given environment variables.

        Delegates to the sync command's ``with_env`` (which produces a
        ``BoundEnvCommand``) and re-wraps it in the same async type.
        """
        return self.__class__(self._base_cmd.with_env(**env))

    def with_cwd(self, path: Any) -> Self:
        """Return a new async command with the given working directory.

        Delegates to the sync command's ``with_cwd`` and re-wraps the result.
        """
        return self.__class__(self._base_cmd.with_cwd(path))

    def __call__(self, *args: Any, **kwargs: Any) -> Coroutine[Any, Any, str]:
        """Execute the command asynchronously and return stdout.

        This is a shortcut for run() that returns only stdout, matching the
        behavior of the sync API's __call__ method.
        """

        async def _run() -> str:
            result = await self.run(args, **kwargs)
            return result.stdout

        return _run()

    def __or__(self, other: AsyncCommandMixin) -> Self:
        """Create a pipeline using the base command's logic.

        This delegates to the sync command's __or__ method to create a sync
        Pipeline, then wraps it with the same type as self.
        """
        sync_pipeline = self._base_cmd | other._base_cmd
        return self.__class__(sync_pipeline)

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        """Delegate formulation to the base command.

        This reuses the sync command's formulation logic, which handles:
        - Converting the command to an argv list
        - Proper shell quoting based on nesting level
        - Handling of bound arguments
        - Support for nested commands
        """
        return self._base_cmd.formulate(level, args)

    async def run(
        self,
        args: Sequence[Any] = (),
        retcode: int | Container[int] | None = 0,
        timeout: float | None = None,
        cwd: str | None = None,
        env: dict[str, str] | None = None,
    ) -> AsyncResult:
        """Run the command asynchronously.

        Args:
            args: Additional arguments to pass to the command
            retcode: Expected return code(s). None to disable checking.
            timeout: Maximum time to wait for command completion
            cwd: Working directory for the command
            env: Environment variables for the command

        Returns:
            AsyncResult with returncode, stdout, and stderr

        Raises:
            ProcessExecutionError: If return code doesn't match expected
            asyncio.TimeoutError: If timeout is exceeded
        """
        loop = asyncio.get_running_loop()

        def _run_sync() -> tuple[int, str, str]:
            retcode_val, stdout, stderr = self._base_cmd.run(
                args, retcode=retcode, timeout=timeout, cwd=cwd, env=env
            )
            return retcode_val, stdout, stderr

        try:
            retcode_val, stdout, stderr = await loop.run_in_executor(None, _run_sync)
        except ProcessTimedOut as e:
            raise asyncio.TimeoutError() from e

        return AsyncResult(retcode_val, stdout, stderr)

    async def popen(
        self,
        args: Sequence[Any] = (),
        cwd: str | None = None,
        env: dict[str, str] | None = None,
    ) -> asyncio.subprocess.Process | AsyncPipelineProcess:
        """Create an async subprocess without waiting for it to complete.

        This is useful for long-running processes or when you need to
        interact with stdin/stdout/stderr.

        .. note::
            Streaming is only supported for **local** commands. Remote commands
            raise :class:`NotImplementedError` here -- run them with ``.run()``
            or by calling the command directly (those execute the underlying
            sync command in a thread).

        Args:
            args: Additional arguments to pass to the command
            cwd: Working directory for the command
            env: Environment variables for the command

        Returns:
            An :class:`asyncio.subprocess.Process` for a plain command. For a
            pipeline, an :class:`AsyncPipelineProcess` proxy that exposes the same
            interface (``stdout``/``stderr`` from the last stage, ``stdin`` to
            the first) while reaping every stage on ``wait``/``communicate``.
        """
        return await self._popen(args, cwd=cwd, env=env)

    async def _popen(
        self,
        args: Sequence[Any] = (),
        *,
        cwd: str | None = None,
        env: dict[str, str] | None = None,
        stdin: int | None = asyncio.subprocess.PIPE,
        stdout: int | None = asyncio.subprocess.PIPE,
        stderr: int | None = asyncio.subprocess.PIPE,
    ) -> asyncio.subprocess.Process | AsyncPipelineProcess:
        """Spawn the subprocess, threading stdin/stdout through pipeline stages.

        The public ``popen`` always uses pipes for all three streams, but
        pipelines need to connect one stage's stdout to the next stage's stdin,
        so this internal helper exposes those handles.
        """
        # Streaming via ``asyncio.create_subprocess_exec`` only works for local
        # commands -- a RemoteCommand would be formulated without its ssh wrapper
        # and silently run on *this* machine. Remote commands are still supported
        # through ``run()``/``__call__`` (which delegate to the sync command in a
        # thread); fail loudly here rather than executing the wrong thing.
        if self._base_cmd.machine is not local:
            raise NotImplementedError(
                "Async popen/streaming (popen, AsyncTEE) is only supported for "
                "local commands. Use .run() or call the command directly for "
                "remote commands -- those execute the sync command in a thread."
            )

        # Binding args/env/cwd wraps the command in BoundCommand/BoundEnvCommand
        # (possibly around a Pipeline -- e.g. (a | b)["--flag"]). Unwrap those to
        # reach the concrete command (or Pipeline) and thread the bound
        # args/env/cwd inward, so they take effect for plain commands too.
        base: BaseCommand = self._base_cmd
        cmd_args = list(args)
        cmd_env = env
        cmd_cwd = cwd
        while isinstance(base, (BoundCommand, BoundEnvCommand)):
            if isinstance(base, BoundCommand):
                cmd_args = [*base.args, *cmd_args]
            else:
                cmd_env = {**base.env, **(cmd_env or {})}
                cmd_cwd = base.cwd if cmd_cwd is None else cmd_cwd
            base = base.cmd

        # A pipeline has no single argv; connect the two stages with an OS pipe,
        # mirroring the synchronous Pipeline.popen (extra args go to the head
        # stage, as in that implementation).
        if isinstance(base, Pipeline):
            read_fd, write_fd = os.pipe()
            # The parent closes each fd once the child has inherited it. The
            # nested try ensures both fds are closed even if either spawn raises
            # (write_fd by the inner finally, read_fd by the outer one).
            try:
                try:
                    srcproc = await self.__class__(base.srccmd)._popen(
                        cmd_args,
                        cwd=cmd_cwd,
                        env=cmd_env,
                        stdin=stdin,
                        stdout=write_fd,
                        stderr=stderr,
                    )
                finally:
                    os.close(write_fd)

                try:
                    dstproc = await self.__class__(base.dstcmd)._popen(
                        cwd=cmd_cwd,
                        env=cmd_env,
                        stdin=read_fd,
                        stdout=stdout,
                        stderr=stderr,
                    )
                except BaseException:
                    # The upstream stage is already running; reap it rather than
                    # leaking a child process and its open pipe ends.
                    with contextlib.suppress(ProcessLookupError):
                        srcproc.kill()
                    await srcproc.wait()
                    raise
            finally:
                os.close(read_fd)

            # Waiting on the pipeline reaps the upstream process too, and the
            # return code is the downstream stage's, or the upstream stage's if
            # the downstream one succeeded (a pipefail-like combination, as in
            # the synchronous Pipeline.popen).
            return AsyncPipelineProcess(dstproc, srcproc)

        # Formulate at level 0 (no shell-quoting -- we exec the argv directly),
        # matching the synchronous LocalCommand.popen.
        argv = base.formulate(0, cmd_args)

        full_env = dict(local.env.getdict())
        if cmd_env:
            full_env.update(cmd_env)

        working_dir = cmd_cwd or str(local.cwd)

        return await asyncio.create_subprocess_exec(
            *argv,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            cwd=working_dir,
            env=full_env,
        )


class AsyncCommand(AsyncCommandMixin):
    """Async wrapper for BaseCommand.

    This class wraps any BaseCommand and provides async execution capabilities.
    It reuses all the formulation, binding, and pipeline logic from the base command.

    Example::

        # The sync command is looked up and wrapped
        async_cmd = async_local["ls"]

        # Binding works via delegation to sync command
        bound_cmd = async_cmd["-la"]

        # Execution is async
        result = await bound_cmd.run()
    """

    __slots__ = ()


class AsyncLocalCommand(AsyncCommand):
    """Async version of LocalCommand.

    This class wraps a LocalCommand and provides async execution methods.
    It reuses all the LocalCommand logic for formulation, binding, etc.
    """

    __slots__ = ()

    @property
    def executable(self) -> Any:
        """The path to the executable."""
        # Unwrap any bound layers; the concrete LocalCommand carries `executable`.
        return self._concrete.executable  # type: ignore[attr-defined]


class AsyncRemoteCommand(AsyncCommand):
    """Async wrapper for RemoteCommand.

    This class wraps a RemoteCommand and provides async execution capabilities.
    It reuses all the RemoteCommand logic for formulation, binding, etc.

    Example::

        async with AsyncSshMachine("host") as rem:
            ls = rem["ls"]
            result = await ls("-la")

    .. versionadded:: 2.0
    """

    __slots__ = ()

    @property
    def remote(self) -> Any:
        """The remote machine this command belongs to."""
        # Unwrap any bound layers; the concrete RemoteCommand carries `remote`.
        return self._concrete.remote  # type: ignore[attr-defined]


# ===================================================================================================
# Async execution modifiers
# ===================================================================================================


class AsyncExecutionModifier:
    """Base class for async execution modifiers."""

    __slots__ = ("__weakref__",)

    def __repr__(self) -> str:
        """Automatically creates a representation for given subclass with slots."""
        slots = {}
        for cls in self.__class__.__mro__:
            slots_list = getattr(cls, "__slots__", ())
            if isinstance(slots_list, str):
                slots_list = (slots_list,)
            for prop in slots_list:
                if prop[0] != "_":
                    slots[prop] = getattr(self, prop)
        mystrs = (f"{name} = {value}" for name, value in slots.items())
        mystrs_str = ", ".join(mystrs)
        return f"{self.__class__.__name__}({mystrs_str})"

    @classmethod
    def __call__(cls, *args: Any, **kwargs: Any) -> Self:
        return cls(*args, **kwargs)


class _AsyncTF(AsyncExecutionModifier):
    """Async execution modifier that returns True/False based on return code.

    This is the async equivalent of the sync TF modifier. It runs the command
    and returns True if the exit code matches the expected value, False otherwise.

    Unlike the sync version, there is no FG parameter because async commands
    don't have a concept of foreground/background execution - they're all
    non-blocking by nature.

    Example::

        # Check if a file exists
        exists = await (async_local["test"]["-f", "file.txt"] & AsyncTF)

        # Check for specific exit code
        result = await (async_local["grep"]["pattern", "file.txt"] & AsyncTF(retcode=(0, 1)))

    .. versionadded:: 2.0
    """

    __slots__ = ("retcode", "timeout")

    def __init__(
        self,
        retcode: int | Container[int] = 0,
        timeout: float | None = None,
    ) -> None:
        """Initialize AsyncTF modifier.

        Args:
            retcode: Expected return code(s). Default is 0.
            timeout: Maximum time to wait for command completion.
        """
        self.retcode = retcode
        self.timeout = timeout

    async def __rand__(self, cmd: AsyncCommandMixin) -> bool:
        """Execute command and return True/False based on return code."""
        try:
            await cmd.run(retcode=self.retcode, timeout=self.timeout)
        except ProcessExecutionError:
            return False
        return True


class _AsyncRETCODE(AsyncExecutionModifier):
    """Async execution modifier that returns only the exit code.

    This is the async equivalent of the sync RETCODE modifier. It runs the
    command and returns only the exit code, ignoring stdout/stderr.

    Unlike the sync version, there is no FG parameter because async commands
    don't have a concept of foreground/background execution.

    Example::

        # Get exit code
        code = await (async_local["ls"]["/nonexistent"] & AsyncRETCODE)
        print(f"Exit code: {code}")

    .. versionadded:: 2.0
    """

    __slots__ = ("timeout",)

    def __init__(self, timeout: float | None = None) -> None:
        """Initialize AsyncRETCODE modifier.

        Args:
            timeout: Maximum time to wait for command completion.
        """
        self.timeout = timeout

    async def __rand__(self, cmd: AsyncCommandMixin) -> int:
        """Execute command and return exit code."""
        result = await cmd.run(retcode=None, timeout=self.timeout)
        return result.returncode


class _AsyncTEE(AsyncExecutionModifier):
    """Async execution modifier that displays output in real-time and returns it.

    This is the async equivalent of the sync TEE modifier. It runs the command,
    displays stdout/stderr to the console in real-time, and also returns them.

    Unlike the sync version, buffering is always enabled because async I/O
    handles buffering differently.

    Example::

        # Run command and see output in real-time
        retcode, stdout, stderr = await (async_local["npm"]["install"] & AsyncTEE)

        # With custom expected return code
        retcode, stdout, stderr = await (async_local["grep"]["pattern"] & AsyncTEE(retcode=(0, 1)))

    .. versionadded:: 2.0
    """

    __slots__ = ("retcode", "timeout")

    def __init__(
        self,
        retcode: int | Container[int] = 0,
        timeout: float | None = None,
    ) -> None:
        """Initialize AsyncTEE modifier.

        Args:
            retcode: Expected return code(s). Default is 0.
            timeout: Maximum time to wait for command completion.
        """
        self.retcode = retcode
        self.timeout = timeout

    async def __rand__(self, cmd: AsyncCommandMixin) -> tuple[int, str, str]:
        """Execute command, display output, and return (retcode, stdout, stderr)."""
        encoding = cmd._base_cmd._get_encoding() or local.custom_encoding

        # Reuse ``_popen`` so this works for pipelines too (its own
        # create_subprocess_exec path only handled single commands), and so the
        # local-only guard applies. stdin is inherited from the parent process.
        proc = await cmd._popen((), stdin=None)

        # Collect output while displaying it
        stdout_lines: list[str] = []
        stderr_lines: list[str] = []

        async def read_stream(
            stream: asyncio.StreamReader | None,
            target: TextIO,
            output_list: list[str] | None = None,
        ) -> None:
            """Read a stream in fixed-size chunks, display, and optionally collect.

            Reads by chunk rather than ``readline``: ``StreamReader.readline``
            raises ``LimitOverrunError``/``ValueError`` once a line exceeds the
            stream's buffer limit (64KiB by default). A single incremental
            decoder per stream keeps multibyte characters that straddle a chunk
            boundary intact.
            """
            if stream is None:
                return

            decoder = codecs.getincrementaldecoder(encoding)(errors="ignore")
            while True:
                chunk = await stream.read(4096)
                text = decoder.decode(chunk, final=not chunk)
                if text:
                    if output_list is not None:
                        output_list.append(text)
                    target.write(text)
                    target.flush()
                if not chunk:
                    break

        # ``proc.stdout``/``proc.stderr`` are the pipeline's *final* stage. Each
        # upstream stage of a pipeline has its own stderr pipe that must also be
        # drained -- otherwise a chatty upstream stderr fills its buffer and
        # deadlocks the wait() below (only communicate() drains every stage). We
        # display these too, like a shell pipeline; the returned stderr stays the
        # final stage's, matching ``popen().communicate()``. The loop adds nothing
        # for a plain command (no upstream stages).
        readers = [
            read_stream(proc.stdout, sys.stdout, stdout_lines),
            read_stream(proc.stderr, sys.stderr, stderr_lines),
        ]
        node: Any = proc
        while isinstance(node, AsyncPipelineProcess):
            node = node.srcproc
            readers.append(read_stream(node.stderr, sys.stderr))

        # Read every stage's streams concurrently
        try:
            await asyncio.wait_for(
                asyncio.gather(*readers),
                timeout=self.timeout,
            )
        except BaseException:
            # On timeout or any other failure (e.g. a decode/read error or
            # cancellation), make sure the subprocess is killed and reaped
            # instead of being left running with open pipes.
            with contextlib.suppress(Exception):
                proc.kill()
            with contextlib.suppress(Exception):
                await proc.wait()
            raise

        # Wait for process to complete (reaps every pipeline stage)
        await proc.wait()

        # Combine output
        stdout = "".join(stdout_lines)
        stderr = "".join(stderr_lines)
        retcode = proc.returncode

        # Check return code
        if self.retcode is not None:
            expected_codes: set[int] = (
                {self.retcode} if isinstance(self.retcode, int) else set(self.retcode)  # type: ignore[call-overload]
            )

            if retcode not in expected_codes:
                raise ProcessExecutionError(
                    argv=cmd.formulate(0, ()),
                    retcode=retcode,
                    stdout=stdout,
                    stderr=stderr,
                )

        return (retcode or 0), stdout, stderr


# Singleton instances
AsyncTF = _AsyncTF()
AsyncRETCODE = _AsyncRETCODE()
AsyncTEE = _AsyncTEE()


__all__ = (
    "AsyncCommand",
    "AsyncLocalCommand",
    "AsyncPipelineProcess",
    "AsyncRETCODE",
    "AsyncRemoteCommand",
    "AsyncResult",
    "AsyncTEE",
    "AsyncTF",
)


def __dir__() -> list[str]:
    return list(__all__)
