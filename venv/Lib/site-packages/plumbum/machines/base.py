from __future__ import annotations

__lazy_modules__ = {"plumbum.commands", "plumbum.commands.processes"}

import abc
import os
import typing

from plumbum.commands.processes import (
    CommandNotFound,
    ProcessExecutionError,
    ProcessTimedOut,
)

if typing.TYPE_CHECKING:
    import subprocess
    from collections.abc import Container, Sequence

    from plumbum.commands.base import BaseCommand, ConcreteCommand


StrOrBytesPath = typing.Union[str, bytes, os.PathLike[str], os.PathLike[bytes]]


class PopenAddons:
    """This adds a verify to popen objects to that the correct command is attributed when
    an error is thrown."""

    _proc: subprocess.Popen[bytes]
    custom_encoding: str | None
    returncode: int | None

    def verify(
        self,
        retcode: int | Container[int] | None,
        timeout: float | None,
        stdout: str | bytes,
        stderr: str | bytes,
    ) -> None:
        """This verifies that the correct command is attributed."""
        if getattr(self, "_timed_out", False):
            raise ProcessTimedOut(
                f"Process did not terminate within {timeout} seconds",
                getattr(self, "argv", None),
            )

        if retcode is not None:
            # TODO: this will break if argv is not set, as it doesn't handle None
            if hasattr(retcode, "__contains__"):
                if self.returncode not in retcode:
                    raise ProcessExecutionError(
                        getattr(self, "argv", None),  # type: ignore[arg-type]
                        self.returncode,
                        stdout,
                        stderr,
                    )
            elif self.returncode != retcode:
                raise ProcessExecutionError(
                    getattr(self, "argv", None),  # type: ignore[arg-type]
                    self.returncode,
                    stdout,
                    stderr,
                )


AnyStr = typing.TypeVar("AnyStr", str, bytes)


class PopenWithAddons(typing.Protocol[AnyStr]):
    def verify(
        self,
        retcode: int | Container[int] | None,
        timeout: float | None,
        stdout: AnyStr,
        stderr: AnyStr,
    ) -> None: ...

    custom_encoding: str | None

    args: StrOrBytesPath | Sequence[StrOrBytesPath]
    stdout: typing.IO[AnyStr] | None
    stderr: typing.IO[AnyStr] | None
    stdin: typing.IO[AnyStr] | None
    returncode: int | None
    pid: int | None

    def poll(self) -> int | None: ...
    def communicate(self) -> tuple[AnyStr, AnyStr]: ...
    def send_signal(self, sig: int) -> None: ...
    def terminate(self) -> None: ...
    def wait(self) -> int: ...
    def kill(self) -> None: ...


class MachineCmd:
    __slots__ = ("_machine",)

    def __init__(self, machine: BaseMachine) -> None:
        self._machine = machine

    def __getattr__(self, name: str) -> ConcreteCommand:
        try:
            return self._machine[name]
        except CommandNotFound:
            raise AttributeError(name) from None


class BaseMachine(metaclass=abc.ABCMeta):
    """This is a base class for other machines. It contains common code to
    all machines in Plumbum."""

    __slots__ = ("custom_encoding",)

    custom_encoding: str

    @abc.abstractmethod
    def __getitem__(self, cmd: str) -> ConcreteCommand:
        pass

    def get(self, cmd: str, *othercommands: str) -> ConcreteCommand:
        """This works a little like the ``.get`` method with dict's, only
        it supports an unlimited number of arguments, since later arguments
        are tried as commands and could also fail. It
        will try to call the first command, and if that is not found,
        it will call the next, etc. Will raise if no file named for the
        executable if a path is given, unlike ``[]`` access.

        Usage::

            best_zip = local.get('pigz','gzip')
        """
        try:
            command = self[cmd]
        except CommandNotFound:
            if othercommands:
                return self.get(*othercommands)
            raise

        if not command.executable.exists():
            if othercommands:
                return self.get(*othercommands)
            raise CommandNotFound(cmd, command.executable)

        return command

    def __contains__(self, cmd: str) -> bool:
        """Tests for the existence of the command, e.g., ``"ls" in plumbum.local``.
        ``cmd`` can be anything acceptable by ``__getitem__``.
        """
        try:
            self[cmd]
        except CommandNotFound:
            return False
        return True

    @property
    def encoding(self) -> str:
        "This is a wrapper for custom_encoding"
        return self.custom_encoding

    @encoding.setter
    def encoding(self, value: str) -> None:
        self.custom_encoding = value

    def daemonic_popen(
        self,
        command: BaseCommand,
        cwd: str = "/",
        stdout: str | None = None,
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        raise NotImplementedError("This is not implemented on this machine!")

    Cmd = MachineCmd

    @property
    def cmd(self) -> MachineCmd:
        return self.Cmd(self)

    @abc.abstractmethod
    def clear_program_cache(self) -> None:
        """
        Clear the program cache, which is populated via ``machine.which(progname)`` calls.

        This cache speeds up the lookup of a program in the machines PATH, and is particularly
        effective for RemoteMachines.
        """


__all__ = [
    "BaseMachine",
    "MachineCmd",
    "PopenAddons",
    "PopenWithAddons",
    "StrOrBytesPath",
]


def __dir__() -> list[str]:
    return list(__all__)
