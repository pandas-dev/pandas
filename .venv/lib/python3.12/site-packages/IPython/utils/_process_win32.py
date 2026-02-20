"""Windows-specific implementation of process utilities.

This file is only meant to be imported by process.py, not by end-users.
"""

import ctypes
import os
import subprocess
import sys
import time
from ctypes import POINTER, c_int
from ctypes.wintypes import HLOCAL, LPCWSTR
from subprocess import STDOUT
from threading import Thread
from types import TracebackType
from typing import List, Optional

from . import py3compat
from ._process_common import arg_split as py_arg_split

from ._process_common import process_handler, read_no_interrupt
from .encoding import DEFAULT_ENCODING


class AvoidUNCPath:
    """A context manager to protect command execution from UNC paths.

    In the Win32 API, commands can't be invoked with the cwd being a UNC path.
    This context manager temporarily changes directory to the 'C:' drive on
    entering, and restores the original working directory on exit.

    The context manager returns the starting working directory *if* it made a
    change and None otherwise, so that users can apply the necessary adjustment
    to their system calls in the event of a change.

    Examples
    --------
    ::
        cmd = 'dir'
        with AvoidUNCPath() as path:
            if path is not None:
                cmd = '"pushd %s &&"%s' % (path, cmd)
            os.system(cmd)
    """

    def __enter__(self) -> Optional[str]:
        self.path = os.getcwd()
        self.is_unc_path = self.path.startswith(r"\\")
        if self.is_unc_path:
            # change to c drive (as cmd.exe cannot handle UNC addresses)
            os.chdir("C:")
            return self.path
        else:
            # We return None to signal that there was no change in the working
            # directory
            return None

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: TracebackType,
    ) -> None:
        if self.is_unc_path:
            os.chdir(self.path)


def _system_body(p: subprocess.Popen) -> int:
    """Callback for _system."""
    enc = DEFAULT_ENCODING

    # Dec 2024: in both of these functions, I'm not sure why we .splitlines()
    # the bytes and then decode each line individually instead of just decoding
    # the whole thing at once.
    def stdout_read() -> None:
        try:
            assert p.stdout is not None
            for byte_line in read_no_interrupt(p.stdout).splitlines():
                line = byte_line.decode(enc, "replace")
                print(line, file=sys.stdout)
        except Exception as e:
            print(f"Error reading stdout: {e}", file=sys.stderr)

    def stderr_read() -> None:
        try:
            assert p.stderr is not None
            for byte_line in read_no_interrupt(p.stderr).splitlines():
                line = byte_line.decode(enc, "replace")
                print(line, file=sys.stderr)
        except Exception as e:
            print(f"Error reading stderr: {e}", file=sys.stderr)

    stdout_thread = Thread(target=stdout_read)
    stderr_thread = Thread(target=stderr_read)

    stdout_thread.start()
    stderr_thread.start()

    # Wait to finish for returncode. Unfortunately, Python has a bug where
    # wait() isn't interruptible (https://bugs.python.org/issue28168) so poll in
    # a loop instead of just doing `return p.wait()`
    while True:
        result = p.poll()
        if result is None:
            time.sleep(0.01)
        else:
            break

    # Join the threads to ensure they complete before returning
    stdout_thread.join()
    stderr_thread.join()

    return result


def system(cmd: str) -> Optional[int]:
    """Win32 version of os.system() that works with network shares.

    Note that this implementation returns None, as meant for use in IPython.

    Parameters
    ----------
    cmd : str or list
        A command to be executed in the system shell.

    Returns
    -------
    int : child process' exit code.
    """
    # The controller provides interactivity with both
    # stdin and stdout
    # import _process_win32_controller
    # _process_win32_controller.system(cmd)

    with AvoidUNCPath() as path:
        if path is not None:
            cmd = '"pushd %s &&"%s' % (path, cmd)
        res = process_handler(cmd, _system_body)
        assert isinstance(res, int | type(None))
        return res
    return None


def getoutput(cmd: str) -> str:
    """Return standard output of executing cmd in a shell.

    Accepts the same arguments as os.system().

    Parameters
    ----------
    cmd : str or list
        A command to be executed in the system shell.

    Returns
    -------
    stdout : str
    """

    with AvoidUNCPath() as path:
        if path is not None:
            cmd = '"pushd %s &&"%s' % (path, cmd)
        out = process_handler(cmd, lambda p: p.communicate()[0], STDOUT)

    if out is None:
        out = b""
    return py3compat.decode(out)


try:
    windll = ctypes.windll  # type: ignore [attr-defined]
    CommandLineToArgvW = windll.shell32.CommandLineToArgvW
    CommandLineToArgvW.arg_types = [LPCWSTR, POINTER(c_int)]
    CommandLineToArgvW.restype = POINTER(LPCWSTR)
    LocalFree = windll.kernel32.LocalFree
    LocalFree.res_type = HLOCAL
    LocalFree.arg_types = [HLOCAL]

    def arg_split(
        commandline: str, posix: bool = False, strict: bool = True
    ) -> List[str]:
        """Split a command line's arguments in a shell-like manner.

        This is a special version for windows that use a ctypes call to CommandLineToArgvW
        to do the argv splitting. The posix parameter is ignored.

        If strict=False, process_common.arg_split(...strict=False) is used instead.
        """
        # CommandLineToArgvW returns path to executable if called with empty string.
        if commandline.strip() == "":
            return []
        if not strict:
            # not really a cl-arg, fallback on _process_common
            return py_arg_split(commandline, posix=posix, strict=strict)
        argvn = c_int()
        result_pointer = CommandLineToArgvW(commandline.lstrip(), ctypes.byref(argvn))
        try:
            result_array_type = LPCWSTR * argvn.value
            result = [
                arg
                for arg in result_array_type.from_address(
                    ctypes.addressof(result_pointer.contents)
                )
                if arg is not None
            ]
        finally:
            # for side effects
            _ = LocalFree(result_pointer)
        return result
except AttributeError:
    arg_split = py_arg_split


def check_pid(pid: int) -> bool:
    # OpenProcess returns 0 if no such process (of ours) exists
    # positive int otherwise
    return bool(windll.kernel32.OpenProcess(1, 0, pid))
