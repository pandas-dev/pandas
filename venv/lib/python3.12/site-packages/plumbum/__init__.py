r"""
Plumbum Shell Combinators
-------------------------
A wrist-handy library for writing shell-like scripts in Python, that can serve
as a ``Popen`` replacement, and much more::

    >>> from plumbum.cmd import ls, grep, wc, cat
    >>> ls()
    'build.py\ndist\ndocs\nLICENSE\nplumbum\nREADME.rst\nsetup.py\ntests\ntodo.txt\n'
    >>> chain = ls["-a"] | grep["-v", "py"] | wc["-l"]
    >>> print(chain)
    /bin/ls -a | /bin/grep -v py | /usr/bin/wc -l
    >>> chain()
    '12\n'
    >>> ((ls["-a"] | grep["-v", "py"]) > "/tmp/foo.txt")()
    ''
    >>> ((cat < "/tmp/foo.txt") | wc["-l"])()
    '12\n'
    >>> from plumbum import local, FG, BG
    >>> with local.cwd("/tmp"):
    ...     (ls | wc["-l"]) & FG
    ...
    13              # printed directly to the interpreter's stdout
    >>> (ls | wc["-l"]) & BG
    <Future ['/usr/bin/wc', '-l'] (running)>
    >>> f = _
    >>> f.stdout    # will wait for the process to terminate
    '9\n'

Plumbum includes local/remote path abstraction, working directory and environment
manipulation, process execution, remote process execution over SSH, tunneling,
SCP-based upload/download, and a {arg|opt}parse replacement for the easy creation
of command-line interface (CLI) programs.

See https://plumbum.readthedocs.io for full details
"""

from __future__ import annotations

# Avoids a circular import error later
import plumbum.path  # noqa: F401
from plumbum.commands import (
    BG,
    ERROUT,
    FG,
    NOHUP,
    RETCODE,
    TEE,
    TF,
    CommandNotFound,
    ProcessExecutionError,
    ProcessLineTimedOut,
    ProcessTimedOut,
)
from plumbum.machines import BaseRemoteMachine, PuttyMachine, SshMachine, local
from plumbum.path import LocalPath, Path, RemotePath
from plumbum.version import version

__author__ = "Tomer Filiba (tomerfiliba@gmail.com)"
__version__ = version

__all__ = (
    "BG",
    "ERROUT",
    "FG",
    "NOHUP",
    "RETCODE",
    "TEE",
    "TF",
    "CommandNotFound",
    "ProcessExecutionError",
    "ProcessLineTimedOut",
    "ProcessTimedOut",
    "BaseRemoteMachine",
    "PuttyMachine",
    "SshMachine",
    "local",
    "LocalPath",
    "Path",
    "RemotePath",
    "__author__",
    "__version__",
    "cmd",
)

from . import cmd


def __dir__():
    "Support nice tab completion"
    return __all__
