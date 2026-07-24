from __future__ import annotations

__lazy_modules__ = {
    "plumbum.commands.base",
    "plumbum.commands.modifiers",
    "plumbum.commands.processes",
}

from plumbum.commands.base import (
    ERROUT,
    BaseCommand,
    ConcreteCommand,
    shquote,
    shquote_list,
)
from plumbum.commands.modifiers import (
    BG,
    FG,
    NOHUP,
    RETCODE,
    TEE,
    TF,
    ExecutionModifier,
    Future,
)
from plumbum.commands.processes import (
    CommandNotFound,
    ProcessExecutionError,
    ProcessLineTimedOut,
    ProcessTimedOut,
    run_proc,
)

__all__ = (
    "BG",
    "ERROUT",
    "FG",
    "NOHUP",
    "RETCODE",
    "TEE",
    "TF",
    "BaseCommand",
    "CommandNotFound",
    "ConcreteCommand",
    "ExecutionModifier",
    "Future",
    "ProcessExecutionError",
    "ProcessLineTimedOut",
    "ProcessTimedOut",
    "run_proc",
    "shquote",
    "shquote_list",
)


def __dir__() -> list[str]:
    return list(__all__)
