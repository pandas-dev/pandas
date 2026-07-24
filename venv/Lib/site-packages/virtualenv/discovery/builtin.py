"""Virtualenv-specific Builtin discovery wrapping py_discovery."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from python_discovery import get_interpreter as _get_interpreter

from .discover import Discover

if TYPE_CHECKING:
    from argparse import ArgumentParser
    from collections.abc import Iterable, Mapping, Sequence

    from python_discovery import PyInfoCache, PythonInfo

    from virtualenv.config.cli.parser import VirtualEnvOptions


def get_interpreter(
    key: str,
    try_first_with: Iterable[str],
    cache: PyInfoCache | None = None,
    env: Mapping[str, str] | None = None,
    app_data: PyInfoCache | None = None,
) -> PythonInfo | None:
    return _get_interpreter(key, try_first_with, cache or app_data, env)


class Builtin(Discover):
    python_spec: Sequence[str]
    app_data: PyInfoCache
    try_first_with: Sequence[str]

    def __init__(self, options: VirtualEnvOptions) -> None:
        super().__init__(options)
        self.python_spec = options.python or [sys.executable]
        if self._env.get("VIRTUALENV_PYTHON"):
            self.python_spec = self.python_spec[1:] + self.python_spec[:1]
        self.app_data = options.app_data
        self.try_first_with = options.try_first_with

    @classmethod
    def add_parser_arguments(cls, parser: ArgumentParser) -> None:
        parser.add_argument(
            "-p",
            "--python",
            dest="python",
            metavar="py",
            type=str,
            action="append",
            default=[],
            help="interpreter based on what to create environment (path/identifier/version-specifier) "
            "- by default use the interpreter where the tool is installed - first found wins. "
            "Version specifiers (e.g., >=3.12, ~=3.11.0, ==3.10) are also supported",
        )
        parser.add_argument(
            "--try-first-with",
            dest="try_first_with",
            metavar="py_exe",
            type=str,
            action="append",
            default=[],
            help="try first these interpreters before starting the discovery",
        )

    def run(self) -> PythonInfo | None:
        for python_spec in self.python_spec:
            if result := get_interpreter(
                python_spec,
                self.try_first_with,
                app_data=self.app_data,
                env=self._env,
            ):
                return result
        return None

    def __repr__(self) -> str:
        spec = self.python_spec[0] if len(self.python_spec) == 1 else self.python_spec
        return f"{self.__class__.__name__} discover of python_spec={spec!r}"


__all__ = [
    "Builtin",
    "get_interpreter",
]
