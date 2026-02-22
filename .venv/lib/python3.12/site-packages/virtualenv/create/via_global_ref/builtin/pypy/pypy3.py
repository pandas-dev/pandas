from __future__ import annotations

import abc
from pathlib import Path

from virtualenv.create.describe import PosixSupports, Python3Supports, WindowsSupports
from virtualenv.create.via_global_ref.builtin.ref import PathRefToDest

from .common import PyPy


class PyPy3(PyPy, Python3Supports, abc.ABC):
    @classmethod
    def exe_stem(cls):
        return "pypy3"

    @classmethod
    def exe_names(cls, interpreter):
        return super().exe_names(interpreter) | {"pypy"}


class PyPy3Posix(PyPy3, PosixSupports):
    """PyPy 3 on POSIX."""

    @classmethod
    def _shared_libs(cls, python_dir):
        # glob for libpypy3-c.so, libpypy3-c.dylib, libpypy3.9-c.so ...
        return python_dir.glob("libpypy3*.*")

    def to_lib(self, src):
        return self.dest / "lib" / src.name

    @classmethod
    def sources(cls, interpreter):
        yield from super().sources(interpreter)
        # PyPy >= 3.8 supports a standard prefix installation, where older versions always used a portable/development
        # style installation. If this is a standard prefix installation, skip the below:
        if interpreter.system_prefix == "/usr":
            return
        # Also copy/symlink anything under prefix/lib, which, for "portable" PyPy builds, includes the tk,tcl runtime
        # and a number of shared objects. In distro-specific builds or on conda this should be empty (on PyPy3.8+ it
        # will, like on CPython, hold the stdlib).
        host_lib = Path(interpreter.system_prefix) / "lib"
        stdlib = Path(interpreter.system_stdlib)
        if host_lib.exists() and host_lib.is_dir():
            if (deps_file := host_lib / "PYPY_PORTABLE_DEPS.txt").exists():
                for line in deps_file.read_text(encoding="utf-8").splitlines():
                    dep = line.strip()
                    if dep and (path := host_lib / dep).exists():
                        yield PathRefToDest(path, dest=cls.to_lib)
            else:
                for path in host_lib.iterdir():
                    if stdlib == path:
                        continue
                    yield PathRefToDest(path, dest=cls.to_lib)


class Pypy3Windows(PyPy3, WindowsSupports):
    """PyPy 3 on Windows."""

    @property
    def less_v37(self):
        return self.interpreter.version_info.minor < 7  # noqa: PLR2004

    @classmethod
    def _shared_libs(cls, python_dir):
        # PyPy does not use a PEP 397 launcher, so all DLLs from the interpreter directory are needed for the venv
        yield from python_dir.glob("*.dll")


__all__ = [
    "PyPy3",
    "PyPy3Posix",
    "Pypy3Windows",
]
