from __future__ import annotations

import abc
import fnmatch
from itertools import chain
from operator import methodcaller as method
from pathlib import Path
from textwrap import dedent

from virtualenv.create.describe import Python3Supports
from virtualenv.create.via_global_ref.builtin.ref import ExePathRefToDest, PathRefToDest, RefWhen
from virtualenv.create.via_global_ref.store import is_store_python
from virtualenv.util.path import copy as copy_path
from virtualenv.util.path import ensure_dir

from .common import CPython, CPythonPosix, CPythonWindows, is_mac_os_framework, is_macos_brew


class CPython3(CPython, Python3Supports, abc.ABC):
    """CPython 3 or later."""


class CPython3Posix(CPythonPosix, CPython3):
    @classmethod
    def can_describe(cls, interpreter):
        return (
            is_mac_os_framework(interpreter) is False
            and is_macos_brew(interpreter) is False
            and super().can_describe(interpreter)
        )

    @classmethod
    def sources(cls, interpreter):
        yield from super().sources(interpreter)
        if shared_lib := cls._shared_libpython(interpreter):
            yield PathRefToDest(shared_lib, dest=cls._to_lib, when=RefWhen.COPY)

    @classmethod
    def _to_lib(cls, creator, src):
        return creator.dest / "lib" / src.name

    @classmethod
    def _shared_libpython(cls, interpreter):
        if not interpreter.sysconfig_vars.get("Py_ENABLE_SHARED"):
            return None
        if not (instsoname := interpreter.sysconfig_vars.get("INSTSONAME")):
            return None
        if not (libdir := interpreter.sysconfig_vars.get("LIBDIR")):
            return None
        if not (lib_path := Path(libdir) / instsoname).exists():
            return None
        return lib_path

    def install_venv_shared_libs(self, venv_creator):
        if venv_creator.symlinks:
            return
        if not (shared_lib := self._shared_libpython(venv_creator.interpreter)):
            return
        dest = venv_creator.dest / "lib" / shared_lib.name
        ensure_dir(dest.parent)
        copy_path(shared_lib, dest)

    def env_patch_text(self):
        text = super().env_patch_text()
        if self.pyvenv_launch_patch_active(self.interpreter):
            text += dedent(
                """
                # for https://github.com/python/cpython/pull/9516, see https://github.com/pypa/virtualenv/issues/1704
                import os
                if "__PYVENV_LAUNCHER__" in os.environ:
                    del os.environ["__PYVENV_LAUNCHER__"]
                """,
            )
        return text

    @classmethod
    def pyvenv_launch_patch_active(cls, interpreter):
        ver = interpreter.version_info
        return interpreter.platform == "darwin" and ((3, 7, 8) > ver >= (3, 7) or (3, 8, 3) > ver >= (3, 8))


class CPython3Windows(CPythonWindows, CPython3):
    """CPython 3 on Windows."""

    @classmethod
    def setup_meta(cls, interpreter):
        if is_store_python(interpreter):  # store python is not supported here
            return None
        return super().setup_meta(interpreter)

    @classmethod
    def sources(cls, interpreter):
        if cls.has_shim(interpreter):
            refs = cls.executables(interpreter)
        else:
            refs = chain(
                cls.executables(interpreter),
                cls.dll_and_pyd(interpreter),
                cls.python_zip(interpreter),
            )
        yield from refs

    @classmethod
    def executables(cls, interpreter):
        sources = super().sources(interpreter)
        if interpreter.version_info >= (3, 13):
            t_suffix = "t" if interpreter.free_threaded else ""
            updated_sources = []
            for ref in sources:
                if ref.src.name == "python.exe":
                    launcher_path = ref.src.with_name(f"venvlauncher{t_suffix}.exe")
                    if launcher_path.exists():
                        new_ref = ExePathRefToDest(
                            launcher_path, dest=ref.dest, targets=[ref.base, *ref.aliases], must=ref.must, when=ref.when
                        )
                        updated_sources.append(new_ref)
                        continue
                elif ref.src.name == "pythonw.exe":
                    w_launcher_path = ref.src.with_name(f"venvwlauncher{t_suffix}.exe")
                    if w_launcher_path.exists():
                        new_ref = ExePathRefToDest(
                            w_launcher_path,
                            dest=ref.dest,
                            targets=[ref.base, *ref.aliases],
                            must=ref.must,
                            when=ref.when,
                        )
                        updated_sources.append(new_ref)
                        continue
                updated_sources.append(ref)
            return updated_sources
        return sources

    @classmethod
    def has_shim(cls, interpreter):
        return interpreter.version_info.minor >= 7 and cls.shim(interpreter) is not None  # noqa: PLR2004

    @classmethod
    def shim(cls, interpreter):
        root = Path(interpreter.system_stdlib) / "venv" / "scripts" / "nt"
        if interpreter.version_info >= (3, 13):
            # https://github.com/python/cpython/issues/112984
            t_suffix = "t" if interpreter.free_threaded else ""
            exe_name = f"venvlauncher{t_suffix}.exe"
        else:
            exe_name = "python.exe"
        if (shim := root / exe_name).exists():
            return shim
        return None

    @classmethod
    def host_python(cls, interpreter):
        if cls.has_shim(interpreter):
            # starting with CPython 3.7 Windows ships with a venvlauncher.exe that avoids the need for dll/pyd copies
            # it also means the wrapper must be copied to avoid bugs such as https://bugs.python.org/issue42013
            return cls.shim(interpreter)
        return super().host_python(interpreter)

    @classmethod
    def dll_and_pyd(cls, interpreter):
        folders = [Path(interpreter.system_executable).parent]

        # May be missing on some Python hosts.
        # See https://github.com/pypa/virtualenv/issues/2368
        dll_folder = Path(interpreter.system_prefix) / "DLLs"
        if dll_folder.is_dir():
            folders.append(dll_folder)

        for folder in folders:
            for file in folder.iterdir():
                if file.suffix in {".pyd", ".dll"}:
                    # Skip pywin32 DLLs to avoid conflicts with pywin32 installation
                    # pywin32 has its own post-install that places DLLs in site-packages/pywin32_system32
                    # See https://github.com/pypa/virtualenv/issues/2662
                    if cls._is_pywin32_dll(file.name):
                        continue
                    yield PathRefToDest(file, cls.to_bin)

    @classmethod
    def _is_pywin32_dll(cls, filename):
        """Check if a DLL file belongs to pywin32."""
        # pywin32 DLLs follow patterns like: pywintypes39.dll, pythoncom39.dll
        name_lower = filename.lower()
        return name_lower.startswith(("pywintypes", "pythoncom"))

    @classmethod
    def python_zip(cls, interpreter):
        """``python{VERSION}.zip`` contains compiled ``*.pyc`` std lib packages, where ``VERSION`` is ``py_version_nodot`` var from the ``sysconfig`` module.

        See https://docs.python.org/3/using/windows.html#the-embeddable-package, ``discovery.py_info.PythonInfo`` class
        (interpreter), and ``python -m sysconfig`` output.

        The embeddable Python distribution for Windows includes ``python{VERSION}.zip`` and ``python{VERSION}._pth``
        files. User can move/rename the zip file and edit ``sys.path`` by editing the ``_pth`` file. Here the
        ``pattern`` is used only for the default zip file name.

        """
        pattern = f"*python{interpreter.version_nodot}.zip"
        matches = fnmatch.filter(interpreter.path, pattern)
        matched_paths = map(Path, matches)
        existing_paths = filter(method("exists"), matched_paths)
        if (path := next(existing_paths, None)) is not None:
            yield PathRefToDest(path, cls.to_bin)


__all__ = [
    "CPython3",
    "CPython3Posix",
    "CPython3Windows",
]
