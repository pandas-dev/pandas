"""Patches that are applied at runtime to the virtual environment."""

from __future__ import annotations

import contextlib
import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import types
    from collections.abc import Callable
    from importlib.machinery import ModuleSpec

VIRTUALENV_PATCH_FILE = os.path.abspath(__file__)


def patch_dist(dist: types.ModuleType) -> None:
    """Distutils allows user to configure some arguments via a configuration file: https://docs.python.org/3/install/index.html#distutils-configuration-files.

    Some of this arguments though don't make sense in context of the virtual environment files, let's fix them up.

    """
    # we cannot allow some install config as that would get packages installed outside of the virtual environment
    old_parse_config_files = dist.Distribution.parse_config_files

    def parse_config_files(self, *args: object, **kwargs: object) -> object:  # ruff:ignore[missing-type-function-argument]
        result = old_parse_config_files(self, *args, **kwargs)
        install = self.get_option_dict("install")

        if "prefix" in install:  # the prefix governs where to install the libraries
            install["prefix"] = VIRTUALENV_PATCH_FILE, os.path.abspath(sys.prefix)
        for base in ("purelib", "platlib", "headers", "scripts", "data"):
            key = f"install_{base}"
            if key in install:  # do not allow global configs to hijack venv paths
                install.pop(key, None)
        return result

    dist.Distribution.parse_config_files = parse_config_files


# Import hook that patches some modules to ignore configuration values that break package installation in case
# of virtual environments.
_DISTUTILS_PATCH = "distutils.dist", "setuptools.dist"
# https://docs.python.org/3/library/importlib.html#setting-up-an-importer


class _Finder:
    """A meta path finder that allows patching the imported distutils modules."""

    fullname = None

    # lock[0] is threading.Lock(), but initialized lazily to avoid importing threading very early at startup,
    # because there are gevent-based applications that need to be first to import threading by themselves.
    # See https://github.com/pypa/virtualenv/issues/1895 for details.
    lock = []  # ruff:ignore[mutable-class-default]

    def find_spec(self, fullname: str, path: object, target: object = None) -> ModuleSpec | None:  # ruff:ignore[unused-method-argument]
        # Guard against race conditions during file rewrite by checking if _DISTUTILS_PATCH is defined.
        # This can happen when the file is being overwritten while it's being imported by another process.
        # See https://github.com/pypa/virtualenv/issues/2969 for details.
        try:
            distutils_patch = _DISTUTILS_PATCH
        except NameError:
            return None
        if fullname in distutils_patch and self.fullname is None:
            # initialize lock[0] lazily
            if len(self.lock) == 0:
                import threading  # ruff:ignore[import-outside-top-level]

                lock = threading.Lock()
                # there is possibility that two threads T1 and T2 are simultaneously running into find_spec,
                # observing .lock as empty, and further going into hereby initialization. However due to the GIL,
                # list.append() operation is atomic and this way only one of the threads will "win" to put the lock
                # - that every thread will use - into .lock[0].
                # https://docs.python.org/3/faq/library.html#what-kinds-of-global-value-mutation-are-thread-safe
                self.lock.append(lock)

            from functools import partial  # ruff:ignore[import-outside-top-level]
            from importlib.util import find_spec  # ruff:ignore[import-outside-top-level]

            with self.lock[0]:
                self.fullname = fullname
                spec = None
                try:
                    spec = find_spec(fullname, path)  # ty: ignore[invalid-argument-type]
                finally:
                    self.fullname = None
                if spec is not None:
                    return self._patch_spec(spec, partial)
        return None

    def _patch_spec(self, spec: ModuleSpec, partial: Callable[..., object]) -> ModuleSpec:
        old = getattr(spec.loader, "exec_module", None)
        if old is not None and old is not self.exec_module:
            try:  # ruff:ignore[suppressible-exception]
                spec.loader.exec_module = partial(self.exec_module, old)  # ty: ignore[invalid-assignment]
            except AttributeError:
                pass
        return spec

    @staticmethod
    def exec_module(old: Callable[..., object], module: types.ModuleType) -> None:
        old(module)
        try:
            distutils_patch = _DISTUTILS_PATCH
        except NameError:
            return
        if module.__name__ in distutils_patch:
            with contextlib.suppress(NameError):
                patch_dist(module)


sys.meta_path.insert(0, _Finder())
