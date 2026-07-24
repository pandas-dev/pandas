from __future__ import annotations

import logging
import os
import sys
from contextlib import suppress
from pathlib import Path
from typing import TYPE_CHECKING, Final

from platformdirs import user_data_path

from ._compat import fs_path_id
from ._py_info import PythonInfo
from ._py_spec import PythonSpec

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable, Iterator, Mapping, Sequence

    from ._cache import PyInfoCache

_LOGGER: Final[logging.Logger] = logging.getLogger(__name__)
IS_WIN: Final[bool] = sys.platform == "win32"


def get_interpreter(
    key: str | Sequence[str],
    try_first_with: Iterable[str] | None = None,
    cache: PyInfoCache | None = None,
    env: Mapping[str, str] | None = None,
    predicate: Callable[[PythonInfo], bool] | None = None,
) -> PythonInfo | None:
    """
    Find a Python interpreter matching *key*.

    Iterates over one or more specification strings and returns the first interpreter that satisfies the spec and passes
    the optional *predicate*.

    :param key: interpreter specification string(s) — an absolute path, a version (``3.12``), an implementation prefix
        (``cpython3.12``), or a
        `version specifier <https://packaging.python.org/en/latest/specifications/version-specifiers/>`_
        (``>=3.10``). When a sequence is given each entry is tried in order.
    :param try_first_with: executables to probe before the normal discovery search.
    :param cache: interpreter metadata cache; when ``None`` results are not cached.
    :param env: environment mapping for ``PATH`` lookup; defaults to :data:`os.environ`.
    :param predicate: optional callback applied after an interpreter matches the spec. Return ``True`` to accept the
        interpreter, ``False`` to skip it and continue searching.
    :return: the first matching interpreter, or ``None`` if no match is found.
    """
    specs = [key] if isinstance(key, str) else key
    for spec_str in specs:
        if result := _find_interpreter(spec_str, try_first_with or (), cache, env, predicate):
            return result
    return None


def iter_interpreters(
    key: str | Sequence[str] | None = None,
    try_first_with: Iterable[str] | None = None,
    cache: PyInfoCache | None = None,
    env: Mapping[str, str] | None = None,
    predicate: Callable[[PythonInfo], bool] | None = None,
) -> Iterator[PythonInfo]:
    """
    Yield every interpreter on the system that satisfies *key*.

    Iteration order is discovery order: ``try_first_with`` paths first, then the running interpreter, then ``PATH``
    (left to right), then UV-managed installs. Results are deduplicated by the resolved real path of the underlying
    system interpreter, so symlinked aliases (``/bin`` vs ``/usr/bin``) and venvs that symlink to a base interpreter
    collapse to a single entry. Callers that want a different ordering should sort the result.

    :param key: interpreter specification — same syntax as :func:`get_interpreter`. ``None`` enumerates every Python
        implementation python-discovery knows about (see :data:`KNOWN_IMPLEMENTATIONS`).
    :param try_first_with: executables to probe before the normal discovery search.
    :param cache: interpreter metadata cache; when ``None`` results are not cached. Strongly recommended for
        enumeration, which interrogates every candidate as a subprocess on a cold cache.
    :param env: environment mapping for ``PATH`` lookup; defaults to :data:`os.environ`.
    :param predicate: optional filter applied after the spec match; return ``True`` to include the interpreter.
    """
    if key is None:
        keys: tuple[str | None, ...] = (None,)
    elif isinstance(key, str):
        keys = (key,)
    else:
        keys = tuple(key)
    first_with = tuple(try_first_with or ())
    env_map = os.environ if env is None else env
    seen: set[str] = set()
    for spec_str in keys:
        yield from _iter_for_spec(spec_str, first_with, cache, env_map, predicate, seen)


def _iter_for_spec(  # ruff:ignore[too-many-arguments, too-many-positional-arguments]
    spec_str: str | None,
    try_first_with: tuple[str, ...],
    cache: PyInfoCache | None,
    env: Mapping[str, str],
    predicate: Callable[[PythonInfo], bool] | None,
    seen: set[str],
) -> Iterator[PythonInfo]:
    if spec_str is None:
        spec = PythonSpec("", None, None, None, None, None, None)
        wide = True
    else:
        spec = PythonSpec.from_string_spec(spec_str)
        wide = False
    for interpreter, impl_must_match in propose_interpreters(
        spec, try_first_with, cache, env, all_implementations=wide
    ):
        if interpreter is None:
            continue
        if (anchor := interpreter.system_executable or interpreter.executable) is None:
            continue
        if (real_path := os.path.realpath(anchor)) in seen:
            continue
        if not interpreter.satisfies(spec, impl_must_match=impl_must_match):
            continue
        if predicate is not None and not predicate(interpreter):
            continue
        seen.add(real_path)
        yield interpreter


def _find_interpreter(
    key: str,
    try_first_with: Iterable[str],
    cache: PyInfoCache | None = None,
    env: Mapping[str, str] | None = None,
    predicate: Callable[[PythonInfo], bool] | None = None,
) -> PythonInfo | None:
    spec = PythonSpec.from_string_spec(key)
    _LOGGER.info("find interpreter for spec %r", spec)
    proposed_paths: set[tuple[str | None, bool]] = set()
    env = os.environ if env is None else env
    for interpreter, impl_must_match in propose_interpreters(spec, try_first_with, cache, env):
        if interpreter is None:  # pragma: no cover
            continue
        proposed_key = interpreter.system_executable, impl_must_match
        if proposed_key in proposed_paths:
            continue
        _LOGGER.info("proposed %s", interpreter)
        if interpreter.satisfies(spec, impl_must_match=impl_must_match) and (
            predicate is None or predicate(interpreter)
        ):
            _LOGGER.debug("accepted %s", interpreter)
            return interpreter
        proposed_paths.add(proposed_key)
    return None


def _check_exe(path: str, tested_exes: set[str]) -> str | None:
    """Resolve *path* to an absolute path and return it if not yet tested, otherwise ``None``."""
    try:
        os.lstat(path)
    except OSError:
        return None
    resolved = str(Path(path).resolve())
    exe_id = fs_path_id(resolved)
    if exe_id in tested_exes:
        return None
    tested_exes.add(exe_id)
    return str(Path(path).absolute())


def _is_new_exe(exe_raw: str, tested_exes: set[str]) -> bool:
    """Return ``True`` and register *exe_raw* if it hasn't been tested yet."""
    exe_id = fs_path_id(exe_raw)
    if exe_id in tested_exes:
        return False
    tested_exes.add(exe_id)
    return True


def propose_interpreters(
    spec: PythonSpec,
    try_first_with: Iterable[str],
    cache: PyInfoCache | None = None,
    env: Mapping[str, str] | None = None,
    *,
    all_implementations: bool = False,
) -> Generator[tuple[PythonInfo | None, bool], None, None]:
    """
    Yield ``(interpreter, impl_must_match)`` candidates for *spec*.

    :param spec: the parsed interpreter specification to match against.
    :param try_first_with: executable paths to probe before the standard search.
    :param cache: interpreter metadata cache; when ``None`` results are not cached.
    :param env: environment mapping for ``PATH`` lookup; defaults to :data:`os.environ`.
    :param all_implementations: when ``True`` and *spec* does not constrain the implementation, also surface
        non-CPython binaries on ``PATH`` and under UV's install directory. Used by enumeration APIs.
    """
    env = os.environ if env is None else env
    tested_exes: set[str] = set()
    if spec.is_abs and spec.path is not None:
        if exe_raw := _check_exe(spec.path, tested_exes):  # pragma: no branch # first exe always new
            yield PythonInfo.from_exe(exe_raw, cache, env=env), True
        return

    yield from _propose_explicit(spec, try_first_with, cache, env, tested_exes)
    if spec.path is not None and spec.is_abs:  # pragma: no cover # relative spec.path is never abs
        return
    yield from _propose_from_path(spec, cache, env, tested_exes, all_implementations=all_implementations)
    yield from _propose_from_uv(cache, env, all_implementations=all_implementations)


def _propose_explicit(
    spec: PythonSpec,
    try_first_with: Iterable[str],
    cache: PyInfoCache | None,
    env: Mapping[str, str],
    tested_exes: set[str],
) -> Generator[tuple[PythonInfo | None, bool], None, None]:
    for py_exe in try_first_with:
        if exe_raw := _check_exe(str(Path(py_exe).resolve()), tested_exes):
            yield PythonInfo.from_exe(exe_raw, cache, env=env), True

    if spec.path is not None:
        if exe_raw := _check_exe(spec.path, tested_exes):  # pragma: no branch
            yield PythonInfo.from_exe(exe_raw, cache, env=env), True
    else:
        yield from _propose_current_and_windows(spec, cache, env, tested_exes)


def _propose_current_and_windows(
    spec: PythonSpec,
    cache: PyInfoCache | None,
    env: Mapping[str, str],
    tested_exes: set[str],
) -> Generator[tuple[PythonInfo | None, bool], None, None]:
    current_python = PythonInfo.current_system(cache)
    if _is_new_exe(str(current_python.executable), tested_exes):
        yield current_python, True

    if IS_WIN:  # pragma: win32 cover
        from ._windows import propose_interpreters as win_propose  # ruff:ignore[import-outside-top-level]

        for interpreter in win_propose(spec, cache, env):
            if _is_new_exe(str(interpreter.executable), tested_exes):
                yield interpreter, True


def _propose_from_path(
    spec: PythonSpec,
    cache: PyInfoCache | None,
    env: Mapping[str, str],
    tested_exes: set[str],
    *,
    all_implementations: bool = False,
) -> Generator[tuple[PythonInfo | None, bool], None, None]:
    find_candidates = path_exe_finder(spec, all_implementations=all_implementations)
    for pos, path in enumerate(get_paths(env)):
        _LOGGER.debug(LazyPathDump(pos, path, env))
        for exe, impl_must_match in find_candidates(path):
            exe_raw = str(exe)
            if resolved := _resolve_shim(exe_raw, env):
                _LOGGER.debug("resolved shim %s to %s", exe_raw, resolved)
                exe_raw = resolved
            if not _is_new_exe(exe_raw, tested_exes):
                continue
            interpreter = PathPythonInfo.from_exe(exe_raw, cache, raise_on_error=False, env=env)
            if interpreter is not None:
                yield interpreter, impl_must_match


def _propose_from_uv(
    cache: PyInfoCache | None,
    env: Mapping[str, str],
    *,
    all_implementations: bool = False,
) -> Generator[tuple[PythonInfo | None, bool], None, None]:
    if uv_python_dir := os.getenv("UV_PYTHON_INSTALL_DIR"):
        uv_python_path = Path(uv_python_dir).expanduser()
    elif xdg_data_home := os.getenv("XDG_DATA_HOME"):
        uv_python_path = Path(xdg_data_home).expanduser() / "uv" / "python"
    else:
        uv_python_path = user_data_path("uv") / "python"

    patterns: list[str] = ["*/bin/python", "*/python.exe"]
    if all_implementations:
        patterns.extend(("*/bin/pypy*", "*/bin/graalpy", "*/pypy*.exe", "*/bin/graalpy.exe"))
    seen_uv_paths: set[str] = set()
    for pattern in patterns:
        for exe_path in uv_python_path.glob(pattern):
            resolved = str(Path(exe_path).resolve())
            if resolved in seen_uv_paths:
                continue
            seen_uv_paths.add(resolved)
            if interpreter := PathPythonInfo.from_exe(str(exe_path), cache, raise_on_error=False, env=env):
                yield interpreter, True


def get_paths(env: Mapping[str, str]) -> Generator[Path, None, None]:
    path = env.get("PATH", None)
    if path is None:
        try:
            path = os.confstr("CS_PATH")
        except (AttributeError, ValueError):  # pragma: no cover # Windows only (no confstr)
            path = os.defpath
    if path:
        for entry in map(Path, path.split(os.pathsep)):
            with suppress(OSError):
                if entry.is_dir() and next(entry.iterdir(), None):
                    yield entry


class LazyPathDump:
    def __init__(self, pos: int, path: Path, env: Mapping[str, str]) -> None:
        self.pos = pos
        self.path = path
        self.env = env

    def __repr__(self) -> str:
        content = f"discover PATH[{self.pos}]={self.path}"
        if self.env.get("_VIRTUALENV_DEBUG"):
            content += " with =>"
            for file_path in self.path.iterdir():
                try:
                    if not self._is_executable(file_path):
                        continue
                except OSError:
                    pass
                content += " "
                content += file_path.name
        return content

    def _is_executable(self, file_path: Path) -> bool:
        if file_path.is_dir():
            return False
        if IS_WIN:  # pragma: win32 cover
            pathext = self.env.get("PATHEXT", ".COM;.EXE;.BAT;.CMD").split(";")
            return any(file_path.name.upper().endswith(ext) for ext in pathext)
        return bool(file_path.stat().st_mode & os.X_OK)


def path_exe_finder(
    spec: PythonSpec, *, all_implementations: bool = False
) -> Callable[[Path], Generator[tuple[Path, bool], None, None]]:
    """Given a spec, return a function that can be called on a path to find all matching files in it."""
    pat = spec.generate_re(windows=sys.platform == "win32", all_implementations=all_implementations)
    direct = spec.str_spec
    if sys.platform == "win32":  # pragma: win32 cover
        direct = f"{direct}.exe"

    def path_exes(path: Path) -> Generator[tuple[Path, bool], None, None]:
        direct_path = path / direct
        if direct_path.exists():
            yield direct_path, False

        for exe in path.iterdir():
            match = pat.fullmatch(exe.name)
            if match:
                yield exe.absolute(), match["impl"] == "python"

    return path_exes


def _resolve_shim(exe_path: str, env: Mapping[str, str]) -> str | None:
    """Resolve a version-manager shim to the actual Python binary."""
    for shims_dir_env, versions_path in _VERSION_MANAGER_LAYOUTS:
        if root := env.get(shims_dir_env):
            shims_dir = os.path.join(root, "shims")
            if os.path.dirname(exe_path) == shims_dir:
                exe_name = os.path.basename(exe_path)
                versions_dir = os.path.join(root, *versions_path)
                return _resolve_shim_to_binary(exe_name, versions_dir, env)
    return None


_VERSION_MANAGER_LAYOUTS: list[tuple[str, tuple[str, ...]]] = [
    ("PYENV_ROOT", ("versions",)),
    ("MISE_DATA_DIR", ("installs", "python")),
    ("ASDF_DATA_DIR", ("installs", "python")),
]


def _resolve_shim_to_binary(exe_name: str, versions_dir: str, env: Mapping[str, str]) -> str | None:
    for version in _active_versions(env):
        resolved = os.path.join(versions_dir, version, "bin", exe_name)
        if Path(resolved).is_file() and os.access(resolved, os.X_OK):
            return resolved
    return None


def _active_versions(env: Mapping[str, str]) -> Generator[str, None, None]:
    """Yield active Python version strings by reading version-manager configuration."""
    if pyenv_version := env.get("PYENV_VERSION"):
        yield from pyenv_version.split(":")
        return
    if versions := _read_python_version_file(Path.cwd()):
        yield from versions
        return
    if (pyenv_root := env.get("PYENV_ROOT")) and (
        versions := _read_python_version_file(os.path.join(pyenv_root, "version"), search_parents=False)
    ):
        yield from versions


def _read_python_version_file(start: str | Path, *, search_parents: bool = True) -> list[str] | None:
    """Read a ``.python-version`` file, optionally searching parent directories."""
    current = start
    while True:
        candidate = os.path.join(current, ".python-version") if Path(current).is_dir() else current
        if Path(candidate).is_file():
            with Path(candidate).open(encoding="utf-8") as fh:
                if versions := [v for line in fh if (v := line.strip()) and not v.startswith("#")]:
                    return versions
        if not search_parents:
            return None
        parent = Path(current).parent
        if parent == current:
            return None
        current = parent


class PathPythonInfo(PythonInfo):
    """python info from path."""


__all__ = [
    "LazyPathDump",
    "PathPythonInfo",
    "get_interpreter",
    "get_paths",
    "iter_interpreters",
    "propose_interpreters",
]
