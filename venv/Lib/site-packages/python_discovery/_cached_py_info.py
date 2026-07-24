"""Acquire Python information via subprocess interrogation with multi-level caching."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pkgutil
import secrets
import subprocess  # ruff:ignore[suspicious-subprocess-import]
import sys
import tempfile
from collections import OrderedDict
from contextlib import contextmanager
from pathlib import Path
from shlex import quote
from subprocess import Popen, TimeoutExpired  # ruff:ignore[suspicious-subprocess-import]
from typing import TYPE_CHECKING, Final

from ._cache import NoOpCache
from ._py_info import PythonInfo

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping

    from ._cache import ContentStore, PyInfoCache


_CACHE: OrderedDict[Path, PythonInfo | Exception] = OrderedDict()
_CACHE[Path(sys.executable)] = PythonInfo()
_LOGGER: Final[logging.Logger] = logging.getLogger(__name__)


def from_exe(  # ruff:ignore[too-many-arguments]
    cls: type[PythonInfo],
    cache: PyInfoCache | None,
    exe: str,
    env: Mapping[str, str] | None = None,
    *,
    raise_on_error: bool = True,
    ignore_cache: bool = False,
) -> PythonInfo | None:
    env = os.environ if env is None else env
    result = _get_from_cache(cls, cache, exe, env, ignore_cache=ignore_cache)
    if isinstance(result, Exception):
        if raise_on_error:
            raise result
        _LOGGER.info("%s", result)
        result = None
    return result


def _get_from_cache(
    cls: type[PythonInfo],
    cache: PyInfoCache | None,
    exe: str,
    env: Mapping[str, str],
    *,
    ignore_cache: bool = True,
) -> PythonInfo | Exception:
    exe_path = Path(exe)
    if not ignore_cache and exe_path in _CACHE:
        result = _CACHE[exe_path]
    else:
        py_info = _get_via_file_cache(cls, cache, exe_path, exe, env)
        result = _CACHE[exe_path] = py_info
    if isinstance(result, PythonInfo):
        result.executable = exe
    return result


def _get_via_file_cache(
    cls: type[PythonInfo],
    cache: PyInfoCache | None,
    path: Path,
    exe: str,
    env: Mapping[str, str],
) -> PythonInfo | Exception:
    path_text = str(path)
    try:
        path_modified = path.stat().st_mtime
    except OSError:
        path_modified = -1
    py_info_script = Path(Path(__file__).resolve()).parent / "_py_info.py"
    try:
        py_info_hash: str | None = hashlib.sha256(py_info_script.read_bytes()).hexdigest()
    except OSError:
        py_info_hash = None

    resolved_cache = cache if cache is not None else NoOpCache()
    py_info: PythonInfo | None = None
    py_info_store = resolved_cache.py_info(path)
    with py_info_store.locked():
        if py_info_store.exists() and (data := py_info_store.read()) is not None:
            of_path, of_st_mtime = data.get("path"), data.get("st_mtime")
            of_content, of_hash = data.get("content"), data.get("hash")
            if (
                of_path == path_text
                and of_st_mtime == path_modified
                and of_hash == py_info_hash
                and isinstance(of_content, dict)
            ):
                py_info = _load_cached_py_info(cls, py_info_store, of_content)
            else:
                py_info_store.remove()
        if py_info is None:
            failure, py_info = _run_subprocess(cls, exe, env)
            if failure is not None:
                _LOGGER.debug("first subprocess attempt failed for %s (%s), retrying", exe, failure)
                failure, py_info = _run_subprocess(cls, exe, env)
            if failure is not None:
                return failure
            if py_info is not None:
                py_info_store.write({
                    "st_mtime": path_modified,
                    "path": path_text,
                    "content": py_info.to_dict(),
                    "hash": py_info_hash,
                })
    if py_info is None:
        msg = f"{exe} failed to produce interpreter info"
        return RuntimeError(msg)
    return py_info


def _load_cached_py_info(
    cls: type[PythonInfo],
    py_info_store: ContentStore,
    content: dict,
) -> PythonInfo | None:
    try:
        py_info = cls.from_dict(content.copy())
    except (KeyError, TypeError):
        py_info_store.remove()
        return None
    if (sys_exe := py_info.system_executable) is not None and not Path(sys_exe).exists():
        py_info_store.remove()
        return None
    return py_info


COOKIE_LENGTH: Final[int] = 32


def gen_cookie() -> str:
    return secrets.token_hex(COOKIE_LENGTH // 2)


@contextmanager
def _resolve_py_info_script() -> Generator[Path]:
    py_info_script = Path(Path(__file__).resolve()).parent / "_py_info.py"
    if py_info_script.is_file():
        yield py_info_script
    else:
        data = pkgutil.get_data(__package__ or __name__, "_py_info.py")
        if data is None:
            msg = "cannot locate _py_info.py for subprocess interrogation"
            raise FileNotFoundError(msg)
        fd, tmp = tempfile.mkstemp(suffix=".py")
        try:
            os.write(fd, data)
            os.close(fd)
            yield Path(tmp)
        finally:
            Path(tmp).unlink()


def _extract_between_cookies(out: str, start_cookie: str, end_cookie: str) -> tuple[str, str, int, int]:
    """Extract payload between reversed cookie markers, forwarding any surrounding output to stdout."""
    raw_out = out
    out_starts = out.find(start_cookie[::-1])
    if out_starts > -1:
        if pre_cookie := out[:out_starts]:
            sys.stdout.write(pre_cookie)
        out = out[out_starts + COOKIE_LENGTH :]
    out_ends = out.find(end_cookie[::-1])
    if out_ends > -1:
        if post_cookie := out[out_ends + COOKIE_LENGTH :]:
            sys.stdout.write(post_cookie)
        out = out[:out_ends]
    return out, raw_out, out_starts, out_ends


def _run_subprocess(
    cls: type[PythonInfo],
    exe: str,
    env: Mapping[str, str],
) -> tuple[Exception | None, PythonInfo | None]:
    start_cookie = gen_cookie()
    end_cookie = gen_cookie()
    timeout = float(env.get("PY_DISCOVERY_TIMEOUT", "15"))
    with _resolve_py_info_script() as py_info_script:
        cmd = [exe, str(py_info_script), start_cookie, end_cookie]
        env = dict(env)
        env.pop("__PYVENV_LAUNCHER__", None)
        env["PYTHONUTF8"] = "1"
        _LOGGER.debug("get interpreter info via cmd: %s", LogCmd(cmd))
        try:
            process = Popen(  # ruff:ignore[subprocess-without-shell-equals-true]
                cmd,
                universal_newlines=True,
                stdin=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                env=env,
                encoding="utf-8",
                errors="backslashreplace",
            )
            out, err = process.communicate(timeout=timeout)
            code = process.returncode
        except TimeoutExpired:
            process.kill()
            process.communicate()
            out, err, code = "", "timed out", -1
        except OSError as os_error:
            out, err, code = "", os_error.strerror, os_error.errno
    if code != 0:
        msg = f"{exe} with code {code}{f' out: {out!r}' if out else ''}{f' err: {err!r}' if err else ''}"
        return RuntimeError(f"failed to query {msg}"), None
    out, raw_out, out_starts, out_ends = _extract_between_cookies(out, start_cookie, end_cookie)
    try:
        result = cls.from_json(out)
        result.executable = exe
    except json.JSONDecodeError as exc:
        _LOGGER.warning(
            "subprocess %s returned invalid JSON; raw stdout %d chars, start cookie %s, end cookie %s, "
            "parsed output %d chars: %r",
            exe,
            len(raw_out),
            "found" if out_starts > -1 else "missing",
            "found" if out_ends > -1 else "missing",
            len(out),
            out[:200] if out else "<empty>",
        )
        msg = f"{exe} returned invalid JSON (exit code {code}){f', stderr: {err!r}' if err else ''}"
        failure = RuntimeError(msg)
        failure.__cause__ = exc
        return failure, None
    return None, result


class LogCmd:
    def __init__(self, cmd: list[str], env: Mapping[str, str] | None = None) -> None:
        self.cmd = cmd
        self.env = env

    def __repr__(self) -> str:
        cmd_repr = " ".join(quote(str(c)) for c in self.cmd)
        if self.env is not None:
            cmd_repr = f"{cmd_repr} env of {self.env!r}"
        return cmd_repr


def clear(cache: PyInfoCache) -> None:
    cache.py_info_clear()
    _CACHE.clear()


__all__ = [
    "LogCmd",
    "clear",
    "from_exe",
]
