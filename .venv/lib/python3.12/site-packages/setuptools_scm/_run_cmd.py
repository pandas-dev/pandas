from __future__ import annotations

import os
import shlex
import subprocess
import textwrap
import warnings

from typing import TYPE_CHECKING
from typing import Callable
from typing import Final
from typing import Mapping
from typing import Sequence
from typing import TypeVar
from typing import overload

from . import _log
from . import _types as _t

if TYPE_CHECKING:
    BaseCompletedProcess = subprocess.CompletedProcess[str]
else:
    BaseCompletedProcess = subprocess.CompletedProcess

# pick 40 seconds
# unfortunately github CI for windows sometimes needs
# up to 30 seconds to start a command


def _get_timeout(env: Mapping[str, str]) -> int:
    return int(env.get("SETUPTOOLS_SCM_SUBPROCESS_TIMEOUT") or 40)


BROKEN_TIMEOUT: Final[int] = _get_timeout(os.environ)

log = _log.log.getChild("run_cmd")

PARSE_RESULT = TypeVar("PARSE_RESULT")
T = TypeVar("T")


class CompletedProcess(BaseCompletedProcess):
    @classmethod
    def from_raw(
        cls, input: BaseCompletedProcess, strip: bool = True
    ) -> CompletedProcess:
        return cls(
            args=input.args,
            returncode=input.returncode,
            stdout=input.stdout.strip() if strip and input.stdout else input.stdout,
            stderr=input.stderr.strip() if strip and input.stderr else input.stderr,
        )

    @overload
    def parse_success(
        self,
        parse: Callable[[str], PARSE_RESULT],
        default: None = None,
        error_msg: str | None = None,
    ) -> PARSE_RESULT | None: ...

    @overload
    def parse_success(
        self,
        parse: Callable[[str], PARSE_RESULT],
        default: T,
        error_msg: str | None = None,
    ) -> PARSE_RESULT | T: ...

    def parse_success(
        self,
        parse: Callable[[str], PARSE_RESULT],
        default: T | None = None,
        error_msg: str | None = None,
    ) -> PARSE_RESULT | T | None:
        if self.returncode:
            if error_msg:
                log.warning("%s %s", error_msg, self)
            return default
        else:
            return parse(self.stdout)


KEEP_GIT_ENV = (
    "GIT_CEILING_DIRECTORIES",
    "GIT_EXEC_PATH",
    "GIT_SSH",
    "GIT_SSH_COMMAND",
    "GIT_AUTHOR_DATE",
    "GIT_COMMITTER_DATE",
)


def no_git_env(env: Mapping[str, str]) -> dict[str, str]:
    # adapted from pre-commit
    # Too many bugs dealing with environment variables and GIT:
    # https://github.com/pre-commit/pre-commit/issues/300
    # In git 2.6.3 (maybe others), git exports GIT_WORK_TREE while running
    # pre-commit hooks
    # In git 1.9.1 (maybe others), git exports GIT_DIR and GIT_INDEX_FILE
    # while running pre-commit hooks in submodules.
    # GIT_DIR: Causes git clone to clone wrong thing
    # GIT_INDEX_FILE: Causes 'error invalid object ...' during commit
    for k, v in env.items():
        if k.startswith("GIT_"):
            log.debug("%s: %s", k, v)
    return {
        k: v for k, v in env.items() if not k.startswith("GIT_") or k in KEEP_GIT_ENV
    }


def avoid_pip_isolation(env: Mapping[str, str]) -> dict[str, str]:
    """
    pip build isolation can break Mercurial
    (see https://github.com/pypa/pip/issues/10635)

    pip uses PYTHONNOUSERSITE and a path in PYTHONPATH containing "pip-build-env-".
    """
    new_env = {k: v for k, v in env.items() if k != "PYTHONNOUSERSITE"}
    if "PYTHONPATH" not in new_env:
        return new_env

    new_env["PYTHONPATH"] = os.pathsep.join(
        [
            path
            for path in new_env["PYTHONPATH"].split(os.pathsep)
            if "-build-env-" not in path
        ]
    )
    return new_env


def ensure_stripped_str(str_or_bytes: str | bytes) -> str:
    if isinstance(str_or_bytes, str):
        return str_or_bytes.strip()
    else:
        return str_or_bytes.decode("utf-8", "surrogateescape").strip()


def run(
    cmd: _t.CMD_TYPE,
    cwd: _t.PathT,
    *,
    strip: bool = True,
    trace: bool = True,
    timeout: int | None = None,
    check: bool = False,
) -> CompletedProcess:
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    else:
        cmd = [os.fspath(x) for x in cmd]
    cmd_4_trace = " ".join(map(_unsafe_quote_for_display, cmd))
    log.debug("at %s\n    $ %s ", cwd, cmd_4_trace)
    if timeout is None:
        timeout = BROKEN_TIMEOUT
    res = subprocess.run(
        cmd,
        capture_output=True,
        cwd=os.fspath(cwd),
        env=dict(
            avoid_pip_isolation(no_git_env(os.environ)),
            # os.environ,
            # try to disable i18n, but still allow UTF-8 encoded text.
            LC_ALL="C.UTF-8",
            LANGUAGE="",
            HGPLAIN="1",
        ),
        text=True,
        encoding="utf-8",
        timeout=timeout,
    )

    res = CompletedProcess.from_raw(res, strip=strip)
    if trace:
        if res.stdout:
            log.debug("out:\n%s", textwrap.indent(res.stdout, "    "))
        if res.stderr:
            log.debug("err:\n%s", textwrap.indent(res.stderr, "    "))
        if res.returncode:
            log.debug("ret: %s", res.returncode)
    if check:
        res.check_returncode()
    return res


def _unsafe_quote_for_display(item: _t.PathT) -> str:
    # give better results than shlex.join in our cases
    text = os.fspath(item)
    return text if all(c not in text for c in " {[:") else f'"{text}"'


def has_command(
    name: str, args: Sequence[str] = ["version"], warn: bool = True
) -> bool:
    try:
        p = run([name, *args], cwd=".")
        if p.returncode != 0:
            log.error("Command '%s' returned non-zero. This is stderr:", name)
            log.error(p.stderr)
    except OSError as e:
        log.warning("command %s missing: %s", name, e)
        res = False
    except subprocess.TimeoutExpired as e:
        log.warning("command %s timed out %s", name, e)
        res = False

    else:
        res = not p.returncode
    if not res and warn:
        warnings.warn(f"{name!r} was not found", category=RuntimeWarning)
    return res


class CommandNotFoundError(LookupError, FileNotFoundError):
    pass


def require_command(name: str) -> None:
    if not has_command(name, warn=False):
        raise CommandNotFoundError(name)
