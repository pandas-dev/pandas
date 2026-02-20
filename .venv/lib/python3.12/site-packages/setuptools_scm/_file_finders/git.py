from __future__ import annotations

import logging
import os
import subprocess
import tarfile

from typing import IO

from .. import _types as _t
from .._run_cmd import run as _run
from ..integration import data_from_mime
from . import is_toplevel_acceptable
from . import scm_find_files
from .pathtools import norm_real

log = logging.getLogger(__name__)


def _git_toplevel(path: str) -> str | None:
    try:
        cwd = os.path.abspath(path or ".")
        res = _run(["git", "rev-parse", "HEAD"], cwd=cwd)
        if res.returncode:
            # BAIL if there is no commit
            log.error("listing git files failed - pretending there aren't any")
            return None
        res = _run(
            ["git", "rev-parse", "--show-prefix"],
            cwd=cwd,
        )
        if res.returncode:
            return None
        out = res.stdout[:-1]  # remove the trailing pathsep
        if not out:
            out = cwd
        else:
            # Here, ``out`` is a relative path to root of git.
            # ``cwd`` is absolute path to current working directory.
            # the below method removes the length of ``out`` from
            # ``cwd``, which gives the git toplevel
            from .._compat import strip_path_suffix

            out = strip_path_suffix(cwd, out, f"cwd={cwd!r}\nout={out!r}")
        log.debug("find files toplevel %s", out)
        return norm_real(out)
    except subprocess.CalledProcessError:
        # git returned error, we are not in a git repo
        return None
    except OSError:
        # git command not found, probably
        return None


def _git_interpret_archive(fd: IO[bytes], toplevel: str) -> tuple[set[str], set[str]]:
    with tarfile.open(fileobj=fd, mode="r|*") as tf:
        git_files = set()
        git_dirs = {toplevel}
        for member in tf.getmembers():
            name = os.path.normcase(member.name).replace("/", os.path.sep)
            if member.type == tarfile.DIRTYPE:
                git_dirs.add(name)
            else:
                git_files.add(name)
        return git_files, git_dirs


def _git_ls_files_and_dirs(toplevel: str) -> tuple[set[str], set[str]]:
    # use git archive instead of git ls-file to honor
    # export-ignore git attribute

    cmd = ["git", "archive", "--prefix", toplevel + os.path.sep, "HEAD"]
    log.info("running %s", " ".join(str(x) for x in cmd))
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, cwd=toplevel, stderr=subprocess.DEVNULL
    )
    assert proc.stdout is not None
    try:
        try:
            return _git_interpret_archive(proc.stdout, toplevel)
        finally:
            # ensure we avoid resource warnings by cleaning up the process
            proc.stdout.close()
            proc.terminate()
            # Wait for process to actually terminate and be reaped
            try:
                proc.wait(timeout=5)  # Add timeout to avoid hanging
            except subprocess.TimeoutExpired:
                log.warning("git archive process did not terminate gracefully, killing")
                proc.kill()
                proc.wait()
    except Exception:
        # proc.wait() already called in finally block, check if it failed
        if proc.returncode != 0:
            log.error("listing git files failed - pretending there aren't any")
        return set(), set()


def git_find_files(path: _t.PathT = "") -> list[str]:
    toplevel = _git_toplevel(os.fspath(path))
    if not is_toplevel_acceptable(toplevel):
        return []
    fullpath = norm_real(path)
    if not fullpath.startswith(toplevel):
        log.warning("toplevel mismatch computed %s vs resolved %s ", toplevel, fullpath)
    git_files, git_dirs = _git_ls_files_and_dirs(toplevel)
    return scm_find_files(path, git_files, git_dirs)


def git_archive_find_files(path: _t.PathT = "") -> list[str]:
    # This function assumes that ``path`` is obtained from a git archive
    # and therefore all the files that should be ignored were already removed.
    archival = os.path.join(path, ".git_archival.txt")
    if not os.path.exists(archival):
        return []

    data = data_from_mime(archival)

    if "$Format" in data.get("node", ""):
        # Substitutions have not been performed, so not a reliable archive
        return []

    log.warning("git archive detected - fallback to listing all files")
    return scm_find_files(path, set(), set(), force_all_files=True)
