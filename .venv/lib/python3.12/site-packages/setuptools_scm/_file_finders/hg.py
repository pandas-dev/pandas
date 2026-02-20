from __future__ import annotations

import logging
import os
import subprocess

from .. import _types as _t
from .._file_finders import is_toplevel_acceptable
from .._file_finders import scm_find_files
from ..hg import run_hg
from ..integration import data_from_mime
from .pathtools import norm_real

log = logging.getLogger(__name__)


def _hg_toplevel(path: str) -> str | None:
    try:
        return run_hg(
            ["root"],
            cwd=(path or "."),
            check=True,
        ).parse_success(norm_real)
    except subprocess.CalledProcessError:
        # hg returned error, we are not in a mercurial repo
        return None
    except OSError:
        # hg command not found, probably
        return None


def _hg_ls_files_and_dirs(toplevel: str) -> tuple[set[str], set[str]]:
    hg_files: set[str] = set()
    hg_dirs = {toplevel}
    res = run_hg(["files"], cwd=toplevel)
    if res.returncode:
        return set(), set()
    for name in res.stdout.splitlines():
        name = os.path.normcase(name).replace("/", os.path.sep)
        fullname = os.path.join(toplevel, name)
        hg_files.add(fullname)
        dirname = os.path.dirname(fullname)
        while len(dirname) > len(toplevel) and dirname not in hg_dirs:
            hg_dirs.add(dirname)
            dirname = os.path.dirname(dirname)
    return hg_files, hg_dirs


def hg_find_files(path: str = "") -> list[str]:
    toplevel = _hg_toplevel(path)
    if not is_toplevel_acceptable(toplevel):
        return []
    assert toplevel is not None
    hg_files, hg_dirs = _hg_ls_files_and_dirs(toplevel)
    return scm_find_files(path, hg_files, hg_dirs)


def hg_archive_find_files(path: _t.PathT = "") -> list[str]:
    # This function assumes that ``path`` is obtained from a mercurial archive
    # and therefore all the files that should be ignored were already removed.
    archival = os.path.join(path, ".hg_archival.txt")
    if not os.path.exists(archival):
        return []

    data = data_from_mime(archival)

    if "node" not in data:
        # Ensure file is valid
        return []

    log.warning("hg archive detected - fallback to listing all files")
    return scm_find_files(path, set(), set(), force_all_files=True)
