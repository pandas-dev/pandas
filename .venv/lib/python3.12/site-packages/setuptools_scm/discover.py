from __future__ import annotations

import os

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Iterable
from typing import Iterator

from . import _entrypoints
from . import _log
from . import _types as _t
from ._config import Configuration

if TYPE_CHECKING:
    from ._entrypoints import im


log = _log.log.getChild("discover")


def walk_potential_roots(root: _t.PathT, search_parents: bool = True) -> Iterator[Path]:
    """
    Iterate though a path and each of its parents.
    :param root: File path.
    :param search_parents: If ``False`` the parents are not considered.
    """
    root = Path(root)
    yield root
    if search_parents:
        yield from root.parents


def match_entrypoint(root: _t.PathT, name: str) -> bool:
    """
    Consider a ``root`` as entry-point.
    :param root: File path.
    :param name: Subdirectory name.
    :return: ``True`` if a subdirectory ``name`` exits in ``root``.
    """

    if os.path.exists(os.path.join(root, name)):
        if not os.path.isabs(name):
            return True
        log.debug("ignoring bad ep %s", name)

    return False


# blocked entrypints from legacy plugins
_BLOCKED_EP_TARGETS = {"setuptools_scm_git_archive:parse"}


def iter_matching_entrypoints(
    root: _t.PathT, entrypoint: str, config: Configuration
) -> Iterable[im.EntryPoint]:
    """
    Consider different entry-points in ``root`` and optionally its parents.
    :param root: File path.
    :param entrypoint: Entry-point to consider.
    :param config: Configuration,
        read ``search_parent_directories``, write found parent to ``parent``.
    """

    log.debug("looking for ep %s in %s", entrypoint, root)

    for wd in walk_potential_roots(root, config.search_parent_directories):
        for ep in _entrypoints.entry_points(group=entrypoint):
            if ep.value in _BLOCKED_EP_TARGETS:
                continue
            if match_entrypoint(wd, ep.name):
                log.debug("found ep %s in %s", ep, wd)
                config.parent = wd
                yield ep
