from __future__ import annotations

import logging
import os

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import _types as _t
from . import Configuration
from .integration import data_from_mime
from .version import ScmVersion
from .version import meta
from .version import tag_to_version

log = logging.getLogger(__name__)

_UNKNOWN = "UNKNOWN"


def parse_pkginfo(root: _t.PathT, config: Configuration) -> ScmVersion | None:
    pkginfo = Path(root) / "PKG-INFO"
    log.debug("pkginfo %s", pkginfo)
    data = data_from_mime(pkginfo)
    version = data.get("Version", _UNKNOWN)
    if version != _UNKNOWN:
        return meta(version, preformatted=True, config=config)
    else:
        return None


def fallback_version(root: _t.PathT, config: Configuration) -> ScmVersion | None:
    if config.parentdir_prefix_version is not None:
        _, parent_name = os.path.split(os.path.abspath(root))
        if parent_name.startswith(config.parentdir_prefix_version):
            version = tag_to_version(
                parent_name[len(config.parentdir_prefix_version) :], config
            )
            if version is not None:
                return meta(str(version), preformatted=True, config=config)
    if config.fallback_version is not None:
        log.debug("FALLBACK %s", config.fallback_version)
        return meta(config.fallback_version, preformatted=True, config=config)
    return None
