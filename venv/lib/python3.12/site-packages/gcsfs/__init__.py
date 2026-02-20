import logging
import os

from ._version import get_versions

logger = logging.getLogger(__name__)
__version__ = get_versions()["version"]
del get_versions
from .core import GCSFileSystem
from .mapping import GCSMap

if os.getenv("GCSFS_EXPERIMENTAL_ZB_HNS_SUPPORT", "false").lower() in ("true", "1"):
    try:
        from .extended_gcsfs import ExtendedGcsFileSystem as GCSFileSystem

        logger.info(
            "gcsfs experimental features enabled via GCSFS_EXPERIMENTAL_ZB_HNS_SUPPORT."
        )
    except ImportError as e:
        logger.warning(
            f"GCSFS_EXPERIMENTAL_ZB_HNS_SUPPORT is set, but failed to import experimental features: {e}"
        )
        # Fallback to core GCSFileSystem, do not register here

# TODO: GCSMap still refers to the original GCSFileSystem. This will be
# addressed in a future update.
__all__ = ["GCSFileSystem", "GCSMap"]

from . import _version

__version__ = _version.get_versions()["version"]
