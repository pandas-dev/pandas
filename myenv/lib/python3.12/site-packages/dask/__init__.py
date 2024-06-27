from __future__ import annotations

from dask import config, datasets
from dask._version import get_versions
from dask.base import (
    annotate,
    compute,
    get_annotations,
    is_dask_collection,
    optimize,
    persist,
    visualize,
)
from dask.core import istask
from dask.delayed import delayed
from dask.local import get_sync as get

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
