from __future__ import annotations

from dask import config, datasets
from dask._expr import Expr, HLGExpr, LLGExpr, SingletonExpr

try:
    # Backwards compatibility with versioneer
    from dask._version import __commit_id__ as __git_revision__
    from dask._version import __version__
except ImportError:  # pragma: no cover
    __git_revision__ = "unknown"
    __version__ = "unknown"
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
