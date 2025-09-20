from xdist._version import version as __version__
from xdist.plugin import get_xdist_worker_id
from xdist.plugin import is_xdist_controller
from xdist.plugin import is_xdist_master
from xdist.plugin import is_xdist_worker


__all__ = [
    "__version__",
    "get_xdist_worker_id",
    "is_xdist_controller",
    "is_xdist_master",
    "is_xdist_worker",
]
