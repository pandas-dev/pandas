from logging import Logger

from .decorator import cross_origin as cross_origin
from .extension import CORS as CORS
from .version import __version__ as __version__

rootlogger: Logger

__all__ = ["CORS", "__version__", "cross_origin"]
