
import warnings
from distutils.version import LooseVersion

_NUMEXPR_INSTALLED = False

try:
    import numexpr as ne
    ver = ne.__version__
    _NUMEXPR_INSTALLED = ver >= LooseVersion('2.4.6')

    if not _NUMEXPR_INSTALLED:
        warnings.warn(
            "The installed version of numexpr {ver} is not supported "
            "in pandas and will be not be used\nThe minimum supported "
            "version is 2.4.6\n".format(ver=ver), UserWarning)

except ImportError:  # pragma: no cover
    pass

__all__ = ['_NUMEXPR_INSTALLED']
