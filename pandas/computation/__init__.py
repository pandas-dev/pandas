
import warnings
from distutils.version import LooseVersion

_NUMEXPR_INSTALLED = False

try:
    import numexpr as ne
    ver = ne.__version__
    _NUMEXPR_INSTALLED = ver >= LooseVersion('2.1')

    # we specifically disallow 2.4.4 as
    # has some hard-to-diagnose bugs
    if ver == LooseVersion('2.4.4'):
        _NUMEXPR_INSTALLED = False
        warnings.warn(
            "The installed version of numexpr {ver} is not supported "
            "in pandas and will be not be used\n".format(ver=ver),
            UserWarning)

    elif not _NUMEXPR_INSTALLED:
        warnings.warn(
            "The installed version of numexpr {ver} is not supported "
            "in pandas and will be not be used\nThe minimum supported "
            "version is 2.1\n".format(ver=ver), UserWarning)

except ImportError:  # pragma: no cover
    pass

__all__ = ['_NUMEXPR_INSTALLED']
