import sys
from numba.core.utils import _RedirectSubpackage
sys.modules[__name__] = _RedirectSubpackage(locals(), "numba.core.types")
