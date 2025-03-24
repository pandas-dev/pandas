import sys
from numba.core.utils import _RedirectSubpackage
from numba.core import config

if config.USE_LEGACY_TYPE_SYSTEM:
    sys.modules[__name__] = _RedirectSubpackage(locals(),
                                                "numba.cpython.old_tupleobj")
else:
    sys.modules[__name__] = _RedirectSubpackage(locals(),
                                                "numba.cpython.new_tupleobj")
