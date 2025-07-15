import sys
from numba.core.utils import _RedirectSubpackage
from numba.core import config

if config.USE_LEGACY_TYPE_SYSTEM:
    sys.modules[__name__] = _RedirectSubpackage(locals(),
                                                "numba.np.old_arraymath")
else:
    sys.modules[__name__] = _RedirectSubpackage(locals(),
                                                "numba.np.new_arraymath")
