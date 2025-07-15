import sys
from numba.core.utils import _RedirectSubpackage
from numba.core import config

if config.USE_LEGACY_TYPE_SYSTEM:
    sys.modules[__name__] = \
        _RedirectSubpackage(locals(),
                            "numba.np.random.old_random_methods")
else:
    sys.modules[__name__] = \
        _RedirectSubpackage(locals(),
                            "numba.np.random.new_random_methods")
