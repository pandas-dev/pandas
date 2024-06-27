# -*- coding: utf-8 -*-

import warnings
from numba.core.errors import NumbaPendingDeprecationWarning
# The pycc module requires setuptools.
try:
    import setuptools
except ImportError:
    msg = "The 'setuptools' package is required at runtime for pycc support."
    raise ImportError(msg)

# Public API
from .cc import CC
from .decorators import export, exportmany

# If use of anything is attempted through the `pycc` import path this warning
# will be shown.
__pycc_deprecation_doc_url = ("https://numba.readthedocs.io/en/stable/"
                              "reference/deprecation.html"
                              "#deprecation-of-the-numba-pycc-module")
__pycc_pending_deprecation_message = ("The 'pycc' module is pending "
                                      "deprecation. Replacement technology is "
                                      "being developed.\n\n"
                                      "Pending Deprecation in Numba 0.57.0. "
                                      "For more information please see: "
                                      f"{__pycc_deprecation_doc_url}")

_pend_dep = NumbaPendingDeprecationWarning(__pycc_pending_deprecation_message)
warnings.warn(_pend_dep, stacklevel=2)
