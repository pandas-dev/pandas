from numba.cuda.testing import ensure_supported_ccs_initialized
from numba.testing import load_testsuite
import os


def load_tests(loader, tests, pattern):
    ensure_supported_ccs_initialized()
    return load_testsuite(loader, os.path.dirname(__file__))
