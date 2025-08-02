###################################################################
#  Numexpr - Fast numerical array expression evaluator for NumPy.
#
#      License: MIT
#      Author:  See AUTHORS.txt
#
#  See LICENSE.txt and LICENSES/*.txt for details about copyright and
#  rights to use.
####################################################################

import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "thread_unsafe: mark test as unsafe for parallel execution"
    )
