import sys
from importlib.metadata import version

from moto.utilities.distutils_version import LooseVersion

PYTHON_VERSION_INFO = sys.version_info
PYTHON_311 = sys.version_info >= (3, 11)
RESPONSES_VERSION = version("responses")
WERKZEUG_VERSION = version("werkzeug")


def is_responses_0_17_x() -> bool:
    return LooseVersion(RESPONSES_VERSION) >= LooseVersion("0.17.0")


def is_werkzeug_2_0_x_or_older() -> bool:
    return LooseVersion(WERKZEUG_VERSION) < LooseVersion("2.1.0")


def is_werkzeug_2_3_x() -> bool:
    return LooseVersion(WERKZEUG_VERSION) >= LooseVersion("2.3.0")
