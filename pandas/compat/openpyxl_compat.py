"""
Detect incompatible version of OpenPyXL

GH7169
"""

from distutils.version import LooseVersion

start_ver = '1.6.1'
stop_ver = '2.0.0'


def is_compat():
    """Detect whether the installed version of openpyxl is supported.

    Returns
    -------
    compat : bool
        ``True`` if openpyxl is installed and is between versions 1.6.1 and
        2.0.0, ``False`` otherwise.
    """
    import openpyxl
    ver = LooseVersion(openpyxl.__version__)
    return LooseVersion(start_ver) < ver <= LooseVersion(stop_ver)
