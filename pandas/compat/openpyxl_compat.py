"""
Detect incompatible version of OpenPyXL

GH7169
"""

from distutils.version import LooseVersion

start_ver = '1.6.1'
stop_ver = '2.0.0'

def is_compat():
    """
    Detect whether the installed version of openpyxl is supported
    Returns True/False if openpyxl is installed, None otherwise
    """
    try:
        import openpyxl
    except ImportError:
        return None

    ver = LooseVersion(openpyxl.__version__)
    if ver < LooseVersion(start_ver) or LooseVersion(stop_ver) <= ver:
        return False

    return True
