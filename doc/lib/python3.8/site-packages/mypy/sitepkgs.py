from __future__ import print_function
"""This file is used to find the site packages of a Python executable, which may be Python 2.

This file MUST remain compatible with Python 2. Since we cannot make any assumptions about the
Python being executed, this module should not use *any* dependencies outside of the standard
library found in Python 2. This file is run each mypy run, so it should be kept as fast as
possible.
"""

if __name__ == '__main__':
    import sys
    sys.path = sys.path[1:]  # we don't want to pick up mypy.types

import site

MYPY = False
if MYPY:
    from typing import List


def getsitepackages():
    # type: () -> List[str]
    if hasattr(site, 'getusersitepackages') and hasattr(site, 'getsitepackages'):
        user_dir = site.getusersitepackages()
        return site.getsitepackages() + [user_dir]
    else:
        from distutils.sysconfig import get_python_lib
        return [get_python_lib()]


if __name__ == '__main__':
    print(repr(getsitepackages()))
