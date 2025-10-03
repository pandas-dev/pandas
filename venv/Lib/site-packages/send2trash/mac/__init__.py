# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

from platform import mac_ver
from sys import version_info

# NOTE: version of pyobjc only supports python >= 3.6 and 10.9+
macos_ver = tuple(int(part) for part in mac_ver()[0].split("."))
if version_info >= (3, 6) and macos_ver >= (10, 9):
    try:
        from send2trash.mac.modern import send2trash
    except ImportError:
        # Try to fall back to ctypes version, although likely problematic still
        from send2trash.mac.legacy import send2trash
else:
    # Just use the old version otherwise
    from send2trash.mac.legacy import send2trash  # noqa: F401
