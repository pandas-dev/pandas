# Copyright 2013 Hardcoded Software (http://www.hardcoded.net)

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

import sys

from send2trash.exceptions import TrashPermissionError  # noqa: F401

if sys.platform == "darwin":
    from send2trash.mac import send2trash
elif sys.platform == "win32":
    from send2trash.win import send2trash
else:
    try:
        # If we can use gio, let's use it
        from send2trash.plat_gio import send2trash
    except ImportError:
        # Oh well, let's fallback to our own Freedesktop trash implementation
        from send2trash.plat_other import send2trash  # noqa: F401
