import errno


class TrashPermissionError(PermissionError):
    """A permission error specific to a trash directory.

    Raising this error indicates that permissions prevent us efficiently
    trashing a file, although we might still have permission to delete it.
    This is *not* used when permissions prevent removing the file itself:
    that will be raised as a regular PermissionError (OSError on Python 2).

    Application code that catches this may try to simply delete the file,
    or prompt the user to decide, or (on Freedesktop platforms), move it to
    'home trash' as a fallback. This last option probably involves copying the
    data between partitions, devices, or network drives, so we don't do it as
    a fallback.
    """

    def __init__(self, filename):
        PermissionError.__init__(self, errno.EACCES, "Permission denied", filename)
