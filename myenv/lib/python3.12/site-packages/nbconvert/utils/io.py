"""io-related utilities"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import codecs
import errno
import os
import random
import shutil
import sys
from typing import Any, Dict


def unicode_std_stream(stream="stdout"):
    """Get a wrapper to write unicode to stdout/stderr as UTF-8.

    This ignores environment variables and default encodings, to reliably write
    unicode to stdout or stderr.

    ::

        unicode_std_stream().write(u'ł@e¶ŧ←')
    """
    assert stream in ("stdout", "stderr")
    stream = getattr(sys, stream)

    try:
        stream_b = stream.buffer
    except AttributeError:
        # sys.stdout has been replaced - use it directly
        return stream

    return codecs.getwriter("utf-8")(stream_b)


def unicode_stdin_stream():
    """Get a wrapper to read unicode from stdin as UTF-8.

    This ignores environment variables and default encodings, to reliably read unicode from stdin.

    ::

        totreat = unicode_stdin_stream().read()
    """
    stream = sys.stdin
    try:
        stream_b = stream.buffer
    except AttributeError:
        return stream

    return codecs.getreader("utf-8")(stream_b)


class FormatSafeDict(Dict[Any, Any]):
    """Format a dictionary safely."""

    def __missing__(self, key):
        """Handle missing value."""
        return "{" + key + "}"


try:
    ENOLINK = errno.ENOLINK
except AttributeError:
    ENOLINK = 1998


def link(src, dst):
    """Hard links ``src`` to ``dst``, returning 0 or errno.

    Note that the special errno ``ENOLINK`` will be returned if ``os.link`` isn't
    supported by the operating system.
    """

    if not hasattr(os, "link"):
        return ENOLINK
    link_errno = 0
    try:
        os.link(src, dst)
    except OSError as e:
        link_errno = e.errno
    return link_errno


def link_or_copy(src, dst):
    """Attempts to hardlink ``src`` to ``dst``, copying if the link fails.

    Attempts to maintain the semantics of ``shutil.copy``.

    Because ``os.link`` does not overwrite files, a unique temporary file
    will be used if the target already exists, then that file will be moved
    into place.
    """

    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))

    link_errno = link(src, dst)
    if link_errno == errno.EEXIST:
        if os.stat(src).st_ino == os.stat(dst).st_ino:
            # dst is already a hard link to the correct file, so we don't need
            # to do anything else. If we try to link and rename the file
            # anyway, we get duplicate files - see http://bugs.python.org/issue21876
            return

        new_dst = dst + f"-temp-{random.randint(1, 16**4):04X}"  # noqa: S311
        try:
            link_or_copy(src, new_dst)
        except BaseException:
            try:
                os.remove(new_dst)
            except OSError:
                pass
            raise
        os.rename(new_dst, dst)
    elif link_errno != 0:
        # Either link isn't supported, or the filesystem doesn't support
        # linking, or 'src' and 'dst' are on different filesystems.
        shutil.copy(src, dst)
