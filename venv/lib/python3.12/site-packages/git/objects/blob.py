# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

__all__ = ["Blob"]

from mimetypes import guess_type
import sys

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

from . import base


class Blob(base.IndexObject):
    """A Blob encapsulates a git blob object.

    See :manpage:`gitglossary(7)` on "blob":
    https://git-scm.com/docs/gitglossary#def_blob_object
    """

    DEFAULT_MIME_TYPE = "text/plain"
    type: Literal["blob"] = "blob"

    # Valid blob modes
    executable_mode = 0o100755
    file_mode = 0o100644
    link_mode = 0o120000

    __slots__ = ()

    @property
    def mime_type(self) -> str:
        """
        :return:
            String describing the mime type of this file (based on the filename)

        :note:
            Defaults to ``text/plain`` in case the actual file type is unknown.
        """
        guesses = None
        if self.path:
            guesses = guess_type(str(self.path))
        return guesses and guesses[0] or self.DEFAULT_MIME_TYPE
