# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

"""Module with our own gitdb implementation - it uses the git command."""

__all__ = ["GitCmdObjectDB", "GitDB"]

from gitdb.base import OInfo, OStream
from gitdb.db import GitDB, LooseObjectDB
from gitdb.exc import BadObject

from git.util import bin_to_hex, hex_to_bin
from git.exc import GitCommandError

# typing-------------------------------------------------

from typing import TYPE_CHECKING

from git.types import PathLike

if TYPE_CHECKING:
    from git.cmd import Git

# --------------------------------------------------------


class GitCmdObjectDB(LooseObjectDB):
    """A database representing the default git object store, which includes loose
    objects, pack files and an alternates file.

    It will create objects only in the loose object database.
    """

    def __init__(self, root_path: PathLike, git: "Git") -> None:
        """Initialize this instance with the root and a git command."""
        super().__init__(root_path)
        self._git = git

    def info(self, binsha: bytes) -> OInfo:
        """Get a git object header (using git itself)."""
        hexsha, typename, size = self._git.get_object_header(bin_to_hex(binsha))
        return OInfo(hex_to_bin(hexsha), typename, size)

    def stream(self, binsha: bytes) -> OStream:
        """Get git object data as a stream supporting ``read()`` (using git itself)."""
        hexsha, typename, size, stream = self._git.stream_object_data(bin_to_hex(binsha))
        return OStream(hex_to_bin(hexsha), typename, size, stream)

    # { Interface

    def partial_to_complete_sha_hex(self, partial_hexsha: str) -> bytes:
        """
        :return:
            Full binary 20 byte sha from the given partial hexsha

        :raise gitdb.exc.AmbiguousObjectName:

        :raise gitdb.exc.BadObject:

        :note:
            Currently we only raise :exc:`~gitdb.exc.BadObject` as git does not
            communicate ambiguous objects separately.
        """
        try:
            hexsha, _typename, _size = self._git.get_object_header(partial_hexsha)
            return hex_to_bin(hexsha)
        except (GitCommandError, ValueError) as e:
            raise BadObject(partial_hexsha) from e
        # END handle exceptions

    # } END interface
