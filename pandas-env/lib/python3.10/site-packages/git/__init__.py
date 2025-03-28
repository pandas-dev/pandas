# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

# @PydevCodeAnalysisIgnore

__all__ = [
    "Actor",
    "AmbiguousObjectName",
    "BadName",
    "BadObject",
    "BadObjectType",
    "BaseIndexEntry",
    "Blob",
    "BlobFilter",
    "BlockingLockFile",
    "CacheError",
    "CheckoutError",
    "CommandError",
    "Commit",
    "Diff",
    "DiffConstants",
    "DiffIndex",
    "Diffable",
    "FetchInfo",
    "Git",
    "GitCmdObjectDB",
    "GitCommandError",
    "GitCommandNotFound",
    "GitConfigParser",
    "GitDB",
    "GitError",
    "HEAD",
    "Head",
    "HookExecutionError",
    "INDEX",
    "IndexEntry",
    "IndexFile",
    "IndexObject",
    "InvalidDBRoot",
    "InvalidGitRepositoryError",
    "List",  # Deprecated - import this from `typing` instead.
    "LockFile",
    "NULL_TREE",
    "NoSuchPathError",
    "ODBError",
    "Object",
    "Optional",  # Deprecated - import this from `typing` instead.
    "ParseError",
    "PathLike",
    "PushInfo",
    "RefLog",
    "RefLogEntry",
    "Reference",
    "Remote",
    "RemoteProgress",
    "RemoteReference",
    "Repo",
    "RepositoryDirtyError",
    "RootModule",
    "RootUpdateProgress",
    "Sequence",  # Deprecated - import from `typing`, or `collections.abc` in 3.9+.
    "StageType",
    "Stats",
    "Submodule",
    "SymbolicReference",
    "TYPE_CHECKING",  # Deprecated - import this from `typing` instead.
    "Tag",
    "TagObject",
    "TagReference",
    "Tree",
    "TreeModifier",
    "Tuple",  # Deprecated - import this from `typing` instead.
    "Union",  # Deprecated - import this from `typing` instead.
    "UnmergedEntriesError",
    "UnsafeOptionError",
    "UnsafeProtocolError",
    "UnsupportedOperation",
    "UpdateProgress",
    "WorkTreeRepositoryUnsupported",
    "refresh",
    "remove_password_if_present",
    "rmtree",
    "safe_decode",
    "to_hex_sha",
]

__version__ = '3.1.44'

from typing import Any, List, Optional, Sequence, TYPE_CHECKING, Tuple, Union

if TYPE_CHECKING:
    from types import ModuleType

import warnings

from gitdb.util import to_hex_sha

from git.exc import (
    AmbiguousObjectName,
    BadName,
    BadObject,
    BadObjectType,
    CacheError,
    CheckoutError,
    CommandError,
    GitCommandError,
    GitCommandNotFound,
    GitError,
    HookExecutionError,
    InvalidDBRoot,
    InvalidGitRepositoryError,
    NoSuchPathError,
    ODBError,
    ParseError,
    RepositoryDirtyError,
    UnmergedEntriesError,
    UnsafeOptionError,
    UnsafeProtocolError,
    UnsupportedOperation,
    WorkTreeRepositoryUnsupported,
)
from git.types import PathLike

try:
    from git.compat import safe_decode  # @NoMove
    from git.config import GitConfigParser  # @NoMove
    from git.objects import (  # @NoMove
        Blob,
        Commit,
        IndexObject,
        Object,
        RootModule,
        RootUpdateProgress,
        Submodule,
        TagObject,
        Tree,
        TreeModifier,
        UpdateProgress,
    )
    from git.refs import (  # @NoMove
        HEAD,
        Head,
        RefLog,
        RefLogEntry,
        Reference,
        RemoteReference,
        SymbolicReference,
        Tag,
        TagReference,
    )
    from git.diff import (  # @NoMove
        INDEX,
        NULL_TREE,
        Diff,
        DiffConstants,
        DiffIndex,
        Diffable,
    )
    from git.db import GitCmdObjectDB, GitDB  # @NoMove
    from git.cmd import Git  # @NoMove
    from git.repo import Repo  # @NoMove
    from git.remote import FetchInfo, PushInfo, Remote, RemoteProgress  # @NoMove
    from git.index import (  # @NoMove
        BaseIndexEntry,
        BlobFilter,
        CheckoutError,
        IndexEntry,
        IndexFile,
        StageType,
        # NOTE: This tells type checkers what util resolves to. We delete it, and it is
        # really resolved by __getattr__, which warns. See below on what to use instead.
        util,
    )
    from git.util import (  # @NoMove
        Actor,
        BlockingLockFile,
        LockFile,
        Stats,
        remove_password_if_present,
        rmtree,
    )
except GitError as _exc:
    raise ImportError("%s: %s" % (_exc.__class__.__name__, _exc)) from _exc


def _warned_import(message: str, fullname: str) -> "ModuleType":
    import importlib

    warnings.warn(message, DeprecationWarning, stacklevel=3)
    return importlib.import_module(fullname)


def _getattr(name: str) -> Any:
    # TODO: If __version__ is made dynamic and lazily fetched, put that case right here.

    if name == "util":
        return _warned_import(
            "The expression `git.util` and the import `from git import util` actually "
            "reference git.index.util, and not the git.util module accessed in "
            '`from git.util import XYZ` or `sys.modules["git.util"]`. This potentially '
            "confusing behavior is currently preserved for compatibility, but may be "
            "changed in the future and should not be relied on.",
            fullname="git.index.util",
        )

    for names, prefix in (
        ({"head", "log", "reference", "symbolic", "tag"}, "git.refs"),
        ({"base", "fun", "typ"}, "git.index"),
    ):
        if name not in names:
            continue

        fullname = f"{prefix}.{name}"

        return _warned_import(
            f"{__name__}.{name} is a private alias of {fullname} and subject to "
            f"immediate removal. Use {fullname} instead.",
            fullname=fullname,
        )

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if not TYPE_CHECKING:
    # NOTE: The expression `git.util` gives git.index.util and `from git import util`
    # imports git.index.util, NOT git.util. It may not be feasible to change this until
    # the next major version, to avoid breaking code inadvertently relying on it.
    #
    # - If git.index.util *is* what you want, use (or import from) that, to avoid
    #   confusion.
    #
    # - To use the "real" git.util module, write `from git.util import ...`, or if
    #   necessary access it as `sys.modules["git.util"]`.
    #
    # Note also that `import git.util` technically imports the "real" git.util... but
    # the *expression* `git.util` after doing so is still git.index.util!
    #
    # (This situation differs from that of other indirect-submodule imports that are
    # unambiguously non-public and subject to immediate removal. Here, the public
    # git.util module, though different, makes less discoverable that the expression
    # `git.util` refers to a non-public attribute of the git module.)
    #
    # This had originally come about by a wildcard import. Now that all intended imports
    # are explicit, the intuitive but potentially incompatible binding occurs due to the
    # usual rules for Python submodule bindings. So for now we replace that binding with
    # git.index.util, delete that, and let __getattr__ handle it and issue a warning.
    #
    # For the same runtime behavior, it would be enough to forgo importing util, and
    # delete util as created naturally; __getattr__ would behave the same. But type
    # checkers would not know what util refers to when accessed as an attribute of git.
    del util

    # This is "hidden" to preserve static checking for undefined/misspelled attributes.
    __getattr__ = _getattr

# { Initialize git executable path

GIT_OK = None


def refresh(path: Optional[PathLike] = None) -> None:
    """Convenience method for setting the git executable path.

    :param path:
        Optional path to the Git executable. If not absolute, it is resolved
        immediately, relative to the current directory.

    :note:
        The `path` parameter is usually omitted and cannot be used to specify a custom
        command whose location is looked up in a path search on each call. See
        :meth:`Git.refresh <git.cmd.Git.refresh>` for details on how to achieve this.

    :note:
        This calls :meth:`Git.refresh <git.cmd.Git.refresh>` and sets other global
        configuration according to the effect of doing so. As such, this function should
        usually be used instead of using :meth:`Git.refresh <git.cmd.Git.refresh>` or
        :meth:`FetchInfo.refresh <git.remote.FetchInfo.refresh>` directly.

    :note:
        This function is called automatically, with no arguments, at import time.
    """
    global GIT_OK
    GIT_OK = False

    if not Git.refresh(path=path):
        return
    if not FetchInfo.refresh():  # noqa: F405
        return  # type: ignore[unreachable]

    GIT_OK = True


try:
    refresh()
except Exception as _exc:
    raise ImportError("Failed to initialize: {0}".format(_exc)) from _exc

# } END initialize git executable path
