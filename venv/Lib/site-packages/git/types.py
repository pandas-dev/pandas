# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

import os
import sys
from typing import (
    Any,
    Callable,
    Dict,
    List,
    NoReturn,
    Optional,
    Sequence as Sequence,
    Tuple,
    TYPE_CHECKING,
    TypeVar,
    Union,
)
import warnings

if sys.version_info >= (3, 8):
    from typing import (
        Literal,
        Protocol,
        SupportsIndex as SupportsIndex,
        TypedDict,
        runtime_checkable,
    )
else:
    from typing_extensions import (
        Literal,
        Protocol,
        SupportsIndex as SupportsIndex,
        TypedDict,
        runtime_checkable,
    )

if TYPE_CHECKING:
    from git.objects import Commit, Tree, TagObject, Blob
    from git.repo import Repo

PathLike = Union[str, "os.PathLike[str]"]
"""A :class:`str` (Unicode) based file or directory path."""

TBD = Any
"""Alias of :class:`~typing.Any`, when a type hint is meant to become more specific."""

_T = TypeVar("_T")
"""Type variable used internally in GitPython."""

AnyGitObject = Union["Commit", "Tree", "TagObject", "Blob"]
"""Union of the :class:`~git.objects.base.Object`-based types that represent actual git
object types.

As noted in :class:`~git.objects.base.Object`, which has further details, these are:

* :class:`Blob <git.objects.blob.Blob>`
* :class:`Tree <git.objects.tree.Tree>`
* :class:`Commit <git.objects.commit.Commit>`
* :class:`TagObject <git.objects.tag.TagObject>`

Those GitPython classes represent the four git object types, per
:manpage:`gitglossary(7)`:

* "blob": https://git-scm.com/docs/gitglossary#def_blob_object
* "tree object": https://git-scm.com/docs/gitglossary#def_tree_object
* "commit object": https://git-scm.com/docs/gitglossary#def_commit_object
* "tag object": https://git-scm.com/docs/gitglossary#def_tag_object

For more general information on git objects and their types as git understands them:

* "object": https://git-scm.com/docs/gitglossary#def_object
* "object type": https://git-scm.com/docs/gitglossary#def_object_type

:note:
    See also the :class:`Tree_ish` and :class:`Commit_ish` unions.
"""

Tree_ish = Union["Commit", "Tree", "TagObject"]
"""Union of :class:`~git.objects.base.Object`-based types that are typically tree-ish.

See :manpage:`gitglossary(7)` on "tree-ish":
https://git-scm.com/docs/gitglossary#def_tree-ish

:note:
    :class:`~git.objects.tree.Tree` and :class:`~git.objects.commit.Commit` are the
    classes whose instances are all tree-ish. This union includes them, but also
    :class:`~git.objects.tag.TagObject`, only **most** of whose instances are tree-ish.
    Whether a particular :class:`~git.objects.tag.TagObject` peels (recursively
    dereferences) to a tree or commit, rather than a blob, can in general only be known
    at runtime. In practice, git tag objects are nearly always used for tagging commits,
    and such tags are tree-ish because commits are tree-ish.

:note:
    See also the :class:`AnyGitObject` union of all four classes corresponding to git
    object types.
"""

Commit_ish = Union["Commit", "TagObject"]
"""Union of :class:`~git.objects.base.Object`-based types that are typically commit-ish.

See :manpage:`gitglossary(7)` on "commit-ish":
https://git-scm.com/docs/gitglossary#def_commit-ish

:note:
    :class:`~git.objects.commit.Commit` is the only class whose instances are all
    commit-ish. This union type includes :class:`~git.objects.commit.Commit`, but also
    :class:`~git.objects.tag.TagObject`, only **most** of whose instances are
    commit-ish. Whether a particular :class:`~git.objects.tag.TagObject` peels
    (recursively dereferences) to a commit, rather than a tree or blob, can in general
    only be known at runtime. In practice, git tag objects are nearly always used for
    tagging commits, and such tags are of course commit-ish.

:note:
    See also the :class:`AnyGitObject` union of all four classes corresponding to git
    object types.
"""

GitObjectTypeString = Literal["commit", "tag", "blob", "tree"]
"""Literal strings identifying git object types and the
:class:`~git.objects.base.Object`-based types that represent them.

See the :attr:`Object.type <git.objects.base.Object.type>` attribute. These are its
values in :class:`~git.objects.base.Object` subclasses that represent git objects. These
literals therefore correspond to the types in the :class:`AnyGitObject` union.

These are the same strings git itself uses to identify its four object types.
See :manpage:`gitglossary(7)` on "object type":
https://git-scm.com/docs/gitglossary#def_object_type
"""

if TYPE_CHECKING:
    Lit_commit_ish = Literal["commit", "tag"]
"""Deprecated. Type of literal strings identifying typically-commitish git object types.

Prior to a bugfix, this type had been defined more broadly. Any usage is in practice
ambiguous and likely to be incorrect. This type has therefore been made a static type
error to appear in annotations. It is preserved, with a deprecated status, to avoid
introducing runtime errors in code that refers to it, but it should not be used.

Instead of this type:

* For the type of the string literals associated with :class:`Commit_ish`, use
  ``Literal["commit", "tag"]`` or create a new type alias for it. That is equivalent to
  this type as currently defined (but usable in statically checked type annotations).

* For the type of all four string literals associated with :class:`AnyGitObject`, use
  :class:`GitObjectTypeString`. That is equivalent to the old definition of this type
  prior to the bugfix (and is also usable in statically checked type annotations).
"""


def _getattr(name: str) -> Any:
    if name != "Lit_commit_ish":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    warnings.warn(
        "Lit_commit_ish is deprecated. It is currently defined as "
        '`Literal["commit", "tag"]`, which should be used in its place if desired. It '
        'had previously been defined as `Literal["commit", "tag", "blob", "tree"]`, '
        "covering all four git object type strings including those that are never "
        "commit-ish. For that, use the GitObjectTypeString type instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return Literal["commit", "tag"]


if not TYPE_CHECKING:  # Preserve static checking for undefined/misspelled attributes.
    __getattr__ = _getattr


def __dir__() -> List[str]:
    return [*globals(), "Lit_commit_ish"]


# Config_levels ---------------------------------------------------------

Lit_config_levels = Literal["system", "global", "user", "repository"]
"""Type of literal strings naming git configuration levels.

These strings relate to which file a git configuration variable is in.
"""

ConfigLevels_Tup = Tuple[Literal["system"], Literal["user"], Literal["global"], Literal["repository"]]
"""Static type of a tuple of the four strings representing configuration levels."""

# Progress parameter type alias -----------------------------------------

CallableProgress = Optional[Callable[[int, Union[str, float], Union[str, float, None], str], None]]
"""General type of a function or other callable used as a progress reporter for cloning.

This is the type of a function or other callable that reports the progress of a clone,
when passed as a ``progress`` argument to :meth:`Repo.clone <git.repo.base.Repo.clone>`
or :meth:`Repo.clone_from <git.repo.base.Repo.clone_from>`.

:note:
    Those :meth:`~git.repo.base.Repo.clone` and :meth:`~git.repo.base.Repo.clone_from`
    methods also accept :meth:`~git.util.RemoteProgress` instances, including instances
    of its :meth:`~git.util.CallableRemoteProgress` subclass.

:note:
    Unlike objects that match this type, :meth:`~git.util.RemoteProgress` instances are
    not directly callable, not even when they are instances of
    :meth:`~git.util.CallableRemoteProgress`, which wraps a callable and forwards
    information to it but is not itself callable.

:note:
    This type also allows ``None``, for cloning without reporting progress.
"""

# -----------------------------------------------------------------------------------


def assert_never(inp: NoReturn, raise_error: bool = True, exc: Union[Exception, None] = None) -> None:
    """For use in exhaustive checking of a literal or enum in if/else chains.

    A call to this function should only be reached if not all members are handled, or if
    an attempt is made to pass non-members through the chain.

    :param inp:
        If all members are handled, the argument for `inp` will have the
        :class:`~typing.Never`/:class:`~typing.NoReturn` type.
        Otherwise, the type will mismatch and cause a mypy error.

    :param raise_error:
        If ``True``, will also raise :exc:`ValueError` with a general
        "unhandled literal" message, or the exception object passed as `exc`.

    :param exc:
        It not ``None``, this should be an already-constructed exception object, to be
        raised if `raise_error` is ``True``.
    """
    if raise_error:
        if exc is None:
            raise ValueError(f"An unhandled literal ({inp!r}) in an if/else chain was found")
        else:
            raise exc


class Files_TD(TypedDict):
    """Dictionary with stat counts for the diff of a particular file.

    For the :class:`~git.util.Stats.files` attribute of :class:`~git.util.Stats`
    objects.
    """

    insertions: int
    deletions: int
    lines: int
    change_type: str


class Total_TD(TypedDict):
    """Dictionary with total stats from any number of files.

    For the :class:`~git.util.Stats.total` attribute of :class:`~git.util.Stats`
    objects.
    """

    insertions: int
    deletions: int
    lines: int
    files: int


class HSH_TD(TypedDict):
    """Dictionary carrying the same information as a :class:`~git.util.Stats` object."""

    total: Total_TD
    files: Dict[PathLike, Files_TD]


@runtime_checkable
class Has_Repo(Protocol):
    """Protocol for having a :attr:`repo` attribute, the repository to operate on."""

    repo: "Repo"


@runtime_checkable
class Has_id_attribute(Protocol):
    """Protocol for having :attr:`_id_attribute_` used in iteration and traversal."""

    _id_attribute_: str
