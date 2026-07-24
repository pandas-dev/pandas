# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

__all__ = ["Object", "IndexObject"]

import os.path as osp

import gitdb.typ as dbtyp

from git.exc import WorkTreeRepositoryUnsupported
from git.util import LazyMixin, bin_to_hex, join_path_native, stream_copy

from .util import get_object_type_by_name

# typing ------------------------------------------------------------------

from typing import Any, TYPE_CHECKING, Union

from git.types import AnyGitObject, GitObjectTypeString, PathLike

if TYPE_CHECKING:
    from gitdb.base import OStream

    from git.refs.reference import Reference
    from git.repo import Repo

    from .blob import Blob
    from .submodule.base import Submodule
    from .tree import Tree

IndexObjUnion = Union["Tree", "Blob", "Submodule"]

# --------------------------------------------------------------------------


class Object(LazyMixin):
    """Base class for classes representing git object types.

    The following four leaf classes represent specific kinds of git objects:

    * :class:`Blob <git.objects.blob.Blob>`
    * :class:`Tree <git.objects.tree.Tree>`
    * :class:`Commit <git.objects.commit.Commit>`
    * :class:`TagObject <git.objects.tag.TagObject>`

    See :manpage:`gitglossary(7)` on:

    * "object": https://git-scm.com/docs/gitglossary#def_object
    * "object type": https://git-scm.com/docs/gitglossary#def_object_type
    * "blob": https://git-scm.com/docs/gitglossary#def_blob_object
    * "tree object": https://git-scm.com/docs/gitglossary#def_tree_object
    * "commit object": https://git-scm.com/docs/gitglossary#def_commit_object
    * "tag object": https://git-scm.com/docs/gitglossary#def_tag_object

    :note:
        See the :class:`~git.types.AnyGitObject` union type of the four leaf subclasses
        that represent actual git object types.

    :note:
        :class:`~git.objects.submodule.base.Submodule` is defined under the hierarchy
        rooted at this :class:`Object` class, even though submodules are not really a
        type of git object. (This also applies to its
        :class:`~git.objects.submodule.root.RootModule` subclass.)

    :note:
        This :class:`Object` class should not be confused with :class:`object` (the root
        of the class hierarchy in Python).
    """

    NULL_HEX_SHA = "0" * 40
    NULL_BIN_SHA = b"\0" * 20

    TYPES = (
        dbtyp.str_blob_type,
        dbtyp.str_tree_type,
        dbtyp.str_commit_type,
        dbtyp.str_tag_type,
    )

    __slots__ = ("repo", "binsha", "size")

    type: Union[GitObjectTypeString, None] = None
    """String identifying (a concrete :class:`Object` subtype for) a git object type.

    The subtypes that this may name correspond to the kinds of git objects that exist,
    i.e., the objects that may be present in a git repository.

    :note:
        Most subclasses represent specific types of git objects and override this class
        attribute accordingly. This attribute is ``None`` in the :class:`Object` base
        class, as well as the :class:`IndexObject` intermediate subclass, but never
        ``None`` in concrete leaf subclasses representing specific git object types.

    :note:
        See also :class:`~git.types.GitObjectTypeString`.
    """

    def __init__(self, repo: "Repo", binsha: bytes) -> None:
        """Initialize an object by identifying it by its binary sha.

        All keyword arguments will be set on demand if ``None``.

        :param repo:
            Repository this object is located in.

        :param binsha:
            20 byte SHA1
        """
        super().__init__()
        self.repo = repo
        self.binsha = binsha
        assert len(binsha) == 20, "Require 20 byte binary sha, got %r, len = %i" % (
            binsha,
            len(binsha),
        )

    @classmethod
    def new(cls, repo: "Repo", id: Union[str, "Reference"]) -> AnyGitObject:
        """
        :return:
            New :class:`Object` instance of a type appropriate to the object type behind
            `id`. The id of the newly created object will be a binsha even though the
            input id may have been a :class:`~git.refs.reference.Reference` or rev-spec.

        :param id:
            :class:`~git.refs.reference.Reference`, rev-spec, or hexsha.

        :note:
            This cannot be a ``__new__`` method as it would always call :meth:`__init__`
            with the input id which is not necessarily a binsha.
        """
        return repo.rev_parse(str(id))

    @classmethod
    def new_from_sha(cls, repo: "Repo", sha1: bytes) -> AnyGitObject:
        """
        :return:
            New object instance of a type appropriate to represent the given binary sha1

        :param sha1:
            20 byte binary sha1.
        """
        if sha1 == cls.NULL_BIN_SHA:
            # The NULL binsha is always the root commit.
            return get_object_type_by_name(b"commit")(repo, sha1)
        # END handle special case
        oinfo = repo.odb.info(sha1)
        inst = get_object_type_by_name(oinfo.type)(repo, oinfo.binsha)
        inst.size = oinfo.size
        return inst

    def _set_cache_(self, attr: str) -> None:
        """Retrieve object information."""
        if attr == "size":
            oinfo = self.repo.odb.info(self.binsha)
            self.size = oinfo.size  # type: int
        else:
            super()._set_cache_(attr)

    def __eq__(self, other: Any) -> bool:
        """:return: ``True`` if the objects have the same SHA1"""
        if not hasattr(other, "binsha"):
            return False
        return self.binsha == other.binsha

    def __ne__(self, other: Any) -> bool:
        """:return: ``True`` if the objects do not have the same SHA1"""
        if not hasattr(other, "binsha"):
            return True
        return self.binsha != other.binsha

    def __hash__(self) -> int:
        """:return: Hash of our id allowing objects to be used in dicts and sets"""
        return hash(self.binsha)

    def __str__(self) -> str:
        """:return: String of our SHA1 as understood by all git commands"""
        return self.hexsha

    def __repr__(self) -> str:
        """:return: String with pythonic representation of our object"""
        return '<git.%s "%s">' % (self.__class__.__name__, self.hexsha)

    @property
    def hexsha(self) -> str:
        """:return: 40 byte hex version of our 20 byte binary sha"""
        # b2a_hex produces bytes.
        return bin_to_hex(self.binsha).decode("ascii")

    @property
    def data_stream(self) -> "OStream":
        """
        :return:
            File-object compatible stream to the uncompressed raw data of the object

        :note:
            Returned streams must be read in order.
        """
        return self.repo.odb.stream(self.binsha)

    def stream_data(self, ostream: "OStream") -> "Object":
        """Write our data directly to the given output stream.

        :param ostream:
            File-object compatible stream object.

        :return:
            self
        """
        istream = self.repo.odb.stream(self.binsha)
        stream_copy(istream, ostream)
        return self


class IndexObject(Object):
    """Base for all objects that can be part of the index file.

    The classes representing git object types that can be part of the index file are
    :class:`~git.objects.tree.Tree` and :class:`~git.objects.blob.Blob`. In addition,
    :class:`~git.objects.submodule.base.Submodule`, which is not really a git object
    type but can be part of an index file, is also a subclass.
    """

    __slots__ = ("path", "mode")

    # For compatibility with iterable lists.
    _id_attribute_ = "path"

    def __init__(
        self,
        repo: "Repo",
        binsha: bytes,
        mode: Union[None, int] = None,
        path: Union[None, PathLike] = None,
    ) -> None:
        """Initialize a newly instanced :class:`IndexObject`.

        :param repo:
            The :class:`~git.repo.base.Repo` we are located in.

        :param binsha:
            20 byte sha1.

        :param mode:
            The stat-compatible file mode as :class:`int`.
            Use the :mod:`stat` module to evaluate the information.

        :param path:
            The path to the file in the file system, relative to the git repository
            root, like ``file.ext`` or ``folder/other.ext``.

        :note:
            Path may not be set if the index object has been created directly, as it
            cannot be retrieved without knowing the parent tree.
        """
        super().__init__(repo, binsha)
        if mode is not None:
            self.mode = mode
        if path is not None:
            self.path = path

    def __hash__(self) -> int:
        """
        :return:
            Hash of our path as index items are uniquely identifiable by path, not by
            their data!
        """
        return hash(self.path)

    def _set_cache_(self, attr: str) -> None:
        if attr in IndexObject.__slots__:
            # They cannot be retrieved later on (not without searching for them).
            raise AttributeError(
                "Attribute '%s' unset: path and mode attributes must have been set during %s object creation"
                % (attr, type(self).__name__)
            )
        else:
            super()._set_cache_(attr)
        # END handle slot attribute

    @property
    def name(self) -> str:
        """:return: Name portion of the path, effectively being the basename"""
        return osp.basename(self.path)

    @property
    def abspath(self) -> PathLike:
        R"""
        :return:
            Absolute path to this index object in the file system (as opposed to the
            :attr:`path` field which is a path relative to the git repository).

            The returned path will be native to the system and contains ``\`` on
            Windows.
        """
        if self.repo.working_tree_dir is not None:
            return join_path_native(self.repo.working_tree_dir, self.path)
        else:
            raise WorkTreeRepositoryUnsupported("working_tree_dir was None or empty")
