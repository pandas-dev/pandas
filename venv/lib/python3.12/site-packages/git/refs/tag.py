# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

"""Provides a :class:`~git.refs.reference.Reference`-based type for lightweight tags.

This defines the :class:`TagReference` class (and its alias :class:`Tag`), which
represents lightweight tags. For annotated tags (which are git objects), see the
:mod:`git.objects.tag` module.
"""

__all__ = ["TagReference", "Tag"]

from .reference import Reference

# typing ------------------------------------------------------------------

from typing import Any, TYPE_CHECKING, Type, Union

from git.types import AnyGitObject, PathLike

if TYPE_CHECKING:
    from git.objects import Commit, TagObject
    from git.refs import SymbolicReference
    from git.repo import Repo

# ------------------------------------------------------------------------------


class TagReference(Reference):
    """A lightweight tag reference which either points to a commit, a tag object or any
    other object. In the latter case additional information, like the signature or the
    tag-creator, is available.

    This tag object will always point to a commit object, but may carry additional
    information in a tag object::

     tagref = TagReference.list_items(repo)[0]
     print(tagref.commit.message)
     if tagref.tag is not None:
        print(tagref.tag.message)
    """

    __slots__ = ()

    _common_default = "tags"
    _common_path_default = Reference._common_path_default + "/" + _common_default

    @property  # type: ignore[misc]
    def commit(self) -> "Commit":  # LazyMixin has unrelated commit method
        """:return: Commit object the tag ref points to

        :raise ValueError:
            If the tag points to a tree or blob.
        """
        obj = self.object
        while obj.type != "commit":
            if obj.type == "tag":
                # It is a tag object which carries the commit as an object - we can point to anything.
                obj = obj.object
            else:
                raise ValueError(
                    (
                        "Cannot resolve commit as tag %s points to a %s object - "
                        + "use the `.object` property instead to access it"
                    )
                    % (self, obj.type)
                )
        return obj

    @property
    def tag(self) -> Union["TagObject", None]:
        """
        :return:
            Tag object this tag ref points to, or ``None`` in case we are a lightweight
            tag
        """
        obj = self.object
        if obj.type == "tag":
            return obj
        return None

    # Make object read-only. It should be reasonably hard to adjust an existing tag.
    @property  # type: ignore[misc]
    def object(self) -> AnyGitObject:
        return Reference._get_object(self)

    @classmethod
    def create(
        cls: Type["TagReference"],
        repo: "Repo",
        path: PathLike,
        reference: Union[str, "SymbolicReference"] = "HEAD",
        logmsg: Union[str, None] = None,
        force: bool = False,
        **kwargs: Any,
    ) -> "TagReference":
        """Create a new tag reference.

        :param repo:
            The :class:`~git.repo.base.Repo` to create the tag in.

        :param path:
            The name of the tag, e.g. ``1.0`` or ``releases/1.0``.
            The prefix ``refs/tags`` is implied.

        :param reference:
            A reference to the :class:`~git.objects.base.Object` you want to tag.
            The referenced object can be a commit, tree, or blob.

        :param logmsg:
            If not ``None``, the message will be used in your tag object. This will also
            create an additional tag object that allows to obtain that information,
            e.g.::

                tagref.tag.message

        :param message:
            Synonym for the `logmsg` parameter. Included for backwards compatibility.
            `logmsg` takes precedence if both are passed.

        :param force:
            If ``True``, force creation of a tag even though that tag already exists.

        :param kwargs:
            Additional keyword arguments to be passed to :manpage:`git-tag(1)`.

        :return:
            A new :class:`TagReference`.
        """
        if "ref" in kwargs and kwargs["ref"]:
            reference = kwargs["ref"]

        if "message" in kwargs and kwargs["message"]:
            kwargs["m"] = kwargs["message"]
            del kwargs["message"]

        if logmsg:
            kwargs["m"] = logmsg

        if force:
            kwargs["f"] = True

        args = (path, reference)

        repo.git.tag(*args, **kwargs)
        return TagReference(repo, "%s/%s" % (cls._common_path_default, path))

    @classmethod
    def delete(cls, repo: "Repo", *tags: "TagReference") -> None:  # type: ignore[override]
        """Delete the given existing tag or tags."""
        repo.git.tag("-d", *tags)


# Provide an alias.
Tag = TagReference
