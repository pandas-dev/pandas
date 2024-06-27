# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

from __future__ import annotations

__all__ = ["Repo"]

import gc
import logging
import os
import os.path as osp
from pathlib import Path
import re
import shlex
import sys
import warnings

import gitdb
from gitdb.db.loose import LooseObjectDB
from gitdb.exc import BadObject

from git.cmd import Git, handle_process_output
from git.compat import defenc, safe_decode
from git.config import GitConfigParser
from git.db import GitCmdObjectDB
from git.exc import (
    GitCommandError,
    InvalidGitRepositoryError,
    NoSuchPathError,
)
from git.index import IndexFile
from git.objects import Submodule, RootModule, Commit
from git.refs import HEAD, Head, Reference, TagReference
from git.remote import Remote, add_progress, to_progress_instance
from git.util import (
    Actor,
    cygpath,
    expand_path,
    finalize_process,
    hex_to_bin,
    remove_password_if_present,
)

from .fun import (
    find_submodule_git_dir,
    find_worktree_git_dir,
    is_git_dir,
    rev_parse,
    touch,
)

# typing ------------------------------------------------------

from git.types import (
    CallableProgress,
    Commit_ish,
    Lit_config_levels,
    PathLike,
    TBD,
    Tree_ish,
    assert_never,
)
from typing import (
    Any,
    BinaryIO,
    Callable,
    Dict,
    Iterator,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Sequence,
    TYPE_CHECKING,
    TextIO,
    Tuple,
    Type,
    Union,
    cast,
)

from git.types import ConfigLevels_Tup, TypedDict

if TYPE_CHECKING:
    from git.objects import Tree
    from git.objects.submodule.base import UpdateProgress
    from git.refs.symbolic import SymbolicReference
    from git.remote import RemoteProgress
    from git.util import IterableList

# -----------------------------------------------------------

_logger = logging.getLogger(__name__)


class BlameEntry(NamedTuple):
    commit: Dict[str, Commit]
    linenos: range
    orig_path: Optional[str]
    orig_linenos: range


class Repo:
    """Represents a git repository and allows you to query references, create commit
    information, generate diffs, create and clone repositories, and query the log.

    The following attributes are worth using:

    * :attr:`working_dir` is the working directory of the git command, which is the
      working tree directory if available or the ``.git`` directory in case of bare
      repositories.

    * :attr:`working_tree_dir` is the working tree directory, but will return ``None``
      if we are a bare repository.

    * :attr:`git_dir` is the ``.git`` repository directory, which is always set.
    """

    DAEMON_EXPORT_FILE = "git-daemon-export-ok"

    # Must exist, or  __del__  will fail in case we raise on `__init__()`.
    git = cast("Git", None)

    working_dir: PathLike
    """The working directory of the git command."""

    _working_tree_dir: Optional[PathLike] = None

    git_dir: PathLike
    """The ``.git`` repository directory."""

    _common_dir: PathLike = ""

    # Precompiled regex
    re_whitespace = re.compile(r"\s+")
    re_hexsha_only = re.compile(r"^[0-9A-Fa-f]{40}$")
    re_hexsha_shortened = re.compile(r"^[0-9A-Fa-f]{4,40}$")
    re_envvars = re.compile(r"(\$(\{\s?)?[a-zA-Z_]\w*(\}\s?)?|%\s?[a-zA-Z_]\w*\s?%)")
    re_author_committer_start = re.compile(r"^(author|committer)")
    re_tab_full_line = re.compile(r"^\t(.*)$")

    unsafe_git_clone_options = [
        # Executes arbitrary commands:
        "--upload-pack",
        "-u",
        # Can override configuration variables that execute arbitrary commands:
        "--config",
        "-c",
    ]
    """Options to :manpage:`git-clone(1)` that allow arbitrary commands to be executed.

    The ``--upload-pack``/``-u`` option allows users to execute arbitrary commands
    directly:
    https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---upload-packltupload-packgt

    The ``--config``/``-c`` option allows users to override configuration variables like
    ``protocol.allow`` and ``core.gitProxy`` to execute arbitrary commands:
    https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---configltkeygtltvaluegt
    """

    # Invariants
    config_level: ConfigLevels_Tup = ("system", "user", "global", "repository")
    """Represents the configuration level of a configuration file."""

    # Subclass configuration
    GitCommandWrapperType = Git
    """Subclasses may easily bring in their own custom types by placing a constructor or
    type here."""

    def __init__(
        self,
        path: Optional[PathLike] = None,
        odbt: Type[LooseObjectDB] = GitCmdObjectDB,
        search_parent_directories: bool = False,
        expand_vars: bool = True,
    ) -> None:
        R"""Create a new :class:`Repo` instance.

        :param path:
            The path to either the root git directory or the bare git repo::

                repo = Repo("/Users/mtrier/Development/git-python")
                repo = Repo("/Users/mtrier/Development/git-python.git")
                repo = Repo("~/Development/git-python.git")
                repo = Repo("$REPOSITORIES/Development/git-python.git")
                repo = Repo(R"C:\Users\mtrier\Development\git-python\.git")

            - In *Cygwin*, `path` may be a ``cygdrive/...`` prefixed path.
            - If `path` is ``None`` or an empty string, :envvar:`GIT_DIR` is used. If
              that environment variable is absent or empty, the current directory is
              used.

        :param odbt:
            Object DataBase type - a type which is constructed by providing the
            directory containing the database objects, i.e. ``.git/objects``. It will be
            used to access all object data.

        :param search_parent_directories:
            If ``True``, all parent directories will be searched for a valid repo as
            well.

            Please note that this was the default behaviour in older versions of
            GitPython, which is considered a bug though.

        :raise git.exc.InvalidGitRepositoryError:

        :raise git.exc.NoSuchPathError:

        :return:
            :class:`Repo`
        """

        epath = path or os.getenv("GIT_DIR")
        if not epath:
            epath = os.getcwd()
        if Git.is_cygwin():
            # Given how the tests are written, this seems more likely to catch Cygwin
            # git used from Windows than Windows git used from Cygwin. Therefore
            # changing to Cygwin-style paths is the relevant operation.
            epath = cygpath(str(epath))

        epath = epath or path or os.getcwd()
        if not isinstance(epath, str):
            epath = str(epath)
        if expand_vars and re.search(self.re_envvars, epath):
            warnings.warn(
                "The use of environment variables in paths is deprecated"
                + "\nfor security reasons and may be removed in the future!!",
                stacklevel=1,
            )
        epath = expand_path(epath, expand_vars)
        if epath is not None:
            if not os.path.exists(epath):
                raise NoSuchPathError(epath)

        # Walk up the path to find the `.git` dir.
        curpath = epath
        git_dir = None
        while curpath:
            # ABOUT osp.NORMPATH
            # It's important to normalize the paths, as submodules will otherwise
            # initialize their repo instances with paths that depend on path-portions
            # that will not exist after being removed. It's just cleaner.
            if is_git_dir(curpath):
                git_dir = curpath
                # from man git-config : core.worktree
                # Set the path to the root of the working tree. If GIT_COMMON_DIR
                # environment variable is set, core.worktree is ignored and not used for
                # determining the root of working tree. This can be overridden by the
                # GIT_WORK_TREE environment variable. The value can be an absolute path
                # or relative to the path to the .git directory, which is either
                # specified by GIT_DIR, or automatically discovered. If GIT_DIR is
                # specified but none of GIT_WORK_TREE and core.worktree is specified,
                # the current working directory is regarded as the top level of your
                # working tree.
                self._working_tree_dir = os.path.dirname(git_dir)
                if os.environ.get("GIT_COMMON_DIR") is None:
                    gitconf = self._config_reader("repository", git_dir)
                    if gitconf.has_option("core", "worktree"):
                        self._working_tree_dir = gitconf.get("core", "worktree")
                if "GIT_WORK_TREE" in os.environ:
                    self._working_tree_dir = os.getenv("GIT_WORK_TREE")
                break

            dotgit = osp.join(curpath, ".git")
            sm_gitpath = find_submodule_git_dir(dotgit)
            if sm_gitpath is not None:
                git_dir = osp.normpath(sm_gitpath)

            sm_gitpath = find_submodule_git_dir(dotgit)
            if sm_gitpath is None:
                sm_gitpath = find_worktree_git_dir(dotgit)

            if sm_gitpath is not None:
                git_dir = expand_path(sm_gitpath, expand_vars)
                self._working_tree_dir = curpath
                break

            if not search_parent_directories:
                break
            curpath, tail = osp.split(curpath)
            if not tail:
                break
        # END while curpath

        if git_dir is None:
            raise InvalidGitRepositoryError(epath)
        self.git_dir = git_dir

        self._bare = False
        try:
            self._bare = self.config_reader("repository").getboolean("core", "bare")
        except Exception:
            # Let's not assume the option exists, although it should.
            pass

        try:
            common_dir = (Path(self.git_dir) / "commondir").read_text().splitlines()[0].strip()
            self._common_dir = osp.join(self.git_dir, common_dir)
        except OSError:
            self._common_dir = ""

        # Adjust the working directory in case we are actually bare - we didn't know
        # that in the first place.
        if self._bare:
            self._working_tree_dir = None
        # END working dir handling

        self.working_dir: PathLike = self._working_tree_dir or self.common_dir
        self.git = self.GitCommandWrapperType(self.working_dir)

        # Special handling, in special times.
        rootpath = osp.join(self.common_dir, "objects")
        if issubclass(odbt, GitCmdObjectDB):
            self.odb = odbt(rootpath, self.git)
        else:
            self.odb = odbt(rootpath)

    def __enter__(self) -> "Repo":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def close(self) -> None:
        if self.git:
            self.git.clear_cache()
            # Tempfiles objects on Windows are holding references to open files until
            # they are collected by the garbage collector, thus preventing deletion.
            # TODO: Find these references and ensure they are closed and deleted
            # synchronously rather than forcing a gc collection.
            if sys.platform == "win32":
                gc.collect()
            gitdb.util.mman.collect()
            if sys.platform == "win32":
                gc.collect()

    def __eq__(self, rhs: object) -> bool:
        if isinstance(rhs, Repo):
            return self.git_dir == rhs.git_dir
        return False

    def __ne__(self, rhs: object) -> bool:
        return not self.__eq__(rhs)

    def __hash__(self) -> int:
        return hash(self.git_dir)

    # Description property
    def _get_description(self) -> str:
        filename = osp.join(self.git_dir, "description")
        with open(filename, "rb") as fp:
            return fp.read().rstrip().decode(defenc)

    def _set_description(self, descr: str) -> None:
        filename = osp.join(self.git_dir, "description")
        with open(filename, "wb") as fp:
            fp.write((descr + "\n").encode(defenc))

    description = property(_get_description, _set_description, doc="the project's description")
    del _get_description
    del _set_description

    @property
    def working_tree_dir(self) -> Optional[PathLike]:
        """
        :return:
            The working tree directory of our git repository.
            If this is a bare repository, ``None`` is returned.
        """
        return self._working_tree_dir

    @property
    def common_dir(self) -> PathLike:
        """
        :return:
            The git dir that holds everything except possibly HEAD, FETCH_HEAD,
            ORIG_HEAD, COMMIT_EDITMSG, index, and logs/.
        """
        return self._common_dir or self.git_dir

    @property
    def bare(self) -> bool:
        """:return: ``True`` if the repository is bare"""
        return self._bare

    @property
    def heads(self) -> "IterableList[Head]":
        """A list of :class:`~git.refs.head.Head` objects representing the branch heads
        in this repo.

        :return:
            ``git.IterableList(Head, ...)``
        """
        return Head.list_items(self)

    @property
    def references(self) -> "IterableList[Reference]":
        """A list of :class:`~git.refs.reference.Reference` objects representing tags,
        heads and remote references.

        :return:
            ``git.IterableList(Reference, ...)``
        """
        return Reference.list_items(self)

    # Alias for references.
    refs = references

    # Alias for heads.
    branches = heads

    @property
    def index(self) -> "IndexFile":
        """
        :return:
            A :class:`~git.index.base.IndexFile` representing this repository's index.

        :note:
            This property can be expensive, as the returned
            :class:`~git.index.base.IndexFile` will be reinitialized.
            It is recommended to reuse the object.
        """
        return IndexFile(self)

    @property
    def head(self) -> "HEAD":
        """
        :return:
            :class:`~git.refs.head.HEAD` object pointing to the current head reference
        """
        return HEAD(self, "HEAD")

    @property
    def remotes(self) -> "IterableList[Remote]":
        """A list of :class:`~git.remote.Remote` objects allowing to access and
        manipulate remotes.

        :return:
            ``git.IterableList(Remote, ...)``
        """
        return Remote.list_items(self)

    def remote(self, name: str = "origin") -> "Remote":
        """:return: The remote with the specified name

        :raise ValueError:
            If no remote with such a name exists.
        """
        r = Remote(self, name)
        if not r.exists():
            raise ValueError("Remote named '%s' didn't exist" % name)
        return r

    # { Submodules

    @property
    def submodules(self) -> "IterableList[Submodule]":
        """
        :return:
            git.IterableList(Submodule, ...) of direct submodules available from the
            current head
        """
        return Submodule.list_items(self)

    def submodule(self, name: str) -> "Submodule":
        """:return: The submodule with the given name

        :raise ValueError:
            If no such submodule exists.
        """
        try:
            return self.submodules[name]
        except IndexError as e:
            raise ValueError("Didn't find submodule named %r" % name) from e
        # END exception handling

    def create_submodule(self, *args: Any, **kwargs: Any) -> Submodule:
        """Create a new submodule.

        :note:
            For a description of the applicable parameters, see the documentation of
            :meth:`Submodule.add <git.objects.submodule.base.Submodule.add>`.

        :return:
            The created submodule.
        """
        return Submodule.add(self, *args, **kwargs)

    def iter_submodules(self, *args: Any, **kwargs: Any) -> Iterator[Submodule]:
        """An iterator yielding Submodule instances.

        See the `~git.objects.util.Traversable` interface for a description of `args`
        and `kwargs`.

        :return:
            Iterator
        """
        return RootModule(self).traverse(*args, **kwargs)

    def submodule_update(self, *args: Any, **kwargs: Any) -> Iterator[Submodule]:
        """Update the submodules, keeping the repository consistent as it will
        take the previous state into consideration.

        :note:
            For more information, please see the documentation of
            :meth:`RootModule.update <git.objects.submodule.root.RootModule.update>`.
        """
        return RootModule(self).update(*args, **kwargs)

    # }END submodules

    @property
    def tags(self) -> "IterableList[TagReference]":
        """A list of :class:`~git.refs.tag.TagReference` objects that are available in
        this repo.

        :return:
            ``git.IterableList(TagReference, ...)``
        """
        return TagReference.list_items(self)

    def tag(self, path: PathLike) -> TagReference:
        """
        :return:
            :class:`~git.refs.tag.TagReference` object, reference pointing to a
            :class:`~git.objects.commit.Commit` or tag

        :param path:
            Path to the tag reference, e.g. ``0.1.5`` or ``tags/0.1.5``.
        """
        full_path = self._to_full_tag_path(path)
        return TagReference(self, full_path)

    @staticmethod
    def _to_full_tag_path(path: PathLike) -> str:
        path_str = str(path)
        if path_str.startswith(TagReference._common_path_default + "/"):
            return path_str
        if path_str.startswith(TagReference._common_default + "/"):
            return Reference._common_path_default + "/" + path_str
        else:
            return TagReference._common_path_default + "/" + path_str

    def create_head(
        self,
        path: PathLike,
        commit: Union["SymbolicReference", "str"] = "HEAD",
        force: bool = False,
        logmsg: Optional[str] = None,
    ) -> "Head":
        """Create a new head within the repository.

        :note:
            For more documentation, please see the
            :meth:`Head.create <git.refs.head.Head.create>` method.

        :return:
            Newly created :class:`~git.refs.head.Head` Reference.
        """
        return Head.create(self, path, commit, logmsg, force)

    def delete_head(self, *heads: "Union[str, Head]", **kwargs: Any) -> None:
        """Delete the given heads.

        :param kwargs:
            Additional keyword arguments to be passed to :manpage:`git-branch(1)`.
        """
        return Head.delete(self, *heads, **kwargs)

    def create_tag(
        self,
        path: PathLike,
        ref: Union[str, "SymbolicReference"] = "HEAD",
        message: Optional[str] = None,
        force: bool = False,
        **kwargs: Any,
    ) -> TagReference:
        """Create a new tag reference.

        :note:
            For more documentation, please see the
            :meth:`TagReference.create <git.refs.tag.TagReference.create>` method.

        :return:
            :class:`~git.refs.tag.TagReference` object
        """
        return TagReference.create(self, path, ref, message, force, **kwargs)

    def delete_tag(self, *tags: TagReference) -> None:
        """Delete the given tag references."""
        return TagReference.delete(self, *tags)

    def create_remote(self, name: str, url: str, **kwargs: Any) -> Remote:
        """Create a new remote.

        For more information, please see the documentation of the
        :meth:`Remote.create <git.remote.Remote.create>` method.

        :return:
            :class:`~git.remote.Remote` reference
        """
        return Remote.create(self, name, url, **kwargs)

    def delete_remote(self, remote: "Remote") -> str:
        """Delete the given remote."""
        return Remote.remove(self, remote)

    def _get_config_path(self, config_level: Lit_config_levels, git_dir: Optional[PathLike] = None) -> str:
        if git_dir is None:
            git_dir = self.git_dir
        # We do not support an absolute path of the gitconfig on Windows.
        # Use the global config instead.
        if sys.platform == "win32" and config_level == "system":
            config_level = "global"

        if config_level == "system":
            return "/etc/gitconfig"
        elif config_level == "user":
            config_home = os.environ.get("XDG_CONFIG_HOME") or osp.join(os.environ.get("HOME", "~"), ".config")
            return osp.normpath(osp.expanduser(osp.join(config_home, "git", "config")))
        elif config_level == "global":
            return osp.normpath(osp.expanduser("~/.gitconfig"))
        elif config_level == "repository":
            repo_dir = self._common_dir or git_dir
            if not repo_dir:
                raise NotADirectoryError
            else:
                return osp.normpath(osp.join(repo_dir, "config"))
        else:
            assert_never(  # type: ignore[unreachable]
                config_level,
                ValueError(f"Invalid configuration level: {config_level!r}"),
            )

    def config_reader(
        self,
        config_level: Optional[Lit_config_levels] = None,
    ) -> GitConfigParser:
        """
        :return:
            :class:`~git.config.GitConfigParser` allowing to read the full git
            configuration, but not to write it.

            The configuration will include values from the system, user and repository
            configuration files.

        :param config_level:
            For possible values, see the :meth:`config_writer` method. If ``None``, all
            applicable levels will be used. Specify a level in case you know which file
            you wish to read to prevent reading multiple files.

        :note:
            On Windows, system configuration cannot currently be read as the path is
            unknown, instead the global path will be used.
        """
        return self._config_reader(config_level=config_level)

    def _config_reader(
        self,
        config_level: Optional[Lit_config_levels] = None,
        git_dir: Optional[PathLike] = None,
    ) -> GitConfigParser:
        if config_level is None:
            files = [
                self._get_config_path(cast(Lit_config_levels, f), git_dir)
                for f in self.config_level
                if cast(Lit_config_levels, f)
            ]
        else:
            files = [self._get_config_path(config_level, git_dir)]
        return GitConfigParser(files, read_only=True, repo=self)

    def config_writer(self, config_level: Lit_config_levels = "repository") -> GitConfigParser:
        """
        :return:
            A :class:`~git.config.GitConfigParser` allowing to write values of the
            specified configuration file level. Config writers should be retrieved, used
            to change the configuration, and written right away as they will lock the
            configuration file in question and prevent other's to write it.

        :param config_level:
            One of the following values:

            * ``"system"`` = system wide configuration file
            * ``"global"`` = user level configuration file
            * ``"`repository"`` = configuration file for this repository only
        """
        return GitConfigParser(self._get_config_path(config_level), read_only=False, repo=self, merge_includes=False)

    def commit(self, rev: Union[str, Commit_ish, None] = None) -> Commit:
        """The :class:`~git.objects.commit.Commit` object for the specified revision.

        :param rev:
            Revision specifier, see :manpage:`git-rev-parse(1)` for viable options.

        :return:
            :class:`~git.objects.commit.Commit`
        """
        if rev is None:
            return self.head.commit
        return self.rev_parse(str(rev) + "^0")

    def iter_trees(self, *args: Any, **kwargs: Any) -> Iterator["Tree"]:
        """:return: Iterator yielding :class:`~git.objects.tree.Tree` objects

        :note:
            Accepts all arguments known to the :meth:`iter_commits` method.
        """
        return (c.tree for c in self.iter_commits(*args, **kwargs))

    def tree(self, rev: Union[Tree_ish, str, None] = None) -> "Tree":
        """The :class:`~git.objects.tree.Tree` object for the given tree-ish revision.

        Examples::

              repo.tree(repo.heads[0])

        :param rev:
            A revision pointing to a Treeish (being a commit or tree).

        :return:
            :class:`~git.objects.tree.Tree`

        :note:
            If you need a non-root level tree, find it by iterating the root tree.
            Otherwise it cannot know about its path relative to the repository root and
            subsequent operations might have unexpected results.
        """
        if rev is None:
            return self.head.commit.tree
        return self.rev_parse(str(rev) + "^{tree}")

    def iter_commits(
        self,
        rev: Union[str, Commit, "SymbolicReference", None] = None,
        paths: Union[PathLike, Sequence[PathLike]] = "",
        **kwargs: Any,
    ) -> Iterator[Commit]:
        """An iterator of :class:`~git.objects.commit.Commit` objects representing the
        history of a given ref/commit.

        :param rev:
            Revision specifier, see :manpage:`git-rev-parse(1)` for viable options.
            If ``None``, the active branch will be used.

        :param paths:
            An optional path or a list of paths. If set, only commits that include the
            path or paths will be returned.

        :param kwargs:
            Arguments to be passed to :manpage:`git-rev-list(1)`.
            Common ones are ``max_count`` and ``skip``.

        :note:
            To receive only commits between two named revisions, use the
            ``"revA...revB"`` revision specifier.

        :return:
            Iterator of :class:`~git.objects.commit.Commit` objects
        """
        if rev is None:
            rev = self.head.commit

        return Commit.iter_items(self, rev, paths, **kwargs)

    def merge_base(self, *rev: TBD, **kwargs: Any) -> List[Commit]:
        R"""Find the closest common ancestor for the given revision
        (:class:`~git.objects.commit.Commit`\s, :class:`~git.refs.tag.Tag`\s,
        :class:`~git.refs.reference.Reference`\s, etc.).

        :param rev:
            At least two revs to find the common ancestor for.

        :param kwargs:
            Additional arguments to be passed to the ``repo.git.merge_base()`` command
            which does all the work.

        :return:
            A list of :class:`~git.objects.commit.Commit` objects. If ``--all`` was
            not passed as a keyword argument, the list will have at max one
            :class:`~git.objects.commit.Commit`, or is empty if no common merge base
            exists.

        :raise ValueError:
            If fewer than two revisions are provided.
        """
        if len(rev) < 2:
            raise ValueError("Please specify at least two revs, got only %i" % len(rev))
        # END handle input

        res: List[Commit] = []
        try:
            lines: List[str] = self.git.merge_base(*rev, **kwargs).splitlines()
        except GitCommandError as err:
            if err.status == 128:
                raise
            # END handle invalid rev
            # Status code 1 is returned if there is no merge-base.
            # (See: https://github.com/git/git/blob/v2.44.0/builtin/merge-base.c#L19)
            return res
        # END exception handling

        for line in lines:
            res.append(self.commit(line))
        # END for each merge-base

        return res

    def is_ancestor(self, ancestor_rev: Commit, rev: Commit) -> bool:
        """Check if a commit is an ancestor of another.

        :param ancestor_rev:
            Rev which should be an ancestor.

        :param rev:
            Rev to test against `ancestor_rev`.

        :return:
            ``True`` if `ancestor_rev` is an ancestor to `rev`.
        """
        try:
            self.git.merge_base(ancestor_rev, rev, is_ancestor=True)
        except GitCommandError as err:
            if err.status == 1:
                return False
            raise
        return True

    def is_valid_object(self, sha: str, object_type: Union[str, None] = None) -> bool:
        try:
            complete_sha = self.odb.partial_to_complete_sha_hex(sha)
            object_info = self.odb.info(complete_sha)
            if object_type:
                if object_info.type == object_type.encode():
                    return True
                else:
                    _logger.debug(
                        "Commit hash points to an object of type '%s'. Requested were objects of type '%s'",
                        object_info.type.decode(),
                        object_type,
                    )
                    return False
            else:
                return True
        except BadObject:
            _logger.debug("Commit hash is invalid.")
            return False

    def _get_daemon_export(self) -> bool:
        if self.git_dir:
            filename = osp.join(self.git_dir, self.DAEMON_EXPORT_FILE)
        return osp.exists(filename)

    def _set_daemon_export(self, value: object) -> None:
        if self.git_dir:
            filename = osp.join(self.git_dir, self.DAEMON_EXPORT_FILE)
        fileexists = osp.exists(filename)
        if value and not fileexists:
            touch(filename)
        elif not value and fileexists:
            os.unlink(filename)

    daemon_export = property(
        _get_daemon_export,
        _set_daemon_export,
        doc="If True, git-daemon may export this repository",
    )
    del _get_daemon_export
    del _set_daemon_export

    def _get_alternates(self) -> List[str]:
        """The list of alternates for this repo from which objects can be retrieved.

        :return:
            List of strings being pathnames of alternates
        """
        if self.git_dir:
            alternates_path = osp.join(self.git_dir, "objects", "info", "alternates")

        if osp.exists(alternates_path):
            with open(alternates_path, "rb") as f:
                alts = f.read().decode(defenc)
            return alts.strip().splitlines()
        return []

    def _set_alternates(self, alts: List[str]) -> None:
        """Set the alternates.

        :param alts:
            The array of string paths representing the alternates at which git should
            look for objects, i.e. ``/home/user/repo/.git/objects``.

        :raise git.exc.NoSuchPathError:

        :note:
            The method does not check for the existence of the paths in `alts`, as the
            caller is responsible.
        """
        alternates_path = osp.join(self.common_dir, "objects", "info", "alternates")
        if not alts:
            if osp.isfile(alternates_path):
                os.remove(alternates_path)
        else:
            with open(alternates_path, "wb") as f:
                f.write("\n".join(alts).encode(defenc))

    alternates = property(
        _get_alternates,
        _set_alternates,
        doc="Retrieve a list of alternates paths or set a list paths to be used as alternates",
    )

    def is_dirty(
        self,
        index: bool = True,
        working_tree: bool = True,
        untracked_files: bool = False,
        submodules: bool = True,
        path: Optional[PathLike] = None,
    ) -> bool:
        """
        :return:
            ``True`` if the repository is considered dirty. By default it will react
            like a :manpage:`git-status(1)` without untracked files, hence it is dirty
            if the index or the working copy have changes.
        """
        if self._bare:
            # Bare repositories with no associated working directory are
            # always considered to be clean.
            return False

        # Start from the one which is fastest to evaluate.
        default_args = ["--abbrev=40", "--full-index", "--raw"]
        if not submodules:
            default_args.append("--ignore-submodules")
        if path:
            default_args.extend(["--", str(path)])
        if index:
            # diff index against HEAD.
            if osp.isfile(self.index.path) and len(self.git.diff("--cached", *default_args)):
                return True
        # END index handling
        if working_tree:
            # diff index against working tree.
            if len(self.git.diff(*default_args)):
                return True
        # END working tree handling
        if untracked_files:
            if len(self._get_untracked_files(path, ignore_submodules=not submodules)):
                return True
        # END untracked files
        return False

    @property
    def untracked_files(self) -> List[str]:
        """
        :return:
            list(str,...)

            Files currently untracked as they have not been staged yet. Paths are
            relative to the current working directory of the git command.

        :note:
            Ignored files will not appear here, i.e. files mentioned in ``.gitignore``.

        :note:
            This property is expensive, as no cache is involved. To process the result,
            please consider caching it yourself.
        """
        return self._get_untracked_files()

    def _get_untracked_files(self, *args: Any, **kwargs: Any) -> List[str]:
        # Make sure we get all files, not only untracked directories.
        proc = self.git.status(*args, porcelain=True, untracked_files=True, as_process=True, **kwargs)
        # Untracked files prefix in porcelain mode
        prefix = "?? "
        untracked_files = []
        for line in proc.stdout:
            line = line.decode(defenc)
            if not line.startswith(prefix):
                continue
            filename = line[len(prefix) :].rstrip("\n")
            # Special characters are escaped
            if filename[0] == filename[-1] == '"':
                filename = filename[1:-1]
                # WHATEVER ... it's a mess, but works for me
                filename = filename.encode("ascii").decode("unicode_escape").encode("latin1").decode(defenc)
            untracked_files.append(filename)
        finalize_process(proc)
        return untracked_files

    def ignored(self, *paths: PathLike) -> List[str]:
        """Checks if paths are ignored via ``.gitignore``.

        This does so using the :manpage:`git-check-ignore(1)` method.

        :param paths:
            List of paths to check whether they are ignored or not.

        :return:
            Subset of those paths which are ignored
        """
        try:
            proc: str = self.git.check_ignore(*paths)
        except GitCommandError as err:
            if err.status == 1:
                # If return code is 1, this means none of the items in *paths are
                # ignored by Git, so return an empty list.
                return []
            else:
                # Raise the exception on all other return codes.
                raise

        return proc.replace("\\\\", "\\").replace('"', "").split("\n")

    @property
    def active_branch(self) -> Head:
        """The name of the currently active branch.

        :raise TypeError:
            If HEAD is detached.

        :return:
            :class:`~git.refs.head.Head` to the active branch
        """
        # reveal_type(self.head.reference)  # => Reference
        return self.head.reference

    def blame_incremental(self, rev: str | HEAD | None, file: str, **kwargs: Any) -> Iterator["BlameEntry"]:
        """Iterator for blame information for the given file at the given revision.

        Unlike :meth:`blame`, this does not return the actual file's contents, only a
        stream of :class:`BlameEntry` tuples.

        :param rev:
            Revision specifier. If ``None``, the blame will include all the latest
            uncommitted changes. Otherwise, anything successfully parsed by
            :manpage:`git-rev-parse(1)` is a valid option.

        :return:
            Lazy iterator of :class:`BlameEntry` tuples, where the commit indicates the
            commit to blame for the line, and range indicates a span of line numbers in
            the resulting file.

        If you combine all line number ranges outputted by this command, you should get
        a continuous range spanning all line numbers in the file.
        """

        data: bytes = self.git.blame(rev, "--", file, p=True, incremental=True, stdout_as_string=False, **kwargs)
        commits: Dict[bytes, Commit] = {}

        stream = (line for line in data.split(b"\n") if line)
        while True:
            try:
                # When exhausted, causes a StopIteration, terminating this function.
                line = next(stream)
            except StopIteration:
                return
            split_line = line.split()
            hexsha, orig_lineno_b, lineno_b, num_lines_b = split_line
            lineno = int(lineno_b)
            num_lines = int(num_lines_b)
            orig_lineno = int(orig_lineno_b)
            if hexsha not in commits:
                # Now read the next few lines and build up a dict of properties for this
                # commit.
                props: Dict[bytes, bytes] = {}
                while True:
                    try:
                        line = next(stream)
                    except StopIteration:
                        return
                    if line == b"boundary":
                        # "boundary" indicates a root commit and occurs instead of the
                        # "previous" tag.
                        continue

                    tag, value = line.split(b" ", 1)
                    props[tag] = value
                    if tag == b"filename":
                        # "filename" formally terminates the entry for --incremental.
                        orig_filename = value
                        break

                c = Commit(
                    self,
                    hex_to_bin(hexsha),
                    author=Actor(
                        safe_decode(props[b"author"]),
                        safe_decode(props[b"author-mail"].lstrip(b"<").rstrip(b">")),
                    ),
                    authored_date=int(props[b"author-time"]),
                    committer=Actor(
                        safe_decode(props[b"committer"]),
                        safe_decode(props[b"committer-mail"].lstrip(b"<").rstrip(b">")),
                    ),
                    committed_date=int(props[b"committer-time"]),
                )
                commits[hexsha] = c
            else:
                # Discard all lines until we find "filename" which is guaranteed to be
                # the last line.
                while True:
                    try:
                        # Will fail if we reach the EOF unexpectedly.
                        line = next(stream)
                    except StopIteration:
                        return
                    tag, value = line.split(b" ", 1)
                    if tag == b"filename":
                        orig_filename = value
                        break

            yield BlameEntry(
                commits[hexsha],
                range(lineno, lineno + num_lines),
                safe_decode(orig_filename),
                range(orig_lineno, orig_lineno + num_lines),
            )

    def blame(
        self,
        rev: Union[str, HEAD, None],
        file: str,
        incremental: bool = False,
        rev_opts: Optional[List[str]] = None,
        **kwargs: Any,
    ) -> List[List[Commit | List[str | bytes] | None]] | Iterator[BlameEntry] | None:
        """The blame information for the given file at the given revision.

        :param rev:
            Revision specifier. If ``None``, the blame will include all the latest
            uncommitted changes. Otherwise, anything successfully parsed by
            :manpage:`git-rev-parse(1)` is a valid option.

        :return:
            list: [git.Commit, list: [<line>]]

            A list of lists associating a :class:`~git.objects.commit.Commit` object
            with a list of lines that changed within the given commit. The
            :class:`~git.objects.commit.Commit` objects will be given in order of
            appearance.
        """
        if incremental:
            return self.blame_incremental(rev, file, **kwargs)
        rev_opts = rev_opts or []
        data: bytes = self.git.blame(rev, *rev_opts, "--", file, p=True, stdout_as_string=False, **kwargs)
        commits: Dict[str, Commit] = {}
        blames: List[List[Commit | List[str | bytes] | None]] = []

        class InfoTD(TypedDict, total=False):
            sha: str
            id: str
            filename: str
            summary: str
            author: str
            author_email: str
            author_date: int
            committer: str
            committer_email: str
            committer_date: int

        info: InfoTD = {}

        keepends = True
        for line_bytes in data.splitlines(keepends):
            try:
                line_str = line_bytes.rstrip().decode(defenc)
            except UnicodeDecodeError:
                firstpart = ""
                parts = []
                is_binary = True
            else:
                # As we don't have an idea when the binary data ends, as it could
                # contain multiple newlines in the process. So we rely on being able to
                # decode to tell us what it is. This can absolutely fail even on text
                # files, but even if it does, we should be fine treating it as binary
                # instead.
                parts = self.re_whitespace.split(line_str, 1)
                firstpart = parts[0]
                is_binary = False
            # END handle decode of line

            if self.re_hexsha_only.search(firstpart):
                # handles
                # 634396b2f541a9f2d58b00be1a07f0c358b999b3 1 1 7        - indicates blame-data start
                # 634396b2f541a9f2d58b00be1a07f0c358b999b3 2 2          - indicates
                # another line of blame with the same data
                digits = parts[-1].split(" ")
                if len(digits) == 3:
                    info = {"id": firstpart}
                    blames.append([None, []])
                elif info["id"] != firstpart:
                    info = {"id": firstpart}
                    blames.append([commits.get(firstpart), []])
                # END blame data initialization
            else:
                m = self.re_author_committer_start.search(firstpart)
                if m:
                    # handles:
                    # author Tom Preston-Werner
                    # author-mail <tom@mojombo.com>
                    # author-time 1192271832
                    # author-tz -0700
                    # committer Tom Preston-Werner
                    # committer-mail <tom@mojombo.com>
                    # committer-time 1192271832
                    # committer-tz -0700  - IGNORED BY US
                    role = m.group(0)
                    if role == "author":
                        if firstpart.endswith("-mail"):
                            info["author_email"] = parts[-1]
                        elif firstpart.endswith("-time"):
                            info["author_date"] = int(parts[-1])
                        elif role == firstpart:
                            info["author"] = parts[-1]
                    elif role == "committer":
                        if firstpart.endswith("-mail"):
                            info["committer_email"] = parts[-1]
                        elif firstpart.endswith("-time"):
                            info["committer_date"] = int(parts[-1])
                        elif role == firstpart:
                            info["committer"] = parts[-1]
                    # END distinguish mail,time,name
                else:
                    # handle
                    # filename lib/grit.rb
                    # summary add Blob
                    # <and rest>
                    if firstpart.startswith("filename"):
                        info["filename"] = parts[-1]
                    elif firstpart.startswith("summary"):
                        info["summary"] = parts[-1]
                    elif firstpart == "":
                        if info:
                            sha = info["id"]
                            c = commits.get(sha)
                            if c is None:
                                c = Commit(
                                    self,
                                    hex_to_bin(sha),
                                    author=Actor._from_string(f"{info['author']} {info['author_email']}"),
                                    authored_date=info["author_date"],
                                    committer=Actor._from_string(f"{info['committer']} {info['committer_email']}"),
                                    committed_date=info["committer_date"],
                                )
                                commits[sha] = c
                            blames[-1][0] = c
                            # END if commit objects needs initial creation

                            if blames[-1][1] is not None:
                                line: str | bytes
                                if not is_binary:
                                    if line_str and line_str[0] == "\t":
                                        line_str = line_str[1:]
                                    line = line_str
                                else:
                                    line = line_bytes
                                    # NOTE: We are actually parsing lines out of binary
                                    # data, which can lead to the binary being split up
                                    # along the newline separator. We will append this
                                    # to the blame we are currently looking at, even
                                    # though it should be concatenated with the last
                                    # line we have seen.
                                blames[-1][1].append(line)

                            info = {"id": sha}
                        # END if we collected commit info
                    # END distinguish filename,summary,rest
                # END distinguish author|committer vs filename,summary,rest
            # END distinguish hexsha vs other information
        return blames

    @classmethod
    def init(
        cls,
        path: Union[PathLike, None] = None,
        mkdir: bool = True,
        odbt: Type[GitCmdObjectDB] = GitCmdObjectDB,
        expand_vars: bool = True,
        **kwargs: Any,
    ) -> "Repo":
        """Initialize a git repository at the given path if specified.

        :param path:
            The full path to the repo (traditionally ends with ``/<name>.git``). Or
            ``None``, in which case the repository will be created in the current
            working directory.

        :param mkdir:
            If specified, will create the repository directory if it doesn't already
            exist. Creates the directory with a mode=0755.
            Only effective if a path is explicitly given.

        :param odbt:
            Object DataBase type - a type which is constructed by providing the
            directory containing the database objects, i.e. ``.git/objects``. It will be
            used to access all object data.

        :param expand_vars:
            If specified, environment variables will not be escaped. This can lead to
            information disclosure, allowing attackers to access the contents of
            environment variables.

        :param kwargs:
            Keyword arguments serving as additional options to the
            :manpage:`git-init(1)` command.

        :return:
            :class:`Repo` (the newly created repo)
        """
        if path:
            path = expand_path(path, expand_vars)
        if mkdir and path and not osp.exists(path):
            os.makedirs(path, 0o755)

        # git command automatically chdir into the directory
        git = cls.GitCommandWrapperType(path)
        git.init(**kwargs)
        return cls(path, odbt=odbt)

    @classmethod
    def _clone(
        cls,
        git: "Git",
        url: PathLike,
        path: PathLike,
        odb_default_type: Type[GitCmdObjectDB],
        progress: Union["RemoteProgress", "UpdateProgress", Callable[..., "RemoteProgress"], None] = None,
        multi_options: Optional[List[str]] = None,
        allow_unsafe_protocols: bool = False,
        allow_unsafe_options: bool = False,
        **kwargs: Any,
    ) -> "Repo":
        odbt = kwargs.pop("odbt", odb_default_type)

        # When pathlib.Path or other class-based path is passed
        if not isinstance(path, str):
            path = str(path)

        ## A bug win cygwin's Git, when `--bare` or `--separate-git-dir`
        #  it prepends the cwd or(?) the `url` into the `path, so::
        #        git clone --bare  /cygwin/d/foo.git  C:\\Work
        #  becomes::
        #        git clone --bare  /cygwin/d/foo.git  /cygwin/d/C:\\Work
        #
        clone_path = Git.polish_url(path) if Git.is_cygwin() and "bare" in kwargs else path
        sep_dir = kwargs.get("separate_git_dir")
        if sep_dir:
            kwargs["separate_git_dir"] = Git.polish_url(sep_dir)
        multi = None
        if multi_options:
            multi = shlex.split(" ".join(multi_options))

        if not allow_unsafe_protocols:
            Git.check_unsafe_protocols(str(url))
        if not allow_unsafe_options:
            Git.check_unsafe_options(options=list(kwargs.keys()), unsafe_options=cls.unsafe_git_clone_options)
        if not allow_unsafe_options and multi_options:
            Git.check_unsafe_options(options=multi_options, unsafe_options=cls.unsafe_git_clone_options)

        proc = git.clone(
            multi,
            "--",
            Git.polish_url(str(url)),
            clone_path,
            with_extended_output=True,
            as_process=True,
            v=True,
            universal_newlines=True,
            **add_progress(kwargs, git, progress),
        )
        if progress:
            handle_process_output(
                proc,
                None,
                to_progress_instance(progress).new_message_handler(),
                finalize_process,
                decode_streams=False,
            )
        else:
            (stdout, stderr) = proc.communicate()
            cmdline = getattr(proc, "args", "")
            cmdline = remove_password_if_present(cmdline)

            _logger.debug("Cmd(%s)'s unused stdout: %s", cmdline, stdout)
            finalize_process(proc, stderr=stderr)

        # Our git command could have a different working dir than our actual
        # environment, hence we prepend its working dir if required.
        if not osp.isabs(path):
            path = osp.join(git._working_dir, path) if git._working_dir is not None else path

        repo = cls(path, odbt=odbt)

        # Retain env values that were passed to _clone().
        repo.git.update_environment(**git.environment())

        # Adjust remotes - there may be operating systems which use backslashes, These
        # might be given as initial paths, but when handling the config file that
        # contains the remote from which we were clones, git stops liking it as it will
        # escape the backslashes. Hence we undo the escaping just to be sure.
        if repo.remotes:
            with repo.remotes[0].config_writer as writer:
                writer.set_value("url", Git.polish_url(repo.remotes[0].url))
        # END handle remote repo
        return repo

    def clone(
        self,
        path: PathLike,
        progress: Optional[CallableProgress] = None,
        multi_options: Optional[List[str]] = None,
        allow_unsafe_protocols: bool = False,
        allow_unsafe_options: bool = False,
        **kwargs: Any,
    ) -> "Repo":
        """Create a clone from this repository.

        :param path:
            The full path of the new repo (traditionally ends with ``./<name>.git``).

        :param progress:
            See :meth:`Remote.push <git.remote.Remote.push>`.

        :param multi_options:
            A list of :manpage:`git-clone(1)` options that can be provided multiple
            times.

            One option per list item which is passed exactly as specified to clone.
            For example::

                [
                    "--config core.filemode=false",
                    "--config core.ignorecase",
                    "--recurse-submodule=repo1_path",
                    "--recurse-submodule=repo2_path",
                ]

        :param allow_unsafe_protocols:
            Allow unsafe protocols to be used, like ``ext``.

        :param allow_unsafe_options:
            Allow unsafe options to be used, like ``--upload-pack``.

        :param kwargs:
            * ``odbt`` = ObjectDatabase Type, allowing to determine the object database
              implementation used by the returned :class:`Repo` instance.
            * All remaining keyword arguments are given to the :manpage:`git-clone(1)`
              command.

        :return:
            :class:`Repo` (the newly cloned repo)
        """
        return self._clone(
            self.git,
            self.common_dir,
            path,
            type(self.odb),
            progress,
            multi_options,
            allow_unsafe_protocols=allow_unsafe_protocols,
            allow_unsafe_options=allow_unsafe_options,
            **kwargs,
        )

    @classmethod
    def clone_from(
        cls,
        url: PathLike,
        to_path: PathLike,
        progress: CallableProgress = None,
        env: Optional[Mapping[str, str]] = None,
        multi_options: Optional[List[str]] = None,
        allow_unsafe_protocols: bool = False,
        allow_unsafe_options: bool = False,
        **kwargs: Any,
    ) -> "Repo":
        """Create a clone from the given URL.

        :param url:
            Valid git url, see: https://git-scm.com/docs/git-clone#URLS

        :param to_path:
            Path to which the repository should be cloned to.

        :param progress:
            See :meth:`Remote.push <git.remote.Remote.push>`.

        :param env:
            Optional dictionary containing the desired environment variables.

            Note: Provided variables will be used to update the execution environment
            for ``git``. If some variable is not specified in `env` and is defined in
            :attr:`os.environ`, value from :attr:`os.environ` will be used. If you want
            to unset some variable, consider providing empty string as its value.

        :param multi_options:
            See the :meth:`clone` method.

        :param allow_unsafe_protocols:
            Allow unsafe protocols to be used, like ``ext``.

        :param allow_unsafe_options:
            Allow unsafe options to be used, like ``--upload-pack``.

        :param kwargs:
            See the :meth:`clone` method.

        :return:
            :class:`Repo` instance pointing to the cloned directory.
        """
        git = cls.GitCommandWrapperType(os.getcwd())
        if env is not None:
            git.update_environment(**env)
        return cls._clone(
            git,
            url,
            to_path,
            GitCmdObjectDB,
            progress,
            multi_options,
            allow_unsafe_protocols=allow_unsafe_protocols,
            allow_unsafe_options=allow_unsafe_options,
            **kwargs,
        )

    def archive(
        self,
        ostream: Union[TextIO, BinaryIO],
        treeish: Optional[str] = None,
        prefix: Optional[str] = None,
        **kwargs: Any,
    ) -> Repo:
        """Archive the tree at the given revision.

        :param ostream:
            File-compatible stream object to which the archive will be written as bytes.

        :param treeish:
            The treeish name/id, defaults to active branch.

        :param prefix:
            The optional prefix to prepend to each filename in the archive.

        :param kwargs:
            Additional arguments passed to :manpage:`git-archive(1)`:

            * Use the ``format`` argument to define the kind of format. Use specialized
              ostreams to write any format supported by Python.
            * You may specify the special ``path`` keyword, which may either be a
              repository-relative path to a directory or file to place into the archive,
              or a list or tuple of multiple paths.

        :raise git.exc.GitCommandError:
            If something went wrong.

        :return:
            self
        """
        if treeish is None:
            treeish = self.head.commit
        if prefix and "prefix" not in kwargs:
            kwargs["prefix"] = prefix
        kwargs["output_stream"] = ostream
        path = kwargs.pop("path", [])
        path = cast(Union[PathLike, List[PathLike], Tuple[PathLike, ...]], path)
        if not isinstance(path, (tuple, list)):
            path = [path]
        # END ensure paths is list (or tuple)
        self.git.archive("--", treeish, *path, **kwargs)
        return self

    def has_separate_working_tree(self) -> bool:
        """
        :return:
            True if our :attr:`git_dir` is not at the root of our
            :attr:`working_tree_dir`, but a ``.git`` file with a platform-agnostic
            symbolic link. Our :attr:`git_dir` will be wherever the ``.git`` file points
            to.

        :note:
            Bare repositories will always return ``False`` here.
        """
        if self.bare:
            return False
        if self.working_tree_dir:
            return osp.isfile(osp.join(self.working_tree_dir, ".git"))
        else:
            return False  # Or raise Error?

    rev_parse = rev_parse

    def __repr__(self) -> str:
        clazz = self.__class__
        return "<%s.%s %r>" % (clazz.__module__, clazz.__name__, self.git_dir)

    def currently_rebasing_on(self) -> Commit | None:
        """
        :return:
            The commit which is currently being replayed while rebasing.

            ``None`` if we are not currently rebasing.
        """
        if self.git_dir:
            rebase_head_file = osp.join(self.git_dir, "REBASE_HEAD")
        if not osp.isfile(rebase_head_file):
            return None
        with open(rebase_head_file, "rt") as f:
            content = f.readline().strip()
        return self.commit(content)
