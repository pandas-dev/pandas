# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

__all__ = ["Submodule", "UpdateProgress"]

import gc
from io import BytesIO
import logging
import os
import os.path as osp
import stat
import sys
import uuid
import urllib

import git
from git.cmd import Git
from git.compat import defenc
from git.config import GitConfigParser, SectionConstraint, cp
from git.exc import (
    BadName,
    InvalidGitRepositoryError,
    NoSuchPathError,
    RepositoryDirtyError,
)
from git.objects.base import IndexObject, Object
from git.objects.util import TraversableIterableObj
from git.util import (
    IterableList,
    RemoteProgress,
    join_path_native,
    rmtree,
    to_native_path_linux,
    unbare_repo,
)

from .util import (
    SubmoduleConfigParser,
    find_first_remote_branch,
    mkhead,
    sm_name,
    sm_section,
)

# typing ----------------------------------------------------------------------

from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    Mapping,
    Sequence,
    TYPE_CHECKING,
    Union,
    cast,
)

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

from git.types import Commit_ish, PathLike, TBD

if TYPE_CHECKING:
    from git.index import IndexFile
    from git.objects.commit import Commit
    from git.refs import Head, RemoteReference
    from git.repo import Repo

# -----------------------------------------------------------------------------

_logger = logging.getLogger(__name__)


class UpdateProgress(RemoteProgress):
    """Class providing detailed progress information to the caller who should
    derive from it and implement the
    :meth:`update(...) <git.util.RemoteProgress.update>` message."""

    CLONE, FETCH, UPDWKTREE = [1 << x for x in range(RemoteProgress._num_op_codes, RemoteProgress._num_op_codes + 3)]
    _num_op_codes: int = RemoteProgress._num_op_codes + 3

    __slots__ = ()


BEGIN = UpdateProgress.BEGIN
END = UpdateProgress.END
CLONE = UpdateProgress.CLONE
FETCH = UpdateProgress.FETCH
UPDWKTREE = UpdateProgress.UPDWKTREE


# IndexObject comes via the util module. It's a 'hacky' fix thanks to Python's import
# mechanism, which causes plenty of trouble if the only reason for packages and modules
# is refactoring - subpackages shouldn't depend on parent packages.
class Submodule(IndexObject, TraversableIterableObj):
    """Implements access to a git submodule. They are special in that their sha
    represents a commit in the submodule's repository which is to be checked out
    at the path of this instance.

    The submodule type does not have a string type associated with it, as it exists
    solely as a marker in the tree and index.

    All methods work in bare and non-bare repositories.
    """

    _id_attribute_ = "name"
    k_modules_file = ".gitmodules"
    k_head_option = "branch"
    k_head_default = "master"
    k_default_mode = stat.S_IFDIR | stat.S_IFLNK
    """Submodule flags. Submodules are directories with link-status."""

    type: Literal["submodule"] = "submodule"  # type: ignore[assignment]
    """This is a bogus type string for base class compatibility."""

    __slots__ = ("_parent_commit", "_url", "_branch_path", "_name", "__weakref__")

    _cache_attrs = ("path", "_url", "_branch_path")

    def __init__(
        self,
        repo: "Repo",
        binsha: bytes,
        mode: Union[int, None] = None,
        path: Union[PathLike, None] = None,
        name: Union[str, None] = None,
        parent_commit: Union["Commit", None] = None,
        url: Union[str, None] = None,
        branch_path: Union[PathLike, None] = None,
    ) -> None:
        """Initialize this instance with its attributes.

        We only document the parameters that differ from
        :class:`~git.objects.base.IndexObject`.

        :param repo:
            Our parent repository.

        :param binsha:
            Binary sha referring to a commit in the remote repository.
            See the `url` parameter.

        :param parent_commit:
            The :class:`~git.objects.commit.Commit` whose tree is supposed to contain
            the ``.gitmodules`` blob, or ``None`` to always point to the most recent
            commit. See :meth:`set_parent_commit` for details.

        :param url:
            The URL to the remote repository which is the submodule.

        :param branch_path:
            Full repository-relative path to ref to checkout when cloning the remote
            repository.
        """
        super().__init__(repo, binsha, mode, path)
        self.size = 0
        self._parent_commit = parent_commit
        if url is not None:
            self._url = url
        if branch_path is not None:
            self._branch_path = branch_path
        if name is not None:
            self._name = name

    def _set_cache_(self, attr: str) -> None:
        if attr in ("path", "_url", "_branch_path"):
            reader: SectionConstraint = self.config_reader()
            # Default submodule values.
            try:
                self.path = reader.get("path")
            except cp.NoSectionError as e:
                if self.repo.working_tree_dir is not None:
                    raise ValueError(
                        "This submodule instance does not exist anymore in '%s' file"
                        % osp.join(self.repo.working_tree_dir, ".gitmodules")
                    ) from e

            self._url = reader.get("url")
            # GitPython extension values - optional.
            self._branch_path = reader.get_value(self.k_head_option, git.Head.to_full_path(self.k_head_default))
        elif attr == "_name":
            raise AttributeError("Cannot retrieve the name of a submodule if it was not set initially")
        else:
            super()._set_cache_(attr)
        # END handle attribute name

    @classmethod
    def _get_intermediate_items(cls, item: "Submodule") -> IterableList["Submodule"]:
        """:return: All the submodules of our module repository"""
        try:
            return cls.list_items(item.module())
        except InvalidGitRepositoryError:
            return IterableList("")
        # END handle intermediate items

    @classmethod
    def _need_gitfile_submodules(cls, git: Git) -> bool:
        return git.version_info[:3] >= (1, 7, 5)

    def __eq__(self, other: Any) -> bool:
        """Compare with another submodule."""
        # We may only compare by name as this should be the ID they are hashed with.
        # Otherwise this type wouldn't be hashable.
        # return self.path == other.path and self.url == other.url and super().__eq__(other)
        return self._name == other._name

    def __ne__(self, other: object) -> bool:
        """Compare with another submodule for inequality."""
        return not (self == other)

    def __hash__(self) -> int:
        """Hash this instance using its logical id, not the sha."""
        return hash(self._name)

    def __str__(self) -> str:
        return self._name

    def __repr__(self) -> str:
        return "git.%s(name=%s, path=%s, url=%s, branch_path=%s)" % (
            type(self).__name__,
            self._name,
            self.path,
            self.url,
            self.branch_path,
        )

    @classmethod
    def _config_parser(
        cls, repo: "Repo", parent_commit: Union["Commit", None], read_only: bool
    ) -> SubmoduleConfigParser:
        """
        :return:
            Config parser constrained to our submodule in read or write mode

        :raise IOError:
            If the ``.gitmodules`` file cannot be found, either locally or in the
            repository at the given parent commit. Otherwise the exception would be
            delayed until the first access of the config parser.
        """
        parent_matches_head = True
        if parent_commit is not None:
            try:
                parent_matches_head = repo.head.commit == parent_commit
            except ValueError:
                # We are most likely in an empty repository, so the HEAD doesn't point
                # to a valid ref.
                pass
        # END handle parent_commit
        fp_module: Union[str, BytesIO]
        if not repo.bare and parent_matches_head and repo.working_tree_dir:
            fp_module = osp.join(repo.working_tree_dir, cls.k_modules_file)
        else:
            assert parent_commit is not None, "need valid parent_commit in bare repositories"
            try:
                fp_module = cls._sio_modules(parent_commit)
            except KeyError as e:
                raise IOError(
                    "Could not find %s file in the tree of parent commit %s" % (cls.k_modules_file, parent_commit)
                ) from e
            # END handle exceptions
        # END handle non-bare working tree

        if not read_only and (repo.bare or not parent_matches_head):
            raise ValueError("Cannot write blobs of 'historical' submodule configurations")
        # END handle writes of historical submodules

        return SubmoduleConfigParser(fp_module, read_only=read_only)

    def _clear_cache(self) -> None:
        """Clear the possibly changed values."""
        for name in self._cache_attrs:
            try:
                delattr(self, name)
            except AttributeError:
                pass
            # END try attr deletion
        # END for each name to delete

    @classmethod
    def _sio_modules(cls, parent_commit: "Commit") -> BytesIO:
        """
        :return:
            Configuration file as :class:`~io.BytesIO` - we only access it through the
            respective blob's data
        """
        sio = BytesIO(parent_commit.tree[cls.k_modules_file].data_stream.read())
        sio.name = cls.k_modules_file
        return sio

    def _config_parser_constrained(self, read_only: bool) -> SectionConstraint:
        """:return: Config parser constrained to our submodule in read or write mode"""
        try:
            pc = self.parent_commit
        except ValueError:
            pc = None
        # END handle empty parent repository
        parser = self._config_parser(self.repo, pc, read_only)
        parser.set_submodule(self)
        return SectionConstraint(parser, sm_section(self.name))

    @classmethod
    def _module_abspath(cls, parent_repo: "Repo", path: PathLike, name: str) -> PathLike:
        if cls._need_gitfile_submodules(parent_repo.git):
            return osp.join(parent_repo.git_dir, "modules", name)
        if parent_repo.working_tree_dir:
            return osp.join(parent_repo.working_tree_dir, path)
        raise NotADirectoryError()

    @classmethod
    def _clone_repo(
        cls,
        repo: "Repo",
        url: str,
        path: PathLike,
        name: str,
        allow_unsafe_options: bool = False,
        allow_unsafe_protocols: bool = False,
        **kwargs: Any,
    ) -> "Repo":
        """
        :return:
            :class:`~git.repo.base.Repo` instance of newly cloned repository.

        :param repo:
            Our parent repository.

        :param url:
            URL to clone from.

        :param path:
            Repository-relative path to the submodule checkout location.

        :param name:
            Canonical name of the submodule.

        :param allow_unsafe_protocols:
            Allow unsafe protocols to be used, like ``ext``.

        :param allow_unsafe_options:
            Allow unsafe options to be used, like ``--upload-pack``.

        :param kwargs:
            Additional arguments given to :manpage:`git-clone(1)`.
        """
        module_abspath = cls._module_abspath(repo, path, name)
        module_checkout_path = module_abspath
        if cls._need_gitfile_submodules(repo.git):
            kwargs["separate_git_dir"] = module_abspath
            module_abspath_dir = osp.dirname(module_abspath)
            if not osp.isdir(module_abspath_dir):
                os.makedirs(module_abspath_dir)
            module_checkout_path = osp.join(repo.working_tree_dir, path)  # type: ignore[arg-type]

        if url.startswith("../"):
            remote_name = cast("RemoteReference", repo.active_branch.tracking_branch()).remote_name
            repo_remote_url = repo.remote(remote_name).url
            url = os.path.join(repo_remote_url, url)

        clone = git.Repo.clone_from(
            url,
            module_checkout_path,
            allow_unsafe_options=allow_unsafe_options,
            allow_unsafe_protocols=allow_unsafe_protocols,
            **kwargs,
        )
        if cls._need_gitfile_submodules(repo.git):
            cls._write_git_file_and_module_config(module_checkout_path, module_abspath)

        return clone

    @classmethod
    def _to_relative_path(cls, parent_repo: "Repo", path: PathLike) -> PathLike:
        """:return: A path guaranteed to be relative to the given parent repository

        :raise ValueError:
            If path is not contained in the parent repository's working tree.
        """
        path = to_native_path_linux(path)
        if path.endswith("/"):
            path = path[:-1]
        # END handle trailing slash

        if osp.isabs(path) and parent_repo.working_tree_dir:
            working_tree_linux = to_native_path_linux(parent_repo.working_tree_dir)
            if not path.startswith(working_tree_linux):
                raise ValueError(
                    "Submodule checkout path '%s' needs to be within the parents repository at '%s'"
                    % (working_tree_linux, path)
                )
            path = path[len(working_tree_linux.rstrip("/")) + 1 :]
            if not path:
                raise ValueError("Absolute submodule path '%s' didn't yield a valid relative path" % path)
            # END verify converted relative path makes sense
        # END convert to a relative path

        return path

    @classmethod
    def _write_git_file_and_module_config(cls, working_tree_dir: PathLike, module_abspath: PathLike) -> None:
        """Write a ``.git`` file containing a (preferably) relative path to the actual
        git module repository.

        It is an error if the `module_abspath` cannot be made into a relative path,
        relative to the `working_tree_dir`.

        :note:
            This will overwrite existing files!

        :note:
            As we rewrite both the git file as well as the module configuration, we
            might fail on the configuration and will not roll back changes done to the
            git file. This should be a non-issue, but may easily be fixed if it becomes
            one.

        :param working_tree_dir:
            Directory to write the ``.git`` file into.

        :param module_abspath:
            Absolute path to the bare repository.
        """
        git_file = osp.join(working_tree_dir, ".git")
        rela_path = osp.relpath(module_abspath, start=working_tree_dir)
        if sys.platform == "win32" and osp.isfile(git_file):
            os.remove(git_file)
        with open(git_file, "wb") as fp:
            fp.write(("gitdir: %s" % rela_path).encode(defenc))

        with GitConfigParser(osp.join(module_abspath, "config"), read_only=False, merge_includes=False) as writer:
            writer.set_value(
                "core",
                "worktree",
                to_native_path_linux(osp.relpath(working_tree_dir, start=module_abspath)),
            )

    # { Edit Interface

    @classmethod
    def add(
        cls,
        repo: "Repo",
        name: str,
        path: PathLike,
        url: Union[str, None] = None,
        branch: Union[str, None] = None,
        no_checkout: bool = False,
        depth: Union[int, None] = None,
        env: Union[Mapping[str, str], None] = None,
        clone_multi_options: Union[Sequence[TBD], None] = None,
        allow_unsafe_options: bool = False,
        allow_unsafe_protocols: bool = False,
    ) -> "Submodule":
        """Add a new submodule to the given repository. This will alter the index as
        well as the ``.gitmodules`` file, but will not create a new commit. If the
        submodule already exists, no matter if the configuration differs from the one
        provided, the existing submodule will be returned.

        :param repo:
            Repository instance which should receive the submodule.

        :param name:
            The name/identifier for the submodule.

        :param path:
            Repository-relative or absolute path at which the submodule should be
            located.
            It will be created as required during the repository initialization.

        :param url:
            ``git clone ...``-compatible URL. See :manpage:`git-clone(1)` for more
            information. If ``None``, the repository is assumed to exist, and the URL of
            the first remote is taken instead. This is useful if you want to make an
            existing repository a submodule of another one.

        :param branch:
            Name of branch at which the submodule should (later) be checked out. The
            given branch must exist in the remote repository, and will be checked out
            locally as a tracking branch.
            It will only be written into the configuration if it not ``None``, which is
            when the checked out branch will be the one the remote HEAD pointed to.
            The result you get in these situation is somewhat fuzzy, and it is
            recommended to specify at least ``master`` here.
            Examples are ``master`` or ``feature/new``.

        :param no_checkout:
            If ``True``, and if the repository has to be cloned manually, no checkout
            will be performed.

        :param depth:
            Create a shallow clone with a history truncated to the specified number of
            commits.

        :param env:
            Optional dictionary containing the desired environment variables.

            Note: Provided variables will be used to update the execution environment
            for ``git``. If some variable is not specified in `env` and is defined in
            attr:`os.environ`, the value from attr:`os.environ` will be used. If you
            want to unset some variable, consider providing an empty string as its
            value.

        :param clone_multi_options:
            A list of clone options. Please see
            :meth:`Repo.clone <git.repo.base.Repo.clone>` for details.

        :param allow_unsafe_protocols:
            Allow unsafe protocols to be used, like ``ext``.

        :param allow_unsafe_options:
            Allow unsafe options to be used, like ``--upload-pack``.

        :return:
            The newly created :class:`Submodule` instance.

        :note:
            Works atomically, such that no change will be done if, for example, the
            repository update fails.
        """
        if repo.bare:
            raise InvalidGitRepositoryError("Cannot add submodules to bare repositories")
        # END handle bare repos

        path = cls._to_relative_path(repo, path)

        # Ensure we never put backslashes into the URL, as might happen on Windows.
        if url is not None:
            url = to_native_path_linux(url)
        # END ensure URL correctness

        # INSTANTIATE INTERMEDIATE SM
        sm = cls(
            repo,
            cls.NULL_BIN_SHA,
            cls.k_default_mode,
            path,
            name,
            url="invalid-temporary",
        )
        if sm.exists():
            # Reretrieve submodule from tree.
            try:
                sm = repo.head.commit.tree[os.fspath(path)]
                sm._name = name
                return sm
            except KeyError:
                # Could only be in index.
                index = repo.index
                entry = index.entries[index.entry_key(path, 0)]
                sm.binsha = entry.binsha
                return sm
            # END handle exceptions
        # END handle existing

        # fake-repo - we only need the functionality on the branch instance.
        br = git.Head(repo, git.Head.to_full_path(str(branch) or cls.k_head_default))
        has_module = sm.module_exists()
        branch_is_default = branch is None
        if has_module and url is not None:
            if url not in [r.url for r in sm.module().remotes]:
                raise ValueError(
                    "Specified URL '%s' does not match any remote url of the repository at '%s'" % (url, sm.abspath)
                )
            # END check url
        # END verify urls match

        mrepo: Union[Repo, None] = None

        if url is None:
            if not has_module:
                raise ValueError("A URL was not given and a repository did not exist at %s" % path)
            # END check url
            mrepo = sm.module()
            # assert isinstance(mrepo, git.Repo)
            urls = [r.url for r in mrepo.remotes]
            if not urls:
                raise ValueError("Didn't find any remote url in repository at %s" % sm.abspath)
            # END verify we have url
            url = urls[0]
        else:
            # Clone new repo.
            kwargs: Dict[str, Union[bool, int, str, Sequence[TBD]]] = {"n": no_checkout}
            if not branch_is_default:
                kwargs["b"] = br.name
            # END setup checkout-branch

            if depth:
                if isinstance(depth, int):
                    kwargs["depth"] = depth
                else:
                    raise ValueError("depth should be an integer")
            if clone_multi_options:
                kwargs["multi_options"] = clone_multi_options

            # _clone_repo(cls, repo, url, path, name, **kwargs):
            mrepo = cls._clone_repo(
                repo,
                url,
                path,
                name,
                env=env,
                allow_unsafe_options=allow_unsafe_options,
                allow_unsafe_protocols=allow_unsafe_protocols,
                **kwargs,
            )
        # END verify url

        ## See #525 for ensuring git URLs in config-files are valid under Windows.
        url = Git.polish_url(url)

        # It's important to add the URL to the parent config, to let `git submodule` know.
        # Otherwise there is a '-' character in front of the submodule listing:
        #  a38efa84daef914e4de58d1905a500d8d14aaf45 mymodule (v0.9.0-1-ga38efa8)
        # -a38efa84daef914e4de58d1905a500d8d14aaf45 submodules/intermediate/one
        writer: Union[GitConfigParser, SectionConstraint]

        with sm.repo.config_writer() as writer:
            writer.set_value(sm_section(name), "url", url)

        # Update configuration and index.
        index = sm.repo.index
        with sm.config_writer(index=index, write=False) as writer:
            writer.set_value("url", url)
            writer.set_value("path", path)

            sm._url = url
            if not branch_is_default:
                # Store full path.
                writer.set_value(cls.k_head_option, br.path)
                sm._branch_path = br.path

        # We deliberately assume that our head matches our index!
        if mrepo:
            sm.binsha = mrepo.head.commit.binsha
        index.add([sm], write=True)

        return sm

    def update(
        self,
        recursive: bool = False,
        init: bool = True,
        to_latest_revision: bool = False,
        progress: Union["UpdateProgress", None] = None,
        dry_run: bool = False,
        force: bool = False,
        keep_going: bool = False,
        env: Union[Mapping[str, str], None] = None,
        clone_multi_options: Union[Sequence[TBD], None] = None,
        allow_unsafe_options: bool = False,
        allow_unsafe_protocols: bool = False,
    ) -> "Submodule":
        """Update the repository of this submodule to point to the checkout we point at
        with the binsha of this instance.

        :param recursive:
            If ``True``, we will operate recursively and update child modules as well.

        :param init:
            If ``True``, the module repository will be cloned into place if necessary.

        :param to_latest_revision:
            If ``True``, the submodule's sha will be ignored during checkout. Instead,
            the remote will be fetched, and the local tracking branch updated. This only
            works if we have a local tracking branch, which is the case if the remote
            repository had a master branch, or if the ``branch`` option was specified
            for this submodule and the branch existed remotely.

        :param progress:
            :class:`UpdateProgress` instance, or ``None`` if no progress should be
            shown.

        :param dry_run:
            If ``True``, the operation will only be simulated, but not performed.
            All performed operations are read-only.

        :param force:
            If ``True``, we may reset heads even if the repository in question is dirty.
            Additionally we will be allowed to set a tracking branch which is ahead of
            its remote branch back into the past or the location of the remote branch.
            This will essentially 'forget' commits.

            If ``False``, local tracking branches that are in the future of their
            respective remote branches will simply not be moved.

        :param keep_going:
            If ``True``, we will ignore but log all errors, and keep going recursively.
            Unless `dry_run` is set as well, `keep_going` could cause
            subsequent/inherited errors you wouldn't see otherwise.
            In conjunction with `dry_run`, it can be useful to anticipate all errors
            when updating submodules.

        :param env:
            Optional dictionary containing the desired environment variables.

            Note: Provided variables will be used to update the execution environment
            for ``git``. If some variable is not specified in `env` and is defined in
            attr:`os.environ`, value from attr:`os.environ` will be used.

            If you want to unset some variable, consider providing the empty string as
            its value.

        :param clone_multi_options:
            List of :manpage:`git-clone(1)` options.
            Please see :meth:`Repo.clone <git.repo.base.Repo.clone>` for details.
            They only take effect with the `init` option.

        :param allow_unsafe_protocols:
            Allow unsafe protocols to be used, like ``ext``.

        :param allow_unsafe_options:
            Allow unsafe options to be used, like ``--upload-pack``.

        :note:
            Does nothing in bare repositories.

        :note:
            This method is definitely not atomic if `recursive` is ``True``.

        :return:
            self
        """
        if self.repo.bare:
            return self
        # END pass in bare mode

        if progress is None:
            progress = UpdateProgress()
        # END handle progress
        prefix = ""
        if dry_run:
            prefix = "DRY-RUN: "
        # END handle prefix

        # To keep things plausible in dry-run mode.
        if dry_run:
            mrepo = None
        # END init mrepo

        try:
            # ENSURE REPO IS PRESENT AND UP-TO-DATE
            #######################################
            try:
                mrepo = self.module()
                rmts = mrepo.remotes
                len_rmts = len(rmts)
                for i, remote in enumerate(rmts):
                    op = FETCH
                    if i == 0:
                        op |= BEGIN
                    # END handle start

                    progress.update(
                        op,
                        i,
                        len_rmts,
                        prefix + "Fetching remote %s of submodule %r" % (remote, self.name),
                    )
                    # ===============================
                    if not dry_run:
                        remote.fetch(progress=progress)
                    # END handle dry-run
                    # ===============================
                    if i == len_rmts - 1:
                        op |= END
                    # END handle end
                    progress.update(
                        op,
                        i,
                        len_rmts,
                        prefix + "Done fetching remote of submodule %r" % self.name,
                    )
                # END fetch new data
            except InvalidGitRepositoryError:
                mrepo = None
                if not init:
                    return self
                # END early abort if init is not allowed

                # There is no git-repository yet - but delete empty paths.
                checkout_module_abspath = self.abspath
                if not dry_run and osp.isdir(checkout_module_abspath):
                    try:
                        os.rmdir(checkout_module_abspath)
                    except OSError as e:
                        raise OSError(
                            "Module directory at %r does already exist and is non-empty" % checkout_module_abspath
                        ) from e
                    # END handle OSError
                # END handle directory removal

                # Don't check it out at first - nonetheless it will create a local
                # branch according to the remote-HEAD if possible.
                progress.update(
                    BEGIN | CLONE,
                    0,
                    1,
                    prefix
                    + "Cloning url '%s' to '%s' in submodule %r" % (self.url, checkout_module_abspath, self.name),
                )
                if not dry_run:
                    if self.url.startswith("."):
                        url = urllib.parse.urljoin(self.repo.remotes.origin.url + "/", self.url)
                    else:
                        url = self.url
                    mrepo = self._clone_repo(
                        self.repo,
                        url,
                        self.path,
                        self.name,
                        n=True,
                        env=env,
                        multi_options=clone_multi_options,
                        allow_unsafe_options=allow_unsafe_options,
                        allow_unsafe_protocols=allow_unsafe_protocols,
                    )
                # END handle dry-run
                progress.update(
                    END | CLONE,
                    0,
                    1,
                    prefix + "Done cloning to %s" % checkout_module_abspath,
                )

                if not dry_run:
                    # See whether we have a valid branch to check out.
                    try:
                        mrepo = cast("Repo", mrepo)
                        # Find a remote which has our branch - we try to be flexible.
                        remote_branch = find_first_remote_branch(mrepo.remotes, self.branch_name)
                        local_branch = mkhead(mrepo, self.branch_path)

                        # Have a valid branch, but no checkout - make sure we can figure
                        # that out by marking the commit with a null_sha.
                        local_branch.set_object(Object(mrepo, self.NULL_BIN_SHA))
                        # END initial checkout + branch creation

                        # Make sure HEAD is not detached.
                        mrepo.head.set_reference(
                            local_branch,
                            logmsg="submodule: attaching head to %s" % local_branch,
                        )
                        mrepo.head.reference.set_tracking_branch(remote_branch)
                    except (IndexError, InvalidGitRepositoryError):
                        _logger.warning("Failed to checkout tracking branch %s", self.branch_path)
                    # END handle tracking branch

                    # NOTE: Have to write the repo config file as well, otherwise the
                    # default implementation will be offended and not update the
                    # repository. Maybe this is a good way to ensure it doesn't get into
                    # our way, but we want to stay backwards compatible too... It's so
                    # redundant!
                    with self.repo.config_writer() as writer:
                        writer.set_value(sm_section(self.name), "url", self.url)
                # END handle dry_run
            # END handle initialization

            # DETERMINE SHAS TO CHECK OUT
            #############################
            binsha = self.binsha
            hexsha = self.hexsha
            if mrepo is not None:
                # mrepo is only set if we are not in dry-run mode or if the module
                # existed.
                is_detached = mrepo.head.is_detached
            # END handle dry_run

            if mrepo is not None and to_latest_revision:
                msg_base = "Cannot update to latest revision in repository at %r as " % mrepo.working_dir
                if not is_detached:
                    rref = mrepo.head.reference.tracking_branch()
                    if rref is not None:
                        rcommit = rref.commit
                        binsha = rcommit.binsha
                        hexsha = rcommit.hexsha
                    else:
                        _logger.error(
                            "%s a tracking branch was not set for local branch '%s'",
                            msg_base,
                            mrepo.head.reference,
                        )
                    # END handle remote ref
                else:
                    _logger.error("%s there was no local tracking branch", msg_base)
                # END handle detached head
            # END handle to_latest_revision option

            # Update the working tree.
            # Handles dry_run.
            if mrepo is not None and mrepo.head.commit.binsha != binsha:
                # We must ensure that our destination sha (the one to point to) is in
                # the future of our current head. Otherwise, we will reset changes that
                # might have been done on the submodule, but were not yet pushed. We
                # also handle the case that history has been rewritten, leaving no
                # merge-base. In that case we behave conservatively, protecting possible
                # changes the user had done.
                may_reset = True
                if mrepo.head.commit.binsha != self.NULL_BIN_SHA:
                    base_commit = mrepo.merge_base(mrepo.head.commit, hexsha)
                    if len(base_commit) == 0 or (base_commit[0] is not None and base_commit[0].hexsha == hexsha):
                        if force:
                            msg = "Will force checkout or reset on local branch that is possibly in the future of"
                            msg += " the commit it will be checked out to, effectively 'forgetting' new commits"
                            _logger.debug(msg)
                        else:
                            msg = "Skipping %s on branch '%s' of submodule repo '%s' as it contains un-pushed commits"
                            msg %= (
                                is_detached and "checkout" or "reset",
                                mrepo.head,
                                mrepo,
                            )
                            _logger.info(msg)
                            may_reset = False
                        # END handle force
                    # END handle if we are in the future

                    if may_reset and not force and mrepo.is_dirty(index=True, working_tree=True, untracked_files=True):
                        raise RepositoryDirtyError(mrepo, "Cannot reset a dirty repository")
                    # END handle force and dirty state
                # END handle empty repo

                # END verify future/past
                progress.update(
                    BEGIN | UPDWKTREE,
                    0,
                    1,
                    prefix
                    + "Updating working tree at %s for submodule %r to revision %s" % (self.path, self.name, hexsha),
                )

                if not dry_run and may_reset:
                    if is_detached:
                        # NOTE: For now we force. The user is not supposed to change
                        # detached submodules anyway. Maybe at some point this becomes
                        # an option, to properly handle user modifications - see below
                        # for future options regarding rebase and merge.
                        mrepo.git.checkout(hexsha, force=force)
                    else:
                        mrepo.head.reset(hexsha, index=True, working_tree=True)
                    # END handle checkout
                # If we may reset/checkout.
                progress.update(
                    END | UPDWKTREE,
                    0,
                    1,
                    prefix + "Done updating working tree for submodule %r" % self.name,
                )
            # END update to new commit only if needed
        except Exception as err:
            if not keep_going:
                raise
            _logger.error(str(err))
        # END handle keep_going

        # HANDLE RECURSION
        ##################
        if recursive:
            # In dry_run mode, the module might not exist.
            if mrepo is not None:
                for submodule in self.iter_items(self.module()):
                    submodule.update(
                        recursive,
                        init,
                        to_latest_revision,
                        progress=progress,
                        dry_run=dry_run,
                        force=force,
                        keep_going=keep_going,
                    )
                # END handle recursive update
            # END handle dry run
        # END for each submodule

        return self

    @unbare_repo
    def move(self, module_path: PathLike, configuration: bool = True, module: bool = True) -> "Submodule":
        """Move the submodule to a another module path. This involves physically moving
        the repository at our current path, changing the configuration, as well as
        adjusting our index entry accordingly.

        :param module_path:
            The path to which to move our module in the parent repository's working
            tree, given as repository-relative or absolute path. Intermediate
            directories will be created accordingly. If the path already exists, it must
            be empty. Trailing (back)slashes are removed automatically.

        :param configuration:
            If ``True``, the configuration will be adjusted to let the submodule point
            to the given path.

        :param module:
            If ``True``, the repository managed by this submodule will be moved as well.
            If ``False``, we don't move the submodule's checkout, which may leave the
            parent repository in an inconsistent state.

        :return:
            self

        :raise ValueError:
            If the module path existed and was not empty, or was a file.

        :note:
            Currently the method is not atomic, and it could leave the repository in an
            inconsistent state if a sub-step fails for some reason.
        """
        if module + configuration < 1:
            raise ValueError("You must specify to move at least the module or the configuration of the submodule")
        # END handle input

        module_checkout_path = self._to_relative_path(self.repo, module_path)

        # VERIFY DESTINATION
        if module_checkout_path == self.path:
            return self
        # END handle no change

        module_checkout_abspath = join_path_native(str(self.repo.working_tree_dir), module_checkout_path)
        if osp.isfile(module_checkout_abspath):
            raise ValueError("Cannot move repository onto a file: %s" % module_checkout_abspath)
        # END handle target files

        index = self.repo.index
        tekey = index.entry_key(module_checkout_path, 0)
        # if the target item already exists, fail
        if configuration and tekey in index.entries:
            raise ValueError("Index entry for target path did already exist")
        # END handle index key already there

        # Remove existing destination.
        if module:
            if osp.exists(module_checkout_abspath):
                if len(os.listdir(module_checkout_abspath)):
                    raise ValueError("Destination module directory was not empty")
                # END handle non-emptiness

                if osp.islink(module_checkout_abspath):
                    os.remove(module_checkout_abspath)
                else:
                    os.rmdir(module_checkout_abspath)
                # END handle link
            else:
                # Recreate parent directories.
                # NOTE: renames() does that now.
                pass
            # END handle existence
        # END handle module

        # Move the module into place if possible.
        cur_path = self.abspath
        renamed_module = False
        if module and osp.exists(cur_path):
            os.renames(cur_path, module_checkout_abspath)
            renamed_module = True

            if osp.isfile(osp.join(module_checkout_abspath, ".git")):
                module_abspath = self._module_abspath(self.repo, self.path, self.name)
                self._write_git_file_and_module_config(module_checkout_abspath, module_abspath)
            # END handle git file rewrite
        # END move physical module

        # Rename the index entry - we have to manipulate the index directly as git-mv
        # cannot be used on submodules... yeah.
        previous_sm_path = self.path
        try:
            if configuration:
                try:
                    ekey = index.entry_key(self.path, 0)
                    entry = index.entries[ekey]
                    del index.entries[ekey]
                    nentry = git.IndexEntry(entry[:3] + (module_checkout_path,) + entry[4:])
                    index.entries[tekey] = nentry
                except KeyError as e:
                    raise InvalidGitRepositoryError("Submodule's entry at %r did not exist" % (self.path)) from e
                # END handle submodule doesn't exist

                # Update configuration.
                with self.config_writer(index=index) as writer:  # Auto-write.
                    writer.set_value("path", module_checkout_path)
                    self.path = module_checkout_path
            # END handle configuration flag
        except Exception:
            if renamed_module:
                os.renames(module_checkout_abspath, cur_path)
            # END undo module renaming
            raise
        # END handle undo rename

        # Auto-rename submodule if its name was 'default', that is, the checkout
        # directory.
        if previous_sm_path == self.name:
            self.rename(module_checkout_path)

        return self

    @unbare_repo
    def remove(
        self,
        module: bool = True,
        force: bool = False,
        configuration: bool = True,
        dry_run: bool = False,
    ) -> "Submodule":
        """Remove this submodule from the repository. This will remove our entry
        from the ``.gitmodules`` file and the entry in the ``.git/config`` file.

        :param module:
            If ``True``, the checked out module we point to will be deleted as well. If
            that module is currently on a commit outside any branch in the remote, or if
            it is ahead of its tracking branch, or if there are modified or untracked
            files in its working tree, then the removal will fail. In case the removal
            of the repository fails for these reasons, the submodule status will not
            have been altered.

            If this submodule has child modules of its own, these will be deleted prior
            to touching the direct submodule.

        :param force:
            Enforces the deletion of the module even though it contains modifications.
            This basically enforces a brute-force file system based deletion.

        :param configuration:
            If ``True``, the submodule is deleted from the configuration, otherwise it
            isn't. Although this should be enabled most of the time, this flag enables
            you to safely delete the repository of your submodule.

        :param dry_run:
            If ``True``, we will not actually do anything, but throw the errors we would
            usually throw.

        :return:
            self

        :note:
            Doesn't work in bare repositories.

        :note:
            Doesn't work atomically, as failure to remove any part of the submodule will
            leave an inconsistent state.

        :raise git.exc.InvalidGitRepositoryError:
            Thrown if the repository cannot be deleted.

        :raise OSError:
            If directories or files could not be removed.
        """
        if not (module or configuration):
            raise ValueError("Need to specify to delete at least the module, or the configuration")
        # END handle parameters

        # Recursively remove children of this submodule.
        nc = 0
        for csm in self.children():
            nc += 1
            csm.remove(module, force, configuration, dry_run)
            del csm

        if configuration and not dry_run and nc > 0:
            # Ensure we don't leave the parent repository in a dirty state, and commit
            # our changes. It's important for recursive, unforced, deletions to work as
            # expected.
            self.module().index.commit("Removed at least one of child-modules of '%s'" % self.name)
        # END handle recursion

        # DELETE REPOSITORY WORKING TREE
        ################################
        if module and self.module_exists():
            mod = self.module()
            git_dir = mod.git_dir
            if force:
                # Take the fast lane and just delete everything in our module path.
                # TODO: If we run into permission problems, we have a highly
                # inconsistent state. Delete the .git folders last, start with the
                # submodules first.
                mp = self.abspath
                method: Union[None, Callable[[PathLike], None]] = None
                if osp.islink(mp):
                    method = os.remove
                elif osp.isdir(mp):
                    method = rmtree
                elif osp.exists(mp):
                    raise AssertionError("Cannot forcibly delete repository as it was neither a link, nor a directory")
                # END handle brutal deletion
                if not dry_run:
                    assert method
                    method(mp)
                # END apply deletion method
            else:
                # Verify we may delete our module.
                if mod.is_dirty(index=True, working_tree=True, untracked_files=True):
                    raise InvalidGitRepositoryError(
                        "Cannot delete module at %s with any modifications, unless force is specified"
                        % mod.working_tree_dir
                    )
                # END check for dirt

                # Figure out whether we have new commits compared to the remotes.
                # NOTE: If the user pulled all the time, the remote heads might not have
                # been updated, so commits coming from the remote look as if they come
                # from us. But we stay strictly read-only and don't fetch beforehand.
                for remote in mod.remotes:
                    num_branches_with_new_commits = 0
                    rrefs = remote.refs
                    for rref in rrefs:
                        num_branches_with_new_commits += len(mod.git.cherry(rref)) != 0
                    # END for each remote ref
                    # Not a single remote branch contained all our commits.
                    if len(rrefs) and num_branches_with_new_commits == len(rrefs):
                        raise InvalidGitRepositoryError(
                            "Cannot delete module at %s as there are new commits" % mod.working_tree_dir
                        )
                    # END handle new commits
                    # We have to manually delete some references to allow resources to
                    # be cleaned up immediately when we are done with them, because
                    # Python's scoping is no more granular than the whole function (loop
                    # bodies are not scopes). When the objects stay alive longer, they
                    # can keep handles open. On Windows, this is a problem.
                    if len(rrefs):
                        del rref  # skipcq: PYL-W0631
                    # END handle remotes
                    del rrefs
                    del remote
                # END for each remote

                # Finally delete our own submodule.
                if not dry_run:
                    self._clear_cache()
                    wtd = mod.working_tree_dir
                    del mod  # Release file-handles (Windows).
                    gc.collect()
                    rmtree(str(wtd))
                # END delete tree if possible
            # END handle force

            if not dry_run and osp.isdir(git_dir):
                self._clear_cache()
                rmtree(git_dir)
            # END handle separate bare repository
        # END handle module deletion

        # Void our data so as not to delay invalid access.
        if not dry_run:
            self._clear_cache()

        # DELETE CONFIGURATION
        ######################
        if configuration and not dry_run:
            # First the index-entry.
            parent_index = self.repo.index
            try:
                del parent_index.entries[parent_index.entry_key(self.path, 0)]
            except KeyError:
                pass
            # END delete entry
            parent_index.write()

            # Now git config - we need the config intact, otherwise we can't query
            # information anymore.

            with self.repo.config_writer() as gcp_writer:
                gcp_writer.remove_section(sm_section(self.name))

            with self.config_writer() as sc_writer:
                sc_writer.remove_section()
        # END delete configuration

        return self

    def set_parent_commit(self, commit: Union[Commit_ish, str, None], check: bool = True) -> "Submodule":
        """Set this instance to use the given commit whose tree is supposed to
        contain the ``.gitmodules`` blob.

        :param commit:
            Commit-ish reference pointing at the root tree, or ``None`` to always point
            to the most recent commit.

        :param check:
            If ``True``, relatively expensive checks will be performed to verify
            validity of the submodule.

        :raise ValueError:
            If the commit's tree didn't contain the ``.gitmodules`` blob.

        :raise ValueError:
            If the parent commit didn't store this submodule under the current path.

        :return:
            self
        """
        if commit is None:
            self._parent_commit = None
            return self
        # END handle None
        pcommit = self.repo.commit(commit)
        pctree = pcommit.tree
        if self.k_modules_file not in pctree:
            raise ValueError("Tree of commit %s did not contain the %s file" % (commit, self.k_modules_file))
        # END handle exceptions

        prev_pc = self._parent_commit
        self._parent_commit = pcommit

        if check:
            parser = self._config_parser(self.repo, self._parent_commit, read_only=True)
            if not parser.has_section(sm_section(self.name)):
                self._parent_commit = prev_pc
                raise ValueError("Submodule at path %r did not exist in parent commit %s" % (self.path, commit))
            # END handle submodule did not exist
        # END handle checking mode

        # Update our sha, it could have changed.
        # If check is False, we might see a parent-commit that doesn't even contain the
        # submodule anymore. in that case, mark our sha as being NULL.
        try:
            self.binsha = pctree[str(self.path)].binsha
        except KeyError:
            self.binsha = self.NULL_BIN_SHA

        self._clear_cache()
        return self

    @unbare_repo
    def config_writer(
        self, index: Union["IndexFile", None] = None, write: bool = True
    ) -> SectionConstraint["SubmoduleConfigParser"]:
        """
        :return:
            A config writer instance allowing you to read and write the data belonging
            to this submodule into the ``.gitmodules`` file.

        :param index:
            If not ``None``, an :class:`~git.index.base.IndexFile` instance which should
            be written. Defaults to the index of the :class:`Submodule`'s parent
            repository.

        :param write:
            If ``True``, the index will be written each time a configuration value changes.

        :note:
            The parameters allow for a more efficient writing of the index, as you can
            pass in a modified index on your own, prevent automatic writing, and write
            yourself once the whole operation is complete.

        :raise ValueError:
            If trying to get a writer on a parent_commit which does not match the
            current head commit.

        :raise IOError:
            If the ``.gitmodules`` file/blob could not be read.
        """
        writer = self._config_parser_constrained(read_only=False)
        if index is not None:
            writer.config._index = index
        writer.config._auto_write = write
        return writer

    @unbare_repo
    def rename(self, new_name: str) -> "Submodule":
        """Rename this submodule.

        :note:
            This method takes care of renaming the submodule in various places, such as:

            * ``$parent_git_dir / config``
            * ``$working_tree_dir / .gitmodules``
            * (git >= v1.8.0: move submodule repository to new name)

        As ``.gitmodules`` will be changed, you would need to make a commit afterwards.
        The changed ``.gitmodules`` file will already be added to the index.

        :return:
            This :class:`Submodule` instance
        """
        if self.name == new_name:
            return self

        # .git/config
        with self.repo.config_writer() as pw:
            # As we ourselves didn't write anything about submodules into the parent
            # .git/config, we will not require it to exist, and just ignore missing
            # entries.
            if pw.has_section(sm_section(self.name)):
                pw.rename_section(sm_section(self.name), sm_section(new_name))

        # .gitmodules
        with self.config_writer(write=True).config as cw:
            cw.rename_section(sm_section(self.name), sm_section(new_name))

        self._name = new_name

        # .git/modules
        mod = self.module()
        if mod.has_separate_working_tree():
            destination_module_abspath = self._module_abspath(self.repo, self.path, new_name)
            source_dir = mod.git_dir
            # Let's be sure the submodule name is not so obviously tied to a directory.
            if str(destination_module_abspath).startswith(str(mod.git_dir)):
                tmp_dir = self._module_abspath(self.repo, self.path, str(uuid.uuid4()))
                os.renames(source_dir, tmp_dir)
                source_dir = tmp_dir
            # END handle self-containment
            os.renames(source_dir, destination_module_abspath)
            if mod.working_tree_dir:
                self._write_git_file_and_module_config(mod.working_tree_dir, destination_module_abspath)
        # END move separate git repository

        return self

    # } END edit interface

    # { Query Interface

    @unbare_repo
    def module(self) -> "Repo":
        """
        :return:
            :class:`~git.repo.base.Repo` instance initialized from the repository at our
            submodule path

        :raise git.exc.InvalidGitRepositoryError:
            If a repository was not available.
            This could also mean that it was not yet initialized.
        """
        module_checkout_abspath = self.abspath
        try:
            repo = git.Repo(module_checkout_abspath)
            if repo != self.repo:
                return repo
            # END handle repo uninitialized
        except (InvalidGitRepositoryError, NoSuchPathError) as e:
            raise InvalidGitRepositoryError("No valid repository at %s" % module_checkout_abspath) from e
        else:
            raise InvalidGitRepositoryError("Repository at %r was not yet checked out" % module_checkout_abspath)
        # END handle exceptions

    def module_exists(self) -> bool:
        """
        :return:
            ``True`` if our module exists and is a valid git repository.
            See the :meth:`module` method.
        """
        try:
            self.module()
            return True
        except Exception:
            return False
        # END handle exception

    def exists(self) -> bool:
        """
        :return:
            ``True`` if the submodule exists, ``False`` otherwise.
            Please note that a submodule may exist (in the ``.gitmodules`` file) even
            though its module doesn't exist on disk.
        """
        # Keep attributes for later, and restore them if we have no valid data.
        # This way we do not actually alter the state of the object.
        loc = locals()
        for attr in self._cache_attrs:
            try:
                if hasattr(self, attr):
                    loc[attr] = getattr(self, attr)
                # END if we have the attribute cache
            except (cp.NoSectionError, ValueError):
                # On PY3, this can happen apparently... don't know why this doesn't
                # happen on PY2.
                pass
        # END for each attr
        self._clear_cache()

        try:
            try:
                self.path  # noqa: B018
                return True
            except Exception:
                return False
            # END handle exceptions
        finally:
            for attr in self._cache_attrs:
                if attr in loc:
                    setattr(self, attr, loc[attr])
                # END if we have a cache
            # END reapply each attribute
        # END handle object state consistency

    @property
    def branch(self) -> "Head":
        """
        :return:
            The branch instance that we are to checkout

        :raise git.exc.InvalidGitRepositoryError:
            If our module is not yet checked out.
        """
        return mkhead(self.module(), self._branch_path)

    @property
    def branch_path(self) -> PathLike:
        """
        :return:
            Full repository-relative path as string to the branch we would checkout from
            the remote and track
        """
        return self._branch_path

    @property
    def branch_name(self) -> str:
        """
        :return:
            The name of the branch, which is the shortest possible branch name
        """
        # Use an instance method, for this we create a temporary Head instance which
        # uses a repository that is available at least (it makes no difference).
        return git.Head(self.repo, self._branch_path).name

    @property
    def url(self) -> str:
        """:return: The url to the repository our submodule's repository refers to"""
        return self._url

    @property
    def parent_commit(self) -> "Commit":
        """
        :return:
            :class:`~git.objects.commit.Commit` instance with the tree containing the
            ``.gitmodules`` file

        :note:
            Will always point to the current head's commit if it was not set explicitly.
        """
        if self._parent_commit is None:
            return self.repo.commit()
        return self._parent_commit

    @property
    def name(self) -> str:
        """
        :return:
            The name of this submodule. It is used to identify it within the
            ``.gitmodules`` file.

        :note:
            By default, this is the name is the path at which to find the submodule, but
            in GitPython it should be a unique identifier similar to the identifiers
            used for remotes, which allows to change the path of the submodule easily.
        """
        return self._name

    def config_reader(self) -> SectionConstraint[SubmoduleConfigParser]:
        """
        :return:
            ConfigReader instance which allows you to query the configuration values of
            this submodule, as provided by the ``.gitmodules`` file.

        :note:
            The config reader will actually read the data directly from the repository
            and thus does not need nor care about your working tree.

        :note:
            Should be cached by the caller and only kept as long as needed.

        :raise IOError:
            If the ``.gitmodules`` file/blob could not be read.
        """
        return self._config_parser_constrained(read_only=True)

    def children(self) -> IterableList["Submodule"]:
        """
        :return:
            IterableList(Submodule, ...) An iterable list of :class:`Submodule`
            instances which are children of this submodule or 0 if the submodule is not
            checked out.
        """
        return self._get_intermediate_items(self)

    # } END query interface

    # { Iterable Interface

    @classmethod
    def iter_items(
        cls,
        repo: "Repo",
        parent_commit: Union[Commit_ish, str] = "HEAD",
        *args: Any,
        **kwargs: Any,
    ) -> Iterator["Submodule"]:
        """
        :return:
            Iterator yielding :class:`Submodule` instances available in the given
            repository
        """
        try:
            pc = repo.commit(parent_commit)  # Parent commit instance
            parser = cls._config_parser(repo, pc, read_only=True)
        except (IOError, BadName):
            return
        # END handle empty iterator

        for sms in parser.sections():
            n = sm_name(sms)
            p = parser.get(sms, "path")
            u = parser.get(sms, "url")
            b = cls.k_head_default
            if parser.has_option(sms, cls.k_head_option):
                b = str(parser.get(sms, cls.k_head_option))
            # END handle optional information

            # Get the binsha.
            index = repo.index
            try:
                rt = pc.tree  # Root tree
                sm = rt[p]
            except KeyError:
                # Try the index, maybe it was just added.
                try:
                    entry = index.entries[index.entry_key(p, 0)]
                    sm = Submodule(repo, entry.binsha, entry.mode, entry.path)
                except KeyError:
                    # The submodule doesn't exist, probably it wasn't removed from the
                    # .gitmodules file.
                    continue
                # END handle keyerror
            # END handle critical error

            # Make sure we are looking at a submodule object.
            if type(sm) is not git.objects.submodule.base.Submodule:
                continue

            # Fill in remaining info - saves time as it doesn't have to be parsed again.
            sm._name = n
            if pc != repo.commit():
                sm._parent_commit = pc
            # END set only if not most recent!
            sm._branch_path = git.Head.to_full_path(b)
            sm._url = u

            yield sm
        # END for each section

    # } END iterable interface
