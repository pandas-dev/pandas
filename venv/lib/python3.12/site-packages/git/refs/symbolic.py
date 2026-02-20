# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

__all__ = ["SymbolicReference"]

import os
from pathlib import Path

from gitdb.exc import BadName, BadObject

from git.compat import defenc
from git.objects.base import Object
from git.objects.commit import Commit
from git.refs.log import RefLog
from git.util import (
    LockedFD,
    assure_directory_exists,
    hex_to_bin,
    join_path,
    join_path_native,
    to_native_path_linux,
)

# typing ------------------------------------------------------------------

from typing import (
    Any,
    Iterator,
    List,
    TYPE_CHECKING,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

from git.types import AnyGitObject, PathLike

if TYPE_CHECKING:
    from git.config import GitConfigParser
    from git.objects.commit import Actor
    from git.refs.log import RefLogEntry
    from git.refs.reference import Reference
    from git.repo import Repo


T_References = TypeVar("T_References", bound="SymbolicReference")

# ------------------------------------------------------------------------------


def _git_dir(repo: "Repo", path: Union[PathLike, None]) -> PathLike:
    """Find the git dir that is appropriate for the path."""
    name = f"{path}"
    if name in ["HEAD", "ORIG_HEAD", "FETCH_HEAD", "index", "logs"]:
        return repo.git_dir
    return repo.common_dir


class SymbolicReference:
    """Special case of a reference that is symbolic.

    This does not point to a specific commit, but to another
    :class:`~git.refs.head.Head`, which itself specifies a commit.

    A typical example for a symbolic reference is :class:`~git.refs.head.HEAD`.
    """

    __slots__ = ("repo", "path")

    _resolve_ref_on_create = False
    _points_to_commits_only = True
    _common_path_default = ""
    _remote_common_path_default = "refs/remotes"
    _id_attribute_ = "name"

    def __init__(self, repo: "Repo", path: PathLike, check_path: bool = False) -> None:
        self.repo = repo
        self.path: PathLike = path

    def __str__(self) -> str:
        return os.fspath(self.path)

    def __repr__(self) -> str:
        return '<git.%s "%s">' % (self.__class__.__name__, self.path)

    def __eq__(self, other: object) -> bool:
        if hasattr(other, "path"):
            other = cast(SymbolicReference, other)
            return self.path == other.path
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __hash__(self) -> int:
        return hash(self.path)

    @property
    def name(self) -> str:
        """
        :return:
            In case of symbolic references, the shortest assumable name is the path
            itself.
        """
        return os.fspath(self.path)

    @property
    def abspath(self) -> PathLike:
        return join_path_native(_git_dir(self.repo, self.path), self.path)

    @classmethod
    def _get_packed_refs_path(cls, repo: "Repo") -> str:
        return os.path.join(repo.common_dir, "packed-refs")

    @classmethod
    def _iter_packed_refs(cls, repo: "Repo") -> Iterator[Tuple[str, str]]:
        """Return an iterator yielding pairs of sha1/path pairs (as strings) for the
        corresponding refs.

        :note:
            The packed refs file will be kept open as long as we iterate.
        """
        try:
            with open(cls._get_packed_refs_path(repo), "rt", encoding="UTF-8") as fp:
                for line in fp:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith("#"):
                        # "# pack-refs with: peeled fully-peeled sorted"
                        # the git source code shows "peeled",
                        # "fully-peeled" and "sorted" as the keywords
                        # that can go on this line, as per comments in git file
                        # refs/packed-backend.c
                        # I looked at master on 2017-10-11,
                        # commit 111ef79afe, after tag v2.15.0-rc1
                        # from repo https://github.com/git/git.git
                        if line.startswith("# pack-refs with:") and "peeled" not in line:
                            raise TypeError("PackingType of packed-Refs not understood: %r" % line)
                        # END abort if we do not understand the packing scheme
                        continue
                    # END parse comment

                    # Skip dereferenced tag object entries - previous line was actual
                    # tag reference for it.
                    if line[0] == "^":
                        continue

                    yield cast(Tuple[str, str], tuple(line.split(" ", 1)))
                # END for each line
        except OSError:
            return None
        # END no packed-refs file handling

    @classmethod
    def dereference_recursive(cls, repo: "Repo", ref_path: Union[PathLike, None]) -> str:
        """
        :return:
            hexsha stored in the reference at the given `ref_path`, recursively
            dereferencing all intermediate references as required

        :param repo:
            The repository containing the reference at `ref_path`.
        """

        while True:
            hexsha, ref_path = cls._get_ref_info(repo, ref_path)
            if hexsha is not None:
                return hexsha
        # END recursive dereferencing

    @staticmethod
    def _check_ref_name_valid(ref_path: PathLike) -> None:
        """Check a ref name for validity.

        This is based on the rules described in :manpage:`git-check-ref-format(1)`.
        """
        previous: Union[str, None] = None
        one_before_previous: Union[str, None] = None
        for c in os.fspath(ref_path):
            if c in " ~^:?*[\\":
                raise ValueError(
                    f"Invalid reference '{ref_path}': references cannot contain spaces, tildes (~), carets (^),"
                    f" colons (:), question marks (?), asterisks (*), open brackets ([) or backslashes (\\)"
                )
            elif c == ".":
                if previous is None or previous == "/":
                    raise ValueError(
                        f"Invalid reference '{ref_path}': references cannot start with a period (.) or contain '/.'"
                    )
                elif previous == ".":
                    raise ValueError(f"Invalid reference '{ref_path}': references cannot contain '..'")
            elif c == "/":
                if previous == "/":
                    raise ValueError(f"Invalid reference '{ref_path}': references cannot contain '//'")
                elif previous is None:
                    raise ValueError(
                        f"Invalid reference '{ref_path}': references cannot start with forward slashes '/'"
                    )
            elif c == "{" and previous == "@":
                raise ValueError(f"Invalid reference '{ref_path}': references cannot contain '@{{'")
            elif ord(c) < 32 or ord(c) == 127:
                raise ValueError(f"Invalid reference '{ref_path}': references cannot contain ASCII control characters")

            one_before_previous = previous
            previous = c

        if previous == ".":
            raise ValueError(f"Invalid reference '{ref_path}': references cannot end with a period (.)")
        elif previous == "/":
            raise ValueError(f"Invalid reference '{ref_path}': references cannot end with a forward slash (/)")
        elif previous == "@" and one_before_previous is None:
            raise ValueError(f"Invalid reference '{ref_path}': references cannot be '@'")
        elif any(component.endswith(".lock") for component in Path(ref_path).parts):
            raise ValueError(
                f"Invalid reference '{ref_path}': references cannot have slash-separated components that end with"
                " '.lock'"
            )

    @classmethod
    def _get_ref_info_helper(
        cls, repo: "Repo", ref_path: Union[PathLike, None]
    ) -> Union[Tuple[str, None], Tuple[None, str]]:
        """
        :return:
            *(str(sha), str(target_ref_path))*, where:

            * *sha* is of the file at rela_path points to if available, or ``None``.
            * *target_ref_path* is the reference we point to, or ``None``.
        """
        if ref_path:
            cls._check_ref_name_valid(ref_path)

        tokens: Union[None, List[str], Tuple[str, str]] = None
        repodir = _git_dir(repo, ref_path)
        try:
            with open(os.path.join(repodir, ref_path), "rt", encoding="UTF-8") as fp:  # type: ignore[arg-type]
                value = fp.read().rstrip()
            # Don't only split on spaces, but on whitespace, which allows to parse lines like:
            # 60b64ef992065e2600bfef6187a97f92398a9144                branch 'master' of git-server:/path/to/repo
            tokens = value.split()
            assert len(tokens) != 0
        except OSError:
            # Probably we are just packed. Find our entry in the packed refs file.
            # NOTE: We are not a symbolic ref if we are in a packed file, as these
            # are excluded explicitly.
            for sha, path in cls._iter_packed_refs(repo):
                if path != ref_path:
                    continue
                # sha will be used.
                tokens = sha, path
                break
            # END for each packed ref
        # END handle packed refs
        if tokens is None:
            raise ValueError("Reference at %r does not exist" % ref_path)

        # Is it a reference?
        if tokens[0] == "ref:":
            return (None, tokens[1])

        # It's a commit.
        if repo.re_hexsha_only.match(tokens[0]):
            return (tokens[0], None)

        raise ValueError("Failed to parse reference information from %r" % ref_path)

    @classmethod
    def _get_ref_info(cls, repo: "Repo", ref_path: Union[PathLike, None]) -> Union[Tuple[str, None], Tuple[None, str]]:
        """
        :return:
            *(str(sha), str(target_ref_path))*, where:

            * *sha* is of the file at rela_path points to if available, or ``None``.
            * *target_ref_path* is the reference we point to, or ``None``.
        """
        return cls._get_ref_info_helper(repo, ref_path)

    def _get_object(self) -> AnyGitObject:
        """
        :return:
            The object our ref currently refers to. Refs can be cached, they will always
            point to the actual object as it gets re-created on each query.
        """
        # We have to be dynamic here as we may be a tag which can point to anything.
        # Our path will be resolved to the hexsha which will be used accordingly.
        return Object.new_from_sha(self.repo, hex_to_bin(self.dereference_recursive(self.repo, self.path)))

    def _get_commit(self) -> "Commit":
        """
        :return:
            :class:`~git.objects.commit.Commit` object we point to. This works for
            detached and non-detached :class:`SymbolicReference` instances. The symbolic
            reference will be dereferenced recursively.
        """
        obj = self._get_object()
        if obj.type == "tag":
            obj = obj.object
        # END dereference tag

        if obj.type != Commit.type:
            raise TypeError("Symbolic Reference pointed to object %r, commit was required" % obj)
        # END handle type
        return obj

    def set_commit(
        self,
        commit: Union[Commit, "SymbolicReference", str],
        logmsg: Union[str, None] = None,
    ) -> "SymbolicReference":
        """Like :meth:`set_object`, but restricts the type of object to be a
        :class:`~git.objects.commit.Commit`.

        :raise ValueError:
            If `commit` is not a :class:`~git.objects.commit.Commit` object, nor does it
            point to a commit.

        :return:
            self
        """
        # Check the type - assume the best if it is a base-string.
        invalid_type = False
        if isinstance(commit, Object):
            invalid_type = commit.type != Commit.type
        elif isinstance(commit, SymbolicReference):
            invalid_type = commit.object.type != Commit.type
        else:
            try:
                invalid_type = self.repo.rev_parse(commit).type != Commit.type
            except (BadObject, BadName) as e:
                raise ValueError("Invalid object: %s" % commit) from e
            # END handle exception
        # END verify type

        if invalid_type:
            raise ValueError("Need commit, got %r" % commit)
        # END handle raise

        # We leave strings to the rev-parse method below.
        self.set_object(commit, logmsg)

        return self

    def set_object(
        self,
        object: Union[AnyGitObject, "SymbolicReference", str],
        logmsg: Union[str, None] = None,
    ) -> "SymbolicReference":
        """Set the object we point to, possibly dereference our symbolic reference
        first. If the reference does not exist, it will be created.

        :param object:
            A refspec, a :class:`SymbolicReference` or an
            :class:`~git.objects.base.Object` instance.

            * :class:`SymbolicReference` instances will be dereferenced beforehand to
              obtain the git object they point to.
            * :class:`~git.objects.base.Object` instances must represent git objects
              (:class:`~git.types.AnyGitObject`).

        :param logmsg:
            If not ``None``, the message will be used in the reflog entry to be written.
            Otherwise the reflog is not altered.

        :note:
            Plain :class:`SymbolicReference` instances may not actually point to objects
            by convention.

        :return:
            self
        """
        if isinstance(object, SymbolicReference):
            object = object.object  # @ReservedAssignment
        # END resolve references

        is_detached = True
        try:
            is_detached = self.is_detached
        except ValueError:
            pass
        # END handle non-existing ones

        if is_detached:
            return self.set_reference(object, logmsg)

        # set the commit on our reference
        return self._get_reference().set_object(object, logmsg)

    @property
    def commit(self) -> "Commit":
        """Query or set commits directly"""
        return self._get_commit()

    @commit.setter
    def commit(self, commit: Union[Commit, "SymbolicReference", str]) -> "SymbolicReference":
        return self.set_commit(commit)

    @property
    def object(self) -> AnyGitObject:
        """Return the object our ref currently refers to"""
        return self._get_object()

    @object.setter
    def object(self, object: Union[AnyGitObject, "SymbolicReference", str]) -> "SymbolicReference":
        return self.set_object(object)

    def _get_reference(self) -> "Reference":
        """
        :return:
            :class:`~git.refs.reference.Reference` object we point to

        :raise TypeError:
            If this symbolic reference is detached, hence it doesn't point to a
            reference, but to a commit.
        """
        sha, target_ref_path = self._get_ref_info(self.repo, self.path)
        if target_ref_path is None:
            raise TypeError("%s is a detached symbolic reference as it points to %r" % (self, sha))
        return cast("Reference", self.from_path(self.repo, target_ref_path))

    def set_reference(
        self,
        ref: Union[AnyGitObject, "SymbolicReference", str],
        logmsg: Union[str, None] = None,
    ) -> "SymbolicReference":
        """Set ourselves to the given `ref`.

        It will stay a symbol if the `ref` is a :class:`~git.refs.reference.Reference`.

        Otherwise a git object, specified as a :class:`~git.objects.base.Object`
        instance or refspec, is assumed. If it is valid, this reference will be set to
        it, which effectively detaches the reference if it was a purely symbolic one.

        :param ref:
            A :class:`SymbolicReference` instance, an :class:`~git.objects.base.Object`
            instance (specifically an :class:`~git.types.AnyGitObject`), or a refspec
            string. Only if the ref is a :class:`SymbolicReference` instance, we will
            point to it. Everything else is dereferenced to obtain the actual object.

        :param logmsg:
            If set to a string, the message will be used in the reflog.
            Otherwise, a reflog entry is not written for the changed reference.
            The previous commit of the entry will be the commit we point to now.

            See also: :meth:`log_append`

        :return:
            self

        :note:
            This symbolic reference will not be dereferenced. For that, see
            :meth:`set_object`.
        """
        write_value = None
        obj = None
        if isinstance(ref, SymbolicReference):
            write_value = "ref: %s" % ref.path
        elif isinstance(ref, Object):
            obj = ref
            write_value = ref.hexsha
        elif isinstance(ref, str):
            try:
                obj = self.repo.rev_parse(ref + "^{}")  # Optionally dereference tags.
                write_value = obj.hexsha
            except (BadObject, BadName) as e:
                raise ValueError("Could not extract object from %s" % ref) from e
            # END end try string
        else:
            raise ValueError("Unrecognized Value: %r" % ref)
        # END try commit attribute

        # typecheck
        if obj is not None and self._points_to_commits_only and obj.type != Commit.type:
            raise TypeError("Require commit, got %r" % obj)
        # END verify type

        oldbinsha: bytes = b""
        if logmsg is not None:
            try:
                oldbinsha = self.commit.binsha
            except ValueError:
                oldbinsha = Commit.NULL_BIN_SHA
            # END handle non-existing
        # END retrieve old hexsha

        fpath = self.abspath
        assure_directory_exists(fpath, is_file=True)

        lfd = LockedFD(fpath)
        fd = lfd.open(write=True, stream=True)
        try:
            fd.write(write_value.encode("utf-8") + b"\n")
            lfd.commit()
        except BaseException:
            lfd.rollback()
            raise
        # Adjust the reflog
        if logmsg is not None:
            self.log_append(oldbinsha, logmsg)

        return self

    # Aliased reference
    @property
    def reference(self) -> "Reference":
        return self._get_reference()

    @reference.setter
    def reference(self, ref: Union[AnyGitObject, "SymbolicReference", str]) -> "SymbolicReference":
        return self.set_reference(ref)

    ref = reference

    def is_valid(self) -> bool:
        """
        :return:
            ``True`` if the reference is valid, hence it can be read and points to a
            valid object or reference.
        """
        try:
            self.object  # noqa: B018
        except (OSError, ValueError):
            return False
        else:
            return True

    @property
    def is_detached(self) -> bool:
        """
        :return:
            ``True`` if we are a detached reference, hence we point to a specific commit
            instead to another reference.
        """
        try:
            self.ref  # noqa: B018
            return False
        except TypeError:
            return True

    def log(self) -> "RefLog":
        """
        :return:
            :class:`~git.refs.log.RefLog` for this reference.
            Its last entry reflects the latest change applied to this reference.

        :note:
            As the log is parsed every time, its recommended to cache it for use instead
            of calling this method repeatedly. It should be considered read-only.
        """
        return RefLog.from_file(RefLog.path(self))

    def log_append(
        self,
        oldbinsha: bytes,
        message: Union[str, None],
        newbinsha: Union[bytes, None] = None,
    ) -> "RefLogEntry":
        """Append a logentry to the logfile of this ref.

        :param oldbinsha:
            Binary sha this ref used to point to.

        :param message:
            A message describing the change.

        :param newbinsha:
            The sha the ref points to now. If None, our current commit sha will be used.

        :return:
            The added :class:`~git.refs.log.RefLogEntry` instance.
        """
        # NOTE: We use the committer of the currently active commit - this should be
        # correct to allow overriding the committer on a per-commit level.
        # See https://github.com/gitpython-developers/GitPython/pull/146.
        try:
            committer_or_reader: Union["Actor", "GitConfigParser"] = self.commit.committer
        except ValueError:
            committer_or_reader = self.repo.config_reader()
        # END handle newly cloned repositories
        if newbinsha is None:
            newbinsha = self.commit.binsha

        if message is None:
            message = ""

        return RefLog.append_entry(committer_or_reader, RefLog.path(self), oldbinsha, newbinsha, message)

    def log_entry(self, index: int) -> "RefLogEntry":
        """
        :return:
            :class:`~git.refs.log.RefLogEntry` at the given index

        :param index:
            Python list compatible positive or negative index.

        :note:
            This method must read part of the reflog during execution, hence it should
            be used sparingly, or only if you need just one index. In that case, it will
            be faster than the :meth:`log` method.
        """
        return RefLog.entry_at(RefLog.path(self), index)

    @classmethod
    def to_full_path(cls, path: Union[PathLike, "SymbolicReference"]) -> PathLike:
        """
        :return:
            String with a full repository-relative path which can be used to initialize
            a :class:`~git.refs.reference.Reference` instance, for instance by using
            :meth:`Reference.from_path <git.refs.reference.Reference.from_path>`.
        """
        if isinstance(path, SymbolicReference):
            path = path.path
        full_ref_path = path
        if not cls._common_path_default:
            return full_ref_path
        if not os.fspath(path).startswith(cls._common_path_default + "/"):
            full_ref_path = "%s/%s" % (cls._common_path_default, path)
        return full_ref_path

    @classmethod
    def delete(cls, repo: "Repo", path: PathLike) -> None:
        """Delete the reference at the given path.

        :param repo:
            Repository to delete the reference from.

        :param path:
            Short or full path pointing to the reference, e.g. ``refs/myreference`` or
            just ``myreference``, hence ``refs/`` is implied.
            Alternatively the symbolic reference to be deleted.
        """
        full_ref_path = cls.to_full_path(path)
        abs_path = os.path.join(repo.common_dir, full_ref_path)
        if os.path.exists(abs_path):
            os.remove(abs_path)
        else:
            # Check packed refs.
            pack_file_path = cls._get_packed_refs_path(repo)
            try:
                with open(pack_file_path, "rb") as reader:
                    new_lines = []
                    made_change = False
                    dropped_last_line = False
                    for line_bytes in reader:
                        line = line_bytes.decode(defenc)
                        _, _, line_ref = line.partition(" ")
                        line_ref = line_ref.strip()
                        # Keep line if it is a comment or if the ref to delete is not in
                        # the line.
                        # If we deleted the last line and this one is a tag-reference
                        # object, we drop it as well.
                        if (line.startswith("#") or full_ref_path != line_ref) and (
                            not dropped_last_line or dropped_last_line and not line.startswith("^")
                        ):
                            new_lines.append(line)
                            dropped_last_line = False
                            continue
                        # END skip comments and lines without our path

                        # Drop this line.
                        made_change = True
                        dropped_last_line = True

                # Write the new lines.
                if made_change:
                    # Binary writing is required, otherwise Windows will open the file
                    # in text mode and change LF to CRLF!
                    with open(pack_file_path, "wb") as fd:
                        fd.writelines(line.encode(defenc) for line in new_lines)

            except OSError:
                pass  # It didn't exist at all.

        # Delete the reflog.
        reflog_path = RefLog.path(cls(repo, full_ref_path))
        if os.path.isfile(reflog_path):
            os.remove(reflog_path)
        # END remove reflog

    @classmethod
    def _create(
        cls: Type[T_References],
        repo: "Repo",
        path: PathLike,
        resolve: bool,
        reference: Union["SymbolicReference", str],
        force: bool,
        logmsg: Union[str, None] = None,
    ) -> T_References:
        """Internal method used to create a new symbolic reference.

        If `resolve` is ``False``, the reference will be taken as is, creating a proper
        symbolic reference. Otherwise it will be resolved to the corresponding object
        and a detached symbolic reference will be created instead.
        """
        git_dir = _git_dir(repo, path)
        full_ref_path = cls.to_full_path(path)
        abs_ref_path = os.path.join(git_dir, full_ref_path)

        # Figure out target data.
        target = reference
        if resolve:
            target = repo.rev_parse(str(reference))

        if not force and os.path.isfile(abs_ref_path):
            target_data = str(target)
            if isinstance(target, SymbolicReference):
                target_data = os.fspath(target.path)
            if not resolve:
                target_data = "ref: " + target_data
            with open(abs_ref_path, "rb") as fd:
                existing_data = fd.read().decode(defenc).strip()
            if existing_data != target_data:
                raise OSError(
                    "Reference at %r does already exist, pointing to %r, requested was %r"
                    % (full_ref_path, existing_data, target_data)
                )
        # END no force handling

        ref = cls(repo, full_ref_path)
        ref.set_reference(target, logmsg)
        return ref

    @classmethod
    def create(
        cls: Type[T_References],
        repo: "Repo",
        path: PathLike,
        reference: Union["SymbolicReference", str] = "HEAD",
        logmsg: Union[str, None] = None,
        force: bool = False,
        **kwargs: Any,
    ) -> T_References:
        """Create a new symbolic reference: a reference pointing to another reference.

        :param repo:
            Repository to create the reference in.

        :param path:
            Full path at which the new symbolic reference is supposed to be created at,
            e.g. ``NEW_HEAD`` or ``symrefs/my_new_symref``.

        :param reference:
            The reference which the new symbolic reference should point to.
            If it is a commit-ish, the symbolic ref will be detached.

        :param force:
            If ``True``, force creation even if a symbolic reference with that name
            already exists. Raise :exc:`OSError` otherwise.

        :param logmsg:
            If not ``None``, the message to append to the reflog.
            If ``None``, no reflog entry is written.

        :return:
            Newly created symbolic reference

        :raise OSError:
            If a (Symbolic)Reference with the same name but different contents already
            exists.

        :note:
            This does not alter the current HEAD, index or working tree.
        """
        return cls._create(repo, path, cls._resolve_ref_on_create, reference, force, logmsg)

    def rename(self, new_path: PathLike, force: bool = False) -> "SymbolicReference":
        """Rename self to a new path.

        :param new_path:
            Either a simple name or a full path, e.g. ``new_name`` or
            ``features/new_name``.
            The prefix ``refs/`` is implied for references and will be set as needed.
            In case this is a symbolic ref, there is no implied prefix.

        :param force:
            If ``True``, the rename will succeed even if a head with the target name
            already exists. It will be overwritten in that case.

        :return:
            self

        :raise OSError:
            If a file at path but with different contents already exists.
        """
        new_path = self.to_full_path(new_path)
        if self.path == new_path:
            return self

        new_abs_path = os.path.join(_git_dir(self.repo, new_path), new_path)
        cur_abs_path = os.path.join(_git_dir(self.repo, self.path), self.path)
        if os.path.isfile(new_abs_path):
            if not force:
                # If they point to the same file, it's not an error.
                with open(new_abs_path, "rb") as fd1:
                    f1 = fd1.read().strip()
                with open(cur_abs_path, "rb") as fd2:
                    f2 = fd2.read().strip()
                if f1 != f2:
                    raise OSError("File at path %r already exists" % new_abs_path)
                # else: We could remove ourselves and use the other one, but...
                # ...for clarity, we just continue as usual.
            # END not force handling
            os.remove(new_abs_path)
        # END handle existing target file

        dname = os.path.dirname(new_abs_path)
        if not os.path.isdir(dname):
            os.makedirs(dname)
        # END create directory

        os.rename(cur_abs_path, new_abs_path)
        self.path = new_path

        return self

    @classmethod
    def _iter_items(
        cls: Type[T_References], repo: "Repo", common_path: Union[PathLike, None] = None
    ) -> Iterator[T_References]:
        if common_path is None:
            common_path = cls._common_path_default
        rela_paths = set()

        # Walk loose refs.
        # Currently we do not follow links.
        for root, dirs, files in os.walk(join_path_native(repo.common_dir, common_path)):
            if "refs" not in root.split(os.sep):  # Skip non-refs subfolders.
                refs_id = [d for d in dirs if d == "refs"]
                if refs_id:
                    dirs[0:] = ["refs"]
            # END prune non-refs folders

            for f in files:
                if f == "packed-refs":
                    continue
                abs_path = to_native_path_linux(join_path(root, f))
                rela_paths.add(abs_path.replace(to_native_path_linux(repo.common_dir) + "/", ""))
            # END for each file in root directory
        # END for each directory to walk

        # Read packed refs.
        for _sha, rela_path in cls._iter_packed_refs(repo):
            if rela_path.startswith(os.fspath(common_path)):
                rela_paths.add(rela_path)
            # END relative path matches common path
        # END packed refs reading

        # Yield paths in sorted order.
        for path in sorted(rela_paths):
            try:
                yield cls.from_path(repo, path)
            except ValueError:
                continue
        # END for each sorted relative refpath

    @classmethod
    def iter_items(
        cls: Type[T_References],
        repo: "Repo",
        common_path: Union[PathLike, None] = None,
        *args: Any,
        **kwargs: Any,
    ) -> Iterator[T_References]:
        """Find all refs in the repository.

        :param repo:
            The :class:`~git.repo.base.Repo`.

        :param common_path:
            Optional keyword argument to the path which is to be shared by all returned
            Ref objects.
            Defaults to class specific portion if ``None``, ensuring that only refs
            suitable for the actual class are returned.

        :return:
            A list of :class:`SymbolicReference`, each guaranteed to be a symbolic ref
            which is not detached and pointing to a valid ref.

            The list is lexicographically sorted. The returned objects are instances of
            concrete subclasses, such as :class:`~git.refs.head.Head` or
            :class:`~git.refs.tag.TagReference`.
        """
        return (r for r in cls._iter_items(repo, common_path) if r.__class__ is SymbolicReference or not r.is_detached)

    @classmethod
    def from_path(cls: Type[T_References], repo: "Repo", path: PathLike) -> T_References:
        """Make a symbolic reference from a path.

        :param path:
            Full ``.git``-directory-relative path name to the Reference to instantiate.

        :note:
            Use :meth:`to_full_path` if you only have a partial path of a known
            Reference type.

        :return:
            Instance of type :class:`~git.refs.reference.Reference`,
            :class:`~git.refs.head.Head`, or :class:`~git.refs.tag.Tag`, depending on
            the given path.
        """
        if not path:
            raise ValueError("Cannot create Reference from %r" % path)

        # Names like HEAD are inserted after the refs module is imported - we have an
        # import dependency cycle and don't want to import these names in-function.
        from . import HEAD, Head, RemoteReference, TagReference, Reference

        for ref_type in (
            HEAD,
            Head,
            RemoteReference,
            TagReference,
            Reference,
            SymbolicReference,
        ):
            try:
                instance = cast(T_References, ref_type(repo, path))
                if instance.__class__ is SymbolicReference and instance.is_detached:
                    raise ValueError("SymbolicRef was detached, we drop it")
                else:
                    return instance

            except ValueError:
                pass
            # END exception handling
        # END for each type to try
        raise ValueError("Could not find reference type suitable to handle path %r" % path)

    def is_remote(self) -> bool:
        """:return: True if this symbolic reference points to a remote branch"""
        return os.fspath(self.path).startswith(self._remote_common_path_default + "/")
