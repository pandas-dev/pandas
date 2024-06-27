# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

__all__ = ["DiffConstants", "NULL_TREE", "INDEX", "Diffable", "DiffIndex", "Diff"]

import enum
import re
import warnings

from git.cmd import handle_process_output
from git.compat import defenc
from git.objects.blob import Blob
from git.objects.util import mode_str_to_int
from git.util import finalize_process, hex_to_bin

# typing ------------------------------------------------------------------

from typing import (
    Any,
    Iterator,
    List,
    Match,
    Optional,
    Tuple,
    TYPE_CHECKING,
    TypeVar,
    Union,
    cast,
)
from git.types import Literal, PathLike

if TYPE_CHECKING:
    from subprocess import Popen

    from git.cmd import Git
    from git.objects.base import IndexObject
    from git.objects.commit import Commit
    from git.objects.tree import Tree
    from git.repo.base import Repo

Lit_change_type = Literal["A", "D", "C", "M", "R", "T", "U"]

# ------------------------------------------------------------------------


@enum.unique
class DiffConstants(enum.Enum):
    """Special objects for :meth:`Diffable.diff`.

    See the :meth:`Diffable.diff` method's ``other`` parameter, which accepts various
    values including these.

    :note:
        These constants are also available as attributes of the :mod:`git.diff` module,
        the :class:`Diffable` class and its subclasses and instances, and the top-level
        :mod:`git` module.
    """

    NULL_TREE = enum.auto()
    """Stand-in indicating you want to compare against the empty tree in diffs.

    Also accessible as :const:`git.NULL_TREE`, :const:`git.diff.NULL_TREE`, and
    :const:`Diffable.NULL_TREE`.
    """

    INDEX = enum.auto()
    """Stand-in indicating you want to diff against the index.

    Also accessible as :const:`git.INDEX`, :const:`git.diff.INDEX`, and
    :const:`Diffable.INDEX`, as well as :const:`Diffable.Index`. The latter has been
    kept for backward compatibility and made an alias of this, so it may still be used.
    """


NULL_TREE: Literal[DiffConstants.NULL_TREE] = DiffConstants.NULL_TREE
"""Stand-in indicating you want to compare against the empty tree in diffs.

See :meth:`Diffable.diff`, which accepts this as a value of its ``other`` parameter.

This is an alias of :const:`DiffConstants.NULL_TREE`, which may also be accessed as
:const:`git.NULL_TREE` and :const:`Diffable.NULL_TREE`.
"""

INDEX: Literal[DiffConstants.INDEX] = DiffConstants.INDEX
"""Stand-in indicating you want to diff against the index.

See :meth:`Diffable.diff`, which accepts this as a value of its ``other`` parameter.

This is an alias of :const:`DiffConstants.INDEX`, which may also be accessed as
:const:`git.INDEX` and :const:`Diffable.INDEX`, as well as :const:`Diffable.Index`.
"""

_octal_byte_re = re.compile(rb"\\([0-9]{3})")


def _octal_repl(matchobj: Match) -> bytes:
    value = matchobj.group(1)
    value = int(value, 8)
    value = bytes(bytearray((value,)))
    return value


def decode_path(path: bytes, has_ab_prefix: bool = True) -> Optional[bytes]:
    if path == b"/dev/null":
        return None

    if path.startswith(b'"') and path.endswith(b'"'):
        path = path[1:-1].replace(b"\\n", b"\n").replace(b"\\t", b"\t").replace(b'\\"', b'"').replace(b"\\\\", b"\\")

    path = _octal_byte_re.sub(_octal_repl, path)

    if has_ab_prefix:
        assert path.startswith(b"a/") or path.startswith(b"b/")
        path = path[2:]

    return path


class Diffable:
    """Common interface for all objects that can be diffed against another object of
    compatible type.

    :note:
        Subclasses require a :attr:`repo` member, as it is the case for
        :class:`~git.objects.base.Object` instances. For practical reasons we do not
        derive from :class:`~git.objects.base.Object`.
    """

    __slots__ = ()

    repo: "Repo"
    """Repository to operate on. Must be provided by subclass or sibling class."""

    NULL_TREE = NULL_TREE
    """Stand-in indicating you want to compare against the empty tree in diffs.

    See the :meth:`diff` method, which accepts this as a value of its ``other``
    parameter.

    This is the same as :const:`DiffConstants.NULL_TREE`, and may also be accessed as
    :const:`git.NULL_TREE` and :const:`git.diff.NULL_TREE`.
    """

    INDEX = INDEX
    """Stand-in indicating you want to diff against the index.

    See the :meth:`diff` method, which accepts this as a value of its ``other``
    parameter.

    This is the same as :const:`DiffConstants.INDEX`, and may also be accessed as
    :const:`git.INDEX` and :const:`git.diff.INDEX`, as well as :class:`Diffable.INDEX`,
    which is kept for backward compatibility (it is now defined an alias of this).
    """

    Index = INDEX
    """Stand-in indicating you want to diff against the index
    (same as :const:`~Diffable.INDEX`).

    This is an alias of :const:`~Diffable.INDEX`, for backward compatibility. See
    :const:`~Diffable.INDEX` and :meth:`diff` for details.

    :note:
        Although always meant for use as an opaque constant, this was formerly defined
        as a class. Its usage is unchanged, but static type annotations that attempt
        to permit only this object must be changed to avoid new mypy errors. This was
        previously not possible to do, though ``Type[Diffable.Index]`` approximated it.
        It is now possible to do precisely, using ``Literal[DiffConstants.INDEX]``.
    """

    def _process_diff_args(
        self,
        args: List[Union[PathLike, "Diffable"]],
    ) -> List[Union[PathLike, "Diffable"]]:
        """
        :return:
            Possibly altered version of the given args list.
            This method is called right before git command execution.
            Subclasses can use it to alter the behaviour of the superclass.
        """
        return args

    def diff(
        self,
        other: Union[DiffConstants, "Tree", "Commit", str, None] = INDEX,
        paths: Union[PathLike, List[PathLike], Tuple[PathLike, ...], None] = None,
        create_patch: bool = False,
        **kwargs: Any,
    ) -> "DiffIndex":
        """Create diffs between two items being trees, trees and index or an index and
        the working tree. Detects renames automatically.

        :param other:
            This the item to compare us with.

            * If ``None``, we will be compared to the working tree.

            * If a :class:`~git.types.Tree_ish` or string, it will be compared against
              the respective tree.

            * If :const:`INDEX`, it will be compared against the index.

            * If :const:`NULL_TREE`, it will compare against the empty tree.

            This parameter defaults to :const:`INDEX` (rather than ``None``) so that the
            method will not by default fail on bare repositories.

        :param paths:
            This a list of paths or a single path to limit the diff to. It will only
            include at least one of the given path or paths.

        :param create_patch:
            If ``True``, the returned :class:`Diff` contains a detailed patch that if
            applied makes the self to other. Patches are somewhat costly as blobs have
            to be read and diffed.

        :param kwargs:
            Additional arguments passed to :manpage:`git-diff(1)`, such as ``R=True`` to
            swap both sides of the diff.

        :return:
            A :class:`DiffIndex` representing the computed diff.

        :note:
            On a bare repository, `other` needs to be provided as :const:`INDEX`, or as
            an instance of :class:`~git.objects.tree.Tree` or
            :class:`~git.objects.commit.Commit`, or a git command error will occur.
        """
        args: List[Union[PathLike, Diffable]] = []
        args.append("--abbrev=40")  # We need full shas.
        args.append("--full-index")  # Get full index paths, not only filenames.

        # Remove default '-M' arg (check for renames) if user is overriding it.
        if not any(x in kwargs for x in ("find_renames", "no_renames", "M")):
            args.append("-M")

        if create_patch:
            args.append("-p")
            args.append("--no-ext-diff")
        else:
            args.append("--raw")
            args.append("-z")

        # Ensure we never see colored output.
        # Fixes: https://github.com/gitpython-developers/GitPython/issues/172
        args.append("--no-color")

        if paths is not None and not isinstance(paths, (tuple, list)):
            paths = [paths]

        diff_cmd = self.repo.git.diff
        if other is INDEX:
            args.insert(0, "--cached")
        elif other is NULL_TREE:
            args.insert(0, "-r")  # Recursive diff-tree.
            args.insert(0, "--root")
            diff_cmd = self.repo.git.diff_tree
        elif other is not None:
            args.insert(0, "-r")  # Recursive diff-tree.
            args.insert(0, other)
            diff_cmd = self.repo.git.diff_tree

        args.insert(0, self)

        # paths is a list or tuple here, or None.
        if paths:
            args.append("--")
            args.extend(paths)
        # END paths handling

        kwargs["as_process"] = True
        proc = diff_cmd(*self._process_diff_args(args), **kwargs)

        diff_method = Diff._index_from_patch_format if create_patch else Diff._index_from_raw_format
        index = diff_method(self.repo, proc)

        proc.wait()
        return index


T_Diff = TypeVar("T_Diff", bound="Diff")


class DiffIndex(List[T_Diff]):
    R"""An index for diffs, allowing a list of :class:`Diff`\s to be queried by the diff
    properties.

    The class improves the diff handling convenience.
    """

    change_type = ("A", "C", "D", "R", "M", "T")
    """Change type invariant identifying possible ways a blob can have changed:

    * ``A`` = Added
    * ``D`` = Deleted
    * ``R`` = Renamed
    * ``M`` = Modified
    * ``T`` = Changed in the type
    """

    def iter_change_type(self, change_type: Lit_change_type) -> Iterator[T_Diff]:
        """
        :return:
            Iterator yielding :class:`Diff` instances that match the given `change_type`

        :param change_type:
            Member of :attr:`DiffIndex.change_type`, namely:

            * 'A' for added paths
            * 'D' for deleted paths
            * 'R' for renamed paths
            * 'M' for paths with modified data
            * 'T' for changed in the type paths
        """
        if change_type not in self.change_type:
            raise ValueError("Invalid change type: %s" % change_type)

        for diffidx in self:
            if diffidx.change_type == change_type:
                yield diffidx
            elif change_type == "A" and diffidx.new_file:
                yield diffidx
            elif change_type == "D" and diffidx.deleted_file:
                yield diffidx
            elif change_type == "C" and diffidx.copied_file:
                yield diffidx
            elif change_type == "R" and diffidx.renamed:
                yield diffidx
            elif change_type == "M" and diffidx.a_blob and diffidx.b_blob and diffidx.a_blob != diffidx.b_blob:
                yield diffidx
        # END for each diff


class Diff:
    """A Diff contains diff information between two Trees.

    It contains two sides a and b of the diff. Members are prefixed with "a" and "b"
    respectively to indicate that.

    Diffs keep information about the changed blob objects, the file mode, renames,
    deletions and new files.

    There are a few cases where ``None`` has to be expected as member variable value:

    New File::

        a_mode is None
        a_blob is None
        a_path is None

    Deleted File::

        b_mode is None
        b_blob is None
        b_path is None

    Working Tree Blobs:

        When comparing to working trees, the working tree blob will have a null hexsha
        as a corresponding object does not yet exist. The mode will be null as well. The
        path will be available, though.

        If it is listed in a diff, the working tree version of the file must differ from
        the version in the index or tree, and hence has been modified.
    """

    # Precompiled regex.
    re_header = re.compile(
        rb"""
                                ^diff[ ]--git
                                    [ ](?P<a_path_fallback>"?[ab]/.+?"?)[ ](?P<b_path_fallback>"?[ab]/.+?"?)\n
                                (?:^old[ ]mode[ ](?P<old_mode>\d+)\n
                                   ^new[ ]mode[ ](?P<new_mode>\d+)(?:\n|$))?
                                (?:^similarity[ ]index[ ]\d+%\n
                                   ^rename[ ]from[ ](?P<rename_from>.*)\n
                                   ^rename[ ]to[ ](?P<rename_to>.*)(?:\n|$))?
                                (?:^new[ ]file[ ]mode[ ](?P<new_file_mode>.+)(?:\n|$))?
                                (?:^deleted[ ]file[ ]mode[ ](?P<deleted_file_mode>.+)(?:\n|$))?
                                (?:^similarity[ ]index[ ]\d+%\n
                                   ^copy[ ]from[ ].*\n
                                   ^copy[ ]to[ ](?P<copied_file_name>.*)(?:\n|$))?
                                (?:^index[ ](?P<a_blob_id>[0-9A-Fa-f]+)
                                    \.\.(?P<b_blob_id>[0-9A-Fa-f]+)[ ]?(?P<b_mode>.+)?(?:\n|$))?
                                (?:^---[ ](?P<a_path>[^\t\n\r\f\v]*)[\t\r\f\v]*(?:\n|$))?
                                (?:^\+\+\+[ ](?P<b_path>[^\t\n\r\f\v]*)[\t\r\f\v]*(?:\n|$))?
                            """,
        re.VERBOSE | re.MULTILINE,
    )

    # These can be used for comparisons.
    NULL_HEX_SHA = "0" * 40
    NULL_BIN_SHA = b"\0" * 20

    __slots__ = (
        "a_blob",
        "b_blob",
        "a_mode",
        "b_mode",
        "a_rawpath",
        "b_rawpath",
        "new_file",
        "deleted_file",
        "copied_file",
        "raw_rename_from",
        "raw_rename_to",
        "diff",
        "change_type",
        "score",
    )

    def __init__(
        self,
        repo: "Repo",
        a_rawpath: Optional[bytes],
        b_rawpath: Optional[bytes],
        a_blob_id: Union[str, bytes, None],
        b_blob_id: Union[str, bytes, None],
        a_mode: Union[bytes, str, None],
        b_mode: Union[bytes, str, None],
        new_file: bool,
        deleted_file: bool,
        copied_file: bool,
        raw_rename_from: Optional[bytes],
        raw_rename_to: Optional[bytes],
        diff: Union[str, bytes, None],
        change_type: Optional[Lit_change_type],
        score: Optional[int],
    ) -> None:
        assert a_rawpath is None or isinstance(a_rawpath, bytes)
        assert b_rawpath is None or isinstance(b_rawpath, bytes)
        self.a_rawpath = a_rawpath
        self.b_rawpath = b_rawpath

        self.a_mode = mode_str_to_int(a_mode) if a_mode else None
        self.b_mode = mode_str_to_int(b_mode) if b_mode else None

        # Determine whether this diff references a submodule. If it does then
        # we need to overwrite "repo" to the corresponding submodule's repo instead.
        if repo and a_rawpath:
            for submodule in repo.submodules:
                if submodule.path == a_rawpath.decode(defenc, "replace"):
                    if submodule.module_exists():
                        repo = submodule.module()
                    break

        self.a_blob: Union["IndexObject", None]
        if a_blob_id is None or a_blob_id == self.NULL_HEX_SHA:
            self.a_blob = None
        else:
            self.a_blob = Blob(repo, hex_to_bin(a_blob_id), mode=self.a_mode, path=self.a_path)

        self.b_blob: Union["IndexObject", None]
        if b_blob_id is None or b_blob_id == self.NULL_HEX_SHA:
            self.b_blob = None
        else:
            self.b_blob = Blob(repo, hex_to_bin(b_blob_id), mode=self.b_mode, path=self.b_path)

        self.new_file: bool = new_file
        self.deleted_file: bool = deleted_file
        self.copied_file: bool = copied_file

        # Be clear and use None instead of empty strings.
        assert raw_rename_from is None or isinstance(raw_rename_from, bytes)
        assert raw_rename_to is None or isinstance(raw_rename_to, bytes)
        self.raw_rename_from = raw_rename_from or None
        self.raw_rename_to = raw_rename_to or None

        self.diff = diff
        self.change_type: Union[Lit_change_type, None] = change_type
        self.score = score

    def __eq__(self, other: object) -> bool:
        for name in self.__slots__:
            if getattr(self, name) != getattr(other, name):
                return False
        # END for each name
        return True

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __hash__(self) -> int:
        return hash(tuple(getattr(self, n) for n in self.__slots__))

    def __str__(self) -> str:
        h = "%s"
        if self.a_blob:
            h %= self.a_blob.path
        elif self.b_blob:
            h %= self.b_blob.path

        msg = ""
        line = None
        line_length = 0
        for b, n in zip((self.a_blob, self.b_blob), ("lhs", "rhs")):
            if b:
                line = "\n%s: %o | %s" % (n, b.mode, b.hexsha)
            else:
                line = "\n%s: None" % n
            # END if blob is not None
            line_length = max(len(line), line_length)
            msg += line
        # END for each blob

        # Add headline.
        h += "\n" + "=" * line_length

        if self.deleted_file:
            msg += "\nfile deleted in rhs"
        if self.new_file:
            msg += "\nfile added in rhs"
        if self.copied_file:
            msg += "\nfile %r copied from %r" % (self.b_path, self.a_path)
        if self.rename_from:
            msg += "\nfile renamed from %r" % self.rename_from
        if self.rename_to:
            msg += "\nfile renamed to %r" % self.rename_to
        if self.diff:
            msg += "\n---"
            try:
                msg += self.diff.decode(defenc) if isinstance(self.diff, bytes) else self.diff
            except UnicodeDecodeError:
                msg += "OMITTED BINARY DATA"
            # END handle encoding
            msg += "\n---"
        # END diff info

        return h + msg

    @property
    def a_path(self) -> Optional[str]:
        return self.a_rawpath.decode(defenc, "replace") if self.a_rawpath else None

    @property
    def b_path(self) -> Optional[str]:
        return self.b_rawpath.decode(defenc, "replace") if self.b_rawpath else None

    @property
    def rename_from(self) -> Optional[str]:
        return self.raw_rename_from.decode(defenc, "replace") if self.raw_rename_from else None

    @property
    def rename_to(self) -> Optional[str]:
        return self.raw_rename_to.decode(defenc, "replace") if self.raw_rename_to else None

    @property
    def renamed(self) -> bool:
        """Deprecated, use :attr:`renamed_file` instead.

        :return:
            ``True`` if the blob of our diff has been renamed

        :note:
            This property is deprecated.
            Please use the :attr:`renamed_file` property instead.
        """
        warnings.warn(
            "Diff.renamed is deprecated, use Diff.renamed_file instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.renamed_file

    @property
    def renamed_file(self) -> bool:
        """:return: ``True`` if the blob of our diff has been renamed"""
        return self.rename_from != self.rename_to

    @classmethod
    def _pick_best_path(cls, path_match: bytes, rename_match: bytes, path_fallback_match: bytes) -> Optional[bytes]:
        if path_match:
            return decode_path(path_match)

        if rename_match:
            return decode_path(rename_match, has_ab_prefix=False)

        if path_fallback_match:
            return decode_path(path_fallback_match)

        return None

    @classmethod
    def _index_from_patch_format(cls, repo: "Repo", proc: Union["Popen", "Git.AutoInterrupt"]) -> DiffIndex:
        """Create a new :class:`DiffIndex` from the given process output which must be
        in patch format.

        :param repo:
            The repository we are operating on.

        :param proc:
            :manpage:`git-diff(1)` process to read from
            (supports :class:`Git.AutoInterrupt <git.cmd.Git.AutoInterrupt>` wrapper).

        :return:
            :class:`DiffIndex`
        """

        # FIXME: Here SLURPING raw, need to re-phrase header-regexes linewise.
        text_list: List[bytes] = []
        handle_process_output(proc, text_list.append, None, finalize_process, decode_streams=False)

        # For now, we have to bake the stream.
        text = b"".join(text_list)
        index: "DiffIndex" = DiffIndex()
        previous_header: Union[Match[bytes], None] = None
        header: Union[Match[bytes], None] = None
        a_path, b_path = None, None  # For mypy.
        a_mode, b_mode = None, None  # For mypy.
        for _header in cls.re_header.finditer(text):
            (
                a_path_fallback,
                b_path_fallback,
                old_mode,
                new_mode,
                rename_from,
                rename_to,
                new_file_mode,
                deleted_file_mode,
                copied_file_name,
                a_blob_id,
                b_blob_id,
                b_mode,
                a_path,
                b_path,
            ) = _header.groups()

            new_file, deleted_file, copied_file = (
                bool(new_file_mode),
                bool(deleted_file_mode),
                bool(copied_file_name),
            )

            a_path = cls._pick_best_path(a_path, rename_from, a_path_fallback)
            b_path = cls._pick_best_path(b_path, rename_to, b_path_fallback)

            # Our only means to find the actual text is to see what has not been matched
            # by our regex, and then retro-actively assign it to our index.
            if previous_header is not None:
                index[-1].diff = text[previous_header.end() : _header.start()]
            # END assign actual diff

            # Make sure the mode is set if the path is set. Otherwise the resulting blob
            # is invalid. We just use the one mode we should have parsed.
            a_mode = old_mode or deleted_file_mode or (a_path and (b_mode or new_mode or new_file_mode))
            b_mode = b_mode or new_mode or new_file_mode or (b_path and a_mode)
            index.append(
                Diff(
                    repo,
                    a_path,
                    b_path,
                    a_blob_id and a_blob_id.decode(defenc),
                    b_blob_id and b_blob_id.decode(defenc),
                    a_mode and a_mode.decode(defenc),
                    b_mode and b_mode.decode(defenc),
                    new_file,
                    deleted_file,
                    copied_file,
                    rename_from,
                    rename_to,
                    None,
                    None,
                    None,
                )
            )

            previous_header = _header
            header = _header
        # END for each header we parse
        if index and header:
            index[-1].diff = text[header.end() :]
        # END assign last diff

        return index

    @staticmethod
    def _handle_diff_line(lines_bytes: bytes, repo: "Repo", index: DiffIndex) -> None:
        lines = lines_bytes.decode(defenc)

        # Discard everything before the first colon, and the colon itself.
        _, _, lines = lines.partition(":")

        for line in lines.split("\x00:"):
            if not line:
                # The line data is empty, skip.
                continue
            meta, _, path = line.partition("\x00")
            path = path.rstrip("\x00")
            a_blob_id: Optional[str]
            b_blob_id: Optional[str]
            old_mode, new_mode, a_blob_id, b_blob_id, _change_type = meta.split(None, 4)
            # Change type can be R100
            # R: status letter
            # 100: score (in case of copy and rename)
            change_type: Lit_change_type = cast(Lit_change_type, _change_type[0])
            score_str = "".join(_change_type[1:])
            score = int(score_str) if score_str.isdigit() else None
            path = path.strip()
            a_path = path.encode(defenc)
            b_path = path.encode(defenc)
            deleted_file = False
            new_file = False
            copied_file = False
            rename_from = None
            rename_to = None

            # NOTE: We cannot conclude from the existence of a blob to change type,
            # as diffs with the working do not have blobs yet.
            if change_type == "D":
                b_blob_id = None  # Optional[str]
                deleted_file = True
            elif change_type == "A":
                a_blob_id = None
                new_file = True
            elif change_type == "C":
                copied_file = True
                a_path_str, b_path_str = path.split("\x00", 1)
                a_path = a_path_str.encode(defenc)
                b_path = b_path_str.encode(defenc)
            elif change_type == "R":
                a_path_str, b_path_str = path.split("\x00", 1)
                a_path = a_path_str.encode(defenc)
                b_path = b_path_str.encode(defenc)
                rename_from, rename_to = a_path, b_path
            elif change_type == "T":
                # Nothing to do.
                pass
            # END add/remove handling

            diff = Diff(
                repo,
                a_path,
                b_path,
                a_blob_id,
                b_blob_id,
                old_mode,
                new_mode,
                new_file,
                deleted_file,
                copied_file,
                rename_from,
                rename_to,
                "",
                change_type,
                score,
            )
            index.append(diff)

    @classmethod
    def _index_from_raw_format(cls, repo: "Repo", proc: "Popen") -> "DiffIndex":
        """Create a new :class:`DiffIndex` from the given process output which must be
        in raw format.

        :param repo:
            The repository we are operating on.

        :param proc:
            Process to read output from.

        :return:
            :class:`DiffIndex`
        """
        # handles
        # :100644 100644 687099101... 37c5e30c8... M    .gitignore

        index: "DiffIndex" = DiffIndex()
        handle_process_output(
            proc,
            lambda byt: cls._handle_diff_line(byt, repo, index),
            None,
            finalize_process,
            decode_streams=False,
        )

        return index
