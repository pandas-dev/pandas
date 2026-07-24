# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

"""General repository-related functions."""

from __future__ import annotations

__all__ = [
    "rev_parse",
    "is_git_dir",
    "touch",
    "find_submodule_git_dir",
    "name_to_object",
    "short_to_long",
    "deref_tag",
    "to_commit",
    "find_worktree_git_dir",
]

import os
import os.path as osp
from pathlib import Path
import re
import stat
from string import digits

from gitdb.exc import BadName, BadObject

from git.cmd import Git
from git.exc import WorkTreeRepositoryUnsupported
from git.objects import Object
from git.objects.util import parse_date
from git.refs import SymbolicReference
from git.util import cygpath, bin_to_hex, hex_to_bin

# Typing ----------------------------------------------------------------------

from typing import Iterator, Optional, TYPE_CHECKING, Tuple, Union, cast, overload

from git.types import AnyGitObject, Literal, PathLike

if TYPE_CHECKING:
    from git.db import GitCmdObjectDB
    from git.objects import Commit
    from git.refs.reference import Reference
    from git.refs.log import RefLog, RefLogEntry
    from git.refs.tag import Tag

    from .base import Repo

# ----------------------------------------------------------------------------


def touch(filename: str) -> str:
    with open(filename, "ab"):
        pass
    return filename


def is_git_dir(d: PathLike) -> bool:
    """This is taken from the git setup.c:is_git_directory function.

    :raise git.exc.WorkTreeRepositoryUnsupported:
        If it sees a worktree directory. It's quite hacky to do that here, but at least
        clearly indicates that we don't support it. There is the unlikely danger to
        throw if we see directories which just look like a worktree dir, but are none.
    """
    if osp.isdir(d):
        if (osp.isdir(osp.join(d, "objects")) or "GIT_OBJECT_DIRECTORY" in os.environ) and osp.isdir(
            osp.join(d, "refs")
        ):
            headref = osp.join(d, "HEAD")
            return osp.isfile(headref) or (osp.islink(headref) and os.readlink(headref).startswith("refs"))
        elif (
            osp.isfile(osp.join(d, "gitdir"))
            and osp.isfile(osp.join(d, "commondir"))
            and osp.isfile(osp.join(d, "gitfile"))
        ):
            raise WorkTreeRepositoryUnsupported(d)
    return False


def find_worktree_git_dir(dotgit: PathLike) -> Optional[str]:
    """Search for a gitdir for this worktree."""
    try:
        statbuf = os.stat(dotgit)
    except OSError:
        return None
    if not stat.S_ISREG(statbuf.st_mode):
        return None

    try:
        lines = Path(dotgit).read_text().splitlines()
        for key, value in [line.strip().split(": ") for line in lines]:
            if key == "gitdir":
                return value
    except ValueError:
        pass
    return None


def find_submodule_git_dir(d: PathLike) -> Optional[PathLike]:
    """Search for a submodule repo."""
    if is_git_dir(d):
        return d

    try:
        with open(d) as fp:
            content = fp.read().rstrip()
    except IOError:
        # It's probably not a file.
        pass
    else:
        if content.startswith("gitdir: "):
            path = content[8:]

            if Git.is_cygwin():
                # Cygwin creates submodules prefixed with `/cygdrive/...`.
                # Cygwin git understands Cygwin paths much better than Windows ones.
                # Also the Cygwin tests are assuming Cygwin paths.
                path = cygpath(path)
            if not osp.isabs(path):
                path = osp.normpath(osp.join(osp.dirname(d), path))
            return find_submodule_git_dir(path)
    # END handle exception
    return None


def short_to_long(odb: "GitCmdObjectDB", hexsha: str) -> Optional[bytes]:
    """
    :return:
        Long hexadecimal sha1 from the given less than 40 byte hexsha, or ``None`` if no
        candidate could be found.

    :param hexsha:
        hexsha with less than 40 bytes.
    """
    try:
        return bin_to_hex(odb.partial_to_complete_sha_hex(hexsha))
    except BadObject:
        return None
    # END exception handling


def _describe_to_long(repo: "Repo", name: str) -> Optional[bytes]:
    """Resolve git-describe style names to the abbreviated object they contain."""
    match = re.match(r"^.+-\d+-g([0-9A-Fa-f]{4,40})(?:-dirty)?$", name)
    if match is None:
        match = re.match(r"^.+-g([0-9A-Fa-f]{4,40})(?:-dirty)?$", name)
    if match is None:
        match = re.match(r"^([0-9A-Fa-f]{4,40})-dirty$", name)
    if match is None:
        return None
    # END handle match

    hexsha = match.group(1)
    if len(hexsha) == 40:
        return hexsha.encode("ascii")
    return short_to_long(repo.odb, hexsha)


@overload
def name_to_object(repo: "Repo", name: str, return_ref: Literal[False] = ...) -> AnyGitObject: ...


@overload
def name_to_object(repo: "Repo", name: str, return_ref: Literal[True]) -> Union[AnyGitObject, SymbolicReference]: ...


def name_to_object(repo: "Repo", name: str, return_ref: bool = False) -> Union[AnyGitObject, SymbolicReference]:
    """
    :return:
        Object specified by the given name - hexshas (short and long) as well as
        references are supported.

    :param return_ref:
        If ``True``, and name specifies a reference, we will return the reference
        instead of the object. Otherwise it will raise :exc:`~gitdb.exc.BadObject` or
        :exc:`~gitdb.exc.BadName`.
    """
    hexsha: Union[None, str, bytes] = None

    # Is it a hexsha? Try the most common ones, which is 7 to 40.
    if repo.re_hexsha_shortened.match(name):
        if len(name) != 40:
            # Find long sha for short sha.
            hexsha = short_to_long(repo.odb, name)
        else:
            hexsha = name
        # END handle short shas
    # END find sha if it matches

    # If we couldn't find an object for what seemed to be a short hexsha, try to find it
    # as reference anyway, it could be named 'aaa' for instance.
    if hexsha is None:
        for base in (
            "%s",
            "refs/%s",
            "refs/tags/%s",
            "refs/heads/%s",
            "refs/remotes/%s",
            "refs/remotes/%s/HEAD",
        ):
            try:
                hexsha = SymbolicReference.dereference_recursive(repo, base % name)
                if return_ref:
                    return SymbolicReference(repo, base % name)
                # END handle symbolic ref
                break
            except ValueError:
                pass
        # END for each base
    # END handle hexsha

    if hexsha is None:
        hexsha = _describe_to_long(repo, name)
    # END handle describe output

    # Didn't find any ref, this is an error.
    if return_ref:
        raise BadObject("Couldn't find reference named %r" % name)
    # END handle return ref

    # Tried everything ? fail.
    if hexsha is None:
        raise BadName(name)
    # END assert hexsha was found

    return Object.new_from_sha(repo, hex_to_bin(hexsha))


def deref_tag(tag: "Tag") -> AnyGitObject:
    """Recursively dereference a tag and return the resulting object."""
    while True:
        try:
            tag = tag.object
        except AttributeError:
            break
    # END dereference tag
    return tag


def to_commit(obj: Object) -> "Commit":
    """Convert the given object to a commit if possible and return it."""
    if obj.type == "tag":
        obj = deref_tag(obj)

    if obj.type != "commit":
        raise ValueError("Cannot convert object %r to type commit" % obj)
    # END verify type
    return obj


def _object_from_hexsha(repo: "Repo", hexsha: str) -> AnyGitObject:
    return Object.new_from_sha(repo, hex_to_bin(hexsha))


def _current_reflog_ref(repo: "Repo") -> SymbolicReference:
    try:
        return repo.head.ref
    except TypeError:
        return repo.head
    # END handle detached head


def _common_reflog_path(repo: "Repo", ref: SymbolicReference) -> Optional[str]:
    if repo.common_dir == repo.git_dir:
        return None
    # END handle normal repository
    return SymbolicReference._get_validated_path(osp.join(repo.common_dir, "logs"), ref.path)


def _ref_log(repo: "Repo", ref: SymbolicReference) -> "RefLog":
    try:
        return ref.log()
    except FileNotFoundError:
        common_path = _common_reflog_path(repo, ref)
        if common_path and osp.isfile(common_path):
            from git.refs.log import RefLog

            return RefLog.from_file(common_path)
        # END handle linked-worktree branch logs
        try:
            if ref.path == repo.head.ref.path:
                return repo.head.log()
            # END handle linked-worktree current branch logs
        except TypeError:
            pass
        # END handle detached head
        raise
    # END handle missing branch log


def _ref_log_entry(repo: "Repo", ref: SymbolicReference, index: int) -> "RefLogEntry":
    try:
        return ref.log_entry(index)
    except FileNotFoundError:
        common_path = _common_reflog_path(repo, ref)
        if common_path and osp.isfile(common_path):
            from git.refs.log import RefLog

            return RefLog.entry_at(common_path, index)
        # END handle linked-worktree branch logs
        try:
            if ref.path == repo.head.ref.path:
                return repo.head.log_entry(index)
            # END handle linked-worktree current branch logs
        except TypeError:
            pass
        # END handle detached head
        raise
    # END handle missing branch log


def _find_reflog_entry_by_date(repo: "Repo", ref: SymbolicReference, spec: str) -> str:
    try:
        timestamp, _offset = parse_date(spec)
    except ValueError as e:
        raise NotImplementedError("Support for additional @{...} modes not implemented") from e
    # END handle unsupported dates
    log = _ref_log(repo, ref)
    if not log:
        raise IndexError("Invalid revlog date: %s" % spec)
    # END handle empty log

    for entry in reversed(log):
        if entry.time[0] <= timestamp:
            return entry.newhexsha
        # END found candidate
    # END for each entry
    return log[0].newhexsha


def _previous_checked_out_branch(repo: "Repo", nth: int) -> AnyGitObject:
    if nth <= 0:
        raise ValueError("Invalid previous checkout selector: -%i" % nth)
    # END handle invalid input

    seen = 0
    for entry in reversed(_ref_log(repo, repo.head)):
        message = entry.message or ""
        prefix = "checkout: moving from "
        if not message.startswith(prefix):
            continue
        # END skip non-checkouts

        previous_branch = message[len(prefix) :].split(" to ", 1)[0]
        seen += 1
        if seen == nth:
            return name_to_object(repo, previous_branch)
        # END found selector
    # END for each entry
    raise IndexError("Invalid previous checkout selector: -%i" % nth)


def _tracking_branch_object(repo: "Repo", ref: Optional[SymbolicReference]) -> AnyGitObject:
    from git.refs.head import Head

    if ref is None:
        try:
            head = repo.active_branch
        except TypeError as e:
            raise BadName("@{upstream}") from e
    elif isinstance(ref, Head):
        head = ref
    elif os.fspath(ref.path).startswith("refs/heads/"):
        head = Head(repo, ref.path)
    else:
        raise BadName("%s@{upstream}" % ref.name)
    # END handle head

    tracking_branch = head.tracking_branch()
    if tracking_branch is None:
        raise BadName("%s@{upstream}" % head.name)
    # END handle missing upstream
    return tracking_branch.commit


def _apply_reflog(repo: "Repo", ref: Optional[SymbolicReference], content: str) -> AnyGitObject:
    if content.startswith("+"):
        content = content[1:]
    # END handle explicit positive sign

    if content.startswith("-"):
        if ref is not None:
            raise ValueError("Previous checkout selectors do not take an explicit ref")
        if content == "-0":
            raise ValueError("Negative zero is invalid in reflog selector")
        # END handle invalid negative zero
        try:
            return _previous_checked_out_branch(repo, int(content[1:]))
        except ValueError as e:
            raise ValueError("Invalid previous checkout selector: %s" % content) from e
    # END handle previous checkout branch

    content_lower = content.lower()
    if content_lower in ("u", "upstream", "push"):
        return _tracking_branch_object(repo, ref)
    # END handle sibling branches

    ref = ref or _current_reflog_ref(repo)
    try:
        entry_no = int(content)
    except ValueError:
        hexsha = _find_reflog_entry_by_date(repo, ref, content)
    else:
        if entry_no >= 100000000:
            hexsha = _find_reflog_entry_by_date(repo, ref, "%s +0000" % entry_no)
        elif entry_no == 0:
            return ref.commit
        else:
            try:
                entry = _ref_log_entry(repo, ref, -(entry_no + 1))
            except IndexError as e:
                raise IndexError("Invalid revlog index: %i" % entry_no) from e
            # END handle index out of bound
            hexsha = entry.newhexsha
        # END handle offset or date-like timestamp
    # END handle content
    return _object_from_hexsha(repo, hexsha)


def _find_closing_brace(rev: str, start: int) -> int:
    depth = 1
    escaped = False
    for idx in range(start + 1, len(rev)):
        char = rev[idx]
        if escaped:
            escaped = False
        elif char == "\\":
            escaped = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return idx
            # END found end
        # END handle char
    # END for each char
    raise ValueError("Missing closing brace to define type in %s" % rev)


def _parse_search(pattern: str) -> Tuple[str, bool]:
    if not pattern:
        raise ValueError("Revision search requires a pattern")
    # END handle empty pattern

    if pattern.startswith("!-"):
        return pattern[2:], True
    if pattern.startswith("!!"):
        return pattern[1:], False
    if pattern.startswith("!"):
        raise ValueError("Need one character after /!, typically -")
    return pattern, False


def _unescape_braced_regex(pattern: str) -> str:
    out = []
    idx = 0
    while idx < len(pattern):
        char = pattern[idx]
        if char == "\\" and idx + 1 < len(pattern):
            next_char = pattern[idx + 1]
            if next_char in "{}\\":
                out.append(next_char)
            else:
                out.append(char)
                out.append(next_char)
            # END handle escaped char
            idx += 2
            continue
        # END handle backslash
        out.append(char)
        idx += 1
    # END for each char
    return "".join(out)


def _find_commit_by_message(
    repo: "Repo", rev: Optional[AnyGitObject], pattern: str, braced: bool = False
) -> AnyGitObject:
    pattern, negated = _parse_search(_unescape_braced_regex(pattern) if braced else pattern)
    try:
        regex = re.compile(pattern)
    except re.error as e:
        raise ValueError("Invalid commit message regex %r" % pattern) from e
    # END handle invalid regex
    if rev is None:
        commits = _all_ref_commits(repo)
    else:
        commits = _reachable_commits([to_commit(cast(Object, rev))])
    # END handle starting point

    for commit in commits:
        message = commit.message
        if isinstance(message, bytes):
            message = message.decode(commit.encoding, "replace")
        # END handle bytes message
        matches = regex.search(message or "") is not None
        if matches != negated:
            return commit
        # END found commit
    # END for each commit
    raise BadName("No commit found matching message pattern %r" % pattern)


def _all_ref_commits(repo: "Repo") -> Iterator["Commit"]:
    starts = []
    for ref in repo.references:
        try:
            starts.append(to_commit(cast(Object, ref.object)))
        except (BadName, ValueError):
            pass
        # END skip refs that do not point to commits
    # END for each ref
    try:
        starts.append(repo.head.commit)
    except ValueError:
        pass
    # END handle unborn head
    return _reachable_commits(starts)


def _reachable_commits(starts: list["Commit"]) -> Iterator["Commit"]:
    seen = set()
    pending = starts[:]
    while pending:
        pending.sort(key=lambda commit: commit.committed_date, reverse=True)
        commit = pending.pop(0)
        if commit.binsha in seen:
            continue
        # END skip seen commit
        seen.add(commit.binsha)
        yield commit
        pending.extend(commit.parents)
    # END while commits remain


def _index_lookup(repo: "Repo", spec: str) -> AnyGitObject:
    if not spec:
        raise ValueError("':' must be followed by a path")
    # END handle empty lookup

    stage = 0
    path = spec
    if len(spec) >= 2 and spec[1] == ":" and spec[0] in "0123":
        stage = int(spec[0])
        path = spec[2:]
    # END handle stage

    try:
        return repo.index.entries[(path, stage)].to_blob(repo)
    except KeyError as e:
        raise BadName("Path %r did not exist in the index at stage %i" % (path, stage)) from e


def _tree_lookup(obj: AnyGitObject, path: str) -> AnyGitObject:
    if obj.type != "tree":
        obj = to_commit(cast(Object, obj)).tree
    # END get tree
    if not path:
        return obj
    return obj[path]


def _peel(obj: AnyGitObject, output_type: str, repo: "Repo", rev: str) -> AnyGitObject:
    if output_type.startswith("/"):
        return _find_commit_by_message(repo, obj, output_type[1:], braced=True)
    if output_type == "":
        return deref_tag(obj) if obj.type == "tag" else obj
    if output_type == "object":
        return obj
    if output_type == "commit":
        return to_commit(cast(Object, obj))
    if output_type == "tree":
        return to_commit(cast(Object, obj)).tree if obj.type != "tree" else obj
    if output_type == "blob":
        obj = deref_tag(obj) if obj.type == "tag" else obj
        if obj.type == output_type:
            return obj
        # END handle matching type
        raise ValueError("Could not accommodate requested object type %r, got %s" % (output_type, obj.type))
    if output_type == "tag":
        if obj.type == output_type:
            return obj
        # END handle matching type
        raise ValueError("Could not accommodate requested object type %r, got %s" % (output_type, obj.type))
    # END handle known types
    raise ValueError("Invalid output type: %s ( in %s )" % (output_type, rev))


def _first_rev_token(rev: str) -> Optional[int]:
    for idx, char in enumerate(rev):
        if char in "^~:":
            return idx
        if char == "@":
            next_char = rev[idx + 1] if idx + 1 < len(rev) else None
            if idx == 0 and next_char in (None, "^", "~", ":", "{"):
                return idx
            if next_char == "{":
                return idx
            # END handle reflog selector
        # END handle at symbol
    # END for each char
    return None


def rev_parse(repo: "Repo", rev: str) -> AnyGitObject:
    """Parse a revision string. Like :manpage:`git-rev-parse(1)`.

    :return:
        `~git.objects.base.Object` at the given revision.

        This may be any type of git object:

        * :class:`Commit <git.objects.commit.Commit>`
        * :class:`TagObject <git.objects.tag.TagObject>`
        * :class:`Tree <git.objects.tree.Tree>`
        * :class:`Blob <git.objects.blob.Blob>`

    :param rev:
        :manpage:`git-rev-parse(1)`-compatible revision specification as string.
        Please see :manpage:`git-rev-parse(1)` for details.

    :raise gitdb.exc.BadObject:
        If the given revision could not be found.

    :raise ValueError:
        If `rev` couldn't be parsed.

    :raise IndexError:
        If an invalid reflog index is specified.
    """
    if rev.startswith(":/"):
        return _find_commit_by_message(repo, None, rev[2:])
    if rev.startswith(":"):
        return _index_lookup(repo, rev[1:])
    # END handle top-level colon modes

    obj: Optional[AnyGitObject] = None
    ref = None
    lr = len(rev)
    first_token = _first_rev_token(rev)
    if first_token is None:
        return name_to_object(repo, rev)
    # END handle plain name

    if first_token == 0:
        if rev[0] != "@":
            raise ValueError("Revision specifier must start with an object name: %s" % rev)
        # END handle invalid leading token
        ref = _current_reflog_ref(repo)
        obj = ref.commit
        start = 0 if rev.startswith("@{") else 1
    else:
        if rev[first_token] == "@":
            ref = cast("Reference", name_to_object(repo, rev[:first_token], return_ref=True))
            obj = ref.commit
        else:
            obj = name_to_object(repo, rev[:first_token])
        # END handle anchor
        start = first_token
    # END initialize anchor

    while start < lr:
        token = rev[start]

        if token == "@":
            if start + 1 >= lr or rev[start + 1] != "{":
                raise ValueError("Invalid @ token in revision specifier: %s" % rev)
            # END handle invalid @
            end = _find_closing_brace(rev, start + 1)
            obj = _apply_reflog(repo, ref if first_token != 0 and start == first_token else None, rev[start + 2 : end])
            ref = None
            start = end + 1
            continue
        # END handle reflog

        if token == ":":
            return _tree_lookup(obj, rev[start + 1 :])
        # END handle path

        start += 1

        if token == "^" and start < lr and rev[start] == "{":
            end = _find_closing_brace(rev, start)
            obj = _peel(obj, rev[start + 1 : end], repo, rev)
            ref = None
            start = end + 1
            continue
        # END parse type

        num = 0
        found_digit = False
        while start < lr:
            if rev[start] in digits:
                num = num * 10 + int(rev[start])
                start += 1
                found_digit = True
            else:
                break
            # END handle number
        # END number parse loop

        if not found_digit:
            num = 1
        # END set default num

        try:
            if token == "~":
                obj = to_commit(obj)
                for _ in range(num):
                    obj = obj.parents[0]
                # END for each history item to walk
            elif token == "^":
                obj = to_commit(obj)
                if num == 0:
                    pass
                else:
                    obj = obj.parents[num - 1]
                # END handle parent
            else:
                raise ValueError("Invalid token: %r" % token)
            # END end handle tag
        except (IndexError, AttributeError) as e:
            raise BadName(
                f"Invalid revision spec '{rev}' - not enough parent commits to reach '{token}{int(num)}'"
            ) from e
        # END exception handling
    # END parse loop

    if obj is None:
        raise ValueError("Revision specifier could not be parsed: %s" % rev)

    return obj
