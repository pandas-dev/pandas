from __future__ import annotations

import contextlib
import errno
import os
import secrets
import stat
import sys
import tempfile
import time
from dataclasses import dataclass
from errno import EACCES, EEXIST, ENOENT, ENOSYS, EPERM, ESTALE, EXDEV
from pathlib import Path
from typing import TYPE_CHECKING, Final, Literal, cast

from ._api import BaseFileLock, _canonical, _raise_cleanup_errors
from ._error import SoftFileLockProtocolError
from ._identity import host_name, process_start_token
from ._soft_protocol import STRICT_SOFT_SENTINEL_RECORD
from ._util import ensure_directory_exists, write_all

if TYPE_CHECKING:
    from collections.abc import Iterator

StrictSoftFileClaimState = Literal["intent", "held"]

_CLAIM_STATES: Final[frozenset[str]] = frozenset({"intent", "held"})
_COORDINATION_SUFFIX: Final[str] = ".filelock"
_CLAIM_DIRECTORY_NAME: Final[str] = "claims"
_CLAIM_MAGIC: Final[str] = "filelock-strict-v1"
_CLAIM_RECORD_LIMIT: Final[int] = 1024
_CLAIM_NAME_PART_COUNT: Final[int] = 3
_TOKEN_HEX_LENGTH: Final[int] = 32
_PRIVATE_RECORD_MARKER: Final[str] = ".private-v1-"
_PRIVATE_RECORD_SUFFIX: Final[str] = ".tmp"
_PRIVATE_RECORD_RANDOM_HEX_LENGTH: Final[int] = 32
_PRIVATE_RECORD_GRACE: Final[float] = 2.0
_UNLINK_MAX_RETRIES: Final[int] = 10
#: How long a scan waits out a claim held in Windows' delete-pending state before treating it as unreadable.
_CLAIM_READ_GRACE: Final[float] = 0.5
_CLAIM_READ_RETRY: Final[float] = 0.002
#: Windows opens descriptors in text mode by default, which rewrites newlines and truncates a record at a control byte.
#: The claim and sentinel records are exact binary, so every record descriptor must be binary; POSIX ignores the flag.
_O_BINARY: Final[int] = getattr(os, "O_BINARY", 0)
_LEGACY_SENTINEL: Final[bytes] = STRICT_SOFT_SENTINEL_RECORD.encode()
_WINDOWS_HARD_LINK_UNSUPPORTED: Final[frozenset[int]] = frozenset({1, 17, 50})

# Termux/Android CPython ships without os.link (bionic long had only linkat), so the strict backend's whole hard-link
# mechanism is absent there. Probe once and gate every os.link reference on it, so importing filelock still works and
# only an actual StrictSoftFileLock acquire reports the missing capability.
_HAS_LINK: Final[bool] = hasattr(os, "link")

# Probe dir_fd capability once at import. A per-call ``os.unlink in os.supports_dir_fd`` check flips to False the moment
# a test mocks os.unlink, silently diverting the code to a different branch than the one under test.
_OPEN_SUPPORTS_DIR_FD: Final[bool] = os.open in os.supports_dir_fd
_UNLINK_SUPPORTS_DIR_FD: Final[bool] = os.unlink in os.supports_dir_fd
_STAT_SUPPORTS_DIR_FD: Final[bool] = os.stat in os.supports_dir_fd
_LINK_SUPPORTS_DIR_FD: Final[bool] = _HAS_LINK and os.link in os.supports_dir_fd


def _probe_link_follow_symlinks() -> bool:
    # os.supports_follow_symlinks lists os.link on PyPy, but its linkat then rejects follow_symlinks=False with EINVAL,
    # and Windows raises NotImplementedError for the option outright. Link a throwaway file for real so the answer
    # reflects the runtime rather than its advertisement, and treat any failure as "not honored": the option only
    # hardens a source this process created with O_EXCL, so skipping it is safe, and a real environment fault surfaces
    # when the actual link runs.
    if not _HAS_LINK:
        return False
    try:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "probe-source")
            source.touch()
            os.link(source, Path(directory, "probe-link"), follow_symlinks=False)
    except (OSError, NotImplementedError, ValueError):
        return False
    return True


_LINK_HONORS_FOLLOW_SYMLINKS: Final[bool] = _probe_link_follow_symlinks()


def _probe_hard_link_unsupported_errnos() -> frozenset[int]:
    # GraalPy's errno omits ENOTSUP, so importing the name outright breaks every runtime that ships without it. ENOTSUP
    # wins wherever it exists, leaving every runtime that names it unchanged. EOPNOTSUPP only stands in for the ones
    # that do not, and it approximates rather than matches: the two codes agree on Linux but differ on macOS/BSD and on
    # Windows. Where neither exists a runtime that cannot name "operation not supported" cannot raise it either, and
    # ENOSYS/EXDEV still classify the link failures it can raise.
    not_supported = getattr(errno, "ENOTSUP", getattr(errno, "EOPNOTSUPP", None))
    return frozenset({ENOSYS, EXDEV} if not_supported is None else {ENOSYS, EXDEV, not_supported})


_HARD_LINK_UNSUPPORTED_ERRNOS: Final[frozenset[int]] = _probe_hard_link_unsupported_errnos()


class StrictSoftFileLock(BaseFileLock):
    """Portable fail-closed lock based on immutable owner claims."""

    _preserve_lock_file_supported: bool = True
    _on_acquired_supported: bool = False
    #: Age cannot clear a strict claim: expiring one on a clock is the overlap the fail-closed contract exists to rule
    #: out, so only force_break() removes it.
    _lifetime_supported: bool = False
    _lifetime_unsupported_reason: str = "a strict claim is never broken by age, only by force_break()"
    #: The claim doorway publishes an intent and a held record per owner, so a shared instance must serialize them.
    _serialize_transitions: bool = True
    #: Contending processes each publish and rescan several files, so back their retries off across a jittered window
    #: rather than let them collide on every poll. Seconds; keeps a waiter responsive once it wins.
    _poll_backoff_cap: float = 0.05

    def _acquire(self) -> None:
        # Resolve once per acquisition, not per poll: a waiter on a relative path must keep publishing into the
        # directory it started waiting in even when another thread changes the working directory mid-wait.
        if (claim_root := self._context.claim_root) is None:
            claim_root = self._context.claim_root = _canonical(self.lock_file)
        lock_path = Path(claim_root)
        coordination_directory = Path(f"{lock_path}{_COORDINATION_SUFFIX}")
        claim_directory = coordination_directory / _CLAIM_DIRECTORY_NAME
        ensure_directory_exists(os.fspath(lock_path))
        _ensure_protocol_directory(self.lock_file, coordination_directory)
        _ensure_protocol_directory(self.lock_file, claim_directory)
        if (sentinel_fd := _open_or_create_sentinel(self.lock_file, lock_path, self._open_mode())) is None:
            return
        try:
            sentinel_identity = _file_identity(os.fstat(sentinel_fd))
        except BaseException as inspection_error:  # preserve inspection and descriptor cleanup errors
            try:
                os.close(sentinel_fd)
            except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve inspection and descriptor cleanup errors
                _raise_cleanup_errors("strict sentinel inspection cleanup failed", inspection_error, close_error)
            raise
        self._mark_descriptor_pending(sentinel_fd, sentinel_identity)
        try:
            self._attempt_doorway(claim_directory, sentinel_fd, sentinel_identity)
        except BaseException:
            if self._context.pending_lock_file_fd == sentinel_fd:
                self._discard_doorway(sentinel_fd, sentinel_identity)
            raise

    def _attempt_doorway(self, claim_directory: Path, sentinel_fd: int, sentinel_identity: tuple[int, int]) -> None:
        if _read_existing_claims(self.lock_file, claim_directory):
            self._discard_doorway(sentinel_fd, sentinel_identity)
            return

        token = secrets.token_hex(_TOKEN_HEX_LENGTH // 2)
        intent_name = _claim_name("intent", token)
        intent_path = str(claim_directory / intent_name)
        try:
            publication_cleanup_error = _publish_record(intent_path, _claim_record(token), self._open_mode())
        except _PrivateRecordReclaimedError:
            self._discard_doorway(sentinel_fd, sentinel_identity)
            return
        except (NotImplementedError, OSError) as error:
            _raise_if_hard_links_unsupported(self.lock_file, error)
            if isinstance(error, OSError) and error.errno == EEXIST:
                self._discard_doorway(sentinel_fd, sentinel_identity)
                return
            raise
        if publication_cleanup_error is not None:  # pragma: needs dir-fd
            raise publication_cleanup_error
        self._context.owner_claim_paths = (intent_path,)

        claims = _read_existing_claims(self.lock_file, claim_directory)
        if (
            not claims
            or any(claim.state == "held" for claim in claims)
            or min(claim.name for claim in claims) != intent_name
        ):
            self._discard_doorway(sentinel_fd, sentinel_identity)
            return

        held_name = _claim_name("held", token)
        held_path = str(claim_directory / held_name)
        try:
            link_cleanup_error = _link_no_replace(claim_directory, intent_name, held_name)
        except (NotImplementedError, OSError) as error:
            _raise_if_hard_links_unsupported(self.lock_file, error)
            raise
        self._context.owner_claim_paths = (held_path, intent_path)
        if link_cleanup_error is not None:  # pragma: needs dir-fd
            self._context.owner_claim_paths = ()
            raise link_cleanup_error

        claims = _read_existing_claims(self.lock_file, claim_directory)
        if (
            not {intent_name, held_name}.issubset(claim.name for claim in claims)
            or min(_claim_token_key(claim.name) for claim in claims) != f"v1-{token}.claim"
        ):
            self._discard_doorway(sentinel_fd, sentinel_identity)
            return
        # Keep the intent claim for the whole hold rather than unlinking it now. The intent has existed, unchanged,
        # since this owner published it, so a contender's os.scandir is guaranteed to return it (POSIX only leaves the
        # visibility of entries created or removed *during* a scan unspecified). The freshly linked held claim carries
        # no such guarantee: a scan that races its creation can miss it. Were the intent removed here, that scan could
        # observe neither claim and let a larger-token contender win over this owner. The stable intent is the witness
        # that keeps the phase-five min-token decision computed over the true set. Release unlinks both.
        self._mark_descriptor_owned(sentinel_fd, sentinel_identity)

    @property
    def claims(self) -> tuple[StrictSoftFileClaim, ...]:
        """Published claims that block acquisition."""
        return _read_existing_claims(self.lock_file, self._claim_directory)

    def force_break(self, claim_name: str) -> None:
        """Remove one named claim, allowing overlap if its owner still holds the protected resource."""
        _validate_force_break_name(claim_name)
        _require_exact_name(self._claim_directory, claim_name)
        if (
            cleanup_error := _unlink_in_directory(self._claim_directory, claim_name)
        ) is not None:  # pragma: needs dir-fd
            raise cleanup_error

    def _rollback_failed_acquire(self, acquisition_error: BaseException) -> None:
        # _acquire already reconciles a failed doorway through _discard_doorway: it either closes the pending
        # descriptor or, when a held claim cannot be removed, commits it as owned so a later release retries and
        # raises the cleanup errors. A base rollback would release that owned descriptor again and report each
        # failure a second time, so leave the reconciled state alone.
        if self.is_locked:
            return
        super()._rollback_failed_acquire(acquisition_error)

    def _reconcile_failed_acquire(self, canonical: str) -> None:
        # The acquisition is over, so the next one resolves the working directory again rather than reuse this one's.
        if not self.is_locked:
            self._context.claim_root = None
        super()._reconcile_failed_acquire(canonical)

    def _release(self) -> None:
        fd = cast("int", self._context.lock_file_fd)
        self._context.claim_root = None
        remaining, errors = _unlink_owner_paths(self._context.owner_claim_paths)
        self._context.owner_claim_paths = tuple(remaining)
        if remaining:
            _raise_recorded_errors("strict claim release failed", errors)
        self._mark_descriptor_released()
        try:
            self._close_released_fd(fd, default_suppresses=False)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve claim and sentinel cleanup errors
            errors.append(close_error)
        if errors:
            _raise_recorded_errors("strict release cleanup failed", errors)

    def _discard_doorway(self, fd: int, identity: tuple[int, int]) -> None:
        remaining, errors = _unlink_owner_paths(self._context.owner_claim_paths)
        self._context.owner_claim_paths = tuple(remaining)
        if remaining:
            self._mark_descriptor_owned(fd, identity)
            _raise_recorded_errors("strict doorway claim cleanup failed", errors)
        self._mark_descriptor_released()
        try:
            self._close_released_fd(fd, default_suppresses=False)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve claim and sentinel cleanup errors
            errors.append(close_error)
        if errors:
            _raise_recorded_errors("strict doorway cleanup failed", errors)

    @property
    def _claim_directory(self) -> Path:
        if self._context.owner_claim_paths:
            return Path(self._context.owner_claim_paths[0]).parent
        return Path(f"{_canonical(self.lock_file)}{_COORDINATION_SUFFIX}") / _CLAIM_DIRECTORY_NAME


@dataclass(frozen=True)
class StrictSoftFileClaim:
    """One parsed strict soft-lock claim."""

    name: str
    state: StrictSoftFileClaimState
    token: str
    pid: int
    hostname: str
    #: The owner's process start token, or ``None`` when the platform exposes no proven start time. A strict lock never
    #: reclaims a claim on its own, so this identifies the owner for tooling rather than driving any automatic break.
    start: int | None = None


class _PrivateRecordReclaimedError(Exception):
    pass


def _open_or_create_sentinel(lock_file: str, path: Path, mode: int) -> int | None:
    try:
        return _open_sentinel(path)
    except FileNotFoundError:
        pass
    except OSError:
        return None

    _reclaim_sentinel_private_records(path, time.time())
    try:
        publication_cleanup_error = _publish_record(os.fspath(path), _LEGACY_SENTINEL, mode)
    except _PrivateRecordReclaimedError:
        return None
    except (NotImplementedError, OSError) as error:
        _raise_if_hard_links_unsupported(lock_file, error)
        if not isinstance(error, OSError) or error.errno != EEXIST:
            raise
    else:
        if publication_cleanup_error is not None:  # pragma: needs dir-fd
            raise publication_cleanup_error
    try:
        return _open_sentinel(path)
    except OSError:
        return None


def _open_sentinel(path: Path) -> int | None:
    fd, record = _open_record(path, len(_LEGACY_SENTINEL))
    if record == _LEGACY_SENTINEL:
        return fd
    os.close(fd)
    return None


def _read_claims(lock_file: str, directory: Path) -> tuple[StrictSoftFileClaim, ...]:
    try:
        with os.scandir(directory) as entries:
            names = _public_claim_names(directory, entries)
    except OSError as error:
        reason = f"cannot list claim directory: {error.strerror or type(error).__name__}"
        raise SoftFileLockProtocolError(lock_file, None, reason) from error

    claims: list[StrictSoftFileClaim] = []
    for name in names:
        if (name_parts := _parse_claim_name(name)) is None:
            raise SoftFileLockProtocolError(lock_file, name, "unknown claim name or protocol version")
        if (record := _read_claim_record(lock_file, directory, name)) is not None:
            claims.append(_parse_claim(lock_file, name, name_parts, record))
    return tuple(claims)


def _read_claim_record(lock_file: str, directory: Path, name: str) -> bytes | None:
    # A contended scan can list a claim that is not yet cleanly readable, in two ways that both resolve on a brief
    # retry. Windows holds a claim mid-unlink in a delete-pending state that fails an open with EACCES until the unlink
    # completes. On NFS, a peer that unlinks its own claim leaves this client's cached filehandle stale, so the next
    # open returns ESTALE rather than a clean ENOENT until the lookup revalidates against the server. Retrying re-runs
    # the path lookup, which turns the vanished claim into ENOENT (skip) or reads it if it still exists. A genuinely
    # unreadable record (a locked-down file, an EIO fault) still fails closed.
    deadline = time.monotonic() + _CLAIM_READ_GRACE
    delaying = False
    while True:
        if delaying:
            time.sleep(_CLAIM_READ_RETRY)
        record, pending = _attempt_claim_read(lock_file, directory, name)
        if pending is None:
            return record
        if time.monotonic() >= deadline:
            if pending.errno == ESTALE:
                # A stale handle that outlives revalidation is a claim the server no longer has (RFC 1813
                # NFS3ERR_STALE): skip it like ENOENT. Skipping a peer's vanished claim can only overcount
                # contention, never free a held lock.
                return None
            reason = f"cannot read claim: {pending.strerror or str(pending) or type(pending).__name__}"
            raise SoftFileLockProtocolError(lock_file, name, reason) from pending
        delaying = True


def _attempt_claim_read(lock_file: str, directory: Path, name: str) -> tuple[bytes | None, OSError | None]:
    # Return the record, or (None, None) when the claim has already gone, or (None, error) for an open the caller may
    # retry: a Windows delete-pending EACCES that resolves to a clean removal, or an NFS ESTALE from a peer unlinking
    # its own claim out from under this client's cached filehandle. Any other OSError is a real fault and fails closed.
    try:
        return _read_record(directory / name, _CLAIM_RECORD_LIMIT), None
    except FileNotFoundError:
        return None, None
    except PermissionError as error:
        return None, error
    except OSError as error:
        if error.errno == ESTALE:
            return None, error
        reason = f"cannot read claim: {error.strerror or str(error) or type(error).__name__}"
        raise SoftFileLockProtocolError(lock_file, name, reason) from error


def _public_claim_names(directory: Path, entries: Iterator[os.DirEntry[str]]) -> list[str]:
    names: list[str] = []
    now = time.time()
    for entry in entries:
        if not entry.name.startswith("."):
            names.append(entry.name)
        elif (public_name := _private_public_name(entry.name)) is not None and _parse_claim_name(
            public_name
        ) is not None:
            _reclaim_private_record((os.fspath(directory), None), entry.name, now)
    return sorted(names)


def _read_existing_claims(lock_file: str, directory: Path) -> tuple[StrictSoftFileClaim, ...]:
    if not directory.exists():
        return ()
    return _read_claims(lock_file, directory)


def _parse_claim(
    lock_file: str,
    name: str,
    name_parts: tuple[StrictSoftFileClaimState, str],
    record: bytes,
) -> StrictSoftFileClaim:
    try:
        magic, token, pid_text, hostname_hex, start_text, trailing = record.decode("ascii").split("\n")
        pid = int(pid_text)
        start = int(start_text) if start_text else None
        hostname = bytes.fromhex(hostname_hex).decode("utf-8")
    except (UnicodeDecodeError, ValueError) as error:
        raise SoftFileLockProtocolError(lock_file, name, "malformed claim record") from error
    if not all((
        not trailing,
        magic == _CLAIM_MAGIC,
        token == name_parts[1],
        1 <= pid <= 2**31 - 1,
        str(pid) == pid_text,
        start is None or (start >= 0 and str(start) == start_text),
        hostname.encode().hex() == hostname_hex
        and hostname.isprintable()
        and not any(character.isspace() for character in hostname),
    )):
        raise SoftFileLockProtocolError(lock_file, name, "malformed claim record")
    return StrictSoftFileClaim(name=name, state=name_parts[0], token=token, pid=pid, hostname=hostname, start=start)


def _parse_claim_name(name: str) -> tuple[StrictSoftFileClaimState, str] | None:
    if not name.endswith(".claim"):
        return None
    parts = name.removesuffix(".claim").split("-")
    if len(parts) != _CLAIM_NAME_PART_COUNT or parts[0] not in _CLAIM_STATES or parts[1] != "v1":
        return None
    token = parts[2]
    if len(token) != _TOKEN_HEX_LENGTH or any(character not in "0123456789abcdef" for character in token):
        return None
    return cast("StrictSoftFileClaimState", parts[0]), token


def _claim_name(state: StrictSoftFileClaimState, token: str) -> str:
    return f"{state}-v1-{token}.claim"


def _claim_token_key(name: str) -> str:
    return name.removeprefix("held-").removeprefix("intent-")


def _claim_record(token: str) -> bytes:
    hostname_hex = host_name().encode().hex()
    start = process_start_token(os.getpid())
    start_text = "" if start is None else str(start)
    return f"{_CLAIM_MAGIC}\n{token}\n{os.getpid()}\n{hostname_hex}\n{start_text}\n".encode("ascii")


def _publish_record(
    public_path: str,
    record: bytes,
    mode: int,
) -> BaseException | None:
    directory, public_name = os.path.split(public_path)
    directory = directory or os.curdir
    private_name = _private_record_name(public_name)
    directory_fd = _open_directory(directory) if _OPEN_SUPPORTS_DIR_FD else None
    directory_ref = directory, directory_fd
    try:
        _publish_record_in_directory(directory_ref, (private_name, public_name), mode, record)
    except BaseException as publication_error:  # preserve publication and directory cleanup errors
        try:
            if directory_fd is not None:  # pragma: needs dir-fd
                os.close(directory_fd)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # pragma: needs dir-fd  # preserve publication and directory cleanup errors
            _raise_cleanup_errors("strict publication directory cleanup failed", publication_error, close_error)
        raise
    if directory_fd is not None:  # pragma: needs dir-fd
        try:
            os.close(directory_fd)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # caller records the published path before raising
            return close_error
    return None


def _publish_record_in_directory(
    directory_ref: tuple[str, int | None],
    names: tuple[str, str],
    mode: int,
    record: bytes,
) -> None:
    flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | _O_BINARY
    if (o_nofollow := getattr(os, "O_NOFOLLOW", None)) is not None:  # pragma: needs o-nofollow
        flags |= o_nofollow
    private_fd = _open_relative(directory_ref, names[0], flags, mode)
    private_identity: tuple[int, int] | None = None
    try:
        private_identity = _file_identity(os.fstat(private_fd))
        write_all(private_fd, record)
        _link_private_record(directory_ref, names, private_identity)
    except BaseException as publication_error:  # preserve publication and cleanup errors
        close_error, unlink_error = _close_and_unlink_private_record(
            directory_ref,
            names[0],
            private_fd,
            private_identity,
        )
        if close_error is not None or unlink_error is not None:
            _raise_cleanup_errors(
                "strict record publication cleanup failed",
                publication_error,
                close_error,
                unlink_error,
            )
        raise
    close_error, unlink_error = _close_and_unlink_private_record(
        directory_ref,
        names[0],
        private_fd,
        private_identity,
    )
    if close_error is not None or unlink_error is not None:
        _raise_record_finalization_errors(close_error, unlink_error)


def _close_and_unlink_private_record(
    directory_ref: tuple[str, int | None],
    private_name: str,
    private_fd: int,
    private_identity: tuple[int, int] | None,
) -> tuple[BaseException | None, BaseException | None]:
    close_error: BaseException | None = None
    try:
        os.close(private_fd)
    except BaseException as error:  # ruff:ignore[blind-except]  # returned for grouping with unlink failures
        close_error = error
    unlink_error: BaseException | None = None
    try:
        if private_identity is None:
            _unlink_relative(directory_ref, private_name)
        else:
            _unlink_relative_if_identity(directory_ref, private_name, private_identity)
    except FileNotFoundError:
        pass
    except BaseException as error:  # ruff:ignore[blind-except]  # returned for grouping with close failures
        unlink_error = error
    return close_error, unlink_error


def _link_private_record(
    directory_ref: tuple[str, int | None],
    names: tuple[str, str],
    private_identity: tuple[int, int],
) -> None:
    try:
        _link_relative(directory_ref, *names)
    except FileNotFoundError as error:
        if _relative_identity(directory_ref, names[0]) is not None:
            raise
        msg = "private publication record was reclaimed"
        raise _PrivateRecordReclaimedError(msg) from error
    if _relative_identity(directory_ref, names[1]) == private_identity:
        return
    msg_0 = "private publication record was replaced"
    raise _PrivateRecordReclaimedError(msg_0)


def _raise_record_finalization_errors(
    close_error: BaseException | None,
    unlink_error: BaseException | None,
) -> None:
    errors = [error for error in (close_error, unlink_error) if error is not None]
    if len(errors) > 1:
        _raise_cleanup_errors("strict record finalization failed", errors[0], *errors[1:])
    raise errors[0]


def _private_record_name(public_name: str) -> str:
    return f".{public_name}{_PRIVATE_RECORD_MARKER}{secrets.token_hex(_PRIVATE_RECORD_RANDOM_HEX_LENGTH // 2)}.tmp"


def _private_public_name(private_name: str) -> str | None:
    if not private_name.startswith(".") or not private_name.endswith(_PRIVATE_RECORD_SUFFIX):
        return None
    public_name, marker, random_hex = private_name[1 : -len(_PRIVATE_RECORD_SUFFIX)].rpartition(_PRIVATE_RECORD_MARKER)
    if (
        marker != _PRIVATE_RECORD_MARKER
        or len(random_hex) != _PRIVATE_RECORD_RANDOM_HEX_LENGTH
        or any(character not in "0123456789abcdef" for character in random_hex)
    ):
        return None
    return public_name


def _reclaim_sentinel_private_records(path: Path, now: float) -> None:
    directory_ref = os.fspath(path.parent), None
    with os.scandir(path.parent) as entries:
        for entry in entries:
            if _private_public_name(entry.name) == path.name:
                _reclaim_private_record(directory_ref, entry.name, now)


def _reclaim_private_record(directory_ref: tuple[str, int | None], private_name: str, now: float) -> None:
    directory, directory_fd = directory_ref
    try:
        private_stat = (
            os.stat(private_name, dir_fd=directory_fd, follow_symlinks=False)
            if directory_fd is not None and _STAT_SUPPORTS_DIR_FD
            else Path(directory, private_name).lstat()
        )
    except FileNotFoundError:
        return
    if not stat.S_ISREG(private_stat.st_mode):
        msg = f"{Path(directory, private_name)} is not a regular private record"
        raise OSError(msg)
    if private_stat.st_nlink == 1 and now - private_stat.st_mtime < _PRIVATE_RECORD_GRACE:
        return
    try:
        _unlink_private_record_once(directory_ref, private_name)
    except FileNotFoundError:
        pass
    except OSError as error:
        if error.errno not in {EACCES, EPERM}:
            raise


def _unlink_private_record_once(directory_ref: tuple[str, int | None], private_name: str) -> None:
    directory, directory_fd = directory_ref
    if (
        directory_fd is not None and _UNLINK_SUPPORTS_DIR_FD
    ):  # pragma: no cover  # callers always pass directory_fd=None
        os.unlink(private_name, dir_fd=directory_fd)
    else:
        Path(directory, private_name).unlink()


def _relative_identity(directory_ref: tuple[str, int | None], name: str) -> tuple[int, int] | None:
    directory, directory_fd = directory_ref
    try:
        if directory_fd is not None and _STAT_SUPPORTS_DIR_FD:  # pragma: needs dir-fd
            path_stat = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        else:  # pragma: win32 cover
            path_stat = Path(directory, name).lstat()
    except FileNotFoundError:
        return None
    return _file_identity(path_stat)


def _unlink_relative_if_identity(
    directory_ref: tuple[str, int | None],
    name: str,
    identity: tuple[int, int],
) -> None:
    if _relative_identity(directory_ref, name) == identity:
        with contextlib.suppress(FileNotFoundError):
            _unlink_relative(directory_ref, name)


def _open_relative(directory_ref: tuple[str, int | None], name: str, flags: int, mode: int) -> int:
    directory, directory_fd = directory_ref
    return (
        os.open(name, flags, mode, dir_fd=directory_fd)
        if directory_fd is not None
        else os.open(Path(directory, name), flags, mode)
    )


def _link_relative(directory_ref: tuple[str, int | None], source_name: str, destination_name: str) -> None:
    directory, directory_fd = directory_ref
    if directory_fd is not None and _LINK_SUPPORTS_DIR_FD:  # pragma: needs dir-fd
        _link_no_follow(source_name, destination_name, src_dir_fd=directory_fd, dst_dir_fd=directory_fd)
        return
    _link_no_follow(Path(directory, source_name), Path(directory, destination_name))  # pragma: win32 cover


def _link_no_follow(
    source: str | Path,
    destination: str | Path,
    *,
    src_dir_fd: int | None = None,
    dst_dir_fd: int | None = None,
) -> None:
    if not _HAS_LINK:
        # No os.link at all (Termux/Android): report it as unsupported like a filesystem that refuses hard links, so
        # the acquire path raises SoftFileLockProtocolError instead of a bare AttributeError.
        msg = "os.link is unavailable on this platform"
        raise NotImplementedError(msg)
    # The source is a private record this process created with O_CREAT | O_EXCL, so follow_symlinks guards nothing an
    # attacker can reach. Pass the option only when the runtime honors it: PyPy advertises it through
    # os.supports_follow_symlinks yet its linkat rejects it with EINVAL, so probe once rather than trust the set.
    if _LINK_HONORS_FOLLOW_SYMLINKS:
        os.link(source, destination, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd, follow_symlinks=False)
    else:  # pragma: lacks link-follow-symlinks
        os.link(source, destination, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)


def _unlink_relative(directory_ref: tuple[str, int | None], name: str) -> None:
    directory, directory_fd = directory_ref
    if directory_fd is not None and _UNLINK_SUPPORTS_DIR_FD:  # pragma: needs dir-fd
        os.unlink(name, dir_fd=directory_fd)
    elif sys.platform == "win32":  # pragma: win32 cover
        _unlink_in_directory(Path(directory), name)
    else:
        Path(directory, name).unlink()  # pragma: win32 no cover


def _link_no_replace(directory: Path, source_name: str, destination_name: str) -> BaseException | None:
    if _LINK_SUPPORTS_DIR_FD:  # pragma: needs dir-fd
        directory_fd = _open_directory(str(directory))
        try:
            _link_no_follow(source_name, destination_name, src_dir_fd=directory_fd, dst_dir_fd=directory_fd)
        except BaseException as link_error:  # preserve link and directory cleanup errors
            try:
                os.close(directory_fd)
            except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve link and directory cleanup errors
                _raise_cleanup_errors("strict link directory cleanup failed", link_error, close_error)
            raise
        try:
            os.close(directory_fd)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # caller records the held path before raising
            return close_error
        return None
    _link_no_follow(directory / source_name, directory / destination_name)  # pragma: win32 cover
    return None  # pragma: win32 cover


def _unlink_owner_path(path: str) -> BaseException | None:
    directory, name = os.path.split(path)
    return _unlink_in_directory(Path(directory or os.curdir), name)


def _unlink_owner_paths(paths: tuple[str, ...]) -> tuple[list[str], list[BaseException]]:
    results = [(path, _unlink_owner_path_result(path)) for path in paths]
    return [path for path, (removed, _) in results if not removed], [
        error for _, (_, error) in results if error is not None
    ]


def _unlink_owner_path_result(path: str) -> tuple[bool, BaseException | None]:
    try:
        cleanup_error = _unlink_owner_path(path)
    except FileNotFoundError:
        return True, None
    except BaseException as error:  # ruff:ignore[blind-except]  # keep ownership when unlink did not commit
        return False, error
    return True, cleanup_error


def _raise_recorded_errors(message: str, errors: list[BaseException]) -> None:
    if len(errors) > 1:
        _raise_cleanup_errors(message, errors[0], *errors[1:])
    raise errors[0]


def _unlink_in_directory(directory: Path, name: str) -> BaseException | None:
    if _UNLINK_SUPPORTS_DIR_FD:  # pragma: needs dir-fd
        directory_fd = _open_directory(str(directory))
        try:
            os.unlink(name, dir_fd=directory_fd)
        except BaseException as unlink_error:  # preserve unlink and directory cleanup errors
            try:
                os.close(directory_fd)
            except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve unlink and directory cleanup errors
                _raise_cleanup_errors("strict unlink directory cleanup failed", unlink_error, close_error)
            raise
        try:
            os.close(directory_fd)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # caller commits the removed path before raising
            return close_error
        return None
    if sys.platform != "win32":  # pragma: win32 no cover
        Path(directory / name).unlink()
        return None
    retry_delay = 0.001  # pragma: win32 cover
    for attempt in range(_UNLINK_MAX_RETRIES):  # pragma: win32 cover
        if (error := _unlink_error(directory / name)) is None:
            return None
        if error.errno not in {EACCES, EPERM} or attempt == _UNLINK_MAX_RETRIES - 1:
            raise error
        time.sleep(retry_delay)
        retry_delay *= 2
    return None  # pragma: no cover  # the final retry returns or raises


def _unlink_error(path: Path) -> OSError | None:  # pragma: win32 cover
    try:
        path.unlink()
    except OSError as error:
        return error
    return None


def _open_directory(directory: str) -> int:  # pragma: needs dir-fd
    return os.open(
        directory,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )


def _read_record(path: Path, limit: int) -> bytes:
    fd, record = _open_record(path, limit)
    os.close(fd)
    return record


def _open_record(path: Path, limit: int) -> tuple[int, bytes]:
    path_stat = path.lstat()
    if not stat.S_ISREG(path_stat.st_mode):
        msg = f"{path} is not a regular file"
        raise OSError(msg)
    flags = os.O_RDONLY | _O_BINARY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    fd = os.open(path, flags)
    try:
        record = _read_opened_record(fd, path, path_stat, limit)
    except BaseException as read_error:  # preserve read and descriptor cleanup errors
        try:
            os.close(fd)
        except BaseException as close_error:  # ruff:ignore[blind-except]  # preserve read and descriptor cleanup errors
            _raise_cleanup_errors("strict record read cleanup failed", read_error, close_error)
        raise
    return fd, record


def _read_opened_record(fd: int, path: Path, path_stat: os.stat_result, limit: int) -> bytes:
    opened_stat = os.fstat(fd)
    if not stat.S_ISREG(opened_stat.st_mode):
        msg = f"{path} is not a regular file"
        raise OSError(msg)
    if _file_identity(opened_stat) != _file_identity(path_stat):
        msg = f"{path} changed while opening"
        raise OSError(msg)
    record = os.read(fd, limit + 1)
    if len(record) > limit:
        msg = f"{path} exceeds {limit} bytes"
        raise OSError(msg)
    return record


def _ensure_protocol_directory(lock_file: str, directory: Path) -> None:
    try:
        directory.mkdir()
    except FileExistsError:
        try:
            mode = directory.lstat().st_mode
        except OSError as error:
            raise SoftFileLockProtocolError(lock_file, None, f"cannot inspect {directory}") from error
        if stat.S_ISDIR(mode) and not stat.S_ISLNK(mode):
            return
        raise SoftFileLockProtocolError(lock_file, None, f"{directory} is not a real directory") from None


def _validate_force_break_name(name: str) -> None:
    invalid_component = not name or name.startswith(".") or name in {os.curdir, os.pardir}
    if invalid_component or any(separator in name for separator in ("/", "\\", "\x00")):
        msg = "claim_name must be one public claim basename"
        raise ValueError(msg)


def _require_exact_name(directory: Path, name: str) -> None:
    with os.scandir(directory) as entries:
        if not any(entry.name == name for entry in entries):
            raise FileNotFoundError(ENOENT, os.strerror(ENOENT), name)


def _raise_if_hard_links_unsupported(lock_file: str, error: NotImplementedError | OSError) -> None:
    unsupported = (
        isinstance(error, NotImplementedError)
        or error.errno in _HARD_LINK_UNSUPPORTED_ERRNOS
        or (getattr(error, "winerror", None) in _WINDOWS_HARD_LINK_UNSUPPORTED)
    )
    if unsupported:
        reason = "filesystem does not support atomic no-replace hard-link publication"
        raise SoftFileLockProtocolError(lock_file, None, reason) from error


def _file_identity(st: os.stat_result) -> tuple[int, int]:
    return st.st_dev, st.st_ino


__all__ = [
    "StrictSoftFileClaim",
    "StrictSoftFileClaimState",
    "StrictSoftFileLock",
]
