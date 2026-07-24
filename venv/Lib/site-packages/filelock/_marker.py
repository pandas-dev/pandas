from __future__ import annotations

import math
import os
from contextlib import suppress
from typing import Final, Literal, NamedTuple

from ._identity import host_name, process_start_token
from ._soft import SoftFileLock, _read_lock_file
from ._util import write_all

#: Protocol 1 is the legacy ``<pid>\n<hostname>\n[<start_token>\n]`` marker that :class:`SoftFileLock` still writes.
#: Protocol 2 carries the owner mode and the lease claim. A protocol 1 reader treats a protocol 2 marker as malformed
#: and evicts it after its grace period, so the two never guarantee mutual exclusion against each other.
_PROTOCOL: Final[str] = "filelock/2"

_MAX_PID: Final[int] = 2**31 - 1

#: ``unknown`` is never published: it names a mode some other filelock wrote that this version cannot interpret. Such a
#: record still identifies a live owner, so it is parsed rather than read as malformed and aged out.
OwnerMode = Literal["lease", "unknown"]


class OwnerRecord(NamedTuple):
    """The owner published in a protocol 2 marker."""

    pid: int
    hostname: str
    mode: OwnerMode
    token: str | None = None
    lease_duration: float | None = None
    start: int | None = None


class MarkerSoftFileLock(SoftFileLock):
    """An existence lock whose marker carries a protocol 2 owner record."""

    #: Filled in by each mode so the published record states the contract its holder acquired under.
    _owner_mode: OwnerMode

    @property
    def owner(self) -> OwnerRecord | None:
        """
        The owner named by the marker on disk.

        :returns: the published record, or ``None`` when no marker exists or its record is malformed or protocol 1

        """
        return self._read_owner()

    @property
    def pid(self) -> int | None:
        """
        The PID of the process holding this lock, read from the marker.

        :returns: the PID, or ``None`` when no marker exists or its record is unreadable

        """
        return None if (owner := self._read_owner()) is None else owner.pid

    @property
    def is_lock_held_by_us(self) -> bool:
        """
        Whether the marker on disk names this process.

        :returns: ``True`` when the marker's PID and hostname match this process

        """
        owner = self._read_owner()
        return owner is not None and owner.pid == os.getpid() and owner.hostname == host_name()

    def force_break(self) -> None:
        """
        Remove the marker whoever holds it, so a later contender can acquire.

        Forced breaking voids mutual exclusion: the previous holder keeps running and keeps using whatever the lock
        protects. Reserve it for an operator clearing a marker whose holder is known to be gone.
        """
        self.break_lock()

    def _read_owner(self) -> OwnerRecord | None:
        with suppress(OSError, ValueError):
            return parse_marker(_read_lock_file(self.lock_file)[0])
        return None

    def _write_lock_info(self, fd: int) -> None:
        write_all(fd, encode_marker(self._published_record()))

    def _published_record(self) -> OwnerRecord:
        return OwnerRecord(
            pid=os.getpid(),
            hostname=host_name(),
            mode=self._owner_mode,
            start=process_start_token(os.getpid()),
        )


def encode_marker(record: OwnerRecord) -> bytes:
    """Render an owner record as the bytes a protocol 2 marker holds."""
    lines = [_PROTOCOL, f"pid={record.pid}", f"host={record.hostname}", f"mode={record.mode}"]
    if record.token is not None:
        lines.append(f"token={record.token}")
    if record.lease_duration is not None:
        lines.append(f"duration={record.lease_duration!r}")
    if record.start is not None:
        lines.append(f"start={record.start}")
    return "".join(f"{line}\n" for line in lines).encode()


def parse_marker(content: str | None) -> OwnerRecord | None:
    """Return the owner a protocol 2 marker names, or ``None`` when the record is malformed or protocol 1."""
    if not content or not (lines := content.strip().splitlines()) or lines[0] != _PROTOCOL:
        return None
    fields: dict[str, str] = {}
    for line in lines[1:]:
        key, separator, value = line.partition("=")
        if not separator:
            return None
        fields[key] = value
    return _build_record(fields)


def _build_record(fields: dict[str, str]) -> OwnerRecord | None:
    # An unknown key is a field a newer filelock published, so ignore it rather than read the record as malformed. An
    # unrecognized mode is the same story one level up: a contract this version does not implement. Reading it as
    # malformed would age the marker out of a live owner's hands, so keep it and let the caller refuse to reclaim it.
    # A record naming no mode at all states no contract and stays malformed.
    if (published := fields.get("mode")) is None:
        return None
    mode: OwnerMode = "lease" if published == "lease" else "unknown"
    hostname = fields.get("host")
    if not hostname or "pid" not in fields:
        return None
    try:
        pid = int(fields["pid"])
        duration = float(fields["duration"]) if "duration" in fields else None
        start = int(fields["start"]) if "start" in fields else None
    except ValueError:
        return None
    if not 1 <= pid <= _MAX_PID:
        return None
    token = fields.get("token")
    # float() accepts "nan" and "inf", and neither is non-positive, so a duration <= 0 guard alone would read such a
    # marker as a valid lease. A nan duration mismatches every configured duration and so wedges reclaim, where a
    # malformed marker ages out through the grace window.
    if mode == "lease" and (token is None or duration is None or not (math.isfinite(duration) and duration > 0)):
        return None
    return OwnerRecord(pid=pid, hostname=hostname, mode=mode, token=token, lease_duration=duration, start=start)


__all__ = [
    "MarkerSoftFileLock",
    "OwnerMode",
    "OwnerRecord",
    "encode_marker",
    "parse_marker",
]
