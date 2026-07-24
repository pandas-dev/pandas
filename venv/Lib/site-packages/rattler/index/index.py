from __future__ import annotations

import datetime
from dataclasses import dataclass
import os
import sys
from typing import TYPE_CHECKING, Any, Literal, Mapping, Optional, TypedDict

if TYPE_CHECKING:
    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

from rattler.platform import Platform
from rattler.rattler import py_index_fs, py_index_s3


@dataclass
class S3Credentials:
    """Credentials for accessing an S3 backend."""

    # The endpoint URL of the S3 backend
    endpoint_url: str

    # The region of the S3 backend
    region: str

    # The access key ID for the S3 bucket.
    access_key_id: Optional[str] = None

    # The secret access key for the S3 bucket.
    secret_access_key: Optional[str] = None

    # The session token for the S3 bucket.
    session_token: Optional[str] = None

    # Defines how to address the bucket, either using virtual-hosted-style or path-style.
    addressing_style: Literal["path", "virtual-host"] = "virtual-host"


class RepodataRevisionMetadata(TypedDict, total=False):
    """Metadata for a single revision in the `vN`-keyed dictionary form of
    `repodata_revisions`. The revision itself is the dictionary key."""

    n_packages: Optional[int]
    oldest: Optional[datetime.datetime]
    newest: Optional[datetime.datetime]


# Repodata revisions to advertise: a `vN`-keyed dictionary mapping each revision
# to its metadata (the CEP shape), e.g. `{"v3": {"n_packages": 1}}`.
RepodataRevisions: TypeAlias = Mapping[str, Optional[RepodataRevisionMetadata]]


def _revision_to_wire(revision: str) -> int:
    if revision == "v3":
        return 3
    raise ValueError(f"unsupported repodata revision {revision!r}, expected 'v3'")


def _epoch_ms(timestamp: datetime.datetime) -> int:
    # repodata stores timestamps as Unix milliseconds.
    return int(round(timestamp.timestamp() * 1000))


def _repodata_revisions_to_dicts(
    revisions: Optional[RepodataRevisions],
) -> Optional[list[dict[str, Any]]]:
    if revisions is None:
        return None

    result = []
    for revision, metadata in revisions.items():
        revision_dict: dict[str, Any] = {"revision": _revision_to_wire(revision)}
        if metadata is not None:
            n_packages = metadata.get("n_packages")
            if n_packages is not None:
                revision_dict["n_packages"] = n_packages
            oldest = metadata.get("oldest")
            if oldest is not None:
                revision_dict["oldest"] = _epoch_ms(oldest)
            newest = metadata.get("newest")
            if newest is not None:
                revision_dict["newest"] = _epoch_ms(newest)
        result.append(revision_dict)
    return result


async def index_fs(
    channel_directory: os.PathLike[str],
    target_platform: Optional[Platform] = None,
    repodata_patch: Optional[str] = None,
    write_zst: bool = True,
    write_shards: bool = True,
    repodata_revisions: Optional[RepodataRevisions] = None,
    package_revision_assignment: Literal["from-index-json", "latest"] = "from-index-json",
    force: bool = False,
    max_parallel: int | None = None,
) -> None:
    """
    Indexes dependencies in the `channel_directory` for one or more subdirectories within said directory.
    Will generate repodata.json files in each subdirectory containing metadata about each present package,
    or if `target_platform` is specified will only consider the subdirectory corresponding to this platform.
    Will always index the "noarch" subdirectory, and thus this subdirectory should always be present, because
    conda channels at a minimum must include this subdirectory.

    Arguments:
        channel_directory: A `os.PathLike[str]` that is the directory containing subdirectories
                           of dependencies to index.
        target_platform: A `Platform` to index dependencies for.
        repodata_patch: The name of the conda package (expected to be in the `noarch` subdir) that should be used for repodata patching.
        write_zst: Whether to write repodata.json.zst.
        write_shards: Whether to write sharded repodata.
        repodata_revisions: Repodata revisions to advertise, with optional `n_packages`, `oldest`, and `newest` metadata.
                            Either a sequence of revisions or a `vN`-keyed dictionary, e.g. `{"v3": {"n_packages": 1}}`.
        package_revision_assignment: Whether to assign packages to the revision required by their `index.json`, or to the latest advertised revision.
        force: Whether to forcefully re-index all subdirs.
        max_parallel: The maximum number of packages to process in-memory simultaneously.
    """
    await py_index_fs(
        channel_directory,
        target_platform._inner if target_platform else target_platform,
        repodata_patch,
        write_zst,
        write_shards,
        _repodata_revisions_to_dicts(repodata_revisions),
        package_revision_assignment,
        force,
        max_parallel,
    )


async def index_s3(
    channel_url: str,
    credentials: Optional[S3Credentials] = None,
    target_platform: Optional[Platform] = None,
    repodata_patch: Optional[str] = None,
    write_zst: bool = True,
    write_shards: bool = True,
    repodata_revisions: Optional[RepodataRevisions] = None,
    package_revision_assignment: Literal["from-index-json", "latest"] = "from-index-json",
    force: bool = False,
    max_parallel: int | None = None,
    precondition_checks: bool = True,
) -> None:
    """
    Indexes dependencies in the `channel_url` for one or more subdirectories in the S3 directory.
    Will generate repodata.json files in each subdirectory containing metadata about each present package,
    or if `target_platform` is specified will only consider the subdirectory corresponding to this platform.
    Will always index the "noarch" subdirectory, and thus this subdirectory should always be present, because
    conda channels at a minimum must include this subdirectory.

    Arguments:
        channel_url: An S3 URL (e.g., s3://my-bucket/my-channel that containins the subdirectories
                     of dependencies to index.
        credentials: The credentials to use for accessing the S3 bucket. If not provided, will use the default
                     credentials from the environment.
        target_platform: A `Platform` to index dependencies for.
        repodata_patch: The name of the conda package (expected to be in the `noarch` subdir) that should be used for repodata patching.
        write_zst: Whether to write repodata.json.zst.
        write_shards: Whether to write sharded repodata.
        repodata_revisions: Repodata revisions to advertise, with optional `n_packages`, `oldest`, and `newest` metadata.
                            Either a sequence of revisions or a `vN`-keyed dictionary, e.g. `{"v3": {"n_packages": 1}}`.
        package_revision_assignment: Whether to assign packages to the revision required by their `index.json`, or to the latest advertised revision.
        force: Whether to forcefully re-index all subdirs.
        max_parallel: The maximum number of packages to process in-memory simultaneously.
        precondition_checks: Whether to perform precondition checks before indexing on S3 buckets which helps to prevent data corruption when indexing with multiple processes at the same time.  Defaults to True.
    """
    await py_index_s3(
        channel_url,
        credentials,
        target_platform._inner if target_platform else target_platform,
        repodata_patch,
        write_zst,
        write_shards,
        _repodata_revisions_to_dicts(repodata_revisions),
        package_revision_assignment,
        force,
        max_parallel,
        precondition_checks,
    )
