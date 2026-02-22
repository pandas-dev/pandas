from __future__ import annotations

import dataclasses
import logging
import os
import re
import warnings

from datetime import date
from datetime import datetime
from datetime import timezone
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable
from typing import Match

from . import _entrypoints
from . import _modify_version
from ._node_utils import _format_node_for_output

if TYPE_CHECKING:
    import sys

    if sys.version_info >= (3, 10):
        from typing import Concatenate
        from typing import ParamSpec
    else:
        from typing_extensions import Concatenate
        from typing_extensions import ParamSpec

    _P = ParamSpec("_P")

from typing import TypedDict

from . import _config
from . import _version_cls as _v
from ._version_cls import Version as PkgVersion
from ._version_cls import _VersionT

log = logging.getLogger(__name__)


SEMVER_MINOR = 2
SEMVER_PATCH = 3
SEMVER_LEN = 3


class _TagDict(TypedDict):
    version: str
    prefix: str
    suffix: str


def _parse_version_tag(
    tag: str | object, config: _config.Configuration
) -> _TagDict | None:
    match = config.tag_regex.match(str(tag))

    if match:
        key: str | int = 1 if len(match.groups()) == 1 else "version"
        full = match.group(0)
        log.debug("%r %r %s", tag, config.tag_regex, match)
        log.debug(
            "key %s data %s, %s, %r", key, match.groupdict(), match.groups(), full
        )

        if version := match.group(key):
            result = _TagDict(
                version=version,
                prefix=full[: match.start(key)],
                suffix=full[match.end(key) :],
            )

            log.debug("tag %r parsed to %r", tag, result)
            return result

        raise ValueError(
            f'The tag_regex "{config.tag_regex.pattern}" matched tag "{tag}", '
            "however the matched group has no value."
        )
    else:
        log.debug("tag %r did not parse", tag)

        return None


def callable_or_entrypoint(group: str, callable_or_name: str | Any) -> Any:
    log.debug("ep %r %r", group, callable_or_name)

    if callable(callable_or_name):
        return callable_or_name

    from ._entrypoints import _get_ep

    return _get_ep(group, callable_or_name)


def tag_to_version(
    tag: _VersionT | str, config: _config.Configuration
) -> _VersionT | None:
    """
    take a tag that might be prefixed with a keyword and return only the version part
    """
    log.debug("tag %s", tag)

    tag_dict = _parse_version_tag(tag, config)
    if tag_dict is None or not tag_dict.get("version", None):
        warnings.warn(f"tag {tag!r} no version found")
        return None

    version_str = tag_dict["version"]
    log.debug("version pre parse %s", version_str)

    # Try to create version from base version first
    try:
        version: _VersionT = config.version_cls(version_str)
        log.debug("version=%r", version)
    except Exception:
        warnings.warn(
            f"tag {tag!r} will be stripped of its suffix {tag_dict.get('suffix', '')!r}"
        )
        # Fall back to trying without any suffix
        version = config.version_cls(version_str)
        log.debug("version=%r", version)
        return version

    # If base version is valid, check if we can preserve the suffix
    if suffix := tag_dict.get("suffix", ""):
        log.debug("tag %r includes local build data %r, preserving it", tag, suffix)
        # Try creating version with suffix - if it fails, we'll use the base version
        try:
            version_with_suffix = config.version_cls(version_str + suffix)
            log.debug("version with suffix=%r", version_with_suffix)
            return version_with_suffix
        except Exception:
            warnings.warn(f"tag {tag!r} will be stripped of its suffix {suffix!r}")
            # Return the base version without suffix
            return version

    return version


def _source_epoch_or_utc_now() -> datetime:
    if "SOURCE_DATE_EPOCH" in os.environ:
        date_epoch = int(os.environ["SOURCE_DATE_EPOCH"])
        return datetime.fromtimestamp(date_epoch, timezone.utc)
    else:
        return datetime.now(timezone.utc)


@dataclasses.dataclass
class ScmVersion:
    """represents a parsed version from scm"""

    tag: _v.Version | _v.NonNormalizedVersion
    """the related tag or preformatted version"""
    config: _config.Configuration
    """the configuration used to parse the version"""
    distance: int = 0
    """the number of commits since the tag"""
    node: str | None = None
    """the shortened node id"""
    dirty: bool = False
    """whether the working copy had uncommitted changes"""
    preformatted: bool = False
    """whether the version string was preformatted"""
    branch: str | None = None
    """the branch name if any"""
    node_date: date | None = None
    """the date of the commit if available"""
    time: datetime = dataclasses.field(default_factory=_source_epoch_or_utc_now)
    """the current time or source epoch time
    only set for unit-testing version schemes
    for real usage it must be `now(utc)` or `SOURCE_EPOCH`
    """

    @property
    def exact(self) -> bool:
        """returns true checked out exactly on a tag and no local changes apply"""
        return self.distance == 0 and not self.dirty

    @property
    def short_node(self) -> str | None:
        """Return the node formatted for output."""
        return _format_node_for_output(self.node)

    def __repr__(self) -> str:
        return (
            f"<ScmVersion {self.tag} dist={self.distance} "
            f"node={self.node} dirty={self.dirty} branch={self.branch}>"
        )

    def format_with(self, fmt: str, **kw: object) -> str:
        """format a given format string with attributes of this object"""
        return fmt.format(
            time=self.time,
            tag=self.tag,
            distance=self.distance,
            node=_format_node_for_output(self.node),
            dirty=self.dirty,
            branch=self.branch,
            node_date=self.node_date,
            **kw,
        )

    def format_choice(self, clean_format: str, dirty_format: str, **kw: object) -> str:
        """given `clean_format` and `dirty_format`

        choose one based on `self.dirty` and format it using `self.format_with`"""

        return self.format_with(dirty_format if self.dirty else clean_format, **kw)

    def format_next_version(
        self,
        guess_next: Callable[Concatenate[ScmVersion, _P], str],
        fmt: str = "{guessed}.dev{distance}",
        *k: _P.args,
        **kw: _P.kwargs,
    ) -> str:
        guessed = guess_next(self, *k, **kw)
        return self.format_with(fmt, guessed=guessed)


def _parse_tag(
    tag: _VersionT | str, preformatted: bool, config: _config.Configuration
) -> _VersionT:
    if preformatted:
        # For preformatted versions, tag should already be validated as a version object
        # String validation is handled in meta function before calling this
        if isinstance(tag, str):
            # This should not happen with enhanced meta, but kept for safety
            return _v.NonNormalizedVersion(tag)
        else:
            # Already a version object (including test mocks), return as-is
            return tag
    elif not isinstance(tag, config.version_cls):
        version = tag_to_version(tag, config)
        assert version is not None
        return version
    else:
        return tag


def meta(
    tag: str | _VersionT,
    *,
    distance: int = 0,
    dirty: bool = False,
    node: str | None = None,
    preformatted: bool = False,
    branch: str | None = None,
    config: _config.Configuration,
    node_date: date | None = None,
    time: datetime | None = None,
) -> ScmVersion:
    parsed_version: _VersionT
    # Enhanced string validation for preformatted versions
    if preformatted and isinstance(tag, str):
        # Validate PEP 440 compliance using NonNormalizedVersion
        # Let validation errors bubble up to the caller
        parsed_version = _v.NonNormalizedVersion(tag)
    else:
        # Use existing _parse_tag logic for non-preformatted or already validated inputs
        parsed_version = _parse_tag(tag, preformatted, config)

    log.info("version %s -> %s", tag, parsed_version)
    assert parsed_version is not None, f"Can't parse version {tag}"
    scm_version = ScmVersion(
        parsed_version,
        distance=distance,
        node=node,
        dirty=dirty,
        preformatted=preformatted,
        branch=branch,
        config=config,
        node_date=node_date,
    )
    if time is not None:
        scm_version = dataclasses.replace(scm_version, time=time)
    return scm_version


def guess_next_version(tag_version: ScmVersion) -> str:
    version = _modify_version.strip_local(str(tag_version.tag))
    return _modify_version._bump_dev(version) or _modify_version._bump_regex(version)


def guess_next_dev_version(version: ScmVersion) -> str:
    if version.exact:
        return version.format_with("{tag}")
    else:
        return version.format_next_version(guess_next_version)


def guess_next_simple_semver(
    version: ScmVersion, retain: int, increment: bool = True
) -> str:
    if isinstance(version.tag, _v.Version):
        parts = list(version.tag.release[:retain])
    else:
        try:
            parts = [int(i) for i in str(version.tag).split(".")[:retain]]
        except ValueError:
            raise ValueError(f"{version} can't be parsed as numeric version") from None
    while len(parts) < retain:
        parts.append(0)
    if increment:
        parts[-1] += 1
    while len(parts) < SEMVER_LEN:
        parts.append(0)
    return ".".join(str(i) for i in parts)


def simplified_semver_version(version: ScmVersion) -> str:
    if version.exact:
        return guess_next_simple_semver(version, retain=SEMVER_LEN, increment=False)
    elif version.branch is not None and "feature" in version.branch:
        return version.format_next_version(
            guess_next_simple_semver, retain=SEMVER_MINOR
        )
    else:
        return version.format_next_version(
            guess_next_simple_semver, retain=SEMVER_PATCH
        )


def release_branch_semver_version(version: ScmVersion) -> str:
    if version.exact:
        return version.format_with("{tag}")
    if version.branch is not None:
        # Does the branch name (stripped of namespace) parse as a version?
        branch_ver_data = _parse_version_tag(
            version.branch.split("/")[-1], version.config
        )
        if branch_ver_data is not None:
            branch_ver = branch_ver_data["version"]
            if branch_ver[0] == "v":
                # Allow branches that start with 'v', similar to Version.
                branch_ver = branch_ver[1:]
            # Does the branch version up to the minor part match the tag? If not it
            # might be like, an issue number or something and not a version number, so
            # we only want to use it if it matches.
            tag_ver_up_to_minor = str(version.tag).split(".")[:SEMVER_MINOR]
            branch_ver_up_to_minor = branch_ver.split(".")[:SEMVER_MINOR]
            if branch_ver_up_to_minor == tag_ver_up_to_minor:
                # We're in a release/maintenance branch, next is a patch/rc/beta bump:
                return version.format_next_version(guess_next_version)
    # We're in a development branch, next is a minor bump:
    return version.format_next_version(guess_next_simple_semver, retain=SEMVER_MINOR)


def release_branch_semver(version: ScmVersion) -> str:
    warnings.warn(
        "release_branch_semver is deprecated and will be removed in the future. "
        "Use release_branch_semver_version instead",
        category=DeprecationWarning,
        stacklevel=2,
    )
    return release_branch_semver_version(version)


def only_version(version: ScmVersion) -> str:
    return version.format_with("{tag}")


def no_guess_dev_version(version: ScmVersion) -> str:
    if version.exact:
        return version.format_with("{tag}")
    else:
        return version.format_next_version(_modify_version._dont_guess_next_version)


_DATE_REGEX = re.compile(
    r"""
    ^(?P<date>
        (?P<prefix>[vV]?)
        (?P<year>\d{2}|\d{4})(?:\.\d{1,2}){2})
        (?:\.(?P<patch>\d*))?$
    """,
    re.VERBOSE,
)


def date_ver_match(ver: str) -> Match[str] | None:
    return _DATE_REGEX.match(ver)


def guess_next_date_ver(
    version: ScmVersion,
    node_date: date | None = None,
    date_fmt: str | None = None,
    version_cls: type | None = None,
) -> str:
    """
    same-day -> patch +1
    other-day -> today

    distance is always added as .devX
    """
    match = date_ver_match(str(version.tag))
    if match is None:
        warnings.warn(
            f"{version} does not correspond to a valid versioning date, "
            "assuming legacy version"
        )
        if date_fmt is None:
            date_fmt = "%y.%m.%d"
    else:
        # deduct date format if not provided
        if date_fmt is None:
            date_fmt = "%Y.%m.%d" if len(match.group("year")) == 4 else "%y.%m.%d"
        if prefix := match.group("prefix"):
            if not date_fmt.startswith(prefix):
                date_fmt = prefix + date_fmt

    today = version.time.date()
    head_date = node_date or today
    # compute patch
    if match is None:
        # For legacy non-date tags, always use patch=0 (treat as "other day")
        # Use yesterday to ensure tag_date != head_date
        from datetime import timedelta

        tag_date = head_date - timedelta(days=1)
    else:
        tag_date = (
            datetime.strptime(match.group("date"), date_fmt)
            .replace(tzinfo=timezone.utc)
            .date()
        )
    if tag_date == head_date:
        assert match is not None
        # Same day as existing date tag - increment patch
        patch = int(match.group("patch") or "0") + 1
    else:
        # Different day or legacy non-date tag - use patch 0
        if tag_date > head_date and match is not None:
            # warn on future times (only for actual date tags, not legacy)
            warnings.warn(
                f"your previous tag  ({tag_date}) is ahead your node date ({head_date})"
            )
        patch = 0
    next_version = "{node_date:{date_fmt}}.{patch}".format(
        node_date=head_date, date_fmt=date_fmt, patch=patch
    )
    # rely on the Version object to ensure consistency (e.g. remove leading 0s)
    if version_cls is None:
        version_cls = PkgVersion
    next_version = str(version_cls(next_version))
    return next_version


def calver_by_date(version: ScmVersion) -> str:
    if version.exact and not version.dirty:
        return version.format_with("{tag}")
    # TODO: move the release-X check to a new scheme
    if version.branch is not None and version.branch.startswith("release-"):
        branch_ver = _parse_version_tag(version.branch.split("-")[-1], version.config)
        if branch_ver is not None:
            ver = branch_ver["version"]
            match = date_ver_match(ver)
            if match:
                return ver
    return version.format_next_version(
        guess_next_date_ver,
        node_date=version.node_date,
        version_cls=version.config.version_cls,
    )


def get_local_node_and_date(version: ScmVersion) -> str:
    return _modify_version._format_local_with_time(version, time_format="%Y%m%d")


def get_local_node_and_timestamp(version: ScmVersion) -> str:
    return _modify_version._format_local_with_time(version, time_format="%Y%m%d%H%M%S")


def get_local_dirty_tag(version: ScmVersion) -> str:
    return version.format_choice("", "+dirty")


def get_no_local_node(version: ScmVersion) -> str:
    return ""


def postrelease_version(version: ScmVersion) -> str:
    if version.exact:
        return version.format_with("{tag}")
    else:
        return version.format_with("{tag}.post{distance}")


def _combine_version_with_local_parts(
    main_version: str, *local_parts: str | None
) -> str:
    """
    Combine a main version with multiple local parts into a valid PEP 440 version string.
    Handles deduplication of local parts to avoid adding the same local data twice.

    Args:
        main_version: The main version string (e.g., "1.2.0", "1.2.dev3")
        *local_parts: Variable number of local version parts, can be None or empty

    Returns:
        A valid PEP 440 version string

    Examples:
        _combine_version_with_local_parts("1.2.0", "build.123", "d20090213") -> "1.2.0+build.123.d20090213"
        _combine_version_with_local_parts("1.2.0", "build.123", None) -> "1.2.0+build.123"
        _combine_version_with_local_parts("1.2.0+build.123", "d20090213") -> "1.2.0+build.123.d20090213"
        _combine_version_with_local_parts("1.2.0+build.123", "build.123") -> "1.2.0+build.123"  # no duplication
        _combine_version_with_local_parts("1.2.0", None, None) -> "1.2.0"
    """
    # Split main version into base and existing local parts
    if "+" in main_version:
        main_part, existing_local = main_version.split("+", 1)
        all_local_parts = existing_local.split(".")
    else:
        main_part = main_version
        all_local_parts = []

    # Process each new local part
    for part in local_parts:
        if not part or not part.strip():
            continue

        # Strip any leading + and split into segments
        clean_part = part.strip("+")
        if not clean_part:
            continue

        # Split multi-part local identifiers (e.g., "build.123" -> ["build", "123"])
        part_segments = clean_part.split(".")

        # Add each segment if not already present
        for segment in part_segments:
            if segment and segment not in all_local_parts:
                all_local_parts.append(segment)

    # Return combined result
    if all_local_parts:
        return main_part + "+" + ".".join(all_local_parts)
    else:
        return main_part


def format_version(version: ScmVersion) -> str:
    log.debug("scm version %s", version)
    log.debug("config %s", version.config)
    if version.preformatted:
        return str(version.tag)

    # Extract original tag's local data for later combination
    original_local = ""
    if hasattr(version.tag, "local") and version.tag.local is not None:
        original_local = str(version.tag.local)

    # Create a patched ScmVersion with only the base version (no local data) for version schemes
    from dataclasses import replace

    # Extract the base version (public part) from the tag using config's version_cls
    base_version_str = str(version.tag.public)
    base_tag = version.config.version_cls(base_version_str)
    version_for_scheme = replace(version, tag=base_tag)

    main_version = _entrypoints._call_version_scheme(
        version_for_scheme,
        "setuptools_scm.version_scheme",
        version.config.version_scheme,
    )
    log.debug("version %s", main_version)
    assert main_version is not None

    local_version = _entrypoints._call_version_scheme(
        version, "setuptools_scm.local_scheme", version.config.local_scheme, "+unknown"
    )
    log.debug("local_version %s", local_version)

    # Combine main version with original local data and new local scheme data
    return _combine_version_with_local_parts(
        str(main_version), original_local, local_version
    )
