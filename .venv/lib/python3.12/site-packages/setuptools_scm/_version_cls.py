from __future__ import annotations

from typing import Type
from typing import Union
from typing import cast

try:
    from packaging.version import InvalidVersion
    from packaging.version import Version as Version
except ImportError:
    from setuptools.extern.packaging.version import (  # type: ignore[import-not-found, no-redef]
        InvalidVersion,
    )
    from setuptools.extern.packaging.version import (  # type: ignore[no-redef]
        Version as Version,
    )
from . import _log

log = _log.log.getChild("version_cls")


class NonNormalizedVersion(Version):
    """A non-normalizing version handler.

    You can use this class to preserve version verification but skip normalization.
    For example you can use this to avoid git release candidate version tags
    ("1.0.0-rc1") to be normalized to "1.0.0rc1". Only use this if you fully
    trust the version tags.
    """

    def __init__(self, version: str) -> None:
        # parse and validate using parent
        super().__init__(version)

        # store raw for str
        self._raw_version = version

    def __str__(self) -> str:
        # return the non-normalized version (parent returns the normalized)
        return self._raw_version

    def __repr__(self) -> str:
        # same pattern as parent
        return f"<NonNormalizedVersion({self._raw_version!r})>"


def _version_as_tuple(version_str: str) -> tuple[int | str, ...]:
    try:
        parsed_version = Version(version_str)
    except InvalidVersion as e:
        log.error("failed to parse version %s: %s", e, version_str)
        return (version_str,)
    else:
        version_fields: tuple[int | str, ...] = parsed_version.release
        if parsed_version.epoch:
            version_fields = (f"{parsed_version.epoch}!", *version_fields)
        if parsed_version.pre is not None:
            version_fields += (f"{parsed_version.pre[0]}{parsed_version.pre[1]}",)

        if parsed_version.post is not None:
            version_fields += (f"post{parsed_version.post}",)

        if parsed_version.dev is not None:
            version_fields += (f"dev{parsed_version.dev}",)

        if parsed_version.local is not None:
            version_fields += (parsed_version.local,)
        return version_fields


_VersionT = Union[Version, NonNormalizedVersion]


def import_name(name: str) -> object:
    import importlib

    pkg_name, cls_name = name.rsplit(".", 1)
    pkg = importlib.import_module(pkg_name)
    return getattr(pkg, cls_name)


def _validate_version_cls(
    version_cls: type[_VersionT] | str | None, normalize: bool
) -> type[_VersionT]:
    if not normalize:
        if version_cls is not None:
            raise ValueError(
                "Providing a custom `version_cls` is not permitted when "
                "`normalize=False`"
            )
        return NonNormalizedVersion
    # Use `version_cls` if provided, default to packaging or pkg_resources
    elif version_cls is None:
        return Version
    elif isinstance(version_cls, str):
        try:
            return cast(Type[_VersionT], import_name(version_cls))
        except Exception:
            raise ValueError(f"Unable to import version_cls='{version_cls}'") from None
    else:
        return version_cls
