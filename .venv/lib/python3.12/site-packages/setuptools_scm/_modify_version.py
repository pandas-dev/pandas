from __future__ import annotations

import re

from . import _types as _t


def strip_local(version_string: str) -> str:
    public = version_string.partition("+")[0]
    return public


def _add_post(version: str) -> str:
    if "post" in version:
        raise ValueError(
            f"{version} already is a post release, refusing to guess the update"
        )
    return f"{version}.post1"


def _bump_dev(version: str) -> str | None:
    if ".dev" not in version:
        return None

    prefix, tail = version.rsplit(".dev", 1)
    if tail != "0":
        raise ValueError(
            "choosing custom numbers for the `.devX` distance "
            "is not supported.\n "
            f"The {version} can't be bumped\n"
            "Please drop the tag or create a new supported one ending in .dev0"
        )
    return prefix


def _bump_regex(version: str) -> str:
    match = re.match(r"(.*?)(\d+)$", version)
    if match is None:
        raise ValueError(
            f"{version} does not end with a number to bump, "
            "please correct or use a custom version scheme"
        )
    else:
        prefix, tail = match.groups()
        return f"{prefix}{int(tail) + 1}"


def _format_local_with_time(version: _t.SCMVERSION, time_format: str) -> str:
    if version.exact or version.node is None:
        return version.format_choice(
            "", "+d{time:{time_format}}", time_format=time_format
        )
    else:
        return version.format_choice(
            "+{node}", "+{node}.d{time:{time_format}}", time_format=time_format
        )


def _dont_guess_next_version(tag_version: _t.SCMVERSION) -> str:
    version = strip_local(str(tag_version.tag))
    return _bump_dev(version) or _add_post(version)
