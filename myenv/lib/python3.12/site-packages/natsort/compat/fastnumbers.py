# -*- coding: utf-8 -*-
"""
Interface for natsort to access fastnumbers functions without
having to worry if it is actually installed.
"""
import re
from typing import Callable, Iterable, Iterator, Tuple, Union

StrOrFloat = Union[str, float]
StrOrInt = Union[str, int]

__all__ = ["try_float", "try_int"]


def is_supported_fastnumbers(
    fastnumbers_version: str, minimum: Tuple[int, int, int] = (2, 0, 0)
) -> bool:
    match = re.match(
        r"^(\d+)\.(\d+)(\.(\d+))?([ab](\d+))?$",
        fastnumbers_version,
        flags=re.ASCII,
    )

    if not match:
        raise ValueError(
            "Invalid fastnumbers version number '{}'".format(fastnumbers_version)
        )

    (major, minor, patch) = match.group(1, 2, 4)

    return (int(major), int(minor), int(patch)) >= minimum


# If the user has fastnumbers installed, they will get great speed
# benefits. If not, we use the simulated functions that come with natsort.
try:
    # noinspection PyPackageRequirements
    from fastnumbers import fast_float, fast_int, __version__ as fn_ver

    # Require >= version 2.0.0.
    if not is_supported_fastnumbers(fn_ver):
        raise ImportError  # pragma: no cover

    # For versions of fastnumbers with mapping capability, use that
    if is_supported_fastnumbers(fn_ver, (5, 0, 0)):
        del fast_float, fast_int
        from fastnumbers import try_float, try_int
except ImportError:
    from natsort.compat.fake_fastnumbers import fast_float, fast_int  # type: ignore

# Re-map the old-or-compatibility functions fast_float/fast_int to the
# newer API of try_float/try_int. If we already imported try_float/try_int
# then there is nothing to do.
if "try_float" not in globals():

    def try_float(  # type: ignore[no-redef]  # noqa: F811
        x: Iterable[str],
        map: bool,
        nan: float = float("inf"),
        on_fail: Callable[[str], str] = lambda x: x,
    ) -> Iterator[StrOrFloat]:
        assert map is True
        return (fast_float(y, nan=nan, key=on_fail) for y in x)


if "try_int" not in globals():

    def try_int(  # type: ignore[no-redef]  # noqa: F811
        x: Iterable[str],
        map: bool,
        on_fail: Callable[[str], str] = lambda x: x,
    ) -> Iterator[StrOrInt]:
        assert map is True
        return (fast_int(y, key=on_fail) for y in x)
