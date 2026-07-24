from __future__ import annotations

from typing import TYPE_CHECKING

from virtualenv.seed.wheels.embed import get_embed_wheel

from .periodic_update import periodic_update
from .util import Version, Wheel, discover_wheels

if TYPE_CHECKING:
    from pathlib import Path

    from virtualenv.app_data.base import AppData


def from_bundle(  # ruff:ignore[too-many-arguments]
    distribution: str,
    version: str | None,
    for_py_version: str,
    search_dirs: list[Path],
    app_data: AppData,
    do_periodic_update: bool,
    env: dict[str, str],
) -> Wheel | None:
    """Load the bundled wheel to a cache directory."""
    of_version = Version.of_version(version)
    wheel = load_embed_wheel(app_data, distribution, for_py_version, of_version)

    if version != Version.embed:
        # 2. check if we have upgraded embed
        if app_data.can_update:
            per = do_periodic_update
            wheel = periodic_update(distribution, of_version, for_py_version, wheel, search_dirs, app_data, per, env)

        # 3. acquire from extra search dir
        found_wheel = from_dir(distribution, of_version, for_py_version, search_dirs)
        if found_wheel is not None and (wheel is None or found_wheel.version_tuple > wheel.version_tuple):
            wheel = found_wheel
    return wheel


def load_embed_wheel(app_data: AppData, distribution: str, for_py_version: str, version: str | None) -> Wheel | None:
    wheel = get_embed_wheel(distribution, for_py_version)
    if wheel is not None:
        version_match = version == wheel.version
        if version is None or version_match:
            with app_data.ensure_extracted(wheel.path, lambda: app_data.house) as wheel_path:  # ty: ignore[invalid-argument-type]
                wheel = Wheel(wheel_path)
        else:  # if version does not match ignore
            wheel = None
    return wheel


def from_dir(distribution: str, version: str | None, for_py_version: str, directories: list[Path]) -> Wheel | None:
    """Load a compatible wheel from a given folder."""
    for folder in directories:
        for wheel in discover_wheels(folder, distribution, version, for_py_version):
            return wheel
    return None


__all__ = [
    "from_bundle",
    "load_embed_wheel",
]
