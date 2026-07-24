"""Bootstrap."""

from __future__ import annotations

import logging
import re
import sys
from operator import eq, lt
from pathlib import Path
from subprocess import PIPE, CalledProcessError, Popen
from typing import TYPE_CHECKING

from .bundle import from_bundle
from .periodic_update import add_wheel_to_update_log
from .util import Version, Wheel, discover_wheels

if TYPE_CHECKING:
    from virtualenv.app_data.base import AppData

LOGGER = logging.getLogger(__name__)

# PEP 503 normalized distribution name. Anything outside this character set on the way to ``pip download`` means
# somebody is smuggling pip options or extras, so reject it before we build the command line.
_DISTRIBUTION_RE = re.compile(
    r"""
    ^
    (?P<name>
        [A-Za-z0-9]                 # must start with an alnum
        (?:[A-Za-z0-9._-]*           # inner chars: alnum plus . _ -
           [A-Za-z0-9])?             # must also end with an alnum (unless length is 1)
    )
    $
    """,
    re.VERBOSE,
)

# Version specifier that matches what ``Version.as_version_spec`` emits: either empty, ``==<ver>`` or ``<<ver>`` where
# ``<ver>`` is a subset of PEP 440 public versions. Kept deliberately strict so a crafted version cannot inject pip
# flags.
_VERSION_SPEC_RE = re.compile(
    r"""
    ^
    (?P<operator>==|<)              # only the operators Version.as_version_spec can emit
    (?P<version>[A-Za-z0-9._+!-]+)  # PEP 440 public-version character set, no whitespace
    $
    """,
    re.VERBOSE,
)


def get_wheel(  # ruff:ignore[too-many-arguments]
    distribution: str,
    version: str | None,
    for_py_version: str,
    search_dirs: list[Path],
    download: bool,
    app_data: AppData,
    do_periodic_update: bool,
    env: dict[str, str],
) -> Wheel | None:
    """Get a wheel with the given distribution-version-for_py_version trio, by using the extra search dir + download."""
    # not all wheels are compatible with all python versions, so we need to py version qualify it
    wheel = None

    if not download or version != Version.bundle:
        # 1. acquire from bundle
        wheel = from_bundle(distribution, version, for_py_version, search_dirs, app_data, do_periodic_update, env)

    if download and wheel is None and version != Version.embed:
        # 2. download from the internet
        wheel = download_wheel(
            distribution=distribution,
            version_spec=Version.as_version_spec(version),
            for_py_version=for_py_version,
            search_dirs=search_dirs,
            app_data=app_data,
            to_folder=app_data.house,
            env=env,
        )
        if wheel is not None and app_data.can_update:
            add_wheel_to_update_log(wheel, for_py_version, app_data)

    return wheel


def download_wheel(  # ruff:ignore[too-many-arguments]
    distribution: str,
    version_spec: str | None,
    for_py_version: str,
    search_dirs: list[Path],
    app_data: AppData,
    to_folder: Path,
    env: dict[str, str],
) -> Wheel:
    """Invoke ``pip download`` in a subprocess to fetch a seed wheel.

    :param distribution: PEP 503 normalized project name; rejected if it contains anything other than
        ``[A-Za-z0-9._-]``.
    :param version_spec: optional version specifier of the form ``==<ver>`` or ``<<ver>`` as emitted by
        :func:`Version.as_version_spec`, or ``None``/empty for the latest compatible release.
    :param for_py_version: major.minor Python version to pass through to ``pip --python-version``.
    :param search_dirs: additional directories to treat as a local wheel index when bootstrapping pip.
    :param app_data: application data store used to locate the embedded pip wheel.
    :param to_folder: directory the downloaded wheel is written into.
    :param env: environment mapping passed through to the subprocess.

    :returns: the downloaded :class:`Wheel`.

    :raises ValueError: if ``distribution`` or ``version_spec`` fail the strict allow-list check.
    :raises CalledProcessError: if ``pip download`` exits with a non-zero status.

    """
    _check_distribution(distribution)
    _check_version_spec(version_spec)
    to_download = f"{distribution}{version_spec or ''}"
    LOGGER.debug("download wheel %s %s to %s", to_download, for_py_version, to_folder)
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "download",
        "--progress-bar",
        "off",
        "--disable-pip-version-check",
        "--only-binary=:all:",
        "--no-deps",
        "--python-version",
        for_py_version,
        "-d",
        str(to_folder),
        to_download,
    ]
    # pip has no interface in python - must be a new sub-process
    env = pip_wheel_env_run(search_dirs, app_data, env)
    process = Popen(cmd, env=env, stdout=PIPE, stderr=PIPE, universal_newlines=True, encoding="utf-8")
    out, err = process.communicate()
    if process.returncode != 0:
        kwargs = {"output": out, "stderr": err}
        raise CalledProcessError(process.returncode, cmd, **kwargs)
    result = _find_downloaded_wheel(distribution, version_spec, for_py_version, to_folder, out)
    LOGGER.debug("downloaded wheel %s", result.name)  # ty: ignore[unresolved-attribute]
    return result  # ty: ignore[invalid-return-type]


def _find_downloaded_wheel(
    distribution: str, version_spec: str | None, for_py_version: str, to_folder: Path, out: str
) -> Wheel | None:
    for line in out.splitlines():
        stripped_line = line.lstrip()
        for marker in ("Saved ", "File was already downloaded "):
            if stripped_line.startswith(marker):
                return Wheel(Path(stripped_line[len(marker) :]).absolute())
    # if for some reason the output does not match fallback to the latest version with that spec
    return find_compatible_in_house(distribution, version_spec, for_py_version, to_folder)


def find_compatible_in_house(
    distribution: str, version_spec: str | None, for_py_version: str, in_folder: Path
) -> Wheel | None:
    wheels = discover_wheels(in_folder, distribution, None, for_py_version)
    start, end = 0, len(wheels)
    if version_spec is not None and version_spec:
        if version_spec.startswith("<"):
            from_pos, op = 1, lt
        elif version_spec.startswith("=="):
            from_pos, op = 2, eq
        else:
            raise ValueError(version_spec)
        version = Wheel.as_version_tuple(version_spec[from_pos:])
        start = next((at for at, w in enumerate(wheels) if op(w.version_tuple, version)), len(wheels))

    return None if start == end else wheels[start]


def pip_wheel_env_run(search_dirs: list[Path], app_data: AppData, env: dict[str, str]) -> dict[str, str]:
    env = env.copy()
    env.update({"PIP_USE_WHEEL": "1", "PIP_USER": "0", "PIP_NO_INPUT": "1", "PYTHONIOENCODING": "utf-8"})
    wheel = get_wheel(
        distribution="pip",
        version=None,
        for_py_version=f"{sys.version_info.major}.{sys.version_info.minor}",
        search_dirs=search_dirs,
        download=False,
        app_data=app_data,
        do_periodic_update=False,
        env=env,
    )
    if wheel is None:
        msg = "could not find the embedded pip"
        raise RuntimeError(msg)
    env["PYTHONPATH"] = str(wheel.path)
    return env


def _check_distribution(distribution: str) -> None:
    if not _DISTRIBUTION_RE.fullmatch(distribution):
        msg = f"refusing to download wheel for suspicious distribution name: {distribution!r}"
        raise ValueError(msg)


def _check_version_spec(version_spec: str | None) -> None:
    if not version_spec:
        return
    if not _VERSION_SPEC_RE.fullmatch(version_spec):
        msg = f"refusing to download wheel with suspicious version spec: {version_spec!r}"
        raise ValueError(msg)


__all__ = [
    "download_wheel",
    "get_wheel",
    "pip_wheel_env_run",
]
