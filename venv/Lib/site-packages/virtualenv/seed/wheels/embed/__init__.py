from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

from virtualenv.info import IS_ZIPAPP, ROOT
from virtualenv.seed.wheels.util import Wheel

BUNDLE_FOLDER = Path(__file__).absolute().parent
BUNDLE_SUPPORT = {
    "3.9": {
        "pip": "pip-26.0.1-py3-none-any.whl",
        "setuptools": "setuptools-82.0.1-py3-none-any.whl",
    },
    "3.10": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.11": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.12": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.13": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.14": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.15": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
    "3.16": {
        "pip": "pip-26.1.2-py3-none-any.whl",
        "setuptools": "setuptools-83.0.0-py3-none-any.whl",
    },
}
MAX = next(reversed(BUNDLE_SUPPORT))
MIN = next(iter(BUNDLE_SUPPORT))


def _release_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


# oldest target Python version virtualenv still bundles seed wheels for; anything below this has no embedded pip
OLDEST_SUPPORTED = _release_tuple(MIN)

# SHA-256 of every bundled wheel. Verified on load so a corrupted or tampered wheel on disk fails loud instead of
# being handed to pip. Generated together with ``BUNDLE_SUPPORT`` by ``tasks/upgrade_wheels.py``.
BUNDLE_SHA256 = {
    "pip-26.0.1-py3-none-any.whl": "bdb1b08f4274833d62c1aa29e20907365a2ceb950410df15fc9521bad440122b",
    "pip-26.1.2-py3-none-any.whl": "382ff9f685ee3bc25864f820aa50505825f10f5458ffff07e30a6d96e5715cab",
    "setuptools-82.0.1-py3-none-any.whl": "a59e362652f08dcd477c78bb6e7bd9d80a7995bc73ce773050228a348ce2e5bb",
    "setuptools-83.0.0-py3-none-any.whl": "29b23c360f22f414dc7336bb39178cc7bcbf6021ed2733cde173f09dba19abb3",
}

_VERIFIED_WHEELS: set[str] = set()


def get_embed_wheel(distribution: str, for_py_version: str | None) -> Wheel | None:
    """Return the bundled wheel that ships with virtualenv for a given distribution and Python version.

    :param distribution: project name of the seed package, for example ``pip`` or ``setuptools``.
    :param for_py_version: major.minor Python version string the environment will be created for, or ``None`` to use the
        newest bundle.

    :returns: a :class:`Wheel` pointing at the verified bundled file, or ``None`` when no wheel is bundled for the
        requested combination, including target versions below the oldest bundled one.

    :raises RuntimeError: if the bundled wheel on disk fails SHA-256 verification.

    """
    if for_py_version is None or _release_tuple(for_py_version) > _release_tuple(MAX):
        # no specific target, or a Python newer than anything bundled: reuse the newest bundle
        mapping = BUNDLE_SUPPORT[MAX]
    else:  # versions below the oldest bundled one fall through to None instead of an incompatible newer wheel
        mapping = BUNDLE_SUPPORT.get(for_py_version)
    if not mapping:
        return None
    wheel_file = mapping.get(distribution)
    if wheel_file is None:
        return None
    path = BUNDLE_FOLDER / wheel_file
    _verify_bundled_wheel(path)
    return Wheel.from_path(path)


def _verify_bundled_wheel(path: Path) -> None:
    name = path.name
    if name in _VERIFIED_WHEELS:
        return
    expected = BUNDLE_SHA256.get(name)
    if expected is None:
        msg = f"bundled wheel {name} has no recorded sha256 in BUNDLE_SHA256"
        raise RuntimeError(msg)
    actual = _hash_bundled_wheel(path)
    if actual != expected:
        msg = f"bundled wheel {name} sha256 mismatch: expected {expected}, got {actual}"
        raise RuntimeError(msg)
    _VERIFIED_WHEELS.add(name)


def _hash_bundled_wheel(path: Path) -> str:
    # ``path`` is under the package directory; when virtualenv runs from a zipapp the wheel lives inside the
    # archive and cannot be opened as a regular file, so read the bytes straight from the zipapp entry.
    digest = hashlib.sha256()
    if IS_ZIPAPP:
        entry = path.resolve().relative_to(Path(ROOT).resolve()).as_posix()
        with zipfile.ZipFile(ROOT, "r") as archive, archive.open(entry) as stream:
            for chunk in iter(lambda: stream.read(1 << 20), b""):
                digest.update(chunk)
    else:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1 << 20), b""):
                digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "BUNDLE_FOLDER",
    "BUNDLE_SHA256",
    "BUNDLE_SUPPORT",
    "MAX",
    "MIN",
    "OLDEST_SUPPORTED",
    "get_embed_wheel",
]
