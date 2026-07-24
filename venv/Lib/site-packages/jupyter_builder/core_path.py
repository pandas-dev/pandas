"""Utilities for locating and resolving JupyterLab core package metadata."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import io
import json
import logging
import os
import re
import subprocess
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

from .constants import JPBLD_NPM_URL, JPBLD_RAW_GITHUB_URL

_MAX_CORE_META_BYTES = 5 * 1024 * 1024  # 5 MB — generous upper bound for core.package.json

#: GitHub API endpoint used to resolve wildcard versions to concrete release tags
_GITHUB_TAGS_API_URL = "https://api.github.com/repos/jupyterlab/jupyterlab/tags"

#: Upper bound on tag-list pages fetched when resolving a wildcard (100 tags per page)
_MAX_GITHUB_TAG_PAGES = 10


def _home_dir() -> Path:
    home = os.environ.get("HOME")
    return Path(home) if home else Path.home()


def _http_get(url: str, *, headers: dict[str, str] | None = None, timeout: int = 10) -> bytes:
    if not url.startswith(("http:", "https:")):
        msg = "URL must start with 'http:' or 'https:'"
        raise ValueError(msg)
    request_headers = {"User-Agent": "jupyter-builder"}
    if headers:
        request_headers.update(headers)
    req = urllib.request.Request(url, headers=request_headers)  # noqa: S310
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        return bytes(resp.read())


def get_core_meta(
    version: str | None = None,
    ext_path: str | os.PathLike[str] | None = None,
    logger: logging.Logger | None = None,
) -> str:
    """Return the path to the core package JSON, downloading it if needed."""
    requested_version = version
    if requested_version is None:
        if ext_path is not None:
            installed_core_meta = _get_installed_core_meta(Path(ext_path).resolve())
            if installed_core_meta is not None:
                return installed_core_meta
            if logger:
                logger.warning(
                    "\033[33m@jupyterlab/core-meta was not found in node_modules, "
                    "so a network download will be used as a fallback. To avoid this, "
                    "add @jupyter/builder as a devDependency instead of "
                    "@jupyterlab/builder.\n \033[0m",
                )
        requested_version = "main"
    else:
        # Accept both "vX.Y.Z" and "X.Y.Z" for an explicitly requested version.
        requested_version = _normalize_version(requested_version)

    cache_root = _home_dir() / ".cache" / "jupyterlab_builder" / "core"
    cached_file = _get_cached_core_meta_file(cache_root, requested_version)
    if cached_file is not None:
        return str(cached_file)

    # Try to retrieve core meta from npm first, then fall back to GitHub. If the
    # requested version cannot be found in either source, raise an error.
    try:
        npm_version = _resolve_npm_version(requested_version)
        npm_cache_file = cache_root / npm_version / "core.package.json"
        if npm_cache_file.exists():
            return str(npm_cache_file)
        _download_npm_core_meta(npm_version, npm_cache_file)
        return str(npm_cache_file)
    except urllib.error.URLError as npm_error:
        try:
            # Wildcards (e.g. "4.5.x") have no single git ref, so resolve them
            # to a concrete tag from the GitHub tag list before downloading.
            github_version = (
                _resolve_wildcard_github_version(requested_version)
                if _is_wildcard_version(requested_version)
                else requested_version
            )
            github_cache_file = cache_root / github_version / "core.package.json"
            if github_cache_file.exists():
                return str(github_cache_file)
            _download_github_core_meta(_github_ref(github_version), github_cache_file)
        except urllib.error.URLError as github_error:
            msg = (
                f"Could not resolve @jupyterlab/core-meta for requested version "
                f"{requested_version!r}: not found on the npm registry "
                f"({npm_error}) or in the jupyterlab/jupyterlab GitHub repository "
                f"({github_error}). Verify that the version exists "
                f"(both '4.5.7' and 'v4.5.7' are accepted)."
            )
            raise RuntimeError(msg) from github_error
        return str(github_cache_file)


def _normalize_version(version: str) -> str:
    """Strip a leading 'v' from a numeric version so 'vX.Y.Z' and 'X.Y.Z' are equivalent."""
    return version[1:] if re.match(r"v\d", version) else version


def _github_ref(version: str) -> str:
    """Map a resolved version to its jupyterlab/jupyterlab git ref.

    Numeric releases are published as git tags. Stable releases are tagged like
    'v4.5.7', while npm-style prereleases such as '4.6.0-alpha.4' correspond to
    the PEP 440 tag form JupyterLab uses, 'v4.6.0a4'. Branch names such as
    'main' (and other non-numeric refs) are used verbatim.
    """
    if not re.match(r"\d", version):
        return version
    release, separator, prerelease = version.partition("-")
    if separator:
        # Translate npm prerelease identifiers (alpha.4 / beta.1 / rc.2) to the
        # PEP 440 form JupyterLab tags use (a4 / b1 / rc2).
        prerelease = re.sub(r"alpha\.?", "a", prerelease)
        prerelease = re.sub(r"beta\.?", "b", prerelease)
        prerelease = re.sub(r"rc\.", "rc", prerelease)
        version = release + prerelease
    return f"v{version}"


def _is_wildcard_version(version: str) -> bool:
    """Return True for npm range-style versions like 4.5.x."""
    return bool(re.search(r"\.x(\.|$)|(^|\.)x\.", version, flags=re.IGNORECASE))


def _resolve_npm_version(version: str) -> str:
    """Resolve an abstract version specifier to a concrete npm version string.

    - 'latest'  → fetches the current latest tag from npm
    - '4.5.x'   → fetches all published versions and returns the highest 4.5.x match
    - anything else is returned as-is (assumed to be a concrete version)
    """
    if version == "latest":
        data = _http_get(f"{JPBLD_NPM_URL}/@jupyterlab/core-meta/latest")
        latest_version = json.loads(data).get("version")
        if not isinstance(latest_version, str) or not latest_version:
            msg = "Failed to resolve latest @jupyterlab/core-meta version from npm"
            raise urllib.error.URLError(msg)
        return latest_version

    if _is_wildcard_version(version):
        return _resolve_wildcard_npm_version(version)

    return version  # Concrete version like "4.2.5" — use directly


def _resolve_wildcard_npm_version(version: str) -> str:
    """Fetch the highest published npm version matching a wildcard range like '4.5.x'.

    Raises urllib.error.URLError if no matching version is found.
    """
    data = _http_get(
        f"{JPBLD_NPM_URL}/@jupyterlab/core-meta",
        headers={"Accept": "application/vnd.npm.install-v1+json"},
    )
    all_versions: list[str] = list(json.loads(data).get("versions", {}).keys())

    # Build a regex from the wildcard pattern, e.g. "4.5.x" → r"^4\.5\.\d+(-.+)?$"
    # The (-.+)? suffix matches pre-release identifiers like -alpha.3, -beta.1, -rc.2
    escaped = re.escape(version)
    wildcard_pattern = re.sub(r"x", r"\\d+", escaped, flags=re.IGNORECASE)
    pattern = "^" + wildcard_pattern + r"(-.+)?$"

    matching = [v for v in all_versions if re.match(pattern, v)]
    if not matching:
        msg = f"No published @jupyterlab/core-meta versions match range '{version}'"
        raise urllib.error.URLError(msg)

    return max(matching, key=_semver_key)


def _semver_key(v: str) -> tuple[tuple[int, ...], int, tuple[int, ...]]:
    release, _, prerelease = v.partition("-")
    numeric = tuple(int(p) for p in release.split(".") if p.isdigit())
    # Stable releases sort higher than pre-releases; within pre-releases,
    # order by the numeric identifiers (e.g. alpha.3 < alpha.4).
    pre_numeric = tuple(int(p) for p in prerelease.split(".") if p.isdigit())
    return (numeric, 0 if prerelease else 1, pre_numeric)


def _resolve_wildcard_github_version(version: str) -> str:
    """Resolve a wildcard range like '4.5.x' to the highest matching git tag.

    JupyterLab publishes stable releases as git tags (e.g. 'v4.5.9') that may
    predate @jupyterlab/core-meta on npm, so wildcards that npm cannot satisfy
    are resolved here against the jupyterlab/jupyterlab tag list.

    Raises urllib.error.URLError if no matching tag is found.
    """
    escaped = re.escape(version)
    wildcard_pattern = re.sub(r"x", r"\\d+", escaped, flags=re.IGNORECASE)
    pattern = re.compile("^" + wildcard_pattern + r"$")

    matching: list[str] = []
    for page in range(1, _MAX_GITHUB_TAG_PAGES + 1):
        data = _http_get(f"{_GITHUB_TAGS_API_URL}?per_page=100&page={page}")
        tags = json.loads(data)
        if not tags:
            break
        page_matches = [
            normalized
            for tag in tags
            if (normalized := _normalize_version(tag.get("name", ""))) and pattern.match(normalized)
        ]
        # Tags are returned newest-first, so once matches stop appearing
        # (after they have started) the rest are older and can be skipped.
        if matching and not page_matches:
            break
        matching.extend(page_matches)

    if not matching:
        msg = f"No jupyterlab/jupyterlab git tags match range '{version}'"
        raise urllib.error.URLError(msg)

    return max(matching, key=_semver_key)


def _get_cached_core_meta_file(cache_root: Path, version: str) -> Path | None:
    candidates = [
        cache_root / version / "core.package.json",
        cache_root / version / "package.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _download_npm_core_meta(version: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(_http_get(f"{JPBLD_NPM_URL}/@jupyterlab/core-meta/{version}"))
    try:
        tarball_url = metadata["dist"]["tarball"]
    except (KeyError, TypeError) as exc:
        msg = f"Unexpected registry metadata for {version}: {exc}"
        raise urllib.error.URLError(msg) from exc
    tarball_data = _http_get(tarball_url, timeout=20)
    with tarfile.open(fileobj=io.BytesIO(tarball_data), mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name == "package/core.package.json":
                if not member.isfile():
                    msg = "core.package.json entry is not a regular file"
                    raise urllib.error.URLError(msg)
                if member.size > _MAX_CORE_META_BYTES:
                    msg = f"core.package.json entry exceeds size limit ({member.size} bytes)"
                    raise urllib.error.URLError(msg)
                f = tar.extractfile(member)
                if f:
                    destination.write_bytes(f.read(_MAX_CORE_META_BYTES))
                    return
    msg = f"core.package.json not found in tarball for {version}"
    raise urllib.error.URLError(msg)


def _download_github_core_meta(version: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = f"{JPBLD_RAW_GITHUB_URL}/jupyterlab/jupyterlab/{version}/jupyterlab/staging/package.json"
    destination.write_bytes(_http_get(url))


def _get_installed_core_meta(ext_path: Path) -> str | None:
    if not (ext_path / "node_modules").exists():
        subprocess.check_call(["jlpm"], cwd=ext_path)  # noqa: S607

    target = ext_path
    while True:
        core_meta_path = target / "node_modules" / "@jupyterlab" / "core-meta"
        if (core_meta_path / "core.package.json").exists():
            return str(core_meta_path / "core.package.json")
        if target.parent == target:
            return None
        target = target.parent
