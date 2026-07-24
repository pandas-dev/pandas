# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Extension manager using pip as package manager and PyPi.org as packages source."""

import asyncio
import http.client
import io
import json
import math
import re
import sys
import tempfile
import xmlrpc.client
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from functools import partial
from itertools import groupby
from os import environ
from pathlib import Path
from subprocess import CalledProcessError, run
from tarfile import TarFile
from typing import Any, Optional
from urllib.parse import urlparse
from zipfile import ZipFile

import httpx
import tornado
from async_lru import alru_cache
from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.utils import InvalidName, canonicalize_name
from packaging.version import InvalidVersion, Version
from packaging.version import parse as parse_version
from traitlets import CFloat, CInt, Unicode, config, observe

try:
    from typing import override
except ImportError:
    from typing_extensions import override

from jupyterlab._version import __version__
from jupyterlab.extensions.manager import (
    ActionResult,
    ExtensionManager,
    ExtensionManagerMetadata,
    ExtensionPackage,
)


class ProxiedTransport(xmlrpc.client.Transport):
    def set_proxy(self, host, port=None, headers=None):
        self.proxy = host, port
        self.proxy_headers = headers

    def make_connection(self, host):
        connection = http.client.HTTPConnection(*self.proxy)
        connection.set_tunnel(host, headers=self.proxy_headers)
        self._connection = host, connection
        return connection


all_proxy_url = environ.get("ALL_PROXY")
# For historical reasons, we also support the lowercase environment variables.
# Info: https://about.gitlab.com/blog/2021/01/27/we-need-to-talk-no-proxy/
http_proxy_url = environ.get("http_proxy") or environ.get("HTTP_PROXY") or all_proxy_url
https_proxy_url = (
    environ.get("https_proxy") or environ.get("HTTPS_PROXY") or http_proxy_url or all_proxy_url
)

# sniff ``httpx`` version for version-sensitive API
_httpx_version = Version(httpx.__version__)
_httpx_client_args = {}

xmlrpc_transport_override = None

if http_proxy_url:
    http_proxy = urlparse(http_proxy_url)
    proxy_host, _, proxy_port = http_proxy.netloc.partition(":")

    if _httpx_version >= Version("0.28.0"):
        _httpx_client_args = {
            "mounts": {
                "http://": httpx.AsyncHTTPTransport(proxy=http_proxy_url),
                "https://": httpx.AsyncHTTPTransport(proxy=https_proxy_url),
            }
        }
    else:
        _httpx_client_args = {
            "proxies": {
                "http://": http_proxy_url,
                "https://": https_proxy_url,
            }
        }

    xmlrpc_transport_override = ProxiedTransport()
    xmlrpc_transport_override.set_proxy(proxy_host, proxy_port)


def _check_python_version_compatible(requires_python: str | None) -> tuple[bool, str | None]:
    """Check if the current Python version satisfies the requires_python specifier.

    Args:
        requires_python: The requires_python specifier string from PyPI (e.g., ">=3.10")

    Returns:
        (compatible, explanation)
        compatible: True if compatible or if requires_python is None/empty, False otherwise.
        explanation: A string explaining the mismatch if incompatible, otherwise None.
    """
    if not requires_python:
        return True, None
    try:
        current_version_str = sys.version.split()[0]
        current_version = Version(current_version_str)
        specifier = SpecifierSet(requires_python)
        if current_version in specifier:
            return True, None
        return False, f"Requires Python {requires_python} but detected Python {current_version}"
    except (InvalidSpecifier, InvalidVersion):
        # If parsing fails, assume compatible to avoid false negatives
        return True, None


async def _fetch_package_metadata(
    client: httpx.AsyncClient,
    name: str,
    latest_version: str,
    base_url: str,
) -> dict:
    response = await client.get(
        base_url + f"/{name}/{latest_version}/json",
        headers={"Content-Type": "application/json"},
    )
    if response.status_code < 400:  # noqa PLR2004
        data = json.loads(response.text).get("info")

        # Keep minimal information to limit cache size
        return {
            k: data.get(k)
            for k in [
                "author",
                "bugtrack_url",
                "docs_url",
                "home_page",
                "license",
                "package_url",
                "project_url",
                "project_urls",
                "requires_python",
                "summary",
            ]
        }
    else:
        return {}


# Known language packs from https://github.com/jupyterlab/language-packs
# These are not tagged with the prebuilt extension classifier on PyPI,
# so they are listed explicitly.
LANGUAGE_PACKS = (
    "jupyterlab-language-pack-ar-SA",
    "jupyterlab-language-pack-ca-ES",
    "jupyterlab-language-pack-cs-CZ",
    "jupyterlab-language-pack-da-DK",
    "jupyterlab-language-pack-de-DE",
    "jupyterlab-language-pack-el-GR",
    "jupyterlab-language-pack-es-ES",
    "jupyterlab-language-pack-et-EE",
    "jupyterlab-language-pack-fi-FI",
    "jupyterlab-language-pack-fr-FR",
    "jupyterlab-language-pack-he-IL",
    "jupyterlab-language-pack-hu-HU",
    "jupyterlab-language-pack-hy-AM",
    "jupyterlab-language-pack-id-ID",
    "jupyterlab-language-pack-it-IT",
    "jupyterlab-language-pack-ja-JP",
    "jupyterlab-language-pack-ko-KR",
    "jupyterlab-language-pack-lt-LT",
    "jupyterlab-language-pack-nl-NL",
    "jupyterlab-language-pack-no-NO",
    "jupyterlab-language-pack-pl-PL",
    "jupyterlab-language-pack-pt-BR",
    "jupyterlab-language-pack-ro-RO",
    "jupyterlab-language-pack-ru-RU",
    "jupyterlab-language-pack-tr-TR",
    "jupyterlab-language-pack-uk-UA",
    "jupyterlab-language-pack-vi-VN",
    "jupyterlab-language-pack-zh-CN",
    "jupyterlab-language-pack-zh-TW",
)


class PyPIExtensionManager(ExtensionManager):
    """Extension manager using pip as package manager and PyPi.org as packages source."""

    base_url = Unicode("https://pypi.org/pypi", config=True, help="The base URL of PyPI index.")

    cache_timeout = CFloat(
        5 * 60.0, config=True, help="PyPI extensions list cache timeout in seconds."
    )

    package_metadata_cache_size = CInt(
        1500, config=True, help="The cache size for package metadata."
    )

    rpc_request_throttling = CFloat(
        1.0,
        config=True,
        help="Throttling time in seconds between PyPI requests using the XML-RPC API.",
    )

    def __init__(
        self,
        app_options: dict | None = None,
        ext_options: dict | None = None,
        parent: config.Configurable | None = None,
    ) -> None:
        super().__init__(app_options, ext_options, parent)
        self._httpx_client = httpx.AsyncClient(**_httpx_client_args)
        # Set configurable cache size to fetch function
        self._fetch_package_metadata = partial(_fetch_package_metadata, self._httpx_client)
        self._observe_package_metadata_cache_size({"new": self.package_metadata_cache_size})
        # Combine XML RPC API and JSON API to reduce throttling by PyPI.org
        self._rpc_client = xmlrpc.client.ServerProxy(
            self.base_url, transport=xmlrpc_transport_override
        )
        self.__last_all_packages_request_time = datetime.now(tz=timezone.utc) - timedelta(
            seconds=self.cache_timeout * 1.01
        )
        self.__all_packages_cache = None

        self.log.debug(f"Extensions list will be fetched from {self.base_url}.")
        if xmlrpc_transport_override:
            self.log.info(
                f"Extensions will be fetched using proxy, proxy host and port: {xmlrpc_transport_override.proxy}"
            )

    @property
    def metadata(self) -> ExtensionManagerMetadata:
        """Extension manager metadata."""
        return ExtensionManagerMetadata("PyPI", True, sys.prefix)

    @override
    async def is_install_allowed(self, name: str, version: str | None = None) -> bool:
        try:
            canonicalize_name(name, validate=True)
            if version is not None:
                parse_version(version)
        except InvalidName:
            self.log.warning(f"Invalid extension name: {name!r}")
            return False
        except InvalidVersion:
            self.log.warning(f"Version {version!r} does not comply with PEP 440")
            return False

        allowed = await self._is_allowed_by_listing(name)
        if not allowed:
            self.log.warning(f"Installation denied by allowlist/blocklist for {name}")
        return allowed

    @override
    def _canonicalize_name(self, name: str) -> str:
        """Canonicalize PyPI package names for listing policy comparisons."""
        return canonicalize_name(name)

    async def get_latest_version(self, pkg: str) -> str | None:
        """Return the latest available version for a given extension.

        Args:
            pkg: The extension to search for
        Returns:
            The latest available version
        """
        try:
            response = await self._httpx_client.get(
                self.base_url + f"/{pkg}/json", headers={"Content-Type": "application/json"}
            )

            if response.status_code < 400:  # noqa PLR2004
                data = json.loads(response.content).get("info", {})
            else:
                self.log.debug(f"Failed to get package information on PyPI; {response!s}")
                return None
        except Exception:
            return None
        else:
            return ExtensionManager.get_semver_version(data.get("version", "")) or None

    def get_normalized_name(self, extension: ExtensionPackage) -> str:
        """Normalize extension name.

        Extension have multiple parts, npm package, Python package,...
        Sub-classes may override this method to ensure the name of
        an extension from the service provider and the local installed
        listing is matching.

        Args:
            extension: The extension metadata
        Returns:
            The normalized name
        """
        if extension.install is not None:
            install_metadata = extension.install
            if install_metadata["packageManager"] == "python":
                return self._normalize_name(install_metadata["packageName"])
        return self._normalize_name(extension.name)

    async def __throttleRequest(self, recursive: bool, fn: Callable, *args) -> Any:  # noqa
        """Throttle XMLRPC API request

        Args:
            recursive: Whether to call the throttling recursively once or not.
            fn: API method to call
            *args: API method arguments
        Returns:
            Result of the method
        Raises:
            xmlrpc.client.Fault
        """
        current_loop = tornado.ioloop.IOLoop.current()
        try:
            data = await current_loop.run_in_executor(None, fn, *args)
        except xmlrpc.client.Fault as err:
            if err.faultCode == -32500 and err.faultString.startswith(  # noqa PLR2004
                "HTTPTooManyRequests:"
            ):
                delay = 1.01
                match = re.search(r"Limit may reset in (\d+) seconds.", err.faultString)
                if match is not None:
                    delay = int(match.group(1) or "1")
                self.log.info(
                    f"HTTPTooManyRequests - Perform next call to PyPI XMLRPC API in {delay}s."
                )
                await asyncio.sleep(delay * self.rpc_request_throttling + 0.01)
                if recursive:
                    data = await self.__throttleRequest(False, fn, *args)
                else:
                    data = await current_loop.run_in_executor(None, fn, *args)

        return data

    @observe("package_metadata_cache_size")
    def _observe_package_metadata_cache_size(self, change):
        self._fetch_package_metadata = alru_cache(maxsize=change["new"])(
            partial(_fetch_package_metadata, self._httpx_client)
        )

    async def list_packages(
        self, query: str, page: int, per_page: int
    ) -> tuple[dict[str, ExtensionPackage], int | None]:
        """List the available extensions.

        Note:
            This will list the packages based on the classifier
                Framework :: Jupyter :: JupyterLab :: Extensions :: Prebuilt

            Then it filters it with the query and sorts by organization priority:
            1. Project Jupyter (@jupyter)
            2. JupyterLab Community (@jupyterlab-contrib)
            3. Others

        Args:
            query: The search extension query
            page: The result page
            per_page: The number of results per page
        Returns:
            The available extensions in a mapping {name: metadata}
            The results last page; None if the manager does not support pagination
        """
        matches = await self.__get_all_extensions()

        extensions = {}
        all_matches = []

        for name, group in groupby(filter(lambda m: query in m[0], matches), lambda e: e[0]):
            _, latest_version = list(group)[-1]
            data = await self._fetch_package_metadata(name, latest_version, self.base_url)

            normalized_name = self._normalize_name(name)
            package_urls = data.get("project_urls") or {}

            source_url = package_urls.get("Source Code")
            homepage_url = data.get("home_page") or package_urls.get("Homepage")
            documentation_url = data.get("docs_url") or package_urls.get("Documentation")
            bug_tracker_url = data.get("bugtrack_url") or package_urls.get("Bug Tracker")

            best_guess_home_url = (
                homepage_url
                or data.get("project_url")
                or data.get("package_url")
                or documentation_url
                or source_url
                or bug_tracker_url
            )

            # Check Python version compatibility
            requires_python = data.get("requires_python")
            python_compatible, version_explanation = _check_python_version_compatible(
                requires_python
            )

            description = data.get("summary")
            if version_explanation:
                if description:
                    description += f" ({version_explanation})"
                else:
                    description = version_explanation

            extension = ExtensionPackage(
                name=normalized_name,
                description=description,
                homepage_url=best_guess_home_url,
                author=data.get("author"),
                license=data.get("license"),
                latest_version=ExtensionManager.get_semver_version(latest_version),
                pkg_type="prebuilt",
                allowed=python_compatible,
                bug_tracker_url=bug_tracker_url,
                documentation_url=documentation_url,
                package_manager_url=data.get("package_url"),
                repository_url=source_url,
            )

            # Determine organization priority
            priority = 3  # Default priority for other packages
            urls_to_check = [
                str(url).lower() for url in [source_url, homepage_url, best_guess_home_url] if url
            ]
            exclude = [
                "https://github.com/jupyterlab/jupyterlab_apod",
                "https://github.com/jupyterlab/extension-examples",
            ]

            for url in urls_to_check:
                if url in exclude:
                    priority = 4
                    break
                if any(
                    org in url
                    for org in ["github.com/jupyter/", "jupyter.org", "github.com/jupyterlab/"]
                ):
                    priority = 1
                    break
                elif "github.com/jupyterlab-contrib/" in url:
                    priority = 2
                    break

            all_matches.append((priority, extension))

        sorted_matches = sorted(all_matches, key=lambda x: (x[0], x[1].name))

        # Apply pagination
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        page_matches = sorted_matches[start_idx:end_idx]

        for _, extension in page_matches:
            extensions[extension.name] = extension

        total_pages = math.ceil(len(sorted_matches) / per_page)

        return extensions, total_pages

    async def __get_all_extensions(self) -> list[tuple[str, str]]:
        if self.__all_packages_cache is None or datetime.now(
            tz=timezone.utc
        ) > self.__last_all_packages_request_time + timedelta(seconds=self.cache_timeout):
            self.log.debug("Requesting PyPI.org RPC API for prebuilt JupyterLab extensions.")
            self.__all_packages_cache = await self.__throttleRequest(
                True,
                self._rpc_client.browse,
                ["Framework :: Jupyter :: JupyterLab :: Extensions :: Prebuilt"],
            )

            # Also include known language packs.  They are not tagged with the
            # prebuilt extension classifier so we fetch their latest versions
            # from the JSON API instead.
            extension_names = {p[0] for p in self.__all_packages_cache}
            packs_to_fetch = [name for name in LANGUAGE_PACKS if name not in extension_names]
            language_pack_results = await asyncio.gather(
                *(self.get_latest_version(name) for name in packs_to_fetch),
                return_exceptions=True,
            )
            for name, result in zip(packs_to_fetch, language_pack_results, strict=True):
                if isinstance(result, Exception):
                    self.log.info(
                        "Failed to fetch latest version for language pack %s: %s",
                        name,
                        result,
                    )
                elif result is not None:
                    self.__all_packages_cache.append((name, result))

            self.__last_all_packages_request_time = datetime.now(tz=timezone.utc)

        return self.__all_packages_cache

    async def install(self, name: str, version: Optional[str] = None) -> ActionResult:  # noqa
        """Install the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            name: The extension name
            version: The version to install; default None (i.e. the latest possible)
        Returns:
            The action result
        """
        if not await self.is_install_allowed(name, version):
            # is_install_allowed will log the reason
            return ActionResult(status="error", message="install is not allowed")

        current_loop = tornado.ioloop.IOLoop.current()
        with (
            tempfile.TemporaryDirectory() as ve_dir,
            tempfile.NamedTemporaryFile(mode="w+", dir=ve_dir, delete=False) as fconstraint,
        ):
            fconstraint.write(f"jupyterlab=={__version__}")
            fconstraint.flush()

            cmdline = [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-input",
                "--quiet",
                "--progress-bar",
                "off",
                "--constraint",
                fconstraint.name,
            ]
            if version is not None:
                cmdline.append(f"{name}=={version}")
            else:
                cmdline.append(name)

            pkg_action = {}
            try:
                tmp_cmd = cmdline.copy()
                tmp_cmd.insert(-1, "--dry-run")
                tmp_cmd.insert(-1, "--report")
                tmp_cmd.insert(-1, "-")
                result = await current_loop.run_in_executor(
                    None, partial(run, tmp_cmd, capture_output=True, check=True)
                )

                action_info = json.loads(result.stdout.decode("utf-8"))
                pkg_action = next(
                    filter(
                        lambda p: p.get("metadata", {}).get("name") == name.replace("_", "-"),
                        action_info.get("install", []),
                    )
                )
            except CalledProcessError as e:
                self.log.debug(f"Fail to get installation report: {e.stderr}", exc_info=e)
            except Exception as err:
                self.log.debug("Fail to get installation report.", exc_info=err)
            else:
                self.log.debug(f"Actions to be executed by pip {json.dumps(action_info)}.")

            self.log.debug(f"Executing '{' '.join(cmdline)}'")

            result = await current_loop.run_in_executor(
                None, partial(run, cmdline, capture_output=True)
            )

            self.log.debug(f"return code: {result.returncode}")
            self.log.debug(f"stdout: {result.stdout.decode('utf-8')}")
            error = result.stderr.decode("utf-8")
            if result.returncode == 0:
                self.log.debug(f"stderr: {error}")
                # Figure out if the package has server or kernel parts
                jlab_metadata = None
                try:
                    download_url: str = pkg_action.get("download_info", {}).get("url")
                    if download_url is not None:
                        response = await self._httpx_client.get(download_url)
                        if response.status_code < 400:  # noqa PLR2004
                            if download_url.endswith(".whl"):
                                with ZipFile(io.BytesIO(response.content)) as wheel:
                                    for filename in filter(
                                        lambda f: Path(f).name == "package.json",
                                        wheel.namelist(),
                                    ):
                                        data = json.loads(wheel.read(filename))
                                        jlab_metadata = data.get("jupyterlab")
                                        if jlab_metadata is not None:
                                            break
                            elif download_url.endswith("tar.gz"):
                                with TarFile(io.BytesIO(response.content)) as sdist:
                                    for filename in filter(
                                        lambda f: Path(f).name == "package.json",
                                        sdist.getnames(),
                                    ):
                                        data = json.load(
                                            sdist.extractfile(sdist.getmember(filename))
                                        )
                                        jlab_metadata = data.get("jupyterlab")
                                        if jlab_metadata is not None:
                                            break
                        else:
                            self.log.debug(f"Failed to get '{download_url}'; {response!s}")
                except Exception as e:
                    self.log.debug("Fail to get package.json.", exc_info=e)

                follow_ups = [
                    "frontend",
                ]
                if jlab_metadata is not None:
                    discovery = jlab_metadata.get("discovery", {})
                    if "kernel" in discovery:
                        follow_ups.append("kernel")
                    if "server" in discovery:
                        follow_ups.append("server")

                return ActionResult(status="ok", needs_restart=follow_ups)
            else:
                self.log.error(f"Failed to install {name}: code {result.returncode}\n{error}")
                return ActionResult(status="error", message=error)

    async def uninstall(self, extension: str) -> ActionResult:
        """Uninstall the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            extension: The extension name
        Returns:
            The action result
        """
        current_loop = tornado.ioloop.IOLoop.current()
        cmdline = [
            sys.executable,
            "-m",
            "pip",
            "uninstall",
            "--yes",
            "--no-input",
            extension,
        ]

        # Figure out if the package has server or kernel parts
        jlab_metadata = None
        try:
            tmp_cmd = cmdline.copy()
            tmp_cmd.remove("--yes")
            result = await current_loop.run_in_executor(
                None, partial(run, tmp_cmd, capture_output=True)
            )
            lines = filter(
                lambda line: line.endswith("package.json"),
                map(lambda line: line.strip(), result.stdout.decode("utf-8").splitlines()),  # noqa
            )
            for filepath in filter(
                lambda f: f.name == "package.json",
                map(Path, lines),
            ):
                data = json.loads(filepath.read_bytes())
                jlab_metadata = data.get("jupyterlab")
                if jlab_metadata is not None:
                    break
        except Exception as e:
            self.log.debug("Fail to list files to be uninstalled.", exc_info=e)

        self.log.debug(f"Executing '{' '.join(cmdline)}'")

        result = await current_loop.run_in_executor(
            None, partial(run, cmdline, capture_output=True)
        )

        self.log.debug(f"return code: {result.returncode}")
        self.log.debug(f"stdout: {result.stdout.decode('utf-8')}")
        error = result.stderr.decode("utf-8")
        if result.returncode == 0:
            self.log.debug(f"stderr: {error}")
            follow_ups = [
                "frontend",
            ]
            if jlab_metadata is not None:
                discovery = jlab_metadata.get("discovery", {})
                if "kernel" in discovery:
                    follow_ups.append("kernel")
                if "server" in discovery:
                    follow_ups.append("server")

            return ActionResult(status="ok", needs_restart=follow_ups)
        else:
            self.log.error(f"Failed to installed {extension}: code {result.returncode}\n{error}")
            return ActionResult(status="error", message=error)

    def _normalize_name(self, name: str) -> str:
        """Normalize extension name.

        Remove `@` from npm scope and replace `/` and `_` by `-`.

        Args:
            name: Extension name
        Returns:
            Normalized name
        """
        return name.replace("@", "").replace("/", "-").replace("_", "-")
