# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import json
import sys
from unittest.mock import AsyncMock, Mock, patch

import pytest
from tornado import web
from traitlets.config import Config, Configurable

from jupyterlab.extensions import PyPIExtensionManager, ReadOnlyExtensionManager
from jupyterlab.extensions.manager import (
    ActionResult,
    ExtensionManager,
    ExtensionPackage,
    PluginManager,
)
from jupyterlab.extensions.pypi import _check_python_version_compatible
from jupyterlab.handlers.extension_manager_handler import ExtensionHandler

from . import fake_client_factory


@pytest.mark.parametrize(
    "requires_python, expected",
    (
        (None, True),  # No requirement
        ("", True),  # Empty requirement
        (">=3.6", True),  # Should pass on any modern Python
        (">=3.6,<4", True),  # Common range
        (">=99.0", False),  # Future version
        ("<3.0", False),  # Very old Python
        ("invalid-specifier", True),  # Invalid specifier should default to compatible
    ),
)
def test_check_python_version_compatible(requires_python, expected):
    """Test the Python version compatibility check function."""
    result, _ = _check_python_version_compatible(requires_python)
    assert result == expected


def test_check_python_version_compatible_current_version():
    """Test that current Python version is compatible with its own specifier."""
    current = f">={sys.version_info.major}.{sys.version_info.minor}"
    assert _check_python_version_compatible(current)[0] is True


@pytest.mark.parametrize(
    "requires_python, expected_allowed",
    (
        (">=3.6", True),  # Compatible with current Python
        (">=99.0", False),  # Incompatible - future version
        (None, True),  # No requirement - should be allowed
    ),
)
@patch("jupyterlab.extensions.pypi.LANGUAGE_PACKS", ())
@patch("jupyterlab.extensions.pypi.xmlrpc.client")
async def test_extension_package_allowed_set_from_python_compatibility(
    mocked_rpcclient, requires_python, expected_allowed
):
    """Test that allowed attribute is correctly set based on Python version compatibility."""
    proxy = Mock(browse=Mock(return_value=[["test-extension", "1.0.0"]]))
    mocked_rpcclient.ServerProxy = Mock(return_value=proxy)

    manager = PyPIExtensionManager()

    async def mock_pkg_metadata(name, version, base_url):
        return {
            "name": "test-extension",
            "summary": "A test extension",
            "requires_python": requires_python,
        }

    manager._fetch_package_metadata = mock_pkg_metadata

    extensions, _ = await manager.list_extensions("test")

    assert len(extensions) == 1
    ext = extensions[0]
    assert ext.name == "test-extension"
    assert ext.allowed == expected_allowed
    if not expected_allowed:
        assert "Requires Python" in ext.description
    else:
        assert ext.description == "A test extension"


@pytest.mark.parametrize(
    "version, expected",
    (
        ("1", "1"),
        ("1.0", "1.0"),
        ("1.0.0", "1.0.0"),
        ("1.0.0a52", "1.0.0-alpha.52"),
        ("1.0.0b3", "1.0.0-beta.3"),
        ("1.0.0rc22", "1.0.0-rc.22"),
        ("1.0.0rc23.post2", "1.0.0-rc.23"),
        ("1.0.0rc24.dev2", "1.0.0-rc.24"),
        ("1.0.0rc25.post4.dev2", "1.0.0-rc.25"),
    ),
)
def test_ExtensionManager_get_semver_version(version, expected):
    assert ExtensionManager.get_semver_version(version) == expected


async def test_ExtensionManager_list_extensions_installed(monkeypatch):
    extension1 = ExtensionPackage("extension1", "Extension 1 description", "", "prebuilt")

    async def mock_installed(*args, **kwargs):
        return {"extension1": extension1}

    monkeypatch.setattr(ReadOnlyExtensionManager, "_get_installed_extensions", mock_installed)

    manager = ReadOnlyExtensionManager()

    extensions = await manager.list_extensions()

    assert extensions == ([extension1], 1)


async def test_ExtensionManager_list_extensions_query(monkeypatch):
    extension1 = ExtensionPackage("extension1", "Extension 1 description", "", "prebuilt")
    extension2 = ExtensionPackage("extension2", "Extension 2 description", "", "prebuilt")

    async def mock_list(*args, **kwargs):
        return {"extension1": extension1, "extension2": extension2}, None

    monkeypatch.setattr(ReadOnlyExtensionManager, "list_packages", mock_list)

    manager = ReadOnlyExtensionManager()

    extensions = await manager.list_extensions("ext")

    assert extensions == ([extension1, extension2], 1)


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_list_extensions_query_allow(mock_client, monkeypatch):
    extension1 = ExtensionPackage("extension1", "Extension 1 description", "", "prebuilt")
    extension2 = ExtensionPackage("extension2", "Extension 2 description", "", "prebuilt")

    mock_client.body = json.dumps({"allowed_extensions": [{"name": "extension1"}]}).encode()

    async def mock_list(*args, **kwargs):
        return {"extension1": extension1, "extension2": extension2}, None

    monkeypatch.setattr(ReadOnlyExtensionManager, "list_packages", mock_list)

    manager = ReadOnlyExtensionManager(
        ext_options={"allowed_extensions_uris": {"http://dummy-allowed-extension"}},
    )

    extensions = await manager.list_extensions("ext")

    assert extensions == ([extension1], 1)


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_list_extensions_query_block(mock_client, monkeypatch):
    extension1 = ExtensionPackage("extension1", "Extension 1 description", "", "prebuilt")
    extension2 = ExtensionPackage("extension2", "Extension 2 description", "", "prebuilt")

    mock_client.body = json.dumps({"blocked_extensions": [{"name": "extension1"}]}).encode()

    async def mock_list(*args, **kwargs):
        return {"extension1": extension1, "extension2": extension2}, None

    monkeypatch.setattr(ReadOnlyExtensionManager, "list_packages", mock_list)

    manager = ReadOnlyExtensionManager(
        ext_options={"blocked_extensions_uris": {"http://dummy-blocked-extension"}}
    )

    extensions = await manager.list_extensions("ext")

    assert extensions == ([extension2], 1)


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_list_extensions_query_allow_block(mock_client, monkeypatch):
    extension1 = ExtensionPackage("extension1", "Extension 1 description", "", "prebuilt")
    extension2 = ExtensionPackage("extension2", "Extension 2 description", "", "prebuilt")

    mock_client.body = json.dumps(
        {
            "allowed_extensions": [{"name": "extension1"}],
            "blocked_extensions": [{"name": "extension1"}],
        }
    ).encode()

    async def mock_list(*args, **kwargs):
        return {"extension1": extension1, "extension2": extension2}, None

    monkeypatch.setattr(ReadOnlyExtensionManager, "list_packages", mock_list)

    manager = ReadOnlyExtensionManager(
        ext_options={
            "allowed_extensions_uris": {"http://dummy-allowed-extension"},
            "blocked_extensions_uris": {"http://dummy-blocked-extension"},
        }
    )

    extensions = await manager.list_extensions("ext")

    assert extensions == ([extension1], 1)


async def test_ExtensionManager_is_install_allowed_no_policy():
    manager = PyPIExtensionManager()
    assert await manager.is_install_allowed("any-extension") is True


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_is_install_allowed_allowlist(mock_client):
    mock_client.body = json.dumps({"allowed_extensions": [{"name": "jupyterlab-git"}]}).encode()
    manager = PyPIExtensionManager(
        ext_options={"allowed_extensions_uris": {"http://dummy-allowed-extension"}}
    )
    assert await manager.is_install_allowed("jupyterlab-git") is True
    assert await manager.is_install_allowed("jupyterlab-evil") is False


@pytest.mark.parametrize(
    "name",
    ["jupyterlab-git", "JupyterLab-Git", "JupyterLab.Git", "jupyterlab_git"],
)
async def test_pypi_manager_allows_canonical_package_name_allowlist_matches(name):
    manager = PyPIExtensionManager(
        ext_options={"allowed_extensions_uris": {"http://dummy-allowed-extension"}}
    )
    manager._listings_cache = {"jupyterlab-git": {"name": "jupyterlab-git"}}
    if manager._listing_fetch is not None:
        manager._listing_fetch.stop()

    assert await manager.is_install_allowed(name) is True
    assert await manager.is_install_allowed("jupyterlab-evil") is False


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_is_install_allowed_blocklist(mock_client):
    mock_client.body = json.dumps({"blocked_extensions": [{"name": "jupyterlab-evil"}]}).encode()
    manager = PyPIExtensionManager(
        ext_options={"blocked_extensions_uris": {"http://dummy-blocked-extension"}}
    )
    assert await manager.is_install_allowed("jupyterlab-evil") is False
    assert await manager.is_install_allowed("jupyterlab-git") is True


@pytest.mark.parametrize(
    "name",
    ["jupyterlab-git", "JupyterLab-Git", "JupyterLab.Git", "jupyterlab_git"],
)
async def test_handler_blocks_pypi_canonical_package_name_blocklist_matches(name):
    manager = PyPIExtensionManager(
        ext_options={"blocked_extensions_uris": {"http://dummy-blocked-extension"}}
    )
    manager._listings_cache = {"jupyterlab-git": {"name": "jupyterlab-git"}}
    if manager._listing_fetch is not None:
        manager._listing_fetch.stop()
    manager.install = AsyncMock(return_value=ActionResult(status="ok", needs_restart=[]))

    handler = Mock()
    handler.current_user = "user"
    handler.get_json_body.return_value = {"cmd": "install", "extension_name": name}
    handler.manager = manager
    handler.set_status = Mock()
    handler.finish = Mock()

    with pytest.raises(web.HTTPError) as exc_info:
        await ExtensionHandler.post(handler)
    assert exc_info.value.status_code == 422
    assert "was blocked" in exc_info.value.log_message
    manager.install.assert_not_called()


async def test_pypi_manager_install_blocks_policy_denial_before_pip():
    manager = PyPIExtensionManager(
        ext_options={"blocked_extensions_uris": {"http://dummy-blocked-extension"}}
    )
    manager._listings_cache = {"jupyterlab-git": {"name": "jupyterlab-git"}}
    if manager._listing_fetch is not None:
        manager._listing_fetch.stop()

    with patch("jupyterlab.extensions.pypi.run") as run_mock:
        result = await manager.install("JupyterLab.Git")

    assert result == ActionResult(status="error", message="install is not allowed")
    run_mock.assert_not_called()


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_ExtensionManager_is_install_allowed_allowlist_takes_precedence(mock_client):
    mock_client.body = json.dumps(
        {
            "allowed_extensions": [{"name": "jupyterlab-git"}],
            "blocked_extensions": [{"name": "jupyterlab-git"}],  # allowlist takes precedence
        }
    ).encode()
    manager = PyPIExtensionManager(
        ext_options={
            "allowed_extensions_uris": {"http://dummy-allowed-extension"},
            "blocked_extensions_uris": {"http://dummy-blocked-extension"},
        }
    )
    assert await manager.is_install_allowed("jupyterlab-git") is True
    assert await manager.is_install_allowed("jupyterlab-evil") is False


async def test_ExtensionManager_install():
    manager = ReadOnlyExtensionManager()

    result = await manager.install("extension1")

    assert result.status == "error"
    assert result.message == "Extension installation not supported."


async def test_ExtensionManager_uninstall():
    manager = ReadOnlyExtensionManager()

    result = await manager.uninstall("extension1")

    assert result.status == "error"
    assert result.message == "Extension removal not supported."


@patch("jupyterlab.extensions.pypi.LANGUAGE_PACKS", ())
@patch("jupyterlab.extensions.pypi.xmlrpc.client")
async def test_ExtensionManager_list_extensions_query_sort(mocked_rpcclient):
    extension_data = [
        {
            "name": "jupyterlab-apod",
            "project_urls": {
                "Homepage": "https://github.com/jupyterlab/jupyterlab_apod",
            },
        },
        {
            "name": "jupyterlab-gitlab",
            "project_urls": {
                "Homepage": "https://github.com/jupyterlab-contrib/jupyterlab-gitlab/issues",
            },
        },
        {
            "name": "jupyterlab-git",
            "project_url": "https://github.com/jupyterlab/jupyterlab-git",
        },
        {
            "name": "jupyterlab-rainbow-brackets",
            "project_url": "https://github.com/krassowski/jupyterlab-rainbow-brackets",
        },
        {"name": "nbdime", "home_page": "https://github.com/jupyter/nbdime"},
        {
            "name": "rise",
            "project_urls": {
                "Source Code": "https://github.com/jupyterlab-contrib/rise",
            },
        },
    ]

    proxy = Mock(
        browse=Mock(return_value=[[extension["name"], "1.0.0"] for extension in extension_data]),
    )
    mocked_rpcclient.ServerProxy = Mock(return_value=proxy)

    manager = PyPIExtensionManager()

    extensions = {
        extension["name"]: {"version": "1.0.0", **extension} for extension in extension_data
    }

    async def mock_pkg_metadata(name, l, b):  # noqa
        return extensions[name]

    manager._fetch_package_metadata = mock_pkg_metadata

    first_page, pages_count = await manager.list_extensions("", per_page=3)
    assert [extension.name for extension in first_page] == [
        # jupyter/jupyterlab
        "jupyterlab-git",
        "nbdime",
        # jupyterlab-contrib
        "jupyterlab-gitlab",
    ]
    assert pages_count == 2
    second_page, pages_count = await manager.list_extensions("", page=2, per_page=3)
    assert [extension.name for extension in second_page] == [
        # jupyterlab-contrib
        "rise",
        # other third-party
        "jupyterlab-rainbow-brackets",
        # example extensions
        "jupyterlab-apod",
    ]


@patch("jupyterlab.extensions.pypi.LANGUAGE_PACKS", ())
@patch("jupyterlab.extensions.pypi.xmlrpc.client")
async def test_PyPiExtensionManager_list_extensions_query(mocked_rpcclient):
    extension1 = ExtensionPackage(
        name="jupyterlab-git",
        description="A JupyterLab extension for version control using git",
        homepage_url="https://github.com/jupyterlab/jupyterlab-git",
        pkg_type="prebuilt",
        latest_version="0.37.1",
        author="Jupyter Development Team",
        license="BSD-3-Clause",
        package_manager_url="https://pypi.org/project/jupyterlab-git/",
    )
    extension2 = ExtensionPackage(
        name="jupyterlab-github",
        description="JupyterLab viewer for GitHub repositories",
        homepage_url="https://github.com/jupyterlab/jupyterlab-github/blob/main/README.md",
        pkg_type="prebuilt",
        latest_version="3.0.1",
        author="Ian Rose",
        license="BSD-3-Clause",
        bug_tracker_url="https://github.com/jupyterlab/jupyterlab-github/issues",
        package_manager_url="https://pypi.org/project/jupyterlab-github/",
        repository_url="https://github.com/jupyterlab/jupyterlab-github",
    )

    proxy = Mock(
        browse=Mock(
            return_value=[
                ["jupyterlab-git", "0.33.0"],
                ["jupyterlab-git", "0.34.0"],
                ["jupyterlab-git", "0.34.1"],
                ["jupyterlab-git", "0.37.0"],
                ["jupyterlab-git", "0.37.1"],
                ["jupyterlab-github", "3.0.0"],
                ["jupyterlab-github", "3.0.1"],
            ]
        ),
    )
    mocked_rpcclient.ServerProxy = Mock(return_value=proxy)

    manager = PyPIExtensionManager()

    async def mock_pkg_metadata(n, l, b):  # noqa
        return (
            {
                "name": "jupyterlab-git",
                "version": "0.37.1",
                "stable_version": None,
                "bugtrack_url": None,
                "package_url": "https://pypi.org/project/jupyterlab-git/",
                "release_url": "https://pypi.org/project/jupyterlab-git/0.37.1/",
                "docs_url": None,
                "home_page": "https://github.com/jupyterlab/jupyterlab-git",
                "download_url": "",
                "project_url": "",
                "project_urls": {},
                "author": "Jupyter Development Team",
                "author_email": "",
                "maintainer": "",
                "maintainer_email": "",
                "summary": "A JupyterLab extension for version control using git",
                "license": "BSD-3-Clause",
                "keywords": "Jupyter,JupyterLab,JupyterLab3,jupyterlab-extension,Git",
                "platform": "Linux",
                "classifiers": [
                    "Framework :: Jupyter",
                    "Framework :: Jupyter :: JupyterLab",
                    "Framework :: Jupyter :: JupyterLab :: 3",
                    "Framework :: Jupyter :: JupyterLab :: Extensions",
                    "Framework :: Jupyter :: JupyterLab :: Extensions :: Prebuilt",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: BSD License",
                    "Programming Language :: Python",
                    "Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.10",
                    "Programming Language :: Python :: 3.6",
                    "Programming Language :: Python :: 3.7",
                    "Programming Language :: Python :: 3.8",
                    "Programming Language :: Python :: 3.9",
                ],
                "requires": [],
                "requires_dist": [
                    "jupyter-server",
                    "nbdime (~=3.1)",
                    "nbformat",
                    "packaging",
                    "pexpect",
                    "coverage ; extra == 'dev'",
                    "jupyterlab (~=3.0) ; extra == 'dev'",
                    "pre-commit ; extra == 'dev'",
                    "pytest ; extra == 'dev'",
                    "pytest-asyncio ; extra == 'dev'",
                    "pytest-cov ; extra == 'dev'",
                    "pytest-tornasync ; extra == 'dev'",
                    "coverage ; extra == 'tests'",
                    "jupyterlab (~=3.0) ; extra == 'tests'",
                    "pre-commit ; extra == 'tests'",
                    "pytest ; extra == 'tests'",
                    "pytest-asyncio ; extra == 'tests'",
                    "pytest-cov ; extra == 'tests'",
                    "pytest-tornasync ; extra == 'tests'",
                    "hybridcontents ; extra == 'tests'",
                    "jupytext ; extra == 'tests'",
                ],
                "provides": [],
                "provides_dist": [],
                "obsoletes": [],
                "obsoletes_dist": [],
                "requires_python": "<4,>=3.6",
                "requires_external": [],
                "_pypi_ordering": 55,
                "downloads": {"last_day": -1, "last_week": -1, "last_month": -1},
                "cheesecake_code_kwalitee_id": None,
                "cheesecake_documentation_id": None,
                "cheesecake_installability_id": None,
            }
            if n == "jupyterlab-git"
            else {
                "name": "jupyterlab-github",
                "version": "3.0.1",
                "stable_version": None,
                "bugtrack_url": None,
                "package_url": "https://pypi.org/project/jupyterlab-github/",
                "release_url": "https://pypi.org/project/jupyterlab-github/3.0.1/",
                "docs_url": None,
                "home_page": "",
                "download_url": "",
                "project_url": "",
                "project_urls": {
                    "Homepage": "https://github.com/jupyterlab/jupyterlab-github/blob/main/README.md",
                    "Bug Tracker": "https://github.com/jupyterlab/jupyterlab-github/issues",
                    "Source Code": "https://github.com/jupyterlab/jupyterlab-github",
                },
                "author": "Ian Rose",
                "author_email": "jupyter@googlegroups.com",
                "maintainer": "",
                "maintainer_email": "",
                "summary": "JupyterLab viewer for GitHub repositories",
                "license": "BSD-3-Clause",
                "keywords": "Jupyter,JupyterLab,JupyterLab3",
                "platform": "Linux",
                "classifiers": [
                    "Framework :: Jupyter",
                    "Framework :: Jupyter :: JupyterLab",
                    "Framework :: Jupyter :: JupyterLab :: 3",
                    "Framework :: Jupyter :: JupyterLab :: Extensions",
                    "Framework :: Jupyter :: JupyterLab :: Extensions :: Prebuilt",
                    "License :: OSI Approved :: BSD License",
                    "Programming Language :: Python",
                    "Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.6",
                    "Programming Language :: Python :: 3.7",
                    "Programming Language :: Python :: 3.8",
                    "Programming Language :: Python :: 3.9",
                ],
                "requires": [],
                "requires_dist": ["jupyterlab (~=3.0)"],
                "provides": [],
                "provides_dist": [],
                "obsoletes": [],
                "obsoletes_dist": [],
                "requires_python": ">=3.6",
                "requires_external": [],
                "_pypi_ordering": 12,
                "downloads": {"last_day": -1, "last_week": -1, "last_month": -1},
                "cheesecake_code_kwalitee_id": None,
                "cheesecake_documentation_id": None,
                "cheesecake_installability_id": None,
            }
        )

    manager._fetch_package_metadata = mock_pkg_metadata

    extensions = await manager.list_extensions("git")

    assert extensions == ([extension1, extension2], 1)


async def test_PyPiExtensionManager_custom_server_url():
    BASE_URL = "https://mylocal.pypi.server/pypi"  # noqa

    parent = Configurable(config=Config({"PyPIExtensionManager": {"base_url": BASE_URL}}))

    manager = PyPIExtensionManager(parent=parent)

    assert manager.base_url == BASE_URL


LEVELS = ["user", "sys_prefix", "system"]


@pytest.mark.parametrize("level", LEVELS)
async def test_PyPiExtensionManager_custom_level(level):
    parent = Configurable(config=Config({"PyPIExtensionManager": {"level": level}}))
    manager = PyPIExtensionManager(parent=parent)
    assert manager.level == level


@pytest.mark.parametrize("level", LEVELS)
async def test_PyPiExtensionManager_inherits_custom_level(level):
    parent = Configurable(config=Config({"PluginManager": {"level": level}}))
    manager = PyPIExtensionManager(parent=parent)
    assert manager.level == level


@pytest.mark.parametrize("level", LEVELS)
async def test_PluginManager_custom_level(level):
    parent = Configurable(config=Config({"PluginManager": {"level": level}}))
    manager = PluginManager(parent=parent)
    assert manager.level == level


async def test_PluginManager_default_level():
    manager = PluginManager()
    assert manager.level == "sys_prefix"


@pytest.mark.parametrize(
    "name,expected",
    [
        ("jupyterlab-git", True),
        ("my_extension.pkg", True),
        # GHSA-37w4-hwhx-4rc4: VCS/URL names must be blocked
        ("git+https://github.com/attacker/malicious.git", False),
        ("/tmp/local-pkg", False),  # noqa: S108
        ("http://evil.com/pkg.tar.gz", False),
    ],
)
async def test_pypi_manager_is_install_allowed_rejects_non_pypi_names(name, expected):
    manager = PyPIExtensionManager()
    assert await manager.is_install_allowed(name) is expected


@pytest.mark.parametrize(
    "version,expected",
    [(None, True), ("1.2.3", True), ("not valid version", False)],
)
async def test_pypi_manager_is_install_allowed_rejects_non_pypi_versions(version, expected):
    manager = PyPIExtensionManager()
    assert await manager.is_install_allowed("package", version) is expected


async def test_pypi_manager_install_blocks_when_policy_denies():
    manager = PyPIExtensionManager()
    manager.is_install_allowed = AsyncMock(return_value=False)

    with patch("jupyterlab.extensions.pypi.tornado.ioloop.IOLoop.current") as current_loop:
        result = await manager.install("jupyterlab-evil", "1.0.0")

    assert result.status == "error"
    assert result.message == "install is not allowed"
    manager.is_install_allowed.assert_awaited_once_with("jupyterlab-evil", "1.0.0")
    current_loop.assert_not_called()


async def test_handler_blocks_install_when_policy_denies():
    handler = Mock()
    handler.current_user = "user"
    handler.get_json_body.return_value = {"cmd": "install", "extension_name": "jupyterlab-evil"}
    handler.manager.is_install_allowed = AsyncMock(return_value=False)

    with pytest.raises(web.HTTPError) as exc_info:
        await ExtensionHandler.post(handler)
    assert exc_info.value.status_code == 422
    assert "was blocked" in exc_info.value.log_message
    handler.manager.is_install_allowed.assert_called_once_with("jupyterlab-evil", None)
    handler.manager.install.assert_not_called()
