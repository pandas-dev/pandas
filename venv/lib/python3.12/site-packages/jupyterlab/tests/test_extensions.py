# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import json
from unittest.mock import Mock, patch

import pytest
from traitlets.config import Config, Configurable

from jupyterlab.extensions import PyPIExtensionManager, ReadOnlyExtensionManager
from jupyterlab.extensions.manager import ExtensionManager, ExtensionPackage, PluginManager

from . import fake_client_factory


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
