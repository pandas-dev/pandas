"""Test installation of JupyterLab extensions"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import glob
import json
import logging
import os
import platform
import shutil
import subprocess
import sys
from os.path import join as pjoin
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import pytest
from jupyter_core import paths

from jupyterlab import commands
from jupyterlab.commands import (
    DEV_DIR,
    AppOptions,
    _compare_ranges,
    _test_overlap,
    build,
    build_check,
    check_extension,
    disable_extension,
    enable_extension,
    get_app_info,
    get_app_version,
    install_extension,
    link_package,
    list_extensions,
    uninstall_extension,
    unlink_package,
    update_extension,
)
from jupyterlab.coreconfig import CoreConfig, _get_default_core_data

here = os.path.dirname(os.path.abspath(__file__))


def touch(file, mtime=None):
    """ensure a file exists, and set its modification time

    returns the modification time of the file
    """
    dirname = os.path.dirname(file)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    open(file, "a").close()
    # set explicit mtime
    if mtime:
        atime = os.stat(file).st_atime
        os.utime(file, (atime, mtime))
    return os.stat(file).st_mtime


# @pytest.fixture()
# def resource():
#     print("setup")
#     yield "resource"
#    print("teardown")


class AppHandlerTest(TestCase):
    def tempdir(self):
        td = TemporaryDirectory()
        self.tempdirs.append(td)
        return td.name

    def setUp(self):
        # Any TemporaryDirectory objects appended to this list will be cleaned
        # up at the end of the test run.
        self.tempdirs = []
        self.devnull = open(os.devnull, "w")  # noqa

        @self.addCleanup
        def cleanup_tempdirs():
            for d in self.tempdirs:
                d.cleanup()

        self.test_dir = self.tempdir()

        self.data_dir = pjoin(self.test_dir, "data")
        self.config_dir = pjoin(self.test_dir, "config")
        self.pkg_names = {}

        # Copy in the mock packages.
        for name in ["extension", "incompat", "package", "mimeextension"]:
            src = pjoin(here, "mock_packages", name)

            def ignore(dname, files):
                if "node_modules" in dname:
                    files = []
                if "node_modules" in files:
                    files.remove("node_modules")
                return dname, files

            dest = pjoin(self.test_dir, name)
            shutil.copytree(src, dest, ignore=ignore)

            # Make a node modules folder so npm install is not called.
            if not os.path.exists(pjoin(dest, "node_modules")):
                os.makedirs(pjoin(dest, "node_modules"))

            setattr(self, "mock_" + name, dest)
            with open(pjoin(dest, "package.json")) as fid:
                data = json.load(fid)
            self.pkg_names[name] = data["name"]

        self.patches = []
        p = patch.dict(
            "os.environ",
            {
                "JUPYTER_CONFIG_DIR": self.config_dir,
                "JUPYTER_DATA_DIR": self.data_dir,
                "JUPYTERLAB_DIR": pjoin(self.data_dir, "lab"),
            },
        )
        self.patches.append(p)
        for mod in [paths]:
            if hasattr(mod, "ENV_JUPYTER_PATH"):
                p = patch.object(mod, "ENV_JUPYTER_PATH", [self.data_dir])
                self.patches.append(p)
            if hasattr(mod, "ENV_CONFIG_PATH"):
                p = patch.object(mod, "ENV_CONFIG_PATH", [self.config_dir])
                self.patches.append(p)
            if hasattr(mod, "CONFIG_PATH"):
                p = patch.object(mod, "CONFIG_PATH", self.config_dir)
                self.patches.append(p)
            if hasattr(mod, "BUILD_PATH"):
                p = patch.object(mod, "BUILD_PATH", self.data_dir)
                self.patches.append(p)
        for p in self.patches:
            p.start()
            self.addCleanup(p.stop)

        # verify our patches
        self.assertEqual(paths.ENV_CONFIG_PATH, [self.config_dir])
        self.assertEqual(paths.ENV_JUPYTER_PATH, [self.data_dir])
        self.assertEqual(
            Path(commands.get_app_dir()).resolve(), (Path(self.data_dir) / "lab").resolve()
        )

        self.app_dir = commands.get_app_dir()

        # Set pinned extension names
        self.pinned_packages = ["jupyterlab-test-extension@1.0", "jupyterlab-test-extension@2.0"]


class TestExtension(AppHandlerTest):
    def test_install_extension(self):
        assert install_extension(self.mock_extension) is True
        path = pjoin(self.app_dir, "extensions", "*.tgz")
        assert glob.glob(path)
        extensions = get_app_info()["extensions"]
        name = self.pkg_names["extension"]
        assert name in extensions
        assert check_extension(name)

    def test_install_twice(self):
        assert install_extension(self.mock_extension) is True
        path = pjoin(self.app_dir, "extensions", "*.tgz")
        assert install_extension(self.mock_extension) is True
        assert glob.glob(path)
        extensions = get_app_info()["extensions"]
        name = self.pkg_names["extension"]
        assert name in extensions
        assert check_extension(name)

    def test_install_mime_renderer(self):
        install_extension(self.mock_mimeextension)
        name = self.pkg_names["mimeextension"]
        assert name in get_app_info()["extensions"]
        assert check_extension(name)

        assert uninstall_extension(name) is True
        assert name not in get_app_info()["extensions"]
        assert not check_extension(name)

    def test_install_incompatible(self):
        with pytest.raises(ValueError) as excinfo:
            install_extension(self.mock_incompat)
        assert "Conflicting Dependencies" in str(excinfo.value)
        assert not check_extension(self.pkg_names["incompat"])

    def test_install_failed(self):
        path = self.mock_package
        with pytest.raises(ValueError):
            install_extension(path)
        with open(pjoin(path, "package.json")) as fid:
            data = json.load(fid)
        extensions = get_app_info()["extensions"]
        name = data["name"]
        assert name not in extensions
        assert not check_extension(name)

    def test_validation(self):
        path = self.mock_extension
        os.remove(pjoin(path, "index.js"))
        with pytest.raises(ValueError):
            install_extension(path)
        assert not check_extension(self.pkg_names["extension"])

        path = self.mock_mimeextension
        os.remove(pjoin(path, "index.js"))
        with pytest.raises(ValueError):
            install_extension(path)
        assert not check_extension(self.pkg_names["mimeextension"])

    def test_uninstall_extension(self):
        assert install_extension(self.mock_extension) is True
        name = self.pkg_names["extension"]
        assert check_extension(name)
        assert uninstall_extension(self.pkg_names["extension"]) is True
        path = pjoin(self.app_dir, "extensions", "*.tgz")
        assert not glob.glob(path)
        extensions = get_app_info()["extensions"]
        assert name not in extensions
        assert not check_extension(name)

    def test_uninstall_all_extensions(self):
        install_extension(self.mock_extension)
        install_extension(self.mock_mimeextension)
        ext_name = self.pkg_names["extension"]
        mime_ext_name = self.pkg_names["mimeextension"]
        assert check_extension(ext_name) is True
        assert check_extension(mime_ext_name) is True
        assert uninstall_extension(all_=True) is True
        extensions = get_app_info()["extensions"]
        assert ext_name not in extensions
        assert mime_ext_name not in extensions

    @pytest.mark.slow
    def test_uninstall_core_extension(self):
        assert uninstall_extension("@jupyterlab/console-extension") is True
        app_dir = self.app_dir
        build()
        with open(pjoin(app_dir, "staging", "package.json")) as fid:
            data = json.load(fid)
        extensions = data["jupyterlab"]["extensions"]
        assert "@jupyterlab/console-extension" not in extensions
        assert not check_extension("@jupyterlab/console-extension")

        assert install_extension("@jupyterlab/console-extension") is True
        build()
        with open(pjoin(app_dir, "staging", "package.json")) as fid:
            data = json.load(fid)
        extensions = data["jupyterlab"]["extensions"]
        assert "@jupyterlab/console-extension" in extensions
        assert check_extension("@jupyterlab/console-extension")

    def test_install_and_uninstall_pinned(self):
        """
        You should be able to install different versions of the same extension with different
        pinned names and uninstall them with those names.
        """
        NAMES = ["test-1", "test-2"]  # noqa
        assert install_extension(self.pinned_packages[0], pin=NAMES[0])
        assert install_extension(self.pinned_packages[1], pin=NAMES[1])

        extensions = get_app_info()["extensions"]
        assert NAMES[0] in extensions
        assert NAMES[1] in extensions
        assert check_extension(NAMES[0])
        assert check_extension(NAMES[1])

        # Uninstall
        assert uninstall_extension(NAMES[0])
        assert uninstall_extension(NAMES[1])

        extensions = get_app_info()["extensions"]
        assert NAMES[0] not in extensions
        assert NAMES[1] not in extensions
        assert not check_extension(NAMES[0])
        assert not check_extension(NAMES[1])

    @pytest.mark.skipif(
        platform.system() == "Windows", reason="running npm pack fails on windows CI"
    )
    def test_install_and_uninstall_pinned_folder(self):
        """
        Same as above test, but installs from a local folder instead of from npm.
        """
        # Download each version of the package from NPM:
        base_dir = Path(self.tempdir())

        # The archive file names are printed to stdout when run `npm pack`
        packages = [
            subprocess.run(
                ["npm", "pack", name],  # noqa S603 S607
                stdout=subprocess.PIPE,
                text=True,
                check=True,
                cwd=str(base_dir),
            ).stdout.strip()
            for name in self.pinned_packages
        ]

        shutil.unpack_archive(str(base_dir / packages[0]), str(base_dir / "1"))
        shutil.unpack_archive(str(base_dir / packages[1]), str(base_dir / "2"))
        # Change pinned packages to be these directories now, so we install from these folders
        self.pinned_packages = [str(base_dir / "1" / "package"), str(base_dir / "2" / "package")]
        self.test_install_and_uninstall_pinned()

    def test_link_extension(self):
        path = self.mock_extension
        name = self.pkg_names["extension"]
        link_package(path)
        linked = get_app_info()["linked_packages"]
        assert name not in linked
        assert name in get_app_info()["extensions"]
        assert check_extension(name)
        assert unlink_package(path) is True
        linked = get_app_info()["linked_packages"]
        assert name not in linked
        assert name not in get_app_info()["extensions"]
        assert not check_extension(name)

    def test_link_package(self):
        path = self.mock_package
        name = self.pkg_names["package"]
        assert link_package(path) is True
        linked = get_app_info()["linked_packages"]
        assert name in linked
        assert name not in get_app_info()["extensions"]
        assert check_extension(name)
        assert unlink_package(path)
        linked = get_app_info()["linked_packages"]
        assert name not in linked
        assert not check_extension(name)

    def test_unlink_package(self):
        target = self.mock_package
        assert link_package(target) is True
        assert unlink_package(target) is True
        linked = get_app_info()["linked_packages"]
        name = self.pkg_names["package"]
        assert name not in linked
        assert not check_extension(name)

    def test_list_extensions(self):
        assert install_extension(self.mock_extension) is True
        list_extensions()

    def test_app_dir(self):
        app_dir = self.tempdir()
        options = AppOptions(app_dir=app_dir)

        assert install_extension(self.mock_extension, app_options=options) is True
        path = pjoin(app_dir, "extensions", "*.tgz")
        assert glob.glob(path)
        extensions = get_app_info(app_options=options)["extensions"]
        ext_name = self.pkg_names["extension"]
        assert ext_name in extensions
        assert check_extension(ext_name, app_options=options)

        assert uninstall_extension(self.pkg_names["extension"], app_options=options) is True
        path = pjoin(app_dir, "extensions", "*.tgz")
        assert not glob.glob(path)
        extensions = get_app_info(app_options=options)["extensions"]
        assert ext_name not in extensions
        assert not check_extension(ext_name, app_options=options)

        assert link_package(self.mock_package, app_options=options) is True
        linked = get_app_info(app_options=options)["linked_packages"]
        pkg_name = self.pkg_names["package"]
        assert pkg_name in linked
        assert check_extension(pkg_name, app_options=options)

        assert unlink_package(self.mock_package, app_options=options) is True
        linked = get_app_info(app_options=options)["linked_packages"]
        assert pkg_name not in linked
        assert not check_extension(pkg_name, app_options=options)

    def test_app_dir_use_sys_prefix(self):
        app_dir = self.tempdir()
        options = AppOptions(app_dir=app_dir)
        if os.path.exists(self.app_dir):
            os.removedirs(self.app_dir)

        assert install_extension(self.mock_extension) is True
        path = pjoin(app_dir, "extensions", "*.tgz")
        assert not glob.glob(path)
        extensions = get_app_info(app_options=options)["extensions"]
        ext_name = self.pkg_names["extension"]
        assert ext_name in extensions
        assert check_extension(ext_name, app_options=options)

    def test_app_dir_disable_sys_prefix(self):
        app_dir = self.tempdir()
        options = AppOptions(app_dir=app_dir, use_sys_dir=False)
        if os.path.exists(self.app_dir):
            os.removedirs(self.app_dir)

        assert install_extension(self.mock_extension) is True
        path = pjoin(app_dir, "extensions", "*.tgz")
        assert not glob.glob(path)
        extensions = get_app_info(app_options=options)["extensions"]
        ext_name = self.pkg_names["extension"]
        assert ext_name not in extensions
        assert not check_extension(ext_name, app_options=options)

    def test_app_dir_shadowing(self):
        app_dir = self.tempdir()
        sys_dir = self.app_dir
        app_options = AppOptions(app_dir=app_dir)
        if os.path.exists(sys_dir):
            os.removedirs(sys_dir)

        assert install_extension(self.mock_extension) is True
        sys_path = pjoin(sys_dir, "extensions", "*.tgz")
        assert glob.glob(sys_path)
        app_path = pjoin(app_dir, "extensions", "*.tgz")
        assert not glob.glob(app_path)
        extensions = get_app_info(app_options=app_options)["extensions"]
        ext_name = self.pkg_names["extension"]
        assert ext_name in extensions
        assert check_extension(ext_name, app_options=app_options)

        assert install_extension(self.mock_extension, app_options=app_options) is True
        assert glob.glob(app_path)
        extensions = get_app_info(app_options=app_options)["extensions"]
        assert ext_name in extensions
        assert check_extension(ext_name, app_options=app_options)

        assert uninstall_extension(self.pkg_names["extension"], app_options=app_options) is True
        assert not glob.glob(app_path)
        assert glob.glob(sys_path)
        extensions = get_app_info(app_options=app_options)["extensions"]
        assert ext_name in extensions
        assert check_extension(ext_name, app_options=app_options)

        assert uninstall_extension(self.pkg_names["extension"], app_options=app_options) is True
        assert not glob.glob(app_path)
        assert not glob.glob(sys_path)
        extensions = get_app_info(app_options=app_options)["extensions"]
        assert ext_name not in extensions
        assert not check_extension(ext_name, app_options=app_options)

    @pytest.mark.slow
    def test_build(self):
        assert install_extension(self.mock_extension) is True
        build()
        # check staging directory.
        entry = pjoin(self.app_dir, "staging", "build", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

        # check static directory.
        entry = pjoin(self.app_dir, "static", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

    @pytest.mark.slow
    @pytest.mark.skipif(not os.path.exists(DEV_DIR), reason="Not in git checkout")
    def test_build_splice_packages(self):
        app_options = AppOptions(splice_source=True)
        assert install_extension(self.mock_extension) is True
        build(app_options=app_options)
        assert "-spliced" in get_app_version(app_options)
        # check staging directory.
        entry = pjoin(self.app_dir, "staging", "build", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

        # check static directory.
        entry = pjoin(self.app_dir, "static", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

    @pytest.mark.slow
    def test_build_custom(self):
        assert install_extension(self.mock_extension) is True
        build(name="foo", version="1.0", static_url="bar")

        # check static directory.
        entry = pjoin(self.app_dir, "static", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

        pkg = pjoin(self.app_dir, "static", "package.json")
        with open(pkg) as fid:
            data = json.load(fid)
        assert data["jupyterlab"]["name"] == "foo"
        assert data["jupyterlab"]["version"] == "1.0"
        assert data["jupyterlab"]["staticUrl"] == "bar"

    @pytest.mark.slow
    def test_build_custom_minimal_core_config(self):
        default_config = CoreConfig()
        core_config = CoreConfig()
        core_config.clear_packages()
        logger = logging.getLogger("jupyterlab_test_logger")
        logger.setLevel("DEBUG")
        app_dir = self.tempdir()
        options = AppOptions(
            app_dir=app_dir,
            core_config=core_config,
            logger=logger,
            use_sys_dir=False,
        )

        extensions = (
            "@jupyterlab/application-extension",
            "@jupyterlab/apputils-extension",
        )
        singletons = (
            "@jupyterlab/application",
            "@jupyterlab/apputils",
            "@jupyterlab/coreutils",
            "@jupyterlab/services",
        )
        for name in extensions:
            semver = default_config.extensions[name]
            core_config.add(name, semver, extension=True)
        for name in singletons:
            semver = default_config.singletons[name]
            core_config.add(name, semver)

        assert install_extension(self.mock_extension, app_options=options) is True
        build(app_options=options)

        # check static directory.
        entry = pjoin(app_dir, "static", "index.out.js")
        with open(entry) as fid:
            data = fid.read()
        assert self.pkg_names["extension"] in data

        pkg = pjoin(app_dir, "static", "package.json")
        with open(pkg) as fid:
            data = json.load(fid)
        assert sorted(data["jupyterlab"]["extensions"].keys()) == [
            "@jupyterlab/application-extension",
            "@jupyterlab/apputils-extension",
            "@jupyterlab/mock-extension",
        ]
        assert data["jupyterlab"]["mimeExtensions"] == {}
        for pkg in data["jupyterlab"]["singletonPackages"]:
            if pkg.startswith("@jupyterlab/"):
                assert pkg in singletons

    def test_disable_extension(self):
        options = AppOptions(app_dir=self.tempdir())
        assert install_extension(self.mock_extension, app_options=options) is True
        assert disable_extension(self.pkg_names["extension"], app_options=options) is True
        info = get_app_info(app_options=options)
        name = self.pkg_names["extension"]
        assert info["disabled"].get(name) is True
        assert not check_extension(name, app_options=options)
        assert check_extension(name, installed=True, app_options=options)
        assert disable_extension("@jupyterlab/notebook-extension", app_options=options) is True
        info = get_app_info(app_options=options)
        assert info["disabled"].get("@jupyterlab/notebook-extension") is True
        assert not check_extension("@jupyterlab/notebook-extension", app_options=options)
        assert check_extension(
            "@jupyterlab/notebook-extension", installed=True, app_options=options
        )
        assert info["disabled"].get(name) is True
        assert not check_extension(name, app_options=options)
        assert check_extension(name, installed=True, app_options=options)

    def test_enable_extension(self):
        options = AppOptions(app_dir=self.tempdir())
        assert install_extension(self.mock_extension, app_options=options) is True
        assert disable_extension(self.pkg_names["extension"], app_options=options) is True
        assert enable_extension(self.pkg_names["extension"], app_options=options) is True
        info = get_app_info(app_options=options)
        assert "@jupyterlab/notebook-extension" not in info["disabled"]
        name = self.pkg_names["extension"]
        assert info["disabled"].get(name, False) is False
        assert check_extension(name, app_options=options)
        assert disable_extension("@jupyterlab/notebook-extension", app_options=options) is True
        assert check_extension(name, app_options=options)
        assert not check_extension("@jupyterlab/notebook-extension", app_options=options)

    @pytest.mark.slow
    def test_build_check(self):
        # Do the initial build.
        assert build_check()
        assert install_extension(self.mock_extension) is True
        assert link_package(self.mock_package) is True
        build()
        assert not build_check()

        # Check installed extensions.
        assert install_extension(self.mock_mimeextension) is True
        assert build_check()
        assert uninstall_extension(self.pkg_names["mimeextension"]) is True
        assert not build_check()

        # Check local extensions.
        pkg_path = pjoin(self.mock_extension, "package.json")
        with open(pkg_path) as fid:
            data = json.load(fid)
        with open(pkg_path, "rb") as fid:
            orig = fid.read()
        data["foo"] = "bar"
        with open(pkg_path, "w") as fid:
            json.dump(data, fid)
        assert build_check()
        assert build_check()

        with open(pkg_path, "wb") as fid:
            fid.write(orig)

        assert not build_check()

        # Check linked packages.
        pkg_path = pjoin(self.mock_package, "index.js")
        with open(pkg_path, "rb") as fid:
            orig = fid.read()
        with open(pkg_path, "wb") as fid:
            fid.write(orig + b'\nconsole.log("hello");')
        assert build_check()
        assert build_check()

        with open(pkg_path, "wb") as fid:
            fid.write(orig)
        assert not build_check()

    def test_compatibility(self):
        assert _test_overlap("^0.6.0", "^0.6.1")
        assert _test_overlap(">0.1", "0.6")
        assert _test_overlap("~0.5.0", "~0.5.2")
        assert _test_overlap("0.5.2", "^0.5.0")

        assert not _test_overlap("^0.5.0", "^0.6.0")
        assert not _test_overlap("~1.5.0", "^1.6.0")

        assert _test_overlap("*", "0.6") is None
        assert _test_overlap("<0.6", "0.1") is None

        assert _test_overlap("^1 || ^2", "^1")
        assert _test_overlap("^1 || ^2", "^2")
        assert _test_overlap("^1", "^1 || ^2")
        assert _test_overlap("^2", "^1 || ^2")
        assert _test_overlap("^1 || ^2", "^2 || ^3")
        assert not _test_overlap("^1 || ^2", "^3 || ^4")
        assert not _test_overlap("^2", "^1 || ^3")

    def test_compare_ranges(self):
        assert _compare_ranges("^1 || ^2", "^1") == 0
        assert _compare_ranges("^1 || ^2", "^2 || ^3") == 0
        assert _compare_ranges("^1 || ^2", "^3 || ^4") == 1
        assert _compare_ranges("^3 || ^4", "^1 || ^2") == -1
        assert _compare_ranges("^2 || ^3", "^1 || ^4") is None

    def test_install_compatible(self):
        core_data = _get_default_core_data()
        current_app_dep = core_data["dependencies"]["@jupyterlab/application"]

        def _gen_dep(ver):
            return {"dependencies": {"@jupyterlab/application": ver}}

        def _mock_metadata(registry, name, logger):
            assert name == "mockextension"
            return {
                "name": name,
                "versions": {
                    "0.9.0": _gen_dep(current_app_dep),
                    "1.0.0": _gen_dep(current_app_dep),
                    "1.1.0": _gen_dep(current_app_dep),
                    "2.0.0": _gen_dep("^2000.0.0"),
                    "2.0.0-b0": _gen_dep(current_app_dep),
                    "2.1.0-b0": _gen_dep("^2000.0.0"),
                    "2.1.0": _gen_dep("^2000.0.0"),
                },
            }

        def _mock_extract(self, source, tempdir, *args, **kwargs):
            data = {
                "name": source,
                "version": "2.1.0",
                "jupyterlab": {"extension": True},
                "jupyterlab_extracted_files": ["index.js"],
            }
            data.update(_gen_dep("^2000.0.0"))
            info = {
                "source": source,
                "is_dir": False,
                "data": data,
                "name": source,
                "version": data["version"],
                "filename": "mockextension.tgz",
                "path": pjoin(tempdir, "mockextension.tgz"),
            }
            return info

        class Success(Exception):  # noqa
            pass

        def _mock_install(self, name, *args, **kwargs):
            assert name in ("mockextension", "mockextension@1.1.0")
            if name == "mockextension@1.1.0":
                raise Success()
            return orig_install(self, name, *args, **kwargs)

        p1 = patch.object(commands, "_fetch_package_metadata", _mock_metadata)
        p2 = patch.object(commands._AppHandler, "_extract_package", _mock_extract)
        p3 = patch.object(commands._AppHandler, "_install_extension", _mock_install)
        with p1, p2:
            orig_install = commands._AppHandler._install_extension
            with p3, pytest.raises(Success):
                assert install_extension("mockextension") is True

    def test_update_single(self):
        installed = []

        def _mock_install(self, name, *args, **kwargs):
            installed.append(name[0] + name[1:].split("@")[0])
            return {"name": name, "is_dir": False, "path": "foo/bar/" + name}

        def _mock_latest(self, name):
            return "10000.0.0"

        p1 = patch.object(commands._AppHandler, "_install_extension", _mock_install)
        p2 = patch.object(commands._AppHandler, "_latest_compatible_package_version", _mock_latest)

        assert install_extension(self.mock_extension) is True
        assert install_extension(self.mock_mimeextension) is True

        with p1, p2:
            assert update_extension(self.pkg_names["extension"]) is True
        assert installed == [self.pkg_names["extension"]]

    def test_update_missing_extension(self):
        assert update_extension("foo") is False

    def test_update_multiple(self):
        installed = []

        def _mock_install(self, name, *args, **kwargs):
            installed.append(name[0] + name[1:].split("@")[0])
            return {"name": name, "is_dir": False, "path": "foo/bar/" + name}

        def _mock_latest(self, name):
            return "10000.0.0"

        p1 = patch.object(commands._AppHandler, "_install_extension", _mock_install)
        p2 = patch.object(commands._AppHandler, "_latest_compatible_package_version", _mock_latest)

        install_extension(self.mock_extension)
        install_extension(self.mock_mimeextension)

        with p1, p2:
            assert update_extension(self.pkg_names["extension"]) is True
            assert update_extension(self.pkg_names["mimeextension"]) is True
        assert installed == [self.pkg_names["extension"], self.pkg_names["mimeextension"]]

    def test_update_all(self):
        updated = []

        def _mock_update(self, name, *args, **kwargs):
            updated.append(name[0] + name[1:].split("@")[0])
            return True

        original_app_info = commands._AppHandler._get_app_info

        def _mock_app_info(self):
            info = original_app_info(self)
            info["local_extensions"] = []
            return info

        assert install_extension(self.mock_extension) is True
        assert install_extension(self.mock_mimeextension) is True

        p1 = patch.object(commands._AppHandler, "_update_extension", _mock_update)

        # local packages are not updated, so mock them as non-local:
        p2 = patch.object(commands._AppHandler, "_get_app_info", _mock_app_info)

        with p1, p2:
            assert update_extension(None, all_=True) is True
        assert sorted(updated) == [self.pkg_names["extension"], self.pkg_names["mimeextension"]]


def test_load_extension(jp_serverapp, make_lab_app):
    app = make_lab_app()
    stderr = sys.stderr
    #    sys.stderr = self.devnull
    app._link_jupyter_server_extension(jp_serverapp)
    app.initialize()
    sys.stderr = stderr
