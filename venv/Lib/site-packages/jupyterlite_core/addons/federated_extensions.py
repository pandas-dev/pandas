"""a JupyterLite addon for supporting federated_extensions"""

import json
import re
import sys
import urllib.parse
from pathlib import Path

from traitlets import List, Unicode

from ..constants import (
    ALL_FEDERATED_JSON,
    FEDERATED_EXTENSIONS,
    JSON_FMT,
    JUPYTER_CONFIG_DATA,
    JUPYTERLITE_JSON,
    LAB_EXTENSIONS,
    PACKAGE_JSON,
    SHA256SUMS,
    SHARE_LABEXTENSIONS,
    UTF8,
)
from .base import BaseAddon


class FederatedExtensionAddon(BaseAddon):
    """sync the as-installed federated_extensions and update `jupyter-lite.json`"""

    __all__ = ["pre_build", "post_build", "post_init"]

    labextensions_path = Path(sys.prefix) / SHARE_LABEXTENSIONS

    extra_labextensions_path = List(
        Unicode(), help="""Extra paths to look for federated JupyterLab extensions"""
    ).tag(config=True)

    def env_extensions(self, root):
        """a list of all federated extensions"""
        return [
            p
            for p in [
                *root.glob(f"*/{PACKAGE_JSON}"),
                *root.glob(f"@*/*/{PACKAGE_JSON}"),
            ]
            if self.is_prebuilt(json.loads(p.read_text(**UTF8)))
        ]

    @property
    def ext_cache(self):
        """where extensions will go in the cache"""
        return self.manager.cache_dir / "federated_extensions"

    @property
    def archive_cache(self):
        """where archives will go in the cache"""
        return self.manager.cache_dir / "archives"

    @property
    def output_extensions(self):
        """where labextensions will go in the output folder"""
        return self.manager.output_dir / LAB_EXTENSIONS

    def post_init(self, manager):
        """handle downloading of federated extensions"""
        for path_or_url in manager.federated_extensions:
            yield from self.resolve_one_extension(path_or_url, init=True)

    def pre_build(self, manager):
        """yield a doit task to copy each federated extension into the output_dir"""
        if not self.is_sys_prefix_ignored():
            for pkg_json in self.env_extensions(self.labextensions_path):
                yield from self.copy_one_extension(pkg_json)

        for p in self.extra_labextensions_path:
            for pkg_json in self.env_extensions(Path(p)):
                yield from self.copy_one_extension(pkg_json)

        for path_or_url in manager.federated_extensions:
            yield from self.resolve_one_extension(path_or_url, init=False)

    def build(self, manager):
        """yield a doit task to copy each local extension into the output_dir"""
        root = self.manager.lite_dir / LAB_EXTENSIONS

        for pkg_json in self.env_extensions(root):
            yield self.copy_one_extension(pkg_json)

    def copy_one_env_extension(self, pkg_json):
        """yield tasks to copy one unpacked on-disk extension from sys.prefix into the output dir"""
        yield from self.copy_one_extension(pkg_json)

    def copy_one_extension(self, pkg_json):
        """yield a task to copy one unpacked on-disk extension from anywhere into the output dir"""
        pkg_path = pkg_json.parent
        stem = json.loads(pkg_json.read_text(**UTF8))["name"]
        dest = self.output_extensions / stem
        file_dep = [
            p for p in pkg_path.rglob("*") if not (p.is_dir() or self.is_ignored_sourcemap(p.name))
        ]
        targets = [dest / p.relative_to(pkg_path) for p in file_dep]

        yield self.task(
            name=f"copy:ext:{stem}",
            file_dep=file_dep,
            targets=targets,
            actions=[(self.copy_one, [pkg_path, dest])],
        )

    def resolve_one_extension(self, path_or_url, init):
        """yield tasks try to resolve one URL or local folder/archive
        as a (set of) federated_extension(s)"""
        if re.findall(r"^https?://", path_or_url):
            url = urllib.parse.urlparse(path_or_url)
            name = url.path.split("/")[-1]
            dest = self.ext_cache / name
            if init:
                if not dest.exists():
                    yield self.task(
                        name=f"fetch:{name}",
                        actions=[(self.fetch_one, [path_or_url, dest])],
                        targets=[dest],
                    )
                return
            # if not initializing, assume path is now local
            path_or_url = dest.resolve()

        if init:
            # nothing to do for local files during this phase
            return

        local_path = (self.manager.lite_dir / path_or_url).resolve()

        if local_path.is_dir():
            yield from self.copy_one_folder_extension(local_path)
        elif local_path.exists():
            suffix = local_path.suffix

            if local_path.name.endswith((".whl", ".tar.bz2")):
                yield from self.copy_simple_archive_extensions(local_path)
            elif local_path.name.endswith(".conda"):
                yield from self.copy_conda2_extensions(local_path)
            else:  # pragma: no cover
                raise NotImplementedError(f"unknown archive {suffix}")
        else:  # pragma: no cover
            raise FileNotFoundError(path_or_url)

    def copy_one_folder_extension(self, path):
        """yield a task to copy one extension from the given path"""
        pkg_json = path / PACKAGE_JSON

        if not pkg_json.exists():
            raise ValueError(f"[lite][federated_extensions] No package.json in {path}")

        yield from self.copy_one_extension(pkg_json)

    def copy_simple_archive_extensions(self, archive: Path):
        """yield tasks to extract and copy the labextensions from a local simple archive"""
        unarchived = self.archive_cache / archive.name
        hashfile = self.archive_cache / f"{archive.name}.{SHA256SUMS}"

        yield dict(
            name=f"extract:{archive.name}",
            actions=[
                (self.delete_one, [hashfile]),
                (self.extract_one, [archive, unarchived]),
                lambda: self.hash_all(
                    hashfile,
                    unarchived,
                    [p for p in unarchived.rglob("*") if not p.is_dir()],
                ),
            ],
            file_dep=[archive],
        )

        yield dict(
            name=f"copy:{archive.name}",
            file_dep=[hashfile],
            actions=[(self.copy_all_federated_extensions, [unarchived])],
        )

    def copy_all_federated_extensions(self, unarchived):
        """actually copy all federated extensions found in a folder."""
        for simple_pkg_json in unarchived.rglob(f"{SHARE_LABEXTENSIONS}/*/package.json"):
            self.copy_one_federated_extension(simple_pkg_json)
        for org_pkg_json in unarchived.rglob(f"{SHARE_LABEXTENSIONS}/@*/*/package.json"):
            self.copy_one_federated_extension(org_pkg_json)

    def copy_one_federated_extension(self, pkg_json):
        """actually copy one labextension from an extracted archive"""
        pkg_data = json.loads(pkg_json.read_text(**UTF8))

        if self.is_prebuilt(pkg_data):
            pkg_name = pkg_data["name"]
            output_pkg = self.output_extensions / pkg_name
            self.copy_one(pkg_json.parent, output_pkg)

    def is_prebuilt(self, pkg_json):
        """verify this is an actual pre-built extension, containing load information"""
        return pkg_json.get("jupyterlab", {}).get("_build", {}).get("load") is not None

    def copy_conda2_extensions(self, conda_pkg):
        """copy the labextensions from a local, nested ``.conda`` package"""
        if self.manager.no_libarchive:
            raise RuntimeError(
                "`.conda` packages are not supported by python's stdlib. Please:\n\n"
                "\tconda install python-libarchive-c\n\nor:\n\n"
                "\tpip install libarchive-c"
            )
        unarchived = self.archive_cache / conda_pkg.name
        inner_archive = unarchived / f"pkg-{conda_pkg.stem}.tar.zst"
        inner_unarchived = self.archive_cache / inner_archive.name
        hashfile = self.archive_cache / f"{inner_archive.name}.{SHA256SUMS}"

        yield dict(
            name=f"extract:{conda_pkg.name}",
            actions=[
                (self.extract_one, [conda_pkg, unarchived]),
            ],
            file_dep=[conda_pkg],
            targets=[inner_archive],
        )

        yield dict(
            name=f"extract:{inner_archive.name}",
            actions=[
                (self.delete_one, [hashfile]),
                (self.extract_one, [inner_archive, inner_unarchived]),
                lambda: self.hash_all(
                    hashfile,
                    inner_unarchived,
                    [p for p in inner_unarchived.rglob("*") if not p.is_dir()],
                ),
            ],
            file_dep=[inner_archive],
        )

        yield dict(
            name=f"copy:{inner_archive.name}",
            file_dep=[hashfile],
            actions=[(self.copy_all_federated_extensions, [inner_unarchived])],
        )

    def post_build(self, manager):
        """update the root jupyter-lite.json, and copy each output theme to each app

        .. todo::

            the latter per-app steps should be at least cut in half, if not
            avoided altogether.
            See https://github.com/jupyterlite/jupyterlite/issues/118
        """
        jupyterlite_json = manager.output_dir / JUPYTERLITE_JSON
        lab_extensions_root = manager.output_dir / LAB_EXTENSIONS
        lab_extensions = self.env_extensions(lab_extensions_root)

        yield self.task(
            name="patch",
            doc=f"ensure {JUPYTERLITE_JSON} includes the federated_extensions",
            file_dep=[*lab_extensions, jupyterlite_json],
            actions=[(self.patch_jupyterlite_json, [jupyterlite_json])],
        )

        stems = [p.parent.relative_to(lab_extensions_root) for p in lab_extensions]

        app_themes = manager.output_dir / "build/themes"
        for stem in stems:
            pkg = lab_extensions_root / stem
            # this pattern appears to be canonical
            theme_dir = pkg / "themes" / stem
            if not theme_dir.is_dir():
                continue
            # this may be a package or an @org/package... same result
            file_dep = sorted(
                [
                    p
                    for p in theme_dir.rglob("*")
                    if not (p.is_dir() or self.is_ignored_sourcemap(p.name))
                ]
            )
            dest = app_themes / stem
            targets = [dest / p.relative_to(theme_dir) for p in file_dep]
            yield self.task(
                name=f"copy:theme:{stem}",
                doc=f"copy theme asset for {pkg}",
                file_dep=file_dep,
                targets=targets,
                actions=[(self.copy_one, [theme_dir, dest])],
            )

        app_schemas = manager.output_dir / "build" / "schemas"
        all_federated_json = app_schemas / ALL_FEDERATED_JSON

        if app_schemas.is_dir():
            yield self.task(
                name="settings",
                doc=f"ensure {ALL_FEDERATED_JSON} includes the settings of federated extensions",
                file_dep=[*lab_extensions],
                targets=[all_federated_json],
                actions=[
                    (self.ensure_federated_settings, [manager, lab_extensions, all_federated_json])
                ],
            )

    def ensure_federated_settings(self, manager, lab_extensions, all_federated_json):
        """ensure settings from federated extensions are aggregated in a single file"""
        all_federated_settings = [
            setting for p in lab_extensions for setting in self.get_federated_settings(p.parent)
        ]
        all_federated_json.write_text(json.dumps(all_federated_settings), **UTF8)

    def get_federated_settings(self, extension):
        """get the settings for a federated extension"""
        pkg_json = extension / PACKAGE_JSON
        pkg_data = json.loads(pkg_json.read_text(**UTF8))
        settings_dir = extension / "schemas"
        if not settings_dir.is_dir():
            # bail if there is no settings for that extension
            return []
        setting_files = sorted(settings_dir.rglob("*.json"))

        pkg_name = pkg_data["name"]
        pkg_version = pkg_data["version"]
        all_settings = []
        for setting_file in setting_files:
            plugin_id = f"{pkg_name}:{setting_file.stem}"
            schema = json.loads(setting_file.read_text(**UTF8))
            setting = {
                "id": plugin_id,
                "raw": "{}",
                "schema": schema,
                "settings": {},
                "version": pkg_version,
            }
            all_settings.append(setting)

        return all_settings

    def patch_jupyterlite_json(self, jupyterlite_json):
        """add the federated_extensions to jupyter-lite.json

        .. todo::

            it _really_ doesn't like duplicate ids, probably need to catch it
            earlier... not possible with "pure" schema (but perhaps SHACL?)
        """
        config = json.loads(jupyterlite_json.read_text(**UTF8))

        extensions = config[JUPYTER_CONFIG_DATA].get(FEDERATED_EXTENSIONS, [])
        lab_extensions_root = self.manager.output_dir / LAB_EXTENSIONS

        for pkg_json in self.env_extensions(lab_extensions_root):
            pkg_data = json.loads(pkg_json.read_text(**UTF8))
            is_lite = pkg_data.get("jupyterlite", {}).get("liteExtension", False)
            extension_data = {
                **pkg_data["jupyterlab"]["_build"],
                "liteExtension": is_lite,
            }
            extensions += [dict(name=pkg_data["name"], **extension_data)]

        self.dedupe_federated_extensions(config[JUPYTER_CONFIG_DATA])

        jupyterlite_json.write_text(json.dumps(config, **JSON_FMT), **UTF8)

        self.maybe_timestamp(jupyterlite_json)
