import email.utils
import json
import os
import shutil
import tarfile
import tempfile
import time
import zipfile
from collections.abc import Generator
from pathlib import Path
from typing import Any

from traitlets import Bool, Instance
from traitlets.config import LoggingConfigurable

from ..constants import (
    DISABLED_EXTENSIONS,
    EXTENSION_TAR,
    EXTENSION_ZIP,
    FEDERATED_EXTENSIONS,
    JSON_FMT,
    JUPYTER_CONFIG_DATA,
    JUPYTERLITE_IPYNB,
    JUPYTERLITE_JSON,
    JUPYTERLITE_METADATA,
    LITE_PLUGIN_SETTINGS,
    MOD_DIRECTORY,
    MOD_FILE,
    SETTINGS_OVERRIDES,
    SOURCEMAP_IGNORE_PATTERNS,
    SOURCEMAPS,
    UTF8,
)
from ..manager import LiteManager
from ..optional import has_optional_dependency


class BaseAddon(LoggingConfigurable):
    """A base class for addons to the JupyterLite build chain

    Subclassing this is optional, but provides some useful utilities
    """

    manager: LiteManager = Instance(LiteManager)

    ignore_sys_prefix: bool = Bool(False)

    def __init__(self, manager, *args, **kwargs):
        kwargs["parent"] = manager
        kwargs["manager"] = manager
        super().__init__(*args, **kwargs)

    @property
    def log(self):
        return self.manager.log

    def task(self, **task):
        """Ensure a ``doit`` task is well-formed.

        - Task names may not include the ``=`` character
          - arbitrarily replace with ``--``
        """
        task["name"] = task["name"].replace("=", "--")
        return task

    def copy_one(self, src, dest):
        """copy one Path (a file or folder)"""
        if self.manager.no_sourcemaps and self.is_ignored_sourcemap(src.name):
            return

        if dest.is_dir():
            shutil.rmtree(dest)
        elif dest.exists():
            dest.unlink()

        if not dest.parent.exists():
            self.log.debug(f"creating folder {dest.parent}")
            dest.parent.mkdir(parents=True)

        self.maybe_timestamp(dest.parent)

        copytree_kwargs = {}

        if self.manager.no_sourcemaps:
            copytree_kwargs["ignore"] = SOURCEMAP_IGNORE_PATTERNS

        if src.is_dir():
            shutil.copytree(src, dest, **copytree_kwargs)
        else:
            shutil.copy2(src, dest)

        self.maybe_timestamp(dest)

    def fetch_one(self, url, dest):
        """fetch one file

        TODO: enable other backends, auth, etc.
        """
        import urllib.request

        if dest.exists():
            self.log.info(f"[lite][fetch] already downloaded {dest.name}, skipping...")
            return

        if not dest.parent.exists():
            dest.parent.mkdir(parents=True)

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            with urllib.request.urlopen(url) as response:  # noqa: S310
                tmp_dest = tdp / dest.name
                with tmp_dest.open("wb") as fd:
                    shutil.copyfileobj(response, fd)
                last_modified = response.headers.get("Last-Modified")
                if last_modified:
                    epoch_time = time.mktime(email.utils.parsedate(last_modified))
                    os.utime(tmp_dest, (epoch_time, epoch_time))
            shutil.copy2(tmp_dest, dest)

    def maybe_timestamp(self, path):
        if not path.exists() or self.manager.source_date_epoch is None:
            return

        if path.is_dir():
            for p in path.rglob("*"):
                self.timestamp_one(p)

        self.timestamp_one(path)

    def timestamp_one(self, path):
        """adjust the timestamp to be --source-date-epoch for files newer than then

        see https://reproducible-builds.org/specs/source-date-epoch
        """
        stat = path.stat()
        sde = self.manager.source_date_epoch
        if stat.st_mtime > sde:
            cls = self.__class__.__name__
            self.log.debug(f"[lite][base] <{cls}> set time to source_date_epoch {sde} on {path}")
            os.utime(path, (sde, sde))
            return
        return

    def delete_one(self, *src):
        """delete... somethings"""
        for src_dir in src:
            if src_dir.is_dir():
                shutil.rmtree(src_dir)
            elif src_dir.exists():
                src_dir.unlink()

    def validate_one_json_file(self, validator, path=None, data=None, selector=None):
        loaded = json.loads(path.read_text(**UTF8)) if path else data

        if selector:
            for sel in selector:
                selected = loaded.get(sel, {})
        else:
            selected = loaded

        if validator is None:
            # just checking if the JSON is well-formed
            return True

        if isinstance(validator, Path):
            validator = self.get_validator(validator)
            if validator is None:
                return True

        validator.validate(selected)

    def get_validator(self, schema_path, klass=None):
        if klass is None:
            if not has_optional_dependency(
                "jsonschema", "only checking JSON well-formedness: {error}"
            ):
                return None
            from jsonschema import Draft7Validator

            klass = Draft7Validator

        schema = json.loads(schema_path.read_text(**UTF8))
        return klass(schema)

    def merge_one_jupyterlite(self, out_path, in_paths):  # noqa: C901, PLR0912
        """write the ``out_path`` with the merge content of ``in_paths``, where
        all are valid ``jupyter-lite.*`` files.
        """
        self.log.debug(f"[lite][config][merge] {out_path}")
        config = None

        for in_path in in_paths:
            self.log.debug(f"[lite][config][merge] . {in_path}")
            in_config = None
            try:
                in_config = json.loads(in_path.read_text(**UTF8))
                if out_path.name == JUPYTERLITE_IPYNB:
                    in_config = in_config["metadata"].get(JUPYTERLITE_METADATA)
            except:  # noqa: E722, S110
                pass

            if not in_config:
                continue

            if not config:
                config = in_config
                continue

            for k, v in in_config.items():
                self.log.debug(f"""[lite][config] ... updating {k} => {v}?""")
                if k == JUPYTER_CONFIG_DATA:
                    config[k] = self.merge_jupyter_config_data(config.get(k) or {}, v)
                elif config.get(k) != v:
                    self.log.debug(f"""[lite][config] ..... {k} updated""")
                    config[k] = v

        if config and JUPYTER_CONFIG_DATA in config:
            self.dedupe_federated_extensions(config[JUPYTER_CONFIG_DATA])

        if out_path.name == JUPYTERLITE_IPYNB:
            if out_path.exists():
                doc_path = out_path
            else:
                for in_path in in_paths:
                    if in_path.name == JUPYTERLITE_IPYNB and in_path.exists():
                        doc_path = in_path
                        break

            doc = json.loads(doc_path.read_text(**UTF8))

            doc["metadata"][JUPYTERLITE_METADATA] = config

            out_path.write_text(json.dumps(doc, **JSON_FMT), **UTF8)
        else:
            out_path.write_text(json.dumps(config, **JSON_FMT), **UTF8)

        print("MERGED", out_path, "from", in_paths)

    def merge_jupyter_config_data(self, config, in_config):
        """merge well-known ``jupyter-config-data`` fields"""
        self.log.debug(f"""[lite][config][merge] ..... {config}""")
        self.log.debug(f"""[lite][config][merge] ..... {in_config}""")

        config = config or {}
        in_config = in_config or {}

        for k, v in in_config.items():
            if k in [DISABLED_EXTENSIONS, FEDERATED_EXTENSIONS]:
                config[k] = [*config.get(k, []), *v]
            elif k in [SETTINGS_OVERRIDES]:
                config[k] = config.get(k, {})
                for pkg, pkg_config in v.items():
                    config[k][pkg] = config[k].get(pkg, {})
                    config[k][pkg].update(pkg_config)
            else:
                config[k] = v
        self.log.debug(f"""[lite][config][merge] ..... {config}""")
        return config

    def dedupe_federated_extensions(self, config):
        """update a federated_extension list in-place, ensuring unique names.

        .. todo::

            best we can do, for now.
        """
        if FEDERATED_EXTENSIONS not in config:
            return

        named = {}

        for ext in config[FEDERATED_EXTENSIONS]:
            named[ext["name"]] = ext

        config[FEDERATED_EXTENSIONS] = sorted(named.values(), key=lambda x: x["name"])

    def is_ignored_sourcemap(self, path_name: str):
        is_ignored = False
        if self.manager.no_sourcemaps:
            for map_ext in SOURCEMAPS:
                if path_name.endswith(map_ext):
                    is_ignored = True
                    break
        return is_ignored

    def is_sys_prefix_ignored(self):
        return self.ignore_sys_prefix

    @property
    def should_use_libarchive_c(self):
        """should libarchive-c be used (if available)?"""
        if self.manager.no_libarchive:
            return False

        return has_optional_dependency(
            "libarchive",
            "install libarchive-c for better performance when working with archives: {error}",
        )

    def extract_one(self, archive: Path, dest: Path):
        """extract the contents of an archive to a path."""
        if dest.exists():
            shutil.rmtree(dest)

        dest.mkdir(parents=True)

        if self.should_use_libarchive_c:
            import libarchive

            old_cwd = os.getcwd()
            os.chdir(str(dest))
            try:
                libarchive.extract_file(str(archive))
            finally:
                os.chdir(old_cwd)
            return

        if archive.name.endswith(EXTENSION_ZIP):
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(dest)  # noqa: S202
        elif archive.name.endswith(EXTENSION_TAR):
            mode = "r:bz2" if archive.name.endswith(".bz2") else "r:gz"
            with tarfile.open(archive, mode) as tf:
                self.safe_extract_all(tf, dest)
        else:
            raise ValueError(f"Unrecognized archive format {archive.name}")

        for path in [dest, *dest.rglob("*")]:
            path.chmod(MOD_DIRECTORY if path.is_dir() else MOD_FILE)

    def is_within_directory(self, directory, target):
        abs_directory = os.path.abspath(directory)
        abs_target = os.path.abspath(target)
        prefix = os.path.commonprefix([abs_directory, abs_target])
        return prefix == abs_directory

    def safe_extract_all(self, tar, path=".", members=None, *, numeric_owner=False):
        for member in tar.getmembers():
            member_path = os.path.join(path, member.name)
            if not self.is_within_directory(path, member_path):
                raise Exception("Attempted Path Traversal in Tar File")
        tar.extractall(path, members, numeric_owner=numeric_owner)  # noqa: S202

    def hash_all(self, hashfile: Path, root: Path, paths: list[Path]):
        from hashlib import sha256

        lines = [
            "  ".join([sha256(p.read_bytes()).hexdigest(), p.relative_to(root).as_posix()])
            for p in sorted(paths)
        ]
        hashfile.write_text("\n".join(lines))

    def get_lite_config_paths(self) -> Generator[Path, None, None]:
        """Yield all config paths that exist in the ``lite_dir``."""
        for app in [None, *self.manager.apps]:
            app_dir = self.manager.lite_dir / app if app else self.manager.output_dir
            for path_name in [JUPYTERLITE_JSON, JUPYTERLITE_IPYNB]:
                config_path = app_dir / path_name
                if config_path.exists():
                    yield config_path

    def get_output_config_paths(self) -> Generator[Path, None, None]:
        """Yield all config paths that _might_ exist in the ``output_dir``."""
        for app in [None, *self.manager.apps]:
            app_dir = self.manager.output_dir / app if app else self.manager.output_dir
            for path_name in [JUPYTERLITE_JSON, JUPYTERLITE_IPYNB]:
                config_path = app_dir / path_name
                yield config_path

    def get_lite_plugin_settings(self, config_path: Path, plugin_id: str) -> dict[str, Any]:
        """Get the plugin settings from a config path.

        The keys follow the JupyterLab settings naming convention, of module and
        identifier e.g.

            @jupyterlite/contents:plugin
        """
        if not config_path.exists():
            return {}

        config = json.loads(config_path.read_text(**UTF8))

        # if a notebook, look in the top-level metadata (which must exist)
        if config_path.name == JUPYTERLITE_IPYNB:
            config = config["metadata"].get(JUPYTERLITE_METADATA, {})

        return config.get(JUPYTER_CONFIG_DATA, {}).get(LITE_PLUGIN_SETTINGS, {}).get(plugin_id, {})

    def set_lite_plugin_settings(
        self, config_path: Path, plugin_id: str, settings: dict[str, Any]
    ) -> None:
        """Overwrite the plugin settings for a single plugin in a config path."""
        whole_file = config = json.loads(config_path.read_text(**UTF8))
        if config_path.name == JUPYTERLITE_IPYNB:
            config = whole_file["metadata"][JUPYTERLITE_METADATA]

        config.setdefault(JUPYTER_CONFIG_DATA, {}).setdefault(LITE_PLUGIN_SETTINGS, {}).update(
            {plugin_id: settings}
        )

        config_path.write_text(json.dumps(whole_file, **JSON_FMT), **UTF8)
        self.log.debug("%s wrote settings in %s: %s", plugin_id, config_path, settings)
        self.maybe_timestamp(config_path)
