"""a JupyterLite addon for supporting extension settings"""

import json

from ..constants import (
    JSON_FMT,
    JUPYTER_CONFIG_DATA,
    JUPYTERLITE_IPYNB,
    JUPYTERLITE_JSON,
    JUPYTERLITE_METADATA,
    LAB_EXTENSIONS,
    OVERRIDES_JSON,
    SETTINGS_OVERRIDES,
    UTF8,
)
from .base import BaseAddon


class SettingsAddon(BaseAddon):
    """sync the as-installed `overrides.json` and update `jupyter-lite.json`"""

    __all__ = ["pre_build", "post_build", "check", "status"]

    def status(self, manager):
        apps = []
        for app in [None, *manager.apps]:
            app_dir = manager.lite_dir / app if app else manager.lite_dir
            overrides_json = app_dir / OVERRIDES_JSON
            if not overrides_json.exists():
                continue
            apps.append(overrides_json)

        yield self.task(
            name="overrides",
            actions=[lambda: print(f"""    {OVERRIDES_JSON}: {len(apps)}""")],
        )

    def pre_build(self, manager):
        for app in [None, *manager.apps]:
            app_dir = manager.lite_dir / app if app else manager.lite_dir
            app_output_dir = manager.output_dir / app if app else manager.output_dir
            overrides_json = app_dir / OVERRIDES_JSON
            if not overrides_json.exists():
                continue
            dest = app_output_dir / overrides_json.name
            yield self.task(
                name=f"""copy:{app or "root"}/""",
                file_dep=[overrides_json],
                targets=[dest],
                actions=[(self.copy_one, [overrides_json, dest])],
            )

    def post_build(self, manager):
        """add settings from `overrides.json` per app where they exist"""
        for app in [None, *manager.apps]:
            app_output_dir = manager.output_dir / app if app else manager.output_dir
            jupyterlite_json = app_output_dir / JUPYTERLITE_JSON
            overrides_json = app_output_dir / OVERRIDES_JSON
            if not overrides_json.exists():
                continue

            yield self.task(
                name=f"patch:overrides:{app}",
                file_dep=[overrides_json, jupyterlite_json],
                actions=[(self.patch_one_overrides, [jupyterlite_json, overrides_json])],
            )

    def check(self, manager):
        for lite_file in [
            *manager.output_dir.rglob(JUPYTERLITE_JSON),
            *manager.output_dir.rglob(JUPYTERLITE_IPYNB),
        ]:
            yield from self.check_one_lite_file(lite_file)

    def check_one_lite_file(self, lite_file):
        config = json.loads(lite_file.read_text(**UTF8))

        if lite_file.name == JUPYTERLITE_IPYNB:
            config = config["metadata"][JUPYTERLITE_METADATA]

        overrides = config.get(JUPYTER_CONFIG_DATA, {}).get(SETTINGS_OVERRIDES, {})

        for plugin_id, defaults in overrides.items():
            ext, plugin = plugin_id.split(":")
            plugin_stem = f"schemas/{ext}/{plugin}.json"
            schema = self.output_extensions / ext / plugin_stem
            if not schema.exists():
                core_schema = self.manager.output_dir / "build" / plugin_stem
                if core_schema.exists():
                    schema = core_schema
                else:
                    self.log.debug(f"[lite] [settings] Missing {plugin} (probably harmless)")
                    continue

            validator = self.get_validator(schema)

            yield self.task(
                name=f"overrides:{plugin_id}",
                file_dep=[lite_file, schema],
                actions=[(self.validate_one_json_file, [validator, None, defaults])],
            )

    def patch_one_overrides(self, jupyterlite_json, overrides_json):
        """update and normalize settingsOverrides"""
        try:
            config = json.loads(jupyterlite_json.read_text(**UTF8))
        except (FileNotFoundError, json.JSONDecodeError):
            self.log.debug(f"[lite] [settings] Initializing {jupyterlite_json}")
            config = {JUPYTER_CONFIG_DATA: {}}

        overrides = config[JUPYTER_CONFIG_DATA].get(SETTINGS_OVERRIDES, {})

        from_json = json.loads(overrides_json.read_text(**UTF8))
        for k, v in from_json.items():
            if k in overrides:
                overrides[k].update(v)
            else:
                overrides[k] = v

        config[JUPYTER_CONFIG_DATA][SETTINGS_OVERRIDES] = overrides

        jupyterlite_json.write_text(json.dumps(config, **JSON_FMT), **UTF8)

        self.maybe_timestamp(jupyterlite_json)
        self.log.debug(f"[lite] [settings] Updated {jupyterlite_json}")

    @property
    def output_extensions(self):
        """where labextensions will go in the output folder"""
        return self.manager.output_dir / LAB_EXTENSIONS
