"""common addon features for ``jupyterlite-pyodide-kernel``

This should not be considered part of the public API, and much will disappear
when these features are added upstream:

    https://github.com/jupyterlite/jupyterlite/issues/996
"""

import json
from pathlib import Path
from typing import Generator, Dict, Any
from jupyterlite_core.addons.base import BaseAddon
from jupyterlite_core.constants import (
    JUPYTERLITE_IPYNB,
    JUPYTERLITE_JSON,
    UTF8,
    JUPYTERLITE_METADATA,
    JUPYTER_CONFIG_DATA,
    LITE_PLUGIN_SETTINGS,
    JSON_FMT,
)

from ..constants import PYODIDE_KERNEL_PLUGIN_ID

__all__ = ["_BaseAddon"]


class _BaseAddon(BaseAddon):
    def get_pyodide_settings(self, config_path: Path):
        """Get the settings for the client-side Pyodide kernel."""
        return self.get_lite_plugin_settings(config_path, PYODIDE_KERNEL_PLUGIN_ID)

    def set_pyodide_settings(self, config_path: Path, settings: Dict[str, Any]) -> None:
        """Update the settings for the client-side Pyodide kernel."""
        return self.set_lite_plugin_settings(
            config_path, PYODIDE_KERNEL_PLUGIN_ID, settings
        )

    def get_output_config_paths(self) -> Generator[Path, None, None]:
        """Yield an iterator of all config paths that _might_ exist in the
        ``output_dir``.

        This will likely move upstream.
        """
        for app in [None, *self.manager.apps]:
            app_dir = self.manager.output_dir / app if app else self.manager.output_dir
            for path_name in [JUPYTERLITE_JSON, JUPYTERLITE_IPYNB]:
                config_path = app_dir / path_name
                yield config_path

    def get_lite_plugin_settings(
        self, config_path: Path, plugin_id: str
    ) -> Dict[str, Any]:
        """Get the plugin settings from a config path.

        The keys follow the JupyterLab settings naming convention, of module and
        identifier e.g.

            @jupyterlite/contents:plugin

        This will likely move upstream.
        """
        if not config_path.exists():
            return {}

        config = json.loads(config_path.read_text(**UTF8))

        # if a notebook, look in the top-level metadata (which must exist)
        if config_path.name == JUPYTERLITE_IPYNB:
            config = config["metadata"].get(JUPYTERLITE_METADATA, {})

        return (
            config.get(JUPYTER_CONFIG_DATA, {})
            .get(LITE_PLUGIN_SETTINGS, {})
            .get(plugin_id, {})
        )

    def set_lite_plugin_settings(
        self, config_path: Path, plugin_id: str, settings: Dict[str, Any]
    ) -> None:
        """Overwrite the plugin settings for a single plugin in a config path.

        This will likely move upstream.
        """
        whole_file = config = json.loads(config_path.read_text(**UTF8))
        if config_path.name == JUPYTERLITE_IPYNB:
            config = whole_file["metadata"][JUPYTERLITE_METADATA]

        config.setdefault(JUPYTER_CONFIG_DATA, {}).setdefault(
            LITE_PLUGIN_SETTINGS, {}
        ).update({plugin_id: settings})

        config_path.write_text(json.dumps(whole_file, **JSON_FMT), **UTF8)
        self.log.debug("%s wrote settings in %s: %s", plugin_id, config_path, settings)
        self.maybe_timestamp(config_path)
