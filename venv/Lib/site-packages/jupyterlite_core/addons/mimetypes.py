"""a JupyterLite addon for customizing mime types"""

import json

import doit

from ..constants import (
    JSON_FMT,
    JUPYTER_CONFIG_DATA,
    JUPYTERLITE_JSON,
    SETTINGS_FILE_TYPES,
    UTF8,
)
from .base import BaseAddon


class MimetypesAddon(BaseAddon):
    """Handle custom MIME types."""

    __all__ = ["post_build", "status"]

    def status(self, manager):
        """Yield status about file types."""
        yield self.task(
            name=JUPYTERLITE_JSON,
            actions=[
                lambda: print(f"""    filetypes:         {len(self.file_types)} """),
            ],
        )

    @property
    def file_types(self):
        """A merged view of all configured file types."""
        file_types = dict()
        file_types.update(self.manager.file_types)
        file_types.update(self.manager.extra_file_types)
        return file_types

    def post_build(self, manager):
        """Yield ``doit`` tasks to update with file type config."""
        jupyterlite_json = manager.output_dir / JUPYTERLITE_JSON

        yield self.task(
            name="patch",
            uptodate=[
                doit.tools.config_changed(
                    dict(
                        file_types=self.manager.file_types,
                        extra_file_types=self.manager.file_types,
                    )
                )
            ],
            doc=f"ensure {jupyterlite_json} includes the file_types",
            file_dep=[jupyterlite_json],
            actions=[(self.patch_jupyterlite_json, [jupyterlite_json])],
        )

    def patch_jupyterlite_json(self, jupyterlite_json):
        """add the file_types to the base"""
        config = json.loads(jupyterlite_json.read_text(**UTF8))

        file_types = config[JUPYTER_CONFIG_DATA].get(SETTINGS_FILE_TYPES, {})
        file_types.update(self.file_types)

        config = json.loads(jupyterlite_json.read_text(**UTF8))
        config[JUPYTER_CONFIG_DATA][SETTINGS_FILE_TYPES] = file_types
        jupyterlite_json.write_text(json.dumps(config, **JSON_FMT), **UTF8)

        self.maybe_timestamp(jupyterlite_json)
