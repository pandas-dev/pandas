"""a JupyterLite addon for jupyterlite-specific tasks"""

from ..constants import (
    JUPYTERLITE_IPYNB,
    JUPYTERLITE_JSON,
    JUPYTERLITE_METADATA,
    JUPYTERLITE_SCHEMA,
)
from .base import BaseAddon


class LiteAddon(BaseAddon):
    """ensure jupyterlite files have been merged, and validate them"""

    __all__ = ["build", "check", "status"]

    def status(self, manager):
        yield self.task(
            name=JUPYTERLITE_JSON,
            actions=[
                lambda: self.log.debug(f"""    jupyter-lite.(json|ipynb): {self.lite_files}"""),
                lambda: self.log.info(
                    f"""    jupyter-lite.(json|ipynb): {len(self.lite_files)} files"""
                ),
            ],
        )

    def build(self, manager):
        """merge jupyter-lite.json into the output_dir"""
        lite_dir = manager.lite_dir
        output_dir = manager.output_dir

        for jupyterlite_file in self.lite_files:
            rel = jupyterlite_file.relative_to(lite_dir)
            dest = output_dir / rel
            yield self.task(
                name=f"patch:{rel}",
                file_dep=[jupyterlite_file],
                actions=[
                    (self.merge_one_jupyterlite, [dest, [dest, jupyterlite_file]]),
                    (self.maybe_timestamp, [dest]),
                ],
            )

    def check(self, manager):
        """apply schema validation to all `jupyter-lite.json` in the `output_dir`"""
        schema = manager.output_dir / JUPYTERLITE_SCHEMA
        file_dep = []

        if not schema.exists():  # pragma: no cover
            return

        file_dep += [schema]

        for lite_file in [
            *manager.output_dir.rglob(JUPYTERLITE_JSON),
            *manager.output_dir.rglob(JUPYTERLITE_IPYNB),
        ]:
            stem = lite_file.relative_to(manager.output_dir)
            selector = (
                None if lite_file.name == JUPYTERLITE_JSON else ["metadata", JUPYTERLITE_METADATA]
            )
            yield self.task(
                name=f"validate:{stem}",
                file_dep=[schema, lite_file],
                actions=[
                    (
                        self.validate_one_json_file,
                        [schema, lite_file, None, selector],
                    )
                ],
            )

    @property
    def lite_files(self):
        """all the source `jupyter-lite.*` files"""
        lite_dir = self.manager.lite_dir
        all_lite_files = [
            *lite_dir.rglob(JUPYTERLITE_JSON),
            *lite_dir.rglob(JUPYTERLITE_IPYNB),
        ]
        return [p for p in all_lite_files if not str(p).startswith(str(self.manager.output_dir))]
