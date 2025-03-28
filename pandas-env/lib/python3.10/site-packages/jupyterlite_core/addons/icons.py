"""a JupyterLite addon for customizing favicons

TODO:
  - this should do some best-effort image conversions.
  - may also impact `serviceworker.py`
"""

import pprint
from typing import TYPE_CHECKING

from ..optional import has_optional_dependency
from .base import BaseAddon

if TYPE_CHECKING:  # pragma: no cover
    from ..manager import LiteManager

from importlib.resources import files


class IconsAddon(BaseAddon):
    """Copy Jupyter Server favicons to /static"""

    __all__ = ["build", "status"]

    def status(self, manager: "LiteManager"):
        """yield some status information about the icons"""
        yield self.task(
            name="icons",
            actions=[
                lambda: self.log.debug(
                    "[lite] [icons] All favicons %s",
                    pprint.pformat([str(p) for p in self.favicon_files]),
                ),
                lambda: print(f"""    favicon files: {len(list(self.favicon_files))} files"""),
            ],
        )

    def build(self, manager: "LiteManager"):
        if not self.is_sys_prefix_ignored() and has_optional_dependency(
            "jupyter_server",
            "[lite] [icons] install `jupyter_server` to copy notebook favicons: {error}",
        ):
            # get the favicon files from the jupyter_server package
            src_favicons = files("jupyter_server") / "static" / "favicons"
            dest_favicons = self.favicon_dir

            yield self.task(
                name="copy",
                doc="copy the favicons",
                actions=[
                    (self.copy_one, [src_favicons, dest_favicons]),
                    (self.maybe_timestamp, [dest_favicons]),
                ],
            )

    @property
    def favicon_dir(self):
        return self.manager.output_dir / "static" / "favicons"

    @property
    def favicon_files(self):
        return sorted(self.favicon_dir.glob("*.ico"))
