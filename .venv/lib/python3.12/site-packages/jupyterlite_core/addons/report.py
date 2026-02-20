"""a JupyterLite addon for generating hashes"""

from typing import TYPE_CHECKING

from ..constants import SHA256SUMS
from .base import BaseAddon

if TYPE_CHECKING:  # pragma: no cover
    from ..manager import LiteManager


class ReportAddon(BaseAddon):
    """update static listings of the site contents in various formats

    having these in various formats down the line can be handy for various publishing
    tasks
    """

    __all__ = ["pre_archive"]

    def pre_archive(self, manager: "LiteManager"):
        """generate a hash file of all files in the distribution.

        As this is relatively expensive for hundreds of files, this is performed
        as late as possible, while still providing some useful publishing / QA
        features.

        TODO: develop some contract with the frontend in relation to this file,
              or a derivative, as it has precisely the right information for certain
              cache tasks.
        """
        sha256sums = self.sha256sums

        all_output_files = self.all_output_files

        yield self.task(
            name=SHA256SUMS,
            doc="hash all of the files",
            actions=[
                (self.hash_all, [sha256sums, manager.output_dir, all_output_files]),
            ],
            file_dep=all_output_files,
            targets=[sha256sums],
        )

    @property
    def sha256sums(self):
        """The location of the hashfile."""
        return self.manager.output_dir / SHA256SUMS

    @property
    def all_output_files(self):
        return [
            p
            for p in sorted(self.manager.output_dir.rglob("*"))
            if not p.is_dir() and p not in [self.sha256sums, self.manager.output_archive]
        ]
