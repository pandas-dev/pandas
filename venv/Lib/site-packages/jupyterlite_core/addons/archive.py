"""a JupyterLite addon for generating app archives which can be used as input"""

import contextlib
import gzip
import locale
import os
import tarfile
import tempfile
from hashlib import sha256
from pathlib import Path

from ..constants import C_LOCALE, MOD_FILE, NPM_SOURCE_DATE_EPOCH
from .base import BaseAddon


class ArchiveAddon(BaseAddon):
    """Adds contents from the ``lite_dir`` to the ``output_dir``, creates API output

    If ``--source-date-epoch`` (SDE) is set, a number of features
    will be enabled to improve reproducibility of the final artifact. In addition
    to timestamps newer than SDE being "clamped" to SDE, this will also adjust some
    permissions inside the tarball
    """

    __all__ = ["archive", "status"]

    def status(self, manager):
        tarball = manager.output_archive
        yield self.task(
            name="archive",
            actions=[(self.log_archive, [tarball])],
        )

    def archive(self, manager):
        """add all files created prior to ``pre_archive`` to an archive"""
        output_dir = manager.output_dir

        tarball = self.manager.output_archive

        file_dep = [p for p in output_dir.rglob("*") if not p.is_dir() and p not in [tarball]]

        yield self.task(
            name=f"archive:{tarball.name}",
            doc="generate a new app archive",
            file_dep=file_dep,
            actions=[
                (self.make_archive_stdlib, [tarball, output_dir, file_dep]),
                (self.log_archive, [tarball, "[lite] [archive] "]),
            ],
            targets=[tarball],
        )

    def filter_tarinfo(self, tarinfo: tarfile.TarInfo):
        """apply best-effort entropy fixes to give more reproducible archives"""
        tarinfo.uid = tarinfo.gid = 0
        tarinfo.uname = tarinfo.gname = "root"
        tarinfo.mode = MOD_FILE

        norm_path = str(Path(tarinfo.name).as_posix())

        if norm_path.startswith("package/files") or norm_path == "package/files":
            tarinfo.mtime = NPM_SOURCE_DATE_EPOCH
        elif self.manager.source_date_epoch is not None:
            tarinfo.mtime = self.manager.source_date_epoch

        return tarinfo

    @contextlib.contextmanager
    def setlocale(self, name):
        """Context manager for changing the current locale"""
        saved_locale = locale.setlocale(locale.LC_ALL)
        try:
            yield locale.setlocale(locale.LC_ALL, name)
        finally:
            locale.setlocale(locale.LC_ALL, saved_locale)

    def make_archive_stdlib(self, tarball, root, members):
        """actually build the archive.

        * this takes longer than any other hook
            * while this pure-python implementation needs to be maintained,
              a ``libarchive``-based build might be preferable for e.g. CI performance.
        * an npm-compatible ``.tgz`` is the only supported archive format, as this
          is compatible with the upstream ``webpack`` build and its native packaged format.
        """

        # if the command fails, but this still exists, it can cause problems
        if tarball.exists():
            tarball.unlink()

        # best-effort stable sorting
        with self.setlocale(C_LOCALE):
            members = sorted(root.rglob("*"), key=lambda p: locale.strxfrm(str(p)))

        len_members = str(len(members))
        rjust = len(len_members)

        with tempfile.TemporaryDirectory() as td:
            temp_ball = Path(td) / tarball.name
            with (
                os.fdopen(os.open(temp_ball, os.O_WRONLY | os.O_CREAT, MOD_FILE), "wb") as tar_gz,
                gzip.GzipFile(fileobj=tar_gz, mode="wb", mtime=0) as gz,
                tarfile.open(fileobj=gz, mode="w:") as tar,
            ):
                for i, path in enumerate(members):
                    if path.is_dir():
                        continue
                    if i == 0:
                        self.log.info(f"""[lite] [archive] files: {len_members}""")
                    if not (i % 100):
                        self.log.info(
                            """[lite] [archive] """
                            f"""... {str(i + 1).rjust(rjust)} """
                            f"""of {len_members}"""
                        )
                    tar.add(
                        path,
                        arcname=f"package/{path.relative_to(root)}",
                        filter=self.filter_tarinfo,
                        recursive=False,
                    )

            self.copy_one(temp_ball, tarball)

    def log_archive(self, tarball, prefix=""):
        """print some information about an archive"""
        sde = self.manager.source_date_epoch
        if sde is not None:
            self.log.info(f"{prefix}SOURCE_DATE_EPOCH: {sde}")
        if not tarball.exists():
            self.log.info(f"{prefix}No archive (yet): {tarball.name}")
        else:
            stat = tarball.stat()
            size = stat.st_size / (1024 * 1024)
            self.log.info(f"{prefix}filename:   {tarball.name}")
            shasum = sha256(tarball.read_bytes()).hexdigest()
            self.log.info(f"{prefix}size:       {size} Mb")
            # extra details, for the curious
            self.log.debug(f"{prefix}created:  {int(stat.st_mtime)}")
            self.log.debug(f"{prefix}modified: {int(stat.st_mtime)}")
            self.log.debug(f"{prefix}SHA256:   {shasum}")
