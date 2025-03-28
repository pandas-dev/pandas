"""a JupyterLite addon for Jupyter Server-compatible contents"""

import datetime
import json
import pprint
import re
from pathlib import Path

from ..constants import ALL_JSON, API_CONTENTS, JSON_FMT, UTF8
from ..optional import has_optional_dependency
from .base import BaseAddon


class ContentsAddon(BaseAddon):
    """Adds contents from the ``lite_dir`` to the ``output_dir``, creates API output"""

    __all__ = ["build", "post_build", "check", "status"]

    def status(self, manager):
        """yield some status information about the state of contents"""
        yield self.task(
            name="contents",
            actions=[
                lambda: self.log.debug(
                    "[lite] [contents] All Contents %s",
                    pprint.pformat([str(p[0]) for p in self.file_src_dest]),
                ),
                lambda: print(f"""    contents: {len(list(self.file_src_dest))} files"""),
            ],
        )

    def build(self, manager):
        """perform the main user build of pre-populating ``/files/``"""
        contents = sorted(self.file_src_dest)
        output_files_dir = self.output_files_dir
        all_dest_files = []

        for src_file, dest_file in contents:
            all_dest_files += [dest_file]
            rel = dest_file.relative_to(output_files_dir)
            yield self.task(
                name=f"copy:{rel}",
                doc=f"copy {src_file} to {rel}",
                file_dep=[src_file],
                targets=[dest_file],
                actions=[
                    (self.copy_one, [src_file, dest_file]),
                ],
            )

        if manager.source_date_epoch is not None:
            yield self.task(
                name="timestamp",
                file_dep=all_dest_files,
                actions=[
                    (self.maybe_timestamp, [self.output_files_dir]),
                ],
            )

    def post_build(self, manager):
        """create a Contents API index for each subdirectory in ``/files/``"""
        if not self.output_files_dir.exists():
            return

        output_file_dirs = [d for d in self.output_files_dir.rglob("*") if d.is_dir()] + [
            self.output_files_dir
        ]
        for output_file_dir in output_file_dirs:
            stem = output_file_dir.relative_to(self.output_files_dir)
            api_path = self.api_dir / stem / ALL_JSON

            yield self.task(
                name=f"contents:{stem}",
                doc=f"create a Jupyter Contents API response for {stem}",
                actions=[
                    (self.one_contents_path, [output_file_dir, api_path]),
                    (self.maybe_timestamp, [api_path]),
                ],
                file_dep=[p for p in output_file_dir.rglob("*") if not p.is_dir()],
                targets=[api_path],
            )

    def check(self, manager):
        """verify that all Contents API is valid (sorta)"""
        for all_json in self.api_dir.rglob(ALL_JSON):
            stem = all_json.relative_to(self.api_dir)
            yield self.task(
                name=f"validate:{stem}",
                doc=f"(eventually) validate {stem} with the Jupyter Contents API",
                file_dep=[all_json],
                actions=[(self.validate_one_json_file, [None, all_json])],
            )

    @property
    def api_dir(self):
        return self.manager.output_dir / API_CONTENTS

    @property
    def output_files_dir(self):
        return self.manager.output_dir / "files"

    @property
    def file_src_dest(self):
        """the pairs of contents that will be copied

        these are processed in `reverse` order, such that only the last path
        wins
        """
        yielded_dests = []
        for mgr_file in reversed(self.manager.contents):
            path = Path(mgr_file)
            for from_path in self.maybe_add_one_path(path):
                stem = from_path.relative_to(path) if path.is_dir() else path.name
                to_path = self.output_files_dir / stem
                resolved = str(to_path.resolve())
                if resolved in yielded_dests:  # pragma: no cover
                    self.log.debug("Already populated %s", resolved)
                    continue
                yielded_dests += [resolved]
                yield from_path, to_path

    def maybe_add_one_path(self, path, root=None):
        """add a file or folder's contents (if not ignored)"""

        if root is not None:
            rel_posix_path = f"/{path.relative_to(root).as_posix()}"

            for ignore in [
                *self.manager.ignore_contents,
                *self.manager.extra_ignore_contents,
            ]:
                if re.findall(ignore, rel_posix_path):
                    return

        if path.is_dir():
            for child in path.glob("*"):
                yield from self.maybe_add_one_path(child, root or path)
        else:
            yield path.resolve()

    def one_contents_path(self, output_file_dir, api_path):
        """A lazy reuse of a ``jupyter_server`` Contents API generator

        .. todo::

            Ideally we'd have a fallback, schema-verified generator, which we could
            later port to e.g. JS
        """
        if not has_optional_dependency(
            "jupyter_server",
            "[lite] [contents] install `jupyter_server` to index contents: {error}",
        ):
            return

        if not self.output_files_dir.exists():
            return

        self.maybe_timestamp(self.output_files_dir)

        from jupyter_server.services.contents.filemanager import FileContentsManager

        fm = FileContentsManager(root_dir=str(self.output_files_dir), parent=self)

        listing_path = str(output_file_dir.relative_to(self.output_files_dir))
        # normalize the root folder to avoid adding a `./` prefix to the
        # path field in the generated listing
        if listing_path == ".":
            listing_path = ""

        try:
            listing = fm.get(listing_path)
        except Exception as error:
            print(
                f"""Couldn't fetch {listing_path} as Jupyter contents.  {error}
                If this folder, or one of its parents, starts with a `.`, you can
                enable indexing hidden files with a `jupyter_lite_config.json` such as:

                    "ContentsManager": {{
                        "allow_hidden": true
                    }}

                Alternately, to skip it:

                    "LiteBuildConfig": {{
                        "extra_ignore_contents": [
                            "/\\.<the offendings path name>"
                        ]
                    }}
                """
            )
            return False

        if self.manager.source_date_epoch is not None:
            listing = self.patch_listing_timestamps(listing)

        api_path.parent.mkdir(parents=True, exist_ok=True)

        api_path.write_text(
            json.dumps(listing, **JSON_FMT, cls=DateTimeEncoder),
            **UTF8,
        )

        self.maybe_timestamp(api_path.parent)

    def patch_listing_timestamps(self, listing, sde=None):
        """clamp a contents listing's times to ``SOURCE_DATE_EPOCH``

        .. todo::

            pre-validated this structure with the ``jupyter_server`` API spec
        """
        sde = datetime.datetime.fromtimestamp(
            self.manager.source_date_epoch, tz=datetime.timezone.utc
        )

        if isinstance(listing, dict):
            for field in ["created", "last_modified"]:
                if field not in listing:  # pragma: no cover
                    continue
                value = listing[field]
                if isoformat(value) > isoformat(sde):
                    self.log.info(f"""[lite][contents][patch] {field} on {listing["name"]}""")
                    listing[field] = sde
            if listing["type"] == "directory":
                for child in listing.get("content") or []:
                    self.patch_listing_timestamps(child, sde)

        else:  # pragma: no cover
            self.log.error(f"[lite][contents] Don't know how to patch {listing}")
            return None

        return listing


class DateTimeEncoder(json.JSONEncoder):
    """A custom date-aware JSON encoder"""

    def default(self, o):
        if isinstance(o, datetime.datetime):
            return isoformat(o)

        return json.JSONEncoder.default(self, o)


def isoformat(dt):
    """a small helper to user ``Z`` for UTC ISO strings"""
    return dt.isoformat().replace("+00:00", "Z")
