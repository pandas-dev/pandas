# flake8: noqa: W503
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from shutil import rmtree
from typing import List

from tornado.concurrent import run_on_executor
from tornado.gen import convert_yielded

from .manager import lsp_message_listener
from .paths import file_uri_to_path, is_relative
from .types import LanguageServerManagerAPI

# TODO: make configurable
MAX_WORKERS = 4


def extract_or_none(obj, path):
    for crumb in path:
        try:
            obj = obj[crumb]
        except (KeyError, TypeError):
            return None
    return obj


class EditableFile:
    executor = ThreadPoolExecutor(max_workers=MAX_WORKERS)

    def __init__(self, path):
        # Python 3.5 relict:
        self.path = Path(path) if isinstance(path, str) else path

    async def read(self):
        self.lines = await convert_yielded(self.read_lines())

    async def write(self):
        return await convert_yielded(self.write_lines())

    @run_on_executor
    def read_lines(self):
        # empty string required by the assumptions of the gluing algorithm
        lines = [""]
        try:
            # TODO: what to do about bad encoding reads?
            lines = self.path.read_text(encoding="utf-8").splitlines()
        except FileNotFoundError:
            pass
        return lines

    @run_on_executor
    def write_lines(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text("\n".join(self.lines), encoding="utf-8")

    @staticmethod
    def trim(lines: list, character: int, side: int):
        needs_glue = False
        if lines:
            trimmed = lines[side][character:]
            if lines[side] != trimmed:
                needs_glue = True
            lines[side] = trimmed
        return needs_glue

    @staticmethod
    def join(left, right, glue: bool):
        if not glue:
            return []
        return [(left[-1] if left else "") + (right[0] if right else "")]

    def apply_change(self, text: str, start, end):
        before = self.lines[: start["line"]]
        after = self.lines[end["line"] :]

        needs_glue_left = self.trim(lines=before, character=start["character"], side=0)
        needs_glue_right = self.trim(lines=after, character=end["character"], side=-1)

        inner = text.split("\n")

        self.lines = (
            before[: -1 if needs_glue_left else None]
            + self.join(before, inner, needs_glue_left)
            + inner[1 if needs_glue_left else None : -1 if needs_glue_right else None]
            + self.join(inner, after, needs_glue_right)
            + after[1 if needs_glue_right else None :]
        ) or [""]

    @property
    def full_range(self):
        start = {"line": 0, "character": 0}
        end = {
            "line": len(self.lines),
            "character": len(self.lines[-1]) if self.lines else 0,
        }
        return {"start": start, "end": end}


WRITE_ONE = ["textDocument/didOpen", "textDocument/didChange", "textDocument/didSave"]


class ShadowFilesystemError(ValueError):
    """Error in the shadow file system."""


def setup_shadow_filesystem(virtual_documents_uri: str):
    if not virtual_documents_uri.startswith("file:/"):
        raise ShadowFilesystemError(  # pragma: no cover
            'Virtual documents URI has to start with "file:/", got '
            + virtual_documents_uri
        )

    initialized = False
    failures: List[Exception] = []

    shadow_filesystem = Path(file_uri_to_path(virtual_documents_uri))

    @lsp_message_listener("client")
    async def shadow_virtual_documents(scope, message, language_server, manager):
        """Intercept a message with document contents creating a shadow file for it.

        Only create the shadow file if the URI matches the virtual documents URI.
        Returns the path on filesystem where the content was stored.
        """
        nonlocal initialized

        # short-circut if language server does not require documents on disk
        server_spec = manager.language_servers[language_server]
        if not server_spec.get("requires_documents_on_disk", True):
            return

        if not message.get("method") in WRITE_ONE:
            return

        document = extract_or_none(message, ["params", "textDocument"])
        if document is None:
            raise ShadowFilesystemError(
                "Could not get textDocument from: {}".format(message)
            )

        uri = extract_or_none(document, ["uri"])
        if not uri:
            raise ShadowFilesystemError("Could not get URI from: {}".format(message))

        if not uri.startswith(virtual_documents_uri):
            return

        # initialization (/any file system operations) delayed until needed
        if not initialized:
            if len(failures) == 3:
                return
            try:
                # create if does no exist (so that removal does not raise)
                shadow_filesystem.mkdir(parents=True, exist_ok=True)
                # remove with contents
                rmtree(str(shadow_filesystem))
                # create again
                shadow_filesystem.mkdir(parents=True, exist_ok=True)
            except (OSError, PermissionError, FileNotFoundError) as e:
                failures.append(e)
                if len(failures) == 3:
                    manager.log.warn(
                        "[lsp] initialization of shadow filesystem failed three times"
                        " check if the path set by `LanguageServerManager.virtual_documents_dir`"
                        " or `JP_LSP_VIRTUAL_DIR` is correct; if this is happening with a server"
                        " for which you control (or wish to override) jupyter-lsp specification"
                        " you can try switching `requires_documents_on_disk` off. The errors were: %s",
                        failures,
                    )
                return
            initialized = True

        path = file_uri_to_path(uri)
        if not is_relative(shadow_filesystem, path):
            raise ShadowFilesystemError(
                f"Path {path} is not relative to shadow filesystem root"
            )

        editable_file = EditableFile(path)

        await editable_file.read()

        text = extract_or_none(document, ["text"])

        if text is not None:
            # didOpen and didSave may provide text within the document
            changes = [{"text": text}]
        else:
            # didChange is the only one which can also provide it in params (as contentChanges)
            if message["method"] != "textDocument/didChange":
                return
            if "contentChanges" not in message["params"]:
                raise ShadowFilesystemError(
                    "textDocument/didChange is missing contentChanges"
                )
            changes = message["params"]["contentChanges"]

        if len(changes) > 1:
            manager.log.warn(  # pragma: no cover
                "LSP warning: up to one change supported for textDocument/didChange"
            )

        for change in changes[:1]:
            change_range = change.get("range", editable_file.full_range)
            editable_file.apply_change(change["text"], **change_range)

        await editable_file.write()

        return path

    return shadow_virtual_documents
