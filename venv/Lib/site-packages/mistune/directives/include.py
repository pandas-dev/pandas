import os
from typing import TYPE_CHECKING, Any, Dict, List, Match, Union

from ..util import escape as escape_text
from ._base import BaseDirective, DirectivePlugin

if TYPE_CHECKING:
    from ..block_parser import BlockParser
    from ..core import BaseRenderer, BlockState
    from ..markdown import Markdown


class Include(DirectivePlugin):
    def parse(
        self, block: "BlockParser", m: Match[str], state: "BlockState"
    ) -> Union[Dict[str, Any], List[Dict[str, Any]]]:
        source_file = state.env.get("__file__")
        if not source_file:
            return {"type": "block_error", "raw": "Missing source file"}

        encoding = "utf-8"
        options = self.parse_options(m)
        if options:
            attrs = dict(options)
            if "encoding" in attrs:
                encoding = attrs["encoding"]
        else:
            attrs = {}

        relpath = self.parse_title(m)
        source_file = os.path.realpath(source_file)
        source_dir = os.path.dirname(source_file)
        dest = os.path.realpath(os.path.join(source_dir, relpath))

        if os.path.isabs(relpath) or os.path.commonpath([source_dir, dest]) != source_dir:
            return {
                "type": "block_error",
                "raw": "Could not include outside source dir: " + relpath,
            }

        if dest == source_file:
            return {
                "type": "block_error",
                "raw": "Could not include self: " + relpath,
            }

        include_stack = state.env.setdefault("__include_stack__", [])
        source_added = False
        if source_file not in include_stack:
            include_stack.append(source_file)
            source_added = True
        if dest in include_stack:
            if source_added:
                include_stack.pop()
            return {
                "type": "block_error",
                "raw": "Could not include circular reference: " + relpath,
            }

        if not os.path.isfile(dest):
            if source_added:
                include_stack.pop()
            return {
                "type": "block_error",
                "raw": "Could not find file: " + relpath,
            }

        include_stack.append(dest)
        try:
            with open(dest, "rb") as f:
                content = f.read().decode(encoding)

            ext = os.path.splitext(dest)[1]
            if ext in {".md", ".markdown", ".mkd"}:
                content = content.replace("\r\n", "\n").replace("\r", "\n")
                new_state = state.child_state(content)
                previous_file = new_state.env.get("__file__")
                new_state.env["__file__"] = dest
                try:
                    block.parse(new_state)
                finally:
                    if previous_file is None:
                        new_state.env.pop("__file__", None)
                    else:
                        new_state.env["__file__"] = previous_file
                return new_state.tokens

            elif ext in {".html", ".xhtml", ".htm"}:
                return {"type": "block_html", "raw": content}

            attrs["filepath"] = dest
            return {
                "type": "include",
                "raw": content,
                "attrs": attrs,
            }
        finally:
            include_stack.pop()
            if source_added:
                include_stack.pop()

    def __call__(self, directive: BaseDirective, md: "Markdown") -> None:
        directive.register("include", self.parse)
        if md.renderer and md.renderer.NAME == "html":
            md.renderer.register("include", render_html_include)


def render_html_include(renderer: "BaseRenderer", text: str, **attrs: Any) -> str:
    if getattr(renderer, "_escape", True):
        text = escape_text(text)
    return '<pre class="directive-include">\n' + text + "</pre>\n"
