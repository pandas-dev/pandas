import re
from textwrap import indent
from typing import Any, Dict, Iterable, cast

from ..core import BaseRenderer, BlockState
from ..util import strip_end
from ._list import render_list, render_list_item

fenced_re = re.compile(r"^[`~]+", re.M)
_backtick_run_re = re.compile(r"`+")

#: leading markers that would be parsed as a new block (list, heading, block
#: quote) if they appear unescaped at the start of a line.
_block_prefix_re = re.compile(r"^(\s*)(>|[-+*]|#{1,6}|\d{1,9}[.)])(\s|$)")


class MarkdownRenderer(BaseRenderer):
    """A renderer to re-format Markdown text."""

    NAME = "markdown"

    def __call__(self, tokens: Iterable[Dict[str, Any]], state: BlockState) -> str:
        out = self.render_tokens(tokens, state)
        # special handle for line breaks
        out += "\n\n".join(self.render_referrences(state)) + "\n"
        return strip_end(out)

    def render_referrences(self, state: BlockState) -> Iterable[str]:
        ref_links = state.env["ref_links"]
        for key in ref_links:
            attrs = ref_links[key]
            text = "[" + attrs["label"] + "]: " + attrs["url"]
            title = attrs.get("title")
            if title:
                text += ' "' + _escape_title(title) + '"'
            yield text

    def render_children(self, token: Dict[str, Any], state: BlockState) -> str:
        children = token["children"]
        return self.render_tokens(children, state)

    def text(self, token: Dict[str, Any], state: BlockState) -> str:
        raw = cast(str, token["raw"])
        # a text token that is made up entirely of "*"/"_" is a literal
        # emphasis delimiter -- either an escaped marker from the source
        # (``\*``) or an unmatched leftover -- so every character must stay
        # escaped, or it would re-parse as emphasis on the round-trip. Prose
        # punctuation such as ``2 * 3`` or ``snake_case`` arrives mixed with
        # other characters and is left untouched.
        if raw and all(c in "*_" for c in raw):
            return "".join("\\" + c for c in raw)
        # a backtick always opens a code span, so it must stay escaped to
        # survive a re-parse as literal text.
        return raw.replace("`", "\\`")

    def emphasis(self, token: Dict[str, Any], state: BlockState) -> str:
        return "*" + self.render_children(token, state) + "*"

    def strong(self, token: Dict[str, Any], state: BlockState) -> str:
        return "**" + self.render_children(token, state) + "**"

    def link(self, token: Dict[str, Any], state: BlockState) -> str:
        label = cast(str, token.get("label"))
        text = self.render_children(token, state)
        out = "[" + text + "]"
        if label:
            return out + "[" + label + "]"

        attrs = token["attrs"]
        url: str = attrs["url"]
        title = attrs.get("title")
        if text == url and not title:
            return "<" + text + ">"
        elif "mailto:" + text == url and not title:
            return "<" + text + ">"

        out += "("
        if "(" in url or ")" in url:
            out += "<" + url + ">"
        else:
            out += url
        if title:
            out += ' "' + _escape_title(title) + '"'
        return out + ")"

    def image(self, token: Dict[str, Any], state: BlockState) -> str:
        return "!" + self.link(token, state)

    def codespan(self, token: Dict[str, Any], state: BlockState) -> str:
        code = cast(str, token["raw"])
        # The delimiter must be a run of backticks longer than any backtick run
        # inside the content, otherwise the content would close the span early.
        longest = max((len(run) for run in _backtick_run_re.findall(code)), default=0)
        fence = "`" * (longest + 1)
        # A space on each side keeps the delimiter from merging with a leading or
        # trailing backtick; the parser strips this padding back off.
        if code.startswith("`") or code.endswith("`"):
            return fence + " " + code + " " + fence
        return fence + code + fence

    def linebreak(self, token: Dict[str, Any], state: BlockState) -> str:
        return "  \n"

    def softbreak(self, token: Dict[str, Any], state: BlockState) -> str:
        return "\n"

    def blank_line(self, token: Dict[str, Any], state: BlockState) -> str:
        return ""

    def inline_html(self, token: Dict[str, Any], state: BlockState) -> str:
        return cast(str, token["raw"])

    def paragraph(self, token: Dict[str, Any], state: BlockState) -> str:
        text = self.render_children(token, state)
        return _escape_block_prefix(text) + "\n\n"

    def heading(self, token: Dict[str, Any], state: BlockState) -> str:
        level = cast(int, token["attrs"]["level"])
        text = self.render_children(token, state)
        # An ATX heading ("# ...") occupies a single line, so a heading whose
        # rendered text spans several lines -- a multi-line setext heading --
        # must be re-emitted in setext form. As ATX the continuation lines
        # would fall out of the heading and become a paragraph on a re-parse.
        # Only levels 1 and 2 reach this branch, since an ATX heading never
        # contains a line break.
        if "\n" in text and level in (1, 2):
            underline = "=" if level == 1 else "-"
            return text + "\n" + underline * 3 + "\n\n"
        marker = "#" * level
        return marker + " " + text + "\n\n"

    def thematic_break(self, token: Dict[str, Any], state: BlockState) -> str:
        return "***\n\n"

    def block_text(self, token: Dict[str, Any], state: BlockState) -> str:
        return _escape_block_prefix(self.render_children(token, state)) + "\n"

    def block_code(self, token: Dict[str, Any], state: BlockState) -> str:
        attrs = token.get("attrs", {})
        info = cast(str, attrs.get("info", ""))
        code = cast(str, token["raw"])
        if code and code[-1] != "\n":
            code += "\n"

        marker = token.get("marker")
        if not marker:
            marker = _get_fenced_marker(code)
        marker2 = cast(str, marker)
        return marker2 + info + "\n" + code + marker2 + "\n\n"

    def block_quote(self, token: Dict[str, Any], state: BlockState) -> str:
        # strip the children's trailing blank lines first so the quote marker is
        # not added to a dangling empty line; stripping it back off afterwards
        # would also eat a ">" that ends the content (an autolink or HTML tag).
        text = self.render_children(token, state).rstrip("\n")
        text = indent(text, "> ", lambda _: True)
        return text + "\n\n"

    def block_html(self, token: Dict[str, Any], state: BlockState) -> str:
        return cast(str, token["raw"]) + "\n\n"

    def block_error(self, token: Dict[str, Any], state: BlockState) -> str:
        return ""

    def list(self, token: Dict[str, Any], state: BlockState) -> str:
        return render_list(self, token, state)

    def list_item(self, token: Dict[str, Any], state: BlockState) -> str:
        return render_list_item(self, token, state)

    def task_list_item(self, token: Dict[str, Any], state: BlockState) -> str:
        checked = token.get("attrs", {}).get("checked")
        marker = "[x] " if checked else "[ ] "
        return render_list_item(self, token, state, marker)

    def table(self, token: Dict[str, Any], state: BlockState) -> str:
        children = token.get("children", [])
        if not children:
            return "\n"

        head = children[0]
        body = children[1] if len(children) > 1 else None
        head_cells = head.get("children", [])
        align = [_table_cell_align(cell) for cell in head_cells]
        lines = [
            _render_table_row(self, head_cells, state),
            _render_table_delimiter(align),
        ]
        if body:
            for row in body.get("children", []):
                lines.append(_render_table_row(self, row.get("children", []), state))
        return "\n".join(lines) + "\n\n"

    def table_head(self, token: Dict[str, Any], state: BlockState) -> str:
        cells = token.get("children", [])
        return (
            _render_table_row(self, cells, state)
            + "\n"
            + _render_table_delimiter([_table_cell_align(c) for c in cells])
        )

    def table_body(self, token: Dict[str, Any], state: BlockState) -> str:
        return "\n".join(self.render_token(row, state).rstrip("\n") for row in token.get("children", []))

    def table_row(self, token: Dict[str, Any], state: BlockState) -> str:
        return _render_table_row(self, token.get("children", []), state) + "\n"

    def table_cell(self, token: Dict[str, Any], state: BlockState) -> str:
        return _render_table_cell(self, token, state)


def _escape_title(title: str) -> str:
    """Escape a link/image title for emission inside double quotes. The closing
    quote would otherwise end the title early on a re-parse; a backslash is
    escaped first so it can't combine with the following character."""
    return title.replace("\\", "\\\\").replace('"', '\\"')


def _escape_block_prefix(text: str) -> str:
    """Backslash-escape a leading block marker on each line so that literal
    text is not re-parsed as a list, heading or block quote."""
    return "\n".join(_escape_line_prefix(line) for line in text.split("\n"))


def _escape_line_prefix(line: str) -> str:
    m = _block_prefix_re.match(line)
    if not m:
        return line
    indent_, marker = m.group(1), m.group(2)
    return indent_ + marker[:-1] + "\\" + marker[-1] + line[m.end(2) :]


def _get_fenced_marker(code: str) -> str:
    found = fenced_re.findall(code)
    if not found:
        return "```"

    ticks = []  # `
    waves = []  # ~
    for s in found:
        if s[0] == "`":
            ticks.append(len(s))
        else:
            waves.append(len(s))

    if not ticks:
        return "```"

    if not waves:
        return "~~~"
    return "`" * (max(ticks) + 1)


def _render_table_row(renderer: MarkdownRenderer, cells: Iterable[Dict[str, Any]], state: BlockState) -> str:
    return "| " + " | ".join(_render_table_cell(renderer, cell, state) for cell in cells) + " |"


def _render_table_delimiter(aligns: Iterable[Any]) -> str:
    cells = []
    for align in aligns:
        if align == "left":
            cells.append(":---")
        elif align == "center":
            cells.append(":---:")
        elif align == "right":
            cells.append("---:")
        else:
            cells.append("---")
    return "| " + " | ".join(cells) + " |"


def _render_table_cell(renderer: MarkdownRenderer, token: Dict[str, Any], state: BlockState) -> str:
    if "children" in token:
        text = renderer.render_children(token, state)
    else:
        text = cast(str, token.get("raw", ""))
    return text.replace("\n", " ").replace("|", "\\|").strip()


def _table_cell_align(token: Dict[str, Any]) -> Any:
    return token.get("attrs", {}).get("align")
