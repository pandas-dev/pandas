from typing import TYPE_CHECKING, Match
from ..util import escape as escape_text

if TYPE_CHECKING:
    from ..block_parser import BlockParser
    from ..core import BaseRenderer, BlockState, InlineState
    from ..inline_parser import InlineParser
    from ..markdown import Markdown

__all__ = ["math", "math_in_quote", "math_in_list"]

BLOCK_MATH_PATTERN = (
    r"^ {0,3}\$\$(?P<math_text_single>[^\n]*?)\$\$[ \t]*(?:\n|$)|"
    r"^ {0,3}\$\$[ \t]*\n(?P<math_text_multi>[\s\S]*?)\n\$\$[ \t]*(?:\n|$)"
)
INLINE_MATH_PATTERN = (
    r"\$\$(?P<display_math_text>(?:[^$\\]|\\.)*?)\$\$|"
    r"\$(?P<backtick_math_marker>`+)(?P<backtick_math_text>[\s\S]*?)(?P=backtick_math_marker)\$|"
    r"\$(?!\$)(?!\s)(?P<math_text>(?:[^$\\\n]|\\.)+?)\$(?!\d)"
)


def parse_block_math(block: "BlockParser", m: Match[str], state: "BlockState") -> int:
    text = m.group("math_text_single")
    if text is None:
        text = m.group("math_text_multi")
    assert text is not None
    state.append_token({"type": "block_math", "raw": text})
    return m.end()


def parse_inline_math(inline: "InlineParser", m: Match[str], state: "InlineState") -> int:
    display_text = m.group("display_math_text")
    if display_text is not None:
        state.append_token({"type": "block_math", "raw": display_text})
        return m.end()

    text = m.group("backtick_math_text")
    if text is None:
        text = m.group("math_text")
    assert text is not None
    state.append_token({"type": "inline_math", "raw": text})
    return m.end()


def render_block_math(renderer: "BaseRenderer", text: str) -> str:
    return '<div class="math">$$\n' + escape_text(text) + "\n$$</div>\n"


def render_inline_math(renderer: "BaseRenderer", text: str) -> str:
    return r'<span class="math">\(' + escape_text(text) + r"\)</span>"


def math(md: "Markdown") -> None:
    """A mistune plugin to support math. The syntax is used
    by many markdown extensions:

    .. code-block:: text

        Block math is surrounded by $$:

        $$
        f(a)=f(b)
        $$

        Inline math is surrounded by `$`, such as $f(a)=f(b)$

    :param md: Markdown instance
    """
    md.block.register("block_math", BLOCK_MATH_PATTERN, parse_block_math, before="list")
    md.inline.register("inline_math", INLINE_MATH_PATTERN, parse_inline_math, before="codespan")
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("block_math", render_block_math)
        md.renderer.register("inline_math", render_inline_math)


def math_in_quote(md: "Markdown") -> None:
    """Enable block math plugin in block quote."""
    md.block.insert_rule(md.block.block_quote_rules, "block_math", before="list")


def math_in_list(md: "Markdown") -> None:
    """Enable block math plugin in list."""
    md.block.insert_rule(md.block.list_rules, "block_math", before="list")
