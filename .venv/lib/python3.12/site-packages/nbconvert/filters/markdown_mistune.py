"""Markdown filters with mistune

Used from markdown.py
"""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import base64
import mimetypes
import os
from collections.abc import Iterable
from html import escape
from re import Match
from typing import TYPE_CHECKING, Any, ClassVar, Optional, Protocol

import bs4  # type: ignore[import-not-found]
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexer import Lexer
from pygments.lexers import get_lexer_by_name
from pygments.util import ClassNotFound

from nbconvert.filters.strings import add_anchor

if TYPE_CHECKING:
    try:
        from mistune.plugins import Plugin
    except ImportError:

        class Plugin(Protocol):  # type: ignore[no-redef]
            """Mistune plugin interface."""

            def __call__(self, markdown: "Markdown") -> None:
                """Apply the plugin on the markdown document."""
                ...


try:  # for Mistune >= 3.0
    from mistune import (  # type:ignore[attr-defined]
        BlockParser,
        BlockState,
        HTMLRenderer,
        InlineParser,
        InlineState,
        Markdown,
        import_plugin,
    )

    MISTUNE_V3 = True
    MISTUNE_V3_ATX = "atx_heading" in BlockParser.SPECIFICATION

except ImportError:  # for Mistune >= 2.0
    import re

    from mistune import (  # type: ignore[attr-defined]
        PLUGINS,
        BlockParser,
        HTMLRenderer,
        InlineParser,
        Markdown,
    )

    MISTUNE_V3 = False
    MISTUNE_V3_ATX = False

    def import_plugin(name: str) -> "Plugin":  # type: ignore[misc]
        """Simple implementation of Mistune V3's import_plugin for V2."""
        return PLUGINS[name]  # type: ignore[no-any-return]


class InvalidNotebook(Exception):
    """An invalid notebook model."""


def _dotall(pattern: str) -> str:
    """Makes the '.' special character match any character inside the pattern, including a newline.

    This is implemented with the inline flag `(?s:...)` and is equivalent to using `re.DOTALL`.
    It is useful for LaTeX environments, where line breaks may be present.
    """
    return f"(?s:{pattern})"


if MISTUNE_V3:  # Parsers for Mistune >= 3.0.0

    class MathBlockParser(BlockParser):
        """This acts as a pass-through to the MathInlineParser. It is needed in
        order to avoid other block level rules splitting math sections apart.

        It works by matching each multiline math environment as a single paragraph,
        so that other rules don't think each section is its own paragraph. Inline
        is ignored here.
        """

        ATX_HEADING_WITHOUT_LEADING_SPACES = (
            r"^ {0,3}(?P<atx_1>#{1,6})(?!#+)(?P<atx_2>[ \t]*(.*?)?)$"
            if MISTUNE_V3_ATX
            else r"^ {0,3}(?P<axt_1>#{1,6})(?!#+)(?P<axt_2>[ \t]*(.*?)?)$"
        )

        MULTILINE_MATH = _dotall(
            # Display math mode, old TeX delimiter: $$ \sqrt{2} $$
            r"(?<!\\)[$]{2}.*?(?<!\\)[$]{2}"
            "|"
            # Display math mode, new LaTeX delimiter: \[ \sqrt{2} \]
            r"\\\\\[.*?\\\\\]"
            "|"
            # LaTeX environment: \begin{equation} \sqrt{2} \end{equation}
            r"\\begin\{(?P<math_env_name>[a-z]*\*?)\}.*?\\end\{(?P=math_env_name)\}"
        )

        SPECIFICATION = {
            **BlockParser.SPECIFICATION,
            (
                "atx_heading" if MISTUNE_V3_ATX else "axt_heading"
            ): ATX_HEADING_WITHOUT_LEADING_SPACES,
            "multiline_math": MULTILINE_MATH,
        }

        # Multiline math must be searched before other rules
        DEFAULT_RULES: ClassVar[Iterable[str]] = ("multiline_math", *BlockParser.DEFAULT_RULES)  # type: ignore[assignment]

        def parse_multiline_math(self, m: Match[str], state: BlockState) -> int:
            """Send mutiline math as a single paragraph to MathInlineParser."""
            matched_text = m[0]
            state.add_paragraph(matched_text)
            return m.end()

    class MathInlineParser(InlineParser):
        r"""This interprets the content of LaTeX style math objects.

        In particular this grabs ``$$...$$``, ``\\[...\\]``, ``\\(...\\)``, ``$...$``,
        and ``\begin{foo}...\end{foo}`` styles for declaring mathematics. It strips
        delimiters from all these varieties, and extracts the type of environment
        in the last case (``foo`` in this example).
        """

        # Display math mode, using older TeX delimiter: $$ \pi $$
        BLOCK_MATH_TEX = _dotall(r"(?<!\\)\$\$(?P<math_block_tex>.*?)(?<!\\)\$\$")
        # Display math mode, using newer LaTeX delimiter: \[ \pi \]
        BLOCK_MATH_LATEX = _dotall(r"(?<!\\)\\\\\[(?P<math_block_latex>.*?)(?<!\\)\\\\\]")
        # Inline math mode, using older TeX delimiter: $ \pi $  (cannot be empty!)
        INLINE_MATH_TEX = _dotall(r"(?<![$\\])\$(?P<math_inline_tex>.+?)(?<![$\\])\$")
        # Inline math mode, using newer LaTeX delimiter: \( \pi \)
        INLINE_MATH_LATEX = _dotall(r"(?<!\\)\\\\\((?P<math_inline_latex>.*?)(?<!\\)\\\\\)")
        # LaTeX math environment: \begin{equation} \pi \end{equation}
        LATEX_ENVIRONMENT = _dotall(
            r"\\begin\{(?P<math_env_name>[a-z]*\*?)\}"
            r"(?P<math_env_body>.*?)"
            r"\\end\{(?P=math_env_name)\}"
        )

        SPECIFICATION = {
            **InlineParser.SPECIFICATION,
            "block_math_tex": BLOCK_MATH_TEX,
            "block_math_latex": BLOCK_MATH_LATEX,
            "inline_math_tex": INLINE_MATH_TEX,
            "inline_math_latex": INLINE_MATH_LATEX,
            "latex_environment": LATEX_ENVIRONMENT,
        }

        # Block math must be matched first, and all math must come before text
        DEFAULT_RULES: ClassVar[Iterable[str]] = (
            "block_math_tex",
            "block_math_latex",
            "inline_math_tex",
            "inline_math_latex",
            "latex_environment",
            *InlineParser.DEFAULT_RULES,
        )  # type: ignore[assignment]

        def parse_block_math_tex(self, m: Match[str], state: InlineState) -> int:
            """Parse older TeX-style display math."""
            body = m.group("math_block_tex")
            state.append_token({"type": "block_math", "raw": body})
            return m.end()

        def parse_block_math_latex(self, m: Match[str], state: InlineState) -> int:
            """Parse newer LaTeX-style display math."""
            body = m.group("math_block_latex")
            state.append_token({"type": "block_math", "raw": body})
            return m.end()

        def parse_inline_math_tex(self, m: Match[str], state: InlineState) -> int:
            """Parse older TeX-style inline math."""
            body = m.group("math_inline_tex")
            state.append_token({"type": "inline_math", "raw": body})
            return m.end()

        def parse_inline_math_latex(self, m: Match[str], state: InlineState) -> int:
            """Parse newer LaTeX-style inline math."""
            body = m.group("math_inline_latex")
            state.append_token({"type": "inline_math", "raw": body})
            return m.end()

        def parse_latex_environment(self, m: Match[str], state: InlineState) -> int:
            """Parse a latex environment."""
            attrs = {"name": m.group("math_env_name"), "body": m.group("math_env_body")}
            state.append_token({"type": "latex_environment", "attrs": attrs})
            return m.end()

else:  # Parsers for Mistune >= 2.0.0 < 3.0.0

    class MathBlockParser(BlockParser):  # type: ignore[no-redef]
        """This acts as a pass-through to the MathInlineParser. It is needed in
        order to avoid other block level rules splitting math sections apart.
        """

        MULTILINE_MATH = re.compile(
            # Display math mode, old TeX delimiter: $$ \sqrt{2} $$
            r"(?<!\\)[$]{2}.*?(?<!\\)[$]{2}|"
            # Display math mode, new LaTeX delimiter: \[ \sqrt{2} \]
            r"\\\\\[.*?\\\\\]|"
            # LaTeX environment: \begin{equation} \sqrt{2} \end{equation}
            r"\\begin\{([a-z]*\*?)\}.*?\\end\{\1\}",
            re.DOTALL,
        )

        # Regex for header that doesn't require space after '#'
        AXT_HEADING = re.compile(r" {0,3}(#{1,6})(?!#+)(?: *\n+|([^\n]*?)(?:\n+|\s+?#+\s*\n+))")

        # Multiline math must be searched before other rules
        RULE_NAMES = ("multiline_math", *BlockParser.RULE_NAMES)  # type: ignore[attr-defined]

        def parse_multiline_math(self, m: Match[str], state: Any) -> dict[str, str]:
            """Pass token through mutiline math."""
            return {"type": "multiline_math", "text": m.group(0)}

    class MathInlineParser(InlineParser):  # type: ignore[no-redef]
        r"""This interprets the content of LaTeX style math objects.

        In particular this grabs ``$$...$$``, ``\\[...\\]``, ``\\(...\\)``, ``$...$``,
        and ``\begin{foo}...\end{foo}`` styles for declaring mathematics. It strips
        delimiters from all these varieties, and extracts the type of environment
        in the last case (``foo`` in this example).
        """

        # Display math mode, using older TeX delimiter: $$ \pi $$
        BLOCK_MATH_TEX = _dotall(r"(?<!\\)\$\$(.*?)(?<!\\)\$\$")
        # Display math mode, using newer LaTeX delimiter: \[ \pi \]
        BLOCK_MATH_LATEX = _dotall(r"(?<!\\)\\\\\[(.*?)(?<!\\)\\\\\]")
        # Inline math mode, using older TeX delimiter: $ \pi $  (cannot be empty!)
        INLINE_MATH_TEX = _dotall(r"(?<![$\\])\$(.+?)(?<![$\\])\$")
        # Inline math mode, using newer LaTeX delimiter: \( \pi \)
        INLINE_MATH_LATEX = _dotall(r"(?<!\\)\\\\\((.*?)(?<!\\)\\\\\)")
        # LaTeX math environment: \begin{equation} \pi \end{equation}
        LATEX_ENVIRONMENT = _dotall(r"\\begin\{([a-z]*\*?)\}(.*?)\\end\{\1\}")

        RULE_NAMES = (
            "block_math_tex",
            "block_math_latex",
            "inline_math_tex",
            "inline_math_latex",
            "latex_environment",
            *InlineParser.RULE_NAMES,  # type: ignore[attr-defined]
        )

        def parse_block_math_tex(self, m: Match[str], state: Any) -> tuple[str, str]:
            """Parse block text math."""
            # sometimes the Scanner keeps the final '$$', so we use the
            # full matched string and remove the math markers
            text = m.group(0)[2:-2]
            return "block_math", text

        def parse_block_math_latex(self, m: Match[str], state: Any) -> tuple[str, str]:
            """Parse block latex math ."""
            text = m.group(1)
            return "block_math", text

        def parse_inline_math_tex(self, m: Match[str], state: Any) -> tuple[str, str]:
            """Parse inline tex math."""
            text = m.group(1)
            return "inline_math", text

        def parse_inline_math_latex(self, m: Match[str], state: Any) -> tuple[str, str]:
            """Parse inline latex math."""
            text = m.group(1)
            return "inline_math", text

        def parse_latex_environment(self, m: Match[str], state: Any) -> tuple[str, str, str]:
            """Parse a latex environment."""
            name, text = m.group(1), m.group(2)
            return "latex_environment", name, text


class IPythonRenderer(HTMLRenderer):
    """An ipython html renderer."""

    def __init__(
        self,
        escape: bool = True,
        allow_harmful_protocols: bool = True,
        embed_images: bool = False,
        exclude_anchor_links: bool = False,
        anchor_link_text: str = "¶",
        path: str = "",
        attachments: Optional[dict[str, dict[str, str]]] = None,
        **lexer_options,
    ):
        """Initialize the renderer."""
        super().__init__(escape, allow_harmful_protocols)
        self.embed_images = embed_images
        self.exclude_anchor_links = exclude_anchor_links
        self.anchor_link_text = anchor_link_text
        self.path = path
        self.lexer_options = lexer_options
        if attachments is not None:
            self.attachments = attachments
        else:
            self.attachments = {}

    def block_code(self, code: str, info: Optional[str] = None) -> str:
        """Handle block code."""
        lang: Optional[str] = ""
        lexer: Optional[Lexer] = None

        if info:
            if info.startswith("mermaid"):
                return self.block_mermaidjs(code)

            try:
                if info.strip().split(None, 1):
                    lang = info.strip().split(maxsplit=1)[0]
                    lexer = get_lexer_by_name(lang, **self.lexer_options)
            except ClassNotFound:
                code = f"{lang}\n{code}"
                lang = None

        if not lang:
            return super().block_code(code, info=info)

        formatter = HtmlFormatter()
        return highlight(code, lexer, formatter)

    def block_mermaidjs(self, code: str) -> str:
        """Handle mermaid syntax."""
        return (
            """<div class="jp-Mermaid"><pre class="mermaid">\n"""
            f"""{code.strip()}"""
            """\n</pre></div>"""
        )

    def block_html(self, html: str) -> str:
        """Handle block html."""
        if self.embed_images:
            html = self._html_embed_images(html)

        return super().block_html(html)

    def inline_html(self, html: str) -> str:
        """Handle inline html."""
        if self.embed_images:
            html = self._html_embed_images(html)

        return super().inline_html(html)

    def heading(self, text: str, level: int, **attrs: dict[str, Any]) -> str:
        """Handle a heading."""
        html = super().heading(text, level, **attrs)
        if self.exclude_anchor_links:
            return html
        return str(add_anchor(html, anchor_link_text=self.anchor_link_text))

    def escape_html(self, text: str) -> str:
        """Escape html content."""
        return escape(text, quote=False)

    def block_math(self, body: str) -> str:
        """Handle block math."""
        return f"$${self.escape_html(body)}$$"

    def multiline_math(self, text: str) -> str:
        """Handle mulitline math for older mistune versions."""
        return text

    def latex_environment(self, name: str, body: str) -> str:
        """Handle a latex environment."""
        name, body = self.escape_html(name), self.escape_html(body)
        return f"\\begin{{{name}}}{body}\\end{{{name}}}"

    def inline_math(self, body: str) -> str:
        """Handle inline math."""
        return f"${self.escape_html(body)}$"

    def image(self, text: str, url: str, title: Optional[str] = None) -> str:
        """Rendering a image with title and text.

        :param text: alt text of the image.
        :param url: source link of the image.
        :param title: title text of the image.

        :note: The parameters `text` and `url` are swapped in older versions
            of mistune.
        """
        if MISTUNE_V3:
            url = self._embed_image_or_attachment(url)
        else:  # for mistune v2, the first argument is the URL
            text = self._embed_image_or_attachment(text)

        return super().image(text, url, title)

    def _embed_image_or_attachment(self, src: str) -> str:
        """Embed an image or attachment, depending on the configuration.
        If neither is possible, returns the original URL.
        """

        attachment_prefix = "attachment:"
        if src.startswith(attachment_prefix):
            name = src[len(attachment_prefix) :]

            if name not in self.attachments:
                msg = f"missing attachment: {name}"
                raise InvalidNotebook(msg)

            attachment = self.attachments[name]
            # we choose vector over raster, and lossless over lossy
            preferred_mime_types = ("image/svg+xml", "image/png", "image/jpeg")
            for mime_type in preferred_mime_types:
                if mime_type in attachment:
                    return f"data:{mime_type};base64,{attachment[mime_type]}"
            # otherwise we choose the first mimetype we can find
            default_mime_type = next(iter(attachment.keys()))
            return f"data:{default_mime_type};base64,{attachment[default_mime_type]}"

        if self.embed_images:
            base64_url = self._src_to_base64(src)
            if base64_url is not None:
                return base64_url

        return src

    def _src_to_base64(self, src: str) -> Optional[str]:
        """Turn the source file into a base64 url.

        :param src: source link of the file.
        :return: the base64 url or None if the file was not found.
        """
        src_path = os.path.join(self.path, src)

        if not os.path.exists(src_path):
            return None

        with open(src_path, "rb") as fobj:
            mime_type, _ = mimetypes.guess_type(src_path)

            base64_data = base64.b64encode(fobj.read())
            base64_str = base64_data.replace(b"\n", b"").decode("ascii")

            return f"data:{mime_type};base64,{base64_str}"

    def _html_embed_images(self, html: str) -> str:
        parsed_html = bs4.BeautifulSoup(html, features="html.parser")
        imgs: bs4.ResultSet[bs4.Tag] = parsed_html.find_all("img")

        # Replace img tags's sources by base64 dataurls
        for img in imgs:
            src = img.attrs.get("src")
            if src is None:
                continue

            base64_url = self._src_to_base64(img.attrs["src"])
            if base64_url is not None:
                img.attrs["src"] = base64_url

        return str(parsed_html)


class MarkdownWithMath(Markdown):
    """Markdown text with math enabled."""

    DEFAULT_PLUGINS = (
        # "abbr",  (see https://github.com/jupyter/nbconvert/pull/1853)
        # "footnotes",
        "strikethrough",
        "table",
        "url",
        "task_lists",
        "def_list",
    )

    def __init__(
        self,
        renderer: HTMLRenderer,
        block: Optional[BlockParser] = None,
        inline: Optional[InlineParser] = None,
        plugins: Optional[Iterable["Plugin"]] = None,
    ):
        """Initialize the parser."""
        if block is None:
            block = MathBlockParser()
        if inline is None:
            if MISTUNE_V3:
                inline = MathInlineParser(hard_wrap=False)
            else:
                inline = MathInlineParser(renderer, hard_wrap=False)  # type: ignore[arg-type,misc]
        if plugins is None:
            plugins = (import_plugin(p) for p in self.DEFAULT_PLUGINS)

        super().__init__(renderer, block, inline, plugins)

    def render(self, source: str) -> str:
        """Render the HTML output for a Markdown source."""
        return str(super().__call__(source))


def markdown2html_mistune(source: str) -> str:
    """Convert a markdown string to HTML using mistune"""
    return MarkdownWithMath(renderer=IPythonRenderer(escape=False)).render(source)
