from __future__ import annotations

import re
from typing import (
    Any,
    Dict,
    List,
    Match,
    MutableMapping,
    Optional,
    Set,
    Tuple,
)

from ._inline.emphasis import finalize_emphasis_tokens, is_entity_boundary
from ._inline.links import parse_link as parse_inline_link
from .core import InlineState, Parser
from .helpers import (
    HTML_ATTRIBUTES,
    HTML_TAGNAME,
    PUNCTUATION,
    unescape_char,
)
from .util import escape_url

_REGEX_META_CHARS = set(r"()[]{}?*+|.^$")
DEFAULT_MAX_EMPHASIS_DEPTH = 20
DEFAULT_MAX_IMAGE_DEPTH = 20

AUTO_EMAIL = (
    r"""<[a-zA-Z0-9.!#$%&'*+\/=?^_`{|}~-]+@[a-zA-Z0-9]"""
    r"(?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?"
    r"(?:\.[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?)*>"
)

INLINE_HTML = (
    r"<" + HTML_TAGNAME + HTML_ATTRIBUTES + r"\s*/?>|"  # open tag
    r"</" + HTML_TAGNAME + r"\s*>|"  # close tag
    r"<!--(?!>|->)(?:(?!--)[\s\S])+?(?<!-)-->|"  # comment
    r"<\?[\s\S]+?\?>|"  # script like <?php?>
    r"<![A-Z][\s\S]+?>|"  # doctype
    r"<!\[CDATA[\s\S]+?\]\]>"  # cdata
)


class InlineParser(Parser[InlineState]):
    sc_flag = 0
    state_cls = InlineState

    #: linebreak leaves two spaces at the end of line
    STD_LINEBREAK = r"(?:\\| {2,})\n\s*"

    #: every new line becomes <br>
    HARD_LINEBREAK = r" *\n\s*"

    # we only need to find the start pattern of an inline token
    SPECIFICATION = {
        # e.g. \`, \$
        "escape": r"(?:\\" + PUNCTUATION + ")+",
        # `code, ```code
        "codespan": r"`{1,}",
        # *w, **w, _w, __w
        "emphasis": r"\*{1,3}(?=[^\s*])|\b_{1,3}(?=[^\s_])",
        # [link], ![img]
        "link": r"!?\[",
        # <https://example.com>. regex copied from commonmark.js
        "auto_link": r"<[A-Za-z][A-Za-z0-9.+-]{1,31}:[^<>\x00-\x20]*>",
        "auto_email": AUTO_EMAIL,
        "inline_html": INLINE_HTML,
        "linebreak": STD_LINEBREAK,
        "softbreak": HARD_LINEBREAK,
        "prec_auto_link": r"<[A-Za-z][A-Za-z\d.+-]{1,31}:",
        "prec_inline_html": r"</?" + HTML_TAGNAME + r"|<!|<\?",
    }
    DEFAULT_RULES = (
        "escape",
        "codespan",
        "emphasis",
        "link",
        "auto_link",
        "auto_email",
        "inline_html",
        "linebreak",
    )

    def __init__(
        self,
        hard_wrap: bool = False,
        max_emphasis_depth: int = DEFAULT_MAX_EMPHASIS_DEPTH,
        max_image_depth: int = DEFAULT_MAX_IMAGE_DEPTH,
    ) -> None:
        super(InlineParser, self).__init__()

        self.hard_wrap = hard_wrap
        self.max_emphasis_depth = max_emphasis_depth
        self.max_image_depth = max_image_depth
        self._fast_trigger_chars: Optional[Set[str]] = None
        self._fast_trigger_re: Optional[re.Pattern[str]] = None
        self._fast_trigger_re_chars: Optional[Tuple[str, ...]] = None
        # lazy add linebreak
        if hard_wrap:
            self.specification["linebreak"] = self.HARD_LINEBREAK
        else:
            self.rules.append("softbreak")

        self._methods = {name: getattr(self, "parse_" + name) for name in self.rules}

    def register(
        self,
        name: str,
        pattern: Optional[str],
        func: Any,
        before: Optional[str] = None,
    ) -> None:
        super().register(name, pattern, func, before=before)
        self._fast_trigger_chars = None
        self._fast_trigger_re = None
        self._fast_trigger_re_chars = None

    def parse_escape(self, m: Match[str], state: InlineState) -> int:
        text = m.group(0)
        text = unescape_char(text)
        self.process_text(text, state, parse_emphasis=False)
        return m.end()

    def parse_link(self, m: Match[str], state: InlineState) -> Optional[int]:
        return parse_inline_link(self, m, state)

    def parse_auto_link(self, m: Match[str], state: InlineState) -> int:
        text = m.group(0)
        pos = m.end()
        if state.in_link:
            self.process_text(text, state)
            return pos

        text = text[1:-1]
        self._add_auto_link(text, text, state)
        return pos

    def parse_auto_email(self, m: Match[str], state: InlineState) -> int:
        text = m.group(0)
        pos = m.end()
        if state.in_link:
            self.process_text(text, state)
            return pos

        text = text[1:-1]
        url = "mailto:" + text
        self._add_auto_link(url, text, state)
        return pos

    def _add_auto_link(self, url: str, text: str, state: InlineState) -> None:
        state.append_token(
            {
                "type": "link",
                "children": [{"type": "text", "raw": text}],
                "attrs": {"url": escape_url(url)},
            }
        )

    def parse_emphasis(self, m: Match[str], state: InlineState) -> int:
        # Keep a delimiter run separate from the preceding text token.  The
        # emphasis finalizer needs to inspect these markers later, and merging
        # one-character markers with a growing text token makes inputs such as
        # ``*a*a*a`` repeatedly copy the whole accumulated string.
        marker = m.group(0)
        if len(marker) == 1:
            state.append_token({"type": "text", "raw": marker})
        else:
            self.process_text(marker, state)
        return m.end()

    def parse_codespan(self, m: Match[str], state: InlineState) -> int:
        marker = m.group(0)
        # require same marker with same length at end

        pattern = re.compile(r"(.*?[^`])" + marker + r"(?!`)", re.S)

        pos = m.end()
        m2 = pattern.match(state.src, pos)
        if m2:
            end_pos = m2.end()
            code = m2.group(1)
            # Line endings are treated like spaces
            code = code.replace("\n", " ")
            if len(code.strip()):
                if code.startswith(" ") and code.endswith(" "):
                    code = code[1:-1]
            state.append_token({"type": "codespan", "raw": code})
            return end_pos
        else:
            state.append_token({"type": "text", "raw": marker})
            return pos

    def parse_linebreak(self, m: Match[str], state: InlineState) -> int:
        state.append_token({"type": "linebreak"})
        return m.end()

    def parse_softbreak(self, m: Match[str], state: InlineState) -> int:
        state.append_token({"type": "softbreak"})
        return m.end()

    def parse_inline_html(self, m: Match[str], state: InlineState) -> int:
        end_pos = m.end()
        html = m.group(0)
        state.append_token({"type": "inline_html", "raw": html})
        if html.startswith(("<a ", "<a>", "<A ", "<A>")):
            state.in_link = True
        elif html.startswith(("</a ", "</a>", "</A ", "</A>")):
            state.in_link = False
        return end_pos

    def process_text(self, text: str, state: InlineState, parse_emphasis: bool = True) -> None:
        if (
            parse_emphasis
            and state.tokens
            and state.tokens[-1]["type"] == "text"
            and state.tokens[-1].get("_emphasis", True)
            and not is_entity_boundary(state.tokens[-1]["raw"], text)
        ):
            state.tokens[-1]["raw"] += text
        else:
            token: Dict[str, Any] = {"type": "text", "raw": text}
            if not parse_emphasis:
                token["_emphasis"] = False
            state.append_token(token)

    def parse(self, state: InlineState) -> List[Dict[str, Any]]:
        pos = 0
        sc = self.compile_sc()
        while pos < len(state.src):
            fast_end = self._find_fast_text_end(state.src, pos)
            if fast_end is None:
                m = sc.search(state.src, pos)
            else:
                if fast_end > pos:
                    self.process_text(state.src[pos:fast_end], state)
                    pos = fast_end
                if pos >= len(state.src):
                    break
                m = sc.match(state.src, pos)

            if not m:
                if fast_end is not None:
                    self.process_text(state.src[pos : pos + 1], state)
                    pos += 1
                    continue
                break

            end_pos = m.start()
            if end_pos > pos:
                hole = state.src[pos:end_pos]
                self.process_text(hole, state)

            new_pos = self.parse_method(m, state)
            if not new_pos:
                # move cursor 1 character forward
                pos = end_pos + 1
                hole = state.src[end_pos:pos]
                self.process_text(hole, state)
            else:
                pos = new_pos

        if pos == 0:
            # special case, just pure text
            self.process_text(state.src, state)
        elif pos < len(state.src):
            self.process_text(state.src[pos:], state)
        state.tokens = finalize_emphasis_tokens(
            state.tokens,
            "emphasis" in self.rules,
            self.max_emphasis_depth,
        )
        return state.tokens

    def _find_fast_text_end(self, src: str, pos: int) -> Optional[int]:
        chars = self._get_fast_trigger_chars()
        if chars is None:
            return None

        trigger_re = self._get_fast_trigger_re(chars)
        m = trigger_re.search(src, pos)
        if m is None:
            return len(src)

        if m.group(0) == "\n":
            return self._find_linebreak_start(src, pos, m.start())
        return m.start()

    def _get_fast_trigger_re(self, chars: Set[str]) -> re.Pattern[str]:
        key = tuple(sorted(chars))
        if self._fast_trigger_re is None or self._fast_trigger_re_chars != key:
            pattern = "[" + re.escape("".join(key)) + "]"
            self._fast_trigger_re = re.compile(pattern)
            self._fast_trigger_re_chars = key
        assert self._fast_trigger_re is not None
        return self._fast_trigger_re

    def _find_linebreak_start(self, src: str, min_pos: int, newline_pos: int) -> int:
        pos = newline_pos
        while pos > min_pos and src[pos - 1] == " ":
            pos -= 1
        if pos == newline_pos and pos > min_pos and src[pos - 1] == "\\":
            return pos - 1
        return pos

    def _get_fast_trigger_chars(self) -> Optional[Set[str]]:
        chars = self._fast_trigger_chars
        if chars is not None:
            return chars

        chars = set()
        for name in self.rules:
            pattern = self.specification.get(name)
            rule_chars = _get_rule_start_chars(name, pattern)
            if rule_chars is None:
                self._fast_trigger_chars = None
                return None
            chars.update(rule_chars)
        self._fast_trigger_chars = chars
        return chars

    def precedence_scan(
        self,
        m: Match[str],
        state: InlineState,
        end_pos: int,
        rules: Optional[List[str]] = None,
    ) -> Optional[int]:
        if rules is None:
            rules = ["codespan", "link", "prec_auto_link", "prec_inline_html"]

        mark_pos = m.end()
        sc = self.compile_sc(rules)
        m1 = sc.search(state.src, mark_pos, end_pos)
        if not m1:
            return None

        lastgroup = m1.lastgroup
        if not lastgroup:
            return None
        rule_name = lastgroup.replace("prec_", "")
        sc = self.compile_sc([rule_name])
        m2 = sc.match(state.src, m1.start())
        if not m2:
            return None

        func = self._methods[rule_name]
        new_state = state.copy()
        new_state.src = state.src
        m2_pos = func(m2, new_state)
        if not m2_pos or m2_pos < end_pos:
            return None

        raw_text = state.src[m.start() : m2.start()]
        state.append_token({"type": "text", "raw": raw_text})
        for token in new_state.tokens:
            state.append_token(token)
        return m2_pos

    def render(self, state: InlineState) -> List[Dict[str, Any]]:
        self.parse(state)
        return state.tokens

    def __call__(self, s: str, env: MutableMapping[str, Any]) -> List[Dict[str, Any]]:
        state = self.state_cls(env)
        state.src = s
        return self.render(state)


def _get_rule_start_chars(name: str, pattern: Optional[str]) -> Optional[Set[str]]:
    known = {
        "escape": {"\\"},
        "codespan": {"`"},
        "emphasis": {"*", "_"},
        "link": {"!", "["},
        "auto_link": {"<"},
        "auto_email": {"<"},
        "inline_html": {"<"},
        "linebreak": {"\n"},
        "softbreak": {"\n"},
        "prec_auto_link": {"<"},
        "prec_inline_html": {"<"},
        # built-in plugins
        "url_link": {"h"},
        "strikethrough": {"~"},
        "mark": {"="},
        "insert": {"^"},
        "superscript": {"^"},
        "subscript": {"~"},
        "footnote": {"["},
        "inline_math": {"$"},
        "ruby": {"["},
        "inline_spoiler": {">"},
    }
    if name in known:
        return known[name]
    if not pattern:
        return set()
    return _guess_pattern_start_chars(pattern)


def _guess_pattern_start_chars(pattern: str) -> Optional[Set[str]]:
    if not pattern:
        return set()

    if pattern.startswith("\\") and len(pattern) > 1:
        c = pattern[1]
        if c in _REGEX_META_CHARS or c in PUNCTUATION:
            return {c}
        return None

    c = pattern[0]
    if c in _REGEX_META_CHARS or c.isspace():
        return None
    return {c}
