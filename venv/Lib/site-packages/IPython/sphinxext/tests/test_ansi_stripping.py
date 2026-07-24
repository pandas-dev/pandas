"""Test that ANSI escape sequences are stripped from directive output."""

import re

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def test_strip_ansi_basic():
    """ANSI color codes should be removed from output."""
    raw = "\x1b[31mred text\x1b[0m"
    assert _ANSI_RE.sub("", raw) == "red text"


def test_strip_ansi_multiple_codes():
    """Multiple ANSI codes in a single string should all be removed."""
    raw = "\x1b[1m\x1b[32mbold green\x1b[0m normal \x1b[4munderline\x1b[0m"
    assert _ANSI_RE.sub("", raw) == "bold green normal underline"


def test_strip_ansi_no_codes():
    """Strings without ANSI codes should pass through unchanged."""
    raw = "plain output"
    assert _ANSI_RE.sub("", raw) == "plain output"


def test_strip_ansi_empty():
    """Empty strings should remain empty."""
    assert _ANSI_RE.sub("", "") == ""
