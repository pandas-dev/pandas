"""String filters.

Contains a collection of useful string manipulation filters for use in Jinja
templates.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import base64
import os
import re
import textwrap
import warnings
from urllib.parse import quote
from xml.etree.ElementTree import Element

import bleach

# defusedxml does safe(r) parsing of untrusted XML data
from defusedxml import ElementTree

from nbconvert.preprocessors.sanitize import _get_default_css_sanitizer

__all__ = [
    "add_anchor",
    "add_prompts",
    "ascii_only",
    "clean_html",
    "comment_lines",
    "get_lines",
    "html2text",
    "ipython2python",
    "path2url",
    "posix_path",
    "prevent_list_blocks",
    "strip_dollars",
    "strip_files_prefix",
    "strip_trailing_newline",
    "text_base64",
    "wrap_text",
]


def wrap_text(text, width=100):
    """
    Intelligently wrap text.
    Wrap text without breaking words if possible.

    Parameters
    ----------
    text : str
        Text to wrap.
    width : int, optional
        Number of characters to wrap to, default 100.
    """

    split_text = text.split("\n")
    wrp = map(lambda x: textwrap.wrap(x, width), split_text)  # noqa: C417
    wrpd = map("\n".join, wrp)
    return "\n".join(wrpd)


def html2text(element):
    """extract inner text from html

    Analog of jQuery's $(element).text()
    """
    if isinstance(element, (str,)):
        try:
            element = ElementTree.fromstring(element)
        except Exception:
            # failed to parse, just return it unmodified
            return element

    text = element.text or ""
    for child in element:
        text += html2text(child)
    text += element.tail or ""
    return text


def clean_html(element):
    """Clean an html element."""
    element = element.decode() if isinstance(element, bytes) else str(element)
    kwargs = {}
    css_sanitizer = _get_default_css_sanitizer()
    if css_sanitizer:
        kwargs["css_sanitizer"] = css_sanitizer
    return bleach.clean(
        element,
        tags=[*bleach.ALLOWED_TAGS, "div", "pre", "code", "span", "table", "tr", "td"],
        attributes={
            **bleach.ALLOWED_ATTRIBUTES,
            "*": ["class", "id"],
        },
        **kwargs,
    )


def _convert_header_id(header_contents):
    """Convert header contents to valid id value. Takes string as input, returns string.

    Note: this may be subject to change in the case of changes to how we wish to generate ids.

    For use on markdown headings.
    """
    # Valid IDs need to be non-empty and contain no space characters, but are otherwise arbitrary.
    # However, these IDs are also used in URL fragments, which are more restrictive, so we URL
    # encode any characters that are not valid in URL fragments.
    return quote(header_contents.replace(" ", "-"), safe="?/:@!$&'()*+,;=")


def add_anchor(html, anchor_link_text="¶"):
    """Add an id and an anchor-link to an html header

    For use on markdown headings
    """
    try:
        h = ElementTree.fromstring(html)
    except Exception:
        # failed to parse, just return it unmodified
        return html
    link = _convert_header_id(html2text(h))
    h.set("id", link)
    a = Element("a", {"class": "anchor-link", "href": "#" + link})
    try:
        # Test if the anchor link text is HTML (e.g. an image)
        a.append(ElementTree.fromstring(anchor_link_text))
    except Exception:
        # If we fail to parse, assume we've just got regular text
        a.text = anchor_link_text
    h.append(a)

    return ElementTree.tostring(h).decode(encoding="utf-8")


def add_prompts(code, first=">>> ", cont="... "):
    """Add prompts to code snippets"""
    new_code = []
    code_list = code.split("\n")
    new_code.append(first + code_list[0])
    for line in code_list[1:]:
        new_code.append(cont + line)
    return "\n".join(new_code)


def strip_dollars(text):
    """
    Remove all dollar symbols from text

    Parameters
    ----------
    text : str
        Text to remove dollars from
    """

    return text.strip("$")


files_url_pattern = re.compile(r'(src|href)\=([\'"]?)/?files/')
markdown_url_pattern = re.compile(r"(!?)\[(?P<caption>.*?)\]\(/?files/(?P<location>.*?)\)")


def strip_files_prefix(text):
    """
    Fix all fake URLs that start with ``files/``, stripping out the ``files/`` prefix.
    Applies to both urls (for html) and relative paths (for markdown paths).

    Parameters
    ----------
    text : str
        Text in which to replace 'src="files/real...' with 'src="real...'
    """
    cleaned_text = files_url_pattern.sub(r"\1=\2", text)
    cleaned_text = markdown_url_pattern.sub(r"\1[\2](\3)", cleaned_text)
    return cleaned_text  # noqa: RET504


def comment_lines(text, prefix="# "):
    """
    Build a Python comment line from input text.

    Parameters
    ----------
    text : str
        Text to comment out.
    prefix : str
        Character to append to the start of each line.
    """

    # Replace line breaks with line breaks and comment symbols.
    # Also add a comment symbol at the beginning to comment out
    # the first line.
    return prefix + ("\n" + prefix).join(text.split("\n"))


def get_lines(text, start=None, end=None):
    """
    Split the input text into separate lines and then return the
    lines that the caller is interested in.

    Parameters
    ----------
    text : str
        Text to parse lines from.
    start : int, optional
        First line to grab from.
    end : int, optional
        Last line to grab from.
    """

    # Split the input into lines.
    lines = text.split("\n")

    # Return the right lines.
    return "\n".join(lines[start:end])  # re-join


def ipython2python(code):
    """Transform IPython syntax to pure Python syntax

    Parameters
    ----------
    code : str
        IPython code, to be transformed to pure Python
    """
    try:
        from IPython.core.inputtransformer2 import TransformerManager  # noqa: PLC0415
    except ImportError:
        warnings.warn(
            "IPython is needed to transform IPython syntax to pure Python."
            " Install ipython if you need this functionality.",
            stacklevel=2,
        )
        return code
    else:
        isp = TransformerManager()
        return isp.transform_cell(code)


def posix_path(path):
    """Turn a path into posix-style path/to/etc

    Mainly for use in latex on Windows,
    where native Windows paths are not allowed.
    """
    if os.path.sep != "/":
        return path.replace(os.path.sep, "/")
    return path


def path2url(path):
    """Turn a file path into a URL"""
    parts = path.split(os.path.sep)
    return "/".join(quote(part) for part in parts)


def ascii_only(s):
    """ensure a string is ascii"""
    return s.encode("ascii", "replace").decode("ascii")


def prevent_list_blocks(s):
    """
    Prevent presence of enumerate or itemize blocks in latex headings cells
    """
    out = re.sub(r"(^\s*\d*)\.", r"\1\.", s)
    out = re.sub(r"(^\s*)\-", r"\1\-", out)
    out = re.sub(r"(^\s*)\+", r"\1\+", out)
    out = re.sub(r"(^\s*)\*", r"\1\*", out)
    return out  # noqa: RET504


def strip_trailing_newline(text):
    """
    Strips a newline from the end of text.
    """
    if text.endswith("\n"):
        text = text[:-1]
    return text


def text_base64(text):
    """
    Encode base64 text
    """
    return base64.b64encode(text.encode()).decode()
