# Abbreviation Extension for Python-Markdown
# ==========================================

# This extension adds abbreviation handling to Python-Markdown.

# See https://Python-Markdown.github.io/extensions/abbreviations
# for documentation.

# Original code Copyright 2007-2008 [Waylan Limberg](http://achinghead.com/)
# and [Seemant Kulleen](http://www.kulleen.org/)

# All changes Copyright 2008-2014 The Python Markdown Project

# License: [BSD](https://opensource.org/licenses/bsd-license.php)

"""
This extension adds abbreviation handling to Python-Markdown.

See the [documentation](https://Python-Markdown.github.io/extensions/abbreviations)
for details.
"""

from __future__ import annotations

from . import Extension
from ..blockprocessors import BlockProcessor
from ..inlinepatterns import InlineProcessor
from ..util import AtomicString
import re
import xml.etree.ElementTree as etree


class AbbrExtension(Extension):
    """ Abbreviation Extension for Python-Markdown. """

    def extendMarkdown(self, md):
        """ Insert `AbbrPreprocessor` before `ReferencePreprocessor`. """
        md.parser.blockprocessors.register(AbbrPreprocessor(md.parser), 'abbr', 16)


class AbbrPreprocessor(BlockProcessor):
    """ Abbreviation Preprocessor - parse text for abbr references. """

    RE = re.compile(r'^[*]\[(?P<abbr>[^\\]*?)\][ ]?:[ ]*\n?[ ]*(?P<title>.*)$', re.MULTILINE)

    def test(self, parent: etree.Element, block: str) -> bool:
        return True

    def run(self, parent: etree.Element, blocks: list[str]) -> bool:
        """
        Find and remove all Abbreviation references from the text.
        Each reference is set as a new `AbbrPattern` in the markdown instance.

        """
        block = blocks.pop(0)
        m = self.RE.search(block)
        if m:
            abbr = m.group('abbr').strip()
            title = m.group('title').strip()
            self.parser.md.inlinePatterns.register(
                AbbrInlineProcessor(self._generate_pattern(abbr), title), 'abbr-%s' % abbr, 2
            )
            if block[m.end():].strip():
                # Add any content after match back to blocks as separate block
                blocks.insert(0, block[m.end():].lstrip('\n'))
            if block[:m.start()].strip():
                # Add any content before match back to blocks as separate block
                blocks.insert(0, block[:m.start()].rstrip('\n'))
            return True
        # No match. Restore block.
        blocks.insert(0, block)
        return False

    def _generate_pattern(self, text: str) -> str:
        """ Given a string, returns a regex pattern to match that string. """
        return f"(?P<abbr>\\b{ re.escape(text) }\\b)"


class AbbrInlineProcessor(InlineProcessor):
    """ Abbreviation inline pattern. """

    def __init__(self, pattern: str, title: str):
        super().__init__(pattern)
        self.title = title

    def handleMatch(self, m: re.Match[str], data: str) -> tuple[etree.Element, int, int]:
        abbr = etree.Element('abbr')
        abbr.text = AtomicString(m.group('abbr'))
        abbr.set('title', self.title)
        return abbr, m.start(0), m.end(0)


def makeExtension(**kwargs):  # pragma: no cover
    return AbbrExtension(**kwargs)
