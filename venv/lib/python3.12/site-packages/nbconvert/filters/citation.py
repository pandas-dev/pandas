"""Citation handling for LaTeX output."""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
from html.parser import HTMLParser

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

__all__ = ["citation2latex"]


def citation2latex(s):
    """Parse citations in Markdown cells.

    This looks for HTML tags having a data attribute names ``data-cite``
    and replaces it by the call to LaTeX cite command. The transformation
    looks like this::

        <cite data-cite="granger">(Granger, 2013)</cite>

    Becomes ::

        \\cite{granger}

    Any HTML tag can be used, which allows the citations to be formatted
    in HTML in any manner.
    """
    parser = CitationParser()
    parser.feed(s)
    parser.close()
    outtext = ""
    startpos = 0
    for citation in parser.citelist:
        outtext += s[startpos : citation[1]]
        outtext += "\\cite{%s}" % citation[0]
        startpos = citation[2] if len(citation) == 3 else -1
    outtext += s[startpos:] if startpos != -1 else ""
    return outtext


# -----------------------------------------------------------------------------
# Classes
# -----------------------------------------------------------------------------
class CitationParser(HTMLParser):
    """Citation Parser

    Replaces html tags with data-cite attribute with respective latex \\cite.

    Inherites from HTMLParser, overrides:
     - handle_starttag
     - handle_endtag
    """

    # number of open tags
    opentags = None
    # list of found citations
    citelist = None  # type:ignore[var-annotated]
    # active citation tag
    citetag = None

    def __init__(self):
        """Initialize the parser."""
        self.citelist = []
        self.opentags = 0
        HTMLParser.__init__(self)

    def get_offset(self):
        """Get the offset position."""
        # Compute startposition in source
        lin, offset = self.getpos()
        pos = 0
        for _ in range(lin - 1):
            pos = self.data.find("\n", pos) + 1
        return pos + offset

    def handle_starttag(self, tag, attrs):
        """Handle a start tag."""
        # for each tag check if attributes are present and if no citation is active
        if self.opentags == 0 and len(attrs) > 0:
            for atr, data in attrs:
                if atr.lower() == "data-cite":
                    self.citetag = tag
                    self.opentags = 1
                    self.citelist.append([data, self.get_offset()])
                    return

        if tag == self.citetag:
            # found an open citation tag but not the starting one
            self.opentags += 1  # type:ignore[operator]

    def handle_endtag(self, tag):
        """Handle an end tag."""
        if tag == self.citetag:
            # found citation tag check if starting one
            if self.opentags == 1:
                pos = self.get_offset()
                self.citelist[-1].append(pos + len(tag) + 3)
            self.opentags -= 1  # type:ignore[operator]

    def feed(self, data):
        """Handle a feed."""
        self.data = data
        HTMLParser.feed(self, data)
