# Copyright 2010-2025 Kurt McKee <contactme@kurtmckee.org>
# Copyright 2002-2008 Mark Pilgrim
# All rights reserved.
#
# This file is a part of feedparser.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import html.entities
import re

from .sgml import *

_cp1252 = {
    128: '\u20ac',  # euro sign
    130: '\u201a',  # single low-9 quotation mark
    131: '\u0192',  # latin small letter f with hook
    132: '\u201e',  # double low-9 quotation mark
    133: '\u2026',  # horizontal ellipsis
    134: '\u2020',  # dagger
    135: '\u2021',  # double dagger
    136: '\u02c6',  # modifier letter circumflex accent
    137: '\u2030',  # per mille sign
    138: '\u0160',  # latin capital letter s with caron
    139: '\u2039',  # single left-pointing angle quotation mark
    140: '\u0152',  # latin capital ligature oe
    142: '\u017d',  # latin capital letter z with caron
    145: '\u2018',  # left single quotation mark
    146: '\u2019',  # right single quotation mark
    147: '\u201c',  # left double quotation mark
    148: '\u201d',  # right double quotation mark
    149: '\u2022',  # bullet
    150: '\u2013',  # en dash
    151: '\u2014',  # em dash
    152: '\u02dc',  # small tilde
    153: '\u2122',  # trade mark sign
    154: '\u0161',  # latin small letter s with caron
    155: '\u203a',  # single right-pointing angle quotation mark
    156: '\u0153',  # latin small ligature oe
    158: '\u017e',  # latin small letter z with caron
    159: '\u0178',  # latin capital letter y with diaeresis
}


class _BaseHTMLProcessor(sgmllib.SGMLParser, object):
    special = re.compile("""[<>'"]""")
    bare_ampersand = re.compile(r"&(?!#\d+;|#x[0-9a-fA-F]+;|\w+;)")
    elements_no_end_tag = {
        'area',
        'base',
        'basefont',
        'br',
        'col',
        'command',
        'embed',
        'frame',
        'hr',
        'img',
        'input',
        'isindex',
        'keygen',
        'link',
        'meta',
        'param',
        'source',
        'track',
        'wbr',
    }

    def __init__(self, encoding=None, _type='application/xhtml+xml'):
        if encoding:
            self.encoding = encoding
        self._type = _type
        self.pieces = []
        super(_BaseHTMLProcessor, self).__init__()

    def reset(self):
        self.pieces = []
        super(_BaseHTMLProcessor, self).reset()

    def _shorttag_replace(self, match):
        """
        :type match: Match[str]
        :rtype: str
        """

        tag = match.group(1)
        if tag in self.elements_no_end_tag:
            return '<' + tag + ' />'
        else:
            return '<' + tag + '></' + tag + '>'

    # By declaring these methods and overriding their compiled code
    # with the code from sgmllib, the original code will execute in
    # feedparser's scope instead of sgmllib's. This means that the
    # `tagfind` and `charref` regular expressions will be found as
    # they're declared above, not as they're declared in sgmllib.
    def goahead(self, i):
        raise NotImplementedError

    # Replace goahead with SGMLParser's goahead() code object.
    try:
        goahead.__code__ = sgmllib.SGMLParser.goahead.__code__
    except AttributeError:
        # Python 2
        # noinspection PyUnresolvedReferences
        goahead.func_code = sgmllib.SGMLParser.goahead.func_code

    def __parse_starttag(self, i):
        raise NotImplementedError

    # Replace __parse_starttag with SGMLParser's parse_starttag() code object.
    try:
        __parse_starttag.__code__ = sgmllib.SGMLParser.parse_starttag.__code__
    except AttributeError:
        # Python 2
        # noinspection PyUnresolvedReferences
        __parse_starttag.func_code = sgmllib.SGMLParser.parse_starttag.func_code

    def parse_starttag(self, i):
        j = self.__parse_starttag(i)
        if self._type == 'application/xhtml+xml':
            if j > 2 and self.rawdata[j-2:j] == '/>':
                self.unknown_endtag(self.lasttag)
        return j

    def feed(self, data):
        """
        :type data: str
        :rtype: None
        """

        data = re.sub(r"<!((?!DOCTYPE|--|\[))", r"&lt;!\1", data, flags=re.IGNORECASE)
        data = re.sub(r'<([^<>\s]+?)\s*/>', self._shorttag_replace, data)
        data = data.replace('&#39;', "'")
        data = data.replace('&#34;', '"')
        super(_BaseHTMLProcessor, self).feed(data)
        super(_BaseHTMLProcessor, self).close()

    @staticmethod
    def normalize_attrs(attrs):
        """
        :type attrs: List[Tuple[str, str]]
        :rtype: List[Tuple[str, str]]
        """

        if not attrs:
            return attrs
        # utility method to be called by descendants
        # Collapse any duplicate attribute names and values by converting
        # *attrs* into a dictionary, then convert it back to a list.
        attrs_d = {k.lower(): v for k, v in attrs}
        attrs = [
            (k, k in ('rel', 'type') and v.lower() or v)
            for k, v in attrs_d.items()
        ]
        attrs.sort()
        return attrs

    def unknown_starttag(self, tag, attrs):
        """
        :type tag: str
        :type attrs: List[Tuple[str, str]]
        :rtype: None
        """

        # Called for each start tag
        # attrs is a list of (attr, value) tuples
        # e.g. for <pre class='screen'>, tag='pre', attrs=[('class', 'screen')]
        uattrs = []
        strattrs = ''
        if attrs:
            for key, value in attrs:
                value = value.replace('>', '&gt;')
                value = value.replace('<', '&lt;')
                value = value.replace('"', '&quot;')
                value = self.bare_ampersand.sub("&amp;", value)
                uattrs.append((key, value))
            strattrs = ''.join(
                ' %s="%s"' % (key, value)
                for key, value in uattrs
            )
        if tag in self.elements_no_end_tag:
            self.pieces.append('<%s%s />' % (tag, strattrs))
        else:
            self.pieces.append('<%s%s>' % (tag, strattrs))

    def unknown_endtag(self, tag):
        """
        :type tag: str
        :rtype: None
        """

        # Called for each end tag, e.g. for </pre>, tag will be 'pre'
        # Reconstruct the original end tag.
        if tag not in self.elements_no_end_tag:
            self.pieces.append("</%s>" % tag)

    def handle_charref(self, ref):
        """
        :type ref: str
        :rtype: None
        """

        # Called for each character reference, e.g. '&#160;' will extract '160'
        # Reconstruct the original character reference.
        ref = ref.lower()
        if ref.startswith('x'):
            value = int(ref[1:], 16)
        else:
            value = int(ref)

        if value in _cp1252:
            self.pieces.append('&#%s;' % hex(ord(_cp1252[value]))[1:])
        else:
            self.pieces.append('&#%s;' % ref)

    def handle_entityref(self, ref):
        """
        :type ref: str
        :rtype: None
        """

        # Called for each entity reference, e.g. '&copy;' will extract 'copy'
        # Reconstruct the original entity reference.
        if ref in html.entities.name2codepoint or ref == 'apos':
            self.pieces.append('&%s;' % ref)
        else:
            self.pieces.append('&amp;%s' % ref)

    def handle_data(self, text):
        """
        :type text: str
        :rtype: None
        """

        # called for each block of plain text, i.e. outside of any tag and
        # not containing any character or entity references
        # Store the original text verbatim.
        self.pieces.append(text)

    def handle_comment(self, text):
        """
        :type text: str
        :rtype: None
        """

        # Called for HTML comments, e.g. <!-- insert Javascript code here -->
        # Reconstruct the original comment.
        self.pieces.append('<!--%s-->' % text)

    def handle_pi(self, text):
        """
        :type text: str
        :rtype: None
        """

        # Called for each processing instruction, e.g. <?instruction>
        # Reconstruct original processing instruction.
        self.pieces.append('<?%s>' % text)

    def handle_decl(self, text):
        """
        :type text: str
        :rtype: None
        """

        # called for the DOCTYPE, if present, e.g.
        # <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        #     "http://www.w3.org/TR/html4/loose.dtd">
        # Reconstruct original DOCTYPE
        self.pieces.append('<!%s>' % text)

    _new_declname_match = re.compile(r'[a-zA-Z][-_.a-zA-Z0-9:]*\s*').match

    def _scan_name(self, i, declstartpos):
        """
        :type i: int
        :type declstartpos: int
        :rtype: Tuple[Optional[str], int]
        """

        rawdata = self.rawdata
        n = len(rawdata)
        if i == n:
            return None, -1
        m = self._new_declname_match(rawdata, i)
        if m:
            s = m.group()
            name = s.strip()
            if (i + len(s)) == n:
                return None, -1  # end of buffer
            return name.lower(), m.end()
        else:
            self.handle_data(rawdata)
            # self.updatepos(declstartpos, i)
            return None, -1

    @staticmethod
    def convert_charref(name):
        """
        :type name: str
        :rtype: str
        """

        return '&#%s;' % name

    @staticmethod
    def convert_entityref(name):
        """
        :type name: str
        :rtype: str
        """

        return '&%s;' % name

    def output(self):
        """Return processed HTML as a single string.

        :rtype: str
        """

        return ''.join(self.pieces)

    def parse_declaration(self, i):
        """
        :type i: int
        :rtype: int
        """

        try:
            return sgmllib.SGMLParser.parse_declaration(self, i)
        except (AssertionError, sgmllib.SGMLParseError):
            # Escape the doctype declaration and continue parsing.
            self.handle_data('&lt;')
            return i+1
