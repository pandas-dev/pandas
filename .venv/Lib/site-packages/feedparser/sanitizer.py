# Copyright 2010-2023 Kurt McKee <contactme@kurtmckee.org>
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

import re

from .html import _BaseHTMLProcessor
from .urls import make_safe_absolute_uri


class _HTMLSanitizer(_BaseHTMLProcessor):
    acceptable_elements = {
        'a',
        'abbr',
        'acronym',
        'address',
        'area',
        'article',
        'aside',
        'audio',
        'b',
        'big',
        'blockquote',
        'br',
        'button',
        'canvas',
        'caption',
        'center',
        'cite',
        'code',
        'col',
        'colgroup',
        'command',
        'datagrid',
        'datalist',
        'dd',
        'del',
        'details',
        'dfn',
        'dialog',
        'dir',
        'div',
        'dl',
        'dt',
        'em',
        'event-source',
        'fieldset',
        'figcaption',
        'figure',
        'font',
        'footer',
        'form',
        'h1',
        'h2',
        'h3',
        'h4',
        'h5',
        'h6',
        'header',
        'hr',
        'i',
        'img',
        'input',
        'ins',
        'kbd',
        'keygen',
        'label',
        'legend',
        'li',
        'm',
        'map',
        'menu',
        'meter',
        'multicol',
        'nav',
        'nextid',
        'noscript',
        'ol',
        'optgroup',
        'option',
        'output',
        'p',
        'pre',
        'progress',
        'q',
        's',
        'samp',
        'section',
        'select',
        'small',
        'sound',
        'source',
        'spacer',
        'span',
        'strike',
        'strong',
        'sub',
        'sup',
        'table',
        'tbody',
        'td',
        'textarea',
        'tfoot',
        'th',
        'thead',
        'time',
        'tr',
        'tt',
        'u',
        'ul',
        'var',
        'video',
    }

    acceptable_attributes = {
        'abbr',
        'accept',
        'accept-charset',
        'accesskey',
        'action',
        'align',
        'alt',
        'autocomplete',
        'autofocus',
        'axis',
        'background',
        'balance',
        'bgcolor',
        'bgproperties',
        'border',
        'bordercolor',
        'bordercolordark',
        'bordercolorlight',
        'bottompadding',
        'cellpadding',
        'cellspacing',
        'ch',
        'challenge',
        'char',
        'charoff',
        'charset',
        'checked',
        'choff',
        'cite',
        'class',
        'clear',
        'color',
        'cols',
        'colspan',
        'compact',
        'contenteditable',
        'controls',
        'coords',
        'data',
        'datafld',
        'datapagesize',
        'datasrc',
        'datetime',
        'default',
        'delay',
        'dir',
        'disabled',
        'draggable',
        'dynsrc',
        'enctype',
        'end',
        'face',
        'for',
        'form',
        'frame',
        'galleryimg',
        'gutter',
        'headers',
        'height',
        'hidden',
        'hidefocus',
        'high',
        'href',
        'hreflang',
        'hspace',
        'icon',
        'id',
        'inputmode',
        'ismap',
        'keytype',
        'label',
        'lang',
        'leftspacing',
        'list',
        'longdesc',
        'loop',
        'loopcount',
        'loopend',
        'loopstart',
        'low',
        'lowsrc',
        'max',
        'maxlength',
        'media',
        'method',
        'min',
        'multiple',
        'name',
        'nohref',
        'noshade',
        'nowrap',
        'open',
        'optimum',
        'pattern',
        'ping',
        'point-size',
        'poster',
        'pqg',
        'preload',
        'prompt',
        'radiogroup',
        'readonly',
        'rel',
        'repeat-max',
        'repeat-min',
        'replace',
        'required',
        'rev',
        'rightspacing',
        'rows',
        'rowspan',
        'rules',
        'scope',
        'selected',
        'shape',
        'size',
        'span',
        'src',
        'start',
        'step',
        'style',
        'summary',
        'suppress',
        'tabindex',
        'target',
        'template',
        'title',
        'toppadding',
        'type',
        'unselectable',
        'urn',
        'usemap',
        'valign',
        'value',
        'variable',
        'volume',
        'vrml',
        'vspace',
        'width',
        'wrap',
        'xml:lang',
    }

    unacceptable_elements_with_end_tag = {
        'applet',
        'script',
        'style',
    }

    acceptable_css_properties = {
        'azimuth',
        'background-color',
        'border-bottom-color',
        'border-collapse',
        'border-color',
        'border-left-color',
        'border-right-color',
        'border-top-color',
        'clear',
        'color',
        'cursor',
        'direction',
        'display',
        'elevation',
        'float',
        'font',
        'font-family',
        'font-size',
        'font-style',
        'font-variant',
        'font-weight',
        'height',
        'letter-spacing',
        'line-height',
        'overflow',
        'pause',
        'pause-after',
        'pause-before',
        'pitch',
        'pitch-range',
        'richness',
        'speak',
        'speak-header',
        'speak-numeral',
        'speak-punctuation',
        'speech-rate',
        'stress',
        'text-align',
        'text-decoration',
        'text-indent',
        'unicode-bidi',
        'vertical-align',
        'voice-family',
        'volume',
        'white-space',
        'width',
    }

    # survey of common keywords found in feeds
    acceptable_css_keywords = {
        '!important',
        'aqua',
        'auto',
        'black',
        'block',
        'blue',
        'bold',
        'both',
        'bottom',
        'brown',
        'center',
        'collapse',
        'dashed',
        'dotted',
        'fuchsia',
        'gray',
        'green',
        'italic',
        'left',
        'lime',
        'maroon',
        'medium',
        'navy',
        'none',
        'normal',
        'nowrap',
        'olive',
        'pointer',
        'purple',
        'red',
        'right',
        'silver',
        'solid',
        'teal',
        'top',
        'transparent',
        'underline',
        'white',
        'yellow',
    }

    valid_css_values = re.compile(
        r'^('
        r'#[0-9a-f]+'  # Hex values
        r'|rgb\(\d+%?,\d*%?,?\d*%?\)?'  # RGB values
        r'|\d{0,2}\.?\d{0,2}(cm|em|ex|in|mm|pc|pt|px|%|,|\))?'  # Sizes/widths
        r')$'
    )

    mathml_elements = {
        'annotation',
        'annotation-xml',
        'maction',
        'maligngroup',
        'malignmark',
        'math',
        'menclose',
        'merror',
        'mfenced',
        'mfrac',
        'mglyph',
        'mi',
        'mlabeledtr',
        'mlongdiv',
        'mmultiscripts',
        'mn',
        'mo',
        'mover',
        'mpadded',
        'mphantom',
        'mprescripts',
        'mroot',
        'mrow',
        'ms',
        'mscarries',
        'mscarry',
        'msgroup',
        'msline',
        'mspace',
        'msqrt',
        'msrow',
        'mstack',
        'mstyle',
        'msub',
        'msubsup',
        'msup',
        'mtable',
        'mtd',
        'mtext',
        'mtr',
        'munder',
        'munderover',
        'none',
        'semantics',
    }

    mathml_attributes = {
        'accent',
        'accentunder',
        'actiontype',
        'align',
        'alignmentscope',
        'altimg',
        'altimg-height',
        'altimg-valign',
        'altimg-width',
        'alttext',
        'bevelled',
        'charalign',
        'close',
        'columnalign',
        'columnlines',
        'columnspacing',
        'columnspan',
        'columnwidth',
        'crossout',
        'decimalpoint',
        'denomalign',
        'depth',
        'dir',
        'display',
        'displaystyle',
        'edge',
        'encoding',
        'equalcolumns',
        'equalrows',
        'fence',
        'fontstyle',
        'fontweight',
        'form',
        'frame',
        'framespacing',
        'groupalign',
        'height',
        'href',
        'id',
        'indentalign',
        'indentalignfirst',
        'indentalignlast',
        'indentshift',
        'indentshiftfirst',
        'indentshiftlast',
        'indenttarget',
        'infixlinebreakstyle',
        'largeop',
        'length',
        'linebreak',
        'linebreakmultchar',
        'linebreakstyle',
        'lineleading',
        'linethickness',
        'location',
        'longdivstyle',
        'lquote',
        'lspace',
        'mathbackground',
        'mathcolor',
        'mathsize',
        'mathvariant',
        'maxsize',
        'minlabelspacing',
        'minsize',
        'movablelimits',
        'notation',
        'numalign',
        'open',
        'other',
        'overflow',
        'position',
        'rowalign',
        'rowlines',
        'rowspacing',
        'rowspan',
        'rquote',
        'rspace',
        'scriptlevel',
        'scriptminsize',
        'scriptsizemultiplier',
        'selection',
        'separator',
        'separators',
        'shift',
        'side',
        'src',
        'stackalign',
        'stretchy',
        'subscriptshift',
        'superscriptshift',
        'symmetric',
        'voffset',
        'width',
        'xlink:href',
        'xlink:show',
        'xlink:type',
        'xmlns',
        'xmlns:xlink',
    }

    # svgtiny - foreignObject + linearGradient + radialGradient + stop
    svg_elements = {
        'a',
        'animate',
        'animateColor',
        'animateMotion',
        'animateTransform',
        'circle',
        'defs',
        'desc',
        'ellipse',
        'font-face',
        'font-face-name',
        'font-face-src',
        'foreignObject',
        'g',
        'glyph',
        'hkern',
        'line',
        'linearGradient',
        'marker',
        'metadata',
        'missing-glyph',
        'mpath',
        'path',
        'polygon',
        'polyline',
        'radialGradient',
        'rect',
        'set',
        'stop',
        'svg',
        'switch',
        'text',
        'title',
        'tspan',
        'use',
    }

    # svgtiny + class + opacity + offset + xmlns + xmlns:xlink
    svg_attributes = {
        'accent-height',
        'accumulate',
        'additive',
        'alphabetic',
        'arabic-form',
        'ascent',
        'attributeName',
        'attributeType',
        'baseProfile',
        'bbox',
        'begin',
        'by',
        'calcMode',
        'cap-height',
        'class',
        'color',
        'color-rendering',
        'content',
        'cx',
        'cy',
        'd',
        'descent',
        'display',
        'dur',
        'dx',
        'dy',
        'end',
        'fill',
        'fill-opacity',
        'fill-rule',
        'font-family',
        'font-size',
        'font-stretch',
        'font-style',
        'font-variant',
        'font-weight',
        'from',
        'fx',
        'fy',
        'g1',
        'g2',
        'glyph-name',
        'gradientUnits',
        'hanging',
        'height',
        'horiz-adv-x',
        'horiz-origin-x',
        'id',
        'ideographic',
        'k',
        'keyPoints',
        'keySplines',
        'keyTimes',
        'lang',
        'marker-end',
        'marker-mid',
        'marker-start',
        'markerHeight',
        'markerUnits',
        'markerWidth',
        'mathematical',
        'max',
        'min',
        'name',
        'offset',
        'opacity',
        'orient',
        'origin',
        'overline-position',
        'overline-thickness',
        'panose-1',
        'path',
        'pathLength',
        'points',
        'preserveAspectRatio',
        'r',
        'refX',
        'refY',
        'repeatCount',
        'repeatDur',
        'requiredExtensions',
        'requiredFeatures',
        'restart',
        'rotate',
        'rx',
        'ry',
        'slope',
        'stemh',
        'stemv',
        'stop-color',
        'stop-opacity',
        'strikethrough-position',
        'strikethrough-thickness',
        'stroke',
        'stroke-dasharray',
        'stroke-dashoffset',
        'stroke-linecap',
        'stroke-linejoin',
        'stroke-miterlimit',
        'stroke-opacity',
        'stroke-width',
        'systemLanguage',
        'target',
        'text-anchor',
        'to',
        'transform',
        'type',
        'u1',
        'u2',
        'underline-position',
        'underline-thickness',
        'unicode',
        'unicode-range',
        'units-per-em',
        'values',
        'version',
        'viewBox',
        'visibility',
        'width',
        'widths',
        'x',
        'x-height',
        'x1',
        'x2',
        'xlink:actuate',
        'xlink:arcrole',
        'xlink:href',
        'xlink:role',
        'xlink:show',
        'xlink:title',
        'xlink:type',
        'xml:base',
        'xml:lang',
        'xml:space',
        'xmlns',
        'xmlns:xlink',
        'y',
        'y1',
        'y2',
        'zoomAndPan',
    }

    svg_attr_map = None
    svg_elem_map = None

    acceptable_svg_properties = {
        'fill',
        'fill-opacity',
        'fill-rule',
        'stroke',
        'stroke-linecap',
        'stroke-linejoin',
        'stroke-opacity',
        'stroke-width',
    }

    def __init__(self, encoding=None, _type='application/xhtml+xml'):
        super(_HTMLSanitizer, self).__init__(encoding, _type)

        self.unacceptablestack = 0
        self.mathmlOK = 0
        self.svgOK = 0

    def reset(self):
        super(_HTMLSanitizer, self).reset()
        self.unacceptablestack = 0
        self.mathmlOK = 0
        self.svgOK = 0

    def unknown_starttag(self, tag, attrs):
        acceptable_attributes = self.acceptable_attributes
        keymap = {}
        if tag not in self.acceptable_elements or self.svgOK:
            if tag in self.unacceptable_elements_with_end_tag:
                self.unacceptablestack += 1

            # add implicit namespaces to html5 inline svg/mathml
            if self._type.endswith('html'):
                if not dict(attrs).get('xmlns'):
                    if tag == 'svg':
                        attrs.append(('xmlns', 'http://www.w3.org/2000/svg'))
                    if tag == 'math':
                        attrs.append(('xmlns', 'http://www.w3.org/1998/Math/MathML'))

            # not otherwise acceptable, perhaps it is MathML or SVG?
            if tag == 'math' and ('xmlns', 'http://www.w3.org/1998/Math/MathML') in attrs:
                self.mathmlOK += 1
            if tag == 'svg' and ('xmlns', 'http://www.w3.org/2000/svg') in attrs:
                self.svgOK += 1

            # chose acceptable attributes based on tag class, else bail
            if self.mathmlOK and tag in self.mathml_elements:
                acceptable_attributes = self.mathml_attributes
            elif self.svgOK and tag in self.svg_elements:
                # For most vocabularies, lowercasing is a good idea. Many
                # svg elements, however, are camel case.
                if not self.svg_attr_map:
                    lower = [attr.lower() for attr in self.svg_attributes]
                    mix = [a for a in self.svg_attributes if a not in lower]
                    self.svg_attributes = lower
                    self.svg_attr_map = {a.lower(): a for a in mix}

                    lower = [attr.lower() for attr in self.svg_elements]
                    mix = [a for a in self.svg_elements if a not in lower]
                    self.svg_elements = lower
                    self.svg_elem_map = {a.lower(): a for a in mix}
                acceptable_attributes = self.svg_attributes
                tag = self.svg_elem_map.get(tag, tag)
                keymap = self.svg_attr_map
            elif tag not in self.acceptable_elements:
                return

        # declare xlink namespace, if needed
        if self.mathmlOK or self.svgOK:
            if any((a for a in attrs if a[0].startswith('xlink:'))):
                if not ('xmlns:xlink', 'http://www.w3.org/1999/xlink') in attrs:
                    attrs.append(('xmlns:xlink', 'http://www.w3.org/1999/xlink'))

        clean_attrs = []
        for key, value in self.normalize_attrs(attrs):
            if key == 'style' and 'style' in acceptable_attributes:
                clean_value = self.sanitize_style(value)
                if clean_value:
                    clean_attrs.append((key, clean_value))
            elif key in acceptable_attributes:
                key = keymap.get(key, key)
                # make sure the uri uses an acceptable uri scheme
                if key == 'href':
                    value = make_safe_absolute_uri(value)
                clean_attrs.append((key, value))
        super(_HTMLSanitizer, self).unknown_starttag(tag, clean_attrs)

    def unknown_endtag(self, tag):
        if tag not in self.acceptable_elements:
            if tag in self.unacceptable_elements_with_end_tag:
                self.unacceptablestack -= 1
            if self.mathmlOK and tag in self.mathml_elements:
                if tag == 'math' and self.mathmlOK:
                    self.mathmlOK -= 1
            elif self.svgOK and tag in self.svg_elements:
                tag = self.svg_elem_map.get(tag, tag)
                if tag == 'svg' and self.svgOK:
                    self.svgOK -= 1
            else:
                return
        super(_HTMLSanitizer, self).unknown_endtag(tag)

    def handle_pi(self, text):
        pass

    def handle_decl(self, text):
        pass

    def handle_data(self, text):
        if not self.unacceptablestack:
            super(_HTMLSanitizer, self).handle_data(text)

    def sanitize_style(self, style):
        # disallow urls
        style = re.compile(r'url\s*\(\s*[^\s)]+?\s*\)\s*').sub(' ', style)

        # gauntlet
        if not re.match(r"""^([:,;#%.\sa-zA-Z0-9!]|\w-\w|'[\s\w]+'|"[\s\w]+"|\([\d,\s]+\))*$""", style):
            return ''
        # This replaced a regexp that used re.match and was prone to
        # pathological back-tracking.
        if re.sub(r"\s*[-\w]+\s*:\s*[^:;]*;?", '', style).strip():
            return ''

        clean = []
        for prop, value in re.findall(r"([-\w]+)\s*:\s*([^:;]*)", style):
            if not value:
                continue
            if prop.lower() in self.acceptable_css_properties:
                clean.append(prop + ': ' + value + ';')
            elif prop.split('-')[0].lower() in ['background', 'border', 'margin', 'padding']:
                for keyword in value.split():
                    if (
                            keyword not in self.acceptable_css_keywords
                            and not self.valid_css_values.match(keyword)
                    ):
                        break
                else:
                    clean.append(prop + ': ' + value + ';')
            elif self.svgOK and prop.lower() in self.acceptable_svg_properties:
                clean.append(prop + ': ' + value + ';')

        return ' '.join(clean)

    def parse_comment(self, i, report=1):
        ret = super(_HTMLSanitizer, self).parse_comment(i, report)
        if ret >= 0:
            return ret
        # if ret == -1, this may be a malicious attempt to circumvent
        # sanitization, or a page-destroying unclosed comment
        match = re.compile(r'--[^>]*>').search(self.rawdata, i+4)
        if match:
            return match.end()
        # unclosed comment; deliberately fail to handle_data()
        return len(self.rawdata)


def _sanitize_html(html_source, encoding, _type):
    p = _HTMLSanitizer(encoding, _type)
    html_source = html_source.replace('<![CDATA[', '&lt;![CDATA[')
    p.feed(html_source)
    data = p.output()
    data = data.strip().replace('\r\n', '\n')
    return data


# Match XML entity declarations.
# Example: <!ENTITY copyright "(C)">
RE_ENTITY_PATTERN = re.compile(br'^\s*<!ENTITY([^>]*?)>', re.MULTILINE)

# Match XML DOCTYPE declarations.
# Example: <!DOCTYPE feed [ ]>
RE_DOCTYPE_PATTERN = re.compile(br'^\s*<!DOCTYPE([^>]*?)>', re.MULTILINE)

# Match safe entity declarations.
# This will allow hexadecimal character references through,
# as well as text, but not arbitrary nested entities.
# Example: cubed "&#179;"
# Example: copyright "(C)"
# Forbidden: explode1 "&explode2;&explode2;"
RE_SAFE_ENTITY_PATTERN = re.compile(br'\s+(\w+)\s+"(&#\w+;|[^&"]*)"')


def replace_doctype(data):
    """Strips and replaces the DOCTYPE, returns (rss_version, stripped_data)

    rss_version may be 'rss091n' or None
    stripped_data is the same XML document with a replaced DOCTYPE
    """

    # Divide the document into two groups by finding the location
    # of the first element that doesn't begin with '<?' or '<!'.
    start = re.search(br'<\w', data)
    start = start and start.start() or -1
    head, data = data[:start+1], data[start+1:]

    # Save and then remove all of the ENTITY declarations.
    entity_results = RE_ENTITY_PATTERN.findall(head)
    head = RE_ENTITY_PATTERN.sub(b'', head)

    # Find the DOCTYPE declaration and check the feed type.
    doctype_results = RE_DOCTYPE_PATTERN.findall(head)
    doctype = doctype_results and doctype_results[0] or b''
    if b'netscape' in doctype.lower():
        version = 'rss091n'
    else:
        version = None

    # Re-insert the safe ENTITY declarations if a DOCTYPE was found.
    replacement = b''
    if len(doctype_results) == 1 and entity_results:
        safe_entities = [
            e
            for e in entity_results
            if RE_SAFE_ENTITY_PATTERN.match(e)
        ]
        if safe_entities:
            replacement = b'<!DOCTYPE feed [\n<!ENTITY' \
                        + b'>\n<!ENTITY '.join(safe_entities) \
                        + b'>\n]>'
    data = RE_DOCTYPE_PATTERN.sub(replacement, head) + data

    # Precompute the safe entities for the loose parser.
    safe_entities = {
        k.decode('utf-8'): v.decode('utf-8')
        for k, v in RE_SAFE_ENTITY_PATTERN.findall(replacement)
    }
    return version, data, safe_entities
