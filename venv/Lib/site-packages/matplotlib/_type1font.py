"""
A class representing a Type 1 font.

This version reads pfa and pfb files and splits them for embedding in
pdf files. It also supports SlantFont and ExtendFont transformations,
similarly to pdfTeX and friends.

Usage::

    font = Type1Font(filename)
    clear_part, encrypted_part, finale = font.parts
    slanted_font = font.transform({'slant': 0.167})
    extended_font = font.transform({'extend': 1.2})
    subset_font = font.subset([ord(c) for c in 'Hello World'])

Sources:

* Adobe Technical Note #5040, Supporting Downloadable PostScript
  Language Fonts.

* Adobe Type 1 Font Format, Adobe Systems Incorporated, third printing,
  v1.1, 1993. ISBN 0-201-57044-0.
"""

from __future__ import annotations

import binascii
import functools
import itertools
import logging
import re
import string
import struct
import typing as T

import numpy as np

from matplotlib.cbook import _format_approx
from . import _api

_log = logging.getLogger(__name__)


class _Token:
    """
    A token in a PostScript stream.

    Attributes
    ----------
    pos : int
        Position, i.e. offset from the beginning of the data.
    raw : str
        Raw text of the token.
    kind : str
        Description of the token (for debugging or testing).
    """
    __slots__ = ('pos', 'raw')
    kind = '?'

    def __init__(self, pos, raw):
        _log.debug('type1font._Token %s at %d: %r', self.kind, pos, raw)
        self.pos = pos
        self.raw = raw

    def __str__(self):
        return f"<{self.kind} {self.raw} @{self.pos}>"

    def endpos(self):
        """Position one past the end of the token"""
        return self.pos + len(self.raw)

    def is_keyword(self, *names):
        """Is this a name token with one of the names?"""
        return False

    def is_slash_name(self):
        """Is this a name token that starts with a slash?"""
        return False

    def is_delim(self):
        """Is this a delimiter token?"""
        return False

    def is_number(self):
        """Is this a number token?"""
        return False

    def value(self):
        return self.raw


class _NameToken(_Token):
    kind = 'name'

    def is_slash_name(self):
        return self.raw.startswith('/')

    def value(self):
        return self.raw[1:]


class _BooleanToken(_Token):
    kind = 'boolean'

    def value(self):
        return self.raw == 'true'


class _KeywordToken(_Token):
    kind = 'keyword'

    def is_keyword(self, *names):
        return self.raw in names


class _DelimiterToken(_Token):
    kind = 'delimiter'

    def is_delim(self):
        return True

    def opposite(self):
        return {'[': ']', ']': '[',
                '{': '}', '}': '{',
                '<<': '>>', '>>': '<<'
                }[self.raw]


class _WhitespaceToken(_Token):
    kind = 'whitespace'


class _StringToken(_Token):
    kind = 'string'
    _escapes_re = re.compile(r'\\([\\()nrtbf]|[0-7]{1,3})')
    _replacements = {'\\': '\\', '(': '(', ')': ')', 'n': '\n',
                     'r': '\r', 't': '\t', 'b': '\b', 'f': '\f'}
    _ws_re = re.compile('[\0\t\r\f\n ]')

    @classmethod
    def _escape(cls, match):
        group = match.group(1)
        try:
            return cls._replacements[group]
        except KeyError:
            return chr(int(group, 8))

    @functools.lru_cache
    def value(self):
        if self.raw[0] == '(':
            return self._escapes_re.sub(self._escape, self.raw[1:-1])
        else:
            data = self._ws_re.sub('', self.raw[1:-1])
            if len(data) % 2 == 1:
                data += '0'
            return binascii.unhexlify(data)


class _BinaryToken(_Token):
    kind = 'binary'

    def value(self):
        return self.raw[1:]


class _NumberToken(_Token):
    kind = 'number'

    def is_number(self):
        return True

    def value(self):
        if '.' not in self.raw:
            return int(self.raw)
        else:
            return float(self.raw)


def _tokenize(data: bytes, skip_ws: bool) -> T.Generator[_Token, int, None]:
    """
    A generator that produces _Token instances from Type-1 font code.

    The consumer of the generator may send an integer to the tokenizer to
    indicate that the next token should be _BinaryToken of the given length.

    Parameters
    ----------
    data : bytes
        The data of the font to tokenize.

    skip_ws : bool
        If true, the generator will drop any _WhitespaceTokens from the output.
    """

    text = data.decode('ascii', 'replace')
    whitespace_or_comment_re = re.compile(r'[\0\t\r\f\n ]+|%[^\r\n]*')
    token_re = re.compile(r'/{0,2}[^]\0\t\r\f\n ()<>{}/%[]+')
    instring_re = re.compile(r'[()\\]')
    hex_re = re.compile(r'^<[0-9a-fA-F\0\t\r\f\n ]*>$')
    oct_re = re.compile(r'[0-7]{1,3}')
    pos = 0
    next_binary: int | None = None

    while pos < len(text):
        if next_binary is not None:
            n = next_binary
            next_binary = (yield _BinaryToken(pos, data[pos:pos+n]))
            pos += n
            continue
        match = whitespace_or_comment_re.match(text, pos)
        if match:
            if not skip_ws:
                next_binary = (yield _WhitespaceToken(pos, match.group()))
            pos = match.end()
        elif text[pos] == '(':
            # PostScript string rules:
            # - parentheses must be balanced
            # - backslashes escape backslashes and parens
            # - also codes \n\r\t\b\f and octal escapes are recognized
            # - other backslashes do not escape anything
            start = pos
            pos += 1
            depth = 1
            while depth:
                match = instring_re.search(text, pos)
                if match is None:
                    raise ValueError(
                        f'Unterminated string starting at {start}')
                pos = match.end()
                if match.group() == '(':
                    depth += 1
                elif match.group() == ')':
                    depth -= 1
                else:  # a backslash
                    char = text[pos]
                    if char in r'\()nrtbf':
                        pos += 1
                    else:
                        octal = oct_re.match(text, pos)
                        if octal:
                            pos = octal.end()
                        else:
                            pass  # non-escaping backslash
            next_binary = (yield _StringToken(start, text[start:pos]))
        elif text[pos:pos + 2] in ('<<', '>>'):
            next_binary = (yield _DelimiterToken(pos, text[pos:pos + 2]))
            pos += 2
        elif text[pos] == '<':
            start = pos
            try:
                pos = text.index('>', pos) + 1
            except ValueError as e:
                raise ValueError(f'Unterminated hex string starting at {start}'
                                 ) from e
            if not hex_re.match(text[start:pos]):
                raise ValueError(f'Malformed hex string starting at {start}')
            next_binary = (yield _StringToken(pos, text[start:pos]))
        else:
            match = token_re.match(text, pos)
            if match:
                raw = match.group()
                if raw.startswith('/'):
                    next_binary = (yield _NameToken(pos, raw))
                elif match.group() in ('true', 'false'):
                    next_binary = (yield _BooleanToken(pos, raw))
                else:
                    try:
                        float(raw)
                        next_binary = (yield _NumberToken(pos, raw))
                    except ValueError:
                        next_binary = (yield _KeywordToken(pos, raw))
                pos = match.end()
            else:
                next_binary = (yield _DelimiterToken(pos, text[pos]))
                pos += 1


class _BalancedExpression(_Token):
    pass


def _expression(initial, tokens, data):
    """
    Consume some number of tokens and return a balanced PostScript expression.

    Parameters
    ----------
    initial : _Token
        The token that triggered parsing a balanced expression.
    tokens : iterator of _Token
        Following tokens.
    data : bytes
        Underlying data that the token positions point to.

    Returns
    -------
    _BalancedExpression
    """
    delim_stack = []
    token = initial
    while True:
        if token.is_delim():
            if token.raw in ('[', '{'):
                delim_stack.append(token)
            elif token.raw in (']', '}'):
                if not delim_stack:
                    raise RuntimeError(f"unmatched closing token {token}")
                match = delim_stack.pop()
                if match.raw != token.opposite():
                    raise RuntimeError(
                        f"opening token {match} closed by {token}"
                    )
                if not delim_stack:
                    break
            else:
                raise RuntimeError(f'unknown delimiter {token}')
        elif not delim_stack:
            break
        token = next(tokens)
    return _BalancedExpression(
        initial.pos,
        data[initial.pos:token.endpos()].decode('ascii', 'replace')
    )


class Type1Font:
    """
    A class representing a Type-1 font, for use by backends.

    Attributes
    ----------
    parts : tuple
        A 3-tuple of the cleartext part, the encrypted part, and the finale of
        zeros.

    decrypted : bytes
        The decrypted form of ``parts[1]``.

    prop : dict[str, Any]
        A dictionary of font properties. Noteworthy keys include:

        - FontName: PostScript name of the font
        - Encoding: dict from numeric codes to glyph names
        - FontMatrix: bytes object encoding a matrix
        - UniqueID: optional font identifier, dropped when modifying the font
        - CharStrings: dict from glyph names to byte code
        - Subrs: array of byte code subroutines
        - OtherSubrs: bytes object encoding some PostScript code
    """
    __slots__ = ('parts', 'decrypted', 'prop', '_pos', '_abbr')
    # the _pos dict contains (begin, end) indices to parts[0] + decrypted
    # so that they can be replaced when transforming the font;
    # but since sometimes a definition appears in both parts[0] and decrypted,
    # _pos[name] is an array of such pairs
    #
    # _abbr maps three standard abbreviations to their particular names in
    # this font (e.g. 'RD' is named '-|' in some fonts)

    def __init__(self, input):
        """
        Initialize a Type-1 font.

        Parameters
        ----------
        input : str or 3-tuple
            Either a pfb file name, or a 3-tuple of already-decoded Type-1
            font `~.Type1Font.parts`.
        """
        if isinstance(input, tuple) and len(input) == 3:
            self.parts = input
        else:
            with open(input, 'rb') as file:
                data = self._read(file)
            self.parts = self._split(data)

        self.decrypted = self._decrypt(self.parts[1], 'eexec')
        self._abbr = {'RD': 'RD', 'ND': 'ND', 'NP': 'NP'}
        self._parse()

    def _read(self, file):
        """Read the font from a file, decoding into usable parts."""
        rawdata = file.read()
        if not rawdata.startswith(b'\x80'):
            return rawdata

        data = b''
        while rawdata:
            if not rawdata.startswith(b'\x80'):
                raise RuntimeError('Broken pfb file (expected byte 128, '
                                   'got %d)' % rawdata[0])
            type = rawdata[1]
            if type in (1, 2):
                length, = struct.unpack('<i', rawdata[2:6])
                segment = rawdata[6:6 + length]
                rawdata = rawdata[6 + length:]

            if type == 1:       # ASCII text: include verbatim
                data += segment
            elif type == 2:     # binary data: encode in hexadecimal
                data += binascii.hexlify(segment)
            elif type == 3:     # end of file
                break
            else:
                raise RuntimeError('Unknown segment type %d in pfb file' % type)

        return data

    def _split(self, data):
        """
        Split the Type 1 font into its three main parts.

        The three parts are: (1) the cleartext part, which ends in a
        eexec operator; (2) the encrypted part; (3) the fixed part,
        which contains 512 ASCII zeros possibly divided on various
        lines, a cleartomark operator, and possibly something else.
        """

        # Cleartext part: just find the eexec and skip whitespace
        idx = data.index(b'eexec')
        idx += len(b'eexec')
        while data[idx] in b' \t\r\n':
            idx += 1
        len1 = idx

        # Encrypted part: find the cleartomark operator and count
        # zeros backward
        idx = data.rindex(b'cleartomark') - 1
        zeros = 512
        while zeros and data[idx] in b'0' or data[idx] in b'\r\n':
            if data[idx] in b'0':
                zeros -= 1
            idx -= 1
        if zeros:
            # this may have been a problem on old implementations that
            # used the zeros as necessary padding
            _log.info('Insufficiently many zeros in Type 1 font')

        # Convert encrypted part to binary (if we read a pfb file, we may end
        # up converting binary to hexadecimal to binary again; but if we read
        # a pfa file, this part is already in hex, and I am not quite sure if
        # even the pfb format guarantees that it will be in binary).
        idx1 = len1 + ((idx - len1 + 2) & ~1)  # ensure an even number of bytes
        binary = binascii.unhexlify(data[len1:idx1])

        return data[:len1], binary, data[idx+1:]

    @staticmethod
    def _decrypt(ciphertext, key, ndiscard=4):
        """
        Decrypt ciphertext using the Type-1 font algorithm.

        The algorithm is described in Adobe's "Adobe Type 1 Font Format".
        The key argument can be an integer, or one of the strings
        'eexec' and 'charstring', which map to the key specified for the
        corresponding part of Type-1 fonts.

        The ndiscard argument should be an integer, usually 4.
        That number of bytes is discarded from the beginning of plaintext.
        """

        key = _api.getitem_checked({'eexec': 55665, 'charstring': 4330}, key=key)
        plaintext = []
        for byte in ciphertext:
            plaintext.append(byte ^ (key >> 8))
            key = ((key+byte) * 52845 + 22719) & 0xffff

        return bytes(plaintext[ndiscard:])

    @staticmethod
    def _encrypt(plaintext, key, ndiscard=4):
        """
        Encrypt plaintext using the Type-1 font algorithm.

        The algorithm is described in Adobe's "Adobe Type 1 Font Format".
        The key argument can be an integer, or one of the strings
        'eexec' and 'charstring', which map to the key specified for the
        corresponding part of Type-1 fonts.

        The ndiscard argument should be an integer, usually 4. That
        number of bytes is prepended to the plaintext before encryption.
        This function prepends NUL bytes for reproducibility, even though
        the original algorithm uses random bytes, presumably to avoid
        cryptanalysis.
        """

        key = _api.getitem_checked({'eexec': 55665, 'charstring': 4330}, key=key)
        ciphertext = []
        for byte in b'\0' * ndiscard + plaintext:
            c = byte ^ (key >> 8)
            ciphertext.append(c)
            key = ((key + c) * 52845 + 22719) & 0xffff

        return bytes(ciphertext)

    def _parse(self):
        """
        Find the values of various font properties. This limited kind
        of parsing is described in Chapter 10 "Adobe Type Manager
        Compatibility" of the Type-1 spec.
        """
        # Start with reasonable defaults
        prop = {'Weight': 'Regular', 'ItalicAngle': 0.0, 'isFixedPitch': False,
                'UnderlinePosition': -100, 'UnderlineThickness': 50}
        pos = {}
        data = self.parts[0] + self.decrypted

        source = _tokenize(data, True)
        while True:
            # See if there is a key to be assigned a value
            # e.g. /FontName in /FontName /Helvetica def
            try:
                token = next(source)
            except StopIteration:
                break
            if token.is_delim():
                # skip over this - we want top-level keys only
                _expression(token, source, data)
            if token.is_slash_name():
                key = token.value()
                keypos = token.pos
            else:
                continue

            # Some values need special parsing
            if key in ('Subrs', 'CharStrings', 'Encoding', 'OtherSubrs'):
                prop[key], endpos = {
                    'Subrs': self._parse_subrs,
                    'CharStrings': self._parse_charstrings,
                    'Encoding': self._parse_encoding,
                    'OtherSubrs': self._parse_othersubrs
                }[key](source, data)
                pos.setdefault(key, []).append((keypos, endpos))
                continue

            try:
                token = next(source)
            except StopIteration:
                break

            if isinstance(token, _KeywordToken):
                # constructs like
                # FontDirectory /Helvetica known {...} {...} ifelse
                # mean the key was not really a key
                continue

            if token.is_delim():
                value = _expression(token, source, data).raw
            else:
                value = token.value()

            # look for a 'def' possibly preceded by access modifiers
            try:
                kw = next(
                    kw for kw in source
                    if not kw.is_keyword('readonly', 'noaccess', 'executeonly')
                )
            except StopIteration:
                break

            # sometimes noaccess def and readonly def are abbreviated
            if kw.is_keyword('def', self._abbr['ND'], self._abbr['NP']):
                prop[key] = value
                pos.setdefault(key, []).append((keypos, kw.endpos()))

            # detect the standard abbreviations
            if value == '{noaccess def}':
                self._abbr['ND'] = key
            elif value == '{noaccess put}':
                self._abbr['NP'] = key
            elif value == '{string currentfile exch readstring pop}':
                self._abbr['RD'] = key

        # Fill in the various *Name properties
        if 'FontName' not in prop:
            prop['FontName'] = (prop.get('FullName') or
                                prop.get('FamilyName') or
                                'Unknown')
        if 'FullName' not in prop:
            prop['FullName'] = prop['FontName']
        if 'FamilyName' not in prop:
            extras = ('(?i)([ -](regular|plain|italic|oblique|(semi)?bold|'
                      '(ultra)?light|extra|condensed))+$')
            prop['FamilyName'] = re.sub(extras, '', prop['FullName'])

        # Parse FontBBox
        toks = [*_tokenize(prop['FontBBox'].encode('ascii'), True)]
        if ([tok.kind for tok in toks]
                != ['delimiter', 'number', 'number', 'number', 'number', 'delimiter']
                or toks[-1].raw != toks[0].opposite()):
            raise RuntimeError(
                f"FontBBox should be a size-4 array, was {prop['FontBBox']}")
        prop['FontBBox'] = [tok.value() for tok in toks[1:-1]]

        # Decrypt the encrypted parts
        ndiscard = prop.get('lenIV', 4)
        cs = prop['CharStrings']
        for key, value in cs.items():
            cs[key] = self._decrypt(value, 'charstring', ndiscard)
        if 'Subrs' in prop:
            prop['Subrs'] = [
                self._decrypt(value, 'charstring', ndiscard)
                for value in prop['Subrs']
            ]

        self.prop = prop
        self._pos = pos

    def _parse_subrs(self, tokens, _data):
        count_token = next(tokens)
        if not count_token.is_number():
            raise RuntimeError(
                f"Token following /Subrs must be a number, was {count_token}"
            )
        count = count_token.value()
        array = [None] * count
        next(t for t in tokens if t.is_keyword('array'))
        for _ in range(count):
            next(t for t in tokens if t.is_keyword('dup'))
            index_token = next(tokens)
            if not index_token.is_number():
                raise RuntimeError(
                    "Token following dup in Subrs definition must be a "
                    f"number, was {index_token}"
                )
            nbytes_token = next(tokens)
            if not nbytes_token.is_number():
                raise RuntimeError(
                    "Second token following dup in Subrs definition must "
                    f"be a number, was {nbytes_token}"
                )
            token = next(tokens)
            if not token.is_keyword(self._abbr['RD']):
                raise RuntimeError(
                    f"Token preceding subr must be {self._abbr['RD']}, "
                    f"was {token}"
                )
            binary_token = tokens.send(1+nbytes_token.value())
            array[index_token.value()] = binary_token.value()

        return array, next(tokens).endpos()

    def _parse_charstrings(self, tokens, _data):
        count_token = next(tokens)
        if not count_token.is_number():
            raise RuntimeError(
                "Token following /CharStrings must be a number, "
                f"was {count_token}"
            )
        count = count_token.value()
        charstrings = {}
        next(t for t in tokens if t.is_keyword('begin'))
        while True:
            token = next(t for t in tokens
                         if t.is_keyword('end') or t.is_slash_name())
            if token.raw == 'end':
                return charstrings, token.endpos()
            glyphname = token.value()
            nbytes_token = next(tokens)
            if not nbytes_token.is_number():
                raise RuntimeError(
                    f"Token following /{glyphname} in CharStrings definition "
                    f"must be a number, was {nbytes_token}"
                )
            token = next(tokens)
            if not token.is_keyword(self._abbr['RD']):
                raise RuntimeError(
                    f"Token preceding charstring must be {self._abbr['RD']}, "
                    f"was {token}"
                )
            binary_token = tokens.send(1+nbytes_token.value())
            charstrings[glyphname] = binary_token.value()

    @staticmethod
    def _parse_encoding(tokens, _data):
        # this only works for encodings that follow the Adobe manual
        # but some old fonts include non-compliant data - we log a warning
        # and return a possibly incomplete encoding
        encoding = {}
        while True:
            token = next(t for t in tokens
                         if t.is_keyword('StandardEncoding', 'dup', 'def'))
            if token.is_keyword('StandardEncoding'):
                return _StandardEncoding, token.endpos()
            if token.is_keyword('def'):
                return encoding, token.endpos()
            index_token = next(tokens)
            if not index_token.is_number():
                _log.warning(
                    f"Parsing encoding: expected number, got {index_token}"
                )
                continue
            name_token = next(tokens)
            if not name_token.is_slash_name():
                _log.warning(
                    f"Parsing encoding: expected slash-name, got {name_token}"
                )
                continue
            encoding[index_token.value()] = name_token.value()

    def _parse_othersubrs(self, tokens, data):
        init_pos = None
        while True:
            token = next(tokens)
            if init_pos is None:
                init_pos = token.pos
            if token.is_delim():
                _expression(token, tokens, data)
            elif token.is_keyword('def', self._abbr['ND']):
                return data[init_pos:token.endpos()], token.endpos()

    def transform(self, effects):
        """
        Return a new font that is slanted and/or extended.

        Parameters
        ----------
        effects : dict
            A dict with optional entries:

            - 'slant' : float, default: 0
                Tangent of the angle that the font is to be slanted to the
                right. Negative values slant to the left.
            - 'extend' : float, default: 1
                Scaling factor for the font width. Values less than 1 condense
                the glyphs.

        Returns
        -------
        `Type1Font`
        """
        fontname = self.prop['FontName']
        italicangle = self.prop['ItalicAngle']

        array = [
            float(x) for x in (self.prop['FontMatrix']
                               .lstrip('[').rstrip(']').split())
        ]
        oldmatrix = np.eye(3, 3)
        oldmatrix[0:3, 0] = array[::2]
        oldmatrix[0:3, 1] = array[1::2]
        modifier = np.eye(3, 3)

        if 'slant' in effects:
            slant = effects['slant']
            fontname += f'_Slant_{int(1000 * slant)}'
            italicangle = round(
                float(italicangle) - np.arctan(slant) / np.pi * 180,
                5
            )
            modifier[1, 0] = slant

        if 'extend' in effects:
            extend = effects['extend']
            fontname += f'_Extend_{int(1000 * extend)}'
            modifier[0, 0] = extend

        newmatrix = np.dot(modifier, oldmatrix)
        array[::2] = newmatrix[0:3, 0]
        array[1::2] = newmatrix[0:3, 1]
        fontmatrix = (
            f"[{' '.join(_format_approx(x, 6) for x in array)}]"
        )
        newparts = self._replace(
            [(x, f'/FontName/{fontname} def')
             for x in self._pos['FontName']]
            + [(x, f'/ItalicAngle {italicangle} def')
               for x in self._pos['ItalicAngle']]
            + [(x, f'/FontMatrix {fontmatrix} readonly def')
               for x in self._pos['FontMatrix']]
            + [(x, '') for x in self._pos.get('UniqueID', [])]
        )

        return Type1Font((
            newparts[0],
            self._encrypt(newparts[1], 'eexec'),
            self.parts[2]
        ))

    def with_encoding(self, encoding):
        """
        Change the encoding of the font.

        Parameters
        ----------
        encoding : dict
            A dictionary mapping character codes to glyph names.

        Returns
        -------
        `Type1Font`
        """
        newparts = self._replace(
            [(x, '') for x in self._pos.get('UniqueID', [])]
            + [(self._pos['Encoding'][0], self._postscript_encoding(encoding))]
        )
        return Type1Font((
            newparts[0],
            self._encrypt(newparts[1], 'eexec'),
            self.parts[2]
        ))

    def _replace(self, replacements):
        """
        Change the font according to `replacements`

        Parameters
        ----------
        replacements : list of ((int, int), str)
            Each element is ((pos0, pos1), replacement) where pos0 and
            pos1 are indices to the original font data (parts[0] and the
            decrypted part concatenated). The data in the interval
            pos0:pos1 will be replaced by the replacement text. To
            accommodate binary data, the replacement is taken to be in
            Latin-1 encoding.

            The case where pos0 is inside parts[0] and pos1 inside
            the decrypted part is not supported.

        Returns
        -------
        (bytes, bytes)
            The new parts[0] and decrypted part (which needs to be
            encrypted in the transformed font).
        """
        data = bytearray(self.parts[0])
        data.extend(self.decrypted)
        len0 = len(self.parts[0])
        for (pos0, pos1), value in sorted(replacements, reverse=True):
            data[pos0:pos1] = value.encode('latin-1')
            if pos0 < len(self.parts[0]):
                if pos1 >= len(self.parts[0]):
                    raise RuntimeError(
                        f"text to be replaced with {value} spans "
                        "the eexec boundary"
                    )
                len0 += len(value) - pos1 + pos0

        return bytes(data[:len0]), bytes(data[len0:])

    def subset(self, characters, name_prefix):
        """
        Return a new font that only defines the given characters.

        Parameters
        ----------
        characters : sequence of bytes
            The subset of characters to include. These are indices into the
            font's encoding array. The encoding array of a Type-1 font can
            only include 256 characters, but other glyphs may be accessed
            via the seac operator.
        name_prefix : str
            Prefix to prepend to the font name.

        Returns
        -------
        `Type1Font`
        """
        characters = frozenset(characters)
        if _log.isEnabledFor(logging.DEBUG):
            _log.debug(
                "Subsetting font %s to characters %s = %s",
                self.prop['FontName'],
                sorted(characters),
                [self.prop['Encoding'].get(code) for code in sorted(characters)],
            )
        encoding = {code: glyph
                    for code, glyph in self.prop['Encoding'].items()
                    if code in characters}
        # todo and done include strings (glyph names)
        todo = {'.notdef', *encoding.values()}
        done = set()
        seen_subrs = {0, 1, 2, 3}
        while todo:
            glyph = todo.pop()
            called_glyphs, called_subrs = _CharstringSimulator(self).run(glyph)
            done.add(glyph)
            todo.update(called_glyphs - done)
            seen_subrs.update(called_subrs)

        charstrings = self._subset_charstrings(done)
        subrs = self._subset_subrs(seen_subrs)
        newparts = self._replace(
            [(x, f'/FontName /{name_prefix}{self.prop["FontName"]} def')
             for x in self._pos['FontName']]
            + [(self._pos['CharStrings'][0], charstrings),
               (self._pos['Subrs'][0], subrs),
               (self._pos['Encoding'][0], self._postscript_encoding(encoding))
               ] + [(x, '') for x in self._pos.get('UniqueID', [])]
        )
        return type(self)((
            newparts[0],
            self._encrypt(newparts[1], 'eexec'),
            self.parts[2]
        ))

    @staticmethod
    def _charstring_tokens(data):
        """Parse a Type-1 charstring

        Yield opcode names and integer parameters.
        """
        data = iter(data)
        for byte in data:
            if 32 <= byte <= 246:
                yield byte - 139
            elif 247 <= byte <= 250:
                byte2 = next(data)
                yield (byte-247) * 256 + byte2 + 108
            elif 251 <= byte <= 254:
                byte2 = next(data)
                yield -(byte-251)*256 - byte2 - 108
            elif byte == 255:
                bs = bytes(itertools.islice(data, 4))
                yield struct.unpack('>i', bs)[0]
            elif byte == 12:
                byte1 = next(data)
                yield {
                    0: 'dotsection',
                    1: 'vstem3',
                    2: 'hstem3',
                    6: 'seac',
                    7: 'sbw',
                    12: 'div',
                    16: 'callothersubr',
                    17: 'pop',
                    33: 'setcurrentpoint'
                }[byte1]
            else:
                yield {
                    1: 'hstem',
                    3: 'vstem',
                    4: 'vmoveto',
                    5: 'rlineto',
                    6: 'hlineto',
                    7: 'vlineto',
                    8: 'rrcurveto',
                    9: 'closepath',
                    10: 'callsubr',
                    11: 'return',
                    13: 'hsbw',
                    14: 'endchar',
                    21: 'rmoveto',
                    22: 'hmoveto',
                    30: 'vhcurveto',
                    31: 'hvcurveto'
                }[byte]

    def _postscript_encoding(self, encoding):
        """Return a PostScript encoding array for the encoding."""
        return '\n'.join([
            '/Encoding 256 array',
            '0 1 255 { 1 index exch /.notdef put} for',
            *(
                f'dup {i} /{glyph} put'
                for i, glyph in sorted(encoding.items())
                if glyph != '.notdef'
            ),
            'readonly def\n',
        ])

    def _subset_charstrings(self, glyphs):
        """Return a PostScript CharStrings array for the glyphs."""
        charstrings = self.prop['CharStrings']
        lenIV = self.prop.get('lenIV', 4)
        ordered = sorted(glyphs)
        encrypted = [
            self._encrypt(charstrings[glyph], 'charstring', lenIV).decode('latin-1')
            for glyph in ordered
        ]
        RD, ND = self._abbr['RD'], self._abbr['ND']
        return '\n'.join([
            f'/CharStrings {len(ordered)} dict dup begin',
            *(
                f'/{glyph} {len(enc)} {RD} {enc} {ND}'
                for glyph, enc in zip(ordered, encrypted)
            ),
            'end\n',
        ])

    def _subset_subrs(self, indices):
        """Return a PostScript Subrs array for the subroutines."""
        # we can't remove subroutines, we just replace unused ones with a stub
        subrs = self.prop['Subrs']
        n_subrs = len(subrs)
        lenIV = self.prop.get('lenIV', 4)
        stub = self._encrypt(b'\x0b', 'charstring', lenIV).decode('latin-1')
        encrypted = [
            self._encrypt(subrs[i], 'charstring', lenIV).decode('latin-1')
            if i in indices
            else stub
            for i in range(n_subrs)
        ]
        RD, ND, NP = self._abbr['RD'], self._abbr['ND'], self._abbr['NP']
        return '\n'.join([
            f'/Subrs {n_subrs} array',
            *(
                f'dup {i} {len(enc)} {RD} {enc} {NP}'
                for i, enc in enumerate(encrypted)
            ),
        ])


class _CharstringSimulator:
    __slots__ = ('font', 'buildchar_stack', 'postscript_stack', 'glyphs', 'subrs')

    def __init__(self, font):
        self.font = font
        self.buildchar_stack = []
        self.postscript_stack = []
        self.glyphs = set()
        self.subrs = set()

    def run(self, glyph_or_subr):
        """Run the charstring interpreter on a glyph or subroutine.

        This does not actually execute the code but simulates it to find out
        which subroutines get called when executing the glyph or subroutine.

        Parameters
        ----------
        glyph_or_subr : str or int
            The name of the glyph or the index of the subroutine to simulate.

        Returns
        -------
        glyphs : set[str]
            The set of glyph names called by the glyph or subroutine.
        subrs : set[int]
            The set of subroutines called by the glyph or subroutine.
        """
        if isinstance(glyph_or_subr, str):
            program = self.font.prop['CharStrings'][glyph_or_subr]
            self.glyphs.add(glyph_or_subr)
        else:
            program = self.font.prop['Subrs'][glyph_or_subr]
            self.subrs.add(glyph_or_subr)
        for opcode in self.font._charstring_tokens(program):
            if opcode in ('return', 'endchar'):
                return self.glyphs, self.subrs
            self._step(opcode)
        else:
            font_name = self.font.prop.get('FontName', '(unknown)')
            _log.info(
                f"Glyph or subr {glyph_or_subr} in font {font_name} does not end "
                "with return or endchar"
            )
            return self.glyphs, self.subrs

    def _step(self, opcode):
        """Run one step in the charstring interpreter."""
        match opcode:
            case int():
                self.buildchar_stack.append(opcode)
            case (
                'hsbw' | 'sbw' | 'closepath' | 'hlineto' | 'hmoveto' | 'hcurveto' |
                'hvcurveto' | 'rlineto' | 'rmoveto' | 'rrcurveto' | 'vhcurveto' |
                'vlineto' | 'vmoveto' | 'dotsection' | 'hstem' | 'hstem3' |
                'vstem' | 'vstem3' | 'setcurrentpoint'
            ):
                self.buildchar_stack.clear()
            case 'seac':  # Standard Encoding Accented Character
                codes = self.buildchar_stack[3:5]
                self.glyphs.update(_StandardEncoding[int(x)] for x in codes)
                self.buildchar_stack.clear()
            case 'div':
                num1, num2 = self.buildchar_stack[-2:]
                if num2 == 0:
                    _log.warning(
                        f"Division by zero in font {self.font.prop['FontName']}"
                    )
                    self.buildchar_stack[-2:] = [0]
                else:
                    self.buildchar_stack[-2:] = [num1/num2]
            case 'callothersubr':
                n, othersubr = self.buildchar_stack[-2:]
                if not isinstance(n, int):
                    _log.warning(
                        f"callothersubr {othersubr} with non-integer argument "
                        f"count in font {self.font.prop['FontName']}"
                    )
                    n = int(n)
                args = self.buildchar_stack[-2-n:-2]
                if othersubr == 3:
                    self.postscript_stack.append(args[0])
                else:
                    self.postscript_stack.extend(args[::-1])
                self.buildchar_stack[-2-n:] = []
            case 'callsubr':
                subr = self.buildchar_stack.pop()
                if not isinstance(subr, int):
                    _log.warning(
                        f"callsubr with non-integer argument {subr} in font "
                        f"{self.font.prop['FontName']}"
                    )
                    subr = int(subr)
                self.run(subr)
            case 'pop':
                if not self.postscript_stack:
                    _log.warning(
                        f"pop with empty stack in font {self.font.prop['FontName']}"
                    )
                    self.postscript_stack.append(0)
                self.buildchar_stack.append(self.postscript_stack.pop())
            case _:
                raise RuntimeError(f'opcode {opcode}')


_StandardEncoding = {
    **{ord(letter): letter for letter in string.ascii_letters},
    0: '.notdef',
    32: 'space',
    33: 'exclam',
    34: 'quotedbl',
    35: 'numbersign',
    36: 'dollar',
    37: 'percent',
    38: 'ampersand',
    39: 'quoteright',
    40: 'parenleft',
    41: 'parenright',
    42: 'asterisk',
    43: 'plus',
    44: 'comma',
    45: 'hyphen',
    46: 'period',
    47: 'slash',
    48: 'zero',
    49: 'one',
    50: 'two',
    51: 'three',
    52: 'four',
    53: 'five',
    54: 'six',
    55: 'seven',
    56: 'eight',
    57: 'nine',
    58: 'colon',
    59: 'semicolon',
    60: 'less',
    61: 'equal',
    62: 'greater',
    63: 'question',
    64: 'at',
    91: 'bracketleft',
    92: 'backslash',
    93: 'bracketright',
    94: 'asciicircum',
    95: 'underscore',
    96: 'quoteleft',
    123: 'braceleft',
    124: 'bar',
    125: 'braceright',
    126: 'asciitilde',
    161: 'exclamdown',
    162: 'cent',
    163: 'sterling',
    164: 'fraction',
    165: 'yen',
    166: 'florin',
    167: 'section',
    168: 'currency',
    169: 'quotesingle',
    170: 'quotedblleft',
    171: 'guillemotleft',
    172: 'guilsinglleft',
    173: 'guilsinglright',
    174: 'fi',
    175: 'fl',
    177: 'endash',
    178: 'dagger',
    179: 'daggerdbl',
    180: 'periodcentered',
    182: 'paragraph',
    183: 'bullet',
    184: 'quotesinglbase',
    185: 'quotedblbase',
    186: 'quotedblright',
    187: 'guillemotright',
    188: 'ellipsis',
    189: 'perthousand',
    191: 'questiondown',
    193: 'grave',
    194: 'acute',
    195: 'circumflex',
    196: 'tilde',
    197: 'macron',
    198: 'breve',
    199: 'dotaccent',
    200: 'dieresis',
    202: 'ring',
    203: 'cedilla',
    205: 'hungarumlaut',
    206: 'ogonek',
    207: 'caron',
    208: 'emdash',
    225: 'AE',
    227: 'ordfeminine',
    232: 'Lslash',
    233: 'Oslash',
    234: 'OE',
    235: 'ordmasculine',
    241: 'ae',
    245: 'dotlessi',
    248: 'lslash',
    249: 'oslash',
    250: 'oe',
    251: 'germandbls',
}
