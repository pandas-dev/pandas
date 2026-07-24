"""
A Python interface to Adobe Font Metrics Files.

Although a number of other Python implementations exist, and may be more
complete than this, it was decided not to go with them because they were
either:

1) copyrighted or used a non-BSD compatible license
2) had too many dependencies and a free standing lib was needed
3) did more than needed and it was easier to write afresh rather than
   figure out how to get just what was needed.

It is pretty easy to use, and has no external dependencies:

>>> import matplotlib as mpl
>>> from pathlib import Path
>>> afm_path = Path(mpl.get_data_path(), 'fonts', 'afm', 'ptmr8a.afm')
>>>
>>> from matplotlib._afm import AFM
>>> with afm_path.open('rb') as fh:
...     afm = AFM(fh)
>>> afm.get_fontname()
'Times-Roman'

As in the Adobe Font Metrics File Format Specification, all dimensions
are given in units of 1/1000 of the scale factor (point size) of the font
being used.
"""

import inspect
import logging
import re
from typing import BinaryIO, NamedTuple, TypedDict, cast

from ._mathtext_data import uni2type1
from .ft2font import CharacterCodeType, GlyphIndexType


_log = logging.getLogger(__name__)


def _to_int(x: bytes | str) -> int:
    # Some AFM files have floats where we are expecting ints -- there is
    # probably a better way to handle this (support floats, round rather than
    # truncate).  But I don't know what the best approach is now and this
    # change to _to_int should at least prevent Matplotlib from crashing on
    # these.  JDH (2009-11-06)
    return int(float(x))


def _to_float(x: bytes | str) -> float:
    # Some AFM files use "," instead of "." as decimal separator -- this
    # shouldn't be ambiguous (unless someone is wicked enough to use "," as
    # thousands separator...).
    if isinstance(x, bytes):
        # Encoding doesn't really matter -- if we have codepoints >127 the call
        # to float() will error anyways.
        x = x.decode('latin-1')
    return float(x.replace(',', '.'))


def _to_str(x: bytes) -> str:
    return x.decode('utf8')


def _to_list_of_ints(s: bytes) -> list[int]:
    s = s.replace(b',', b' ')
    return [_to_int(val) for val in s.split()]


def _to_list_of_floats(s: bytes | str) -> list[float]:
    return [_to_float(val) for val in s.split()]


def _to_bool(s: bytes) -> bool:
    if s.lower().strip() in (b'false', b'0', b'no'):
        return False
    else:
        return True


class FontMetricsHeader(TypedDict, total=False):
    StartFontMetrics: float
    FontName: str
    FullName: str
    FamilyName: str
    Weight: str
    ItalicAngle: float
    IsFixedPitch: bool
    FontBBox: list[int]
    UnderlinePosition: float
    UnderlineThickness: float
    Version: str
    # Some AFM files have non-ASCII characters (which are not allowed by the spec).
    # Given that there is actually no public API to even access this field, just return
    # it as straight bytes.
    Notice: bytes
    EncodingScheme: str
    CapHeight: float  # Is the second version a mistake, or
    Capheight: float  # do some AFM files contain 'Capheight'? -JKS
    XHeight: float
    Ascender: float
    Descender: float
    StdHW: float
    StdVW: float
    StartCharMetrics: int
    CharacterSet: str
    Characters: int


def _parse_header(fh: BinaryIO) -> FontMetricsHeader:
    """
    Read the font metrics header (up to the char metrics).

    Returns
    -------
    dict
        A dictionary mapping *key* to *val*. Dictionary keys are:

            StartFontMetrics, FontName, FullName, FamilyName, Weight, ItalicAngle,
            IsFixedPitch, FontBBox, UnderlinePosition, UnderlineThickness, Version,
            Notice, EncodingScheme, CapHeight, XHeight, Ascender, Descender,
            StartCharMetrics

        *val* will be converted to the appropriate Python type as necessary, e.g.,:

            * 'False' -> False
            * '0' -> 0
            * '-168 -218 1000 898' -> [-168, -218, 1000, 898]
    """
    header_converters = {
        bool: _to_bool,
        bytes: lambda x: x,
        float: _to_float,
        int: _to_int,
        list[int]: _to_list_of_ints,
        str: _to_str,
    }
    header_value_types = inspect.get_annotations(FontMetricsHeader)
    d: FontMetricsHeader = {}
    first_line = True
    for line in fh:
        line = line.rstrip()
        if line.startswith(b'Comment'):
            continue
        lst = line.split(b' ', 1)
        key = lst[0]
        if first_line:
            # AFM spec, Section 4: The StartFontMetrics keyword
            # [followed by a version number] must be the first line in
            # the file, and the EndFontMetrics keyword must be the
            # last non-empty line in the file.  We just check the
            # first header entry.
            if key != b'StartFontMetrics':
                raise RuntimeError('Not an AFM file')
            first_line = False
        if len(lst) == 2:
            val = lst[1]
        else:
            val = b''
        try:
            key_str = _to_str(key)
            value_type = header_value_types[key_str]
        except (KeyError, UnicodeDecodeError):
            _log.error("Found an unknown keyword in AFM header (was %r)", key)
            continue
        try:
            converter = header_converters[value_type]
            d[key_str] = converter(val)  # type: ignore[literal-required]
        except ValueError:
            _log.error('Value error parsing header in AFM: %r, %r', key, val)
            continue
        if key == b'StartCharMetrics':
            break
    else:
        raise RuntimeError('Bad parse')
    return d


class CharMetrics(NamedTuple):
    """
    Represents the character metrics of a single character.

    Notes
    -----
    The fields do currently only describe a subset of character metrics
    information defined in the AFM standard.
    """

    width: float
    name: str
    bbox: tuple[int, int, int, int]


CharMetrics.width.__doc__ = """The character width (WX)."""
CharMetrics.name.__doc__ = """The character name (N)."""
CharMetrics.bbox.__doc__ = """
    The bbox of the character (B) as a tuple (*llx*, *lly*, *urx*, *ury*)."""


def _parse_char_metrics(fh: BinaryIO) -> tuple[dict[CharacterCodeType, CharMetrics],
                                               dict[str, CharMetrics]]:
    """
    Parse the given filehandle for character metrics information.

    It is assumed that the file cursor is on the line behind 'StartCharMetrics'.

    Returns
    -------
    ascii_d : dict
         A mapping "ASCII num of the character" to `.CharMetrics`.
    name_d : dict
         A mapping "character name" to `.CharMetrics`.

    Notes
    -----
    This function is incomplete per the standard, but thus far parses
    all the sample afm files tried.
    """
    required_keys = {'C', 'WX', 'N', 'B'}

    ascii_d: dict[CharacterCodeType, CharMetrics] = {}
    name_d: dict[str, CharMetrics] = {}
    for bline in fh:
        # We are defensively letting values be utf8. The spec requires
        # ascii, but there are non-compliant fonts in circulation
        line = _to_str(bline.rstrip())
        if line.startswith('EndCharMetrics'):
            return ascii_d, name_d
        # Split the metric line into a dictionary, keyed by metric identifiers
        vals = dict(s.strip().split(' ', 1) for s in line.split(';') if s)
        # There may be other metrics present, but only these are needed
        if not required_keys.issubset(vals):
            raise RuntimeError('Bad char metrics line: %s' % line)
        num = _to_int(vals['C'])
        wx = _to_float(vals['WX'])
        name = vals['N']
        bbox = tuple(map(int, _to_list_of_floats(vals['B'])))
        if len(bbox) != 4:
            raise RuntimeError(f'Bad parse: bbox has {len(bbox)} elements, should be 4')
        metrics = CharMetrics(wx, name, bbox)
        # Workaround: If the character name is 'Euro', give it the
        # corresponding character code, according to WinAnsiEncoding (see PDF
        # Reference).
        if name == 'Euro':
            num = 128
        elif name == 'minus':
            num = ord("\N{MINUS SIGN}")  # 0x2212
        if num != -1:
            ascii_d[num] = metrics
        name_d[name] = metrics
    raise RuntimeError('Bad parse')


def _parse_kern_pairs(fh: BinaryIO) -> dict[tuple[str, str], float]:
    """
    Return a kern pairs dictionary.

    Returns
    -------
    dict
        Keys are (*char1*, *char2*) tuples and values are the kern pair value. For
        example, a kern pairs line like ``KPX A y -50`` will be represented as::

            d['A', 'y'] = -50
    """
    line = next(fh)
    if not line.startswith(b'StartKernPairs'):
        raise RuntimeError(f'Bad start of kern pairs data: {line!r}')

    d: dict[tuple[str, str], float] = {}
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(b'EndKernPairs'):
            next(fh)  # EndKernData
            return d
        vals = line.split()
        if len(vals) != 4 or vals[0] != b'KPX':
            raise RuntimeError(f'Bad kern pairs line: {line!r}')
        c1, c2, val = _to_str(vals[1]), _to_str(vals[2]), _to_float(vals[3])
        d[(c1, c2)] = val
    raise RuntimeError('Bad kern pairs parse')


class CompositePart(NamedTuple):
    """Represents the information on a composite element of a composite char."""

    name: bytes
    dx: float
    dy: float


CompositePart.name.__doc__ = """Name of the part, e.g. 'acute'."""
CompositePart.dx.__doc__ = """x-displacement of the part from the origin."""
CompositePart.dy.__doc__ = """y-displacement of the part from the origin."""


def _parse_composites(fh: BinaryIO) -> dict[bytes, list[CompositePart]]:
    """
    Parse the given filehandle for composites information.

    It is assumed that the file cursor is on the line behind 'StartComposites'.

    Returns
    -------
    dict
        A dict mapping composite character names to a parts list. The parts
        list is a list of `.CompositePart` entries describing the parts of
        the composite.

    Examples
    --------
    A composite definition line::

      CC Aacute 2 ; PCC A 0 0 ; PCC acute 160 170 ;

    will be represented as::

      composites[b'Aacute'] = [CompositePart(name=b'A', dx=0, dy=0),
                               CompositePart(name=b'acute', dx=160, dy=170)]

    """
    composites: dict[bytes, list[CompositePart]] = {}
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(b'EndComposites'):
            return composites
        vals = line.split(b';')
        cc = vals[0].split()
        name, _num_parts = cc[1], _to_int(cc[2])
        if len(vals) != _num_parts + 2:  # First element is 'CC', last is empty.
            raise RuntimeError(f'Bad composites parse: expected {_num_parts} parts, '
                               f'but got {len(vals) - 2}')
        pccParts = []
        for s in vals[1:-1]:
            pcc = s.split()
            part = CompositePart(pcc[1], _to_float(pcc[2]), _to_float(pcc[3]))
            pccParts.append(part)
        composites[name] = pccParts

    raise RuntimeError('Bad composites parse')


def _parse_optional(fh: BinaryIO) -> tuple[dict[tuple[str, str], float],
                                           dict[bytes, list[CompositePart]]]:
    """
    Parse the optional fields for kern pair data and composites.

    Returns
    -------
    kern_data : dict
        A dict containing kerning information. May be empty.
        See `._parse_kern_pairs`.
    composites : dict
        A dict containing composite information. May be empty.
        See `._parse_composites`.
    """
    kern_data: dict[tuple[str, str], float] = {}
    composites: dict[bytes, list[CompositePart]] = {}
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        match line.split()[0]:
            case b'StartKernData':
                kern_data = _parse_kern_pairs(fh)
            case b'StartComposites':
                composites = _parse_composites(fh)

    return kern_data, composites


class AFM:

    def __init__(self, fh: BinaryIO):
        """Parse the AFM file in file object *fh*."""
        self._header = _parse_header(fh)
        self._metrics, self._metrics_by_name = _parse_char_metrics(fh)
        self._kern, self._composite = _parse_optional(fh)

    def get_str_bbox_and_descent(self, s: str) -> tuple[int, int, float, int, int]:
        """Return the string bounding box and the maximal descent."""
        if not len(s):
            return 0, 0, 0, 0, 0
        total_width = 0.0
        namelast = ''
        miny = 1_000_000_000
        maxy = 0
        left = 0
        for c in s:
            if c == '\n':
                continue
            name = uni2type1.get(ord(c), f"uni{ord(c):04X}")
            try:
                wx, _, bbox = self._metrics_by_name[name]
            except KeyError:
                name = 'question'
                wx, _, bbox = self._metrics_by_name[name]
            total_width += wx + self._kern.get((namelast, name), 0)
            l, b, w, h = bbox
            left = min(left, l)
            miny = min(miny, b)
            maxy = max(maxy, b + h)

            namelast = name

        return left, miny, total_width, maxy - miny, -miny

    def get_glyph_name(self,  # For consistency with FT2Font.
                       glyph_ind: GlyphIndexType) -> str:
        """Get the name of the glyph, i.e., ord(';') is 'semicolon'."""
        return self._metrics[cast(CharacterCodeType, glyph_ind)].name

    def get_char_index(self,  # For consistency with FT2Font.
                       c: CharacterCodeType) -> GlyphIndexType:
        """
        Return the glyph index corresponding to a character code point.

        Note, for AFM fonts, we treat the glyph index the same as the codepoint.
        """
        return cast(GlyphIndexType, c)

    def get_width_char(self, c: CharacterCodeType) -> float:
        """Get the width of the character code from the character metric WX field."""
        return self._metrics[c].width

    def get_width_from_char_name(self, name: str) -> float:
        """Get the width of the character from a type1 character name."""
        return self._metrics_by_name[name].width

    def get_kern_dist_from_name(self, name1: str, name2: str) -> float:
        """
        Return the kerning pair distance (possibly 0) for chars *name1* and *name2*.
        """
        return self._kern.get((name1, name2), 0)

    def get_fontname(self) -> str:
        """Return the font name, e.g., 'Times-Roman'."""
        return self._header['FontName']

    @property
    def postscript_name(self) -> str:  # For consistency with FT2Font.
        return self.get_fontname()

    def get_fullname(self) -> str:
        """Return the font full name, e.g., 'Times-Roman'."""
        name = self._header.get('FullName')
        if name is None:  # use FontName as a substitute
            name = self._header['FontName']
        return name

    def get_familyname(self) -> str:
        """Return the font family name, e.g., 'Times'."""
        name = self._header.get('FamilyName')
        if name is not None:
            return name

        # FamilyName not specified so we'll make a guess
        name = self.get_fullname()
        extras = (r'(?i)([ -](regular|plain|italic|oblique|bold|semibold|'
                  r'light|ultralight|extra|condensed))+$')
        return re.sub(extras, '', name)

    @property
    def family_name(self) -> str:  # For consistency with FT2Font.
        """The font family name, e.g., 'Times'."""
        return self.get_familyname()

    def get_weight(self) -> str:
        """Return the font weight, e.g., 'Bold' or 'Roman'."""
        return self._header['Weight']

    def get_angle(self) -> float:
        """Return the fontangle as float."""
        return self._header['ItalicAngle']

    def get_ascender(self) -> float:
        """Return the ascent as float."""
        return self._header['Ascender']

    def get_capheight(self) -> float:
        """Return the cap height as float."""
        return self._header['CapHeight']

    def get_descender(self) -> float:
        """Return the descent as float."""
        return self._header['Descender']

    def get_xheight(self) -> float:
        """Return the xheight as float."""
        return self._header['XHeight']

    def get_underline_thickness(self) -> float:
        """Return the underline thickness as float."""
        return self._header['UnderlineThickness']
