"""
A PDF Matplotlib backend.

Author: Jouni K Seppänen <jks@iki.fi> and others.
"""

import codecs
from datetime import timezone
from datetime import datetime
from enum import Enum
from functools import total_ordering
from io import BytesIO
import itertools
import logging
import math
import os
import string
import struct
import sys
import time
import types
import zlib

import numpy as np
from PIL import Image

import matplotlib as mpl
from matplotlib import _api, _text_helpers, _type1font, cbook, dviread
from matplotlib._pylab_helpers import Gcf
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, GraphicsContextBase,
    RendererBase)
from matplotlib.backends.backend_mixed import MixedModeRenderer
from matplotlib.figure import Figure
from matplotlib.font_manager import FontPath, get_font, fontManager as _fontManager
from matplotlib._afm import AFM
from matplotlib.ft2font import FT2Font, FaceFlags, LoadFlags, StyleFlags
from matplotlib.transforms import Affine2D, BboxBase
from matplotlib.path import Path
from matplotlib.dates import UTC
from matplotlib import _path
from . import _backend_pdf_ps

_log = logging.getLogger(__name__)

# Overview
#
# The low-level knowledge about pdf syntax lies mainly in the pdfRepr
# function and the classes Reference, Name, Operator, and Stream.  The
# PdfFile class knows about the overall structure of pdf documents.
# It provides a "write" method for writing arbitrary strings in the
# file, and an "output" method that passes objects through the pdfRepr
# function before writing them in the file.  The output method is
# called by the RendererPdf class, which contains the various draw_foo
# methods.  RendererPdf contains a GraphicsContextPdf instance, and
# each draw_foo calls self.check_gc before outputting commands.  This
# method checks whether the pdf graphics state needs to be modified
# and outputs the necessary commands.  GraphicsContextPdf represents
# the graphics state, and its "delta" method returns the commands that
# modify the state.

# Add "pdf.use14corefonts: True" in your configuration file to use only
# the 14 PDF core fonts. These fonts do not need to be embedded; every
# PDF viewing application is required to have them. This results in very
# light PDF files you can use directly in LaTeX or ConTeXt documents
# generated with pdfTeX, without any conversion.

# These fonts are: Helvetica, Helvetica-Bold, Helvetica-Oblique,
# Helvetica-BoldOblique, Courier, Courier-Bold, Courier-Oblique,
# Courier-BoldOblique, Times-Roman, Times-Bold, Times-Italic,
# Times-BoldItalic, Symbol, ZapfDingbats.
#
# Some tricky points:
#
# 1. The clip path can only be widened by popping from the state
# stack.  Thus the state must be pushed onto the stack before narrowing
# the clip path.  This is taken care of by GraphicsContextPdf.
#
# 2. Sometimes it is necessary to refer to something (e.g., font,
# image, or extended graphics state, which contains the alpha value)
# in the page stream by a name that needs to be defined outside the
# stream.  PdfFile provides the methods fontName, imageObject, and
# alphaState for this purpose.  The implementations of these methods
# should perhaps be generalized.

# TODOs:
#
# * encoding of fonts, including mathtext fonts and Unicode support
# * TTF support has lots of small TODOs, e.g., how do you know if a font
#   is serif/sans-serif, or symbolic/non-symbolic?
# * draw_quad_mesh


def _fill(strings, linelen=75):
    """
    Make one string from sequence of strings, with whitespace in between.

    The whitespace is chosen to form lines of at most *linelen* characters,
    if possible.
    """
    currpos = 0
    lasti = 0
    result = []
    for i, s in enumerate(strings):
        length = len(s)
        if currpos + length < linelen:
            currpos += length + 1
        else:
            result.append(b' '.join(strings[lasti:i]))
            lasti = i
            currpos = length
    result.append(b' '.join(strings[lasti:]))
    return b'\n'.join(result)


def _create_pdf_info_dict(backend, metadata):
    """
    Create a PDF infoDict based on user-supplied metadata.

    A default ``Creator``, ``Producer``, and ``CreationDate`` are added, though
    the user metadata may override it. The date may be the current time, or a
    time set by the ``SOURCE_DATE_EPOCH`` environment variable.

    Metadata is verified to have the correct keys and their expected types. Any
    unknown keys/types will raise a warning.

    Parameters
    ----------
    backend : str
        The name of the backend to use in the Producer value.

    metadata : dict[str, Union[str, datetime, Name]]
        A dictionary of metadata supplied by the user with information
        following the PDF specification, also defined in
        `~.backend_pdf.PdfPages` below.

        If any value is *None*, then the key will be removed. This can be used
        to remove any pre-defined values.

    Returns
    -------
    dict[str, Union[str, datetime, Name]]
        A validated dictionary of metadata.
    """

    # get source date from SOURCE_DATE_EPOCH, if set
    # See https://reproducible-builds.org/specs/source-date-epoch/
    source_date_epoch = os.getenv("SOURCE_DATE_EPOCH")
    if source_date_epoch:
        source_date = datetime.fromtimestamp(int(source_date_epoch), timezone.utc)
        source_date = source_date.replace(tzinfo=UTC)
    else:
        source_date = datetime.today()

    info = {
        'Creator': f'Matplotlib v{mpl.__version__}, https://matplotlib.org',
        'Producer': f'Matplotlib {backend} backend v{mpl.__version__}',
        'CreationDate': source_date,
        **metadata
    }
    info = {k: v for (k, v) in info.items() if v is not None}

    def is_string_like(x):
        return isinstance(x, str)
    is_string_like.text_for_warning = "an instance of str"

    def is_date(x):
        return isinstance(x, datetime)
    is_date.text_for_warning = "an instance of datetime.datetime"

    def check_trapped(x):
        if isinstance(x, Name):
            return x.name in (b'True', b'False', b'Unknown')
        else:
            return x in ('True', 'False', 'Unknown')
    check_trapped.text_for_warning = 'one of {"True", "False", "Unknown"}'

    keywords = {
        'Title': is_string_like,
        'Author': is_string_like,
        'Subject': is_string_like,
        'Keywords': is_string_like,
        'Creator': is_string_like,
        'Producer': is_string_like,
        'CreationDate': is_date,
        'ModDate': is_date,
        'Trapped': check_trapped,
    }
    for k in info:
        if k not in keywords:
            _api.warn_external(f'Unknown infodict keyword: {k!r}. '
                               f'Must be one of {set(keywords)!r}.')
        elif not keywords[k](info[k]):
            _api.warn_external(f'Bad value for infodict keyword {k}. '
                               f'Got {info[k]!r} which is not '
                               f'{keywords[k].text_for_warning}.')
    if 'Trapped' in info:
        info['Trapped'] = Name(info['Trapped'])

    return info


def _datetime_to_pdf(d):
    """
    Convert a datetime to a PDF string representing it.

    Used for PDF and PGF.
    """
    r = d.strftime('D:%Y%m%d%H%M%S')
    z = d.utcoffset()
    if z is not None:
        z = z.seconds
    else:
        if time.daylight:
            z = time.altzone
        else:
            z = time.timezone
    if z == 0:
        r += 'Z'
    elif z < 0:
        r += "+%02d'%02d'" % ((-z) // 3600, (-z) % 3600)
    else:
        r += "-%02d'%02d'" % (z // 3600, z % 3600)
    return r


def _calculate_quad_point_coordinates(x, y, width, height, angle=0):
    """
    Calculate the coordinates of rectangle when rotated by angle around x, y
    """

    angle = math.radians(-angle)
    sin_angle = math.sin(angle)
    cos_angle = math.cos(angle)
    a = x + height * sin_angle
    b = y + height * cos_angle
    c = x + width * cos_angle + height * sin_angle
    d = y - width * sin_angle + height * cos_angle
    e = x + width * cos_angle
    f = y - width * sin_angle
    return ((x, y), (e, f), (c, d), (a, b))


def _get_coordinates_of_block(x, y, width, height, angle=0):
    """
    Get the coordinates of rotated rectangle and rectangle that covers the
    rotated rectangle.
    """

    vertices = _calculate_quad_point_coordinates(x, y, width,
                                                 height, angle)

    # Find min and max values for rectangle
    # adjust so that QuadPoints is inside Rect
    # PDF docs says that QuadPoints should be ignored if any point lies
    # outside Rect, but for Acrobat it is enough that QuadPoints is on the
    # border of Rect.

    pad = 0.00001 if angle % 90 else 0
    min_x = min(v[0] for v in vertices) - pad
    min_y = min(v[1] for v in vertices) - pad
    max_x = max(v[0] for v in vertices) + pad
    max_y = max(v[1] for v in vertices) + pad
    return (tuple(itertools.chain.from_iterable(vertices)),
            (min_x, min_y, max_x, max_y))


def _get_link_annotation(gc, x, y, width, height, angle=0):
    """
    Create a link annotation object for embedding URLs.
    """
    quadpoints, rect = _get_coordinates_of_block(x, y, width, height, angle)
    link_annotation = {
        'Type': Name('Annot'),
        'Subtype': Name('Link'),
        'Rect': rect,
        'Border': [0, 0, 0],
        'A': {
            'S': Name('URI'),
            'URI': gc.get_url(),
        },
    }
    if angle % 90:
        # Add QuadPoints
        link_annotation['QuadPoints'] = quadpoints
    return link_annotation


# PDF strings are supposed to be able to include any eight-bit data, except
# that unbalanced parens and backslashes must be escaped by a backslash.
# However, sf bug #2708559 shows that the carriage return character may get
# read as a newline; these characters correspond to \gamma and \Omega in TeX's
# math font encoding. Escaping them fixes the bug.
_str_escapes = str.maketrans({
    '\\': '\\\\', '(': '\\(', ')': '\\)', '\n': '\\n', '\r': '\\r'})


def pdfRepr(obj):
    """Map Python objects to PDF syntax."""

    # Some objects defined later have their own pdfRepr method.
    if hasattr(obj, 'pdfRepr'):
        return obj.pdfRepr()

    # Floats. PDF does not have exponential notation (1.0e-10) so we
    # need to use %f with some precision.  Perhaps the precision
    # should adapt to the magnitude of the number?
    elif isinstance(obj, (float, np.floating)):
        if not np.isfinite(obj):
            raise ValueError("Can only output finite numbers in PDF")
        r = b"%.10f" % obj
        return r.rstrip(b'0').rstrip(b'.')

    # Booleans. Needs to be tested before integers since
    # isinstance(True, int) is true.
    elif isinstance(obj, bool):
        return [b'false', b'true'][obj]

    # Integers are written as such.
    elif isinstance(obj, (int, np.integer)):
        return b"%d" % obj

    # Non-ASCII Unicode strings are encoded in UTF-16BE with byte-order mark.
    elif isinstance(obj, str):
        return pdfRepr(obj.encode('ascii') if obj.isascii()
                       else codecs.BOM_UTF16_BE + obj.encode('UTF-16BE'))

    # Strings are written in parentheses, with backslashes and parens
    # escaped. Actually balanced parens are allowed, but it is
    # simpler to escape them all. TODO: cut long strings into lines;
    # I believe there is some maximum line length in PDF.
    # Despite the extra decode/encode, translate is faster than regex.
    elif isinstance(obj, bytes):
        return (
            b'(' +
            obj.decode('latin-1').translate(_str_escapes).encode('latin-1')
            + b')')

    # Dictionaries. The keys must be PDF names, so if we find strings
    # there, we make Name objects from them. The values may be
    # anything, so the caller must ensure that PDF names are
    # represented as Name objects.
    elif isinstance(obj, dict):
        return _fill([
            b"<<",
            *[Name(k).pdfRepr() + b" " + pdfRepr(v) for k, v in obj.items()],
            b">>",
        ])

    # Lists.
    elif isinstance(obj, (list, tuple)):
        return _fill([b"[", *[pdfRepr(val) for val in obj], b"]"])

    # The null keyword.
    elif obj is None:
        return b'null'

    # A date.
    elif isinstance(obj, datetime):
        return pdfRepr(_datetime_to_pdf(obj))

    # A bounding box
    elif isinstance(obj, BboxBase):
        return _fill([pdfRepr(val) for val in obj.bounds])

    else:
        raise TypeError(f"Don't know a PDF representation for {type(obj)} "
                        "objects")


class Reference:
    """
    PDF reference object.

    Use PdfFile.reserveObject() to create References.
    """

    def __init__(self, id):
        self.id = id

    def __repr__(self):
        return "<Reference %d>" % self.id

    def pdfRepr(self):
        return b"%d 0 R" % self.id

    def write(self, contents, file):
        write = file.write
        write(b"%d 0 obj\n" % self.id)
        write(pdfRepr(contents))
        write(b"\nendobj\n")


@total_ordering
class Name:
    """PDF name object."""
    __slots__ = ('name',)
    _hexify = {c: '#%02x' % c
               for c in {*range(256)} - {*range(ord('!'), ord('~') + 1)}}

    def __init__(self, name):
        if isinstance(name, Name):
            self.name = name.name
        else:
            if isinstance(name, bytes):
                name = name.decode('ascii')
            self.name = name.translate(self._hexify).encode('ascii')

    def __repr__(self):
        return "<Name %s>" % self.name

    def __str__(self):
        return '/' + self.name.decode('ascii')

    def __eq__(self, other):
        return isinstance(other, Name) and self.name == other.name

    def __lt__(self, other):
        return isinstance(other, Name) and self.name < other.name

    def __hash__(self):
        return hash(self.name)

    def pdfRepr(self):
        return b'/' + self.name


class Verbatim:
    """Store verbatim PDF command content for later inclusion in the stream."""
    def __init__(self, x):
        self._x = x

    def pdfRepr(self):
        return self._x


class Op(Enum):
    """PDF operators (not an exhaustive list)."""

    close_fill_stroke = b'b'
    fill_stroke = b'B'
    fill = b'f'
    closepath = b'h'
    close_stroke = b's'
    stroke = b'S'
    endpath = b'n'
    begin_text = b'BT'
    end_text = b'ET'
    curveto = b'c'
    rectangle = b're'
    lineto = b'l'
    moveto = b'm'
    concat_matrix = b'cm'
    use_xobject = b'Do'
    setgray_stroke = b'G'
    setgray_nonstroke = b'g'
    setrgb_stroke = b'RG'
    setrgb_nonstroke = b'rg'
    setcolorspace_stroke = b'CS'
    setcolorspace_nonstroke = b'cs'
    setcolor_stroke = b'SCN'
    setcolor_nonstroke = b'scn'
    setdash = b'd'
    setlinejoin = b'j'
    setlinecap = b'J'
    setgstate = b'gs'
    gsave = b'q'
    grestore = b'Q'
    textpos = b'Td'
    selectfont = b'Tf'
    textmatrix = b'Tm'
    textrise = b'Ts'
    show = b'Tj'
    showkern = b'TJ'
    setlinewidth = b'w'
    clip = b'W'
    shading = b'sh'

    def pdfRepr(self):
        return self.value

    @classmethod
    def paint_path(cls, fill, stroke):
        """
        Return the PDF operator to paint a path.

        Parameters
        ----------
        fill : bool
            Fill the path with the fill color.
        stroke : bool
            Stroke the outline of the path with the line color.
        """
        if stroke:
            if fill:
                return cls.fill_stroke
            else:
                return cls.stroke
        else:
            if fill:
                return cls.fill
            else:
                return cls.endpath


class Stream:
    """
    PDF stream object.

    This has no pdfRepr method. Instead, call begin(), then output the
    contents of the stream by calling write(), and finally call end().
    """
    __slots__ = ('id', 'len', 'pdfFile', 'file', 'compressobj', 'extra', 'pos')

    def __init__(self, id, len, file, extra=None, png=None):
        """
        Parameters
        ----------
        id : int
            Object id of the stream.
        len : Reference or None
            An unused Reference object for the length of the stream;
            None means to use a memory buffer so the length can be inlined.
        file : PdfFile
            The underlying object to write the stream to.
        extra : dict from Name to anything, or None
            Extra key-value pairs to include in the stream header.
        png : dict or None
            If the data is already png encoded, the decode parameters.
        """
        self.id = id            # object id
        self.len = len          # id of length object
        self.pdfFile = file
        self.file = file.fh      # file to which the stream is written
        self.compressobj = None  # compression object
        if extra is None:
            self.extra = dict()
        else:
            self.extra = extra.copy()
        if png is not None:
            self.extra.update({'Filter':      Name('FlateDecode'),
                               'DecodeParms': png})

        self.pdfFile.recordXref(self.id)
        if mpl.rcParams['pdf.compression'] and not png:
            self.compressobj = zlib.compressobj(
                mpl.rcParams['pdf.compression'])
        if self.len is None:
            self.file = BytesIO()
        else:
            self._writeHeader()
            self.pos = self.file.tell()

    def _writeHeader(self):
        write = self.file.write
        write(b"%d 0 obj\n" % self.id)
        dict = self.extra
        dict['Length'] = self.len
        if mpl.rcParams['pdf.compression']:
            dict['Filter'] = Name('FlateDecode')

        write(pdfRepr(dict))
        write(b"\nstream\n")

    def end(self):
        """Finalize stream."""

        self._flush()
        if self.len is None:
            contents = self.file.getvalue()
            self.len = len(contents)
            self.file = self.pdfFile.fh
            self._writeHeader()
            self.file.write(contents)
            self.file.write(b"\nendstream\nendobj\n")
        else:
            length = self.file.tell() - self.pos
            self.file.write(b"\nendstream\nendobj\n")
            self.pdfFile.writeObject(self.len, length)

    def write(self, data):
        """Write some data on the stream."""

        if self.compressobj is None:
            self.file.write(data)
        else:
            compressed = self.compressobj.compress(data)
            self.file.write(compressed)

    def _flush(self):
        """Flush the compression object."""

        if self.compressobj is not None:
            compressed = self.compressobj.flush()
            self.file.write(compressed)
            self.compressobj = None


def _get_pdf_charprocs(font_path, glyph_indices):
    font = get_font(font_path)
    conv = 1000 / font.units_per_EM  # Conversion to PS units (1/1000's).
    procs = {}
    for glyph_index in glyph_indices:
        g = font.load_glyph(glyph_index, LoadFlags.NO_SCALE)
        d1 = [
            round(g.horiAdvance * conv), 0,
            # Round bbox corners *outwards*, so that they indeed bound the glyph.
            math.floor(g.bbox[0] * conv), math.floor(g.bbox[1] * conv),
            math.ceil(g.bbox[2] * conv), math.ceil(g.bbox[3] * conv),
        ]
        v, c = font.get_path()
        v = (v * 64 * conv).round()  # Back to TrueType's internal units (1/64's).
        procs[font.get_glyph_name(glyph_index)] = (
            " ".join(map(str, d1)).encode("ascii") + b" d1\n"
            + _path.convert_to_string(
                Path(v, c), None, None, False, None, 0,
                # no code for quad Beziers triggers auto-conversion to cubics.
                [b"m", b"l", b"", b"c", b"h"], True)
            + b"f")
    return procs


class PdfFile:
    """PDF file object."""

    def __init__(self, filename, metadata=None):
        """
        Parameters
        ----------
        filename : str or path-like or file-like
            Output target; if a string, a file will be opened for writing.

        metadata : dict from strings to strings and dates
            Information dictionary object (see PDF reference section 10.2.1
            'Document Information Dictionary'), e.g.:
            ``{'Creator': 'My software', 'Author': 'Me', 'Title': 'Awesome'}``.

            The standard keys are 'Title', 'Author', 'Subject', 'Keywords',
            'Creator', 'Producer', 'CreationDate', 'ModDate', and
            'Trapped'. Values have been predefined for 'Creator', 'Producer'
            and 'CreationDate'. They can be removed by setting them to `None`.
        """
        super().__init__()

        self._object_seq = itertools.count(1)  # consumed by reserveObject
        self.xrefTable = [[0, 65535, 'the zero object']]
        self.passed_in_file_object = False
        self.original_file_like = None
        self.tell_base = 0
        fh, opened = cbook.to_filehandle(filename, "wb", return_opened=True)
        if not opened:
            try:
                self.tell_base = filename.tell()
            except OSError:
                fh = BytesIO()
                self.original_file_like = filename
            else:
                fh = filename
                self.passed_in_file_object = True

        self.fh = fh
        self.currentstream = None  # stream object to write to, if any
        fh.write(b"%PDF-1.4\n")    # 1.4 is the first version to have alpha
        # Output some eight-bit chars as a comment so various utilities
        # recognize the file as binary by looking at the first few
        # lines (see note in section 3.4.1 of the PDF reference).
        fh.write(b"%\254\334 \253\272\n")

        self.rootObject = self.reserveObject('root')
        self.pagesObject = self.reserveObject('pages')
        self.pageList = []
        self.fontObject = self.reserveObject('fonts')
        self._extGStateObject = self.reserveObject('extended graphics states')
        self.hatchObject = self.reserveObject('tiling patterns')
        self.gouraudObject = self.reserveObject('Gouraud triangles')
        self.XObjectObject = self.reserveObject('external objects')
        self.resourceObject = self.reserveObject('resources')

        root = {'Type': Name('Catalog'),
                'Pages': self.pagesObject}
        self.writeObject(self.rootObject, root)

        self.infoDict = _create_pdf_info_dict('pdf', metadata or {})

        self._internal_font_seq = (Name(f'F{i}') for i in itertools.count(1))
        self._fontNames = {}     # maps filenames to internal font names
        self._dviFontInfo = {}   # maps pdf names to dvifonts
        self._character_tracker = _backend_pdf_ps.CharacterTracker(
            _backend_pdf_ps._FONT_MAX_GLYPH.get(mpl.rcParams['ps.fonttype'], 0))

        self.alphaStates = {}   # maps alpha values to graphics state objects
        self._alpha_state_seq = (Name(f'A{i}') for i in itertools.count(1))
        self._soft_mask_states = {}
        self._soft_mask_seq = (Name(f'SM{i}') for i in itertools.count(1))
        self._soft_mask_groups = []
        self._hatch_patterns = {}
        self._hatch_pattern_seq = (Name(f'H{i}') for i in itertools.count(1))
        self.gouraudTriangles = []

        self._images = {}
        self._image_seq = (Name(f'I{i}') for i in itertools.count(1))

        self.markers = {}

        self.paths = []

        # A list of annotations for each page. Each entry is a tuple of the
        # overall Annots object reference that's inserted into the page object,
        # followed by a list of the actual annotations.
        self._annotations = []
        # For annotations added before a page is created; mostly for the
        # purpose of newTextnote.
        self.pageAnnotations = []

        # The PDF spec recommends to include every procset
        procsets = [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]

        # Write resource dictionary.
        # Possibly TODO: more general ExtGState (graphics state dictionaries)
        #                ColorSpace Pattern Shading Properties
        resources = {'Font': self.fontObject,
                     'XObject': self.XObjectObject,
                     'ExtGState': self._extGStateObject,
                     'Pattern': self.hatchObject,
                     'Shading': self.gouraudObject,
                     'ProcSet': procsets}
        self.writeObject(self.resourceObject, resources)

    fontNames = _api.deprecated("3.11")(property(lambda self: self._fontNames))
    multi_byte_charprocs = _api.deprecated("3.11")(property(lambda _: {}))
    type1Descriptors = _api.deprecated("3.11")(property(lambda _: {}))

    @_api.deprecated("3.11")
    @property
    def dviFontInfo(self):
        d = {}
        tex_font_map = dviread.PsfontsMap(dviread.find_tex_file('pdftex.map'))
        for pdfname, dvifont in self._dviFontInfo.items():
            psfont = tex_font_map[dvifont.texname]
            if psfont.filename is None:
                raise ValueError(
                    "No usable font file found for {} (TeX: {}); "
                    "the font may lack a Type-1 version"
                    .format(psfont.psname, dvifont.texname))
            d[dvifont.texname] = types.SimpleNamespace(
                dvifont=dvifont,
                pdfname=pdfname,
                fontfile=psfont.filename,
                basefont=psfont.psname,
                encodingfile=psfont.encoding,
                effects=psfont.effects,
            )
        return d

    def newPage(self, width, height):
        self.endStream()

        self.width, self.height = width, height
        contentObject = self.reserveObject('page contents')
        annotsObject = self.reserveObject('annotations')
        thePage = {'Type': Name('Page'),
                   'Parent': self.pagesObject,
                   'Resources': self.resourceObject,
                   'MediaBox': [0, 0, 72 * width, 72 * height],
                   'Contents': contentObject,
                   'Annots': annotsObject,
                   }
        pageObject = self.reserveObject('page')
        self.writeObject(pageObject, thePage)
        self.pageList.append(pageObject)
        self._annotations.append((annotsObject, self.pageAnnotations))

        self.beginStream(contentObject.id,
                         self.reserveObject('length of content stream'))
        # Initialize the pdf graphics state to match the default Matplotlib
        # graphics context (colorspace and joinstyle).
        self.output(Name('DeviceRGB'), Op.setcolorspace_stroke)
        self.output(Name('DeviceRGB'), Op.setcolorspace_nonstroke)
        self.output(GraphicsContextPdf.joinstyles['round'], Op.setlinejoin)

        # Clear the list of annotations for the next page
        self.pageAnnotations = []

    def newTextnote(self, text, positionRect=[-100, -100, 0, 0]):
        # Create a new annotation of type text
        theNote = {'Type': Name('Annot'),
                   'Subtype': Name('Text'),
                   'Contents': text,
                   'Rect': positionRect,
                   }
        self.pageAnnotations.append(theNote)

    @staticmethod
    def _get_subset_prefix(charset):
        """
        Get a prefix for a subsetted font name.

        The prefix is six uppercase letters followed by a plus sign;
        see PDF reference section 5.5.3 Font Subsets.
        """
        def toStr(n, base):
            if n < base:
                return string.ascii_uppercase[n]
            else:
                return (
                    toStr(n // base, base) + string.ascii_uppercase[n % base]
                )

        # encode to string using base 26
        hashed = hash(charset) % ((sys.maxsize + 1) * 2)
        prefix = toStr(hashed, 26)

        # get first 6 characters from prefix
        return prefix[:6] + "+"

    @staticmethod
    def _get_subsetted_psname(ps_name, charmap):
        return PdfFile._get_subset_prefix(frozenset(charmap.values())) + ps_name

    def finalize(self):
        """Write out the various deferred objects and the pdf end matter."""

        self.endStream()
        self._write_annotations()
        self.writeFonts()
        self.writeExtGSTates()
        self._write_soft_mask_groups()
        self.writeHatches()
        self.writeGouraudTriangles()
        xobjects = {
            name: ob for image, name, ob in self._images.values()}
        for tup in self.markers.values():
            xobjects[tup[0]] = tup[1]
        for name, path, trans, ob, join, cap, padding, filled, stroked \
                in self.paths:
            xobjects[name] = ob
        self.writeObject(self.XObjectObject, xobjects)
        self.writeImages()
        self.writeMarkers()
        self.writePathCollectionTemplates()
        self.writeObject(self.pagesObject,
                         {'Type': Name('Pages'),
                          'Kids': self.pageList,
                          'Count': len(self.pageList)})
        self.writeInfoDict()

        # Finalize the file
        self.writeXref()
        self.writeTrailer()

    def close(self):
        """Flush all buffers and free all resources."""

        self.endStream()
        if self.passed_in_file_object:
            self.fh.flush()
        else:
            if self.original_file_like is not None:
                self.original_file_like.write(self.fh.getvalue())
            self.fh.close()

    def write(self, data):
        if self.currentstream is None:
            self.fh.write(data)
        else:
            self.currentstream.write(data)

    def output(self, *data):
        self.write(_fill([pdfRepr(x) for x in data]))
        self.write(b'\n')

    def beginStream(self, id, len, extra=None, png=None):
        assert self.currentstream is None
        self.currentstream = Stream(id, len, self, extra, png)

    def endStream(self):
        if self.currentstream is not None:
            self.currentstream.end()
            self.currentstream = None

    def outputStream(self, ref, data, *, extra=None):
        self.beginStream(ref.id, None, extra)
        self.currentstream.write(data)
        self.endStream()

    def _write_annotations(self):
        for annotsObject, annotations in self._annotations:
            self.writeObject(annotsObject, annotations)

    def fontName(self, fontprop, subset=0):
        """
        Select a font based on fontprop and return a name suitable for
        ``Op.selectfont``. If fontprop is a string, it will be interpreted
        as the filename of the font.
        """

        if isinstance(fontprop, FontPath):
            filenames = [fontprop]
        elif isinstance(fontprop, str):
            filenames = [FontPath(fontprop, 0)]
        elif mpl.rcParams['pdf.use14corefonts']:
            filenames = _fontManager._find_fonts_by_props(
                fontprop, fontext='afm', directory=RendererPdf._afm_font_dir
            )
        else:
            filenames = _fontManager._find_fonts_by_props(fontprop)
        first_Fx = None
        for fname in filenames:
            Fx = self._fontNames.get((fname, subset))
            if not first_Fx:
                first_Fx = Fx
            if Fx is None:
                Fx = next(self._internal_font_seq)
                self._fontNames[(fname, subset)] = Fx
                _log.debug('Assigning font %s (subset %d) = %r', Fx, subset, fname)
                if not first_Fx:
                    first_Fx = Fx

        # find_fontsprop's first value always adheres to
        # findfont's value, so technically no behaviour change
        return first_Fx

    def dviFontName(self, dvifont):
        """
        Given a dvi font object, return a name suitable for Op.selectfont.

        Register the font internally (in ``_dviFontInfo``) if not yet registered.
        """
        pdfname = Name(f"F-{dvifont.texname.decode('ascii')}")
        _log.debug('Assigning font %s = %s (dvi)', pdfname, dvifont.texname)
        self._dviFontInfo[pdfname] = dvifont
        return Name(pdfname)

    def writeFonts(self):
        fonts = {}
        for pdfname, dvifont in sorted(self._dviFontInfo.items()):
            _log.debug('Embedding Type-1 font %s from dvi.', dvifont.texname)
            fonts[pdfname] = self._embedTeXFont(dvifont)
        for (filename, subset), Fx in sorted(self._fontNames.items()):
            _log.debug('Embedding font %r:%d.', filename, subset)
            if filename.endswith('.afm'):
                # from pdf.use14corefonts
                _log.debug('Writing AFM font.')
                fonts[Fx] = self._write_afm_font(filename)
            else:
                # a normal TrueType font
                _log.debug('Writing TrueType font.')
                charmap = self._character_tracker.used[filename][subset]
                fonts[Fx] = self.embedTTF(filename, subset, charmap)
        self.writeObject(self.fontObject, fonts)

    def _write_afm_font(self, filename):
        with open(filename, 'rb') as fh:
            font = AFM(fh)
        fontname = font.get_fontname()
        fontdict = {'Type': Name('Font'),
                    'Subtype': Name('Type1'),
                    'BaseFont': Name(fontname),
                    'Encoding': Name('WinAnsiEncoding')}
        fontdictObject = self.reserveObject('font dictionary')
        self.writeObject(fontdictObject, fontdict)
        return fontdictObject

    def _embedTeXFont(self, dvifont):
        tex_font_map = dviread.PsfontsMap(dviread.find_tex_file('pdftex.map'))
        psfont = tex_font_map[dvifont.texname]
        if psfont.filename is None:
            raise ValueError(
                "No usable font file found for {} (TeX: {}); "
                "the font may lack a Type-1 version"
                .format(psfont.psname, dvifont.texname))

        # The font dictionary is the top-level object describing a font
        fontdictObject = self.reserveObject('font dictionary')
        fontdict = {
            'Type':      Name('Font'),
            'Subtype':   Name('Type1'),
        }

        # Read the font file and apply any encoding changes and effects
        t1font = _type1font.Type1Font(psfont.filename)
        if psfont.encoding is not None:
            t1font = t1font.with_encoding(
                {i: c for i, c in enumerate(dviread._parse_enc(psfont.encoding))}
            )
        if psfont.effects:
            t1font = t1font.transform(psfont.effects)

        # Reduce the font to only the glyphs used in the document, get the encoding
        # for that subset, and compute various properties based on the encoding.
        font_path = FontPath(dvifont.fname, dvifont.face_index)
        charmap = self._character_tracker.used[font_path][0]
        chars = {
            # DVI type 1 fonts always map single glyph to single character.
            ord(self._character_tracker.subset_to_unicode(font_path, 0, ccode))
            for ccode in charmap
        }
        t1font = t1font.subset(chars, self._get_subset_prefix(charmap.values()))
        fontdict['BaseFont'] = Name(t1font.prop['FontName'])
        # createType1Descriptor writes the font data as a side effect
        fontdict['FontDescriptor'] = self.createType1Descriptor(t1font)
        encoding = t1font.prop['Encoding']
        fontdict['Encoding'] = self._generate_encoding(encoding)
        fc = fontdict['FirstChar'] = min(encoding.keys(), default=0)
        lc = fontdict['LastChar'] = max(encoding.keys(), default=255)
        # Convert glyph widths from TeX 12.20 fixed point to 1/1000 text space units
        font_metrics = dvifont._metrics
        widths = [(1000 * glyph_metrics.tex_width) >> 20
                  if (glyph_metrics := font_metrics.get_metrics(char)) else 0
                  for char in range(fc, lc + 1)]
        fontdict['Widths'] = widthsObject = self.reserveObject('glyph widths')
        self.writeObject(widthsObject, widths)
        self.writeObject(fontdictObject, fontdict)
        return fontdictObject

    def _generate_encoding(self, encoding):
        prev = -2
        result = []
        for code, name in sorted(encoding.items()):
            if code != prev + 1:
                result.append(code)
            prev = code
            result.append(Name(name))
        return {
            'Type': Name('Encoding'),
            'Differences': result
        }

    @_api.delete_parameter("3.11", "fontfile")
    def createType1Descriptor(self, t1font, fontfile=None):
        # Create and write the font descriptor and the font file
        # of a Type-1 font
        fontdescObject = self.reserveObject('font descriptor')
        fontfileObject = self.reserveObject('font file')

        italic_angle = t1font.prop['ItalicAngle']
        fixed_pitch = t1font.prop['isFixedPitch']

        flags = 0
        # fixed width
        if fixed_pitch:
            flags |= 1 << 0
        # TODO: serif
        if 0:
            flags |= 1 << 1
        # TODO: symbolic (most TeX fonts are)
        if 1:
            flags |= 1 << 2
        # non-symbolic
        else:
            flags |= 1 << 5
        # italic
        if italic_angle:
            flags |= 1 << 6
        # TODO: all caps
        if 0:
            flags |= 1 << 16
        # TODO: small caps
        if 0:
            flags |= 1 << 17
        # TODO: force bold
        if 0:
            flags |= 1 << 18

        encoding = t1font.prop['Encoding']
        charset = ''.join(
            sorted(
                f'/{c}' for c in encoding.values()
                if c != '.notdef'
            )
        )

        descriptor = {
            'Type':        Name('FontDescriptor'),
            'FontName':    Name(t1font.prop['FontName']),
            'Flags':       flags,
            'FontBBox':    t1font.prop['FontBBox'],
            'ItalicAngle': italic_angle,
            'Ascent':      t1font.prop['FontBBox'][3],
            'Descent':     t1font.prop['FontBBox'][1],
            'CapHeight':   1000,  # TODO: find this out
            'XHeight':     500,  # TODO: this one too
            'FontFile':    fontfileObject,
            'FontFamily':  t1font.prop['FamilyName'],
            'StemV':       50,  # TODO
            'CharSet':     charset,
            # (see also revision 3874; but not all TeX distros have AFM files!)
            # 'FontWeight': a number where 400 = Regular, 700 = Bold
        }

        self.writeObject(fontdescObject, descriptor)

        self.outputStream(fontfileObject, b"".join(t1font.parts[:2]),
                          extra={'Length1': len(t1font.parts[0]),
                                 'Length2': len(t1font.parts[1]),
                                 'Length3': 0})

        return fontdescObject

    _identityToUnicodeCMap = b"""/CIDInit /ProcSet findresource begin
12 dict begin
begincmap
/CIDSystemInfo
<< /Registry (Adobe)
   /Ordering (UCS)
   /Supplement 0
>> def
/CMapName /Adobe-Identity-UCS def
/CMapType 2 def
1 begincodespacerange
<0000> <ffff>
endcodespacerange
%d beginbfrange
%s
endbfrange
endcmap
CMapName currentdict /CMap defineresource pop
end
end"""

    def embedTTF(self, filename, subset_index, charmap):
        """Embed the TTF font from the named file into the document."""
        font = get_font(filename)
        fonttype = mpl.rcParams['pdf.fonttype']

        def cvt(length, upe=font.units_per_EM, nearest=True):
            """Convert font coordinates to PDF glyph coordinates."""
            value = length / upe * 1000
            if nearest:
                return round(value)
            # Best(?) to round away from zero for bounding boxes and the like.
            if value < 0:
                return math.floor(value)
            else:
                return math.ceil(value)

        def generate_unicode_cmap(subset_index, charmap):
            # Make the ToUnicode CMap.
            last_ccode = -2
            unicode_groups = []
            for ccode in sorted(charmap.keys()):
                if ccode != last_ccode + 1:
                    unicode_groups.append([ccode, ccode])
                else:
                    unicode_groups[-1][1] = ccode
                last_ccode = ccode

            def _to_unicode(ccode):
                chars = self._character_tracker.subset_to_unicode(
                    filename, subset_index, ccode)
                hexstr = chars.encode('utf-16be').hex()
                return f'<{hexstr}>'

            width = 2 if fonttype == 3 else 4
            unicode_bfrange = []
            for start, end in unicode_groups:
                real_values = ' '.join(_to_unicode(x) for x in range(start, end+1))
                unicode_bfrange.append(
                    f'<{start:0{width}x}> <{end:0{width}x}> [{real_values}]')
            unicode_cmap = (self._identityToUnicodeCMap %
                            (len(unicode_groups),
                             '\n'.join(unicode_bfrange).encode('ascii')))

            return unicode_cmap

        def embedTTFType3(font, subset_index, charmap, descriptor):
            """The Type 3-specific part of embedding a Truetype font"""
            widthsObject = self.reserveObject('font widths')
            fontdescObject = self.reserveObject('font descriptor')
            fontdictObject = self.reserveObject('font dictionary')
            charprocsObject = self.reserveObject('character procs')
            toUnicodeMapObject = self.reserveObject('ToUnicode map')
            differencesArray = []
            firstchar, lastchar = min(charmap), max(charmap)
            bbox = [cvt(x, nearest=False) for x in font.bbox]

            fontdict = {
                'Type': Name('Font'),
                'BaseFont': ps_name,
                'FirstChar': firstchar,
                'LastChar': lastchar,
                'FontDescriptor': fontdescObject,
                'Subtype': Name('Type3'),
                'Name': descriptor['FontName'],
                'FontBBox': bbox,
                'FontMatrix': [.001, 0, 0, .001, 0, 0],
                'CharProcs': charprocsObject,
                'Encoding': {
                    'Type': Name('Encoding'),
                    'Differences': differencesArray},
                'Widths': widthsObject,
                'ToUnicode': toUnicodeMapObject,
            }

            # Make the "Widths" array
            def get_char_width(charcode):
                width = font.load_glyph(
                    charmap.get(charcode, 0),
                    flags=LoadFlags.NO_SCALE | LoadFlags.NO_HINTING).horiAdvance
                return cvt(width)
            widths = [get_char_width(charcode)
                      for charcode in range(firstchar, lastchar+1)]
            descriptor['MaxWidth'] = max(widths)

            # Make the "Differences" array with the whole set of character codes that we
            # need from this font.
            differences = sorted([
                (ccode, font.get_glyph_name(gind)) for ccode, gind in charmap.items()
            ])
            last_c = -2
            for c, name in differences:
                if c != last_c + 1:
                    differencesArray.append(c)
                differencesArray.append(Name(name))
                last_c = c

            # Make the charprocs array.
            rawcharprocs = _get_pdf_charprocs(filename, charmap.values())
            charprocs = {}
            for charname in sorted(rawcharprocs):
                stream = rawcharprocs[charname]
                charprocObject = self.reserveObject('charProc')
                self.outputStream(charprocObject, stream)
                charprocs[charname] = charprocObject

            unicode_cmap = generate_unicode_cmap(subset_index, charmap)

            # Write everything out
            self.writeObject(fontdictObject, fontdict)
            self.writeObject(fontdescObject, descriptor)
            self.writeObject(widthsObject, widths)
            self.writeObject(charprocsObject, charprocs)
            self.outputStream(toUnicodeMapObject, unicode_cmap)

            return fontdictObject

        def embedTTFType42(font, subset_index, charmap, descriptor):
            """The Type 42-specific part of embedding a Truetype font"""
            fontdescObject = self.reserveObject('font descriptor')
            cidFontDictObject = self.reserveObject('CID font dictionary')
            type0FontDictObject = self.reserveObject('Type 0 font dictionary')
            cidToGidMapObject = self.reserveObject('CIDToGIDMap stream')
            fontfileObject = self.reserveObject('font file stream')
            wObject = self.reserveObject('Type 0 widths')
            toUnicodeMapObject = self.reserveObject('ToUnicode map')

            _log.debug("SUBSET %r:%d characters: %s", filename, subset_index, charmap)
            with _backend_pdf_ps.get_glyphs_subset(filename,
                                                   charmap.values()) as subset:
                fontdata = _backend_pdf_ps.font_as_file(subset.font)
            _log.debug(
                "SUBSET %r:%d %d -> %d", filename, subset_index,
                os.stat(filename).st_size, fontdata.getbuffer().nbytes
            )

            # reload the font object from the subset
            # (all the necessary data could probably be obtained directly
            # using fontLib.ttLib)
            font = FT2Font(fontdata)

            cidFontDict = {
                'Type': Name('Font'),
                'Subtype': Name('CIDFontType2'),
                'BaseFont': ps_name,
                'CIDSystemInfo': {
                    'Registry': 'Adobe',
                    'Ordering': 'Identity',
                    'Supplement': 0},
                'FontDescriptor': fontdescObject,
                'W': wObject,
                'CIDToGIDMap': cidToGidMapObject
                }

            type0FontDict = {
                'Type': Name('Font'),
                'Subtype': Name('Type0'),
                'BaseFont': ps_name,
                'Encoding': Name('Identity-H'),
                'DescendantFonts': [cidFontDictObject],
                'ToUnicode': toUnicodeMapObject
                }

            # Make fontfile stream
            descriptor['FontFile2'] = fontfileObject
            self.outputStream(
                fontfileObject, fontdata.getvalue(),
                extra={'Length1': fontdata.getbuffer().nbytes})

            # Make the 'W' (Widths) array and CidToGidMap at the same time.
            cid_to_gid_map = ['\0'] * 65536
            widths = []
            max_ccode = 0
            for ccode, orig_glyph_index in charmap.items():
                gind = subset.glyph_index_map[orig_glyph_index]
                glyph = font.load_glyph(gind,
                                        flags=LoadFlags.NO_SCALE | LoadFlags.NO_HINTING)
                widths.append((ccode, cvt(glyph.horiAdvance)))
                cid_to_gid_map[ccode] = chr(gind)
                max_ccode = max(ccode, max_ccode)
            widths.sort()
            cid_to_gid_map = cid_to_gid_map[:max_ccode + 1]

            last_ccode = -2
            w = []
            max_width = 0
            for ccode, width in widths:
                if ccode != last_ccode + 1:
                    w.append(ccode)
                    w.append([width])
                else:
                    w[-1].append(width)
                max_width = max(max_width, width)
                last_ccode = ccode

            # CIDToGIDMap stream
            cid_to_gid_map = "".join(cid_to_gid_map).encode("utf-16be")
            self.outputStream(cidToGidMapObject, cid_to_gid_map)

            # ToUnicode CMap
            unicode_cmap = generate_unicode_cmap(subset_index, charmap)
            self.outputStream(toUnicodeMapObject, unicode_cmap)

            descriptor['MaxWidth'] = max_width

            # Write everything out
            self.writeObject(cidFontDictObject, cidFontDict)
            self.writeObject(type0FontDictObject, type0FontDict)
            self.writeObject(fontdescObject, descriptor)
            self.writeObject(wObject, w)

            return type0FontDictObject

        # Beginning of main embedTTF function...

        ps_name = self._get_subsetted_psname(font.postscript_name, charmap)
        ps_name = ps_name.encode('ascii', 'replace')
        ps_name = Name(ps_name)
        pclt = font.get_sfnt_table('pclt') or {'capHeight': 0, 'xHeight': 0}
        post = font.get_sfnt_table('post') or {'italicAngle': (0, 0)}
        ff = font.face_flags
        sf = font.style_flags

        flags = 0
        symbolic = False  # ps_name.name in ('Cmsy10', 'Cmmi10', 'Cmex10')
        if FaceFlags.FIXED_WIDTH in ff:
            flags |= 1 << 0
        if 0:  # TODO: serif
            flags |= 1 << 1
        if symbolic:
            flags |= 1 << 2
        else:
            flags |= 1 << 5
        if StyleFlags.ITALIC in sf:
            flags |= 1 << 6
        if 0:  # TODO: all caps
            flags |= 1 << 16
        if 0:  # TODO: small caps
            flags |= 1 << 17
        if 0:  # TODO: force bold
            flags |= 1 << 18

        descriptor = {
            'Type': Name('FontDescriptor'),
            'FontName': ps_name,
            'Flags': flags,
            'FontBBox': [cvt(x, nearest=False) for x in font.bbox],
            'Ascent': cvt(font.ascender, nearest=False),
            'Descent': cvt(font.descender, nearest=False),
            'CapHeight': cvt(pclt['capHeight'], nearest=False),
            'XHeight': cvt(pclt['xHeight']),
            'ItalicAngle': post['italicAngle'][1],  # ???
            'StemV': 0  # ???
            }

        if fonttype == 3:
            return embedTTFType3(font, subset_index, charmap, descriptor)
        elif fonttype == 42:
            return embedTTFType42(font, subset_index, charmap, descriptor)

    def alphaState(self, alpha):
        """Return name of an ExtGState that sets alpha to the given value."""

        state = self.alphaStates.get(alpha, None)
        if state is not None:
            return state[0]

        name = next(self._alpha_state_seq)
        self.alphaStates[alpha] = \
            (name, {'Type': Name('ExtGState'),
                    'CA': alpha[0], 'ca': alpha[1]})
        return name

    def _soft_mask_state(self, smask):
        """
        Return an ExtGState that sets the soft mask to the given shading.

        Parameters
        ----------
        smask : Reference
            Reference to a shading in DeviceGray color space, whose luminosity
            is to be used as the alpha channel.

        Returns
        -------
        Name
        """

        state = self._soft_mask_states.get(smask, None)
        if state is not None:
            return state[0]

        name = next(self._soft_mask_seq)
        groupOb = self.reserveObject('transparency group for soft mask')
        self._soft_mask_states[smask] = (
            name,
            {
                'Type': Name('ExtGState'),
                'AIS': False,
                'SMask': {
                    'Type': Name('Mask'),
                    'S': Name('Luminosity'),
                    'BC': [1],
                    'G': groupOb
                }
            }
        )
        self._soft_mask_groups.append((
            groupOb,
            {
                'Type': Name('XObject'),
                'Subtype': Name('Form'),
                'FormType': 1,
                'Group': {
                    'S': Name('Transparency'),
                    'CS': Name('DeviceGray')
                },
                'Matrix': [1, 0, 0, 1, 0, 0],
                'Resources': {'Shading': {'S': smask}},
                'BBox': [0, 0, 1, 1]
            },
            [Name('S'), Op.shading]
        ))
        return name

    def writeExtGSTates(self):
        self.writeObject(
            self._extGStateObject,
            dict([
                *self.alphaStates.values(),
                *self._soft_mask_states.values()
            ])
        )

    def _write_soft_mask_groups(self):
        for ob, attributes, content in self._soft_mask_groups:
            self.beginStream(ob.id, None, attributes)
            self.output(*content)
            self.endStream()

    def hatchPattern(self, hatch_style):
        # The colors may come in as numpy arrays, which aren't hashable
        edge, face, hatch, lw = hatch_style
        if edge is not None:
            edge = tuple(edge)
        if face is not None:
            face = tuple(face)
        hatch_style = (edge, face, hatch, lw)

        pattern = self._hatch_patterns.get(hatch_style, None)
        if pattern is not None:
            return pattern

        name = next(self._hatch_pattern_seq)
        self._hatch_patterns[hatch_style] = name
        return name

    hatchPatterns = _api.deprecated("3.10")(property(lambda self: {
        k: (e, f, h) for k, (e, f, h, l) in self._hatch_patterns.items()
    }))

    def writeHatches(self):
        hatchDict = dict()
        sidelen = 72.0
        for hatch_style, name in self._hatch_patterns.items():
            ob = self.reserveObject('hatch pattern')
            hatchDict[name] = ob
            res = {'Procsets':
                   [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]}
            self.beginStream(
                ob.id, None,
                {'Type': Name('Pattern'),
                 'PatternType': 1, 'PaintType': 1, 'TilingType': 1,
                 'BBox': [0, 0, sidelen, sidelen],
                 'XStep': sidelen, 'YStep': sidelen,
                 'Resources': res,
                 # Change origin to match Agg at top-left.
                 'Matrix': [1, 0, 0, 1, 0, self.height * 72]})

            stroke_rgb, fill_rgb, hatch, lw = hatch_style
            self.output(stroke_rgb[0], stroke_rgb[1], stroke_rgb[2],
                        Op.setrgb_stroke)
            if fill_rgb is not None:
                self.output(fill_rgb[0], fill_rgb[1], fill_rgb[2],
                            Op.setrgb_nonstroke,
                            0, 0, sidelen, sidelen, Op.rectangle,
                            Op.fill)
            self.output(stroke_rgb[0], stroke_rgb[1], stroke_rgb[2],
                        Op.setrgb_nonstroke)

            self.output(lw, Op.setlinewidth)

            self.output(*self.pathOperations(
                Path.hatch(hatch),
                Affine2D().scale(sidelen),
                simplify=False))
            self.output(Op.fill_stroke)

            self.endStream()
        self.writeObject(self.hatchObject, hatchDict)

    def addGouraudTriangles(self, points, colors):
        """
        Add a Gouraud triangle shading.

        Parameters
        ----------
        points : np.ndarray
            Triangle vertices, shape (n, 3, 2)
            where n = number of triangles, 3 = vertices, 2 = x, y.
        colors : np.ndarray
            Vertex colors, shape (n, 3, 1) or (n, 3, 4)
            as with points, but last dimension is either (gray,)
            or (r, g, b, alpha).

        Returns
        -------
        Name, Reference
        """
        name = Name('GT%d' % len(self.gouraudTriangles))
        ob = self.reserveObject(f'Gouraud triangle {name}')
        self.gouraudTriangles.append((name, ob, points, colors))
        return name, ob

    def writeGouraudTriangles(self):
        gouraudDict = dict()
        for name, ob, points, colors in self.gouraudTriangles:
            gouraudDict[name] = ob
            shape = points.shape
            flat_points = points.reshape((shape[0] * shape[1], 2))
            colordim = colors.shape[2]
            assert colordim in (1, 4)
            flat_colors = colors.reshape((shape[0] * shape[1], colordim))
            if colordim == 4:
                # strip the alpha channel
                colordim = 3
            points_min = np.min(flat_points, axis=0) - (1 << 8)
            points_max = np.max(flat_points, axis=0) + (1 << 8)
            factor = 0xffffffff / (points_max - points_min)

            self.beginStream(
                ob.id, None,
                {'ShadingType': 4,
                 'BitsPerCoordinate': 32,
                 'BitsPerComponent': 8,
                 'BitsPerFlag': 8,
                 'ColorSpace': Name(
                     'DeviceRGB' if colordim == 3 else 'DeviceGray'
                 ),
                 'AntiAlias': False,
                 'Decode': ([points_min[0], points_max[0],
                             points_min[1], points_max[1]]
                            + [0, 1] * colordim),
                 })

            streamarr = np.empty(
                (shape[0] * shape[1],),
                dtype=[('flags', 'u1'),
                       ('points', '>u4', (2,)),
                       ('colors', 'u1', (colordim,))])
            streamarr['flags'] = 0
            streamarr['points'] = (flat_points - points_min) * factor
            streamarr['colors'] = flat_colors[:, :colordim] * 255.0

            self.write(streamarr.tobytes())
            self.endStream()
        self.writeObject(self.gouraudObject, gouraudDict)

    def imageObject(self, image):
        """Return name of an image XObject representing the given image."""

        entry = self._images.get(id(image), None)
        if entry is not None:
            return entry[1]

        name = next(self._image_seq)
        ob = self.reserveObject(f'image {name}')
        self._images[id(image)] = (image, name, ob)
        return name

    def _unpack(self, im):
        """
        Unpack image array *im* into ``(data, alpha)``, which have shape
        ``(height, width, 3)`` (RGB) or ``(height, width, 1)`` (grayscale or
        alpha), except that alpha is None if the image is fully opaque.
        """
        im = im[::-1]
        if im.ndim == 2:
            return im, None
        else:
            rgb = im[:, :, :3]
            rgb = np.array(rgb, order='C')
            # PDF needs a separate alpha image
            if im.shape[2] == 4:
                alpha = im[:, :, 3][..., None]
                if np.all(alpha == 255):
                    alpha = None
                else:
                    alpha = np.array(alpha, order='C')
            else:
                alpha = None
            return rgb, alpha

    def _writePng(self, img):
        """
        Write the image *img* into the pdf file using png
        predictors with Flate compression.
        """
        buffer = BytesIO()
        img.save(buffer, format="png")
        buffer.seek(8)
        png_data = b''
        bit_depth = palette = None
        while True:
            length, type = struct.unpack(b'!L4s', buffer.read(8))
            if type in [b'IHDR', b'PLTE', b'IDAT']:
                data = buffer.read(length)
                if len(data) != length:
                    raise RuntimeError("truncated data")
                if type == b'IHDR':
                    bit_depth = int(data[8])
                elif type == b'PLTE':
                    palette = data
                elif type == b'IDAT':
                    png_data += data
            elif type == b'IEND':
                break
            else:
                buffer.seek(length, 1)
            buffer.seek(4, 1)   # skip CRC
        return png_data, bit_depth, palette

    def _writeImg(self, data, id, smask=None):
        """
        Write the image *data*, of shape ``(height, width, 1)`` (grayscale) or
        ``(height, width, 3)`` (RGB), as pdf object *id* and with the soft mask
        (alpha channel) *smask*, which should be either None or a ``(height,
        width, 1)`` array.
        """
        height, width, color_channels = data.shape
        obj = {'Type': Name('XObject'),
               'Subtype': Name('Image'),
               'Width': width,
               'Height': height,
               'ColorSpace': Name({1: 'DeviceGray', 3: 'DeviceRGB'}[color_channels]),
               'BitsPerComponent': 8}
        if smask:
            obj['SMask'] = smask
        if mpl.rcParams['pdf.compression']:
            if data.shape[-1] == 1:
                data = data.squeeze(axis=-1)
            png = {'Predictor': 10, 'Colors': color_channels, 'Columns': width}
            img = Image.fromarray(data)
            img_colors = img.getcolors(maxcolors=256)
            if color_channels == 3 and img_colors is not None:
                # Convert to indexed color if there are 256 colors or fewer. This can
                # significantly reduce the file size.
                num_colors = len(img_colors)
                palette = np.array([comp for _, color in img_colors for comp in color],
                                   dtype=np.uint8)
                palette24 = ((palette[0::3].astype(np.uint32) << 16) |
                             (palette[1::3].astype(np.uint32) << 8) |
                             palette[2::3])
                rgb24 = ((data[:, :, 0].astype(np.uint32) << 16) |
                         (data[:, :, 1].astype(np.uint32) << 8) |
                         data[:, :, 2])
                indices = np.argsort(palette24).astype(np.uint8)
                rgb8 = indices[np.searchsorted(palette24, rgb24, sorter=indices)]
                img = Image.fromarray(rgb8).convert("P")
                img.putpalette(palette)
                png_data, bit_depth, palette = self._writePng(img)
                if bit_depth is None or palette is None:
                    raise RuntimeError("invalid PNG header")
                palette = palette[:num_colors * 3]  # Trim padding; remove for Pillow>=9
                obj['ColorSpace'] = [Name('Indexed'), Name('DeviceRGB'),
                                     num_colors - 1, palette]
                obj['BitsPerComponent'] = bit_depth
                png['Colors'] = 1
                png['BitsPerComponent'] = bit_depth
            else:
                png_data, _, _ = self._writePng(img)
        else:
            png = None
        self.beginStream(
            id,
            self.reserveObject('length of image stream'),
            obj,
            png=png
            )
        if png:
            self.currentstream.write(png_data)
        else:
            self.currentstream.write(data.tobytes())
        self.endStream()

    def writeImages(self):
        for img, name, ob in self._images.values():
            data, adata = self._unpack(img)
            if adata is not None:
                smaskObject = self.reserveObject("smask")
                self._writeImg(adata, smaskObject.id)
            else:
                smaskObject = None
            self._writeImg(data, ob.id, smaskObject)

    def markerObject(self, path, trans, fill, stroke, lw, joinstyle,
                     capstyle):
        """Return name of a marker XObject representing the given path."""
        # self.markers used by markerObject, writeMarkers, close:
        # mapping from (path operations, fill?, stroke?) to
        #   [name, object reference, bounding box, linewidth]
        # This enables different draw_markers calls to share the XObject
        # if the gc is sufficiently similar: colors etc can vary, but
        # the choices of whether to fill and whether to stroke cannot.
        # We need a bounding box enclosing all of the XObject path,
        # but since line width may vary, we store the maximum of all
        # occurring line widths in self.markers.
        # close() is somewhat tightly coupled in that it expects the
        # first two components of each value in self.markers to be the
        # name and object reference.
        pathops = self.pathOperations(path, trans, simplify=False)
        key = (tuple(pathops), bool(fill), bool(stroke), joinstyle, capstyle)
        result = self.markers.get(key)
        if result is None:
            name = Name('M%d' % len(self.markers))
            ob = self.reserveObject('marker %d' % len(self.markers))
            bbox = path.get_extents(trans)
            self.markers[key] = [name, ob, bbox, lw]
        else:
            if result[-1] < lw:
                result[-1] = lw
            name = result[0]
        return name

    def writeMarkers(self):
        for ((pathops, fill, stroke, joinstyle, capstyle),
             (name, ob, bbox, lw)) in self.markers.items():
            # bbox wraps the exact limits of the control points, so half a line
            # will appear outside it. If the join style is miter and the line
            # is not parallel to the edge, then the line will extend even
            # further. From the PDF specification, Section 8.4.3.5, the miter
            # limit is miterLength / lineWidth and from Table 52, the default
            # is 10. With half the miter length outside, that works out to the
            # following padding:
            bbox = bbox.padded(lw * 5)
            self.beginStream(
                ob.id, None,
                {'Type': Name('XObject'), 'Subtype': Name('Form'),
                 'BBox': list(bbox.extents)})
            self.output(GraphicsContextPdf.joinstyles[joinstyle],
                        Op.setlinejoin)
            self.output(GraphicsContextPdf.capstyles[capstyle], Op.setlinecap)
            self.output(*pathops)
            self.output(Op.paint_path(fill, stroke))
            self.endStream()

    def pathCollectionObject(self, gc, path, trans, padding, filled, stroked):
        name = Name('P%d' % len(self.paths))
        ob = self.reserveObject('path %d' % len(self.paths))
        self.paths.append(
            (name, path, trans, ob, gc.get_joinstyle(), gc.get_capstyle(),
             padding, filled, stroked))
        return name

    def writePathCollectionTemplates(self):
        for (name, path, trans, ob, joinstyle, capstyle, padding, filled,
             stroked) in self.paths:
            pathops = self.pathOperations(path, trans, simplify=False)
            bbox = path.get_extents(trans)
            if not np.all(np.isfinite(bbox.extents)):
                extents = [0, 0, 0, 0]
            else:
                bbox = bbox.padded(padding)
                extents = list(bbox.extents)
            self.beginStream(
                ob.id, None,
                {'Type': Name('XObject'), 'Subtype': Name('Form'),
                 'BBox': extents})
            self.output(GraphicsContextPdf.joinstyles[joinstyle],
                        Op.setlinejoin)
            self.output(GraphicsContextPdf.capstyles[capstyle], Op.setlinecap)
            self.output(*pathops)
            self.output(Op.paint_path(filled, stroked))
            self.endStream()

    @staticmethod
    def pathOperations(path, transform, clip=None, simplify=None, sketch=None):
        return [Verbatim(_path.convert_to_string(
            path, transform, clip, simplify, sketch,
            6,
            [Op.moveto.value, Op.lineto.value, b'', Op.curveto.value,
             Op.closepath.value],
            True))]

    def writePath(self, path, transform, clip=False, sketch=None):
        if clip:
            clip = (0.0, 0.0, self.width * 72, self.height * 72)
            simplify = path.should_simplify
        else:
            clip = None
            simplify = False
        cmds = self.pathOperations(path, transform, clip, simplify=simplify,
                                   sketch=sketch)
        self.output(*cmds)

    def reserveObject(self, name=''):
        """
        Reserve an ID for an indirect object.

        The name is used for debugging in case we forget to print out
        the object with writeObject.
        """
        id = next(self._object_seq)
        self.xrefTable.append([None, 0, name])
        return Reference(id)

    def recordXref(self, id):
        self.xrefTable[id][0] = self.fh.tell() - self.tell_base

    def writeObject(self, object, contents):
        self.recordXref(object.id)
        object.write(contents, self)

    def writeXref(self):
        """Write out the xref table."""
        self.startxref = self.fh.tell() - self.tell_base
        self.write(b"xref\n0 %d\n" % len(self.xrefTable))
        for i, (offset, generation, name) in enumerate(self.xrefTable):
            if offset is None:
                raise AssertionError(
                    'No offset for object %d (%s)' % (i, name))
            else:
                key = b"f" if name == 'the zero object' else b"n"
                text = b"%010d %05d %b \n" % (offset, generation, key)
                self.write(text)

    def writeInfoDict(self):
        """Write out the info dictionary, checking it for good form"""

        self.infoObject = self.reserveObject('info')
        self.writeObject(self.infoObject, self.infoDict)

    def writeTrailer(self):
        """Write out the PDF trailer."""

        self.write(b"trailer\n")
        self.write(pdfRepr(
            {'Size': len(self.xrefTable),
             'Root': self.rootObject,
             'Info': self.infoObject}))
        # Could add 'ID'
        self.write(b"\nstartxref\n%d\n%%%%EOF\n" % self.startxref)


class RendererPdf(_backend_pdf_ps.RendererPDFPSBase):

    _afm_font_dir = cbook._get_data_path("fonts/pdfcorefonts")
    _use_afm_rc_name = "pdf.use14corefonts"

    def __init__(self, file, image_dpi, height, width):
        super().__init__(width, height)
        self.file = file
        self.gc = self.new_gc()
        self.image_dpi = image_dpi

    def finalize(self):
        self.file.output(*self.gc.finalize())

    def check_gc(self, gc, fillcolor=None):
        orig_fill = getattr(gc, '_fillcolor', (0., 0., 0.))
        gc._fillcolor = fillcolor

        orig_alphas = getattr(gc, '_effective_alphas', (1.0, 1.0))

        if gc.get_rgb() is None:
            # It should not matter what color here since linewidth should be
            # 0 unless affected by global settings in rcParams, hence setting
            # zero alpha just in case.
            gc.set_foreground((0, 0, 0, 0), isRGBA=True)

        if gc._forced_alpha:
            gc._effective_alphas = (gc._alpha, gc._alpha)
        elif fillcolor is None or len(fillcolor) < 4:
            gc._effective_alphas = (gc._rgb[3], 1.0)
        else:
            gc._effective_alphas = (gc._rgb[3], fillcolor[3])

        delta = self.gc.delta(gc)
        if delta:
            self.file.output(*delta)

        # Restore gc to avoid unwanted side effects
        gc._fillcolor = orig_fill
        gc._effective_alphas = orig_alphas

    def get_image_magnification(self):
        return self.image_dpi/72.0

    def draw_image(self, gc, x, y, im, transform=None):
        # docstring inherited

        h, w = im.shape[:2]
        if w == 0 or h == 0:
            return

        if transform is None:
            # If there's no transform, alpha has already been applied
            gc.set_alpha(1.0)

        self.check_gc(gc)

        w = 72.0 * w / self.image_dpi
        h = 72.0 * h / self.image_dpi

        imob = self.file.imageObject(im)

        if transform is None:
            self.file.output(Op.gsave,
                             w, 0, 0, h, x, y, Op.concat_matrix,
                             imob, Op.use_xobject, Op.grestore)
        else:
            tr1, tr2, tr3, tr4, tr5, tr6 = transform.frozen().to_values()

            self.file.output(Op.gsave,
                             1, 0, 0, 1, x, y, Op.concat_matrix,
                             tr1, tr2, tr3, tr4, tr5, tr6, Op.concat_matrix,
                             imob, Op.use_xobject, Op.grestore)

    def draw_path(self, gc, path, transform, rgbFace=None):
        # docstring inherited
        self.check_gc(gc, rgbFace)
        self.file.writePath(
            path, transform,
            rgbFace is None and gc.get_hatch_path() is None,
            gc.get_sketch_params())
        self.file.output(self.gc.paint())

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offset_trans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position, *, hatchcolors=None):
        # We can only reuse the objects if the presence of fill and
        # stroke (and the amount of alpha for each) is the same for
        # all of them
        can_do_optimization = True
        facecolors = np.asarray(facecolors)
        edgecolors = np.asarray(edgecolors)

        if hatchcolors is None:
            hatchcolors = []

        if not len(facecolors):
            filled = False
            can_do_optimization = not gc.get_hatch()
        else:
            if np.all(facecolors[:, 3] == facecolors[0, 3]):
                filled = facecolors[0, 3] != 0.0
            else:
                can_do_optimization = False

        if not len(edgecolors):
            stroked = False
        else:
            if np.all(np.asarray(linewidths) == 0.0):
                stroked = False
            elif np.all(edgecolors[:, 3] == edgecolors[0, 3]):
                stroked = edgecolors[0, 3] != 0.0
            else:
                can_do_optimization = False

        # Is the optimization worth it? Rough calculation:
        # cost of emitting a path in-line is len_path * uses_per_path
        # cost of XObject is len_path + 5 for the definition,
        #    uses_per_path for the uses
        len_path = len(paths[0].vertices) if len(paths) > 0 else 0
        uses_per_path = self._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        should_do_optimization = \
            len_path + uses_per_path + 5 < len_path * uses_per_path

        if (not can_do_optimization) or (not should_do_optimization):
            return RendererBase.draw_path_collection(
                self, gc, master_transform, paths, all_transforms,
                offsets, offset_trans, facecolors, edgecolors,
                linewidths, linestyles, antialiaseds, urls,
                offset_position, hatchcolors=hatchcolors)

        padding = np.max(linewidths)
        path_codes = []
        path_extents = []
        for i, (path, transform) in enumerate(self._iter_collection_raw_paths(
                master_transform, paths, all_transforms)):
            name = self.file.pathCollectionObject(
                gc, path, transform, padding, filled, stroked)
            path_codes.append(name)
            # Compute the extent of each marker path to enable per-marker
            # bounds checking. This allows us to skip markers that are
            # completely outside the visible canvas while preserving markers
            # that are partially visible.
            if len(path.vertices):
                bbox = path.get_extents(transform)
                # Store half-width and half-height for efficient bounds checking
                path_extents.append((bbox.width / 2, bbox.height / 2))
            else:
                path_extents.append((0, 0))

        # Create a mapping from path_id to extent for efficient lookup
        path_extent_map = dict(zip(path_codes, path_extents))

        canvas_width = self.file.width * 72
        canvas_height = self.file.height * 72

        output = self.file.output
        output(*self.gc.push())
        lastx, lasty = 0, 0
        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
                gc, path_codes, offsets, offset_trans,
                facecolors, edgecolors, linewidths, linestyles,
                antialiaseds, urls, offset_position, hatchcolors=hatchcolors):

            # Optimization: Fast path for markers with centers inside canvas.
            # This avoids the dictionary lookup for the common case where
            # markers are visible, improving performance for large scatter plots.
            if 0 <= xo <= canvas_width and 0 <= yo <= canvas_height:
                # Marker center is inside canvas - definitely render it
                self.check_gc(gc0, rgbFace)
                dx, dy = xo - lastx, yo - lasty
                output(1, 0, 0, 1, dx, dy, Op.concat_matrix, path_id,
                       Op.use_xobject)
                lastx, lasty = xo, yo
                continue

            # Marker center is outside canvas - check if partially visible.
            # Skip markers completely outside visible canvas bounds to reduce
            # PDF file size. Use per-marker extents to handle large markers
            # correctly: only skip if the marker's bounding box doesn't
            # intersect the canvas at all.
            extent_x, extent_y = path_extent_map[path_id]
            if not (-extent_x <= xo <= canvas_width + extent_x
                    and -extent_y <= yo <= canvas_height + extent_y):
                continue

            self.check_gc(gc0, rgbFace)
            dx, dy = xo - lastx, yo - lasty
            output(1, 0, 0, 1, dx, dy, Op.concat_matrix, path_id,
                   Op.use_xobject)
            lastx, lasty = xo, yo
        output(*self.gc.pop())

    def draw_markers(self, gc, marker_path, marker_trans, path, trans,
                     rgbFace=None):
        # docstring inherited

        # Same logic as in draw_path_collection
        len_marker_path = len(marker_path)
        uses = len(path)
        if len_marker_path * uses < len_marker_path + uses + 5:
            RendererBase.draw_markers(self, gc, marker_path, marker_trans,
                                      path, trans, rgbFace)
            return

        self.check_gc(gc, rgbFace)
        fill = gc.fill(rgbFace)
        stroke = gc.stroke()

        output = self.file.output
        marker = self.file.markerObject(
            marker_path, marker_trans, fill, stroke, self.gc._linewidth,
            gc.get_joinstyle(), gc.get_capstyle())

        output(Op.gsave)
        lastx, lasty = 0, 0
        for vertices, code in path.iter_segments(
                trans,
                clip=(0, 0, self.file.width*72, self.file.height*72),
                simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                if not (0 <= x <= self.file.width * 72
                        and 0 <= y <= self.file.height * 72):
                    continue
                dx, dy = x - lastx, y - lasty
                output(1, 0, 0, 1, dx, dy, Op.concat_matrix,
                       marker, Op.use_xobject)
                lastx, lasty = x, y
        output(Op.grestore)

    def draw_gouraud_triangles(self, gc, points, colors, trans):
        assert len(points) == len(colors)
        if len(points) == 0:
            return
        assert points.ndim == 3
        assert points.shape[1] == 3
        assert points.shape[2] == 2
        assert colors.ndim == 3
        assert colors.shape[1] == 3
        assert colors.shape[2] in (1, 4)

        shape = points.shape
        points = points.reshape((shape[0] * shape[1], 2))
        tpoints = trans.transform(points)
        tpoints = tpoints.reshape(shape)
        name, _ = self.file.addGouraudTriangles(tpoints, colors)
        output = self.file.output

        if colors.shape[2] == 1:
            # grayscale
            gc.set_alpha(1.0)
            self.check_gc(gc)
            output(name, Op.shading)
            return

        alpha = colors[0, 0, 3]
        if np.allclose(alpha, colors[:, :, 3]):
            # single alpha value
            gc.set_alpha(alpha)
            self.check_gc(gc)
            output(name, Op.shading)
        else:
            # varying alpha: use a soft mask
            alpha = colors[:, :, 3][:, :, None]
            _, smask_ob = self.file.addGouraudTriangles(tpoints, alpha)
            gstate = self.file._soft_mask_state(smask_ob)
            output(Op.gsave, gstate, Op.setgstate,
                   name, Op.shading,
                   Op.grestore)

    def _setup_textpos(self, x, y, angle, oldx=0, oldy=0, oldangle=0):
        if angle == oldangle == 0:
            self.file.output(x - oldx, y - oldy, Op.textpos)
        else:
            angle = math.radians(angle)
            self.file.output(math.cos(angle), math.sin(angle),
                             -math.sin(angle), math.cos(angle),
                             x, y, Op.textmatrix)
            self.file.output(0, 0, Op.textpos)

    def draw_mathtext(self, gc, x, y, s, prop, angle):
        # TODO: fix positioning and encoding
        width, height, descent, glyphs, rects = \
            self._text2path.mathtext_parser.parse(s, 72, prop)

        if gc.get_url() is not None:
            self.file._annotations[-1][1].append(_get_link_annotation(
                gc, x, y, width, height, angle))

        fonttype = mpl.rcParams['pdf.fonttype']

        # Set up a global transformation matrix for the whole math expression
        a = math.radians(angle)
        self.file.output(Op.gsave)
        self.file.output(math.cos(a), math.sin(a),
                         -math.sin(a), math.cos(a),
                         x, y, Op.concat_matrix)

        self.check_gc(gc, gc._rgb)
        prev_font = None, None
        oldx, oldy = 0, 0

        self.file.output(Op.begin_text)
        for font, fontsize, ccode, glyph_index, ox, oy in glyphs:
            subset_index, subset_charcode = self.file._character_tracker.track_glyph(
                font, ccode, glyph_index)
            font_path = FontPath(font.fname, font.face_index)
            self._setup_textpos(ox, oy, 0, oldx, oldy)
            oldx, oldy = ox, oy
            if (font_path, subset_index, fontsize) != prev_font:
                self.file.output(self.file.fontName(font_path, subset_index), fontsize,
                                 Op.selectfont)
                prev_font = font_path, subset_index, fontsize
            self.file.output(self._encode_glyphs([subset_charcode], fonttype),
                             Op.show)
        self.file.output(Op.end_text)

        # Draw any horizontal lines in the math layout
        for ox, oy, width, height in rects:
            self.file.output(Op.gsave, ox, oy, width, height,
                             Op.rectangle, Op.fill, Op.grestore)

        # Pop off the global transformation
        self.file.output(Op.grestore)

    def draw_tex(self, gc, x, y, s, prop, angle, *, mtext=None):
        # docstring inherited
        texmanager = self.get_texmanager()
        fontsize = prop.get_size_in_points()
        dvifile = texmanager.make_dvi(s, fontsize)
        with dviread.Dvi(dvifile, 72) as dvi:
            page, = dvi

        if gc.get_url() is not None:
            self.file._annotations[-1][1].append(_get_link_annotation(
                gc, x, y, page.width, page.height, angle))

        # Gather font information and do some setup for combining
        # characters into strings. The variable seq will contain a
        # sequence of font and text entries. A font entry is a list
        # ['font', name, size] where name is a Name object for the
        # font. A text entry is ['text', x, y, glyphs, x+w] where x
        # and y are the starting coordinates, w is the width, and
        # glyphs is a list; in this phase it will always contain just
        # one single-character string, but later it may have longer
        # strings interspersed with kern amounts.
        oldfont, seq = None, []
        for text in page.text:
            if text.font != oldfont:
                pdfname = self.file.dviFontName(text.font)
                seq += [['font', pdfname, text.font.size]]
                oldfont = text.font
            seq += [['text', text.x, text.y, [bytes([text.glyph])], text.x+text.width]]
            self.file._character_tracker.track_glyph(text.font, text.glyph, text.index)

        # Find consecutive text strings with constant y coordinate and
        # combine into a sequence of strings and kerns, or just one
        # string (if any kerns would be less than 0.1 points).
        i, curx, fontsize = 0, 0, None
        while i < len(seq)-1:
            elt, nxt = seq[i:i+2]
            if elt[0] == 'font':
                fontsize = elt[2]
            elif elt[0] == nxt[0] == 'text' and elt[2] == nxt[2]:
                offset = elt[4] - nxt[1]
                if abs(offset) < 0.1:
                    elt[3][-1] += nxt[3][0]
                    elt[4] += nxt[4]-nxt[1]
                else:
                    elt[3] += [offset*1000.0/fontsize, nxt[3][0]]
                    elt[4] = nxt[4]
                del seq[i+1]
                continue
            i += 1

        # Create a transform to map the dvi contents to the canvas.
        mytrans = Affine2D().rotate_deg(angle).translate(x, y)

        # Output the text.
        self.check_gc(gc, gc._rgb)
        self.file.output(Op.begin_text)
        curx, cury, oldx, oldy = 0, 0, 0, 0
        for elt in seq:
            if elt[0] == 'font':
                self.file.output(elt[1], elt[2], Op.selectfont)
            elif elt[0] == 'text':
                curx, cury = mytrans.transform((elt[1], elt[2]))
                self._setup_textpos(curx, cury, angle, oldx, oldy)
                oldx, oldy = curx, cury
                if len(elt[3]) == 1:
                    self.file.output(elt[3][0], Op.show)
                else:
                    self.file.output(elt[3], Op.showkern)
            else:
                assert False
        self.file.output(Op.end_text)

        # Then output the boxes (e.g., variable-length lines of square
        # roots).
        boxgc = self.new_gc()
        boxgc.copy_properties(gc)
        boxgc.set_linewidth(0)
        pathops = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                   Path.CLOSEPOLY]
        for x1, y1, h, w in page.boxes:
            path = Path([[x1, y1], [x1+w, y1], [x1+w, y1+h], [x1, y1+h],
                         [0, 0]], pathops)
            self.draw_path(boxgc, path, mytrans, gc._rgb)

    def _encode_glyphs(self, subset, fonttype):
        if fonttype in (1, 3):
            return bytes(subset)
        return b''.join(glyph.to_bytes(2, 'big') for glyph in subset)

    def encode_string(self, s, fonttype):
        encoding = {1: 'cp1252', 3: 'latin-1', 42: 'utf-16be'}[fonttype]
        return s.encode(encoding, 'replace')

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        # docstring inherited

        # TODO: combine consecutive texts into one BT/ET delimited section

        self.check_gc(gc, gc._rgb)
        if ismath:
            return self.draw_mathtext(gc, x, y, s, prop, angle)

        fontsize = prop.get_size_in_points()
        if mtext is not None:
            features = mtext.get_fontfeatures()
            language = mtext.get_language()
        else:
            features = language = None

        # For Type-1 fonts, emit the whole string at once without manual kerning.
        if mpl.rcParams['pdf.use14corefonts']:
            font = self._get_font_afm(prop)
            self.file.output(Op.begin_text,
                             self.file.fontName(prop), fontsize, Op.selectfont)
            self._setup_textpos(x, y, angle)
            self.file.output(self.encode_string(s, fonttype=1),
                             Op.show, Op.end_text)

        # A sequence of characters is broken into multiple chunks. The chunking
        # serves two purposes:
        #   - For Type 3 fonts, there is no way to access multibyte characters, as they
        #     cannot have a CIDMap. Therefore, in this case we break the string into
        #     chunks, where each chunk contains a string of consecutive 1-byte
        #     characters in a 256-character subset of the font. A distinct version of
        #     the original font is created for each 256-character subset.
        #   - A sequence of characters is split into chunks to allow for kerning
        #     adjustments between consecutive chunks.
        #
        # Each chunk is emitted with the regular text show command (TJ) with appropriate
        # kerning between chunks.
        else:
            font = self._get_font_ttf(prop)
            fonttype = mpl.rcParams['pdf.fonttype']

            def output_singlebyte_chunk(kerns_or_chars):
                if not kerns_or_chars:
                    return
                self.file.output(
                    # See pdf spec "Text space details" for the 1000/fontsize
                    # (aka. 1000/T_fs) factor.
                    [(-1000 * next(group) / fontsize) if tp == float  # a kern
                     else self._encode_glyphs(group, fonttype)
                     for tp, group in itertools.groupby(kerns_or_chars, type)],
                    Op.showkern)
                kerns_or_chars.clear()
            # Do the rotation and global translation as a single matrix
            # concatenation up front
            self.file.output(Op.gsave)
            a = math.radians(angle)
            self.file.output(math.cos(a), math.sin(a),
                             -math.sin(a), math.cos(a),
                             x, y, Op.concat_matrix)
            # List of [prev_kern, char, char, ...] w/o zero kerns.
            singlebyte_chunk = []
            prev_font = None
            prev_start_x = 0
            # Emit all the characters in a BT/ET group.
            self.file.output(Op.begin_text)
            for item in _text_helpers.layout(s, font, features=features,
                                             language=language):
                subset, charcode = self.file._character_tracker.track_glyph(
                    item.ft_object, item.char, item.glyph_index)
                if (item.ft_object, subset) != prev_font:
                    output_singlebyte_chunk(singlebyte_chunk)
                    font_path = FontPath(item.ft_object.fname,
                                         item.ft_object.face_index)
                    ft_name = self.file.fontName(font_path, subset)
                    self.file.output(ft_name, fontsize, Op.selectfont)
                    self._setup_textpos(item.x, 0, 0, prev_start_x, 0, 0)
                    prev_font = (item.ft_object, subset)
                    prev_start_x = item.x
                if item.y:
                    output_singlebyte_chunk(singlebyte_chunk)
                    self.file.output(item.y, Op.textrise)
                if item.prev_kern:
                    singlebyte_chunk.append(item.prev_kern)
                singlebyte_chunk.append(charcode)
                if item.y:
                    output_singlebyte_chunk(singlebyte_chunk)
                    self.file.output(0, Op.textrise)
            output_singlebyte_chunk(singlebyte_chunk)
            self.file.output(Op.end_text)
            self.file.output(Op.grestore)

        if gc.get_url() is not None:
            font.set_text(s, features=features, language=language)
            width, height = font.get_width_height()
            self.file._annotations[-1][1].append(_get_link_annotation(
                gc, x, y, width / 64, height / 64, angle))

    def new_gc(self):
        # docstring inherited
        return GraphicsContextPdf(self.file)


class GraphicsContextPdf(GraphicsContextBase):

    def __init__(self, file):
        super().__init__()
        self._fillcolor = (0.0, 0.0, 0.0)
        self._effective_alphas = (1.0, 1.0)
        self.file = file
        self.parent = None

    def __repr__(self):
        d = dict(self.__dict__)
        del d['file']
        del d['parent']
        return repr(d)

    def stroke(self):
        """
        Predicate: does the path need to be stroked (its outline drawn)?
        This tests for the various conditions that disable stroking
        the path, in which case it would presumably be filled.
        """
        # _linewidth > 0: in pdf a line of width 0 is drawn at minimum
        #   possible device width, but e.g., agg doesn't draw at all
        return (self._linewidth > 0 and self._alpha > 0 and
                (len(self._rgb) <= 3 or self._rgb[3] != 0.0))

    def fill(self, *args):
        """
        Predicate: does the path need to be filled?

        An optional argument can be used to specify an alternative
        _fillcolor, as needed by RendererPdf.draw_markers.
        """
        if len(args):
            _fillcolor = args[0]
        else:
            _fillcolor = self._fillcolor
        return (self._hatch or
                (_fillcolor is not None and
                 (len(_fillcolor) <= 3 or _fillcolor[3] != 0.0)))

    def paint(self):
        """
        Return the appropriate pdf operator to cause the path to be
        stroked, filled, or both.
        """
        return Op.paint_path(self.fill(), self.stroke())

    capstyles = {'butt': 0, 'round': 1, 'projecting': 2}
    joinstyles = {'miter': 0, 'round': 1, 'bevel': 2}

    def capstyle_cmd(self, style):
        return [self.capstyles[style], Op.setlinecap]

    def joinstyle_cmd(self, style):
        return [self.joinstyles[style], Op.setlinejoin]

    def linewidth_cmd(self, width):
        return [width, Op.setlinewidth]

    def dash_cmd(self, dashes):
        offset, dash = dashes
        if dash is None:
            dash = []
            offset = 0
        return [list(dash), offset, Op.setdash]

    def alpha_cmd(self, alpha, forced, effective_alphas):
        name = self.file.alphaState(effective_alphas)
        return [name, Op.setgstate]

    def hatch_cmd(self, hatch, hatch_color, hatch_linewidth):
        if not hatch:
            if self._fillcolor is not None:
                return self.fillcolor_cmd(self._fillcolor)
            else:
                return [Name('DeviceRGB'), Op.setcolorspace_nonstroke]
        else:
            hatch_style = (hatch_color, self._fillcolor, hatch, hatch_linewidth)
            name = self.file.hatchPattern(hatch_style)
            return [Name('Pattern'), Op.setcolorspace_nonstroke,
                    name, Op.setcolor_nonstroke]

    def rgb_cmd(self, rgb):
        if mpl.rcParams['pdf.inheritcolor']:
            return []
        if rgb[0] == rgb[1] == rgb[2]:
            return [rgb[0], Op.setgray_stroke]
        else:
            return [*rgb[:3], Op.setrgb_stroke]

    def fillcolor_cmd(self, rgb):
        if rgb is None or mpl.rcParams['pdf.inheritcolor']:
            return []
        elif rgb[0] == rgb[1] == rgb[2]:
            return [rgb[0], Op.setgray_nonstroke]
        else:
            return [*rgb[:3], Op.setrgb_nonstroke]

    def push(self):
        parent = GraphicsContextPdf(self.file)
        parent.copy_properties(self)
        parent.parent = self.parent
        self.parent = parent
        return [Op.gsave]

    def pop(self):
        assert self.parent is not None
        self.copy_properties(self.parent)
        self.parent = self.parent.parent
        return [Op.grestore]

    def clip_cmd(self, cliprect, clippath):
        """Set clip rectangle. Calls `.pop()` and `.push()`."""
        cmds = []
        # Pop graphics state until we hit the right one or the stack is empty
        while ((self._cliprect, self._clippath) != (cliprect, clippath)
                and self.parent is not None):
            cmds.extend(self.pop())
        # Unless we hit the right one, set the clip polygon
        if ((self._cliprect, self._clippath) != (cliprect, clippath) or
                self.parent is None):
            cmds.extend(self.push())
            if self._cliprect != cliprect:
                cmds.extend([cliprect, Op.rectangle, Op.clip, Op.endpath])
            if self._clippath != clippath:
                path, affine = clippath.get_transformed_path_and_affine()
                cmds.extend(
                    PdfFile.pathOperations(path, affine, simplify=False) +
                    [Op.clip, Op.endpath])
        return cmds

    commands = (
        # must come first since may pop
        (('_cliprect', '_clippath'), clip_cmd),
        (('_alpha', '_forced_alpha', '_effective_alphas'), alpha_cmd),
        (('_capstyle',), capstyle_cmd),
        (('_fillcolor',), fillcolor_cmd),
        (('_joinstyle',), joinstyle_cmd),
        (('_linewidth',), linewidth_cmd),
        (('_dashes',), dash_cmd),
        (('_rgb',), rgb_cmd),
        # must come after fillcolor and rgb
        (('_hatch', '_hatch_color', '_hatch_linewidth'), hatch_cmd),
    )

    def delta(self, other):
        """
        Copy properties of other into self and return PDF commands
        needed to transform *self* into *other*.
        """
        cmds = []
        fill_performed = False
        for params, cmd in self.commands:
            different = False
            for p in params:
                ours = getattr(self, p)
                theirs = getattr(other, p)
                try:
                    if ours is None or theirs is None:
                        different = ours is not theirs
                    else:
                        different = bool(ours != theirs)
                except ValueError:
                    ours = np.asarray(ours)
                    theirs = np.asarray(theirs)
                    different = (ours.shape != theirs.shape or
                                 np.any(ours != theirs))
                if different:
                    break

            # Need to update hatching if we also updated fillcolor
            if cmd.__name__ == 'hatch_cmd' and fill_performed:
                different = True

            if different:
                if cmd.__name__ == 'fillcolor_cmd':
                    fill_performed = True
                theirs = [getattr(other, p) for p in params]
                cmds.extend(cmd(self, *theirs))
                for p in params:
                    setattr(self, p, getattr(other, p))
        return cmds

    def copy_properties(self, other):
        """
        Copy properties of other into self.
        """
        super().copy_properties(other)
        fillcolor = getattr(other, '_fillcolor', self._fillcolor)
        effective_alphas = getattr(other, '_effective_alphas',
                                   self._effective_alphas)
        self._fillcolor = fillcolor
        self._effective_alphas = effective_alphas

    def finalize(self):
        """
        Make sure every pushed graphics state is popped.
        """
        cmds = []
        while self.parent is not None:
            cmds.extend(self.pop())
        return cmds


class PdfPages:
    """
    A multi-page PDF file.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> # Initialize:
    >>> with PdfPages('foo.pdf') as pdf:
    ...     # As many times as you like, create a figure fig and save it:
    ...     fig = plt.figure()
    ...     pdf.savefig(fig)
    ...     # When no figure is specified the current figure is saved
    ...     pdf.savefig()

    Notes
    -----
    In reality `PdfPages` is a thin wrapper around `PdfFile`, in order to avoid
    confusion when using `~.pyplot.savefig` and forgetting the format argument.
    """

    @_api.delete_parameter("3.10", "keep_empty",
                           addendum="This parameter does nothing.")
    def __init__(self, filename, keep_empty=None, metadata=None):
        """
        Create a new PdfPages object.

        Parameters
        ----------
        filename : str or path-like or file-like
            Plots using `PdfPages.savefig` will be written to a file at this location.
            The file is opened when a figure is saved for the first time (overwriting
            any older file with the same name).

        metadata : dict, optional
            Information dictionary object (see PDF reference section 10.2.1
            'Document Information Dictionary'), e.g.:
            ``{'Creator': 'My software', 'Author': 'Me', 'Title': 'Awesome'}``.

            The standard keys are 'Title', 'Author', 'Subject', 'Keywords',
            'Creator', 'Producer', 'CreationDate', 'ModDate', and
            'Trapped'. Values have been predefined for 'Creator', 'Producer'
            and 'CreationDate'. They can be removed by setting them to `None`.
        """
        self._filename = filename
        self._metadata = metadata
        self._file = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _ensure_file(self):
        if self._file is None:
            self._file = PdfFile(self._filename, metadata=self._metadata)  # init.
        return self._file

    def close(self):
        """
        Finalize this object, making the underlying file a complete
        PDF file.
        """
        if self._file is not None:
            self._file.finalize()
            self._file.close()
            self._file = None

    def infodict(self):
        """
        Return a modifiable information dictionary object
        (see PDF reference section 10.2.1 'Document Information
        Dictionary').
        """
        return self._ensure_file().infoDict

    def savefig(self, figure=None, **kwargs):
        """
        Save a `.Figure` to this file as a new page.

        Any other keyword arguments are passed to `~.Figure.savefig`.

        Parameters
        ----------
        figure : `.Figure` or int, default: the active figure
            The figure, or index of the figure, that is saved to the file.
        """
        if not isinstance(figure, Figure):
            if figure is None:
                manager = Gcf.get_active()
            else:
                manager = Gcf.get_fig_manager(figure)
            if manager is None:
                raise ValueError(f"No figure {figure}")
            figure = manager.canvas.figure
        # Force use of pdf backend, as PdfPages is tightly coupled with it.
        figure.savefig(self, format="pdf", backend="pdf", **kwargs)

    def get_pagecount(self):
        """Return the current number of pages in the multipage pdf file."""
        return len(self._ensure_file().pageList)

    def attach_note(self, text, positionRect=[-100, -100, 0, 0]):
        """
        Add a new text note to the page to be saved next. The optional
        positionRect specifies the position of the new note on the
        page. It is outside the page per default to make sure it is
        invisible on printouts.
        """
        self._ensure_file().newTextnote(text, positionRect)


class FigureCanvasPdf(FigureCanvasBase):
    # docstring inherited

    fixed_dpi = 72
    filetypes = {'pdf': 'Portable Document Format'}

    def get_default_filetype(self):
        return 'pdf'

    def print_pdf(self, filename, *,
                  bbox_inches_restore=None, metadata=None):

        dpi = self.figure.dpi
        self.figure.dpi = 72  # there are 72 pdf points to an inch
        width, height = self.figure.get_size_inches()
        if isinstance(filename, PdfPages):
            file = filename._ensure_file()
        else:
            file = PdfFile(filename, metadata=metadata)
        try:
            file.newPage(width, height)
            renderer = MixedModeRenderer(
                self.figure, width, height, dpi,
                RendererPdf(file, dpi, height, width),
                bbox_inches_restore=bbox_inches_restore)
            self.figure.draw(renderer)
            renderer.finalize()
            if not isinstance(filename, PdfPages):
                file.finalize()
        finally:
            if isinstance(filename, PdfPages):  # finish off this page
                file.endStream()
            else:            # we opened the file above; now finish it off
                file.close()

    def draw(self):
        self.figure.draw_without_rendering()
        return super().draw()


FigureManagerPdf = FigureManagerBase


@_Backend.export
class _BackendPdf(_Backend):
    FigureCanvas = FigureCanvasPdf
