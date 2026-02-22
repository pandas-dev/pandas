from __future__ import annotations

import logging
import os
import traceback
from io import BytesIO, StringIO, UnsupportedOperation
from typing import TYPE_CHECKING, TypedDict, TypeVar, overload

from fontTools.config import Config
from fontTools.misc import xmlWriter
from fontTools.misc.configTools import AbstractConfig
from fontTools.misc.loggingTools import deprecateArgument
from fontTools.misc.textTools import Tag, byteord, tostr
from fontTools.ttLib import TTLibError
from fontTools.ttLib.sfnt import SFNTReader, SFNTWriter
from fontTools.ttLib.ttGlyphSet import (
    _TTGlyph,  # noqa: F401
    _TTGlyphSet,
    _TTGlyphSetCFF,
    _TTGlyphSetGlyf,
    _TTGlyphSetVARC,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, MutableMapping
    from types import ModuleType, TracebackType
    from typing import Any, BinaryIO, Literal, Sequence, TextIO

    from typing_extensions import Self, Unpack

    from fontTools.ttLib.tables import (
        B_A_S_E_,
        C_B_D_T_,
        C_B_L_C_,
        C_F_F_,
        C_F_F__2,
        C_O_L_R_,
        C_P_A_L_,
        D_S_I_G_,
        E_B_D_T_,
        E_B_L_C_,
        F_F_T_M_,
        G_D_E_F_,
        G_M_A_P_,
        G_P_K_G_,
        G_P_O_S_,
        G_S_U_B_,
        G_V_A_R_,
        H_V_A_R_,
        J_S_T_F_,
        L_T_S_H_,
        M_A_T_H_,
        M_E_T_A_,
        M_V_A_R_,
        S_I_N_G_,
        S_T_A_T_,
        S_V_G_,
        T_S_I__0,
        T_S_I__1,
        T_S_I__2,
        T_S_I__3,
        T_S_I__5,
        T_S_I_B_,
        T_S_I_C_,
        T_S_I_D_,
        T_S_I_J_,
        T_S_I_P_,
        T_S_I_S_,
        T_S_I_V_,
        T_T_F_A_,
        V_A_R_C_,
        V_D_M_X_,
        V_O_R_G_,
        V_V_A_R_,
        D__e_b_g,
        F__e_a_t,
        G__l_a_t,
        G__l_o_c,
        O_S_2f_2,
        S__i_l_f,
        S__i_l_l,
        _a_n_k_r,
        _a_v_a_r,
        _b_s_l_n,
        _c_i_d_g,
        _c_m_a_p,
        _c_v_a_r,
        _c_v_t,
        _f_e_a_t,
        _f_p_g_m,
        _f_v_a_r,
        _g_a_s_p,
        _g_c_i_d,
        _g_l_y_f,
        _g_v_a_r,
        _h_d_m_x,
        _h_e_a_d,
        _h_h_e_a,
        _h_m_t_x,
        _k_e_r_n,
        _l_c_a_r,
        _l_o_c_a,
        _l_t_a_g,
        _m_a_x_p,
        _m_e_t_a,
        _m_o_r_t,
        _m_o_r_x,
        _n_a_m_e,
        _o_p_b_d,
        _p_o_s_t,
        _p_r_e_p,
        _p_r_o_p,
        _s_b_i_x,
        _t_r_a_k,
        _v_h_e_a,
        _v_m_t_x,
    )
    from fontTools.ttLib.tables.DefaultTable import DefaultTable

    _VT_co = TypeVar("_VT_co", covariant=True)  # Value type covariant containers.

log = logging.getLogger(__name__)


_NumberT = TypeVar("_NumberT", bound=float)


class TTFont(object):
    """Represents a TrueType font.

    The object manages file input and output, and offers a convenient way of
    accessing tables. Tables will be only decompiled when necessary, ie. when
    they're actually accessed. This means that simple operations can be extremely fast.

    Example usage:

    .. code-block:: pycon

        >>>
        >> from fontTools import ttLib
        >> tt = ttLib.TTFont("afont.ttf") # Load an existing font file
        >> tt['maxp'].numGlyphs
        242
        >> tt['OS/2'].achVendID
        'B&H\000'
        >> tt['head'].unitsPerEm
        2048

    For details of the objects returned when accessing each table, see the
    :doc:`tables </ttLib/tables>` documentation.
    To add a table to the font, use the :py:func:`newTable` function:

    .. code-block:: pycon

        >>>
        >> os2 = newTable("OS/2")
        >> os2.version = 4
        >> # set other attributes
        >> font["OS/2"] = os2

    TrueType fonts can also be serialized to and from XML format (see also the
    :doc:`ttx </ttx>` binary):

    .. code-block:: pycon

        >>
        >> tt.saveXML("afont.ttx")
        Dumping 'LTSH' table...
        Dumping 'OS/2' table...
        [...]

        >> tt2 = ttLib.TTFont() # Create a new font object
        >> tt2.importXML("afont.ttx")
        >> tt2['maxp'].numGlyphs
        242

    The TTFont object may be used as a context manager; this will cause the file
    reader to be closed after the context ``with`` block is exited::

            with TTFont(filename) as f:
                    # Do stuff

    Args:
            file: When reading a font from disk, either a pathname pointing to a file,
                    or a readable file object.
            res_name_or_index: If running on a Macintosh, either a sfnt resource name or
                    an sfnt resource index number. If the index number is zero, TTLib will
                    autodetect whether the file is a flat file or a suitcase. (If it is a suitcase,
                    only the first 'sfnt' resource will be read.)
            sfntVersion (str): When constructing a font object from scratch, sets the four-byte
                    sfnt magic number to be used. Defaults to ``\0\1\0\0`` (TrueType). To create
                    an OpenType file, use ``OTTO``.
            flavor (str): Set this to ``woff`` when creating a WOFF file or ``woff2`` for a WOFF2
                    file.
            checkChecksums (int): How checksum data should be treated. Default is 0
                    (no checking). Set to 1 to check and warn on wrong checksums; set to 2 to
                    raise an exception if any wrong checksums are found.
            recalcBBoxes (bool): If true (the default), recalculates ``glyf``, ``CFF ``,
                    ``head`` bounding box values and ``hhea``/``vhea`` min/max values on save.
                    Also compiles the glyphs on importing, which saves memory consumption and
                    time.
            ignoreDecompileErrors (bool): If true, exceptions raised during table decompilation
                    will be ignored, and the binary data will be returned for those tables instead.
            recalcTimestamp (bool): If true (the default), sets the ``modified`` timestamp in
                    the ``head`` table on save.
            fontNumber (int): The index of the font in a TrueType Collection file.
            lazy (bool): If lazy is set to True, many data structures are loaded lazily, upon
                    access only. If it is set to False, many data structures are loaded immediately.
                    The default is ``lazy=None`` which is somewhere in between.
    """

    tables: dict[Tag, DefaultTable | GlyphOrder]
    reader: SFNTReader | None
    sfntVersion: str
    flavor: str | None
    flavorData: Any | None
    lazy: bool | None
    recalcBBoxes: bool
    recalcTimestamp: bool
    ignoreDecompileErrors: bool
    cfg: AbstractConfig
    glyphOrder: list[str]
    _reverseGlyphOrderDict: dict[str, int]
    _tableCache: MutableMapping[tuple[Tag, bytes], DefaultTable] | None
    disassembleInstructions: bool
    bitmapGlyphDataFormat: str
    # Deprecated attributes
    verbose: bool | None
    quiet: bool | None

    def __init__(
        self,
        file: str | os.PathLike[str] | BinaryIO | None = None,
        res_name_or_index: str | int | None = None,
        sfntVersion: str = "\000\001\000\000",
        flavor: str | None = None,
        checkChecksums: int = 0,
        verbose: bool | None = None,  # Deprecated
        recalcBBoxes: bool = True,
        allowVID: Any = NotImplemented,  # Deprecated/Unused
        ignoreDecompileErrors: bool = False,
        recalcTimestamp: bool = True,
        fontNumber: int = -1,
        lazy: bool | None = None,
        quiet: bool | None = None,  # Deprecated
        _tableCache: MutableMapping[tuple[Tag, bytes], DefaultTable] | None = None,
        cfg: Mapping[str, Any] | AbstractConfig = {},
    ) -> None:
        # Set deprecated attributes
        for name in ("verbose", "quiet"):
            val = locals().get(name)
            if val is not None:
                deprecateArgument(name, "configure logging instead")
            setattr(self, name, val)

        self.lazy = lazy
        self.recalcBBoxes = recalcBBoxes
        self.recalcTimestamp = recalcTimestamp
        self.tables = {}
        self.reader = None
        self.cfg = cfg.copy() if isinstance(cfg, AbstractConfig) else Config(cfg)
        self.ignoreDecompileErrors = ignoreDecompileErrors

        if not file:
            self.sfntVersion = sfntVersion
            self.flavor = flavor
            self.flavorData = None
            return
        seekable = True
        if not hasattr(file, "read"):
            if not isinstance(file, (str, os.PathLike)):
                raise TypeError(
                    "fileOrPath must be a file path (str or PathLike) if it isn't an object with a `read` method."
                )
            closeStream = True
            # assume file is a string
            if res_name_or_index is not None:
                # see if it contains 'sfnt' resources in the resource or data fork
                from . import macUtils

                if res_name_or_index == 0:
                    if macUtils.getSFNTResIndices(file):
                        # get the first available sfnt font.
                        file = macUtils.SFNTResourceReader(file, 1)
                    else:
                        file = open(file, "rb")
                else:
                    file = macUtils.SFNTResourceReader(file, res_name_or_index)
            else:
                file = open(file, "rb")
        else:
            # assume "file" is a readable file object
            assert not isinstance(file, (str, os.PathLike))
            closeStream = False
            # SFNTReader wants the input file to be seekable.
            # SpooledTemporaryFile has no seekable() on < 3.11, but still can seek:
            # https://github.com/fonttools/fonttools/issues/3052
            if hasattr(file, "seekable"):
                seekable = file.seekable()
            elif hasattr(file, "seek"):
                try:
                    file.seek(0)
                except UnsupportedOperation:
                    seekable = False

        if not self.lazy:
            # read input file in memory and wrap a stream around it to allow overwriting
            if seekable:
                file.seek(0)
            tmp = BytesIO(file.read())
            if hasattr(file, "name"):
                # save reference to input file name
                tmp.name = file.name
            if closeStream:
                file.close()
            file = tmp
        elif not seekable:
            raise TTLibError("Input file must be seekable when lazy=True")
        self._tableCache = _tableCache
        self.reader = SFNTReader(file, checkChecksums, fontNumber=fontNumber)
        self.sfntVersion = self.reader.sfntVersion
        self.flavor = self.reader.flavor
        self.flavorData = self.reader.flavorData

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        """If we still have a reader object, close it."""
        if self.reader is not None:
            self.reader.close()
            self.reader = None

    def save(
        self, file: str | os.PathLike[str] | BinaryIO, reorderTables: bool | None = True
    ) -> None:
        """Save the font to disk.

        Args:
                file: Similarly to the constructor, can be either a pathname or a writable
                        binary file object.
                reorderTables (Option[bool]): If true (the default), reorder the tables,
                        sorting them by tag (recommended by the OpenType specification). If
                        false, retain the original font order. If None, reorder by table
                        dependency (fastest).
        """
        if not hasattr(file, "write"):
            if self.lazy and self.reader.file.name == file:
                raise TTLibError("Can't overwrite TTFont when 'lazy' attribute is True")
            createStream = True
        else:
            # assume "file" is a writable file object
            createStream = False

        tmp = BytesIO()

        writer_reordersTables = self._save(tmp)

        if not (
            reorderTables is None
            or writer_reordersTables
            or (reorderTables is False and self.reader is None)
        ):
            if reorderTables is False:
                # sort tables using the original font's order
                if self.reader is None:
                    raise TTLibError(
                        "The original table order is unavailable because there isn't a font to read it from."
                    )
                tableOrder = list(self.reader.keys())
            else:
                # use the recommended order from the OpenType specification
                tableOrder = None
            tmp.flush()
            tmp2 = BytesIO()
            reorderFontTables(tmp, tmp2, tableOrder)
            tmp.close()
            tmp = tmp2

        if createStream:
            # "file" is a path
            assert isinstance(file, (str, os.PathLike))
            with open(file, "wb") as file:
                file.write(tmp.getvalue())
        else:
            assert not isinstance(file, (str, os.PathLike))
            file.write(tmp.getvalue())

        tmp.close()

    def _save(
        self,
        file: BinaryIO,
        tableCache: MutableMapping[tuple[Tag, bytes], Any] | None = None,
    ) -> bool:
        """Internal function, to be shared by save() and TTCollection.save()"""

        if self.recalcTimestamp and "head" in self:
            # make sure 'head' is loaded so the recalculation is actually done
            self["head"]

        tags = self.keys()
        tags.pop(0)  # skip GlyphOrder tag
        numTables = len(tags)
        # write to a temporary stream to allow saving to unseekable streams
        writer = SFNTWriter(
            file, numTables, self.sfntVersion, self.flavor, self.flavorData
        )

        done = []
        for tag in tags:
            self._writeTable(tag, writer, done, tableCache)

        writer.close()

        return writer.reordersTables()

    class XMLSavingOptions(TypedDict):
        writeVersion: bool
        quiet: bool | None
        tables: Sequence[str | bytes] | None
        skipTables: Sequence[str] | None
        splitTables: bool
        splitGlyphs: bool
        disassembleInstructions: bool
        bitmapGlyphDataFormat: str

    def saveXML(
        self,
        fileOrPath: str | os.PathLike[str] | BinaryIO | TextIO,
        newlinestr: str = "\n",
        **kwargs: Unpack[XMLSavingOptions],
    ) -> None:
        """Export the font as TTX (an XML-based text file), or as a series of text
        files when splitTables is true. In the latter case, the 'fileOrPath'
        argument should be a path to a directory.
        The 'tables' argument must either be false (dump all tables) or a
        list of tables to dump. The 'skipTables' argument may be a list of tables
        to skip, but only when the 'tables' argument is false.
        """

        writer = xmlWriter.XMLWriter(fileOrPath, newlinestr=newlinestr)
        self._saveXML(writer, **kwargs)
        writer.close()

    def _saveXML(
        self,
        writer: xmlWriter.XMLWriter,
        writeVersion: bool = True,
        quiet: bool | None = None,  # Deprecated
        tables: Sequence[str | bytes] | None = None,
        skipTables: Sequence[str] | None = None,
        splitTables: bool = False,
        splitGlyphs: bool = False,
        disassembleInstructions: bool = True,
        bitmapGlyphDataFormat: str = "raw",
    ) -> None:
        if quiet is not None:
            deprecateArgument("quiet", "configure logging instead")

        self.disassembleInstructions = disassembleInstructions
        self.bitmapGlyphDataFormat = bitmapGlyphDataFormat
        if not tables:
            tables = self.keys()
            if skipTables:
                tables = [tag for tag in tables if tag not in skipTables]

        if writeVersion:
            from fontTools import version

            version = ".".join(version.split(".")[:2])
            writer.begintag(
                "ttFont",
                sfntVersion=repr(tostr(self.sfntVersion))[1:-1],
                ttLibVersion=version,
            )
        else:
            writer.begintag("ttFont", sfntVersion=repr(tostr(self.sfntVersion))[1:-1])
        writer.newline()

        # always splitTables if splitGlyphs is enabled
        splitTables = splitTables or splitGlyphs

        if not splitTables:
            writer.newline()
        else:
            if writer.filename is None:
                raise TTLibError(
                    "splitTables requires the file name to be a file system path, not a stream."
                )
            path, ext = os.path.splitext(writer.filename)

        for tag in tables:
            if splitTables:
                tablePath = path + "." + tagToIdentifier(tag) + ext
                tableWriter = xmlWriter.XMLWriter(
                    tablePath, newlinestr=writer.newlinestr
                )
                tableWriter.begintag("ttFont", ttLibVersion=version)
                tableWriter.newline()
                tableWriter.newline()
                writer.simpletag(tagToXML(tag), src=os.path.basename(tablePath))
                writer.newline()
            else:
                tableWriter = writer
            self._tableToXML(tableWriter, tag, splitGlyphs=splitGlyphs)
            if splitTables:
                tableWriter.endtag("ttFont")
                tableWriter.newline()
                tableWriter.close()
        writer.endtag("ttFont")
        writer.newline()

    def _tableToXML(
        self,
        writer: xmlWriter.XMLWriter,
        tag: str | bytes,
        quiet: bool | None = None,
        splitGlyphs: bool = False,
    ) -> None:
        if quiet is not None:
            deprecateArgument("quiet", "configure logging instead")
        if tag in self:
            table = self[tag]
            report = "Dumping '%s' table..." % tag
        else:
            report = "No '%s' table found." % tag
        log.info(report)
        if tag not in self:
            return
        xmlTag = tagToXML(tag)
        attrs: dict[str, Any] = {}
        if hasattr(table, "ERROR"):
            attrs["ERROR"] = "decompilation error"
        from .tables.DefaultTable import DefaultTable

        if table.__class__ == DefaultTable:
            attrs["raw"] = True
        writer.begintag(xmlTag, **attrs)
        writer.newline()
        if tag == "glyf":
            table.toXML(writer, self, splitGlyphs=splitGlyphs)
        else:
            table.toXML(writer, self)
        writer.endtag(xmlTag)
        writer.newline()
        writer.newline()

    def importXML(
        self, fileOrPath: str | os.PathLike[str] | BinaryIO, quiet: bool | None = None
    ) -> None:
        """Import a TTX file (an XML-based text format), so as to recreate
        a font object.
        """
        if quiet is not None:
            deprecateArgument("quiet", "configure logging instead")

        if "maxp" in self and "post" in self:
            # Make sure the glyph order is loaded, as it otherwise gets
            # lost if the XML doesn't contain the glyph order, yet does
            # contain the table which was originally used to extract the
            # glyph names from (ie. 'post', 'cmap' or 'CFF ').
            self.getGlyphOrder()

        from fontTools.misc import xmlReader

        reader = xmlReader.XMLReader(fileOrPath, self)
        reader.read()

    def isLoaded(self, tag: str | bytes) -> bool:
        """Return true if the table identified by ``tag`` has been
        decompiled and loaded into memory."""
        return tag in self.tables

    def has_key(self, tag: str | bytes) -> bool:
        """Test if the table identified by ``tag`` is present in the font.

        As well as this method, ``tag in font`` can also be used to determine the
        presence of the table."""
        if self.isLoaded(tag):
            return True
        elif self.reader and tag in self.reader:
            return True
        elif tag == "GlyphOrder":
            return True
        else:
            return False

    __contains__ = has_key

    def keys(self) -> list[str]:
        """Returns the list of tables in the font, along with the ``GlyphOrder`` pseudo-table."""
        keys = list(self.tables.keys())
        if self.reader:
            for key in list(self.reader.keys()):
                if key not in keys:
                    keys.append(key)

        if "GlyphOrder" in keys:
            keys.remove("GlyphOrder")
        keys = sortedTagList(keys)
        return ["GlyphOrder"] + keys

    def ensureDecompiled(self, recurse: bool | None = None) -> None:
        """Decompile all the tables, even if a TTFont was opened in 'lazy' mode."""
        for tag in self.keys():
            table = self[tag]
            if recurse is None:
                recurse = self.lazy is not False
            if recurse and hasattr(table, "ensureDecompiled"):
                table.ensureDecompiled(recurse=recurse)
        self.lazy = False

    def __len__(self) -> int:
        return len(list(self.keys()))

    @overload
    def __getitem__(self, tag: Literal["BASE"]) -> B_A_S_E_.table_B_A_S_E_: ...
    @overload
    def __getitem__(self, tag: Literal["CBDT"]) -> C_B_D_T_.table_C_B_D_T_: ...
    @overload
    def __getitem__(self, tag: Literal["CBLC"]) -> C_B_L_C_.table_C_B_L_C_: ...
    @overload
    def __getitem__(self, tag: Literal["CFF "]) -> C_F_F_.table_C_F_F_: ...
    @overload
    def __getitem__(self, tag: Literal["CFF2"]) -> C_F_F__2.table_C_F_F__2: ...
    @overload
    def __getitem__(self, tag: Literal["COLR"]) -> C_O_L_R_.table_C_O_L_R_: ...
    @overload
    def __getitem__(self, tag: Literal["CPAL"]) -> C_P_A_L_.table_C_P_A_L_: ...
    @overload
    def __getitem__(self, tag: Literal["DSIG"]) -> D_S_I_G_.table_D_S_I_G_: ...
    @overload
    def __getitem__(self, tag: Literal["EBDT"]) -> E_B_D_T_.table_E_B_D_T_: ...
    @overload
    def __getitem__(self, tag: Literal["EBLC"]) -> E_B_L_C_.table_E_B_L_C_: ...
    @overload
    def __getitem__(self, tag: Literal["FFTM"]) -> F_F_T_M_.table_F_F_T_M_: ...
    @overload
    def __getitem__(self, tag: Literal["GDEF"]) -> G_D_E_F_.table_G_D_E_F_: ...
    @overload
    def __getitem__(self, tag: Literal["GMAP"]) -> G_M_A_P_.table_G_M_A_P_: ...
    @overload
    def __getitem__(self, tag: Literal["GPKG"]) -> G_P_K_G_.table_G_P_K_G_: ...
    @overload
    def __getitem__(self, tag: Literal["GPOS"]) -> G_P_O_S_.table_G_P_O_S_: ...
    @overload
    def __getitem__(self, tag: Literal["GSUB"]) -> G_S_U_B_.table_G_S_U_B_: ...
    @overload
    def __getitem__(self, tag: Literal["GVAR"]) -> G_V_A_R_.table_G_V_A_R_: ...
    @overload
    def __getitem__(self, tag: Literal["HVAR"]) -> H_V_A_R_.table_H_V_A_R_: ...
    @overload
    def __getitem__(self, tag: Literal["JSTF"]) -> J_S_T_F_.table_J_S_T_F_: ...
    @overload
    def __getitem__(self, tag: Literal["LTSH"]) -> L_T_S_H_.table_L_T_S_H_: ...
    @overload
    def __getitem__(self, tag: Literal["MATH"]) -> M_A_T_H_.table_M_A_T_H_: ...
    @overload
    def __getitem__(self, tag: Literal["META"]) -> M_E_T_A_.table_M_E_T_A_: ...
    @overload
    def __getitem__(self, tag: Literal["MVAR"]) -> M_V_A_R_.table_M_V_A_R_: ...
    @overload
    def __getitem__(self, tag: Literal["SING"]) -> S_I_N_G_.table_S_I_N_G_: ...
    @overload
    def __getitem__(self, tag: Literal["STAT"]) -> S_T_A_T_.table_S_T_A_T_: ...
    @overload
    def __getitem__(self, tag: Literal["SVG "]) -> S_V_G_.table_S_V_G_: ...
    @overload
    def __getitem__(self, tag: Literal["TSI0"]) -> T_S_I__0.table_T_S_I__0: ...
    @overload
    def __getitem__(self, tag: Literal["TSI1"]) -> T_S_I__1.table_T_S_I__1: ...
    @overload
    def __getitem__(self, tag: Literal["TSI2"]) -> T_S_I__2.table_T_S_I__2: ...
    @overload
    def __getitem__(self, tag: Literal["TSI3"]) -> T_S_I__3.table_T_S_I__3: ...
    @overload
    def __getitem__(self, tag: Literal["TSI5"]) -> T_S_I__5.table_T_S_I__5: ...
    @overload
    def __getitem__(self, tag: Literal["TSIB"]) -> T_S_I_B_.table_T_S_I_B_: ...
    @overload
    def __getitem__(self, tag: Literal["TSIC"]) -> T_S_I_C_.table_T_S_I_C_: ...
    @overload
    def __getitem__(self, tag: Literal["TSID"]) -> T_S_I_D_.table_T_S_I_D_: ...
    @overload
    def __getitem__(self, tag: Literal["TSIJ"]) -> T_S_I_J_.table_T_S_I_J_: ...
    @overload
    def __getitem__(self, tag: Literal["TSIP"]) -> T_S_I_P_.table_T_S_I_P_: ...
    @overload
    def __getitem__(self, tag: Literal["TSIS"]) -> T_S_I_S_.table_T_S_I_S_: ...
    @overload
    def __getitem__(self, tag: Literal["TSIV"]) -> T_S_I_V_.table_T_S_I_V_: ...
    @overload
    def __getitem__(self, tag: Literal["TTFA"]) -> T_T_F_A_.table_T_T_F_A_: ...
    @overload
    def __getitem__(self, tag: Literal["VARC"]) -> V_A_R_C_.table_V_A_R_C_: ...
    @overload
    def __getitem__(self, tag: Literal["VDMX"]) -> V_D_M_X_.table_V_D_M_X_: ...
    @overload
    def __getitem__(self, tag: Literal["VORG"]) -> V_O_R_G_.table_V_O_R_G_: ...
    @overload
    def __getitem__(self, tag: Literal["VVAR"]) -> V_V_A_R_.table_V_V_A_R_: ...
    @overload
    def __getitem__(self, tag: Literal["Debg"]) -> D__e_b_g.table_D__e_b_g: ...
    @overload
    def __getitem__(self, tag: Literal["Feat"]) -> F__e_a_t.table_F__e_a_t: ...
    @overload
    def __getitem__(self, tag: Literal["Glat"]) -> G__l_a_t.table_G__l_a_t: ...
    @overload
    def __getitem__(self, tag: Literal["Gloc"]) -> G__l_o_c.table_G__l_o_c: ...
    @overload
    def __getitem__(self, tag: Literal["OS/2"]) -> O_S_2f_2.table_O_S_2f_2: ...
    @overload
    def __getitem__(self, tag: Literal["Silf"]) -> S__i_l_f.table_S__i_l_f: ...
    @overload
    def __getitem__(self, tag: Literal["Sill"]) -> S__i_l_l.table_S__i_l_l: ...
    @overload
    def __getitem__(self, tag: Literal["ankr"]) -> _a_n_k_r.table__a_n_k_r: ...
    @overload
    def __getitem__(self, tag: Literal["avar"]) -> _a_v_a_r.table__a_v_a_r: ...
    @overload
    def __getitem__(self, tag: Literal["bsln"]) -> _b_s_l_n.table__b_s_l_n: ...
    @overload
    def __getitem__(self, tag: Literal["cidg"]) -> _c_i_d_g.table__c_i_d_g: ...
    @overload
    def __getitem__(self, tag: Literal["cmap"]) -> _c_m_a_p.table__c_m_a_p: ...
    @overload
    def __getitem__(self, tag: Literal["cvar"]) -> _c_v_a_r.table__c_v_a_r: ...
    @overload
    def __getitem__(self, tag: Literal["cvt "]) -> _c_v_t.table__c_v_t: ...
    @overload
    def __getitem__(self, tag: Literal["feat"]) -> _f_e_a_t.table__f_e_a_t: ...
    @overload
    def __getitem__(self, tag: Literal["fpgm"]) -> _f_p_g_m.table__f_p_g_m: ...
    @overload
    def __getitem__(self, tag: Literal["fvar"]) -> _f_v_a_r.table__f_v_a_r: ...
    @overload
    def __getitem__(self, tag: Literal["gasp"]) -> _g_a_s_p.table__g_a_s_p: ...
    @overload
    def __getitem__(self, tag: Literal["gcid"]) -> _g_c_i_d.table__g_c_i_d: ...
    @overload
    def __getitem__(self, tag: Literal["glyf"]) -> _g_l_y_f.table__g_l_y_f: ...
    @overload
    def __getitem__(self, tag: Literal["gvar"]) -> _g_v_a_r.table__g_v_a_r: ...
    @overload
    def __getitem__(self, tag: Literal["hdmx"]) -> _h_d_m_x.table__h_d_m_x: ...
    @overload
    def __getitem__(self, tag: Literal["head"]) -> _h_e_a_d.table__h_e_a_d: ...
    @overload
    def __getitem__(self, tag: Literal["hhea"]) -> _h_h_e_a.table__h_h_e_a: ...
    @overload
    def __getitem__(self, tag: Literal["hmtx"]) -> _h_m_t_x.table__h_m_t_x: ...
    @overload
    def __getitem__(self, tag: Literal["kern"]) -> _k_e_r_n.table__k_e_r_n: ...
    @overload
    def __getitem__(self, tag: Literal["lcar"]) -> _l_c_a_r.table__l_c_a_r: ...
    @overload
    def __getitem__(self, tag: Literal["loca"]) -> _l_o_c_a.table__l_o_c_a: ...
    @overload
    def __getitem__(self, tag: Literal["ltag"]) -> _l_t_a_g.table__l_t_a_g: ...
    @overload
    def __getitem__(self, tag: Literal["maxp"]) -> _m_a_x_p.table__m_a_x_p: ...
    @overload
    def __getitem__(self, tag: Literal["meta"]) -> _m_e_t_a.table__m_e_t_a: ...
    @overload
    def __getitem__(self, tag: Literal["mort"]) -> _m_o_r_t.table__m_o_r_t: ...
    @overload
    def __getitem__(self, tag: Literal["morx"]) -> _m_o_r_x.table__m_o_r_x: ...
    @overload
    def __getitem__(self, tag: Literal["name"]) -> _n_a_m_e.table__n_a_m_e: ...
    @overload
    def __getitem__(self, tag: Literal["opbd"]) -> _o_p_b_d.table__o_p_b_d: ...
    @overload
    def __getitem__(self, tag: Literal["post"]) -> _p_o_s_t.table__p_o_s_t: ...
    @overload
    def __getitem__(self, tag: Literal["prep"]) -> _p_r_e_p.table__p_r_e_p: ...
    @overload
    def __getitem__(self, tag: Literal["prop"]) -> _p_r_o_p.table__p_r_o_p: ...
    @overload
    def __getitem__(self, tag: Literal["sbix"]) -> _s_b_i_x.table__s_b_i_x: ...
    @overload
    def __getitem__(self, tag: Literal["trak"]) -> _t_r_a_k.table__t_r_a_k: ...
    @overload
    def __getitem__(self, tag: Literal["vhea"]) -> _v_h_e_a.table__v_h_e_a: ...
    @overload
    def __getitem__(self, tag: Literal["vmtx"]) -> _v_m_t_x.table__v_m_t_x: ...
    @overload
    def __getitem__(self, tag: Literal["GlyphOrder"]) -> GlyphOrder: ...
    @overload
    def __getitem__(self, tag: str | bytes) -> DefaultTable | GlyphOrder: ...

    def __getitem__(self, tag: str | bytes) -> DefaultTable | GlyphOrder:
        tag = Tag(tag)
        table = self.tables.get(tag)
        if table is None:
            if tag == "GlyphOrder":
                table = GlyphOrder(tag)
                self.tables[tag] = table
            elif self.reader is not None:
                table = self._readTable(tag)
            else:
                raise KeyError("'%s' table not found" % tag)
        return table

    def _readTable(self, tag: Tag) -> DefaultTable:
        log.debug("Reading '%s' table from disk", tag)
        assert self.reader is not None
        data = self.reader[tag]
        if self._tableCache is not None:
            table = self._tableCache.get((tag, data))
            if table is not None:
                return table
        tableClass = getTableClass(tag)
        table = tableClass(tag)
        self.tables[tag] = table
        log.debug("Decompiling '%s' table", tag)
        try:
            table.decompile(data, self)
        except Exception:
            if not self.ignoreDecompileErrors:
                raise
            # fall back to DefaultTable, retaining the binary table data
            log.exception(
                "An exception occurred during the decompilation of the '%s' table", tag
            )
            from .tables.DefaultTable import DefaultTable

            file = StringIO()
            traceback.print_exc(file=file)
            table = DefaultTable(tag)
            table.ERROR = file.getvalue()
            self.tables[tag] = table
            table.decompile(data, self)
        if self._tableCache is not None:
            self._tableCache[(tag, data)] = table
        return table

    def __setitem__(self, tag: str | bytes, table: DefaultTable) -> None:
        self.tables[Tag(tag)] = table

    def __delitem__(self, tag: str | bytes) -> None:
        if tag not in self:
            raise KeyError("'%s' table not found" % tag)
        if tag in self.tables:
            del self.tables[tag]
        if self.reader and tag in self.reader:
            del self.reader[tag]

    @overload
    def get(self, tag: Literal["BASE"]) -> B_A_S_E_.table_B_A_S_E_ | None: ...
    @overload
    def get(self, tag: Literal["CBDT"]) -> C_B_D_T_.table_C_B_D_T_ | None: ...
    @overload
    def get(self, tag: Literal["CBLC"]) -> C_B_L_C_.table_C_B_L_C_ | None: ...
    @overload
    def get(self, tag: Literal["CFF "]) -> C_F_F_.table_C_F_F_ | None: ...
    @overload
    def get(self, tag: Literal["CFF2"]) -> C_F_F__2.table_C_F_F__2 | None: ...
    @overload
    def get(self, tag: Literal["COLR"]) -> C_O_L_R_.table_C_O_L_R_ | None: ...
    @overload
    def get(self, tag: Literal["CPAL"]) -> C_P_A_L_.table_C_P_A_L_ | None: ...
    @overload
    def get(self, tag: Literal["DSIG"]) -> D_S_I_G_.table_D_S_I_G_ | None: ...
    @overload
    def get(self, tag: Literal["EBDT"]) -> E_B_D_T_.table_E_B_D_T_ | None: ...
    @overload
    def get(self, tag: Literal["EBLC"]) -> E_B_L_C_.table_E_B_L_C_ | None: ...
    @overload
    def get(self, tag: Literal["FFTM"]) -> F_F_T_M_.table_F_F_T_M_ | None: ...
    @overload
    def get(self, tag: Literal["GDEF"]) -> G_D_E_F_.table_G_D_E_F_ | None: ...
    @overload
    def get(self, tag: Literal["GMAP"]) -> G_M_A_P_.table_G_M_A_P_ | None: ...
    @overload
    def get(self, tag: Literal["GPKG"]) -> G_P_K_G_.table_G_P_K_G_ | None: ...
    @overload
    def get(self, tag: Literal["GPOS"]) -> G_P_O_S_.table_G_P_O_S_ | None: ...
    @overload
    def get(self, tag: Literal["GSUB"]) -> G_S_U_B_.table_G_S_U_B_ | None: ...
    @overload
    def get(self, tag: Literal["GVAR"]) -> G_V_A_R_.table_G_V_A_R_ | None: ...
    @overload
    def get(self, tag: Literal["HVAR"]) -> H_V_A_R_.table_H_V_A_R_ | None: ...
    @overload
    def get(self, tag: Literal["JSTF"]) -> J_S_T_F_.table_J_S_T_F_ | None: ...
    @overload
    def get(self, tag: Literal["LTSH"]) -> L_T_S_H_.table_L_T_S_H_ | None: ...
    @overload
    def get(self, tag: Literal["MATH"]) -> M_A_T_H_.table_M_A_T_H_ | None: ...
    @overload
    def get(self, tag: Literal["META"]) -> M_E_T_A_.table_M_E_T_A_ | None: ...
    @overload
    def get(self, tag: Literal["MVAR"]) -> M_V_A_R_.table_M_V_A_R_ | None: ...
    @overload
    def get(self, tag: Literal["SING"]) -> S_I_N_G_.table_S_I_N_G_ | None: ...
    @overload
    def get(self, tag: Literal["STAT"]) -> S_T_A_T_.table_S_T_A_T_ | None: ...
    @overload
    def get(self, tag: Literal["SVG "]) -> S_V_G_.table_S_V_G_ | None: ...
    @overload
    def get(self, tag: Literal["TSI0"]) -> T_S_I__0.table_T_S_I__0 | None: ...
    @overload
    def get(self, tag: Literal["TSI1"]) -> T_S_I__1.table_T_S_I__1 | None: ...
    @overload
    def get(self, tag: Literal["TSI2"]) -> T_S_I__2.table_T_S_I__2 | None: ...
    @overload
    def get(self, tag: Literal["TSI3"]) -> T_S_I__3.table_T_S_I__3 | None: ...
    @overload
    def get(self, tag: Literal["TSI5"]) -> T_S_I__5.table_T_S_I__5 | None: ...
    @overload
    def get(self, tag: Literal["TSIB"]) -> T_S_I_B_.table_T_S_I_B_ | None: ...
    @overload
    def get(self, tag: Literal["TSIC"]) -> T_S_I_C_.table_T_S_I_C_ | None: ...
    @overload
    def get(self, tag: Literal["TSID"]) -> T_S_I_D_.table_T_S_I_D_ | None: ...
    @overload
    def get(self, tag: Literal["TSIJ"]) -> T_S_I_J_.table_T_S_I_J_ | None: ...
    @overload
    def get(self, tag: Literal["TSIP"]) -> T_S_I_P_.table_T_S_I_P_ | None: ...
    @overload
    def get(self, tag: Literal["TSIS"]) -> T_S_I_S_.table_T_S_I_S_ | None: ...
    @overload
    def get(self, tag: Literal["TSIV"]) -> T_S_I_V_.table_T_S_I_V_ | None: ...
    @overload
    def get(self, tag: Literal["TTFA"]) -> T_T_F_A_.table_T_T_F_A_ | None: ...
    @overload
    def get(self, tag: Literal["VARC"]) -> V_A_R_C_.table_V_A_R_C_ | None: ...
    @overload
    def get(self, tag: Literal["VDMX"]) -> V_D_M_X_.table_V_D_M_X_ | None: ...
    @overload
    def get(self, tag: Literal["VORG"]) -> V_O_R_G_.table_V_O_R_G_ | None: ...
    @overload
    def get(self, tag: Literal["VVAR"]) -> V_V_A_R_.table_V_V_A_R_ | None: ...
    @overload
    def get(self, tag: Literal["Debg"]) -> D__e_b_g.table_D__e_b_g | None: ...
    @overload
    def get(self, tag: Literal["Feat"]) -> F__e_a_t.table_F__e_a_t | None: ...
    @overload
    def get(self, tag: Literal["Glat"]) -> G__l_a_t.table_G__l_a_t | None: ...
    @overload
    def get(self, tag: Literal["Gloc"]) -> G__l_o_c.table_G__l_o_c | None: ...
    @overload
    def get(self, tag: Literal["OS/2"]) -> O_S_2f_2.table_O_S_2f_2 | None: ...
    @overload
    def get(self, tag: Literal["Silf"]) -> S__i_l_f.table_S__i_l_f | None: ...
    @overload
    def get(self, tag: Literal["Sill"]) -> S__i_l_l.table_S__i_l_l | None: ...
    @overload
    def get(self, tag: Literal["ankr"]) -> _a_n_k_r.table__a_n_k_r | None: ...
    @overload
    def get(self, tag: Literal["avar"]) -> _a_v_a_r.table__a_v_a_r | None: ...
    @overload
    def get(self, tag: Literal["bsln"]) -> _b_s_l_n.table__b_s_l_n | None: ...
    @overload
    def get(self, tag: Literal["cidg"]) -> _c_i_d_g.table__c_i_d_g | None: ...
    @overload
    def get(self, tag: Literal["cmap"]) -> _c_m_a_p.table__c_m_a_p | None: ...
    @overload
    def get(self, tag: Literal["cvar"]) -> _c_v_a_r.table__c_v_a_r | None: ...
    @overload
    def get(self, tag: Literal["cvt "]) -> _c_v_t.table__c_v_t | None: ...
    @overload
    def get(self, tag: Literal["feat"]) -> _f_e_a_t.table__f_e_a_t | None: ...
    @overload
    def get(self, tag: Literal["fpgm"]) -> _f_p_g_m.table__f_p_g_m | None: ...
    @overload
    def get(self, tag: Literal["fvar"]) -> _f_v_a_r.table__f_v_a_r | None: ...
    @overload
    def get(self, tag: Literal["gasp"]) -> _g_a_s_p.table__g_a_s_p | None: ...
    @overload
    def get(self, tag: Literal["gcid"]) -> _g_c_i_d.table__g_c_i_d | None: ...
    @overload
    def get(self, tag: Literal["glyf"]) -> _g_l_y_f.table__g_l_y_f | None: ...
    @overload
    def get(self, tag: Literal["gvar"]) -> _g_v_a_r.table__g_v_a_r | None: ...
    @overload
    def get(self, tag: Literal["hdmx"]) -> _h_d_m_x.table__h_d_m_x | None: ...
    @overload
    def get(self, tag: Literal["head"]) -> _h_e_a_d.table__h_e_a_d | None: ...
    @overload
    def get(self, tag: Literal["hhea"]) -> _h_h_e_a.table__h_h_e_a | None: ...
    @overload
    def get(self, tag: Literal["hmtx"]) -> _h_m_t_x.table__h_m_t_x | None: ...
    @overload
    def get(self, tag: Literal["kern"]) -> _k_e_r_n.table__k_e_r_n | None: ...
    @overload
    def get(self, tag: Literal["lcar"]) -> _l_c_a_r.table__l_c_a_r | None: ...
    @overload
    def get(self, tag: Literal["loca"]) -> _l_o_c_a.table__l_o_c_a | None: ...
    @overload
    def get(self, tag: Literal["ltag"]) -> _l_t_a_g.table__l_t_a_g | None: ...
    @overload
    def get(self, tag: Literal["maxp"]) -> _m_a_x_p.table__m_a_x_p | None: ...
    @overload
    def get(self, tag: Literal["meta"]) -> _m_e_t_a.table__m_e_t_a | None: ...
    @overload
    def get(self, tag: Literal["mort"]) -> _m_o_r_t.table__m_o_r_t | None: ...
    @overload
    def get(self, tag: Literal["morx"]) -> _m_o_r_x.table__m_o_r_x | None: ...
    @overload
    def get(self, tag: Literal["name"]) -> _n_a_m_e.table__n_a_m_e | None: ...
    @overload
    def get(self, tag: Literal["opbd"]) -> _o_p_b_d.table__o_p_b_d | None: ...
    @overload
    def get(self, tag: Literal["post"]) -> _p_o_s_t.table__p_o_s_t | None: ...
    @overload
    def get(self, tag: Literal["prep"]) -> _p_r_e_p.table__p_r_e_p | None: ...
    @overload
    def get(self, tag: Literal["prop"]) -> _p_r_o_p.table__p_r_o_p | None: ...
    @overload
    def get(self, tag: Literal["sbix"]) -> _s_b_i_x.table__s_b_i_x | None: ...
    @overload
    def get(self, tag: Literal["trak"]) -> _t_r_a_k.table__t_r_a_k | None: ...
    @overload
    def get(self, tag: Literal["vhea"]) -> _v_h_e_a.table__v_h_e_a | None: ...
    @overload
    def get(self, tag: Literal["vmtx"]) -> _v_m_t_x.table__v_m_t_x | None: ...
    @overload
    def get(self, tag: Literal["GlyphOrder"]) -> GlyphOrder: ...
    @overload
    def get(self, tag: str | bytes) -> DefaultTable | GlyphOrder | Any | None: ...
    @overload
    def get(
        self, tag: str | bytes, default: _VT_co
    ) -> DefaultTable | GlyphOrder | Any | _VT_co: ...

    def get(
        self, tag: str | bytes, default: Any | None = None
    ) -> DefaultTable | GlyphOrder | Any | None:
        """Returns the table if it exists or (optionally) a default if it doesn't."""
        try:
            return self[tag]
        except KeyError:
            return default

    def setGlyphOrder(self, glyphOrder: list[str]) -> None:
        """Set the glyph order

        Args:
                glyphOrder ([str]): List of glyph names in order.
        """
        self.glyphOrder = glyphOrder
        if hasattr(self, "_reverseGlyphOrderDict"):
            del self._reverseGlyphOrderDict
        if self.isLoaded("glyf"):
            self["glyf"].setGlyphOrder(glyphOrder)

    def getGlyphOrder(self) -> list[str]:
        """Returns a list of glyph names ordered by their position in the font."""
        try:
            return self.glyphOrder
        except AttributeError:
            pass
        if "CFF " in self:
            cff = self["CFF "]
            self.glyphOrder = cff.getGlyphOrder()
        elif "post" in self:
            # TrueType font
            glyphOrder = self["post"].getGlyphOrder()
            if glyphOrder is None:
                #
                # No names found in the 'post' table.
                # Try to create glyph names from the unicode cmap (if available)
                # in combination with the Adobe Glyph List (AGL).
                #
                self._getGlyphNamesFromCmap()
            elif len(glyphOrder) < self["maxp"].numGlyphs:
                #
                # Not enough names found in the 'post' table.
                # Can happen when 'post' format 1 is improperly used on a font that
                # has more than 258 glyphs (the length of 'standardGlyphOrder').
                #
                log.warning(
                    "Not enough names found in the 'post' table, generating them from cmap instead"
                )
                self._getGlyphNamesFromCmap()
            else:
                self.glyphOrder = glyphOrder
        else:
            self._getGlyphNamesFromCmap()
        return self.glyphOrder

    def _getGlyphNamesFromCmap(self) -> None:
        #
        # This is rather convoluted, but then again, it's an interesting problem:
        # - we need to use the unicode values found in the cmap table to
        #   build glyph names (eg. because there is only a minimal post table,
        #   or none at all).
        # - but the cmap parser also needs glyph names to work with...
        # So here's what we do:
        # - make up glyph names based on glyphID
        # - load a temporary cmap table based on those names
        # - extract the unicode values, build the "real" glyph names
        # - unload the temporary cmap table
        #
        if self.isLoaded("cmap"):
            # Bootstrapping: we're getting called by the cmap parser
            # itself. This means self.tables['cmap'] contains a partially
            # loaded cmap, making it impossible to get at a unicode
            # subtable here. We remove the partially loaded cmap and
            # restore it later.
            # This only happens if the cmap table is loaded before any
            # other table that does f.getGlyphOrder()  or f.getGlyphName().
            cmapLoading = self.tables["cmap"]
            del self.tables["cmap"]
        else:
            cmapLoading = None
        # Make up glyph names based on glyphID, which will be used by the
        # temporary cmap and by the real cmap in case we don't find a unicode
        # cmap.
        numGlyphs = int(self["maxp"].numGlyphs)
        glyphOrder = ["glyph%.5d" % i for i in range(numGlyphs)]
        glyphOrder[0] = ".notdef"
        # Set the glyph order, so the cmap parser has something
        # to work with (so we don't get called recursively).
        self.glyphOrder = glyphOrder

        # Make up glyph names based on the reversed cmap table. Because some
        # glyphs (eg. ligatures or alternates) may not be reachable via cmap,
        # this naming table will usually not cover all glyphs in the font.
        # If the font has no Unicode cmap table, reversecmap will be empty.
        if "cmap" in self:
            reversecmap = self["cmap"].buildReversedMin()
        else:
            reversecmap = {}
        useCount = {}
        for i, tempName in enumerate(glyphOrder):
            if tempName in reversecmap:
                # If a font maps both U+0041 LATIN CAPITAL LETTER A and
                # U+0391 GREEK CAPITAL LETTER ALPHA to the same glyph,
                # we prefer naming the glyph as "A".
                glyphName = self._makeGlyphName(reversecmap[tempName])
                numUses = useCount[glyphName] = useCount.get(glyphName, 0) + 1
                if numUses > 1:
                    glyphName = "%s.alt%d" % (glyphName, numUses - 1)
                glyphOrder[i] = glyphName

        if "cmap" in self:
            # Delete the temporary cmap table from the cache, so it can
            # be parsed again with the right names.
            del self.tables["cmap"]
            self.glyphOrder = glyphOrder
            if cmapLoading:
                # restore partially loaded cmap, so it can continue loading
                # using the proper names.
                self.tables["cmap"] = cmapLoading

    @staticmethod
    def _makeGlyphName(codepoint: int) -> str:
        from fontTools import agl  # Adobe Glyph List

        if codepoint in agl.UV2AGL:
            return agl.UV2AGL[codepoint]
        elif codepoint <= 0xFFFF:
            return "uni%04X" % codepoint
        else:
            return "u%X" % codepoint

    def getGlyphNames(self) -> list[str]:
        """Get a list of glyph names, sorted alphabetically."""
        glyphNames = sorted(self.getGlyphOrder())
        return glyphNames

    def getGlyphNames2(self) -> list[str]:
        """Get a list of glyph names, sorted alphabetically,
        but not case sensitive.
        """
        from fontTools.misc import textTools

        return textTools.caselessSort(self.getGlyphOrder())

    def getGlyphName(self, glyphID: int) -> str:
        """Returns the name for the glyph with the given ID.

        If no name is available, synthesises one with the form ``glyphXXXXX``` where
        ```XXXXX`` is the zero-padded glyph ID.
        """
        try:
            return self.getGlyphOrder()[glyphID]
        except IndexError:
            return "glyph%.5d" % glyphID

    def getGlyphNameMany(self, lst: Sequence[int]) -> list[str]:
        """Converts a list of glyph IDs into a list of glyph names."""
        glyphOrder = self.getGlyphOrder()
        cnt = len(glyphOrder)
        return [glyphOrder[gid] if gid < cnt else "glyph%.5d" % gid for gid in lst]

    def getGlyphID(self, glyphName: str) -> int:
        """Returns the ID of the glyph with the given name."""
        try:
            return self.getReverseGlyphMap()[glyphName]
        except KeyError:
            if glyphName[:5] == "glyph":
                try:
                    return int(glyphName[5:])
                except (NameError, ValueError):
                    raise KeyError(glyphName)
            raise

    def getGlyphIDMany(self, lst: Sequence[str]) -> list[int]:
        """Converts a list of glyph names into a list of glyph IDs."""
        d = self.getReverseGlyphMap()
        try:
            return [d[glyphName] for glyphName in lst]
        except KeyError:
            getGlyphID = self.getGlyphID
            return [getGlyphID(glyphName) for glyphName in lst]

    def getReverseGlyphMap(self, rebuild: bool = False) -> dict[str, int]:
        """Returns a mapping of glyph names to glyph IDs."""
        if rebuild or not hasattr(self, "_reverseGlyphOrderDict"):
            self._buildReverseGlyphOrderDict()
        return self._reverseGlyphOrderDict

    def _buildReverseGlyphOrderDict(self) -> dict[str, int]:
        self._reverseGlyphOrderDict = d = {}
        for glyphID, glyphName in enumerate(self.getGlyphOrder()):
            d[glyphName] = glyphID
        return d

    def _writeTable(
        self,
        tag: str | bytes,
        writer: SFNTWriter,
        done: list[str | bytes],  # Use list as original
        tableCache: MutableMapping[tuple[Tag, bytes], DefaultTable] | None = None,
    ) -> None:
        """Internal helper function for self.save(). Keeps track of
        inter-table dependencies.
        """
        if tag in done:
            return
        tableClass = getTableClass(tag)
        for masterTable in tableClass.dependencies:
            if masterTable not in done:
                if masterTable in self:
                    self._writeTable(masterTable, writer, done, tableCache)
                else:
                    done.append(masterTable)
        done.append(tag)
        tabledata = self.getTableData(tag)
        if tableCache is not None:
            entry = tableCache.get((Tag(tag), tabledata))
            if entry is not None:
                log.debug("reusing '%s' table", tag)
                writer.setEntry(tag, entry)
                return
        log.debug("Writing '%s' table to disk", tag)
        writer[tag] = tabledata
        if tableCache is not None:
            tableCache[(Tag(tag), tabledata)] = writer[tag]

    def getTableData(self, tag: str | bytes) -> bytes:
        """Returns the binary representation of a table.

        If the table is currently loaded and in memory, the data is compiled to
        binary and returned; if it is not currently loaded, the binary data is
        read from the font file and returned.
        """
        tag = Tag(tag)
        if self.isLoaded(tag):
            log.debug("Compiling '%s' table", tag)
            return self.tables[tag].compile(self)
        elif self.reader and tag in self.reader:
            log.debug("Reading '%s' table from disk", tag)
            return self.reader[tag]
        else:
            raise KeyError(tag)

    def getGlyphSet(
        self,
        preferCFF: bool = True,
        location: Mapping[str, _NumberT] | None = None,
        normalized: bool = False,
        recalcBounds: bool = True,
    ) -> _TTGlyphSet:
        """Return a generic GlyphSet, which is a dict-like object
        mapping glyph names to glyph objects. The returned glyph objects
        have a ``.draw()`` method that supports the Pen protocol, and will
        have an attribute named 'width'.

        If the font is CFF-based, the outlines will be taken from the ``CFF ``
        or ``CFF2`` tables. Otherwise the outlines will be taken from the
        ``glyf`` table.

        If the font contains both a ``CFF ``/``CFF2`` and a ``glyf`` table, you
        can use the ``preferCFF`` argument to specify which one should be taken.
        If the font contains both a ``CFF `` and a ``CFF2`` table, the latter is
        taken.

        If the ``location`` parameter is set, it should be a dictionary mapping
        four-letter variation tags to their float values, and the returned
        glyph-set will represent an instance of a variable font at that
        location.

        If the ``normalized`` variable is set to True, that location is
        interpreted as in the normalized (-1..+1) space, otherwise it is in the
        font's defined axes space.
        """
        if location and "fvar" not in self:
            location = None
        if location and not normalized:
            location = self.normalizeLocation(location)
        glyphSet = None
        if ("CFF " in self or "CFF2" in self) and (preferCFF or "glyf" not in self):
            glyphSet = _TTGlyphSetCFF(self, location)
        elif "glyf" in self:
            glyphSet = _TTGlyphSetGlyf(self, location, recalcBounds=recalcBounds)
        else:
            raise TTLibError("Font contains no outlines")
        if "VARC" in self:
            glyphSet = _TTGlyphSetVARC(self, location, glyphSet)
        return glyphSet

    def normalizeLocation(self, location: Mapping[str, float]) -> dict[str, float]:
        """Normalize a ``location`` from the font's defined axes space (also
        known as user space) into the normalized (-1..+1) space. It applies
        ``avar`` mapping if the font contains an ``avar`` table.

        The ``location`` parameter should be a dictionary mapping four-letter
        variation tags to their float values.

        Raises ``TTLibError`` if the font is not a variable font.
        """
        from fontTools.varLib.models import normalizeLocation

        if "fvar" not in self:
            raise TTLibError("Not a variable font")

        axes = self["fvar"].getAxes()
        location = normalizeLocation(location, axes)
        if "avar" in self:
            location = self["avar"].renormalizeLocation(location, self)
        return location

    def getBestCmap(
        self,
        cmapPreferences: Sequence[tuple[int, int]] = (
            (3, 10),
            (0, 6),
            (0, 4),
            (3, 1),
            (0, 3),
            (0, 2),
            (0, 1),
            (0, 0),
        ),
    ) -> dict[int, str] | None:
        """Returns the 'best' Unicode cmap dictionary available in the font
        or ``None``, if no Unicode cmap subtable is available.

        By default it will search for the following (platformID, platEncID)
        pairs in order::

                        (3, 10), # Windows Unicode full repertoire
                        (0, 6),  # Unicode full repertoire (format 13 subtable)
                        (0, 4),  # Unicode 2.0 full repertoire
                        (3, 1),  # Windows Unicode BMP
                        (0, 3),  # Unicode 2.0 BMP
                        (0, 2),  # Unicode ISO/IEC 10646
                        (0, 1),  # Unicode 1.1
                        (0, 0)   # Unicode 1.0

        This particular order matches what HarfBuzz uses to choose what
        subtable to use by default. This order prefers the largest-repertoire
        subtable, and among those, prefers the Windows-platform over the
        Unicode-platform as the former has wider support.

        This order can be customized via the ``cmapPreferences`` argument.
        """
        return self["cmap"].getBestCmap(cmapPreferences=cmapPreferences)

    def reorderGlyphs(self, new_glyph_order: list[str]) -> None:
        from .reorderGlyphs import reorderGlyphs

        reorderGlyphs(self, new_glyph_order)


class GlyphOrder(object):
    """A pseudo table. The glyph order isn't in the font as a separate
    table, but it's nice to present it as such in the TTX format.
    """

    def __init__(self, tag: str | None = None) -> None:
        pass

    def toXML(self, writer: xmlWriter.XMLWriter, ttFont: TTFont) -> None:
        glyphOrder = ttFont.getGlyphOrder()
        writer.comment(
            "The 'id' attribute is only for humans; it is ignored when parsed."
        )
        writer.newline()
        for i, glyphName in enumerate(glyphOrder):
            writer.simpletag("GlyphID", id=i, name=glyphName)
            writer.newline()

    def fromXML(
        self, name: str, attrs: dict[str, str], content: list[Any], ttFont: TTFont
    ) -> None:
        if not hasattr(self, "glyphOrder"):
            self.glyphOrder = []
        if name == "GlyphID":
            self.glyphOrder.append(attrs["name"])
        ttFont.setGlyphOrder(self.glyphOrder)


def getTableModule(tag: str | bytes) -> ModuleType | None:
    """Fetch the packer/unpacker module for a table.
    Return None when no module is found.
    """
    from . import tables

    pyTag = tagToIdentifier(tag)
    try:
        __import__("fontTools.ttLib.tables." + pyTag)
    except ImportError as err:
        # If pyTag is found in the ImportError message,
        # means table is not implemented.  If it's not
        # there, then some other module is missing, don't
        # suppress the error.
        if str(err).find(pyTag) >= 0:
            return None
        else:
            raise err
    else:
        return getattr(tables, pyTag)


# Registry for custom table packer/unpacker classes. Keys are table
# tags, values are (moduleName, className) tuples.
# See registerCustomTableClass() and getCustomTableClass()
_customTableRegistry: dict[str | bytes, tuple[str, str]] = {}


def registerCustomTableClass(
    tag: str | bytes, moduleName: str, className: str | None = None
) -> None:
    """Register a custom packer/unpacker class for a table.

    The 'moduleName' must be an importable module. If no 'className'
    is given, it is derived from the tag, for example it will be
    ``table_C_U_S_T_`` for a 'CUST' tag.

    The registered table class should be a subclass of
    :py:class:`fontTools.ttLib.tables.DefaultTable.DefaultTable`
    """
    if className is None:
        className = "table_" + tagToIdentifier(tag)
    _customTableRegistry[tag] = (moduleName, className)


def unregisterCustomTableClass(tag: str | bytes) -> None:
    """Unregister the custom packer/unpacker class for a table."""
    del _customTableRegistry[tag]


def getCustomTableClass(tag: str | bytes) -> type[DefaultTable] | None:
    """Return the custom table class for tag, if one has been registered
    with 'registerCustomTableClass()'. Else return None.
    """
    if tag not in _customTableRegistry:
        return None
    import importlib

    moduleName, className = _customTableRegistry[tag]
    module = importlib.import_module(moduleName)
    return getattr(module, className)


def getTableClass(tag: str | bytes) -> type[DefaultTable]:
    """Fetch the packer/unpacker class for a table."""
    tableClass = getCustomTableClass(tag)
    if tableClass is not None:
        return tableClass
    module = getTableModule(tag)
    if module is None:
        from .tables.DefaultTable import DefaultTable

        return DefaultTable
    pyTag = tagToIdentifier(tag)
    tableClass = getattr(module, "table_" + pyTag)
    return tableClass


def getClassTag(klass: type[DefaultTable]) -> str | bytes:
    """Fetch the table tag for a class object."""
    name = klass.__name__
    assert name[:6] == "table_"
    name = name[6:]  # Chop 'table_'
    return identifierToTag(name)


def newTable(tag: str | bytes) -> DefaultTable:
    """Return a new instance of a table."""
    tableClass = getTableClass(tag)
    return tableClass(tag)


def _escapechar(c: str) -> str:
    """Helper function for tagToIdentifier()"""
    import re

    if re.match("[a-z0-9]", c):
        return "_" + c
    elif re.match("[A-Z]", c):
        return c + "_"
    else:
        return hex(byteord(c))[2:]


def tagToIdentifier(tag: str | bytes) -> str:
    """Convert a table tag to a valid (but UGLY) python identifier,
    as well as a filename that's guaranteed to be unique even on a
    caseless file system. Each character is mapped to two characters.
    Lowercase letters get an underscore before the letter, uppercase
    letters get an underscore after the letter. Trailing spaces are
    trimmed. Illegal characters are escaped as two hex bytes. If the
    result starts with a number (as the result of a hex escape), an
    extra underscore is prepended. Examples:
    .. code-block:: pycon

        >>>
        >> tagToIdentifier('glyf')
        '_g_l_y_f'
        >> tagToIdentifier('cvt ')
        '_c_v_t'
        >> tagToIdentifier('OS/2')
        'O_S_2f_2'
    """
    import re

    tag = Tag(tag)
    if tag == "GlyphOrder":
        return tag
    assert len(tag) == 4, "tag should be 4 characters long"
    while len(tag) > 1 and tag[-1] == " ":
        tag = tag[:-1]
    ident = ""
    for c in tag:
        ident = ident + _escapechar(c)
    if re.match("[0-9]", ident):
        ident = "_" + ident
    return ident


def identifierToTag(ident: str) -> str:
    """the opposite of tagToIdentifier()"""
    if ident == "GlyphOrder":
        return ident
    if len(ident) % 2 and ident[0] == "_":
        ident = ident[1:]
    assert not (len(ident) % 2)
    tag = ""
    for i in range(0, len(ident), 2):
        if ident[i] == "_":
            tag = tag + ident[i + 1]
        elif ident[i + 1] == "_":
            tag = tag + ident[i]
        else:
            # assume hex
            tag = tag + chr(int(ident[i : i + 2], 16))
    # append trailing spaces
    tag = tag + (4 - len(tag)) * " "
    return Tag(tag)


def tagToXML(tag: str | bytes) -> str:
    """Similarly to tagToIdentifier(), this converts a TT tag
    to a valid XML element name. Since XML element names are
    case sensitive, this is a fairly simple/readable translation.
    """
    import re

    tag = Tag(tag)
    if tag == "OS/2":
        return "OS_2"
    elif tag == "GlyphOrder":
        return tag
    if re.match("[A-Za-z_][A-Za-z_0-9]* *$", tag):
        return tag.strip()
    else:
        return tagToIdentifier(tag)


def xmlToTag(tag: str) -> str:
    """The opposite of tagToXML()"""
    if tag == "OS_2":
        return Tag("OS/2")
    if len(tag) == 8:
        return identifierToTag(tag)
    else:
        return Tag(tag + " " * (4 - len(tag)))


# Table order as recommended in the OpenType specification 1.4
TTFTableOrder = [
    "head",
    "hhea",
    "maxp",
    "OS/2",
    "hmtx",
    "LTSH",
    "VDMX",
    "hdmx",
    "cmap",
    "fpgm",
    "prep",
    "cvt ",
    "loca",
    "glyf",
    "kern",
    "name",
    "post",
    "gasp",
    "PCLT",
]

OTFTableOrder = ["head", "hhea", "maxp", "OS/2", "name", "cmap", "post", "CFF "]


def sortedTagList(
    tagList: Sequence[str], tableOrder: Sequence[str] | None = None
) -> list[str]:
    """Return a sorted copy of tagList, sorted according to the OpenType
    specification, or according to a custom tableOrder. If given and not
    None, tableOrder needs to be a list of tag names.
    """
    tagList = sorted(tagList)
    if tableOrder is None:
        if "DSIG" in tagList:
            # DSIG should be last (XXX spec reference?)
            tagList.remove("DSIG")
            tagList.append("DSIG")
        if "CFF " in tagList:
            tableOrder = OTFTableOrder
        else:
            tableOrder = TTFTableOrder
    orderedTables = []
    for tag in tableOrder:
        if tag in tagList:
            orderedTables.append(tag)
            tagList.remove(tag)
    orderedTables.extend(tagList)
    return orderedTables


def reorderFontTables(
    inFile: BinaryIO,  # Takes file-like object as per original
    outFile: BinaryIO,  # Takes file-like object
    tableOrder: Sequence[str] | None = None,
    checkChecksums: bool = False,  # Keep param even if reader handles it
) -> None:
    """Rewrite a font file, ordering the tables as recommended by the
    OpenType specification 1.4.
    """
    inFile.seek(0)
    outFile.seek(0)
    reader = SFNTReader(inFile, checkChecksums=checkChecksums)
    writer = SFNTWriter(
        outFile,
        len(reader.tables),
        reader.sfntVersion,
        reader.flavor,
        reader.flavorData,
    )
    tables = list(reader.keys())
    for tag in sortedTagList(tables, tableOrder):
        writer[tag] = reader[tag]
    writer.close()


def maxPowerOfTwo(x: int) -> int:
    """Return the highest exponent of two, so that
    (2 ** exponent) <= x.  Return 0 if x is 0.
    """
    exponent = 0
    while x:
        x = x >> 1
        exponent = exponent + 1
    return max(exponent - 1, 0)


def getSearchRange(n: int, itemSize: int = 16) -> tuple[int, int, int]:
    """Calculate searchRange, entrySelector, rangeShift."""
    # itemSize defaults to 16, for backward compatibility
    # with upstream fonttools.
    exponent = maxPowerOfTwo(n)
    searchRange = (2**exponent) * itemSize
    entrySelector = exponent
    rangeShift = max(0, n * itemSize - searchRange)
    return searchRange, entrySelector, rangeShift
