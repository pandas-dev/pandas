"""Utilities for conversion to writer-agnostic Excel representation
"""

from functools import reduce
import itertools
import re
from typing import Callable, Dict, List, Optional, Sequence, Union
import warnings

import numpy as np

from pandas.core.dtypes import missing
from pandas.core.dtypes.common import is_float, is_scalar
from pandas.core.dtypes.generic import ABCMultiIndex, ABCPeriodIndex

from pandas import Index
import pandas.core.common as com

from pandas.io.common import stringify_path
from pandas.io.formats.css import CSSResolver, CSSWarning
from pandas.io.formats.format import get_level_lengths
from pandas.io.formats.printing import pprint_thing


class ExcelCell:
    __fields__ = ("row", "col", "val", "style", "mergestart", "mergeend")
    __slots__ = __fields__

    def __init__(
        self, row: int, col: int, val, style=None, mergestart=None, mergeend=None
    ):
        self.row = row
        self.col = col
        self.val = val
        self.style = style
        self.mergestart = mergestart
        self.mergeend = mergeend


class CSSToExcelConverter:
    """A callable for converting CSS declarations to ExcelWriter styles

    Supports parts of CSS 2.2, with minimal CSS 3.0 support (e.g. text-shadow),
    focusing on font styling, backgrounds, borders and alignment.

    Operates by first computing CSS styles in a fairly generic
    way (see :meth:`compute_css`) then determining Excel style
    properties from CSS properties (see :meth:`build_xlstyle`).

    Parameters
    ----------
    inherited : str, optional
        CSS declarations understood to be the containing scope for the
        CSS processed by :meth:`__call__`.
    """

    # NB: Most of the methods here could be classmethods, as only __init__
    #     and __call__ make use of instance attributes.  We leave them as
    #     instancemethods so that users can easily experiment with extensions
    #     without monkey-patching.

    def __init__(self, inherited: Optional[str] = None):
        if inherited is not None:
            inherited = self.compute_css(inherited)

        self.inherited = inherited

    compute_css = CSSResolver()

    def __call__(self, declarations_str: str) -> Dict[str, Dict[str, str]]:
        """
        Convert CSS declarations to ExcelWriter style.

        Parameters
        ----------
        declarations_str : str
            List of CSS declarations.
            e.g. "font-weight: bold; background: blue"

        Returns
        -------
        xlstyle : dict
            A style as interpreted by ExcelWriter when found in
            ExcelCell.style.
        """
        # TODO: memoize?
        properties = self.compute_css(declarations_str, self.inherited)
        return self.build_xlstyle(properties)

    def build_xlstyle(self, props: Dict[str, str]) -> Dict[str, Dict[str, str]]:
        out = {
            "alignment": self.build_alignment(props),
            "border": self.build_border(props),
            "fill": self.build_fill(props),
            "font": self.build_font(props),
            "number_format": self.build_number_format(props),
        }

        # TODO: handle cell width and height: needs support in pandas.io.excel

        def remove_none(d: Dict[str, str]) -> None:
            """Remove key where value is None, through nested dicts"""
            for k, v in list(d.items()):
                if v is None:
                    del d[k]
                elif isinstance(v, dict):
                    remove_none(v)
                    if not v:
                        del d[k]

        remove_none(out)
        return out

    VERTICAL_MAP = {
        "top": "top",
        "text-top": "top",
        "middle": "center",
        "baseline": "bottom",
        "bottom": "bottom",
        "text-bottom": "bottom",
        # OpenXML also has 'justify', 'distributed'
    }

    def build_alignment(self, props) -> Dict[str, Optional[Union[bool, str]]]:
        # TODO: text-indent, padding-left -> alignment.indent
        return {
            "horizontal": props.get("text-align"),
            "vertical": self.VERTICAL_MAP.get(props.get("vertical-align")),
            "wrap_text": (
                None
                if props.get("white-space") is None
                else props["white-space"] not in ("nowrap", "pre", "pre-line")
            ),
        }

    def build_border(self, props: Dict) -> Dict[str, Dict[str, str]]:
        return {
            side: {
                "style": self._border_style(
                    props.get(f"border-{side}-style"),
                    props.get(f"border-{side}-width"),
                ),
                "color": self.color_to_excel(props.get(f"border-{side}-color")),
            }
            for side in ["top", "right", "bottom", "left"]
        }

    def _border_style(self, style: Optional[str], width):
        # convert styles and widths to openxml, one of:
        #       'dashDot'
        #       'dashDotDot'
        #       'dashed'
        #       'dotted'
        #       'double'
        #       'hair'
        #       'medium'
        #       'mediumDashDot'
        #       'mediumDashDotDot'
        #       'mediumDashed'
        #       'slantDashDot'
        #       'thick'
        #       'thin'
        if width is None and style is None:
            return None
        if style == "none" or style == "hidden":
            return None

        if width is None:
            width = "2pt"
        width = float(width[:-2])
        if width < 1e-5:
            return None
        elif width < 1.3:
            width_name = "thin"
        elif width < 2.8:
            width_name = "medium"
        else:
            width_name = "thick"

        if style in (None, "groove", "ridge", "inset", "outset"):
            # not handled
            style = "solid"

        if style == "double":
            return "double"
        if style == "solid":
            return width_name
        if style == "dotted":
            if width_name in ("hair", "thin"):
                return "dotted"
            return "mediumDashDotDot"
        if style == "dashed":
            if width_name in ("hair", "thin"):
                return "dashed"
            return "mediumDashed"

    def build_fill(self, props: Dict[str, str]):
        # TODO: perhaps allow for special properties
        #       -excel-pattern-bgcolor and -excel-pattern-type
        fill_color = props.get("background-color")
        if fill_color not in (None, "transparent", "none"):
            return {"fgColor": self.color_to_excel(fill_color), "patternType": "solid"}

    BOLD_MAP = {
        "bold": True,
        "bolder": True,
        "600": True,
        "700": True,
        "800": True,
        "900": True,
        "normal": False,
        "lighter": False,
        "100": False,
        "200": False,
        "300": False,
        "400": False,
        "500": False,
    }
    ITALIC_MAP = {"normal": False, "italic": True, "oblique": True}

    def build_font(self, props) -> Dict[str, Optional[Union[bool, int, str]]]:
        size = props.get("font-size")
        if size is not None:
            assert size.endswith("pt")
            size = float(size[:-2])

        font_names_tmp = re.findall(
            r"""(?x)
            (
            "(?:[^"]|\\")+"
            |
            '(?:[^']|\\')+'
            |
            [^'",]+
            )(?=,|\s*$)
        """,
            props.get("font-family", ""),
        )
        font_names = []
        for name in font_names_tmp:
            if name[:1] == '"':
                name = name[1:-1].replace('\\"', '"')
            elif name[:1] == "'":
                name = name[1:-1].replace("\\'", "'")
            else:
                name = name.strip()
            if name:
                font_names.append(name)

        family = None
        for name in font_names:
            if name == "serif":
                family = 1  # roman
                break
            elif name == "sans-serif":
                family = 2  # swiss
                break
            elif name == "cursive":
                family = 4  # script
                break
            elif name == "fantasy":
                family = 5  # decorative
                break

        decoration = props.get("text-decoration")
        if decoration is not None:
            decoration = decoration.split()
        else:
            decoration = ()

        return {
            "name": font_names[0] if font_names else None,
            "family": family,
            "size": size,
            "bold": self.BOLD_MAP.get(props.get("font-weight")),
            "italic": self.ITALIC_MAP.get(props.get("font-style")),
            "underline": ("single" if "underline" in decoration else None),
            "strike": ("line-through" in decoration) or None,
            "color": self.color_to_excel(props.get("color")),
            # shadow if nonzero digit before shadow color
            "shadow": (
                bool(re.search("^[^#(]*[1-9]", props["text-shadow"]))
                if "text-shadow" in props
                else None
            ),
            # FIXME: dont leave commented-out
            # 'vertAlign':,
            # 'charset': ,
            # 'scheme': ,
            # 'outline': ,
            # 'condense': ,
        }

    NAMED_COLORS = {
        "maroon": "800000",
        "brown": "A52A2A",
        "red": "FF0000",
        "pink": "FFC0CB",
        "orange": "FFA500",
        "yellow": "FFFF00",
        "olive": "808000",
        "green": "008000",
        "purple": "800080",
        "fuchsia": "FF00FF",
        "lime": "00FF00",
        "teal": "008080",
        "aqua": "00FFFF",
        "blue": "0000FF",
        "navy": "000080",
        "black": "000000",
        "gray": "808080",
        "grey": "808080",
        "silver": "C0C0C0",
        "white": "FFFFFF",
    }

    def color_to_excel(self, val: Optional[str]):
        if val is None:
            return None
        if val.startswith("#") and len(val) == 7:
            return val[1:].upper()
        if val.startswith("#") and len(val) == 4:
            return (val[1] * 2 + val[2] * 2 + val[3] * 2).upper()
        try:
            return self.NAMED_COLORS[val]
        except KeyError:
            warnings.warn(f"Unhandled color format: {repr(val)}", CSSWarning)

    def build_number_format(self, props: Dict) -> Dict[str, Optional[str]]:
        return {"format_code": props.get("number-format")}


class ExcelFormatter:
    """
    Class for formatting a DataFrame to a list of ExcelCells,

    Parameters
    ----------
    df : DataFrame or Styler
    na_rep: na representation
    float_format : string, default None
            Format string for floating point numbers
    cols : sequence, optional
        Columns to write
    header : boolean or list of string, default True
        Write out column names. If a list of string is given it is
        assumed to be aliases for the column names
    index : boolean, default True
        output row names (index)
    index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
    merge_cells : boolean, default False
            Format MultiIndex and Hierarchical Rows as merged cells.
    inf_rep : string, default `'inf'`
        representation for np.inf values (which aren't representable in Excel)
        A `'-'` sign will be added in front of -inf.
    style_converter : callable, optional
        This translates Styler styles (CSS) into ExcelWriter styles.
        Defaults to ``CSSToExcelConverter()``.
        It should have signature css_declarations string -> excel style.
        This is only called for body cells.
    """

    max_rows = 2 ** 20
    max_cols = 2 ** 14

    def __init__(
        self,
        df,
        na_rep: str = "",
        float_format: Optional[str] = None,
        cols: Optional[Sequence] = None,
        header: Union[bool, List[str]] = True,
        index: bool = True,
        index_label: Union[str, Sequence, None] = None,
        merge_cells: bool = False,
        inf_rep: str = "inf",
        style_converter: Optional[Callable] = None,
    ):
        self.rowcounter = 0
        self.na_rep = na_rep
        if hasattr(df, "render"):
            self.styler = df
            df = df.data
            if style_converter is None:
                style_converter = CSSToExcelConverter()
            self.style_converter = style_converter
        else:
            self.styler = None
        self.df = df
        if cols is not None:

            # all missing, raise
            if not len(Index(cols) & df.columns):
                raise KeyError("passes columns are not ALL present dataframe")

            if len(Index(cols) & df.columns) != len(cols):
                # Deprecated in GH#17295, enforced in 1.0.0
                raise KeyError("Not all names specified in 'columns' are found")

            self.df = df.reindex(columns=cols)

        self.columns = self.df.columns
        self.float_format = float_format
        self.index = index
        self.index_label = index_label
        self.header = header
        self.merge_cells = merge_cells
        self.inf_rep = inf_rep

    @property
    def header_style(self):
        return {
            "font": {"bold": True},
            "borders": {
                "top": "thin",
                "right": "thin",
                "bottom": "thin",
                "left": "thin",
            },
            "alignment": {"horizontal": "center", "vertical": "top"},
        }

    def _format_value(self, val):
        if is_scalar(val) and missing.isna(val):
            val = self.na_rep
        elif is_float(val):
            if missing.isposinf_scalar(val):
                val = self.inf_rep
            elif missing.isneginf_scalar(val):
                val = f"-{self.inf_rep}"
            elif self.float_format is not None:
                val = float(self.float_format % val)
        if getattr(val, "tzinfo", None) is not None:
            raise ValueError(
                "Excel does not support datetimes with "
                "timezones. Please ensure that datetimes "
                "are timezone unaware before writing to Excel."
            )
        return val

    def _format_header_mi(self):
        if self.columns.nlevels > 1:
            if not self.index:
                raise NotImplementedError(
                    "Writing to Excel with MultiIndex columns and no "
                    "index ('index'=False) is not yet implemented."
                )

        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if not (has_aliases or self.header):
            return

        columns = self.columns
        level_strs = columns.format(
            sparsify=self.merge_cells, adjoin=False, names=False
        )
        level_lengths = get_level_lengths(level_strs)
        coloffset = 0
        lnum = 0

        if self.index and isinstance(self.df.index, ABCMultiIndex):
            coloffset = len(self.df.index[0]) - 1

        if self.merge_cells:
            # Format multi-index as a merged cells.
            for lnum in range(len(level_lengths)):
                name = columns.names[lnum]
                yield ExcelCell(lnum, coloffset, name, self.header_style)

            for lnum, (spans, levels, level_codes) in enumerate(
                zip(level_lengths, columns.levels, columns.codes)
            ):
                values = levels.take(level_codes)
                for i in spans:
                    if spans[i] > 1:
                        yield ExcelCell(
                            lnum,
                            coloffset + i + 1,
                            values[i],
                            self.header_style,
                            lnum,
                            coloffset + i + spans[i],
                        )
                    else:
                        yield ExcelCell(
                            lnum, coloffset + i + 1, values[i], self.header_style
                        )
        else:
            # Format in legacy format with dots to indicate levels.
            for i, values in enumerate(zip(*level_strs)):
                v = ".".join(map(pprint_thing, values))
                yield ExcelCell(lnum, coloffset + i + 1, v, self.header_style)

        self.rowcounter = lnum

    def _format_header_regular(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            coloffset = 0

            if self.index:
                coloffset = 1
                if isinstance(self.df.index, ABCMultiIndex):
                    coloffset = len(self.df.index[0])

            colnames = self.columns
            if has_aliases:
                if len(self.header) != len(self.columns):
                    raise ValueError(
                        f"Writing {len(self.columns)} cols but got {len(self.header)} "
                        "aliases"
                    )
                else:
                    colnames = self.header

            for colindex, colname in enumerate(colnames):
                yield ExcelCell(
                    self.rowcounter, colindex + coloffset, colname, self.header_style
                )

    def _format_header(self):
        if isinstance(self.columns, ABCMultiIndex):
            gen = self._format_header_mi()
        else:
            gen = self._format_header_regular()

        gen2 = ()
        if self.df.index.names:
            row = [x if x is not None else "" for x in self.df.index.names] + [
                ""
            ] * len(self.columns)
            if reduce(lambda x, y: x and y, map(lambda x: x != "", row)):
                gen2 = (
                    ExcelCell(self.rowcounter, colindex, val, self.header_style)
                    for colindex, val in enumerate(row)
                )
                self.rowcounter += 1
        return itertools.chain(gen, gen2)

    def _format_body(self):
        if isinstance(self.df.index, ABCMultiIndex):
            return self._format_hierarchical_rows()
        else:
            return self._format_regular_rows()

    def _format_regular_rows(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            self.rowcounter += 1

        # output index and index_label?
        if self.index:
            # check aliases
            # if list only take first as this is not a MultiIndex
            if self.index_label and isinstance(
                self.index_label, (list, tuple, np.ndarray, Index)
            ):
                index_label = self.index_label[0]
            # if string good to go
            elif self.index_label and isinstance(self.index_label, str):
                index_label = self.index_label
            else:
                index_label = self.df.index.names[0]

            if isinstance(self.columns, ABCMultiIndex):
                self.rowcounter += 1

            if index_label and self.header is not False:
                yield ExcelCell(self.rowcounter - 1, 0, index_label, self.header_style)

            # write index_values
            index_values = self.df.index
            if isinstance(self.df.index, ABCPeriodIndex):
                index_values = self.df.index.to_timestamp()

            for idx, idxval in enumerate(index_values):
                yield ExcelCell(self.rowcounter + idx, 0, idxval, self.header_style)

            coloffset = 1
        else:
            coloffset = 0

        for cell in self._generate_body(coloffset):
            yield cell

    def _format_hierarchical_rows(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            self.rowcounter += 1

        gcolidx = 0

        if self.index:
            index_labels = self.df.index.names
            # check for aliases
            if self.index_label and isinstance(
                self.index_label, (list, tuple, np.ndarray, Index)
            ):
                index_labels = self.index_label

            # MultiIndex columns require an extra row
            # with index names (blank if None) for
            # unambiguous round-trip, unless not merging,
            # in which case the names all go on one row Issue #11328
            if isinstance(self.columns, ABCMultiIndex) and self.merge_cells:
                self.rowcounter += 1

            # if index labels are not empty go ahead and dump
            if com.any_not_none(*index_labels) and self.header is not False:

                for cidx, name in enumerate(index_labels):
                    yield ExcelCell(self.rowcounter - 1, cidx, name, self.header_style)

            if self.merge_cells:
                # Format hierarchical rows as merged cells.
                level_strs = self.df.index.format(
                    sparsify=True, adjoin=False, names=False
                )
                level_lengths = get_level_lengths(level_strs)

                for spans, levels, level_codes in zip(
                    level_lengths, self.df.index.levels, self.df.index.codes
                ):

                    values = levels.take(
                        level_codes, allow_fill=levels._can_hold_na, fill_value=True
                    )

                    for i in spans:
                        if spans[i] > 1:
                            yield ExcelCell(
                                self.rowcounter + i,
                                gcolidx,
                                values[i],
                                self.header_style,
                                self.rowcounter + i + spans[i] - 1,
                                gcolidx,
                            )
                        else:
                            yield ExcelCell(
                                self.rowcounter + i,
                                gcolidx,
                                values[i],
                                self.header_style,
                            )
                    gcolidx += 1

            else:
                # Format hierarchical rows with non-merged values.
                for indexcolvals in zip(*self.df.index):
                    for idx, indexcolval in enumerate(indexcolvals):
                        yield ExcelCell(
                            self.rowcounter + idx,
                            gcolidx,
                            indexcolval,
                            self.header_style,
                        )
                    gcolidx += 1

        for cell in self._generate_body(gcolidx):
            yield cell

    def _generate_body(self, coloffset: int):
        if self.styler is None:
            styles = None
        else:
            styles = self.styler._compute().ctx
            if not styles:
                styles = None
        xlstyle = None

        # Write the body of the frame data series by series.
        for colidx in range(len(self.columns)):
            series = self.df.iloc[:, colidx]
            for i, val in enumerate(series):
                if styles is not None:
                    xlstyle = self.style_converter(";".join(styles[i, colidx]))
                yield ExcelCell(self.rowcounter + i, colidx + coloffset, val, xlstyle)

    def get_formatted_cells(self):
        for cell in itertools.chain(self._format_header(), self._format_body()):
            cell.val = self._format_value(cell.val)
            yield cell

    def write(
        self,
        writer,
        sheet_name="Sheet1",
        startrow=0,
        startcol=0,
        freeze_panes=None,
        engine=None,
    ):
        """
        writer : string or ExcelWriter object
            File path or existing ExcelWriter
        sheet_name : string, default 'Sheet1'
            Name of sheet which will contain DataFrame
        startrow :
            upper left cell row to dump data frame
        startcol :
            upper left cell column to dump data frame
        freeze_panes : tuple of integer (length 2), default None
            Specifies the one-based bottommost row and rightmost column that
            is to be frozen
        engine : string, default None
            write engine to use if writer is a path - you can also set this
            via the options ``io.excel.xlsx.writer``, ``io.excel.xls.writer``,
            and ``io.excel.xlsm.writer``.
        """
        from pandas.io.excel import ExcelWriter

        num_rows, num_cols = self.df.shape
        if num_rows > self.max_rows or num_cols > self.max_cols:
            raise ValueError(
                f"This sheet is too large! Your sheet size is: {num_rows}, {num_cols} "
                f"Max sheet size is: {self.max_rows}, {self.max_cols}"
            )

        if isinstance(writer, ExcelWriter):
            need_save = False
        else:
            writer = ExcelWriter(stringify_path(writer), engine=engine)
            need_save = True

        formatted_cells = self.get_formatted_cells()
        writer.write_cells(
            formatted_cells,
            sheet_name,
            startrow=startrow,
            startcol=startcol,
            freeze_panes=freeze_panes,
        )
        if need_save:
            writer.save()
