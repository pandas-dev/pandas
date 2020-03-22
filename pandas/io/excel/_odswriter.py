from collections import defaultdict
import datetime

import pandas._libs.json as json

from pandas.io.excel._base import ExcelWriter

# from pandas.io.excel._util import _validate_freeze_panes

from odf.opendocument import OpenDocumentSpreadsheet
from odf.style import Style, TextProperties, TableCellProperties, ParagraphProperties
from odf.table import Table, TableRow, TableCell
from odf.text import P


class _ODSWriter(ExcelWriter):
    engine = "odf"
    supported_extensions = (".ods",)

    def __init__(self, path, engine=None, encoding=None, mode="w", **engine_kwargs):
        engine_kwargs["engine"] = engine

        if mode == "a":
            raise ValueError("Append mode is not supported with odf!")

        super().__init__(path, mode=mode, **engine_kwargs)

        self.book = OpenDocumentSpreadsheet()
        self.style_dict = {}

    def save(self):
        """
        Save workbook to disk.
        """
        for sheet in self.sheets.values():
            self.book.spreadsheet.addElement(sheet)
        return self.book.save(self.path)

    def write_cells(
        self, cells, sheet_name=None, startrow=0, startcol=0, freeze_panes=None
    ):
        # Write the frame cells using odf

        sheet_name = self._get_sheet_name(sheet_name)

        if sheet_name in self.sheets:
            wks = self.sheets[sheet_name]
        else:
            wks = Table(name=sheet_name)
            self.sheets[sheet_name] = wks

        #        if _validate_freeze_panes(freeze_panes):
        #            wks.set_panes_frozen(True)
        #            wks.set_horz_split_pos(freeze_panes[0])
        #            wks.set_vert_split_pos(freeze_panes[1])


        rows = defaultdict(TableRow)
        col_count = defaultdict(int)

        for cell in sorted(cells, key=lambda cell: (cell.row, cell.col)):
            attributes = {}
            style_name = self._process_style(cell.style)
            if style_name is not None:
                attributes["stylename"] = style_name
            print(cell.row, cell.col, cell.val, cell.mergestart, cell.mergeend)
            if cell.mergestart is not None and cell.mergeend is not None:
                attributes["numberrowsspanned"] = max(1, cell.mergestart)
                attributes["numbercolumnsspanned"] = cell.mergeend
            # fill with empty cells if needed
            for _ in range(cell.col - col_count[cell.row]):
                rows[cell.row].addElement(TableCell())
                col_count[cell.row] += 1
            val, fmt = self._value_with_fmt(cell.val)
            # print("type", type(val), "value", val)
            pvalue = value = val
            if isinstance(val, bool):
                value = str(val).lower()
                pvalue = str(val).upper()
            if isinstance(val, datetime.datetime):
                if val.time():
                    value = val.isoformat()
                    pvalue = val.strftime("%c")
                else:
                    value = val.strftime("%Y-%m-%d")
                    pvalue = val.strftime("%x")
                tc = TableCell(valuetype="date", datevalue=value, attributes=attributes)
            elif isinstance(val, datetime.date):
                value = val.strftime("%Y-%m-%d")
                pvalue = val.strftime("%x")
                tc = TableCell(valuetype="date", datevalue=value, attributes=attributes)
            else:
                class_to_cell_type = {
                    str: "string",
                    int: "float",
                    float: "float",
                    bool: "boolean",
                }
                tc = TableCell(
                    valuetype=class_to_cell_type[type(val)],
                    value=value,
                    attributes=attributes,
                )
            rows[cell.row].addElement(tc)
            col_count[cell.row] += 1
            p = P(text=pvalue)
            tc.addElement(p)
        for row_nr in range(max(rows.keys()) + 1):
            wks.addElement(rows[row_nr])

    def _process_style(self, style):
        if style is None:
            return None
        style_key = json.dumps(style)
        if style_key in self.style_dict:
            return self.style_dict[style_key]
        name = f"pd{len(self.style_dict)+1}"
        self.style_dict[style_key] = name
        odf_style = Style(name=name, family="table-cell")
        if "font" in style:
            font = style["font"]
            if font.get("bold", False):
                odf_style.addElement(TextProperties(fontweight="bold"))
        if "borders" in style:
            borders = style["borders"]
            for side, thickness in borders.items():
                thickness_translation = {
                    "thin": "0.75pt solid #000000"
                }
                odf_style.addElement(
                    TableCellProperties(
                        attributes={f"border{side}": thickness_translation[thickness]}))
        if "alignment" in style:
            alignment = style["alignment"]
            horizontal = alignment.get("horizontal")
            if horizontal:
                odf_style.addElement(ParagraphProperties(textalign=horizontal))
            vertical = alignment.get("vertical")
            if vertical:
                odf_style.addElement(TableCellProperties(verticalalign=vertical))
        self.book.styles.addElement(odf_style)
        return name
