from __future__ import annotations

import json
from typing import (
    TYPE_CHECKING,
    Any,
)

if TYPE_CHECKING:
    from pandas._typing import FilePath, StorageOptions, WriteExcelBuffer

from xlsxwriter import Workbook

from pandas.compat._optional import import_optional_dependency

from pandas.io.excel._base import ExcelWriter
from pandas.io.excel._util import validate_freeze_panes


class _XlsxStyler:
    # Map from openpyxl-oriented styles to flatter xlsxwriter representation
    # Ordering necessary for both determinism and because some are keyed by
    # prefixes of others.
    STYLE_MAPPING: dict[str, list[tuple[tuple[str, ...], str]]] = {
        "font": [
            (("name",), "font_name"),
            (("sz",), "font_size"),
            (("size",), "font_size"),
            (("color", "rgb"), "font_color"),
            (("color",), "font_color"),
            (("b",), "bold"),
            (("bold",), "bold"),
            (("i",), "italic"),
            (("italic",), "italic"),
            (("u",), "underline"),
            (("underline",), "underline"),
            (("strike",), "font_strikeout"),
            (("vertAlign",), "font_script"),
            (("vertalign",), "font_script"),
        ],
        "number_format": [(("format_code",), "num_format"), ((), "num_format")],
        "protection": [(("locked",), "locked"), (("hidden",), "hidden")],
        "alignment": [
            (("horizontal",), "align"),
            (("vertical",), "valign"),
            (("text_rotation",), "rotation"),
            (("wrap_text",), "text_wrap"),
            (("indent",), "indent"),
            (("shrink_to_fit",), "shrink"),
        ],
        "fill": [
            (("patternType",), "pattern"),
            (("patterntype",), "pattern"),
            (("fill_type",), "pattern"),
            (("start_color", "rgb"), "fg_color"),
            (("fgColor", "rgb"), "fg_color"),
            (("fgcolor", "rgb"), "fg_color"),
            (("start_color",), "fg_color"),
            (("fgColor",), "fg_color"),
            (("fgcolor",), "fg_color"),
            (("end_color", "rgb"), "bg_color"),
            (("bgColor", "rgb"), "bg_color"),
            (("bgcolor", "rgb"), "bg_color"),
            (("end_color",), "bg_color"),
            (("bgColor",), "bg_color"),
            (("bgcolor",), "bg_color"),
        ],
        "border": [
            (("color", "rgb"), "border_color"),
            (("color",), "border_color"),
            (("style",), "border"),
            (("top", "color", "rgb"), "top_color"),
            (("top", "color"), "top_color"),
            (("top", "style"), "top"),
            (("top",), "top"),
            (("right", "color", "rgb"), "right_color"),
            (("right", "color"), "right_color"),
            (("right", "style"), "right"),
            (("right",), "right"),
            (("bottom", "color", "rgb"), "bottom_color"),
            (("bottom", "color"), "bottom_color"),
            (("bottom", "style"), "bottom"),
            (("bottom",), "bottom"),
            (("left", "color", "rgb"), "left_color"),
            (("left", "color"), "left_color"),
            (("left", "style"), "left"),
            (("left",), "left"),
        ],
    }

    @classmethod
    def convert(
        cls,
        style_dict: dict | None,
        num_format_str: str | None = None,
    ) -> dict[str, Any]:
        """Convert a style_dict to an xlsxwriter format dict."""
        # Normalize and copy to avoid modifying the input
        if style_dict is None:
            style_dict = {}
        else:
            style_dict = style_dict.copy()

        # Map CSS font-weight to xlsxwriter font-weight (bold)
        if style_dict.get("font-weight") in ("bold", "bolder", 700, "700") or (
            isinstance(style_dict.get("font"), dict)
            and style_dict["font"].get("weight") in ("bold", "bolder", 700, "700")
        ):
            # For XLSXWriter, we need to set the font with bold=True
            style_dict = {"font": {"bold": True, "name": "Calibri", "size": 11}}
            # Also set the b property directly as it might be needed
            style_dict["b"] = True

        # Handle font styles
        if "font-style" in style_dict and style_dict["font-style"] == "italic":
            style_dict["italic"] = True
            del style_dict["font-style"]

        # Convert CSS border styles to xlsxwriter format
        #        border_map = {
        #            "border-top": "top",
        #           "border-right": "right",
        #         "border-bottom": "bottom",
        #        "border-left": "left",
        #    }
        if "borders" in style_dict:
            style_dict = style_dict.copy()
            style_dict["border"] = style_dict.pop("borders")

        # Initialize props to track which properties we've processed
        props = {}

        for style_group_key, style_group in style_dict.items():
            for src, dst in cls.STYLE_MAPPING.get(style_group_key, []):
                # src is a sequence of keys into a nested dict
                # dst is a flat key
                if dst in props:
                    continue
                v = style_group
                for k in src:
                    try:
                        v = v[k]
                    except (KeyError, TypeError):
                        break
                else:
                    props[dst] = v

        if isinstance(props.get("pattern"), str):
            # TODO: support other fill patterns
            props["pattern"] = 0 if props["pattern"] == "none" else 1

        for k in ["border", "top", "right", "bottom", "left"]:
            if isinstance(props.get(k), str):
                try:
                    props[k] = [
                        "none",
                        "thin",
                        "medium",
                        "dashed",
                        "dotted",
                        "thick",
                        "double",
                        "hair",
                        "mediumDashed",
                        "dashDot",
                        "mediumDashDot",
                        "dashDotDot",
                        "mediumDashDotDot",
                        "slantDashDot",
                    ].index(props[k])
                except ValueError:
                    props[k] = 2

        if isinstance(props.get("font_script"), str):
            props["font_script"] = ["baseline", "superscript", "subscript"].index(
                props["font_script"]
            )

        if isinstance(props.get("underline"), str):
            props["underline"] = {
                "none": 0,
                "single": 1,
                "double": 2,
                "singleAccounting": 33,
                "doubleAccounting": 34,
            }[props["underline"]]

        # GH 30107 - xlsxwriter uses different name
        if props.get("valign") == "center":
            props["valign"] = "vcenter"

        # Ensure numeric format is applied when provided separately
        if num_format_str and "num_format" not in props:
            props["num_format"] = num_format_str

        return props


class XlsxWriter(ExcelWriter):
    _engine = "xlsxwriter"
    _supported_extensions = (".xlsx",)

    def __init__(  # pyright: ignore[reportInconsistentConstructor]
        self,
        path: FilePath | WriteExcelBuffer | ExcelWriter,
        engine: str | None = None,
        date_format: str | None = None,
        datetime_format: str | None = None,
        mode: str = "w",
        storage_options: StorageOptions | None = None,
        if_sheet_exists: str | None = None,
        engine_kwargs: dict | None = None,
        autofilter: bool = False,
    ) -> None:
        # Use the xlsxwriter module as the Excel writer.
        import_optional_dependency("xlsxwriter")
        # xlsxwriter does not support append; raise before delegating to
        # base init which rewrites mode
        if "a" in (mode or ""):
            raise ValueError("Append mode is not supported with xlsxwriter!")
        super().__init__(
            path,
            mode=mode,
            storage_options=storage_options,
            if_sheet_exists=if_sheet_exists,
            engine_kwargs=engine_kwargs,
        )

        self._engine_kwargs = engine_kwargs or {}
        self.autofilter = autofilter
        self._book = None
        # Let xlsxwriter raise its own TypeError to satisfy tests
        # expecting that error
        self._book = Workbook(self._handles.handle, **self._engine_kwargs)  # type: ignore[arg-type]

    @property
    def book(self):
        """
        Book instance of class xlsxwriter.Workbook.

        This attribute can be used to access engine-specific features.
        """
        return self._book

    @property
    def sheets(self) -> dict[str, Any]:
        result = self.book.sheetnames
        return result

    def _save(self) -> None:
        """
        Save workbook to disk.
        """
        self.book.close()

    def _write_cells(
        self,
        cells,
        sheet_name: str | None = None,
        startrow: int = 0,
        startcol: int = 0,
        freeze_panes: tuple[int, int] | None = None,
    ) -> None:
        # Write the frame cells using xlsxwriter.
        sheet_name = self._get_sheet_name(sheet_name)

        wks = self.book.get_worksheet_by_name(sheet_name)
        if wks is None:
            wks = self.book.add_worksheet(sheet_name)

        style_dict = {"null": None}

        if validate_freeze_panes(freeze_panes):
            wks.freeze_panes(*(freeze_panes))

        # Initialize bounds with first cell
        first_cell = next(cells, None)
        if first_cell is None:
            return

        # Initialize with first cell's position
        min_row = startrow + first_cell.row
        min_col = startcol + first_cell.col
        max_row = min_row
        max_col = min_col

        # Process first cell
        val, fmt = self._value_with_fmt(first_cell.val)
        stylekey = json.dumps(first_cell.style)
        if fmt:
            stylekey += fmt

        if stylekey in style_dict:
            style = style_dict[stylekey]
        else:
            style = self.book.add_format(_XlsxStyler.convert(first_cell.style, fmt))
            style_dict[stylekey] = style

        wks.write(startrow + first_cell.row, startcol + first_cell.col, val, style)

        # Process remaining cells
        for cell in cells:
            val, fmt = self._value_with_fmt(cell.val)

            stylekey = json.dumps(cell.style)
            if fmt:
                stylekey += fmt

            if stylekey in style_dict:
                style = style_dict[stylekey]
            else:
                style = self.book.add_format(_XlsxStyler.convert(cell.style, fmt))
                style_dict[stylekey] = style

            row = startrow + cell.row
            col = startcol + cell.col

            # Write the cell
            wks.write(row, col, val, style)

            # Update bounds
            min_row = min(min_row, row) if min_row is not None else row
            min_col = min(min_col, col) if min_col is not None else col
            max_row = max(max_row, row) if max_row is not None else row
            max_col = max(max_col, col) if max_col is not None else col

            if cell.mergestart is not None and cell.mergeend is not None:
                wks.merge_range(
                    row,
                    col,
                    startrow + cell.mergestart,
                    startcol + cell.mergeend,
                    val,
                    style,
                )
            else:
                wks.write(row, col, val, style)

        # Apply autofilter if requested
        if getattr(self, "autofilter", False):
            wks.autofilter(min_row, min_col, max_row, max_col)

        if hasattr(self, "_engine_kwargs") and bool(
            self._engine_kwargs.get("autofilter_header", False)
        ):
            if (
                min_row is not None
                and min_col is not None
                and max_row is not None
                and max_col is not None
            ):
                try:
                    wks.autofilter(min_row, min_col, max_row, max_col)
                except Exception:
                    pass
