from typing import List

import pandas as pd

from pandas._typing import FilePathOrBuffer, Scalar

from pandas.compat._optional import import_optional_dependency

from pandas.io.excel._base import _BaseExcelReader


class _ODFReader(_BaseExcelReader):
    """Read tables out of OpenDocument formatted files

    Parameters
    ----------
    filepath_or_buffer: string, path to be parsed or
        an open readable stream.
    """
    def __init__(self, filepath_or_buffer: FilePathOrBuffer):
        import_optional_dependency("odf")
        super().__init__(filepath_or_buffer)

    @property
    def _workbook_class(self):
        from odf.opendocument import OpenDocument
        return OpenDocument

    def load_workbook(self, filepath_or_buffer: FilePathOrBuffer):
        from odf.opendocument import load
        return load(filepath_or_buffer)

    @property
    def sheet_names(self) -> List[str]:
        """Return a list of sheet names present in the document"""
        from odf.namespaces import TABLENS
        from odf.table import Table

        tables = self.book.getElementsByType(Table)
        return [t.attributes[(TABLENS, 'name')] for t in tables]

    def get_sheet_by_index(self, index: int):
        from odf.table import Table
        tables = self.book.getElementsByType(Table)
        return tables[index]

    def get_sheet_by_name(self, name: str):
        from odf.namespaces import TABLENS
        from odf.table import Table

        tables = self.book.getElementsByType(Table)

        key = (TABLENS, "name")
        for table in tables:
            if table.attributes[key] == name:
                return table

        raise ValueError("sheet {name} not found".format(name))

    def get_sheet_data(self, sheet, convert_float: bool) -> List[List[Scalar]]:
        """Parse an ODF Table into a list of lists
        """
        from odf.table import TableCell, TableRow

        sheet_rows = sheet.getElementsByType(TableRow)
        table = []
        empty_rows = 0
        max_row_len = 0
        row_spans = {}  # type: Dict[int, int]

        for i, sheet_row in enumerate(sheet_rows):
            sheet_cells = sheet_row.getElementsByType(TableCell)
            empty_cells = 0
            table_row = []

            for j, sheet_cell in enumerate(sheet_cells):
                # Handle vertically merged cells; only works with first column
                if row_spans.get(j, 0) > 1:
                    table_row.append('')
                    row_spans[j] = row_spans[j] - 1

                value = self._get_cell_value(sheet_cell, convert_float)
                column_repeat = self._get_column_repeat(sheet_cell)
                column_span = self._get_column_span(sheet_cell)
                row_span = self._get_row_span(sheet_cell)

                if row_span > 1:
                    if j > 0:
                        raise ValueError(
                            "The odf reader only supports vertical "
                            "merging in the initial column")
                    else:
                        row_spans[j] = row_span

                if len(sheet_cell.childNodes) == 0:
                    empty_cells += column_repeat
                else:
                    if empty_cells > 0:
                        table_row.extend([''] * empty_cells)
                        empty_cells = 0
                    table_row.extend([value] * column_repeat)

                    # horizontally merged cells should only show first value
                    if column_span > 1:
                        table_row.extend([''] * (column_span - 1))

            if max_row_len < len(table_row):
                max_row_len = len(table_row)

            row_repeat = self._get_row_repeat(sheet_row)
            if self._is_empty_row(sheet_row):
                empty_rows += row_repeat
            else:
                if empty_rows > 0:
                    # add blank rows to our table
                    table.extend([['']] * empty_rows)
                    empty_rows = 0
                for _ in range(row_repeat):
                    table.append(table_row)

        # Make our table square
        for row in table:
            if len(row) < max_row_len:
                row.extend([''] * (max_row_len - len(row)))

        return table

    def _get_row_repeat(self, row):
        """Return number of times this row was repeated
        Repeating an empty row appeared to be a common way
        of representing sparse rows in the table.
        """
        from odf.namespaces import TABLENS
        repeat = row.attributes.get((TABLENS, 'number-rows-repeated'))
        if repeat is None:
            return 1
        return int(repeat)

    def _get_column_repeat(self, cell):
        from odf.namespaces import TABLENS
        repeat = cell.attributes.get((TABLENS, 'number-columns-repeated'))
        if repeat is None:
            return 1
        return int(repeat)

    def _get_row_span(self, cell):
        """For handling cells merged vertically."""
        from odf.namespaces import TABLENS
        repeat = cell.attributes.get((TABLENS, 'number-rows-spanned'))
        if repeat is None:
            return 1
        return int(repeat)

    def _get_column_span(self, cell):
        """For handling cells merged horizontally."""
        from odf.namespaces import TABLENS
        repeat = cell.attributes.get((TABLENS, 'number-columns-spanned'))
        if repeat is None:
            return 1
        return int(repeat)

    def _is_empty_row(self, row):
        """Helper function to find empty rows
        """
        for column in row.childNodes:
            if len(column.childNodes) > 0:
                return False

        return True

    def _get_cell_value(self, cell, convert_float: bool) -> Scalar:
        from odf.namespaces import OFFICENS
        cell_type = cell.attributes.get((OFFICENS, 'value-type'))
        if cell_type == 'boolean':
            if str(cell) == "TRUE":
                return True
            return False
        if cell_type is None:
            return ''  # compat with xlrd
        elif cell_type == 'float':
            # GH5394
            cell_value = float(cell.attributes.get((OFFICENS, 'value')))

            if cell_value == 0. and str(cell) != cell_value:  # NA handling
                return str(cell)

            if convert_float:
                val = int(cell_value)
                if val == cell_value:
                    return val
            return cell_value
        elif cell_type == 'percentage':
            cell_value = cell.attributes.get((OFFICENS, 'value'))
            return float(cell_value)
        elif cell_type == 'string':
            return str(cell)
        elif cell_type == 'currency':
            cell_value = cell.attributes.get((OFFICENS, 'value'))
            return float(cell_value)
        elif cell_type == 'date':
            cell_value = cell.attributes.get((OFFICENS, 'date-value'))
            return pd.to_datetime(cell_value)
        elif cell_type == 'time':
            return pd.to_datetime(str(cell)).time()
        else:
            raise ValueError('Unrecognized type {}'.format(cell_type))
