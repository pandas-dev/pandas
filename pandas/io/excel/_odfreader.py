import pandas
from pandas.io.parsers import TextParser


class _ODFReader(object):
    """Read tables out of OpenDocument formatted files

    Parameters
    ----------
    filepath_or_buffer: string, path to be parsed or
        an open readable stream.
    """
    def __init__(self, filepath_or_buffer):
        try:
            from odf.opendocument import load as document_load
            from odf.table import Table
        except ImportError:
            raise ImportError("Install odfpy for OpenDocument support")

        self.filepath_or_buffer = filepath_or_buffer
        self.document = document_load(filepath_or_buffer)
        self.tables = self.document.getElementsByType(Table)

    @property
    def sheet_names(self):
        """Return table names is the document"""
        from odf.namespaces import TABLENS
        return [t.attributes[(TABLENS, 'name')] for t in self.tables]

    def get_sheet_by_index(self, index):
        return self.__get_table(self.tables[index])

    def get_sheet_by_name(self, name):
        i = self.sheet_names.index(name)
        return self.__get_table(self.tables[i])

    def get_sheet(self, name):
        """Given a sheet name or index, return the root ODF Table node
        """
        if isinstance(name, str):
            return self.get_sheet_by_name(name)
        elif isinstance(name, int):
            return self.get_sheet_by_index(name)
        else:
            raise ValueError(
                'Unrecognized sheet identifier type {}. Please use'
                'a string or integer'.format(type(name)))

    def parse(self, sheet_name=0, **kwds):
        data = self.get_sheet(sheet_name)
        parser = TextParser(data, **kwds)
        return parser.read()

    def __get_table(self, sheet):
        """Parse an ODF Table into a list of lists
        """
        from odf.table import TableCell, TableRow

        sheet_rows = sheet.getElementsByType(TableRow)
        table = []
        empty_rows = 0
        max_row_len = 0
        for i, sheet_row in enumerate(sheet_rows):
            sheet_cells = sheet_row.getElementsByType(TableCell)
            empty_cells = 0
            table_row = []
            for j, sheet_cell in enumerate(sheet_cells):
                value = self.__get_cell_value(sheet_cell)
                column_repeat = self.__get_cell_repeat(sheet_cell)

                if len(sheet_cell.childNodes) == 0:
                    empty_cells += column_repeat
                else:
                    if empty_cells > 0:
                        table_row.extend([None] * empty_cells)
                        empty_cells = 0
                    table_row.extend([value] * column_repeat)

            if max_row_len < len(table_row):
                max_row_len = len(table_row)

            row_repeat = self.__get_row_repeat(sheet_row)
            if self.__is_empty_row(sheet_row):
                empty_rows += row_repeat
            else:
                if empty_rows > 0:
                    # add blank rows to our table
                    table.extend([[None]] * empty_rows)
                    empty_rows = 0
                table.append(table_row)

        # Make our table square
        for row in table:
            if len(row) < max_row_len:
                row.extend([None] * (max_row_len - len(row)))

        return table

    def __get_row_repeat(self, row):
        """Return number of times this row was repeated

        Repeating an empty row appeared to be a common way
        of representing sparse rows in the table.
        """
        from odf.namespaces import TABLENS
        repeat = row.attributes.get((TABLENS, 'number-rows-repeated'))
        if repeat is None:
            return 1
        return int(repeat)

    def __get_cell_repeat(self, cell):
        from odf.namespaces import TABLENS
        repeat = cell.attributes.get((TABLENS, 'number-columns-repeated'))
        if repeat is None:
            return 1
        return int(repeat)

    def __is_empty_row(self, row):
        """Helper function to find empty rows
        """
        for column in row.childNodes:
            if len(column.childNodes) > 0:
                return False

        return True

    def __get_cell_value(self, cell):
        from odf.namespaces import OFFICENS
        cell_type = cell.attributes.get((OFFICENS, 'value-type'))
        if cell_type == 'boolean':
            cell_value = cell.attributes.get((OFFICENS, 'boolean'))
            return bool(cell_value)
        elif cell_type in ('float', 'percentage'):
            cell_value = cell.attributes.get((OFFICENS, 'value'))
            return float(cell_value)
        elif cell_type == 'string':
            return str(cell)
        elif cell_type == 'currency':
            cell_value = cell.attributes.get((OFFICENS, 'value'))
            return float(cell_value)
        elif cell_type == 'date':
            cell_value = cell.attributes.get((OFFICENS, 'date-value'))
            return pandas.Timestamp(cell_value)
        elif cell_type == 'time':
            cell_value = cell.attributes.get((OFFICENS, 'time-value'))
            return(pandas_isoduration_compatibility(cell_value))
        elif cell_type is None:
            return None
        else:
            raise ValueError('Unrecognized type {}'.format(cell_type))


def pandas_isoduration_compatibility(duration):
    """Libreoffice returns durations without any day attributes

    For example PT3H45M0S. The current pandas Timedelta
    parse requires the presence of a day component.
    Workaround for https://github.com/pandas-dev/pandas/issues/25422
    """
    if duration.startswith('PT'):
        duration = 'P0DT' + duration[2:]
    return pandas.Timedelta(duration)
