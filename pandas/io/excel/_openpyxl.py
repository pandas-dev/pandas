from collections import OrderedDict
from distutils.version import LooseVersion

from pandas.compat._optional import import_optional_dependency
from pandas.core.dtypes.common import (
    ensure_int_or_float, is_float_dtype, is_integer, is_integer_dtype,
    is_list_like, is_object_dtype)

from pandas.core.internals.construction import get_names_from_index
from pandas.core.frame import DataFrame

from pandas.io.common import _validate_header_arg
from pandas.io.excel._base import (
    ExcelWriter, _BaseExcelReader, _fill_mi_header, _maybe_convert_usecols,
    _pop_header_name)
from pandas.io.excel._util import _validate_freeze_panes
from pandas.io.parsers import _validate_usecols_arg, _validate_usecols_names


class _OpenpyxlWriter(ExcelWriter):
    engine = 'openpyxl'
    supported_extensions = ('.xlsx', '.xlsm')

    def __init__(self, path, engine=None, mode='w', **engine_kwargs):
        # Use the openpyxl module as the Excel writer.
        from openpyxl.workbook import Workbook

        super().__init__(path, mode=mode, **engine_kwargs)

        if self.mode == 'a':  # Load from existing workbook
            from openpyxl import load_workbook
            book = load_workbook(self.path)
            self.book = book
        else:
            # Create workbook object with default optimized_write=True.
            self.book = Workbook()

            if self.book.worksheets:
                try:
                    self.book.remove(self.book.worksheets[0])
                except AttributeError:

                    # compat - for openpyxl <= 2.4
                    self.book.remove_sheet(self.book.worksheets[0])

    def save(self):
        """
        Save workbook to disk.
        """
        return self.book.save(self.path)

    @classmethod
    def _convert_to_style(cls, style_dict):
        """
        converts a style_dict to an openpyxl style object
        Parameters
        ----------
        style_dict : style dictionary to convert
        """

        from openpyxl.style import Style
        xls_style = Style()
        for key, value in style_dict.items():
            for nk, nv in value.items():
                if key == "borders":
                    (xls_style.borders.__getattribute__(nk)
                     .__setattr__('border_style', nv))
                else:
                    xls_style.__getattribute__(key).__setattr__(nk, nv)

        return xls_style

    @classmethod
    def _convert_to_style_kwargs(cls, style_dict):
        """
        Convert a style_dict to a set of kwargs suitable for initializing
        or updating-on-copy an openpyxl v2 style object
        Parameters
        ----------
        style_dict : dict
            A dict with zero or more of the following keys (or their synonyms).
                'font'
                'fill'
                'border' ('borders')
                'alignment'
                'number_format'
                'protection'
        Returns
        -------
        style_kwargs : dict
            A dict with the same, normalized keys as ``style_dict`` but each
            value has been replaced with a native openpyxl style object of the
            appropriate class.
        """

        _style_key_map = {
            'borders': 'border',
        }

        style_kwargs = {}
        for k, v in style_dict.items():
            if k in _style_key_map:
                k = _style_key_map[k]
            _conv_to_x = getattr(cls, '_convert_to_{k}'.format(k=k),
                                 lambda x: None)
            new_v = _conv_to_x(v)
            if new_v:
                style_kwargs[k] = new_v

        return style_kwargs

    @classmethod
    def _convert_to_color(cls, color_spec):
        """
        Convert ``color_spec`` to an openpyxl v2 Color object
        Parameters
        ----------
        color_spec : str, dict
            A 32-bit ARGB hex string, or a dict with zero or more of the
            following keys.
                'rgb'
                'indexed'
                'auto'
                'theme'
                'tint'
                'index'
                'type'
        Returns
        -------
        color : openpyxl.styles.Color
        """

        from openpyxl.styles import Color

        if isinstance(color_spec, str):
            return Color(color_spec)
        else:
            return Color(**color_spec)

    @classmethod
    def _convert_to_font(cls, font_dict):
        """
        Convert ``font_dict`` to an openpyxl v2 Font object
        Parameters
        ----------
        font_dict : dict
            A dict with zero or more of the following keys (or their synonyms).
                'name'
                'size' ('sz')
                'bold' ('b')
                'italic' ('i')
                'underline' ('u')
                'strikethrough' ('strike')
                'color'
                'vertAlign' ('vertalign')
                'charset'
                'scheme'
                'family'
                'outline'
                'shadow'
                'condense'
        Returns
        -------
        font : openpyxl.styles.Font
        """

        from openpyxl.styles import Font

        _font_key_map = {
            'sz': 'size',
            'b': 'bold',
            'i': 'italic',
            'u': 'underline',
            'strike': 'strikethrough',
            'vertalign': 'vertAlign',
        }

        font_kwargs = {}
        for k, v in font_dict.items():
            if k in _font_key_map:
                k = _font_key_map[k]
            if k == 'color':
                v = cls._convert_to_color(v)
            font_kwargs[k] = v

        return Font(**font_kwargs)

    @classmethod
    def _convert_to_stop(cls, stop_seq):
        """
        Convert ``stop_seq`` to a list of openpyxl v2 Color objects,
        suitable for initializing the ``GradientFill`` ``stop`` parameter.
        Parameters
        ----------
        stop_seq : iterable
            An iterable that yields objects suitable for consumption by
            ``_convert_to_color``.
        Returns
        -------
        stop : list of openpyxl.styles.Color
        """

        return map(cls._convert_to_color, stop_seq)

    @classmethod
    def _convert_to_fill(cls, fill_dict):
        """
        Convert ``fill_dict`` to an openpyxl v2 Fill object
        Parameters
        ----------
        fill_dict : dict
            A dict with one or more of the following keys (or their synonyms),
                'fill_type' ('patternType', 'patterntype')
                'start_color' ('fgColor', 'fgcolor')
                'end_color' ('bgColor', 'bgcolor')
            or one or more of the following keys (or their synonyms).
                'type' ('fill_type')
                'degree'
                'left'
                'right'
                'top'
                'bottom'
                'stop'
        Returns
        -------
        fill : openpyxl.styles.Fill
        """

        from openpyxl.styles import PatternFill, GradientFill

        _pattern_fill_key_map = {
            'patternType': 'fill_type',
            'patterntype': 'fill_type',
            'fgColor': 'start_color',
            'fgcolor': 'start_color',
            'bgColor': 'end_color',
            'bgcolor': 'end_color',
        }

        _gradient_fill_key_map = {
            'fill_type': 'type',
        }

        pfill_kwargs = {}
        gfill_kwargs = {}
        for k, v in fill_dict.items():
            pk = gk = None
            if k in _pattern_fill_key_map:
                pk = _pattern_fill_key_map[k]
            if k in _gradient_fill_key_map:
                gk = _gradient_fill_key_map[k]
            if pk in ['start_color', 'end_color']:
                v = cls._convert_to_color(v)
            if gk == 'stop':
                v = cls._convert_to_stop(v)
            if pk:
                pfill_kwargs[pk] = v
            elif gk:
                gfill_kwargs[gk] = v
            else:
                pfill_kwargs[k] = v
                gfill_kwargs[k] = v

        try:
            return PatternFill(**pfill_kwargs)
        except TypeError:
            return GradientFill(**gfill_kwargs)

    @classmethod
    def _convert_to_side(cls, side_spec):
        """
        Convert ``side_spec`` to an openpyxl v2 Side object
        Parameters
        ----------
        side_spec : str, dict
            A string specifying the border style, or a dict with zero or more
            of the following keys (or their synonyms).
                'style' ('border_style')
                'color'
        Returns
        -------
        side : openpyxl.styles.Side
        """

        from openpyxl.styles import Side

        _side_key_map = {
            'border_style': 'style',
        }

        if isinstance(side_spec, str):
            return Side(style=side_spec)

        side_kwargs = {}
        for k, v in side_spec.items():
            if k in _side_key_map:
                k = _side_key_map[k]
            if k == 'color':
                v = cls._convert_to_color(v)
            side_kwargs[k] = v

        return Side(**side_kwargs)

    @classmethod
    def _convert_to_border(cls, border_dict):
        """
        Convert ``border_dict`` to an openpyxl v2 Border object
        Parameters
        ----------
        border_dict : dict
            A dict with zero or more of the following keys (or their synonyms).
                'left'
                'right'
                'top'
                'bottom'
                'diagonal'
                'diagonal_direction'
                'vertical'
                'horizontal'
                'diagonalUp' ('diagonalup')
                'diagonalDown' ('diagonaldown')
                'outline'
        Returns
        -------
        border : openpyxl.styles.Border
        """

        from openpyxl.styles import Border

        _border_key_map = {
            'diagonalup': 'diagonalUp',
            'diagonaldown': 'diagonalDown',
        }

        border_kwargs = {}
        for k, v in border_dict.items():
            if k in _border_key_map:
                k = _border_key_map[k]
            if k == 'color':
                v = cls._convert_to_color(v)
            if k in ['left', 'right', 'top', 'bottom', 'diagonal']:
                v = cls._convert_to_side(v)
            border_kwargs[k] = v

        return Border(**border_kwargs)

    @classmethod
    def _convert_to_alignment(cls, alignment_dict):
        """
        Convert ``alignment_dict`` to an openpyxl v2 Alignment object
        Parameters
        ----------
        alignment_dict : dict
            A dict with zero or more of the following keys (or their synonyms).
                'horizontal'
                'vertical'
                'text_rotation'
                'wrap_text'
                'shrink_to_fit'
                'indent'
        Returns
        -------
        alignment : openpyxl.styles.Alignment
        """

        from openpyxl.styles import Alignment

        return Alignment(**alignment_dict)

    @classmethod
    def _convert_to_number_format(cls, number_format_dict):
        """
        Convert ``number_format_dict`` to an openpyxl v2.1.0 number format
        initializer.
        Parameters
        ----------
        number_format_dict : dict
            A dict with zero or more of the following keys.
                'format_code' : str
        Returns
        -------
        number_format : str
        """
        return number_format_dict['format_code']

    @classmethod
    def _convert_to_protection(cls, protection_dict):
        """
        Convert ``protection_dict`` to an openpyxl v2 Protection object.
        Parameters
        ----------
        protection_dict : dict
            A dict with zero or more of the following keys.
                'locked'
                'hidden'
        Returns
        -------
        """

        from openpyxl.styles import Protection

        return Protection(**protection_dict)

    def write_cells(self, cells, sheet_name=None, startrow=0, startcol=0,
                    freeze_panes=None):
        # Write the frame cells using openpyxl.
        sheet_name = self._get_sheet_name(sheet_name)

        _style_cache = {}

        if sheet_name in self.sheets:
            wks = self.sheets[sheet_name]
        else:
            wks = self.book.create_sheet()
            wks.title = sheet_name
            self.sheets[sheet_name] = wks

        if _validate_freeze_panes(freeze_panes):
            wks.freeze_panes = wks.cell(row=freeze_panes[0] + 1,
                                        column=freeze_panes[1] + 1)

        for cell in cells:
            xcell = wks.cell(
                row=startrow + cell.row + 1,
                column=startcol + cell.col + 1
            )
            xcell.value, fmt = self._value_with_fmt(cell.val)
            if fmt:
                xcell.number_format = fmt

            style_kwargs = {}
            if cell.style:
                key = str(cell.style)
                style_kwargs = _style_cache.get(key)
                if style_kwargs is None:
                    style_kwargs = self._convert_to_style_kwargs(cell.style)
                    _style_cache[key] = style_kwargs

            if style_kwargs:
                for k, v in style_kwargs.items():
                    setattr(xcell, k, v)

            if cell.mergestart is not None and cell.mergeend is not None:

                wks.merge_cells(
                    start_row=startrow + cell.row + 1,
                    start_column=startcol + cell.col + 1,
                    end_column=startcol + cell.mergeend + 1,
                    end_row=startrow + cell.mergestart + 1
                )

                # When cells are merged only the top-left cell is preserved
                # The behaviour of the other cells in a merged range is
                # undefined
                if style_kwargs:
                    first_row = startrow + cell.row + 1
                    last_row = startrow + cell.mergestart + 1
                    first_col = startcol + cell.col + 1
                    last_col = startcol + cell.mergeend + 1

                    for row in range(first_row, last_row + 1):
                        for col in range(first_col, last_col + 1):
                            if row == first_row and col == first_col:
                                # Ignore first cell. It is already handled.
                                continue
                            xcell = wks.cell(column=col, row=row)
                            for k, v in style_kwargs.items():
                                setattr(xcell, k, v)


class _OpenpyxlReader(_BaseExcelReader):

    def __init__(self, filepath_or_buffer):
        """Reader using openpyxl engine.

        Parameters
        ----------
        filepath_or_buffer : string, path object or Workbook
            Object to be parsed.
        """
        import_optional_dependency("openpyxl")
        super().__init__(filepath_or_buffer)

    @property
    def _workbook_class(self):
        from openpyxl import Workbook
        return Workbook

    def load_workbook(self, filepath_or_buffer):
        from openpyxl import load_workbook
        return load_workbook(filepath_or_buffer, data_only=True)

    @property
    def sheet_names(self):
        return self.book.sheetnames

    @staticmethod
    def _handle_usecols(frame, usecols):
        column_names = frame.columns.values
        if usecols:
            _validate_usecols_arg(usecols)
            usecols = sorted(usecols)
            if any(isinstance(i, str) for i in usecols):
                _validate_usecols_names(usecols, column_names)
                frame = frame[usecols]
            else:
                frame = frame.iloc[:, usecols]
        return frame

    def _handle_sheet_name(self, sheet_name):
        """Handle the sheet_name keyword."""
        # Keep sheetname to maintain backwards compatibility.
        if isinstance(sheet_name, list):
            sheets = sheet_name
        elif sheet_name is None:
            sheets = self.sheet_names
        else:
            sheets = [sheet_name]
        return sheets

    @staticmethod
    def _handle_header_keywords(data, header, skiprows, index_col):
        """Handle keywords relating to header parsing."""
        # forward fill and pull out names for MultiIndex column
        header_names = None
        if header is not None and is_list_like(header):
            header_names = []
            control_row = [True] * len(data[0])

            for row in header:
                if is_integer(skiprows):
                    row += skiprows

                data[row], control_row = _fill_mi_header(data[row],
                                                         control_row)

                if index_col is not None:
                    header_name, _ = _pop_header_name(data[row], index_col)
                    header_names.append(header_name)
        return header_names

    @staticmethod
    def _handle_convert_float(series, convert_float):
        """Handle the convert_float keyword."""
        # attempt to convert object columns to integer. Only because this
        # is implicitly done when reading and excel file with xlrd, that
        # behaviour is replicated here.

        if is_object_dtype(series):
            try:
                series = ensure_int_or_float(series)
            except (ValueError):
                return series
        elif (convert_float
                and is_float_dtype(series)
                and all(series % 1 == 0)):
            series = series.astype('int64')
        elif not convert_float:
            if is_integer_dtype(series):
                series = series.astype('float64')
        return series

    @staticmethod
    def _handle_index_col(frame, index_col):
        column_names = frame.columns.values
        if index_col is None:
            return frame
        if is_list_like(index_col):
            if any(isinstance(i, str) for i in index_col):
                # TODO: see if there is already a method for this in
                # pandas.io.parsers
                frame = frame.set_index(index_col)
                if len(index_col) == 1:
                    # TODO: understand why this is needed
                    raise TypeError("list indices must be integers.*, not str")
            else:
                frame = frame.set_index([column_names[i] for i in index_col])
        else:
            if isinstance(index_col, str):
                frame = frame.set_index(index_col)
            else:
                frame = frame.set_index(column_names[index_col])
        return frame

    def get_sheet_by_name(self, name):
        return self.book[name]

    def get_sheet_by_index(self, index):
        return self.book.worksheets[index]

    def _replace_type_error_with_nan(self, rows):
        try:
            from openpyxl.cell.cell import TYPE_ERROR
        except ImportError:  # openpyxl < 2.6
            from openpyxl.cell.cell import Cell
            TYPE_ERROR = Cell.TYPE_ERROR

        for row in rows:
            return [np.nan
                   if cell.data_type == TYPE_ERROR
                   else cell.value
                   for cell in row]

    def get_sheet_data(self, sheet):
        data = self._replace_type_error_with_nan(sheet.rows)
        return list(data)

    def _parse_sheet(self, sheet, convert_float, usecols, header, skiprows,
                     index_col, converters, skipfooter, dtype, squeeze):
        """Parse a single sheet into a dataframe."""

        data = self.get_sheet_data(sheet)
        if not data or data == [[None]]:
            return DataFrame()

        usecols = _maybe_convert_usecols(usecols)

        if is_list_like(header) and len(header) == 1:
            header = header[0]

        header_names = self._handle_header_keywords(data, header, skiprows,
                                                    index_col)

        # TODO: implement whatever this should do
        # has_index_names = is_list_like(header) and len(header) > 1

        if skiprows:
            data = [row for i, row in enumerate(data) if i not in skiprows]

        if skipfooter:
            data = data[:-skipfooter]

        column_names = [cell for cell in data.pop(0)]

        frame = DataFrame(data, columns=column_names)
        frame = self._handle_usecols(frame, usecols)

        if not converters:
            converters = dict()
        if not dtype:
            dtype = dict()

        # handle columns referenced by number so all references are by
        # column name
        handled_converters = {}
        for k, v in converters.items():
            if k not in frame.columns and isinstance(k, int):
                k = frame.columns[k]
            handled_converters[k] = v
        converters = handled_converters

        if len(frame) > 0:
            for column in set(frame) - set(dtype.keys()):
                frame[column] = self._handle_convert_float(frame[column],
                                                           convert_float)

        if converters:
            for k, v in converters.items():
                # for compatibiliy reasons
                if frame[k].dtype == float and convert_float:
                    frame[k] = frame[k].fillna('')
                frame[k] = frame[k].apply(v)

        if dtype:
            for k, v in dtype.items():
                frame[k] = frame[k].astype(v)

        frame = self._handle_index_col(frame, index_col)

        if not squeeze or isinstance(frame, DataFrame):
            if header_names:
                frame = frame.columns.set_names(header_names)

        frame.columns = get_names_from_index(frame.columns)

        return frame

    def parse(self,
              sheet_name=0,
              header=0,
              names=None,
              index_col=None,
              usecols=None,
              squeeze=False,
              converters=None,
              dtype=None,
              true_values=None,
              false_values=None,
              skiprows=None,
              nrows=None,
              na_values=None,
              verbose=False,
              parse_dates=False,
              date_parser=None,
              thousands=None,
              comment=None,
              skipfooter=0,
              convert_float=True,
              mangle_dupe_cols=True,
              **kwds):

        _validate_header_arg(header)

        sheets = self._handle_sheet_name(sheet_name)
        ret_dict = len(sheets) != 1
        output = OrderedDict()

        for asheetname in sheets:
            if verbose:
                print("Reading sheet {sheet}".format(sheet=asheetname))

            if isinstance(asheetname, str):
                sheet = self.get_sheet_by_name(asheetname)
            else:  # assume an integer if not a string
                sheet = self.get_sheet_by_index(asheetname)

            output[asheetname] = self._parse_sheet(
                sheet, convert_float, usecols, header, skiprows, index_col,
                converters, skipfooter, dtype, squeeze)

        if ret_dict:
            return output
        else:
            return output[asheetname]
