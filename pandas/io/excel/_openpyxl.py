from collections import OrderedDict
from io import BytesIO

import pandas.compat as compat
from pandas.core.dtypes.common import is_integer, is_list_like
from pandas.core.frame import DataFrame
from pandas.io.common import (_is_url, _urlopen, _validate_header_arg,
                              get_filepath_or_buffer)
from pandas.io.excel._base import (ExcelFile, ExcelWriter, _BaseExcelReader,
                                   _fill_mi_header, _maybe_convert_to_string,
                                   _maybe_convert_usecols, _pop_header_name)
from pandas.io.excel._util import _validate_freeze_panes
from pandas.io.parsers import _validate_usecols_arg, _validate_usecols_names


class _OpenpyxlWriter(ExcelWriter):
    engine = 'openpyxl'
    supported_extensions = ('.xlsx', '.xlsm')

    def __init__(self, path, engine=None, mode='w', **engine_kwargs):
        # Use the openpyxl module as the Excel writer.
        from openpyxl.workbook import Workbook

        super(_OpenpyxlWriter, self).__init__(path, mode=mode, **engine_kwargs)

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
        err_msg = "Install xlrd >= 1.0.0 for Excel support"

        try:
            import openpyxl
        except ImportError:
            raise ImportError(err_msg)

        # If filepath_or_buffer is a url, want to keep the data as bytes so
        # can't pass to get_filepath_or_buffer()
        if _is_url(filepath_or_buffer):
            filepath_or_buffer = BytesIO(_urlopen(filepath_or_buffer).read())
        elif not isinstance(filepath_or_buffer,
                            (ExcelFile, openpyxl.Workbook)):
            filepath_or_buffer, _, _, _ = get_filepath_or_buffer(
                filepath_or_buffer)

        if isinstance(filepath_or_buffer, openpyxl.Workbook):
            self.book = filepath_or_buffer
        elif hasattr(filepath_or_buffer, "read"):
            if hasattr(filepath_or_buffer, 'seek'):
                filepath_or_buffer.seek(0)
            self.book = openpyxl.load_workbook(
                filepath_or_buffer, data_only=True)
        elif isinstance(filepath_or_buffer, compat.string_types):
            self.book = openpyxl.load_workbook(
                filepath_or_buffer, data_only=True)
        else:
            raise ValueError('Must explicitly set engine if not passing in'
                             ' buffer or path for io.')

    @property
    def sheet_names(self):
        return self.book.sheetnames

    def get_sheet_by_name(self, name):
        return self.book[name]

    def get_sheet_by_index(self, index):
        return self.book.worksheets[index]

    @staticmethod
    def _replace_type_error_with_nan(rows):
        nan = float('nan')
        for row in rows:
            yield [nan
                   if cell.data_type == cell.TYPE_ERROR
                   else cell.value
                   for cell in row]

    def get_sheet_data(self, sheet, convert_float):
        data = self._replace_type_error_with_nan(sheet.rows)
        # TODO: support using iterator
        # TODO: don't make strings out of data
        return list(data)

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

        ret_dict = False

        # Keep sheetname to maintain backwards compatibility.
        if isinstance(sheet_name, list):
            sheets = sheet_name
            ret_dict = True
        elif sheet_name is None:
            sheets = self.sheet_names
            ret_dict = True
        else:
            sheets = [sheet_name]

        # handle same-type duplicates.
        sheets = list(OrderedDict.fromkeys(sheets).keys())

        output = OrderedDict()

        for asheetname in sheets:
            if verbose:
                print("Reading sheet {sheet}".format(sheet=asheetname))

            if isinstance(asheetname, compat.string_types):
                sheet = self.get_sheet_by_name(asheetname)
            else:  # assume an integer if not a string
                sheet = self.get_sheet_by_index(asheetname)

            data = self.get_sheet_data(sheet, convert_float)
            if not data or data == [[None]]:
                output[asheetname] = DataFrame()
                continue

            usecols = _maybe_convert_usecols(usecols)

            if is_list_like(header) and len(header) == 1:
                header = header[0]

            # TODO: scrutinize what is going here
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

            # TODO: implement whatever this should do
            # has_index_names = is_list_like(header) and len(header) > 1

            if skiprows:
                data = [row for i, row in enumerate(data) if i not in skiprows]

            if skipfooter:
                data = data[:-skipfooter]

            column_names = [cell for i, cell in enumerate(data.pop(0))]

            frame = DataFrame(data, columns=column_names)
            if usecols:
                _validate_usecols_arg(usecols)
                usecols = sorted(usecols)
                if any(isinstance(i, str) for i in usecols):
                    _validate_usecols_names(usecols, column_names)
                    frame = frame[usecols]
                else:
                    frame = frame.iloc[:, usecols]

            if not converters:
                converters = dict()
            if not dtype:
                dtype = dict()

            # handle columns referenced by number so all references are by
            #  column name
            handled_converters = {}
            for k, v in converters.items():
                if k not in frame.columns and isinstance(k, int):
                    k = frame.columns[k]
                handled_converters[k] = v
            converters = handled_converters

            # attempt to convert object columns to integer. Only because this
            # is implicitly done when reading and excel file with xlrd
            # TODO: question if this should be default behaviour
            if len(frame) > 0:
                for column in set(frame) - set(dtype.keys()):
                    if frame[column].dtype == object:
                        try:
                            frame[column] = frame[column].astype('int64')
                        except (ValueError, TypeError):
                            try:
                                frame[column] = frame[column].astype('float64')
                            except (ValueError, TypeError):
                                continue
                    elif (convert_float and
                            frame[column].dtype == float and
                            all(frame[column] % 1 == 0)):
                        frame[column] = frame[column].astype('int64')
                    elif not convert_float:
                        if frame[column].dtype == int:
                            frame[column] = frame[column].astype('float64')

            if converters:
                for k, v in converters.items():
                    # for compatibiliy reasons
                    if frame[k].dtype == float and convert_float:
                        frame[k] = frame[k].fillna('')
                    frame[k] = frame[k].apply(v)

            if dtype:
                for k, v in dtype.items():
                    frame[k] = frame[k].astype(v)

            if index_col is not None:
                if is_list_like(index_col):
                    if any(isinstance(i, str) for i in index_col):
                        # TODO: see if there is already a method for this in
                        # pandas.io.parsers
                        frame = frame.set_index(index_col)
                        if len(index_col) == 1:
                            # TODO: understand why this is needed
                            raise TypeError(
                                "list indices must be integers.*, not str")
                    else:
                        frame = frame.set_index(
                            [column_names[i] for i in index_col])
                else:
                    if isinstance(index_col, str):
                        frame = frame.set_index(index_col)
                    else:
                        frame = frame.set_index(column_names[index_col])

            output[asheetname] = frame
            if not squeeze or isinstance(output[asheetname], DataFrame):
                if header_names:
                    output[asheetname].columns = output[
                        asheetname].columns.set_names(header_names)
                elif compat.PY2:
                    output[asheetname].columns = _maybe_convert_to_string(
                        output[asheetname].columns)

            # name unnamed columns
            unnamed = 0
            for i, col_name in enumerate(frame.columns.values):
                if col_name is None:
                    frame.columns.values[i] = "Unnamed: {n}".format(n=unnamed)
                    unnamed += 1

        if ret_dict:
            return output
        else:
            return output[asheetname]
