# -*- coding: utf-8 -*-
"""
Module for formatting output data in HTML.
"""

from collections import OrderedDict
from textwrap import dedent

from pandas._config import get_option

from pandas.compat import lzip

from pandas.core.dtypes.generic import ABCMultiIndex

from pandas import option_context

from pandas.io.common import _is_url
from pandas.io.formats.format import TableFormatter, get_level_lengths
from pandas.io.formats.printing import pprint_thing


class HTMLFormatter(TableFormatter):
    """
    Internal class for formatting output data in html.
    This class is intended for shared functionality between
    DataFrame.to_html() and DataFrame._repr_html_().
    Any logic in common with other output formatting methods
    should ideally be inherited from classes in format.py
    and this class responsible for only producing html markup.
    """

    indent_delta = 2

    def __init__(self, formatter, classes=None, border=None):
        self.fmt = formatter
        self.classes = classes

        self.frame = self.fmt.frame
        self.columns = self.fmt.tr_frame.columns
        self.elements = []
        self.bold_rows = self.fmt.kwds.get('bold_rows', False)
        self.escape = self.fmt.kwds.get('escape', True)
        self.show_dimensions = self.fmt.show_dimensions
        if border is None:
            border = get_option('display.html.border')
        self.border = border
        self.table_id = self.fmt.table_id
        self.render_links = self.fmt.render_links

    @property
    def show_row_idx_names(self):
        return self.fmt.show_row_idx_names

    @property
    def show_col_idx_names(self):
        return self.fmt.show_col_idx_names

    @property
    def row_levels(self):
        if self.fmt.index:
            # showing (row) index
            return self.frame.index.nlevels
        elif self.show_col_idx_names:
            # see gh-22579
            # Column misalignment also occurs for
            # a standard index when the columns index is named.
            # If the row index is not displayed a column of
            # blank cells need to be included before the DataFrame values.
            return 1
        # not showing (row) index
        return 0

    @property
    def is_truncated(self):
        return self.fmt.is_truncated

    @property
    def ncols(self):
        return len(self.fmt.tr_frame.columns)

    def write(self, s, indent=0):
        rs = pprint_thing(s)
        self.elements.append(' ' * indent + rs)

    def write_th(self, s, indent=0, tags=None):
        if self.fmt.col_space is not None and self.fmt.col_space > 0:
            tags = (tags or "")
            tags += ('style="min-width: {colspace};"'
                     .format(colspace=self.fmt.col_space))

        return self._write_cell(s, kind='th', indent=indent, tags=tags)

    def write_td(self, s, indent=0, tags=None):
        return self._write_cell(s, kind='td', indent=indent, tags=tags)

    def _write_cell(self, s, kind='td', indent=0, tags=None):
        if tags is not None:
            start_tag = '<{kind} {tags}>'.format(kind=kind, tags=tags)
        else:
            start_tag = '<{kind}>'.format(kind=kind)

        if self.escape:
            # escape & first to prevent double escaping of &
            esc = OrderedDict([('&', r'&amp;'), ('<', r'&lt;'),
                               ('>', r'&gt;')])
        else:
            esc = {}

        rs = pprint_thing(s, escape_chars=esc).strip()

        if self.render_links and _is_url(rs):
            rs_unescaped = pprint_thing(s, escape_chars={}).strip()
            start_tag += '<a href="{url}" target="_blank">'.format(
                url=rs_unescaped)
            end_a = '</a>'
        else:
            end_a = ''

        self.write('{start}{rs}{end_a}</{kind}>'.format(
            start=start_tag, rs=rs, end_a=end_a, kind=kind), indent)

    def write_tr(self, line, indent=0, indent_delta=0, header=False,
                 align=None, tags=None, nindex_levels=0):
        if tags is None:
            tags = {}

        if align is None:
            self.write('<tr>', indent)
        else:
            self.write('<tr style="text-align: {align};">'
                       .format(align=align), indent)
        indent += indent_delta

        for i, s in enumerate(line):
            val_tag = tags.get(i, None)
            if header or (self.bold_rows and i < nindex_levels):
                self.write_th(s, indent, tags=val_tag)
            else:
                self.write_td(s, indent, tags=val_tag)

        indent -= indent_delta
        self.write('</tr>', indent)

    def render(self):
        self._write_table()

        if self.should_show_dimensions:
            by = chr(215)  # ×
            self.write('<p>{rows} rows {by} {cols} columns</p>'
                       .format(rows=len(self.frame),
                               by=by,
                               cols=len(self.frame.columns)))

        return self.elements

    def _write_table(self, indent=0):
        _classes = ['dataframe']  # Default class.
        use_mathjax = get_option("display.html.use_mathjax")
        if not use_mathjax:
            _classes.append('tex2jax_ignore')
        if self.classes is not None:
            if isinstance(self.classes, str):
                self.classes = self.classes.split()
            if not isinstance(self.classes, (list, tuple)):
                raise TypeError('classes must be a string, list, or tuple, '
                                'not {typ}'.format(typ=type(self.classes)))
            _classes.extend(self.classes)

        if self.table_id is None:
            id_section = ""
        else:
            id_section = ' id="{table_id}"'.format(table_id=self.table_id)

        self.write('<table border="{border}" class="{cls}"{id_section}>'
                   .format(border=self.border, cls=' '.join(_classes),
                           id_section=id_section), indent)

        if self.fmt.header or self.show_row_idx_names:
            self._write_header(indent + self.indent_delta)

        self._write_body(indent + self.indent_delta)

        self.write('</table>', indent)

    def _write_col_header(self, indent):
        truncate_h = self.fmt.truncate_h
        if isinstance(self.columns, ABCMultiIndex):
            template = 'colspan="{span:d}" halign="left"'

            if self.fmt.sparsify:
                # GH3547
                sentinel = object()
            else:
                sentinel = False
            levels = self.columns.format(sparsify=sentinel, adjoin=False,
                                         names=False)
            level_lengths = get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            for lnum, (records, values) in enumerate(zip(level_lengths,
                                                         levels)):
                if truncate_h:
                    # modify the header lines
                    ins_col = self.fmt.tr_col_num
                    if self.fmt.sparsify:
                        recs_new = {}
                        # Increment tags after ... col.
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            elif tag + span > ins_col:
                                recs_new[tag] = span + 1
                                if lnum == inner_lvl:
                                    values = (values[:ins_col] + ('...',) +
                                              values[ins_col:])
                                else:
                                    # sparse col headers do not receive a ...
                                    values = (values[:ins_col] +
                                              (values[ins_col - 1], ) +
                                              values[ins_col:])
                            else:
                                recs_new[tag] = span
                            # if ins_col lies between tags, all col headers
                            # get ...
                            if tag + span == ins_col:
                                recs_new[ins_col] = 1
                                values = (values[:ins_col] + ('...',) +
                                          values[ins_col:])
                        records = recs_new
                        inner_lvl = len(level_lengths) - 1
                        if lnum == inner_lvl:
                            records[ins_col] = 1
                    else:
                        recs_new = {}
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            else:
                                recs_new[tag] = span
                        recs_new[ins_col] = 1
                        records = recs_new
                        values = (values[:ins_col] + ['...'] +
                                  values[ins_col:])

                # see gh-22579
                # Column Offset Bug with to_html(index=False) with
                # MultiIndex Columns and Index.
                # Initially fill row with blank cells before column names.
                # TODO: Refactor to remove code duplication with code
                # block below for standard columns index.
                row = [''] * (self.row_levels - 1)
                if self.fmt.index or self.show_col_idx_names:
                    # see gh-22747
                    # If to_html(index_names=False) do not show columns
                    # index names.
                    # TODO: Refactor to use _get_column_name_list from
                    # DataFrameFormatter class and create a
                    # _get_formatted_column_labels function for code
                    # parity with DataFrameFormatter class.
                    if self.fmt.show_index_names:
                        name = self.columns.names[lnum]
                        row.append(pprint_thing(name or ''))
                    else:
                        row.append('')

                tags = {}
                j = len(row)
                for i, v in enumerate(values):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template.format(span=records[i])
                    else:
                        continue
                    j += 1
                    row.append(v)
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              header=True)
        else:
            # see gh-22579
            # Column misalignment also occurs for
            # a standard index when the columns index is named.
            # Initially fill row with blank cells before column names.
            # TODO: Refactor to remove code duplication with code block
            # above for columns MultiIndex.
            row = [''] * (self.row_levels - 1)
            if self.fmt.index or self.show_col_idx_names:
                # see gh-22747
                # If to_html(index_names=False) do not show columns
                # index names.
                # TODO: Refactor to use _get_column_name_list from
                # DataFrameFormatter class.
                if self.fmt.show_index_names:
                    row.append(self.columns.name or '')
                else:
                    row.append('')
            row.extend(self.columns)
            align = self.fmt.justify

            if truncate_h:
                ins_col = self.row_levels + self.fmt.tr_col_num
                row.insert(ins_col, '...')

            self.write_tr(row, indent, self.indent_delta, header=True,
                          align=align)

    def _write_row_header(self, indent):
        truncate_h = self.fmt.truncate_h
        row = ([x if x is not None else '' for x in self.frame.index.names]
               + [''] * (self.ncols + (1 if truncate_h else 0)))
        self.write_tr(row, indent, self.indent_delta, header=True)

    def _write_header(self, indent):
        self.write('<thead>', indent)

        if self.fmt.header:
            self._write_col_header(indent + self.indent_delta)

        if self.show_row_idx_names:
            self._write_row_header(indent + self.indent_delta)

        self.write('</thead>', indent)

    def _get_formatted_values(self):
        with option_context('display.max_colwidth', 999999):
            fmt_values = {i: self.fmt._format_col(i)
                          for i in range(self.ncols)}
        return fmt_values

    def _write_body(self, indent):
        self.write('<tbody>', indent)
        fmt_values = self._get_formatted_values()

        # write values
        if self.fmt.index and isinstance(self.frame.index, ABCMultiIndex):
            self._write_hierarchical_rows(
                fmt_values, indent + self.indent_delta)
        else:
            self._write_regular_rows(
                fmt_values, indent + self.indent_delta)

        self.write('</tbody>', indent)

    def _write_regular_rows(self, fmt_values, indent):
        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v

        nrows = len(self.fmt.tr_frame)

        if self.fmt.index:
            fmt = self.fmt._get_formatter('__index__')
            if fmt is not None:
                index_values = self.fmt.tr_frame.index.map(fmt)
            else:
                index_values = self.fmt.tr_frame.index.format()

        row = []
        for i in range(nrows):

            if truncate_v and i == (self.fmt.tr_row_num):
                str_sep_row = ['...'] * len(row)
                self.write_tr(str_sep_row, indent, self.indent_delta,
                              tags=None, nindex_levels=self.row_levels)

            row = []
            if self.fmt.index:
                row.append(index_values[i])
            # see gh-22579
            # Column misalignment also occurs for
            # a standard index when the columns index is named.
            # Add blank cell before data cells.
            elif self.show_col_idx_names:
                row.append('')
            row.extend(fmt_values[j][i] for j in range(self.ncols))

            if truncate_h:
                dot_col_ix = self.fmt.tr_col_num + self.row_levels
                row.insert(dot_col_ix, '...')
            self.write_tr(row, indent, self.indent_delta, tags=None,
                          nindex_levels=self.row_levels)

    def _write_hierarchical_rows(self, fmt_values, indent):
        template = 'rowspan="{span}" valign="top"'

        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v
        frame = self.fmt.tr_frame
        nrows = len(frame)

        idx_values = frame.index.format(sparsify=False, adjoin=False,
                                        names=False)
        idx_values = lzip(*idx_values)

        if self.fmt.sparsify:
            # GH3547
            sentinel = object()
            levels = frame.index.format(sparsify=sentinel, adjoin=False,
                                        names=False)

            level_lengths = get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            if truncate_v:
                # Insert ... row and adjust idx_values and
                # level_lengths to take this into account.
                ins_row = self.fmt.tr_row_num
                inserted = False
                for lnum, records in enumerate(level_lengths):
                    rec_new = {}
                    for tag, span in list(records.items()):
                        if tag >= ins_row:
                            rec_new[tag + 1] = span
                        elif tag + span > ins_row:
                            rec_new[tag] = span + 1

                            # GH 14882 - Make sure insertion done once
                            if not inserted:
                                dot_row = list(idx_values[ins_row - 1])
                                dot_row[-1] = '...'
                                idx_values.insert(ins_row, tuple(dot_row))
                                inserted = True
                            else:
                                dot_row = list(idx_values[ins_row])
                                dot_row[inner_lvl - lnum] = '...'
                                idx_values[ins_row] = tuple(dot_row)
                        else:
                            rec_new[tag] = span
                        # If ins_row lies between tags, all cols idx cols
                        # receive ...
                        if tag + span == ins_row:
                            rec_new[ins_row] = 1
                            if lnum == 0:
                                idx_values.insert(ins_row, tuple(
                                    ['...'] * len(level_lengths)))

                            # GH 14882 - Place ... in correct level
                            elif inserted:
                                dot_row = list(idx_values[ins_row])
                                dot_row[inner_lvl - lnum] = '...'
                                idx_values[ins_row] = tuple(dot_row)
                    level_lengths[lnum] = rec_new

                level_lengths[inner_lvl][ins_row] = 1
                for ix_col in range(len(fmt_values)):
                    fmt_values[ix_col].insert(ins_row, '...')
                nrows += 1

            for i in range(nrows):
                row = []
                tags = {}

                sparse_offset = 0
                j = 0
                for records, v in zip(level_lengths, idx_values[i]):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template.format(span=records[i])
                    else:
                        sparse_offset += 1
                        continue

                    j += 1
                    row.append(v)

                row.extend(fmt_values[j][i] for j in range(self.ncols))
                if truncate_h:
                    row.insert(self.row_levels - sparse_offset +
                               self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              nindex_levels=len(levels) - sparse_offset)
        else:
            row = []
            for i in range(len(frame)):
                if truncate_v and i == (self.fmt.tr_row_num):
                    str_sep_row = ['...'] * len(row)
                    self.write_tr(str_sep_row, indent, self.indent_delta,
                                  tags=None, nindex_levels=self.row_levels)

                idx_values = list(zip(*frame.index.format(
                    sparsify=False, adjoin=False, names=False)))
                row = []
                row.extend(idx_values[i])
                row.extend(fmt_values[j][i] for j in range(self.ncols))
                if truncate_h:
                    row.insert(self.row_levels + self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=None,
                              nindex_levels=frame.index.nlevels)


class NotebookFormatter(HTMLFormatter):
    """
    Internal class for formatting output data in html for display in Jupyter
    Notebooks. This class is intended for functionality specific to
    DataFrame._repr_html_() and DataFrame.to_html(notebook=True)
    """

    def _get_formatted_values(self):
        return {i: self.fmt._format_col(i) for i in range(self.ncols)}

    def write_style(self):
        # We use the "scoped" attribute here so that the desired
        # style properties for the data frame are not then applied
        # throughout the entire notebook.
        template_first = """\
            <style scoped>"""
        template_last = """\
            </style>"""
        template_select = """\
                .dataframe %s {
                    %s: %s;
                }"""
        element_props = [('tbody tr th:only-of-type',
                          'vertical-align',
                          'middle'),
                         ('tbody tr th',
                          'vertical-align',
                          'top')]
        if isinstance(self.columns, ABCMultiIndex):
            element_props.append(('thead tr th',
                                  'text-align',
                                  'left'))
            if self.show_row_idx_names:
                element_props.append(('thead tr:last-of-type th',
                                      'text-align',
                                      'right'))
        else:
            element_props.append(('thead th',
                                  'text-align',
                                  'right'))
        template_mid = '\n\n'.join(map(lambda t: template_select % t,
                                       element_props))
        template = dedent('\n'.join((template_first,
                                     template_mid,
                                     template_last)))
        self.write(template)

    def render(self):
        self.write('<div>')
        self.write_style()
        super(NotebookFormatter, self).render()
        self.write('</div>')
        return self.elements
