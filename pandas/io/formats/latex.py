# -*- coding: utf-8 -*-
"""
Module for formatting output data in Latex.
"""

from __future__ import print_function

from pandas.core.index import MultiIndex
from pandas import compat
from pandas.compat import range, map, zip, u
from pandas.io.formats.format import TableFormatter
import numpy as np


class LatexFormatter(TableFormatter):
    """ Used to render a DataFrame to a LaTeX tabular/longtable environment
    output.

    Parameters
    ----------
    formatter : `DataFrameFormatter`
    column_format : str, default None
        The columns format as specified in `LaTeX table format
        <https://en.wikibooks.org/wiki/LaTeX/Tables>`__ e.g 'rcl' for 3 columns
    longtable : boolean, default False
        Use a longtable environment instead of tabular.

    See Also
    --------
    HTMLFormatter
    """

    def __init__(self, formatter, column_format=None, longtable=False,
                 multicolumn=False, multicolumn_format=None, multirow=False):
        self.fmt = formatter
        self.frame = self.fmt.frame
        self.bold_rows = self.fmt.kwds.get('bold_rows', False)
        self.column_format = column_format
        self.longtable = longtable
        self.multicolumn = multicolumn
        self.multicolumn_format = multicolumn_format
        self.multirow = multirow

    def write_result(self, buf):
        """
        Render a DataFrame to a LaTeX tabular/longtable environment output.
        """

        # string representation of the columns
        if len(self.frame.columns) == 0 or len(self.frame.index) == 0:
            info_line = (u('Empty {name}\nColumns: {col}\nIndex: {idx}')
                         .format(name=type(self.frame).__name__,
                                 col=self.frame.columns,
                                 idx=self.frame.index))
            strcols = [[info_line]]
        else:
            strcols = self.fmt._to_str_columns()

        def get_col_type(dtype):
            if issubclass(dtype.type, np.number):
                return 'r'
            else:
                return 'l'

        # reestablish the MultiIndex that has been joined by _to_str_column
        if self.fmt.index and isinstance(self.frame.index, MultiIndex):
            clevels = self.frame.columns.nlevels
            strcols.pop(0)
            name = any(self.frame.index.names)
            cname = any(self.frame.columns.names)
            lastcol = self.frame.index.nlevels - 1
            previous_lev3 = None
            for i, lev in enumerate(self.frame.index.levels):
                lev2 = lev.format()
                blank = (' ' * len(lev2[0]) if lev2 else
                         ' ' * len(self.fmt.na_rep))
                # display column names in last index-column
                if cname and i == lastcol:
                    lev3 = [x if x else '{}' for x in self.frame.columns.names]
                else:
                    lev3 = [blank] * clevels
                if name:
                    lev3.append(lev.name)
                current_idx_val = None
                for level_idx in self.frame.index.labels[i]:
                    idx_val = (lev2[level_idx] if level_idx >= 0 else
                               self.fmt.na_rep)
                    if ((previous_lev3 is None or
                         previous_lev3[len(lev3)].isspace()) and
                            idx_val == current_idx_val):
                        # same index as above row and left index was the same
                        lev3.append(blank)
                    else:
                        lev3.append(idx_val)
                        current_idx_val = idx_val
                strcols.insert(i, lev3)
                previous_lev3 = lev3

        column_format = self.column_format
        if column_format is None:
            dtypes = self.frame.dtypes._values
            column_format = ''.join(map(get_col_type, dtypes))
            if self.fmt.index:
                index_format = 'l' * self.frame.index.nlevels
                column_format = index_format + column_format
        elif not isinstance(column_format,
                            compat.string_types):  # pragma: no cover
            raise AssertionError('column_format must be str or unicode, '
                                 'not {typ}'.format(typ=type(column_format)))

        if not self.longtable:
            buf.write('\\begin{{tabular}}{{{fmt}}}\n'
                      .format(fmt=column_format))
            buf.write('\\toprule\n')
        else:
            buf.write('\\begin{{longtable}}{{{fmt}}}\n'
                      .format(fmt=column_format))
            buf.write('\\toprule\n')

        ilevels = self.frame.index.nlevels
        clevels = self.frame.columns.nlevels
        nlevels = clevels
        if any(self.frame.index.names):
            nlevels += 1
        strrows = list(zip(*strcols))
        self.clinebuf = []

        for i, row in enumerate(strrows):
            if i == nlevels and self.fmt.header:
                buf.write('\\midrule\n')  # End of header
                if self.longtable:
                    buf.write('\\endhead\n')
                    buf.write('\\midrule\n')
                    buf.write('\\multicolumn{{{n}}}{{r}}{{{{Continued on next '
                              'page}}}} \\\\\n'.format(n=len(row)))
                    buf.write('\\midrule\n')
                    buf.write('\\endfoot\n\n')
                    buf.write('\\bottomrule\n')
                    buf.write('\\endlastfoot\n')
            if self.fmt.kwds.get('escape', True):
                # escape backslashes first
                crow = [(x.replace('\\', '\\textbackslash').replace('_', '\\_')
                         .replace('%', '\\%').replace('$', '\\$')
                         .replace('#', '\\#').replace('{', '\\{')
                         .replace('}', '\\}').replace('~', '\\textasciitilde')
                         .replace('^', '\\textasciicircum').replace('&', '\\&')
                         if (x and x != '{}') else '{}') for x in row]
            else:
                crow = [x if x else '{}' for x in row]
            if self.bold_rows and self.fmt.index:
                # bold row labels
                crow = ['\\textbf{{{x}}}'.format(x=x)
                        if j < ilevels and x.strip() not in ['', '{}'] else x
                        for j, x in enumerate(crow)]
            if i < clevels and self.fmt.header and self.multicolumn:
                # sum up columns to multicolumns
                crow = self._format_multicolumn(crow, ilevels)
            if (i >= nlevels and self.fmt.index and self.multirow and
                    ilevels > 1):
                # sum up rows to multirows
                crow = self._format_multirow(crow, ilevels, i, strrows)
            buf.write(' & '.join(crow))
            buf.write(' \\\\\n')
            if self.multirow and i < len(strrows) - 1:
                self._print_cline(buf, i, len(strcols))

        if not self.longtable:
            buf.write('\\bottomrule\n')
            buf.write('\\end{tabular}\n')
        else:
            buf.write('\\end{longtable}\n')

    def _format_multicolumn(self, row, ilevels):
        r"""
        Combine columns belonging to a group to a single multicolumn entry
        according to self.multicolumn_format

        e.g.:
        a &  &  & b & c &
        will become
        \multicolumn{3}{l}{a} & b & \multicolumn{2}{l}{c}
        """
        row2 = list(row[:ilevels])
        ncol = 1
        coltext = ''

        def append_col():
            # write multicolumn if needed
            if ncol > 1:
                row2.append('\\multicolumn{{{ncol:d}}}{{{fmt:s}}}{{{txt:s}}}'
                            .format(ncol=ncol, fmt=self.multicolumn_format,
                                    txt=coltext.strip()))
            # don't modify where not needed
            else:
                row2.append(coltext)
        for c in row[ilevels:]:
            # if next col has text, write the previous
            if c.strip():
                if coltext:
                    append_col()
                coltext = c
                ncol = 1
            # if not, add it to the previous multicolumn
            else:
                ncol += 1
        # write last column name
        if coltext:
            append_col()
        return row2

    def _format_multirow(self, row, ilevels, i, rows):
        r"""
        Check following rows, whether row should be a multirow

        e.g.:     becomes:
        a & 0 &   \multirow{2}{*}{a} & 0 &
          & 1 &     & 1 &
        b & 0 &   \cline{1-2}
                  b & 0 &
        """
        for j in range(ilevels):
            if row[j].strip():
                nrow = 1
                for r in rows[i + 1:]:
                    if not r[j].strip():
                        nrow += 1
                    else:
                        break
                if nrow > 1:
                    # overwrite non-multirow entry
                    row[j] = '\\multirow{{{nrow:d}}}{{*}}{{{row:s}}}'.format(
                        nrow=nrow, row=row[j].strip())
                    # save when to end the current block with \cline
                    self.clinebuf.append([i + nrow - 1, j + 1])
        return row

    def _print_cline(self, buf, i, icol):
        """
        Print clines after multirow-blocks are finished
        """
        for cl in self.clinebuf:
            if cl[0] == i:
                buf.write('\\cline{{{cl:d}-{icol:d}}}\n'
                          .format(cl=cl[1], icol=icol))
        # remove entries that have been written to buffer
        self.clinebuf = [x for x in self.clinebuf if x[0] != i]
