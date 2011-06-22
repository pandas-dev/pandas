"""
High level interface to PyTables for reading and writing pandas data structures
to disk
"""

# pylint: disable-msg=E1101,W0613

from datetime import datetime
import cPickle as pickle
import os.path
import time

import numpy as np
import tables

from pandas import (Series, TimeSeries, DataFrame, DataMatrix, WidePanel,
                    LongPanel)
from pandas.core.pytools import adjoin
import pandas.core.internals as internals

# reading and writing the full object in one go
_TYPE_MAP = {
    Series     : 'series',
    TimeSeries : 'time_series',
    DataFrame  : 'frame',
    WidePanel  : 'wide',
    LongPanel  : 'long'
}

_DEFAULT_TABLE = {
    'wide' : True
}

_NAME_MAP = {
    'series' : 'Series',
    'time_series' : 'TimeSeries',
    'frame' : 'DataFrame',
    'frame_table' : 'DataFrame',
    'wide' : 'WidePanel',
    'wide_table' : 'WidePanel',
    'long' : 'LongPanel'
}

from line_profiler import LineProfiler
prof = LineProfiler()

class HDFStore(object):
    """
    dict-like IO interface for storing pandas objects in PyTables
    format

    Parameters
    ----------
    path : string
        File path to HDF5 file
    mode : {'a', 'w', 'r', 'r+'}, default 'a'

        ``'r'``
            Read-only; no data can be modified.
        ``'w``'
            Write; a new file is created (an existing file with the same
            name would be deleted).
        ``'a'``
            Append; an existing file is opened for reading and writing,
            and if the file does not exist it is created.
        ``'r+'``
            It is similar to ``'a'``, but the file must already exist.

    Examples
    --------
    >>> store = HDFStore('test.h5')
    >>> store['foo'] = bar   # write to HDF5
    >>> bar = store['foo']   # retrieve
    >>> store.close()
    """
    def __init__(self, path, mode='a'):
        self.path = path
        self.mode = mode
        self.handle = None
        self.open(mode=mode, warn=False)

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.put(key, value)

    def __len__(self):
        return len(self.handle.root._v_children)

    def __repr__(self):
        output = '%s\nFile path: %s\n' % (type(self), self.path)

        if len(self) > 0:
            keys = []
            values = []
            for k, v in sorted(self.handle.root._v_children.iteritems()):
                kind = v._v_attrs.pandas_type

                keys.append(str(k))
                values.append(_NAME_MAP[kind])

            output += adjoin(5, keys, values)
        else:
            output += 'Empty'

        return output

    def open(self, mode='a', warn=True):
        """
        Open the file in the specified mode

        Parameters
        ----------
        mode : {'a', 'w', 'r', 'r+'}, default 'a'
            See HDFStore docstring or tables.openFile for info about modes
        """
        self.mode = mode
        if warn and mode == 'w': # pragma: no cover
            while True:
                response = raw_input("Re-opening as mode='w' will delete the "
                                     "current file. Continue (y/n)?")
                if response == 'y':
                    break
                elif response == 'n':
                    return
        if self.handle is not None and self.handle.isopen:
            self.handle.close()
        self.handle = tables.openFile(self.path, self.mode)

    def close(self):
        """
        Close the PyTables file handle
        """
        self.handle.close()

    def flush(self):
        """
        Force all buffered modifications to be written to disk
        """
        self.handle.flush()

    def get(self, key):
        """
        Retrieve pandas object stored in file

        Parameters
        ----------
        key : object

        Returns
        -------
        obj : type of object stored in file
        """
        group = getattr(self.handle.root, key, None)
        if group is not None:
            return self._read_group(group)

    def select(self, key, where=None):
        """
        Retrieve pandas object stored in file, optionally based on where
        criteria

        Parameters
        ----------
        key : object
        """
        group = getattr(self.handle.root, key, None)
        if group is not None:
            return self._read_group(group, where)

    def put(self, key, value, table=None, append=False,
            compression=None):
        """
        Store object in HDFStore

        Parameters
        ----------
        key : object
        value : {Series, DataFrame, WidePanel, LongPanel}
        table : boolean, default False
            Write as a PyTables Table structure which may perform worse but
            allow more flexible operations like searching / selecting subsets of
            the data
        append : boolean, default False
            For table data structures, append the input data to the existing
            table
        compression : {None, 'blosc', 'lzo', 'zlib'}, default None
            Use a compression algorithm to compress the data
        """
        self._write_to_group(key, value, table=table, append=append,
                             comp=compression)

    def _get_handler(self, op, kind):
        return getattr(self,'_%s_%s' % (op, kind))

    def remove(self, key, where):
        """
        Remove pandas object partially by specifying the where condition

        Parameters
        ----------
        key : object
        """
        group = getattr(self.handle.root, key, None)
        if group is not None:
            return self._delete_group(group, where)

    def append(self, key, value):
        """
        Store object in file

        Parameters
        ----------
        key : object
        <value : {Series, DataFrame, WidePanel, LongPanel}
        pandas data structure
        """
        self._write_group(key, value, True)

    def _write_to_group(self, key, value, table=None, append=False,
                        comp=None):
        # TODO WIDEPANEL DEFAULT TO TABLE

        root = self.handle.root
        if key not in root._v_children:
            group = self.handle.createGroup(root, key)
        else:
            group = getattr(root, key)

        kind = _TYPE_MAP[type(value)]
        if table is None:
            table = _DEFAULT_TABLE.get(kind, False)

        if table:
            kind = '%s_table' % kind
            handler = self._get_handler(op='write', kind=kind)
            wrapper = lambda value: handler(group, value, append=append)
        else:
            if append:
                raise ValueError('Can only append to Tables')
            handler = self._get_handler(op='write', kind=kind)
            wrapper = lambda value: handler(group, value)

        wrapper(value)
        group._v_attrs.pandas_type = kind

    def _write_series(self, group, series):
        self._write_index(group, 'index', series.index)
        self._write_array(group, 'values', series.values)

    def _write_frame(self, group, df):
        self._write_index(group, 'index', df.index)
        self._write_index(group, 'columns', df.columns)
        self._write_array(group, 'values', df.values)

    def _write_frame_table(self, group, df, append=False, comp=None):
        self._write_table(group, index=df.index, columns=df.columns,
                          values_dict=dict(value=df.values),
                          append=append, compression=comp)

    def _write_wide(self, group, panel, append=False, comp=None):
        value_d = dict((k, v.values) for k, v in panel.iteritems())
        self._write_table(group, index=panel.major_axis,
                          columns=panel.minor_axis, values_dict=value_d,
                          append=append, compression=comp)

    def _write_wide_table(self, group, panel, append=False, comp=None):
        value_d = dict((k, v.values) for k, v in panel.iteritems())
        self._write_table(group, index=panel.major_axis,
                          columns=panel.minor_axis, values_dict=value_d,
                          append=append, compression=comp)

    def _write_long(self, group, value, append=False):
        pass

    def _write_index(self, group, key, value):
        # don't care about type here
        converted, kind, _ = _convert_index(value)
        self._write_array(group, key, converted)
        node = getattr(group, key)
        node._v_attrs.kind = kind

    def _write_array(self, group, key, value):
        if key in group:
            self.handle.removeNode(group, key)

        self.handle.createArray(group, key, value)

    def _write_table(self, group, index, columns=None, values_dict=None,
                     append=False, compression=None):
        """ need to check for conform to the existing table:
            e.g. columns should match """
        # create dict of types
        index_converted, index_kind, index_t = _convert_index(index)
        columns_converted, cols_kind, col_t = _convert_index(columns)

        # create the table if it doesn't exist (or get it if it does)
        number_fields = len(values_dict)
        if not append:
            if 'table' in group:
                self.handle.removeNode(group, 'table')

        if 'table' not in group:
            # create the table
            desc = {'index'  : index_t,
                    'column' : col_t,
                    'values' : tables.FloatCol(shape=(number_fields))}

            options = {'name' : 'table',
                       'description' : desc}

            if compression:
                options['filters'] = tables.Filters(complevel=9,
                                                    complib=compression)

            table = self.handle.createTable(group, **options)
        else:
            # the table must already exist
            table = getattr(group, 'table', None)

        # add kinds
        fields = values_dict.keys()
        table._v_attrs.index_kind = index_kind
        table._v_attrs.columns_kind = cols_kind
        table._v_attrs.fields = fields

        # add the rows
        values = values_dict.values()

        try:
            for i, index in enumerate(index_converted):
                for c, col in enumerate(columns_converted):
                    v = np.asarray([v[i, c] for v in values])

                    # don't store the row if all values are np.nan
                    if np.isnan(v).all():
                        continue

                    row = table.row
                    row['index'] = index
                    row['column'] = col

                    # create the values array
                    row['values'] = v
                    row.append()
                    self.handle.flush()
        except (ValueError), detail:
            print "value_error in _write_table -> %s" % str(detail)
            raise

    def _read_group(self, group, where=None):
        kind = group._v_attrs.pandas_type
        handler = self._get_handler(op='read', kind=kind)
        try:
            return handler(group, where)
        except Exception:
            raise

    def _delete_group(self,group, where):
        return self._delete_from_table(group, where)

    def _read_series(self, group, where=None):
        index = self._read_index(group, 'index')
        values = group.values[:]
        return Series(values, index=index)

    def _read_frame(self, group, where=None):
        index = self._read_index(group, 'index')
        columns = self._read_index(group, 'columns')
        values = group.values[:]
        return DataFrame(values, index=index, columns=columns)

    def _read_frame_table(self, group, where=None):
        # items is a place holder since the data is stored panel-like
        return self._read_table3(group, where)

        # _, index, columns, values = self._read_table2(group, where)
        # return DataFrame(values.squeeze(), index=index, columns=columns)

    def _read_wide_table(self, group, where=None):
        return self._read_table3(group, where)

        # items, index, columns, values = self._read_table2(group, where)
        # if values is None: return None
        # return WidePanel(values, items, index, columns)

    def _read_long(self, group, where=None):
        index = self._read_index(group, 'index')
        values = group.values[:]

        return Series(values, index=index)

    def _read_index(self, group, key):
        node = getattr(group, key)
        data = node[:]
        kind = node._v_attrs.kind

        return _unconvert_index(data, kind)

    def _read_table2(self, group, where=None):
        table = getattr(group, 'table')

        # create the selection
        s = Selection(table, where)
        s.select()

        # reconstruct
        column_filter = s.column_filter
        col_cache = {}
        index_cache = {}
        fields = table._v_attrs.fields

        # dict - indexed by column
        d = {}

        for r in s.values:
            c, i, v = r

            try:
                c_converted = col_cache[c]
            except (KeyError):
                c_converted = col_cache[c] = _unconvert_value(c, table._v_attrs.columns_kind)

            if column_filter and c_converted not in column_filter:
                continue

            try:
                i_converted = index_cache[i]
            except (KeyError):
                i_converted = index_cache[i] = _unconvert_value(i, table._v_attrs.index_kind)

            d.setdefault(c_converted, {})[i_converted] = v

        # stack the values
        nan = np.asarray([np.nan] * len(fields))

        columns_set = set(col_cache.values())
        columns = sorted(list(column_filter & columns_set
                              if column_filter else columns_set))
        indicies = sorted(index_cache.values())

        # create the 3d stack (items x columns x indicies)
        if len(columns) and len(indicies):
            data = np.dstack([np.asarray([d[c].get(i,nan)
                                          for i in indicies ]).transpose()
                              for c in columns ])
        else:
            data = None

        return fields, indicies, columns, data

    def _read_table3(self, group, where=None):
        from pandas.core.panel import _make_long_index
        table = getattr(group, 'table')

        # create the selection
        sel = Selection(table, where)
        sel.select()
        fields = table._v_attrs.fields

        columns = _maybe_convert(sel.values['column'],
                                 table._v_attrs.columns_kind)
        index = _maybe_convert(sel.values['index'],
                               table._v_attrs.index_kind)
        # reconstruct
        long_index = _make_long_index(np.asarray(index),
                                      np.asarray(columns))
        lp = LongPanel(sel.values['values'], fields, long_index)
        wp = lp.to_wide()

        # TODO
        column_filter = sel.column_filter

        return wp

    def _delete_from_table(self, group, where = None):
        table = getattr(group, 'table')

        # create the selection
        s = Selection(table,where)
        s.select_coords()

        # delete the rows in reverse order
        l = list(s.values)
        l.reverse()
        for c in l:
            table.removeRows(c)
            self.handle.flush()
        return len(s.values)

    _write_time_series = _write_series
    _read_time_series = _read_series

def _convert_index(index):
    # Let's assume the index is homogeneous
    values = np.asarray(index)

    import time
    if isinstance(values[0], datetime):
        converted = np.array([time.mktime(v.timetuple())
                              for v in values], dtype=np.int64)
        return converted, 'datetime', tables.Time64Col()
    elif isinstance(values[0], basestring):
        converted = np.array(list(values), dtype=np.str_)
        itemsize = converted.dtype.itemsize
        return converted, 'string', tables.StringCol(itemsize)
    else:
        # take a guess for now, hope the values fit
        return np.array(list(values)), 'other', tables.StringCol(20)

def _unconvert_index(data, kind):
    if kind == 'datetime':
        index = np.array([datetime.fromtimestamp(v) for v in data],
                         dtype=object)
    elif kind == 'string':
        index = np.array(data, dtype=object)
    else:
        index = data

    return index

def _unconvert_value(item, kind):
    if kind == 'datetime':
        return datetime.fromtimestamp(item)
    return item

def _maybe_convert(values, val_kind):
    if _need_convert(val_kind):
        conv = _get_converter(val_kind)
        conv = np.frompyfunc(conv, 1, 1)
        values = conv(values)
    return values

def _get_converter(kind):
    if kind == 'datetime':
        return datetime.fromtimestamp
    else: # pragma: no cover
        raise ValueError('invalid kind %s' % kind)

def _need_convert(kind):
    if kind == 'datetime':
        return True
    return False

class Selection(object):
    """
    apply selection criteria
    where of the form:  [ list of and_ conditions ]
    condition of the form dict(field = 'index', 'op' = '<=',
                               value=datetime(2010,6,14))
    """
    def __init__(self, table, where=None):
        self.table = table
        self.where = where
        self.column_filter = None
        self.the_condition = None
        self.conditions = []

        if where:
            self.generate(where)

    def generate(self, where):
        # and condictions
        for c in where:
            op = c.get('op',None)
            value = c['value']
            field = c['field']

            if field == 'index' and isinstance(value, datetime):
                value = time.mktime(value.timetuple())
                self.conditions.append('(%s %s %s)' % (field,op,value))
            else:
                self.generate_multiple_conditions(op,value,field)

        if len(self.conditions):
            self.the_condition = '(' + ' & '.join(self.conditions) + ')'

    def generate_multiple_conditions(self, op, value, field):

        if op and op == 'in' or isinstance(value,list):
            if len(value) <= 61:
                l = '(' + ' | '.join([ "(%s == '%s')" % (field,v) for v in value ]) + ')'
                self.conditions.append(l)
            else:
                self.column_filter = set(value)
        else:
            if op is None:
                op = '=='
            self.conditions.append('(%s %s "%s")' % (field,op,value))

    def select(self):
        """ generate the selection """

        if self.the_condition:
            self.values = self.table.readWhere(self.the_condition)

        else:
            self.values = self.table.read()

    def select_coords(self):
        """ generate the selection """
        self.values = self.table.getWhereList(self.the_condition)

