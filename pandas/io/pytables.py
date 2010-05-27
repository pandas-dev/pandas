from datetime import datetime
import numpy as np

from pandas import (Series, TimeSeries, DataFrame, DataMatrix, WidePanel,
                    LongPanel)
from pandas.core.pytools import adjoin
import pandas.lib.tseries as tseries

try:
    import tables
except ImportError:
    pass

class HDFStore(object):
    """
    dict-like IO interface for storing pandas objects in PyTables
    format

    Parameters
    ----------
    path : string
        File path to HDF5 file

    Examples
    --------
    >>> store = HDFStore('test.h5')
    >>> store['foo'] = bar   # write to HDF5
    >>> bar = store['foo']   # retrieve
    >>> store.close()
    """
    def __init__(self, path):
        self.handle = tables.openFile(path, 'a')

    def close(self):
        self.handle.close()

    def flush(self):
        self.handle.flush()

    def __repr__(self):
        output = str(self.__class__) + '\n'

        keys = []
        values = []
        for k, v in sorted(self.handle.root._v_children.iteritems()):
            kind = v._v_attrs.pandas_type

            keys.append(str(k))
            values.append(kind)

        output += adjoin(5, keys, values)
        return output

    def get(self, key):
        """

        """

        return self[key]

    def put(self, key, value):
        """

        """
        self[key] = value

    def __getitem__(self, key):
        group = getattr(self.handle.root, key)
        return _read_group(group)

    def __setitem__(self, key, value):
        self._write_group(key, value)

    def _write_group(self, key, value):
        root = self.handle.root

        if key not in root._v_children:
            group = self.handle.createGroup(root, key)
        else:
            group = getattr(root, key)

        kind = type(value)
        handler = self._get_write_handler(kind)

        try:
            handler(group, value)
        except Exception:
            raise

        group._v_attrs.pandas_type = kind.__name__
        return True

    def _write_series(self, group, series):
        self._write_index(group, 'index', series.index)
        self._write_array(group, 'values', np.asarray(series))

    def _write_frame(self, group, df):
        self._write_index(group, 'index', df.index)
        self._write_index(group, 'columns', df.columns)
        self._write_array(group, 'values', df.asMatrix(df.columns))

    def _write_matrix(self, group, dm):
        self._write_index(group, 'index', dm.index)
        self._write_index(group, 'columns', dm.columns)
        self._write_array(group, 'values', dm.values)

        if dm.objects is not None:
            self._write_index(group, 'obj_columns', dm.objects.columns)
            self._write_array(group, 'obj_values', dm.objects.values)

    def _write_wide(self, group, value):
        pass

    def _write_long(self, group, value):
        pass

    def _write_index(self, group, key, value):
        converted, kind = _convert_index(value)
        self._write_array(group, key, converted)
        node = getattr(group, key)

        node._v_attrs.kind = kind

    def _write_array(self, group, key, value):
        if key in group:
            self.handle.removeNode(group, key)

        self.handle.createArray(group, key, value)

    def _get_write_handler(self, kind):
        handlers = {
            Series : self._write_series,
            TimeSeries : self._write_series,
            DataFrame : self._write_frame,
            DataMatrix : self._write_matrix,
            WidePanel : self._write_wide,
            LongPanel : self._write_long
        }
        return handlers[kind]

def _read_group(group):
    kind = group._v_attrs.pandas_type

    if kind in ('Series', 'TimeSeries'):
        return _read_series(group)
    elif kind == 'DataFrame':
        return _read_frame(group)
    elif kind == 'DataMatrix':
        return _read_matrix(group)
    elif kind == 'WidePanel':
        return _read_wide(group)
    elif kind == 'LongPanel':
        return _read_long(group)


def _read_series(group):
    index = _read_index(group, 'index')
    values = group.values[:]

    return Series(values, index=index)

def _read_frame(group):
    index = _read_index(group, 'index')
    columns = _read_index(group, 'columns')
    values = group.values[:]

    return DataFrame(values, index=index, columns=columns)

def _read_matrix(group):
    index = _read_index(group, 'index')
    columns = _read_index(group, 'columns')
    values = group.values[:]
    objects = None

    if hasattr(group, 'obj_columns'):
        obj_columns = _read_index(group, 'columns')
        obj_values = group.obj_values[:]
        objects = DataMatrix(obj_values, index=index, columns=obj_columns)

    return DataMatrix(values, index=index, columns=columns,
                      objects=objects)

def _read_wide(group):
    index = _read_index(group, 'index')
    values = group.values[:]

    return Series(values, index=index)

def _read_long(group):
    index = _read_index(group, 'index')
    values = group.values[:]

    return Series(values, index=index)

def _read_index(group, key):
    node = getattr(group, key)
    data = node[:]
    kind = node._v_attrs.kind

    return _unconvert_index(data, kind)

def _convert_index(index):
    # Let's assume the index is homogeneous
    values = np.asarray(index)

    if isinstance(values[0], datetime):
        converted = tseries.array_to_timestamp(values)
        return converted, 'datetime'
    elif isinstance(values[0], basestring):
        converted = np.array(list(values), dtype=np.str_)
        return converted, 'string'
    else:
        return values, 'other'

def _unconvert_index(data, kind):
    if kind == 'datetime':
        index = tseries.array_to_datetime(data)
    elif kind == 'string':
        index = np.array(data, dtype=object)
    else:
        index = data

    return index
