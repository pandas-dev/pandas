"""
High level interface to PyTables for reading and writing pandas data structures
to disk
"""

# pylint: disable-msg=E1101,W0613,W0603
from datetime import datetime, date
import time
import re
import copy
import itertools
import warnings

import numpy as np
from pandas import (Series, TimeSeries, DataFrame, Panel, Panel4D, Index,
                    MultiIndex, Int64Index, Timestamp)
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparsePanel
from pandas.sparse.array import BlockIndex, IntIndex
from pandas.tseries.api import PeriodIndex, DatetimeIndex
from pandas.core.base import StringMixin
from pandas.core.common import adjoin, is_list_like, pprint_thing
from pandas.core.algorithms import match, unique
from pandas.core.categorical import Categorical
from pandas.core.common import _asarray_tuplesafe
from pandas.core.internals import BlockManager, make_block
from pandas.core.reshape import block2d_to_blocknd, factor_indexer
from pandas.core.index import _ensure_index
import pandas.core.common as com
from pandas.tools.merge import concat
from pandas.util import py3compat
from pandas.io.common import PerformanceWarning

import pandas.lib as lib
import pandas.algos as algos
import pandas.tslib as tslib

from contextlib import contextmanager

# versioning attribute
_version = '0.10.1'

# PY3 encoding if we don't specify
_default_encoding = 'UTF-8'


def _ensure_decoded(s):
    """ if we have bytes, decode them to unicde """
    if isinstance(s, np.bytes_):
        s = s.decode('UTF-8')
    return s


def _ensure_encoding(encoding):
    # set the encoding if we need
    if encoding is None:
        if py3compat.PY3:
            encoding = _default_encoding
    return encoding


class IncompatibilityWarning(Warning):
    pass


incompatibility_doc = """
where criteria is being ignored as this version [%s] is too old (or
not-defined), read the file in and write it out to a new file to upgrade (with
the copy_to method)
"""


class AttributeConflictWarning(Warning):
    pass


attribute_conflict_doc = """
the [%s] attribute of the existing index is [%s] which conflicts with the new
[%s], resetting the attribute to None
"""


performance_doc = """
your performance may suffer as PyTables will pickle object types that it cannot
map directly to c-types [inferred_type->%s,key->%s] [items->%s]
"""

# map object types
_TYPE_MAP = {

    Series          : u'series',
    SparseSeries    : u'sparse_series',
    TimeSeries      : u'series',
    DataFrame       : u'frame',
    SparseDataFrame : u'sparse_frame',
    Panel           : u'wide',
    Panel4D         : u'ndim',
    SparsePanel     : u'sparse_panel'
}

# storer class map
_STORER_MAP = {
    u'TimeSeries'    : 'LegacySeriesStorer',
    u'Series'        : 'LegacySeriesStorer',
    u'DataFrame'     : 'LegacyFrameStorer',
    u'DataMatrix'    : 'LegacyFrameStorer',
    u'series'        : 'SeriesStorer',
    u'sparse_series' : 'SparseSeriesStorer',
    u'frame'         : 'FrameStorer',
    u'sparse_frame'  : 'SparseFrameStorer',
    u'wide'          : 'PanelStorer',
    u'sparse_panel'  : 'SparsePanelStorer',
}

# table class map
_TABLE_MAP = {
    u'generic_table'    : 'GenericTable',
    u'appendable_frame'      : 'AppendableFrameTable',
    u'appendable_multiframe' : 'AppendableMultiFrameTable',
    u'appendable_panel' : 'AppendablePanelTable',
    u'appendable_ndim'  : 'AppendableNDimTable',
    u'worm'             : 'WORMTable',
    u'legacy_frame'     : 'LegacyFrameTable',
    u'legacy_panel'     : 'LegacyPanelTable',
}

# axes map
_AXES_MAP = {
    DataFrame: [0],
    Panel: [1, 2],
    Panel4D: [1, 2, 3],
}

# oh the troubles to reduce import time
_table_mod = None
_table_supports_index = False


def _tables():
    global _table_mod
    global _table_supports_index
    if _table_mod is None:
        import tables
        from distutils.version import LooseVersion
        _table_mod = tables

        # version requirements
        ver = tables.__version__
        _table_supports_index = LooseVersion(ver) >= '2.3'

    return _table_mod

def h5_open(path, mode):
    tables = _tables()
    return tables.openFile(path, mode)


@contextmanager
def get_store(path, mode='a', complevel=None, complib=None,
              fletcher32=False):
    """
    Creates an HDFStore instance. This function can be used in a with statement

    Parameters
    ----------
    same as HDFStore

    Examples
    --------
    >>> with get_store('test.h5') as store:
    >>>     store['foo'] = bar   # write to HDF5
    >>>     bar = store['foo']   # retrieve
    """
    store = None
    try:
        store = HDFStore(path, mode=mode, complevel=complevel,
                         complib=complib, fletcher32=False)
        yield store
    finally:
        if store is not None:
            store.close()


### interface to/from ###

def to_hdf(path_or_buf, key, value, mode=None, complevel=None, complib=None, append=None, **kwargs):
    """ store this object, close it if we opened it """
    if append:
        f = lambda store: store.append(key, value, **kwargs)
    else:
        f = lambda store: store.put(key, value, **kwargs)

    if isinstance(path_or_buf, basestring):
        with get_store(path_or_buf, mode=mode, complevel=complevel, complib=complib) as store:
            f(store)
    else:
        f(path_or_buf)

def read_hdf(path_or_buf, key, **kwargs):
    """ read from the store, closeit if we opened it """
    f = lambda store, auto_close: store.select(key, auto_close=auto_close, **kwargs)

    if isinstance(path_or_buf, basestring):

        # can't auto open/close if we are using an iterator
        # so delegate to the iterator
        store = HDFStore(path_or_buf)
        try:
            return f(store, True)
        except:

            # if there is an error, close the store
            try:
                store.close()
            except:
                pass

            raise

    # a passed store; user controls open/close
    f(path_or_buf, False)

class HDFStore(StringMixin):
    """
    dict-like IO interface for storing pandas objects in PyTables
    format.

    DataFrame and Panel can be stored in Table format, which is slower to
    read and write but can be searched and manipulated more like an SQL
    table. See HDFStore.put for more information

    Parameters
    ----------
    path : string
        File path to HDF5 file
    mode : {'a', 'w', 'r', 'r+'}, default 'a'

        ``'r'``
            Read-only; no data can be modified.
        ``'w'``
            Write; a new file is created (an existing file with the same
            name would be deleted).
        ``'a'``
            Append; an existing file is opened for reading and writing,
            and if the file does not exist it is created.
        ``'r+'``
            It is similar to ``'a'``, but the file must already exist.
    complevel : int, 1-9, default 0
            If a complib is specified compression will be applied
            where possible
    complib : {'zlib', 'bzip2', 'lzo', 'blosc', None}, default None
            If complevel is > 0 apply compression to objects written
            in the store wherever possible
    fletcher32 : bool, default False
            If applying compression use the fletcher32 checksum

    Examples
    --------
    >>> store = HDFStore('test.h5')
    >>> store['foo'] = bar   # write to HDF5
    >>> bar = store['foo']   # retrieve
    >>> store.close()
    """
    _quiet = False

    def __init__(self, path, mode=None, complevel=None, complib=None,
                 fletcher32=False):
        try:
            import tables as _
        except ImportError:  # pragma: no cover
            raise Exception('HDFStore requires PyTables')

        self._path = path
        if mode is None:
            mode = 'a'
        self._mode = mode
        self._handle = None
        self._complevel = complevel
        self._complib = complib
        self._fletcher32 = fletcher32
        self._filters = None
        self.open(mode=mode, warn=False)

    @property
    def root(self):
        """ return the root node """
        return self._handle.root

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.put(key, value)

    def __delitem__(self, key):
        return self.remove(key)

    def __getattr__(self, name):
        """ allow attribute access to get stores """
        try:
            return self.get(name)
        except:
            pass
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __contains__(self, key):
        """ check for existance of this key
              can match the exact pathname or the pathnm w/o the leading '/'
              """
        node = self.get_node(key)
        if node is not None:
            name = node._v_pathname
            if name == key or name[1:] == key: return True
        return False

    def __len__(self):
        return len(self.groups())

    def __unicode__(self):
        output = '%s\nFile path: %s\n' % (type(self), pprint_thing(self._path))

        if len(self.keys()):
            keys   = []
            values = []

            for k in self.keys():
                try:
                    s = self.get_storer(k)
                    if s is not None:
                        keys.append(pprint_thing(s.pathname or k))
                        values.append(pprint_thing(s or 'invalid_HDFStore node'))
                except Exception as detail:
                    keys.append(k)
                    values.append("[invalid_HDFStore node: %s]" % pprint_thing(detail))

            output += adjoin(12, keys, values)
        else:
            output += 'Empty'

        return output

    def keys(self):
        """
        Return a (potentially unordered) list of the keys corresponding to the
        objects stored in the HDFStore. These are ABSOLUTE path-names (e.g. have the leading '/'
        """
        return [n._v_pathname for n in self.groups()]

    def items(self):
        """
        iterate on key->group
        """
        for g in self.groups():
            yield g._v_pathname, g

    iteritems = items

    def open(self, mode='a', warn=True):
        """
        Open the file in the specified mode

        Parameters
        ----------
        mode : {'a', 'w', 'r', 'r+'}, default 'a'
            See HDFStore docstring or tables.openFile for info about modes
        """
        self._mode = mode
        if warn and mode == 'w':  # pragma: no cover
            while True:
                response = raw_input("Re-opening as mode='w' will delete the "
                                     "current file. Continue (y/n)?")
                if response == 'y':
                    break
                elif response == 'n':
                    return
        if self._handle is not None and self._handle.isopen:
            self._handle.close()

        if self._complib is not None:
            if self._complevel is None:
                self._complevel = 9
            self._filters = _tables().Filters(self._complevel,
                                              self._complib,
                                              fletcher32=self._fletcher32)

        try:
            self._handle = h5_open(self._path, self._mode)
        except IOError, e:  # pragma: no cover
            if 'can not be written' in str(e):
                print ('Opening %s in read-only mode' % self._path)
                self._handle = h5_open(self._path, 'r')
            else:
                raise

    def close(self):
        """
        Close the PyTables file handle
        """
        self._handle.close()

    def flush(self):
        """
        Force all buffered modifications to be written to disk
        """
        self._handle.flush()

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
        group = self.get_node(key)
        if group is None:
            raise KeyError('No object named %s in the file' % key)
        return self._read_group(group)

    def select(self, key, where=None, start=None, stop=None, columns=None, iterator=False, chunksize=None, auto_close=False, **kwargs):
        """
        Retrieve pandas object stored in file, optionally based on where
        criteria

        Parameters
        ----------
        key : object
        where : list of Term (or convertable) objects, optional
        start : integer (defaults to None), row number to start selection
        stop  : integer (defaults to None), row number to stop selection
        columns : a list of columns that if not None, will limit the return columns
        iterator : boolean, return an iterator, default False
        chunksize : nrows to include in iteration, return an iterator
        auto_close : boolean, should automatically close the store when finished, default is False

        """
        group = self.get_node(key)
        if group is None:
            raise KeyError('No object named %s in the file' % key)

        # create the storer and axes
        s = self._create_storer(group)
        s.infer_axes()

        # what we are actually going to do for a chunk
        def func(_start, _stop):
            return s.read(where=where, start=_start, stop=_stop, columns=columns, **kwargs)

        if iterator or chunksize is not None:
            if not s.is_table:
                raise TypeError("can only use an iterator or chunksize on a table")
            return TableIterator(self, func, nrows=s.nrows, start=start, stop=stop, chunksize=chunksize, auto_close=auto_close)

        return TableIterator(self, func, nrows=s.nrows, start=start, stop=stop, auto_close=auto_close).get_values()

    def select_as_coordinates(self, key, where=None, start=None, stop=None, **kwargs):
        """
        return the selection as a Coordinates.

        Parameters
        ----------
        key : object
        where : list of Term (or convertable) objects, optional
        start : integer (defaults to None), row number to start selection
        stop  : integer (defaults to None), row number to stop selection
        """
        return self.get_storer(key).read_coordinates(where=where, start=start, stop=stop, **kwargs)

    def unique(self, key, column, **kwargs):
        warnings.warn("unique(key,column) is deprecated\n"
                      "use select_column(key,column).unique() instead")
        return self.get_storer(key).read_column(column = column, **kwargs).unique()

    def select_column(self, key, column, **kwargs):
        """
        return a single column from the table. This is generally only useful to select an indexable

        Parameters
        ----------
        key : object
        column: the column of interest

        Exceptions
        ----------
        raises KeyError if the column is not found (or key is not a valid store)
        raises ValueError if the column can not be extracted indivually (it is part of a data block)

        """
        return self.get_storer(key).read_column(column = column, **kwargs)

    def select_as_multiple(self, keys, where=None, selector=None, columns=None, start=None, stop=None, iterator=False, chunksize=None, auto_close=False, **kwargs):
        """ Retrieve pandas objects from multiple tables

        Parameters
        ----------
        keys : a list of the tables
        selector : the table to apply the where criteria (defaults to keys[0] if not supplied)
        columns : the columns I want back
        start : integer (defaults to None), row number to start selection
        stop  : integer (defaults to None), row number to stop selection
        iterator : boolean, return an iterator, default False
        chunksize : nrows to include in iteration, return an iterator

        Exceptions
        ----------
        raise if any of the keys don't refer to tables or if they are not ALL THE SAME DIMENSIONS
        """

        # default to single select
        if isinstance(keys, (list, tuple)) and len(keys) == 1:
            keys = keys[0]
        if isinstance(keys, basestring):
            return self.select(key=keys, where=where, columns=columns, start=start, stop=stop, iterator=iterator, chunksize=chunksize, **kwargs)

        if not isinstance(keys, (list, tuple)):
            raise Exception("keys must be a list/tuple")

        if len(keys) == 0:
            raise Exception("keys must have a non-zero length")

        if selector is None:
            selector = keys[0]

        # collect the tables
        tbls = [ self.get_storer(k) for k in keys ]

        # validate rows
        nrows = None
        for t, k in zip(tbls, keys):
            if t is None:
                raise TypeError("Invalid table [%s]" % k)
            if not t.is_table:
                raise TypeError("object [%s] is not a table, and cannot be used in all select as multiple" % t.pathname)

            if nrows is None:
                nrows = t.nrows
            elif t.nrows != nrows:
                raise ValueError("all tables must have exactly the same nrows!")

        # select coordinates from the selector table
        try:
            c = self.select_as_coordinates(selector, where, start=start, stop=stop)
            nrows = len(c)
        except (Exception), detail:
            raise ValueError("invalid selector [%s]" % selector)

        def func(_start, _stop):

            # collect the returns objs
            objs = [t.read(where=c[_start:_stop], columns=columns) for t in tbls]

            # axis is the concentation axes
            axis = list(set([t.non_index_axes[0][0] for t in tbls]))[0]

            # concat and return
            return concat(objs, axis=axis, verify_integrity=True)

        if iterator or chunksize is not None:
            return TableIterator(self, func, nrows=nrows, start=start, stop=stop, chunksize=chunksize, auto_close=auto_close)

        return TableIterator(self, func, nrows=nrows, start=start, stop=stop, auto_close=auto_close).get_values()


    def put(self, key, value, table=None, append=False, **kwargs):
        """
        Store object in HDFStore

        Parameters
        ----------
        key      : object
        value    : {Series, DataFrame, Panel}
        table    : boolean, default False
            Write as a PyTables Table structure which may perform worse but
            allow more flexible operations like searching / selecting subsets
            of the data
        append   : boolean, default False
            For table data structures, append the input data to the existing
            table
        encoding : default None, provide an encoding for strings
        """
        self._write_to_group(key, value, table=table, append=append, **kwargs)

    def remove(self, key, where=None, start=None, stop=None):
        """
        Remove pandas object partially by specifying the where condition

        Parameters
        ----------
        key : string
            Node to remove or delete rows from
        where : list of Term (or convertable) objects, optional
        start : integer (defaults to None), row number to start selection
        stop  : integer (defaults to None), row number to stop selection

        Returns
        -------
        number of rows removed (or None if not a Table)

        Exceptions
        ----------
        raises KeyError if key is not a valid store

        """
        try:
            s = self.get_storer(key)
        except:

            if where is not None:
                raise ValueError("trying to remove a node with a non-None where clause!")

            # we are actually trying to remove a node (with children)
            s = self.get_node(key)
            if s is not None:
                s._f_remove(recursive=True)
                return None

        if s is None:
            raise KeyError('No object named %s in the file' % key)

        # remove the node
        if where is None:
            s.group._f_remove(recursive=True)

        # delete from the table
        else:
            if not s.is_table:
                raise ValueError('can only remove with where on objects written as tables')
            return s.delete(where = where, start=start, stop=stop)

    def append(self, key, value, columns=None, **kwargs):
        """
        Append to Table in file. Node must already exist and be Table
        format.

        Parameters
        ----------
        key : object
        value : {Series, DataFrame, Panel, Panel4D}
        data_columns : list of columns to create as data columns, or True to use all columns
        min_itemsize : dict of columns that specify minimum string sizes
        nan_rep      : string to use as string nan represenation
        chunksize    : size to chunk the writing
        expectedrows : expected TOTAL row size of this table
        encoding     : default None, provide an encoding for strings

        Notes
        -----
        Does *not* check if data being appended overlaps with existing
        data in the table, so be careful
        """
        if columns is not None:
            raise Exception("columns is not a supported keyword in append, try data_columns")

        self._write_to_group(key, value, table=True, append=True, **kwargs)

    def append_to_multiple(self, d, value, selector, data_columns=None, axes=None, **kwargs):
        """
        Append to multiple tables

        Parameters
        ----------
        d : a dict of table_name to table_columns, None is acceptable as the values of
            one node (this will get all the remaining columns)
        value : a pandas object
        selector : a string that designates the indexable table; all of its columns will
                   be designed as data_columns, unless data_columns is passed, in which
                   case these are used

        Notes
        -----
        axes parameter is currently not accepted

        """
        if axes is not None:
            raise Exception("axes is currently not accepted as a paremter to append_to_multiple; you can create the tables indepdently instead")

        if not isinstance(d, dict):
            raise ValueError("append_to_multiple must have a dictionary specified as the way to split the value")

        if selector not in d:
            raise ValueError("append_to_multiple requires a selector that is in passed dict")

        # figure out the splitting axis (the non_index_axis)
        axis = list(set(range(value.ndim)) - set(_AXES_MAP[type(value)]))[0]

        # figure out how to split the value
        remain_key = None
        remain_values = []
        for k, v in d.items():
            if v is None:
                if remain_key is not None:
                    raise ValueError("append_to_multiple can only have one value in d that is None")
                remain_key = k
            else:
                remain_values.extend(v)
        if remain_key is not None:
            ordered = value.axes[axis]
            ordd = ordered - Index(remain_values)
            ordd = sorted(ordered.get_indexer(ordd))
            d[remain_key] = ordered.take(ordd)

        # data_columns
        if data_columns is None:
            data_columns = d[selector]

        # append
        for k, v in d.items():
            dc = data_columns if k == selector else None

            # compute the val
            val = value.reindex_axis(v, axis=axis, copy=False)

            self.append(k, val, data_columns=dc, **kwargs)

    def create_table_index(self, key, **kwargs):
        """ Create a pytables index on the table
        Paramaters
        ----------
        key : object (the node to index)

        Exceptions
        ----------
        raises if the node is not a table

        """

        # version requirements
        _tables()
        if not _table_supports_index:
            raise Exception("PyTables >= 2.3 is required for table indexing")

        s = self.get_storer(key)
        if s is None: return

        if not s.is_table:
            raise TypeError("cannot create table index on a non-table")
        s.create_index(**kwargs)

    def groups(self):
        """ return a list of all the top-level nodes (that are not themselves a pandas storage object) """
        _tables()
        return [ g for g in self._handle.walkNodes() if getattr(g._v_attrs,'pandas_type',None) or getattr(
            g,'table',None) or (isinstance(g,_table_mod.table.Table) and g._v_name != u'table') ]

    def get_node(self, key):
        """ return the node with the key or None if it does not exist """
        try:
            if not key.startswith('/'):
                key = '/' + key
            return self._handle.getNode(self.root, key)
        except:
            return None

    def get_storer(self, key):
        """ return the storer object for a key, raise if not in the file """
        group = self.get_node(key)
        if group is None:
            return None
        s = self._create_storer(group)
        s.infer_axes()
        return s

    def copy(self, file, mode = 'w', propindexes = True, keys = None, complib = None, complevel = None,
             fletcher32 = False, overwrite = True):
        """ copy the existing store to a new file, upgrading in place

            Parameters
            ----------
            propindexes: restore indexes in copied file (defaults to True)
            keys       : list of keys to include in the copy (defaults to all)
            overwrite  : overwrite (remove and replace) existing nodes in the new store (default is True)
            mode, complib, complevel, fletcher32 same as in HDFStore.__init__

            Returns
            -------
            open file handle of the new store

        """
        new_store = HDFStore(file, mode = mode, complib = complib, complevel = complevel, fletcher32 = fletcher32)
        if keys is None:
            keys = self.keys()
        if not isinstance(keys, (tuple,list)):
            keys = [ keys ]
        for k in keys:
            s    = self.get_storer(k)
            if s is not None:

                if k in new_store:
                    if overwrite:
                        new_store.remove(k)

                data = self.select(k)
                if s.is_table:

                    index = False
                    if propindexes:
                        index = [ a.name for a in s.axes if a.is_indexed ]
                    new_store.append(k, data, index=index, data_columns=getattr(s,'data_columns',None), encoding=s.encoding)
                else:
                    new_store.put(k, data, encoding=s.encoding)

        return new_store

    ###### private methods ######

    def _create_storer(self, group, value = None, table = False, append = False, **kwargs):
        """ return a suitable Storer class to operate """

        def error(t):
            raise TypeError("cannot properly create the storer for: [%s] [group->%s,value->%s,table->%s,append->%s,kwargs->%s]" %
                            (t,group,type(value),table,append,kwargs))

        pt = _ensure_decoded(getattr(group._v_attrs,'pandas_type',None))
        tt = _ensure_decoded(getattr(group._v_attrs,'table_type',None))

        # infer the pt from the passed value
        if pt is None:
            if value is None:

                _tables()
                if getattr(group,'table',None) or isinstance(group,_table_mod.table.Table):
                    pt = u'frame_table'
                    tt = u'generic_table'
                else:
                    raise TypeError("cannot create a storer if the object is not existing nor a value are passed")
            else:

                try:
                    pt = _TYPE_MAP[type(value)]
                except:
                    error('_TYPE_MAP')

                # we are actually a table
                if table or append:
                    pt += u'_table'

        # a storer node
        if u'table' not in pt:
            try:
                return globals()[_STORER_MAP[pt]](self, group, **kwargs)
            except:
                error('_STORER_MAP')

        # existing node (and must be a table)
        if tt is None:

            # if we are a writer, determin the tt
            if value is not None:

                if pt == u'frame_table':
                    index = getattr(value,'index',None)
                    if index is not None:
                        if index.nlevels == 1:
                            tt = u'appendable_frame'
                        elif index.nlevels > 1:
                            tt = u'appendable_multiframe'
                elif pt == u'wide_table':
                    tt  = u'appendable_panel'
                elif pt == u'ndim_table':
                    tt = u'appendable_ndim'

            else:

                # distiguish between a frame/table
                tt = u'legacy_panel'
                try:
                    fields = group.table._v_attrs.fields
                    if len(fields) == 1 and fields[0] == u'value':
                        tt = u'legacy_frame'
                except:
                    pass

        try:
            return globals()[_TABLE_MAP[tt]](self, group, **kwargs)
        except:
            error('_TABLE_MAP')

    def _write_to_group(self, key, value, index=True, table=False, append=False,
                        complib=None, encoding=None, **kwargs):
        group = self.get_node(key)

        # remove the node if we are not appending
        if group is not None and not append:
            self._handle.removeNode(group, recursive=True)
            group = None

        if group is None:
            paths = key.split('/')

            # recursively create the groups
            path = '/'
            for p in paths:
                if not len(p):
                    continue
                new_path = path
                if not path.endswith('/'):
                    new_path += '/'
                new_path += p
                group = self.get_node(new_path)
                if group is None:
                    group = self._handle.createGroup(path, p)
                path = new_path

        s = self._create_storer(group, value, table=table, append=append,
                                encoding=encoding, **kwargs)
        if append:
            # raise if we are trying to append to a non-table,
            #       or a table that exists (and we are putting)
            if not s.is_table or (s.is_table and table is None and s.is_exists):
                raise ValueError('Can only append to Tables')
            if not s.is_exists:
                s.set_object_info()
        else:
            s.set_object_info()

        if not s.is_table and complib:
            raise ValueError('Compression not supported on non-table')

        s.write(obj = value, append=append, complib=complib, **kwargs)
        if s.is_table and index:
            s.create_index(columns = index)

    def _read_group(self, group, **kwargs):
        s = self._create_storer(group)
        s.infer_axes()
        return s.read(**kwargs)

class TableIterator(object):
    """ define the iteration interface on a table

        Parameters
        ----------

        store : the reference store
        func  : the function to get results
        nrows : the rows to iterate on
        start : the passed start value (default is None)
        stop  : the passed stop value (default is None)
        chunksize : the passed chunking valeu (default is 50000)
        auto_close : boolean, automatically close the store at the end of iteration,
            default is False
        kwargs : the passed kwargs
        """

    def __init__(self, store, func, nrows, start=None, stop=None, chunksize=None, auto_close=False):
        self.store = store
        self.func  = func
        self.nrows = nrows or 0
        self.start = start or 0

        if stop is None:
            stop = self.nrows
        self.stop  = min(self.nrows,stop)

        if chunksize is None:
            chunksize = 100000

        self.chunksize = chunksize
        self.auto_close = auto_close

    def __iter__(self):
        current = self.start
        while current < self.stop:
            stop = current + self.chunksize
            v = self.func(current, stop)
            current = stop

            if v is None:
                continue

            yield v

        self.close()

    def close(self):
        if self.auto_close:
            self.store.close()

    def get_values(self):
        results = self.func(self.start, self.stop)
        self.close()
        return results

class IndexCol(StringMixin):
    """ an index column description class

        Parameters
        ----------

        axis   : axis which I reference
        values : the ndarray like converted values
        kind   : a string description of this type
        typ    : the pytables type
        pos    : the position in the pytables

        """
    is_an_indexable = True
    is_data_indexable = True
    is_searchable = False
    _info_fields = ['freq','tz','index_name']

    def __init__(self, values=None, kind=None, typ=None, cname=None, itemsize=None,
                 name=None, axis=None, kind_attr=None, pos=None, freq=None, tz=None,
                 index_name=None, **kwargs):
        self.values = values
        self.kind = kind
        self.typ = typ
        self.itemsize = itemsize
        self.name = name
        self.cname = cname
        self.kind_attr = kind_attr
        self.axis = axis
        self.pos = pos
        self.freq = freq
        self.tz = tz
        self.index_name = index_name
        self.table = None

        if name is not None:
            self.set_name(name, kind_attr)
        if pos is not None:
            self.set_pos(pos)

    def set_name(self, name, kind_attr=None):
        """ set the name of this indexer """
        self.name = name
        self.kind_attr = kind_attr or "%s_kind" % name
        if self.cname is None:
            self.cname = name

        return self

    def set_axis(self, axis):
        """ set the axis over which I index """
        self.axis = axis

        return self

    def set_pos(self, pos):
        """ set the position of this column in the Table """
        self.pos = pos
        if pos is not None and self.typ is not None:
            self.typ._v_pos = pos
        return self

    def set_table(self, table):
        self.table = table
        return self

    def __unicode__(self):
        temp = tuple(map(pprint_thing, (self.name, self.cname, self.axis, self.pos, self.kind)))
        return "name->%s,cname->%s,axis->%s,pos->%s,kind->%s" % temp

    def __eq__(self, other):
        """ compare 2 col items """
        return all([getattr(self, a, None) == getattr(other, a, None) for a in ['name', 'cname', 'axis', 'pos']])

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def is_indexed(self):
        """ return whether I am an indexed column """
        try:
            return getattr(self.table.cols,self.cname).is_indexed
        except:
            False

    def copy(self):
        new_self = copy.copy(self)
        return new_self

    def infer(self, table):
        """ infer this column from the table: create and return a new object """
        new_self = self.copy()
        new_self.set_table(table)
        new_self.get_attr()
        return new_self

    def convert(self, values, nan_rep, encoding):
        """ set the values from this selection: take = take ownership """
        try:
            values = values[self.cname]
        except:
            pass

        values =_maybe_convert(values, self.kind, encoding)

        kwargs = dict()
        if self.freq is not None:
            kwargs['freq'] = _ensure_decoded(self.freq)
        if self.index_name is not None:
            kwargs['name'] = _ensure_decoded(self.index_name)
        try:
            self.values = Index(values, **kwargs)
        except:

            # if the output freq is different that what we recorded, then infer it
            if 'freq' in kwargs:
                kwargs['freq'] = 'infer'
            self.values = Index(_maybe_convert(values, self.kind, encoding), **kwargs)

        # set the timezone if indicated
        # we stored in utc, so reverse to local timezone
        if self.tz is not None:
            self.values = self.values.tz_localize('UTC').tz_convert(_ensure_decoded(self.tz))

        return self

    def take_data(self):
        """ return the values & release the memory """
        self.values, values = None, self.values
        return values

    @property
    def attrs(self):
        return self.table._v_attrs

    @property
    def description(self):
        return self.table.description

    @property
    def col(self):
        """ return my current col description """
        return getattr(self.description, self.cname, None)

    @property
    def cvalues(self):
        """ return my cython values """
        return self.values

    def __iter__(self):
        return iter(self.values)

    def maybe_set_size(self, min_itemsize=None, **kwargs):
        """ maybe set a string col itemsize:
               min_itemsize can be an interger or a dict with this columns name with an integer size """
        if _ensure_decoded(self.kind) == u'string':

            if isinstance(min_itemsize, dict):
                min_itemsize = min_itemsize.get(self.name)

            if min_itemsize is not None and self.typ.itemsize < min_itemsize:
                self.typ = _tables(
                ).StringCol(itemsize=min_itemsize, pos=self.pos)

    def validate_and_set(self, table, append, **kwargs):
        self.set_table(table)
        self.validate_col()
        self.validate_attr(append)
        self.set_attr()

    def validate_col(self, itemsize=None):
        """ validate this column: return the compared against itemsize """

        # validate this column for string truncation (or reset to the max size)
        dtype = getattr(self, 'dtype', None)
        if _ensure_decoded(self.kind) == u'string':

            c = self.col
            if c is not None:
                if itemsize is None:
                    itemsize = self.itemsize
                if c.itemsize < itemsize:
                    raise ValueError("Trying to store a string with len [%s] in [%s] column but\n"
                                     "this column has a limit of [%s]!\n"
                                     "Consider using min_itemsize to preset the sizes on these columns"
                                     % (itemsize,self.cname, c.itemsize))
                return c.itemsize

        return None

    def validate_attr(self, append):
        # check for backwards incompatibility
        if append:
            existing_kind = getattr(self.attrs, self.kind_attr, None)
            if existing_kind is not None and existing_kind != self.kind:
                raise TypeError("incompatible kind in col [%s - %s]" %
                                (existing_kind, self.kind))

    def update_info(self, info):
        """ set/update the info for this indexable with the key/value
            if there is a conflict raise/warn as needed """

        for key in self._info_fields:

            value = getattr(self,key,None)

            try:
                idx = info[self.name]
            except:
                idx = info[self.name] = dict()

            existing_value = idx.get(key)
            if key in idx and value is not None and existing_value != value:

                # frequency/name just warn
                if key in ['freq','index_name']:
                    ws = attribute_conflict_doc % (key,existing_value,value)
                    warnings.warn(ws, AttributeConflictWarning)

                    # reset
                    idx[key] = None
                    setattr(self,key,None)

                else:
                    raise ValueError("invalid info for [%s] for [%s]"""
                                     ", existing_value [%s] conflicts with new value [%s]" % (self.name,
                                                                                              key,existing_value,value))
            else:
                if value is not None or existing_value is not None:
                    idx[key] = value

        return self

    def set_info(self, info):
        """ set my state from the passed info """
        idx = info.get(self.name)
        if idx is not None:
            self.__dict__.update(idx)

    def get_attr(self):
        """ set the kind for this colummn """
        self.kind = getattr(self.attrs, self.kind_attr, None)

    def set_attr(self):
        """ set the kind for this colummn """
        setattr(self.attrs, self.kind_attr, self.kind)

class GenericIndexCol(IndexCol):
    """ an index which is not represented in the data of the table """

    @property
    def is_indexed(self):
        return False

    def convert(self, values, nan_rep, encoding):
        """ set the values from this selection: take = take ownership """

        self.values = Int64Index(np.arange(self.table.nrows))
        return self

    def get_attr(self):
        pass

    def set_attr(self):
        pass

class DataCol(IndexCol):
    """ a data holding column, by definition this is not indexable

        Parameters
        ----------

        data   : the actual data
        cname  : the column name in the table to hold the data (typeically values)
        """
    is_an_indexable = False
    is_data_indexable = False
    is_searchable = False
    _info_fields = ['tz']

    @classmethod
    def create_for_block(cls, i=None, name=None, cname=None, version=None, **kwargs):
        """ return a new datacol with the block i """

        if cname is None:
            cname = name or 'values_block_%d' % i
        if name is None:
            name = cname

        # prior to 0.10.1, we named values blocks like: values_block_0 an the
        # name values_0
        try:
            if version[0] == 0 and version[1] <= 10 and version[2] == 0:
                m = re.search("values_block_(\d+)", name)
                if m:
                    name = "values_%s" % m.groups()[0]
        except:
            pass

        return cls(name=name, cname=cname, **kwargs)

    def __init__(self, values=None, kind=None, typ=None, cname=None, data=None, block=None, **kwargs):
        super(DataCol, self).__init__(
            values=values, kind=kind, typ=typ, cname=cname, **kwargs)
        self.dtype = None
        self.dtype_attr = u"%s_dtype" % self.name
        self.set_data(data)

    def __unicode__(self):
        return "name->%s,cname->%s,dtype->%s,shape->%s" % (self.name, self.cname, self.dtype, self.shape)

    def __eq__(self, other):
        """ compare 2 col items """
        return all([getattr(self, a, None) == getattr(other, a, None) for a in ['name', 'cname', 'dtype', 'pos']])

    def set_data(self, data, dtype=None):
        self.data = data
        if data is not None:
            if dtype is not None:
                self.dtype = dtype
                self.set_kind()
            elif self.dtype is None:
                self.dtype = data.dtype.name
                self.set_kind()

    def take_data(self):
        """ return the data & release the memory """
        self.data, data = None, self.data
        return data

    def set_kind(self):
        # set my kind if we can
        if self.dtype is not None:
            dtype = _ensure_decoded(self.dtype)
            if dtype.startswith(u'string') or dtype.startswith(u'bytes'):
                self.kind = 'string'
            elif dtype.startswith(u'float'):
                self.kind = 'float'
            elif dtype.startswith(u'int') or dtype.startswith(u'uint'):
                self.kind = 'integer'
            elif dtype.startswith(u'date'):
                self.kind = 'datetime'
            elif dtype.startswith(u'bool'):
                self.kind = 'bool'
            else:
                raise AssertionError("cannot interpret dtype of [%s] in [%s]" % (dtype,self))

            # set my typ if we need
            if self.typ is None:
                self.typ = getattr(self.description,self.cname,None)

    def set_atom(self, block, existing_col, min_itemsize, nan_rep, info, encoding=None, **kwargs):
        """ create and setup my atom from the block b """

        self.values = list(block.items)
        dtype = block.dtype.name
        rvalues = block.values.ravel()
        inferred_type = lib.infer_dtype(rvalues)

        if inferred_type == 'datetime64':
            self.set_atom_datetime64(block)
        elif inferred_type == 'date':
            raise TypeError(
                "[date] is not implemented as a table column")
        elif inferred_type == 'datetime':
            if getattr(rvalues[0],'tzinfo',None) is not None:

                # if this block has more than one timezone, raise
                if len(set([r.tzinfo for r in rvalues])) != 1:
                    raise TypeError(
                        "too many timezones in this block, create separate data columns")

                # convert this column to datetime64[ns] utc, and save the tz
                index = DatetimeIndex(rvalues)
                tz = getattr(index,'tz',None)
                if tz is None:
                    raise TypeError(
                        "invalid timezone specification")

                values = index.tz_convert('UTC').values.view('i8')

                # store a converted timezone
                zone = tslib.get_timezone(index.tz)
                if zone is None:
                    zone = tslib.tot_seconds(index.tz.utcoffset())
                self.tz = zone

                self.update_info(info)
                self.set_atom_datetime64(block, values.reshape(block.values.shape))

            else:
                raise TypeError(
                    "[datetime] is not implemented as a table column")
        elif inferred_type == 'unicode':
            raise TypeError(
                "[unicode] is not implemented as a table column")

        # this is basically a catchall; if say a datetime64 has nans then will
        # end up here ###
        elif inferred_type == 'string' or dtype == 'object':
            self.set_atom_string(block, existing_col, min_itemsize, nan_rep, encoding)
        else:
            self.set_atom_data(block)

        return self

    def get_atom_string(self, block, itemsize):
        return _tables().StringCol(itemsize=itemsize, shape=block.shape[0])

    def set_atom_string(self, block, existing_col, min_itemsize, nan_rep, encoding):
        # fill nan items with myself
        block = block.fillna(nan_rep)
        data  = block.values

        # see if we have a valid string type
        inferred_type = lib.infer_dtype(data.ravel())
        if inferred_type != 'string':

            # we cannot serialize this data, so report an exception on a column by column basis
            for item in block.items:

                col = block.get(item)
                inferred_type = lib.infer_dtype(col.ravel())
                if inferred_type != 'string':
                    raise TypeError("Cannot serialize the column [%s] because\n"
                                    "its data contents are [%s] object dtype" %
                                    (item,inferred_type))


        # itemsize is the maximum length of a string (along any dimension)
        itemsize = lib.max_len_string_array(com._ensure_object(data.ravel()))

        # specified min_itemsize?
        if isinstance(min_itemsize, dict):
            min_itemsize = int(min_itemsize.get(
                self.name) or min_itemsize.get('values') or 0)
        itemsize = max(min_itemsize or 0, itemsize)

        # check for column in the values conflicts
        if existing_col is not None:
            eci = existing_col.validate_col(itemsize)
            if eci > itemsize:
                itemsize = eci

        self.itemsize = itemsize
        self.kind = 'string'
        self.typ = self.get_atom_string(block, itemsize)
        self.set_data(self.convert_string_data(data, itemsize, encoding))

    def convert_string_data(self, data, itemsize, encoding):
        return _convert_string_array(data, encoding, itemsize)

    def get_atom_coltype(self):
        """ return the PyTables column class for this column """
        if self.kind.startswith('uint'):
            col_name = "UInt%sCol" % self.kind[4:]
        else:
            col_name = "%sCol" % self.kind.capitalize()

        return getattr(_tables(), col_name)

    def get_atom_data(self, block):
        return self.get_atom_coltype()(shape=block.shape[0])

    def set_atom_data(self, block):
        self.kind = block.dtype.name
        self.typ = self.get_atom_data(block)
        self.set_data(block.values.astype(self.typ.type))

    def get_atom_datetime64(self, block):
        return _tables().Int64Col(shape=block.shape[0])

    def set_atom_datetime64(self, block, values = None):
        self.kind = 'datetime64'
        self.typ = self.get_atom_datetime64(block)
        if values is None:
            values = block.values.view('i8')
        self.set_data(values, 'datetime64')

    @property
    def shape(self):
        return getattr(self.data, 'shape', None)

    @property
    def cvalues(self):
        """ return my cython values """
        return self.data

    def validate_attr(self, append):
        """ validate that we have the same order as the existing & same dtype """
        if append:
            existing_fields = getattr(self.attrs, self.kind_attr, None)
            if (existing_fields is not None and
                    existing_fields != list(self.values)):
                raise ValueError("appended items do not match existing items"
                                " in table!")

            existing_dtype = getattr(self.attrs, self.dtype_attr, None)
            if (existing_dtype is not None and
                    existing_dtype != self.dtype):
                raise ValueError("appended items dtype do not match existing items dtype"
                                " in table!")

    def convert(self, values, nan_rep, encoding):
        """ set the data from this selection (and convert to the correct dtype if we can) """
        try:
            values = values[self.cname]
        except:
            pass
        self.set_data(values)

        # convert to the correct dtype
        if self.dtype is not None:
            dtype = _ensure_decoded(self.dtype)

            # reverse converts
            if dtype == u'datetime64':
                # recreate the timezone
                if self.tz is not None:

                    # data should be 2-dim here
                    # we stored as utc, so just set the tz

                    index = DatetimeIndex(self.data.ravel(),tz='UTC').tz_convert(self.tz)
                    self.data = np.array(index.tolist(),dtype=object).reshape(self.data.shape)

                else:
                    self.data = np.asarray(self.data, dtype='M8[ns]')

            elif dtype == u'date':
                self.data = np.array(
                    [date.fromtimestamp(v) for v in self.data], dtype=object)
            elif dtype == u'datetime':
                self.data = np.array(
                    [datetime.fromtimestamp(v) for v in self.data],
                    dtype=object)
            else:

                try:
                    self.data = self.data.astype(dtype)
                except:
                    self.data = self.data.astype('O')

        # convert nans / decode
        if _ensure_decoded(self.kind) == u'string':
            self.data = _unconvert_string_array(self.data, nan_rep=nan_rep, encoding=encoding)

        return self

    def get_attr(self):
        """ get the data for this colummn """
        self.values = getattr(self.attrs, self.kind_attr, None)
        self.dtype  = getattr(self.attrs, self.dtype_attr, None)
        self.set_kind()

    def set_attr(self):
        """ set the data for this colummn """
        setattr(self.attrs, self.kind_attr, self.values)
        if self.dtype is not None:
            setattr(self.attrs, self.dtype_attr, self.dtype)


class DataIndexableCol(DataCol):
    """ represent a data column that can be indexed """
    is_data_indexable = True

    @property
    def is_searchable(self):
        return _ensure_decoded(self.kind) == u'string'

    def get_atom_string(self, block, itemsize):
        return _tables().StringCol(itemsize=itemsize)

    def get_atom_data(self, block):
        return self.get_atom_coltype()()

    def get_atom_datetime64(self, block):
        return _tables().Int64Col()

class GenericDataIndexableCol(DataIndexableCol):
    """ represent a generic pytables data column """

    def get_attr(self):
        pass

class Storer(StringMixin):
    """ represent an object in my store
          facilitate read/write of various types of objects
          this is an abstract base class

        Parameters
        ----------

        parent : my parent HDFStore
        group  : the group node where the table resides
        """
    pandas_kind = None
    obj_type    = None
    ndim        = None
    is_table    = False

    def __init__(self, parent, group, encoding=None, **kwargs):
        self.parent      = parent
        self.group       = group
        self.encoding    = _ensure_encoding(encoding)
        self.set_version()

    @property
    def is_old_version(self):
        return self.version[0] <= 0 and self.version[1] <= 10 and self.version[2] < 1

    def set_version(self):
        """ compute and set our version """
        version = _ensure_decoded(getattr(self.group._v_attrs,'pandas_version',None))
        try:
            self.version = tuple([int(x) for x in version.split('.')])
            if len(self.version) == 2:
                self.version = self.version + (0,)
        except:
            self.version = (0, 0, 0)

    @property
    def pandas_type(self):
        return _ensure_decoded(getattr(self.group._v_attrs, 'pandas_type', None))

    def __unicode__(self):
        """ return a pretty representation of myself """
        self.infer_axes()
        s = self.shape
        if s is not None:
            if isinstance(s, (list,tuple)):
                s = "[%s]" % ','.join([pprint_thing(x) for x in s])
            return "%-12.12s (shape->%s)" % (self.pandas_type,s)
        return self.pandas_type

    def set_object_info(self):
        """ set my pandas type & version """
        self.attrs.pandas_type = str(self.pandas_kind)
        self.attrs.pandas_version = str(_version)
        self.set_version()

    def copy(self):
        new_self = copy.copy(self)
        return new_self

    @property
    def shape(self):
        return self.nrows

    @property
    def pathname(self):
        return self.group._v_pathname

    @property
    def _handle(self):
        return self.parent._handle

    @property
    def _quiet(self):
        return self.parent._quiet

    @property
    def _filters(self):
        return self.parent._filters

    @property
    def _complevel(self):
        return self.parent._complevel

    @property
    def _fletcher32(self):
        return self.parent._fletcher32

    @property
    def _complib(self):
        return self.parent._complib

    @property
    def attrs(self):
        return self.group._v_attrs

    def set_attrs(self):
        """ set our object attributes """
        pass

    def get_attrs(self):
        """ get our object attributes """
        pass

    @property
    def storable(self):
        """ return my storable """
        return self.group

    @property
    def is_exists(self):
        return False

    @property
    def nrows(self):
        return getattr(self.storable,'nrows',None)

    def validate(self, other):
        """ validate against an existing storable """
        if other is None: return
        return True

    def validate_version(self, where = None):
        """ are we trying to operate on an old version? """
        return True

    def infer_axes(self):
        """ infer the axes of my storer
              return a boolean indicating if we have a valid storer or not """

        s = self.storable
        if s is None:
            return False
        self.get_attrs()
        return True

    def read(self, **kwargs):
        raise NotImplementedError("cannot read on an abstract storer: subclasses should implement")

    def write(self, **kwargs):
        raise NotImplementedError("cannot write on an abstract storer: sublcasses should implement")

    def delete(self, where = None, **kwargs):
        """ support fully deleting the node in its entirety (only) - where specification must be None """
        if where is None:
            self._handle.removeNode(self.group, recursive=True)
            return None

        raise TypeError("cannot delete on an abstract storer")

class GenericStorer(Storer):
    """ a generified storer version """
    _index_type_map    = { DatetimeIndex: 'datetime',
                           PeriodIndex: 'period'}
    _reverse_index_map = dict([ (v,k) for k, v in _index_type_map.iteritems() ])
    attributes = []

    # indexer helpders
    def _class_to_alias(self, cls):
        return self._index_type_map.get(cls, '')

    def _alias_to_class(self, alias):
        if isinstance(alias, type):  # pragma: no cover
            return alias  # compat: for a short period of time master stored types
        return self._reverse_index_map.get(alias, Index)

    def _get_index_factory(self, klass):
        if klass == DatetimeIndex:
            def f(values, freq=None, tz=None):
                return DatetimeIndex._simple_new(values, None, freq=freq,
                                                 tz=tz)
            return f
        return klass

    def validate_read(self, kwargs):
        if kwargs.get('columns') is not None:
            raise TypeError("cannot pass a column specification when reading a non-table "
                            "this store must be selected in its entirety")
        if kwargs.get('where') is not None:
            raise TypeError("cannot pass a where specification when reading from a non-table "
                            "this store must be selected in its entirety")

    @property
    def is_exists(self):
        return True

    def set_attrs(self):
        """ set our object attributes """
        self.attrs.encoding = self.encoding

    def get_attrs(self):
        """ retrieve our attributes """
        self.encoding = _ensure_encoding(getattr(self.attrs,'encoding',None))
        for n in self.attributes:
            setattr(self,n,_ensure_decoded(getattr(self.attrs, n, None)))

    def write(self, obj, **kwargs):
        self.set_attrs()

    def read_array(self, key):
        """ read an array for the specified node (off of group """
        import tables
        node = getattr(self.group, key)
        data = node[:]
        attrs = node._v_attrs

        transposed = getattr(attrs, 'transposed', False)

        if isinstance(node, tables.VLArray):
            ret = data[0]
        else:
            dtype = getattr(attrs, 'value_type', None)
            shape = getattr(attrs, 'shape', None)

            if shape is not None:
                # length 0 axis
                ret = np.empty(shape, dtype=dtype)
            else:
                ret = data

            if dtype == u'datetime64':
                ret = np.array(ret, dtype='M8[ns]')

        if transposed:
            return ret.T
        else:
            return ret

    def read_index(self, key):
        variety = _ensure_decoded(getattr(self.attrs, '%s_variety' % key))

        if variety == u'multi':
            return self.read_multi_index(key)
        elif variety == u'block':
            return self.read_block_index(key)
        elif variety == u'sparseint':
            return self.read_sparse_intindex(key)
        elif variety == u'regular':
            _, index = self.read_index_node(getattr(self.group, key))
            return index
        else:  # pragma: no cover
            raise Exception('unrecognized index variety: %s' % variety)

    def write_index(self, key, index):
        if isinstance(index, MultiIndex):
            setattr(self.attrs, '%s_variety' % key, 'multi')
            self.write_multi_index(key, index)
        elif isinstance(index, BlockIndex):
            setattr(self.attrs, '%s_variety' % key, 'block')
            self.write_block_index(key, index)
        elif isinstance(index, IntIndex):
            setattr(self.attrs, '%s_variety' % key, 'sparseint')
            self.write_sparse_intindex(key, index)
        else:
            setattr(self.attrs, '%s_variety' % key, 'regular')
            converted = _convert_index(index,self.encoding).set_name('index')
            self.write_array(key, converted.values)
            node = getattr(self.group, key)
            node._v_attrs.kind = converted.kind
            node._v_attrs.name = index.name

            if isinstance(index, (DatetimeIndex, PeriodIndex)):
                node._v_attrs.index_class = self._class_to_alias(type(index))

            if hasattr(index, 'freq'):
                node._v_attrs.freq = index.freq

            if hasattr(index, 'tz') and index.tz is not None:
                zone = tslib.get_timezone(index.tz)
                if zone is None:
                    zone = tslib.tot_seconds(index.tz.utcoffset())
                node._v_attrs.tz = zone


    def write_block_index(self, key, index):
        self.write_array('%s_blocs' % key, index.blocs)
        self.write_array('%s_blengths' % key, index.blengths)
        setattr(self.attrs, '%s_length' % key, index.length)

    def read_block_index(self, key):
        length = getattr(self.attrs, '%s_length' % key)
        blocs = self.read_array('%s_blocs' % key)
        blengths = self.read_array('%s_blengths' % key)
        return BlockIndex(length, blocs, blengths)

    def write_sparse_intindex(self, key, index):
        self.write_array('%s_indices' % key, index.indices)
        setattr(self.attrs, '%s_length' % key, index.length)

    def read_sparse_intindex(self, key):
        length = getattr(self.attrs, '%s_length' % key)
        indices = self.read_array('%s_indices' % key)
        return IntIndex(length, indices)

    def write_multi_index(self, key, index):
        setattr(self.attrs, '%s_nlevels' % key, index.nlevels)

        for i, (lev, lab, name) in enumerate(zip(index.levels,
                                                 index.labels,
                                                 index.names)):
            # write the level
            level_key = '%s_level%d' % (key, i)
            conv_level = _convert_index(lev, self.encoding).set_name(level_key)
            self.write_array(level_key, conv_level.values)
            node = getattr(self.group, level_key)
            node._v_attrs.kind = conv_level.kind
            node._v_attrs.name = name

            # write the name
            setattr(node._v_attrs, '%s_name%d' % (key, i), name)

            # write the labels
            label_key = '%s_label%d' % (key, i)
            self.write_array(label_key, lab)

    def read_multi_index(self, key):
        nlevels = getattr(self.attrs, '%s_nlevels' % key)

        levels = []
        labels = []
        names = []
        for i in range(nlevels):
            level_key = '%s_level%d' % (key, i)
            name, lev = self.read_index_node(getattr(self.group, level_key))
            levels.append(lev)
            names.append(name)

            label_key = '%s_label%d' % (key, i)
            lab = self.read_array(label_key)
            labels.append(lab)

        return MultiIndex(levels=levels, labels=labels, names=names)

    def read_index_node(self, node):
        data = node[:]
        kind = _ensure_decoded(node._v_attrs.kind)
        name = None

        if 'name' in node._v_attrs:
            name = node._v_attrs.name

        index_class = self._alias_to_class(getattr(node._v_attrs,
                                                   'index_class', ''))
        factory = self._get_index_factory(index_class)

        kwargs = {}
        if u'freq' in node._v_attrs:
            kwargs['freq'] = node._v_attrs['freq']

        if u'tz' in node._v_attrs:
            kwargs['tz'] = node._v_attrs['tz']

        if kind in (u'date', u'datetime'):
            index = factory(_unconvert_index(data, kind, encoding=self.encoding), dtype=object,
                            **kwargs)
        else:
            index = factory(_unconvert_index(data, kind, encoding=self.encoding), **kwargs)

        index.name = name

        return name, index

    def write_array_empty(self, key, value):
        """ write a 0-len array """

        # ugly hack for length 0 axes
        arr = np.empty((1,) * value.ndim)
        self._handle.createArray(self.group, key, arr)
        getattr(self.group, key)._v_attrs.value_type = str(value.dtype)
        getattr(self.group, key)._v_attrs.shape = value.shape

    def write_array(self, key, value, items=None):
        if key in self.group:
            self._handle.removeNode(self.group, key)

        # Transform needed to interface with pytables row/col notation
        empty_array = any(x == 0 for x in value.shape)
        transposed = False

        if not empty_array:
            value = value.T
            transposed = True

        if self._filters is not None:
            atom = None
            try:
                # get the atom for this datatype
                atom = _tables().Atom.from_dtype(value.dtype)
            except ValueError:
                pass

            if atom is not None:
                # create an empty chunked array and fill it from value
                if not empty_array:
                    ca = self._handle.createCArray(self.group, key, atom,
                                                   value.shape,
                                                   filters=self._filters)
                    ca[:] = value
                    getattr(self.group, key)._v_attrs.transposed = transposed

                else:
                    self.write_array_empty(key, value)

                return

        if value.dtype.type == np.object_:

            # infer the type, warn if we have a non-string type here (for performance)
            inferred_type = lib.infer_dtype(value.ravel())
            if empty_array:
                pass
            elif inferred_type == 'string':
                pass
            else:
                try:
                    items = list(items)
                except:
                    pass
                ws = performance_doc % (inferred_type,key,items)
                warnings.warn(ws, PerformanceWarning)

            vlarr = self._handle.createVLArray(self.group, key,
                                              _tables().ObjectAtom())
            vlarr.append(value)
        elif value.dtype.type == np.datetime64:
            self._handle.createArray(self.group, key, value.view('i8'))
            getattr(self.group, key)._v_attrs.value_type = 'datetime64'
        else:
            if empty_array:
                self.write_array_empty(key, value)
            else:
                self._handle.createArray(self.group, key, value)

        getattr(self.group, key)._v_attrs.transposed = transposed

class LegacyStorer(GenericStorer):

    def read_index_legacy(self, key):
        node = getattr(self.group,key)
        data = node[:]
        kind = node._v_attrs.kind
        return _unconvert_index_legacy(data, kind, encoding=self.encoding)

class LegacySeriesStorer(LegacyStorer):

    def read(self, **kwargs):
        self.validate_read(kwargs)
        index = self.read_index_legacy('index')
        values = self.read_array('values')
        return Series(values, index=index)

class LegacyFrameStorer(LegacyStorer):

    def read(self, **kwargs):
        self.validate_read(kwargs)
        index = self.read_index_legacy('index')
        columns = self.read_index_legacy('columns')
        values = self.read_array('values')
        return DataFrame(values, index=index, columns=columns)

class SeriesStorer(GenericStorer):
    pandas_kind = u'series'
    attributes = ['name']

    @property
    def shape(self):
        try:
            return len(getattr(self.group,'values')),
        except:
            return None

    def read(self, **kwargs):
        self.validate_read(kwargs)
        index = self.read_index('index')
        if len(index) > 0:
            values = self.read_array('values')
        else:
            values = []

        return Series(values, index=index, name=self.name)

    def write(self, obj, **kwargs):
        super(SeriesStorer, self).write(obj, **kwargs)
        self.write_index('index', obj.index)
        self.write_array('values', obj.values)
        self.attrs.name = obj.name

class SparseSeriesStorer(GenericStorer):
    pandas_kind = u'sparse_series'
    attributes = ['name','fill_value','kind']

    def read(self, **kwargs):
        self.validate_read(kwargs)
        index = self.read_index('index')
        sp_values = self.read_array('sp_values')
        sp_index = self.read_index('sp_index')
        return SparseSeries(sp_values, index=index, sparse_index=sp_index,
                            kind=self.kind or u'block', fill_value=self.fill_value,
                            name=self.name)

    def write(self, obj, **kwargs):
        super(SparseSeriesStorer, self).write(obj, **kwargs)
        self.write_index('index', obj.index)
        self.write_index('sp_index', obj.sp_index)
        self.write_array('sp_values', obj.sp_values)
        self.attrs.name = obj.name
        self.attrs.fill_value = obj.fill_value
        self.attrs.kind = obj.kind

class SparseFrameStorer(GenericStorer):
    pandas_kind = u'sparse_frame'
    attributes = ['default_kind','default_fill_value']

    def read(self, **kwargs):
        self.validate_read(kwargs)
        columns = self.read_index('columns')
        sdict = {}
        for c in columns:
            key = 'sparse_series_%s' % c
            s = SparseSeriesStorer(self.parent, getattr(self.group,key))
            s.infer_axes()
            sdict[c] = s.read()
        return SparseDataFrame(sdict, columns=columns,
                               default_kind=self.default_kind,
                               default_fill_value=self.default_fill_value)

    def write(self, obj, **kwargs):
        """ write it as a collection of individual sparse series """
        super(SparseFrameStorer, self).write(obj, **kwargs)
        for name, ss in obj.iteritems():
            key = 'sparse_series_%s' % name
            if key not in self.group._v_children:
                node = self._handle.createGroup(self.group, key)
            else:
                node = getattr(self.group, key)
            s = SparseSeriesStorer(self.parent, node)
            s.write(ss)
        self.attrs.default_fill_value = obj.default_fill_value
        self.attrs.default_kind       = obj.default_kind
        self.write_index('columns', obj.columns)

class SparsePanelStorer(GenericStorer):
    pandas_kind = u'sparse_panel'
    attributes = ['default_kind','default_fill_value']

    def read(self, **kwargs):
        self.validate_read(kwargs)
        items = self.read_index('items')

        sdict = {}
        for name in items:
            key = 'sparse_frame_%s' % name
            node = getattr(self.group, key)
            s = SparseFrameStorer(self.parent, getattr(self.group,key))
            s.infer_axes()
            sdict[name] = s.read()
        return SparsePanel(sdict, items=items, default_kind=self.default_kind,
                           default_fill_value=self.default_fill_value)

    def write(self, obj, **kwargs):
        super(SparsePanelStorer, self).write(obj, **kwargs)
        self.attrs.default_fill_value = obj.default_fill_value
        self.attrs.default_kind       = obj.default_kind
        self.write_index('items', obj.items)

        for name, sdf in obj.iterkv():
            key = 'sparse_frame_%s' % name
            if key not in self.group._v_children:
                node = self._handle.createGroup(self.group, key)
            else:
                node = getattr(self.group, key)
            s = SparseFrameStorer(self.parent, node)
            s.write(sdf)

class BlockManagerStorer(GenericStorer):
    attributes = ['ndim','nblocks']
    is_shape_reversed = False

    @property
    def shape(self):
        try:
            ndim = self.ndim

            # items
            items = 0
            for i in range(self.nblocks):
                node = getattr(self.group, 'block%d_items' % i)
                shape = getattr(node,'shape',None)
                if shape is not None:
                    items += shape[0]

            # data shape
            node = getattr(self.group, 'block0_values')
            shape = getattr(node,'shape',None)
            if shape is not None:
                shape = list(shape[0:(ndim-1)])
            else:
                shape = []

            shape.append(items)

            # hacky - this works for frames, but is reversed for panels
            if self.is_shape_reversed:
                shape = shape[::-1]

            return shape
        except:
            return None

    def read(self, **kwargs):
        self.validate_read(kwargs)

        axes = []
        for i in xrange(self.ndim):
            ax = self.read_index('axis%d' % i)
            axes.append(ax)

        items = axes[0]
        blocks = []
        for i in range(self.nblocks):
            blk_items = self.read_index('block%d_items' % i)
            values = self.read_array('block%d_values' % i)
            blk = make_block(values, blk_items, items)
            blocks.append(blk)

        return self.obj_type(BlockManager(blocks, axes))

    def write(self, obj, **kwargs):
        super(BlockManagerStorer, self).write(obj, **kwargs)
        data = obj._data
        if not data.is_consolidated():
            data = data.consolidate()

        self.attrs.ndim = data.ndim
        for i, ax in enumerate(data.axes):
            self.write_index('axis%d' % i, ax)

        # Supporting mixed-type DataFrame objects...nontrivial
        self.attrs.nblocks = nblocks = len(data.blocks)
        for i in range(nblocks):
            blk = data.blocks[i]
            # I have no idea why, but writing values before items fixed #2299
            self.write_array('block%d_values' % i, blk.values, items=blk.items)
            self.write_index('block%d_items' % i, blk.items)

class FrameStorer(BlockManagerStorer):
    pandas_kind = u'frame'
    obj_type    = DataFrame

class PanelStorer(BlockManagerStorer):
    pandas_kind = u'wide'
    obj_type    = Panel
    is_shape_reversed = True

    def write(self, obj, **kwargs):
        obj._consolidate_inplace()
        return super(PanelStorer, self).write(obj, **kwargs)

class Table(Storer):
    """ represent a table:
          facilitate read/write of various types of tables

        Attrs in Table Node
        -------------------
        These are attributes that are store in the main table node, they are necessary
        to recreate these tables when read back in.

        index_axes    : a list of tuples of the (original indexing axis and index column)
        non_index_axes: a list of tuples of the (original index axis and columns on a non-indexing axis)
        values_axes   : a list of the columns which comprise the data of this table
        data_columns  : a list of the columns that we are allowing indexing (these become single columns in values_axes), or True to force all columns
        nan_rep       : the string to use for nan representations for string objects
        levels        : the names of levels

        """
    pandas_kind = u'wide_table'
    table_type  = None
    levels      = 1
    is_table    = True
    is_shape_reversed = False

    def __init__(self, *args, **kwargs):
        super(Table, self).__init__(*args, **kwargs)
        self.index_axes     = []
        self.non_index_axes = []
        self.values_axes    = []
        self.data_columns   = []
        self.info           = dict()
        self.nan_rep        = None
        self.selection      = None

    @property
    def table_type_short(self):
        return self.table_type.split('_')[0]

    def __unicode__(self):
        """ return a pretty representatgion of myself """
        self.infer_axes()
        dc = ",dc->[%s]" % ','.join(self.data_columns) if len(self.data_columns) else ''

        ver = ''
        if self.is_old_version:
            ver = "[%s]" % '.'.join([ str(x) for x in self.version ])

        return "%-12.12s%s (typ->%s,nrows->%s,ncols->%s,indexers->[%s]%s)" % (self.pandas_type,
                                                                              ver,
                                                                              self.table_type_short,
                                                                              self.nrows,
                                                                              self.ncols,
                                                                              ','.join([ a.name for a in self.index_axes ]),
                                                                              dc)

    def __getitem__(self, c):
        """ return the axis for c """
        for a in self.axes:
            if c == a.name:
                return a
        return None

    def validate(self, other):
        """ validate against an existing table """
        if other is None: return

        if other.table_type != self.table_type:
            raise TypeError("incompatible table_type with existing [%s - %s]" %
                            (other.table_type, self.table_type))

        for c in ['index_axes','non_index_axes','values_axes']:
            sv = getattr(self,c,None)
            ov = getattr(other,c,None)
            if sv != ov:

                # show the error for the specific axes
                for i, sax in enumerate(sv):
                    oax = ov[i]
                    if sax != oax:
                        raise ValueError("invalid combinate of [%s] on appending data [%s] vs current table [%s]" % (c,sax,oax))

                # should never get here
                raise Exception("invalid combinate of [%s] on appending data [%s] vs current table [%s]" % (c,sv,ov))

    @property
    def nrows_expected(self):
        """ based on our axes, compute the expected nrows """
        return np.prod([ i.cvalues.shape[0] for i in self.index_axes ])

    @property
    def is_exists(self):
        """ has this table been created """
        return u'table' in self.group

    @property
    def storable(self):
        return getattr(self.group,'table',None)

    @property
    def table(self):
        """ return the table group (this is my storable) """
        return self.storable

    @property
    def dtype(self):
        return self.table.dtype

    @property
    def description(self):
        return self.table.description

    @property
    def axes(self):
        return itertools.chain(self.index_axes, self.values_axes)

    @property
    def ncols(self):
        """ the number of total columns in the values axes """
        return sum([ len(a.values) for a in self.values_axes ])

    @property
    def is_transposed(self):
        return False

    @property
    def data_orientation(self):
        """ return a tuple of my permutated axes, non_indexable at the front """
        return tuple(itertools.chain([int(a[0]) for a in self.non_index_axes], [int(a.axis) for a in self.index_axes]))

    def queryables(self):
        """ return a dict of the kinds allowable columns for this object """

        # compute the values_axes queryables
        return dict([(a.cname, a.kind) for a in self.index_axes] +
                    [(self.obj_type._AXIS_NAMES[axis], None) for axis, values in self.non_index_axes] +
                    [(v.cname, v.kind) for v in self.values_axes if v.name in set(self.data_columns)]
                    )

    def index_cols(self):
        """ return a list of my index cols """
        return [(i.axis, i.cname) for i in self.index_axes]

    def values_cols(self):
        """ return a list of my values cols """
        return [i.cname for i in self.values_axes]

    def set_info(self):
        """ update our table index info """
        self.attrs.info         = self.info

    def set_attrs(self):
        """ set our table type & indexables """
        self.attrs.table_type   = str(self.table_type)
        self.attrs.index_cols   = self.index_cols()
        self.attrs.values_cols  = self.values_cols()
        self.attrs.non_index_axes = self.non_index_axes
        self.attrs.data_columns = self.data_columns
        self.attrs.nan_rep      = self.nan_rep
        self.attrs.encoding     = self.encoding
        self.attrs.levels       = self.levels
        self.set_info()

    def get_attrs(self):
        """ retrieve our attributes """
        self.non_index_axes   = getattr(self.attrs,'non_index_axes',None) or []
        self.data_columns     = getattr(self.attrs,'data_columns',None)   or []
        self.info             = getattr(self.attrs,'info',None)           or dict()
        self.nan_rep          = getattr(self.attrs,'nan_rep',None)
        self.encoding         = _ensure_encoding(getattr(self.attrs,'encoding',None))
        self.levels           = getattr(self.attrs,'levels',None)         or []
        t = self.table
        self.index_axes       = [ a.infer(t) for a in self.indexables if     a.is_an_indexable ]
        self.values_axes      = [ a.infer(t) for a in self.indexables if not a.is_an_indexable ]

    def validate_version(self, where = None):
        """ are we trying to operate on an old version? """
        if where is not None:
            if self.version[0] <= 0 and self.version[1] <= 10 and self.version[2] < 1:
                ws = incompatibility_doc % '.'.join([ str(x) for x in self.version ])
                warnings.warn(ws, IncompatibilityWarning)

    def validate_min_itemsize(self, min_itemsize):
        """ validate the min_itemisze doesn't contain items that are not in the axes
            this needs data_columns to be defined """
        if min_itemsize is None: return
        if not isinstance(min_itemsize, dict): return

        q = self.queryables()
        for k, v in min_itemsize.items():

            # ok, apply generally
            if k == 'values':
                continue
            if k not in q:
                raise ValueError("min_itemsize has the key [%s] which is not an axis or data_column" % k)

    @property
    def indexables(self):
        """ create/cache the indexables if they don't exist """
        if self._indexables is None:

            d = self.description
            self._indexables = []

            # index columns
            self._indexables.extend([ IndexCol(name=name,axis=axis,pos=i) for i, (axis, name) in enumerate(self.attrs.index_cols)])

            # values columns
            dc = set(self.data_columns)
            base_pos = len(self._indexables)

            def f(i, c):
                klass = DataCol
                if c in dc:
                    klass = DataIndexableCol
                return klass.create_for_block(i=i, name=c, pos=base_pos + i, version=self.version)

            self._indexables.extend(
                [f(i, c) for i, c in enumerate(self.attrs.values_cols)])

        return self._indexables

    def create_index(self, columns=None, optlevel=None, kind=None):
        """
        Create a pytables index on the specified columns
          note: cannot index Time64Col() currently; PyTables must be >= 2.3


        Paramaters
        ----------
        columns : False (don't create an index), True (create all columns index), None or list_like (the indexers to index)
        optlevel: optimization level (defaults to 6)
        kind    : kind of index (defaults to 'medium')

        Exceptions
        ----------
        raises if the node is not a table

        """

        if not self.infer_axes():
            return
        if columns is False:
            return

        # index all indexables and data_columns
        if columns is None or columns is True:
            columns = [a.cname for a in self.axes if a.is_data_indexable]
        if not isinstance(columns, (tuple, list)):
            columns = [columns]

        kw = dict()
        if optlevel is not None:
            kw['optlevel'] = optlevel
        if kind is not None:
            kw['kind'] = kind

        table = self.table
        for c in columns:
            v = getattr(table.cols, c, None)
            if v is not None:

                # remove the index if the kind/optlevel have changed
                if v.is_indexed:
                    index = v.index
                    cur_optlevel = index.optlevel
                    cur_kind = index.kind

                    if kind is not None and cur_kind != kind:
                        v.removeIndex()
                    else:
                        kw['kind'] = cur_kind

                    if optlevel is not None and cur_optlevel != optlevel:
                        v.removeIndex()
                    else:
                        kw['optlevel'] = cur_optlevel

                # create the index
                if not v.is_indexed:
                    v.createIndex(**kw)

    def read_axes(self, where, **kwargs):
        """ create and return the axes sniffed from the table: return boolean for success """

        # validate the version
        self.validate_version(where)

        # infer the data kind
        if not self.infer_axes():
            return False

        # create the selection
        self.selection = Selection(self, where=where, **kwargs)
        values = self.selection.select()

        # convert the data
        for a in self.axes:
            a.set_info(self.info)
            a.convert(values, nan_rep=self.nan_rep, encoding=self.encoding)

        return True

    def get_object(self, obj):
        """ return the data for this obj """
        return obj

    def validate_data_columns(self, data_columns, min_itemsize):
        """ take the input data_columns and min_itemize and create a data_columns spec """

        if not len(self.non_index_axes):
            return []

        axis_labels = self.non_index_axes[0][1]

        # evaluate the passed data_columns, True == use all columns
        # take only valide axis labels
        if data_columns is True:
            data_columns = axis_labels
        elif data_columns is None:
            data_columns = []

        # if min_itemsize is a dict, add the keys (exclude 'values')
        if isinstance(min_itemsize,dict):

            existing_data_columns = set(data_columns)
            data_columns.extend([ k for k in min_itemsize.keys() if k != 'values' and k not in existing_data_columns ])

        # return valid columns in the order of our axis
        return [c for c in data_columns if c in axis_labels]

    def create_axes(self, axes, obj, validate=True, nan_rep=None, data_columns=None, min_itemsize=None, **kwargs):
        """ create and return the axes
              leagcy tables create an indexable column, indexable index, non-indexable fields

            Parameters:
            -----------
            axes: a list of the axes in order to create (names or numbers of the axes)
            obj : the object to create axes on
            validate: validate the obj against an existiing object already written
            min_itemsize: a dict of the min size for a column in bytes
            nan_rep : a values to use for string column nan_rep
            encoding : the encoding for string values
            data_columns : a list of columns that we want to create separate to allow indexing (or True will force all colummns)

        """

        # set the default axes if needed
        if axes is None:
            try:
                axes = _AXES_MAP[type(obj)]
            except:
                raise TypeError("cannot properly create the storer for: [group->%s,value->%s]" %
                                (self.group._v_name,type(obj)))

        # map axes to numbers
        axes = [obj._get_axis_number(a) for a in axes]

        # do we have an existing table (if so, use its axes & data_columns)
        if self.infer_axes():
            existing_table = self.copy()
            existing_table.infer_axes()
            axes          = [ a.axis for a in existing_table.index_axes]
            data_columns  = existing_table.data_columns
            nan_rep       = existing_table.nan_rep
            self.encoding = existing_table.encoding
            self.info     = copy.copy(existing_table.info)
        else:
            existing_table = None

        # currently support on ndim-1 axes
        if len(axes) != self.ndim - 1:
            raise ValueError("currently only support ndim-1 indexers in an AppendableTable")

        # create according to the new data
        self.non_index_axes = []
        self.data_columns = []

        # nan_representation
        if nan_rep is None:
            nan_rep = 'nan'

        self.nan_rep = nan_rep

        # create axes to index and non_index
        index_axes_map = dict()
        for i, a in enumerate(obj.axes):

            if i in axes:
                name = obj._AXIS_NAMES[i]
                index_axes_map[i] = _convert_index(
                    a, self.encoding).set_name(name).set_axis(i)
            else:

                # we might be able to change the axes on the appending data if
                # necessary
                append_axis = list(a)
                if existing_table is not None:
                    indexer = len(self.non_index_axes)
                    exist_axis = existing_table.non_index_axes[indexer][1]
                    if append_axis != exist_axis:

                        # ahah! -> reindex
                        if sorted(append_axis) == sorted(exist_axis):
                            append_axis = exist_axis

                self.non_index_axes.append((i, append_axis))

        # set axis positions (based on the axes)
        self.index_axes = [index_axes_map[a].set_pos(j).update_info(self.info) for j,
                           a in enumerate(axes)]
        j = len(self.index_axes)

        # check for column conflicts
        if validate:
            for a in self.axes:
                a.maybe_set_size(min_itemsize=min_itemsize)

        # reindex by our non_index_axes & compute data_columns
        for a in self.non_index_axes:
            obj = obj.reindex_axis(a[1], axis=a[0], copy=False)

        # figure out data_columns and get out blocks
        block_obj = self.get_object(obj).consolidate()
        blocks = block_obj._data.blocks
        if len(self.non_index_axes):
            axis, axis_labels = self.non_index_axes[0]
            data_columns = self.validate_data_columns(data_columns, min_itemsize)
            if len(data_columns):
                blocks = block_obj.reindex_axis(Index(axis_labels) - Index(
                        data_columns), axis=axis, copy=False)._data.blocks
                for c in data_columns:
                    blocks.extend(block_obj.reindex_axis(
                            [c], axis=axis, copy=False)._data.blocks)

        # reorder the blocks in the same order as the existing_table if we can
        if existing_table is not None:
            by_items = dict([ (tuple(b.items.tolist()),b) for b in blocks ])
            new_blocks = []
            for ea in existing_table.values_axes:
                items = tuple(ea.values)
                try:
                    b = by_items.pop(items)
                    new_blocks.append(b)
                except:
                    raise ValueError("cannot match existing table structure for [%s] on appending data" % items)
            blocks = new_blocks

        # add my values
        self.values_axes = []
        for i, b in enumerate(blocks):

            # shape of the data column are the indexable axes
            klass = DataCol
            name = None

            # we have a data_column
            if data_columns and len(b.items) == 1 and b.items[0] in data_columns:
                klass = DataIndexableCol
                name = b.items[0]
                self.data_columns.append(name)

            # make sure that we match up the existing columns
            # if we have an existing table
            if existing_table is not None and validate:
                try:
                    existing_col = existing_table.values_axes[i]
                except:
                    raise ValueError("Incompatible appended table [%s] with existing table [%s]" %
                                    (blocks,existing_table.values_axes))
            else:
                existing_col = None

            try:
                col = klass.create_for_block(
                    i=i, name=name, version=self.version)
                col.set_atom(block=b,
                             existing_col=existing_col,
                             min_itemsize=min_itemsize,
                             nan_rep=nan_rep,
                             encoding=self.encoding,
                             info=self.info,
                             **kwargs)
                col.set_pos(j)

                self.values_axes.append(col)
            except (NotImplementedError, ValueError, TypeError), e:
                raise e
            except (Exception), detail:
                raise Exception("cannot find the correct atom type -> [dtype->%s,items->%s] %s" % (b.dtype.name, b.items, str(detail)))
            j += 1

        # validate our min_itemsize
        self.validate_min_itemsize(min_itemsize)

        # validate the axes if we have an existing table
        if validate:
            self.validate(existing_table)

    def process_axes(self, obj, columns=None):
        """ process axes filters """

        # reorder by any non_index_axes & limit to the select columns
        for axis, labels in self.non_index_axes:
            if columns is not None:
                labels = Index(labels) & Index(columns)
            obj = obj.reindex_axis(labels, axis=axis, copy=False)

        # apply the selection filters (but keep in the same order)
        if self.selection.filter:
            for field, op, filt in self.selection.filter:

                def process_filter(field, filt):

                    for axis_name in obj._AXIS_NAMES.values():
                        axis_number = obj._get_axis_number(axis_name)
                        axis_values = obj._get_axis(axis_name)

                        # see if the field is the name of an axis
                        if field == axis_name:
                            takers = op(axis_values,filt)
                            return obj.ix._getitem_axis(takers,axis=axis_number)

                        # this might be the name of a file IN an axis
                        elif field in axis_values:

                            # we need to filter on this dimension
                            values = _ensure_index(getattr(obj,field).values)
                            filt   = _ensure_index(filt)

                            # hack until we support reversed dim flags
                            if isinstance(obj,DataFrame):
                                axis_number = 1-axis_number
                            takers = op(values,filt)
                            return obj.ix._getitem_axis(takers,axis=axis_number)

                    raise ValueError("cannot find the field [%s] for filtering!" % field)

                obj = process_filter(field, filt)

        return obj

    def create_description(self, complib=None, complevel=None, fletcher32=False, expectedrows=None):
        """ create the description of the table from the axes & values """

        # expected rows estimate
        if expectedrows is None:
            expectedrows = max(self.nrows_expected, 10000)
        d = dict(name='table', expectedrows=expectedrows)

        # description from the axes & values
        d['description'] = dict([(a.cname, a.typ) for a in self.axes])

        if complib:
            if complevel is None:
                complevel = self._complevel or 9
            filters = _tables().Filters(complevel=complevel,
                                        complib=complib,
                                        fletcher32=fletcher32 or self._fletcher32)
            d['filters'] = filters
        elif self._filters is not None:
            d['filters'] = self._filters

        return d

    def read_coordinates(self, where=None, start=None, stop=None, **kwargs):
        """ select coordinates (row numbers) from a table; return the coordinates object """

        # validate the version
        self.validate_version(where)

        # infer the data kind
        if not self.infer_axes():
            return False

        # create the selection
        self.selection = Selection(self, where=where, start=start, stop=stop, **kwargs)
        return Coordinates(self.selection.select_coords(), group=self.group, where=where)

    def read_column(self, column, where = None, **kwargs):
        """ return a single column from the table, generally only indexables are interesting """

        # validate the version
        self.validate_version()

        # infer the data kind
        if not self.infer_axes():
            return False

        if where is not None:
            raise Exception("read_column does not currently accept a where clause")

        # find the axes
        for a in self.axes:
            if column == a.name:

                if not a.is_data_indexable:
                    raise ValueError("column [%s] can not be extracted individually; it is not data indexable" % column)

                # column must be an indexable or a data column
                c = getattr(self.table.cols, column)
                a.set_info(self.info)
                return Series(a.convert(c[:], nan_rep=self.nan_rep, encoding=self.encoding).take_data())

        raise KeyError("column [%s] not found in the table" % column)

class WORMTable(Table):
    """ a write-once read-many table: this format DOES NOT ALLOW appending to a
         table. writing is a one-time operation the data are stored in a format
         that allows for searching the data on disk
         """
    table_type = u'worm'

    def read(self, **kwargs):
        """ read the indicies and the indexing array, calculate offset rows and
        return """
        raise NotImplementedError("WORMTable needs to implement read")

    def write(self, **kwargs):
        """ write in a format that we can search later on (but cannot append
               to): write out the indicies and the values using _write_array
               (e.g. a CArray) create an indexing table so that we can search"""
        raise NotImplementedError("WORKTable needs to implement write")


class LegacyTable(Table):
    """ an appendable table: allow append/query/delete operations to a
          (possibily) already existing appendable table this table ALLOWS
          append (but doesn't require them), and stores the data in a format
          that can be easily searched

        """
    _indexables = [IndexCol(name='index', axis=1, pos=0),
                   IndexCol(name='column', axis=2,
                            pos=1, index_kind='columns_kind'),
                   DataCol(name='fields', cname='values', kind_attr='fields', pos=2)]
    table_type = u'legacy'
    ndim = 3

    def write(self, **kwargs):
        raise Exception("write operations are not allowed on legacy tables!")

    def read(self, where=None, columns=None, **kwargs):
        """ we have n indexable columns, with an arbitrary number of data axes """

        if not self.read_axes(where=where, **kwargs):
            return None

        factors = [Categorical.from_array(a.values) for a in self.index_axes]
        levels = [f.levels for f in factors]
        N = [len(f.levels) for f in factors]
        labels = [f.labels for f in factors]

        # compute the key
        key = factor_indexer(N[1:], labels)

        objs = []
        if len(unique(key)) == len(key):

            sorter, _ = algos.groupsort_indexer(
                com._ensure_int64(key), np.prod(N))
            sorter = com._ensure_platform_int(sorter)

            # create the objs
            for c in self.values_axes:

                # the data need to be sorted
                sorted_values = c.take_data().take(sorter, axis=0)

                take_labels = [l.take(sorter) for l in labels]
                items = Index(c.values)
                block = block2d_to_blocknd(
                    sorted_values, items, tuple(N), take_labels)

                # create the object
                mgr = BlockManager([block], [items] + levels)
                obj = self.obj_type(mgr)

                # permute if needed
                if self.is_transposed:
                    obj = obj.transpose(
                        *tuple(Series(self.data_orientation).argsort()))

                objs.append(obj)

        else:
            if not self._quiet:  # pragma: no cover
                print ('Duplicate entries in table, taking most recently '
                       'appended')

            # reconstruct
            long_index = MultiIndex.from_arrays(
                [i.values for i in self.index_axes])

            for c in self.values_axes:
                lp = DataFrame(c.data, index=long_index, columns=c.values)

                # need a better algorithm
                tuple_index = long_index._tuple_index

                unique_tuples = lib.fast_unique(tuple_index)
                unique_tuples = _asarray_tuplesafe(unique_tuples)

                indexer = match(unique_tuples, tuple_index)
                indexer = com._ensure_platform_int(indexer)

                new_index = long_index.take(indexer)
                new_values = lp.values.take(indexer, axis=0)

                lp = DataFrame(new_values, index=new_index, columns=lp.columns)
                objs.append(lp.to_panel())

        # create the composite object
        if len(objs) == 1:
            wp = objs[0]
        else:
            wp = concat(objs, axis=0, verify_integrity=True)

        # apply the selection filters & axis orderings
        wp = self.process_axes(wp, columns=columns)

        return wp


class LegacyFrameTable(LegacyTable):
    """ support the legacy frame table """
    pandas_kind = u'frame_table'
    table_type = u'legacy_frame'
    obj_type = Panel

    def read(self, *args, **kwargs):
        return super(LegacyFrameTable, self).read(*args, **kwargs)['value']


class LegacyPanelTable(LegacyTable):
    """ support the legacy panel table """
    table_type = u'legacy_panel'
    obj_type = Panel


class AppendableTable(LegacyTable):
    """ suppor the new appendable table formats """
    _indexables = None
    table_type = u'appendable'

    def write(self, obj, axes=None, append=False, complib=None,
              complevel=None, fletcher32=None, min_itemsize=None, chunksize=None,
              expectedrows=None, **kwargs):

        if not append and self.is_exists:
            self._handle.removeNode(self.group, 'table')

        # create the axes
        self.create_axes(axes=axes, obj=obj, validate=append,
                         min_itemsize=min_itemsize,
                         **kwargs)

        if not self.is_exists:

            # create the table
            options = self.create_description(complib=complib,
                                              complevel=complevel,
                                              fletcher32=fletcher32,
                                              expectedrows=expectedrows)

            # set the table attributes
            self.set_attrs()

            # create the table
            table = self._handle.createTable(self.group, **options)

        else:
            table = self.table

        # update my info
        self.set_info()

        # validate the axes and set the kinds
        for a in self.axes:
            a.validate_and_set(table, append)

        # add the rows
        self.write_data(chunksize)

    def write_data(self, chunksize):
        """ fast writing of data: requires specific cython routines each axis shape """

        # create the masks & values
        masks = []
        for a in self.values_axes:

            # figure the mask: only do if we can successfully process this
            # column, otherwise ignore the mask
            mask = com.isnull(a.data).all(axis=0)
            masks.append(mask.astype('u1'))

        # consolidate masks
        mask = masks[0]
        for m in masks[1:]:
            mask = mask & m

        # the arguments
        indexes = [a.cvalues for a in self.index_axes]
        search = np.array(
            [a.is_searchable for a in self.values_axes]).astype('u1')
        values = [a.take_data() for a in self.values_axes]

        # transpose the values so first dimension is last
        values = [ v.transpose(np.roll(np.arange(v.ndim),v.ndim-1)) for v in values ]

        # write the chunks
        if chunksize is None:
            chunksize = 100000

        rows = self.nrows_expected
        chunks = int(rows / chunksize) + 1
        for i in xrange(chunks):
            start_i = i * chunksize
            end_i = min((i + 1) * chunksize, rows)
            if start_i >= end_i:
                break

            self.write_data_chunk(
                indexes=[a[start_i:end_i] for a in indexes],
                mask=mask[start_i:end_i],
                search=search,
                values=[v[start_i:end_i] for v in values])

    def write_data_chunk(self, indexes, mask, search, values):

        # 0 len
        for v in values:
            if not np.prod(v.shape):
                return

        # get our function
        try:
            func = getattr(lib, "create_hdf_rows_%sd" % self.ndim)
            args = list(indexes)
            args.extend([self.dtype, mask, search, values])
            rows = func(*args)
        except (Exception), detail:
            raise Exception("cannot create row-data -> %s" % str(detail))

        try:
            if len(rows):
                self.table.append(rows)
                self.table.flush()
        except (Exception), detail:
            raise Exception("tables cannot write this data -> %s" % str(detail))

    def delete(self, where=None, **kwargs):

        # delete all rows (and return the nrows)
        if where is None or not len(where):
            nrows = self.nrows
            self._handle.removeNode(self.group, recursive=True)
            return nrows

        # infer the data kind
        if not self.infer_axes():
            return None

        # create the selection
        table = self.table
        self.selection = Selection(self, where, **kwargs)
        values = self.selection.select_coords()

        # delete the rows in reverse order
        l = Series(values).order()
        ln = len(l)

        if ln:

            # construct groups of consecutive rows
            diff = l.diff()
            groups = list(diff[diff > 1].index)

            # 1 group
            if not len(groups):
                groups = [0]

            # final element
            if groups[-1] != ln:
                groups.append(ln)

            # initial element
            if groups[0] != 0:
                groups.insert(0, 0)

            # we must remove in reverse order!
            pg = groups.pop()
            for g in reversed(groups):
                rows = l.take(range(g, pg))
                table.removeRows(start=rows[rows.index[0]
                                            ], stop=rows[rows.index[-1]] + 1)
                pg = g

            self.table.flush()

        # return the number of rows removed
        return ln


class AppendableFrameTable(AppendableTable):
    """ suppor the new appendable table formats """
    pandas_kind = u'frame_table'
    table_type = u'appendable_frame'
    ndim = 2
    obj_type = DataFrame

    @property
    def is_transposed(self):
        return self.index_axes[0].axis == 1

    def get_object(self, obj):
        """ these are written transposed """
        if self.is_transposed:
            obj = obj.T
        return obj

    def read(self, where=None, columns=None, **kwargs):

        if not self.read_axes(where=where, **kwargs):
            return None

        index = self.index_axes[0].values
        frames = []
        for a in self.values_axes:
            cols = Index(a.values)

            if self.is_transposed:
                values = a.cvalues
                index_ = cols
                cols_ = Index(index,name=getattr(index,'name',None))
            else:
                values = a.cvalues.T
                index_ = Index(index,name=getattr(index,'name',None))
                cols_ = cols

            # if we have a DataIndexableCol, its shape will only be 1 dim
            if values.ndim == 1:
                values = values.reshape(1, values.shape[0])

            block = make_block(values, cols_, cols_)
            mgr = BlockManager([block], [cols_, index_])
            frames.append(DataFrame(mgr))

        if len(frames) == 1:
            df = frames[0]
        else:
            df = concat(frames, axis=1, verify_integrity=True)

        # apply the selection filters & axis orderings
        df = self.process_axes(df, columns=columns)

        return df


class GenericTable(AppendableFrameTable):
    """ a table that read/writes the generic pytables table format """
    pandas_kind = u'frame_table'
    table_type = u'generic_table'
    ndim = 2
    obj_type = DataFrame

    @property
    def pandas_type(self):
        return self.pandas_kind

    @property
    def storable(self):
        return getattr(self.group,'table',None) or self.group

    def get_attrs(self):
        """ retrieve our attributes """
        self.non_index_axes   = []
        self.nan_rep          = None
        self.levels           = []
        t = self.table
        self.index_axes       = [ a.infer(t) for a in self.indexables if     a.is_an_indexable ]
        self.values_axes      = [ a.infer(t) for a in self.indexables if not a.is_an_indexable ]
        self.data_columns     = [ a.name for a in self.values_axes ]

    @property
    def indexables(self):
        """ create the indexables from the table description """
        if self._indexables is None:

            d = self.description

            # the index columns is just a simple index
            self._indexables = [ GenericIndexCol(name='index',axis=0) ]

            for i, n in enumerate(d._v_names):

                dc = GenericDataIndexableCol(name = n, pos=i, values = [ n ], version = self.version)
                self._indexables.append(dc)

        return self._indexables

    def write(self, **kwargs):
        raise NotImplementedError("cannot write on an generic table")

class AppendableMultiFrameTable(AppendableFrameTable):
    """ a frame with a multi-index """
    table_type = u'appendable_multiframe'
    obj_type = DataFrame
    ndim = 2

    @property
    def table_type_short(self):
        return u'appendable_multi'

    def write(self, obj, data_columns=None, **kwargs):
        if data_columns is None:
            data_columns = []
        elif data_columns is True:
            data_columns = obj.columns[:]
        for n in obj.index.names:
            if n not in data_columns:
                data_columns.insert(0, n)
        self.levels = obj.index.names
        return super(AppendableMultiFrameTable, self).write(obj=obj.reset_index(), data_columns=data_columns, **kwargs)

    def read(self, columns=None, **kwargs):
        if columns is not None:
            for n in self.levels:
                if n not in columns:
                    columns.insert(0, n)
        df = super(AppendableMultiFrameTable, self).read(columns=columns, **kwargs)
        df.set_index(self.levels, inplace=True)
        return df


class AppendablePanelTable(AppendableTable):
    """ suppor the new appendable table formats """
    table_type = u'appendable_panel'
    ndim = 3
    obj_type = Panel

    def get_object(self, obj):
        """ these are written transposed """
        if self.is_transposed:
            obj = obj.transpose(*self.data_orientation)
        return obj

    @property
    def is_transposed(self):
        return self.data_orientation != tuple(range(self.ndim))


class AppendableNDimTable(AppendablePanelTable):
    """ suppor the new appendable table formats """
    table_type = u'appendable_ndim'
    ndim = 4
    obj_type = Panel4D

def _convert_index(index, encoding=None):
    index_name = getattr(index,'name',None)

    if isinstance(index, DatetimeIndex):
        converted = index.asi8
        return IndexCol(converted, 'datetime64', _tables().Int64Col(),
                        freq=getattr(index,'freq',None), tz=getattr(index,'tz',None),
                        index_name=index_name)
    elif isinstance(index, (Int64Index, PeriodIndex)):
        atom = _tables().Int64Col()
        return IndexCol(index.values, 'integer', atom, freq=getattr(index,'freq',None),
                        index_name=index_name)

    if isinstance(index, MultiIndex):
        raise Exception('MultiIndex not supported here!')

    inferred_type = lib.infer_dtype(index)

    values = np.asarray(index)

    if inferred_type == 'datetime64':
        converted = values.view('i8')
        return IndexCol(converted, 'datetime64', _tables().Int64Col(),
                        freq=getattr(index,'freq',None), tz=getattr(index,'tz',None),
                        index_name=index_name)
    elif inferred_type == 'datetime':
        converted = np.array([(time.mktime(v.timetuple()) +
                               v.microsecond / 1E6) for v in values],
                             dtype=np.float64)
        return IndexCol(converted, 'datetime', _tables().Time64Col(),
                        index_name=index_name)
    elif inferred_type == 'date':
        converted = np.array([time.mktime(v.timetuple()) for v in values],
                             dtype=np.int32)
        return IndexCol(converted, 'date', _tables().Time32Col(),
                        index_name=index_name)
    elif inferred_type == 'string':
        # atom = _tables().ObjectAtom()
        # return np.asarray(values, dtype='O'), 'object', atom

        converted = _convert_string_array(values, encoding)
        itemsize = converted.dtype.itemsize
        return IndexCol(converted, 'string', _tables().StringCol(itemsize), itemsize=itemsize,
                        index_name=index_name)
    elif inferred_type == 'unicode':
        atom = _tables().ObjectAtom()
        return IndexCol(np.asarray(values, dtype='O'), 'object', atom,
                        index_name=index_name)
    elif inferred_type == 'integer':
        # take a guess for now, hope the values fit
        atom = _tables().Int64Col()
        return IndexCol(np.asarray(values, dtype=np.int64), 'integer', atom,
                        index_name=index_name)
    elif inferred_type == 'floating':
        atom = _tables().Float64Col()
        return IndexCol(np.asarray(values, dtype=np.float64), 'float', atom,
                        index_name=index_name)
    else:  # pragma: no cover
        atom = _tables().ObjectAtom()
        return IndexCol(np.asarray(values, dtype='O'), 'object', atom,
                        index_name=index_name)

def _unconvert_index(data, kind, encoding=None):
    kind = _ensure_decoded(kind)
    if kind == u'datetime64':
        index = DatetimeIndex(data)
    elif kind == u'datetime':
        index = np.array([datetime.fromtimestamp(v) for v in data],
                         dtype=object)
    elif kind == u'date':
        index = np.array([date.fromtimestamp(v) for v in data], dtype=object)
    elif kind in (u'integer', u'float'):
        index = np.array(data)
    elif kind in (u'string'):
        index = _unconvert_string_array(data, nan_rep=None, encoding=encoding)
    elif kind == u'object':
        index = np.array(data[0])
    else:  # pragma: no cover
        raise ValueError('unrecognized index type %s' % kind)
    return index

def _unconvert_index_legacy(data, kind, legacy=False, encoding=None):
    kind = _ensure_decoded(kind)
    if kind == u'datetime':
        index = lib.time64_to_datetime(data)
    elif kind in (u'integer'):
        index = np.array(data, dtype=object)
    elif kind in (u'string'):
        index = _unconvert_string_array(data, nan_rep=None, encoding=encoding)
    else:  # pragma: no cover
        raise ValueError('unrecognized index type %s' % kind)
    return index

def _convert_string_array(data, encoding, itemsize=None):

    # encode if needed
    if encoding is not None and len(data):
        f = np.vectorize(lambda x: x.encode(encoding), otypes=[np.object])
        data = f(data)

    # create the sized dtype
    if itemsize is None:
        itemsize = lib.max_len_string_array(com._ensure_object(data.ravel()))

    data = np.array(data,dtype="S%d" % itemsize)
    return data

def _unconvert_string_array(data, nan_rep=None, encoding=None):
    """ deserialize a string array, possibly decoding """
    shape = data.shape
    data = np.array(data.ravel(),dtype=object)

    # guard against a None encoding in PY3 (because of a legacy
    # where the passed encoding is actually None)
    encoding = _ensure_encoding(encoding)
    if encoding is not None and len(data):
        f = np.vectorize(lambda x: x.decode(encoding),otypes=[np.object])
        data = f(data)

    if nan_rep is None:
        nan_rep = 'nan'

    data = lib.string_array_replace_from_nan_rep(data, nan_rep)
    return data.reshape(shape)

def _maybe_convert(values, val_kind, encoding):
    if _need_convert(val_kind):
        conv = _get_converter(val_kind, encoding)
        # conv = np.frompyfunc(conv, 1, 1)
        values = conv(values)
    return values

def _get_converter(kind, encoding):
    kind = _ensure_decoded(kind)
    if kind == 'datetime64':
        return lambda x: np.array(x, dtype='M8[ns]')
    elif kind == 'datetime':
        return lib.convert_timestamps
    elif kind == 'string':
        return lambda x: _unconvert_string_array(x,encoding=encoding)
    else:  # pragma: no cover
        raise ValueError('invalid kind %s' % kind)

def _need_convert(kind):
    kind = _ensure_decoded(kind)
    if kind in (u'datetime', u'datetime64', u'string'):
        return True
    return False

class Term(StringMixin):
    """create a term object that holds a field, op, and value

    Parameters
    ----------
    field : dict, string term expression, or the field to operate (must be a valid index/column type of DataFrame/Panel)
    op    : a valid op (defaults to '=') (optional)
            >, >=, <, <=, =, != (not equal) are allowed
    value : a value or list of values (required)
    queryables : a kinds map (dict of column name -> kind), or None i column is non-indexable
    encoding : an encoding that will encode the query terms

    Returns
    -------
    a Term object

    Examples
    --------
    >>> Term(dict(field = 'index', op = '>', value = '20121114'))
    >>> Term('index', '20121114')
    >>> Term('index', '>', '20121114')
    >>> Term('index', ['20121114','20121114'])
    >>> Term('index', datetime(2012,11,14))
    >>> Term('major_axis>20121114')
    >>> Term('minor_axis', ['A','U'])
    """

    _ops = ['<=', '<', '>=', '>', '!=', '==', '=']
    _search = re.compile("^\s*(?P<field>\w+)\s*(?P<op>%s)\s*(?P<value>.+)\s*$" % '|'.join(_ops))
    _max_selectors = 31

    def __init__(self, field, op=None, value=None, queryables=None, encoding=None):
        self.field = None
        self.op = None
        self.value = None
        self.q = queryables or dict()
        self.filter = None
        self.condition = None
        self.encoding = encoding

        # unpack lists/tuples in field
        while(isinstance(field, (tuple, list))):
            f = field
            field = f[0]
            if len(f) > 1:
                op = f[1]
            if len(f) > 2:
                value = f[2]

        # backwards compatible
        if isinstance(field, dict):
            self.field = field.get('field')
            self.op = field.get('op') or '=='
            self.value = field.get('value')

        # passed a term
        elif isinstance(field, Term):
            self.field = field.field
            self.op = field.op
            self.value = field.value

        # a string expression (or just the field)
        elif isinstance(field, basestring):

            # is a term is passed
            s = self._search.match(field)
            if s is not None:
                self.field = s.group('field')
                self.op = s.group('op')
                self.value = s.group('value')

            else:
                self.field = field

                # is an op passed?
                if isinstance(op, basestring) and op in self._ops:
                    self.op = op
                    self.value = value
                else:
                    self.op = '=='
                    self.value = op

        else:
            raise ValueError(
                "Term does not understand the supplied field [%s]" % field)

        # we have valid fields
        if self.field is None or self.op is None or self.value is None:
            raise ValueError("Could not create this term [%s]" % str(self))

         # = vs ==
        if self.op == '=':
            self.op = '=='

        # we have valid conditions
        if self.op in ['>', '>=', '<', '<=']:
            if hasattr(self.value, '__iter__') and len(self.value) > 1 and not isinstance(self.value,basestring):
                raise ValueError("an inequality condition cannot have multiple values [%s]" % str(self))

        if not is_list_like(self.value):
            self.value = [self.value]

        if len(self.q):
            self.eval()

    def __unicode__(self):
        attrs = map(pprint_thing, (self.field, self.op, self.value))
        return "field->%s,op->%s,value->%s" % tuple(attrs)

    @property
    def is_valid(self):
        """ return True if this is a valid field """
        return self.field in self.q

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an actual column in the table) """
        return self.q.get(self.field) is not None

    @property
    def kind(self):
        """ the kind of my field """
        return self.q.get(self.field)

    def generate(self, v):
        """ create and return the op string for this TermValue """
        val = v.tostring(self.encoding)
        return "(%s %s %s)" % (self.field, self.op, val)

    def eval(self):
        """ set the numexpr expression for this term """

        if not self.is_valid:
            raise ValueError("query term is not valid [%s]" % str(self))

        # convert values if we are in the table
        if self.is_in_table:
            values = [self.convert_value(v) for v in self.value]
        else:
            values = [TermValue(v,v,self.kind) for v in self.value]

        # equality conditions
        if self.op in ['==', '!=']:

            # our filter op expression
            if self.op == '!=':
                filter_op = lambda axis, vals: not axis.isin(vals)
            else:
                filter_op = lambda axis, vals: axis.isin(vals)


            if self.is_in_table:

                # too many values to create the expression?
                if len(values) <= self._max_selectors:
                    vs = [ self.generate(v) for v in values ]
                    self.condition = "(%s)" % ' | '.join(vs)

                # use a filter after reading
                else:
                    self.filter = (self.field, filter_op, Index([v.value for v in values]))

            else:

                self.filter = (self.field, filter_op, Index([v.value for v in values]))

        else:

            if self.is_in_table:

                self.condition = self.generate(values[0])

            else:

                raise TypeError("passing a filterable condition to a non-table indexer [%s]" % str(self))

    def convert_value(self, v):
        """ convert the expression that is in the term to something that is accepted by pytables """

        def stringify(value):
            value = str(value)
            if self.encoding is not None:
                value = value.encode(self.encoding)
            return value

        kind = _ensure_decoded(self.kind)
        if kind == u'datetime64' or kind == u'datetime' :
            v = lib.Timestamp(v)
            if v.tz is not None:
                v = v.tz_convert('UTC')
            return TermValue(v,v.value,kind)
        elif isinstance(v, datetime) or hasattr(v, 'timetuple') or kind == u'date':
            v = time.mktime(v.timetuple())
            return TermValue(v,Timestamp(v),kind)
        elif kind == u'integer':
            v = int(float(v))
            return TermValue(v,v,kind)
        elif kind == u'float':
            v = float(v)
            return TermValue(v,v,kind)
        elif kind == u'bool':
            if isinstance(v, basestring):
                v = not v.strip().lower() in [u'false', u'f', u'no', u'n', u'none', u'0', u'[]', u'{}', u'']
            else:
                v = bool(v)
            return TermValue(v,v,kind)
        elif not isinstance(v, basestring):
            v = stringify(v)
            return TermValue(v,stringify(v),u'string')

        # string quoting
        return TermValue(v,stringify(v),u'string')

class TermValue(object):
    """ hold a term value the we use to construct a condition/filter """

    def __init__(self, value, converted, kind):
        self.value = value
        self.converted = converted
        self.kind = kind

    def tostring(self, encoding):
        """ quote the string if not encoded
            else encode and return """
        if self.kind == u'string':
            if encoding is not None:
                return self.converted
            return '"%s"' % self.converted
        return self.converted

class Coordinates(object):
    """ holds a returned coordinates list, useful to select the same rows from different tables

    coordinates : holds the array of coordinates
    group       : the source group
    where       : the source where
    """

    def __init__(self, values, group, where, **kwargs):
        self.values = values
        self.group = group
        self.where = where

    def __len__(self):
        return len(self.values)

    def __getitem__(self, key):
        """ return a new coordinates object, sliced by the key """
        return Coordinates(self.values[key], self.group, self.where)

class Selection(object):
    """
    Carries out a selection operation on a tables.Table object.

    Parameters
    ----------
    table : a Table object
    where : list of Terms (or convertable to)
    start, stop: indicies to start and/or stop selection

    """
    def __init__(self, table, where=None, start=None, stop=None, **kwargs):
        self.table = table
        self.where = where
        self.start = start
        self.stop = stop
        self.condition = None
        self.filter = None
        self.terms = None
        self.coordinates = None

        if isinstance(where, Coordinates):
            self.coordinates = where.values
        else:
            self.terms = self.generate(where)

            # create the numexpr & the filter
            if self.terms:
                terms = [ t for t in self.terms if t.condition is not None ]
                if len(terms):
                    self.condition = "(%s)" % ' & '.join([ t.condition for t in terms ])
                self.filter = []
                for t in self.terms:
                    if t.filter is not None:
                        self.filter.append(t.filter)

    def generate(self, where):
        """ where can be a : dict,list,tuple,string """
        if where is None:
            return None

        if not isinstance(where, (list, tuple)):
            where = [where]
        else:

            # make this a list of we think that we only have a sigle term & no
            # operands inside any terms
            if not any([isinstance(w, (list, tuple, Term)) for w in where]):

                if not any([isinstance(w, basestring) and Term._search.match(w) for w in where]):
                    where = [where]

        queryables = self.table.queryables()
        return [Term(c, queryables=queryables, encoding=self.table.encoding) for c in where]

    def select(self):
        """
        generate the selection
        """
        if self.condition is not None:
            return self.table.table.readWhere(self.condition, start=self.start, stop=self.stop)
        elif self.coordinates is not None:
            return self.table.table.readCoordinates(self.coordinates)
        return self.table.table.read(start=self.start, stop=self.stop)

    def select_coords(self):
        """
        generate the selection
        """
        if self.condition is None:
            return np.arange(self.table.nrows)

        return self.table.table.getWhereList(self.condition, start=self.start, stop=self.stop, sort=True)


### utilities ###

def timeit(key,df,fn=None,remove=True,**kwargs):
    if fn is None:
        fn = 'timeit.h5'
    store = HDFStore(fn,mode='w')
    store.append(key,df,**kwargs)
    store.close()

    if remove:
        import os
        os.remove(fn)
