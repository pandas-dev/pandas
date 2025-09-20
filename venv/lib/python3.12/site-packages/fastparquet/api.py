"""parquet - read parquet files."""
import ast
from collections import OrderedDict, defaultdict
import re
import struct

import numpy as np
import fsspec
import pandas as pd

from fastparquet import core, schema, converted_types, encoding, dataframe, writer
from fastparquet import parquet_thrift
from fastparquet.cencoding import ThriftObject, from_buffer
from fastparquet.json import json_decoder
from fastparquet.util import (default_open, default_remove, ParquetException, val_to_num,
                   ops, ensure_bytes, ensure_str, check_column_names, metadata_from_many,
                   ex_from_sep, _strip_path_tail, get_fs, PANDAS_VERSION, join_path)


# Find in names of partition files the integer matching "**part.*.parquet",
# as 'i'.
PART_ID = re.compile(r'.*part.(?P<i>[\d]+).parquet$')


class ParquetFile(object):
    """The metadata of a parquet file or collection

    Reads the metadata (row-groups and schema definition) and provides
    methods to extract the data from the files.

    Note that when reading parquet files partitioned using directories
    (i.e. using the hive/drill scheme), an attempt is made to coerce
    the partition values to a number, datetime or timedelta. Fastparquet
    cannot read a hive/drill parquet file with partition names which coerce
    to the same value, such as "0.7" and ".7".

    Parameters
    ----------
    fn: path/URL string or list of paths
        Location of the data. If a directory, will attempt to read a file
        "_metadata" within that directory. If a list of paths, will assume
        that they make up a single parquet data set. This parameter can also
        be any file-like object, in which case this must be a single-file
        dataset.
    verify: bool [False]
        test file start/end byte markers
    open_with: function
        With the signature `func(path, mode)`, returns a context which
        evaluated to a file open for reading. Defaults to the built-in `open`.
    root: str
        If passing a list of files, the top directory of the data-set may
        be ambiguous for partitioning where the upmost field has only one
        value. Use this to specify the dataset root directory, if required.
    fs: fsspec-compatible filesystem
        You can use this instead of open_with (otherwise, it will be inferred)
    pandas_nulls: bool (True)
        If True, columns that are int or bool in parquet, but have nulls, will become
        pandas nullale types (Uint, Int, boolean). If False (the only behaviour
        prior to v0.7.0), both kinds will be cast to float, and nulls will be NaN
        unless pandas metadata indicates that the original datatypes were nullable.
        Pandas nullable types were introduces in v1.0.0, but were still marked as
        experimental in v1.3.0.

    Attributes
    ----------
    fn: path/URL
        Of '_metadata' file.
    basepath: path/URL
        Of directory containing files of parquet dataset.
    cats: dict
        Columns derived from hive/drill directory information, with known
        values for each column.
    categories: list
        Columns marked as categorical in the extra metadata (meaning the
        data must have come from pandas).
    columns: list of str
        The data columns available
    count: int
        Total number of rows
    dtypes: dict
        Expected output types for each column
    file_scheme: str
        'simple': all row groups are within the same file; 'hive': all row
        groups are in other files; 'mixed': row groups in this file and others
        too; 'empty': no row groups at all.
    info: dict
        Combination of some of the other attributes
    key_value_metadata: dict
        Additional information about this data's origin, e.g., pandas
        description, and custom metadata defined by user.
    row_groups: list
        Thrift objects for each row group
    schema: schema.SchemaHelper
        print this for a representation of the column structure
    selfmade: bool
        If this file was created by fastparquet
    statistics: dict
        Max/min/count of each column chunk
    fs: fsspec-compatible filesystem
        You can use this instead of open_with (otherwise, it will be inferred)
    """
    _pdm = None
    _kvm = None
    _categories = None

    def __init__(self, fn, verify=False, open_with=default_open, root=False,
                 sep=None, fs=None, pandas_nulls=True, dtypes=None):
        self.pandas_nulls = pandas_nulls
        self._base_dtype = dtypes
        self.tz = None
        self._columns_dtype = None
        if open_with is default_open and fs is None:
            fs = fsspec.filesystem("file")
        elif fs is not None:
            open_with = fs.open
        else:
            fs = getattr(open_with, "__self__", None)
        if fs is None:
            fs, fn, open_with, mkdirs = get_fs(fn, open_with, None)

        if isinstance(fn, (tuple, list)):
            if root and fs is not None:
                root = fs._strip_protocol(root)
            basepath, fmd = metadata_from_many(fn, verify_schema=verify,
                                               open_with=open_with, root=root,
                                               fs=fs)
            writer.consolidate_categories(fmd)
            self.fn = join_path(
                basepath, '_metadata') if basepath else '_metadata'
            self.fmd = fmd
            self._set_attrs()
        elif hasattr(fn, 'read'):
            # file-like
            self.fn = None
            self._parse_header(fn, verify)
            if self.file_scheme not in ['simple', 'empty']:
                raise ValueError('Cannot use file-like input '
                                 'with multi-file data')
            open_with = lambda *args, **kwargs: fn
        elif isinstance(fs, fsspec.AbstractFileSystem):
            if fs.isfile(fn):
                self.fn = join_path(fn)
                with open_with(fn, 'rb') as f:
                    self._parse_header(f, verify)
                if root:
                    paths = [fn.replace(root, "")]
                    self.file_scheme, self.cats = paths_to_cats(paths, None)
            elif "*" in fn or fs.isdir(fn):
                fn2 = join_path(fn, '_metadata')
                if fs.exists(fn2):
                    self.fn = fn2
                    with open_with(fn2, 'rb') as f:
                        self._parse_header(f, verify)
                    fn = fn2
                else:
                    # TODO: get details from fs here, rather than do suffix cat in
                    #  metadata_from_many
                    if "*" in fn:
                        allfiles = fs.glob(fn)
                    else:
                        allfiles = [f for f in fs.find(fn) if
                                    f.endswith(".parquet") or f.endswith(".parq")]
                        root = root or fn
                    if not allfiles:
                        raise ValueError("No files in dir")
                    if root:
                        root = fs._strip_protocol(root)
                    basepath, fmd = metadata_from_many(allfiles, verify_schema=verify,
                                                       open_with=open_with, root=root,
                                                       fs=fs)
                    writer.consolidate_categories(fmd)
                    self.fn = join_path(basepath, '_metadata') if basepath \
                              else '_metadata'
                    self.fmd = fmd
                    self._set_attrs()
                self.fs = fs
            else:
                raise FileNotFoundError(fn)
        else:
            done = False
            try:
                self.fn = fn
                f = open_with(fn)
                self._parse_header(f, verify)
                done = True
            except IOError:
                pass
            if not done:
                # allow this to error with FileNotFound or whatever
                try:
                    self.fn = join_path(fn, "_metadata")
                    f = open_with(self.fn)
                    self._parse_header(f, verify)
                except IOError as e:
                    raise ValueError("Opening directories without a _metadata requires"
                                     "a filesystem compatible with fsspec") from e
        self.open = open_with
        self._statistics = None

    def _parse_header(self, f, verify=True):
        if self.fn and self.fn.endswith("_metadata"):
            #  no point attempting to read footer only for pure metadata
            data = f.read()[4:-8]
            self._head_size = len(data)
        else:
            try:
                f.seek(0)
                if verify:
                    assert f.read(4) == b'PAR1'
                f.seek(-8, 2)
                head_size = struct.unpack('<I', f.read(4))[0]
                if verify:
                    assert f.read() == b'PAR1'
                self._head_size = head_size
                f.seek(-(head_size + 8), 2)
                data = f.read(head_size)
            except (AssertionError, struct.error):
                raise ParquetException('File parse failed: %s' % self.fn)

        try:
            fmd = from_buffer(data, "FileMetaData")
        except Exception:
            raise ParquetException('Metadata parse failed: %s' % self.fn)
        # for rg in fmd.row_groups:
        for rg in fmd[4]:
            # chunks = rg.columns
            chunks = rg[1]    
            if chunks:
                chunk = chunks[0]
                # s = chunk.file_path
                s = chunk.get(1)
                if s:
                    # chunk.file_path = s.decode()
                    chunk[1] = s.decode()
        self.fmd = fmd
        self._set_attrs()

    def _set_attrs(self):
        fmd = self.fmd
        self.version = fmd.version
        self._schema = fmd.schema
        self.row_groups = fmd.row_groups or []
        self.created_by = fmd.created_by
        self.schema = schema.SchemaHelper(self._schema)
        self.selfmade = (
            b"fastparquet" in self.created_by if self.created_by is not None
            else False
        )
        self._read_partitions()
        self._dtypes()

    @property
    def helper(self):
        return self.schema

    @property
    def columns(self):
        """ Column names """
        return [_ for _ in self.dtypes if _ not in self.cats]

    @property
    def statistics(self):
        if not hasattr(self, '_statistics') or self._statistics is None:
            self._statistics = statistics(self)
        return self._statistics

    @property
    def key_value_metadata(self):
        if self._kvm is None:
            self._kvm = {
                ensure_str(k.key, ignore_error=True): ensure_str(k.value, ignore_error=True)
                for k in self.fmd.key_value_metadata or []}
        return self._kvm

    @property
    def partition_meta(self):
        return {col['field_name']: col for col in self.pandas_metadata.get('partition_columns', [])}

    @property
    def basepath(self):
        return re.sub(r'_metadata(/)?$', '', self.fn).rstrip('/')

    def _read_partitions(self):
        # paths = [rg.columns[0].file_path] ... if rg.columns]
        paths = [rg[1][0].get(1, "") for rg in self.row_groups if rg[1]]
        self.file_scheme, self.cats = paths_to_cats(paths, self.partition_meta)

    def head(self, nrows, **kwargs):
        """Get the first nrows of data

        This will load the whole of the first valid row-group for the given
        columns.

        kwargs can include things like columns, filters, etc., with the same
        semantics as to_pandas(). If filters are applied, it may happen that
        data is so reduced that 'nrows' is not ensured (fewer rows). 

        returns: dataframe
        """
        # TODO: implement with truncated assign and early exit
        #  from reading
        total_rows = 0
        for i, rg in enumerate(self.row_groups):
            total_rows += rg.num_rows
            if total_rows >= nrows:
                break
        return self[:i+1].to_pandas(**kwargs).head(nrows)

    def __getitem__(self, item):
        """Select among the row-groups using integer/slicing"""
        import copy
        new_rgs = self.row_groups[item]
        if not isinstance(new_rgs, list):
            new_rgs = [new_rgs]
        new_pf = object.__new__(ParquetFile)
        fmd = copy.copy(self.fmd)
        fmd.row_groups = new_rgs
        new_pf.__setstate__(
            {"fn": self.fn, "open": self.open, "fmd": fmd,
             "pandas_nulls": self.pandas_nulls, "_base_dtype": self._base_dtype,
             "tz": self.tz, "_columns_dtype": self._columns_dtype}
        )
        new_pf._set_attrs()
        return new_pf

    def __len__(self):
        """Return number of row groups."""
        if self.fmd.row_groups:
            return len(self.fmd.row_groups)
        else:
            return 0

    def __bool__(self):
        """Return True, takes precedence over `__len__`."""
        return True

    def row_group_filename(self, rg):
        if rg.columns and rg.columns[0].file_path:
            base = self.basepath
            if base:
                return join_path(base, rg.columns[0].file_path)
            else:
                return rg.columns[0].file_path
        else:
            return self.fn

    def read_row_group_file(self, rg, columns, categories, index=None,
                            assign=None, partition_meta=None, row_filter=False,
                            infile=None):
        """ Open file for reading, and process it as a row-group

        assign is None if this method is called directly (not from to_pandas),
        in which case we return the resultant dataframe

        row_filter can be:

           -  False (don't do row filtering)

           -  a list of filters (do filtering here for this one row-group;
              only makes sense if assign=None

           -  bool array with a size equal to the number of rows in this group
              and the length of the assign arrays
        """
        categories = self.check_categories(categories)
        fn = self.row_group_filename(rg)
        ret = False
        if assign is None:
            if row_filter and isinstance(row_filter, list):
                cs = self._columns_from_filters(row_filter)
                df = self.read_row_group_file(
                    rg, cs, categories, index=False,
                    infile=infile, row_filter=False)
                row_filter = self._column_filter(df, filters=row_filter)
                size = row_filter.sum()
                if size == rg.num_rows:
                    row_filter = False
            else:
                size = rg.num_rows
            df, assign = self.pre_allocate(
                    size, columns, categories, index)
            if "PANDAS_ATTRS" in self.key_value_metadata:
                import json
                df.attrs = json.loads(self.key_value_metadata["PANDAS_ATTRS"])
            ret = True
        f = infile or self.open(fn, mode='rb')

        core.read_row_group(
            f, rg, columns, categories, self.schema, self.cats,
            selfmade=self.selfmade, index=index,
            assign=assign, scheme=self.file_scheme, partition_meta=partition_meta,
            row_filter=row_filter
        )
        if ret:
            return df

    def iter_row_groups(self, filters=None, **kwargs):
        """
        Iterate a dataset by row-groups

        If filters is given, omits row-groups that fail the filer
        (saving execution time)

        Returns
        -------
        Generator yielding one Pandas data-frame per row-group.
        """
        rgs = filter_row_groups(self, filters) if filters else self.row_groups
        for rg in rgs:
            i = self.row_groups.index(rg)
            df = self[i].to_pandas(filters=filters, **kwargs)
            if not df.empty:
                yield df

    def remove_row_groups(self, rgs, sort_pnames:bool=False,
                          write_fmd:bool=True, open_with=default_open,
                          remove_with=None):
        """
        Remove list of row groups from disk. `ParquetFile` metadata are
        updated accordingly. This method can not be applied if file scheme is
        simple.

        Parameter
        ---------
        rgs: row group or list of row groups
            List of row groups to be removed from disk.
        sort_pnames : bool, default False
            Align name of part files with position of the 1st row group they
            contain. Only used if `file_scheme` of parquet file is set to
            `hive` or `drill`.
        write_fmd: bool, True
            Write updated common metadata to disk.
        open_with: function
            When called with f(path, mode), returns an open file-like object.
        remove_with: function
            When called with f(path) removes the file or directory given
            (and any contained files). Not required if this ParquetFile has
            a .fs file system attribute
        """
        if not isinstance(rgs, list):
            if isinstance(rgs, ThriftObject) or isinstance(rgs, dict):
                # Case 'rgs' is a single row group ('ThriftObject' or 'dict'). 
                rgs = [rgs]
            else:
                # Use `list()` here, not `[]`, as the latter does not transform
                # generator or tuple into list but encapsulates them in a list.
                rgs = list(rgs)
        if rgs:
            if self.file_scheme == 'simple':
                raise ValueError("Not possible to remove row groups when file "
                                 "scheme is 'simple'.")
            if remove_with is None:
                if hasattr(self, 'fs'):
                    remove_with = self.fs.rm
                else:
                    remove_with = default_remove
            rgs_to_remove = row_groups_map(rgs)
            if (b"fastparquet" not in self.created_by
                or self.file_scheme == 'flat'):
                # Check if some files contain row groups both to be removed and
                # to be kept.
                all_rgs = row_groups_map(self.fmd.row_groups)
                for file in rgs_to_remove:
                    if len(rgs_to_remove[file]) < len(all_rgs[file]):
                        raise ValueError(
                            f"File {file} contains row groups both to be kept "
                            "and to be removed. Removing row groups partially "
                            "from a file is not possible.")
            if rgs != self.fmd.row_groups:
                rg_new = self.fmd.row_groups
            else:
                # Deep copy required if 'rg_new' and 'rgs' points both to
                # 'self.fmd.row_groups'.
                from copy import deepcopy
                rg_new = deepcopy(self.fmd.row_groups)
            for rg in rgs:
                rg_new.remove(rg)
                self.fmd.num_rows -= rg.num_rows
            self.fmd.row_groups = rg_new
            try:
                basepath = self.basepath
                remove_with([f'{basepath}/{file}' for file in rgs_to_remove])
            except IOError:
                pass
            self._set_attrs()
        if sort_pnames:
            self._sort_part_names(False, open_with)
        if write_fmd:
            self._write_common_metadata(open_with)

    def write_row_groups(self, data, row_group_offsets=None, sort_key=None,
                         sort_pnames:bool=False, compression=None,
                         write_fmd:bool=True, open_with=default_open,
                         mkdirs=None, stats="auto"):
        """Write data as new row groups to disk, with optional sorting.

        Parameters
        ----------
        data : pandas dataframe or iterable of pandas dataframe
            Data to add to existing parquet dataset. Only columns are written
            to disk. Row index is not kept.
            If a dataframe, columns are checked against parquet file schema.
        row_group_offsets: int or list of int
            If int, row-groups will be approximately this many rows, rounded down
            to make row groups about the same size;
            If a list, the explicit index values to start new row groups;
            If `None`, set to 50_000_000.
        sort_key : function, default None
            Sorting function used as `key` parameter for `row_groups.sort()`
            function. This function is called once new row groups have been
            added to list of existing ones.
            If not provided, new row groups are only appended to existing ones
            and the updated list of row groups is not sorted.
        sort_pnames : bool, default False
            Align name of part files with position of the 1st row group they
            contain. Only used if `file_scheme` of parquet file is set to
            `hive` or `drill`.
        compression : str or dict, default None
            Compression to apply to each column, e.g. ``GZIP`` or ``SNAPPY`` or
            a ``dict`` like ``{"col1": "SNAPPY", "col2": None}`` to specify per
            column compression types.
            By default, do not compress.
            Please, review full description of this parameter in `write`
            docstring.
        write_fmd : bool, default True
            Write updated common metadata to disk.
        open_with : function
            When called with a f(path, mode), returns an open file-like object.
        mkdirs : function
            When called with a path/URL, creates any necessary dictionaries to
            make that location writable, e.g., ``os.makedirs``. This is not
            necessary if using the simple file scheme.
        stats : True|False|list of str
            Whether to calculate and write summary statistics.
            If True (default), do it for every column;
            If False, never do;
            If a list of str, do it only for those specified columns.
            "auto" means True for any int/float or timemstamp column, False
            otherwise. This will become the default in a future release.
        """
        from .writer import write_simple, write_multi
        partition_on = list(self.cats)
        if isinstance(data, pd.DataFrame):
            self_cols = sorted(self.columns + partition_on)
            if self_cols != sorted(data.columns):
                diff_cols = set(data.columns) ^ set(self_cols)
                raise ValueError(
                    f'Column names of new data are {sorted(data.columns)}. '
                    f'But column names in existing file are {self_cols}. '
                    f'{diff_cols} are columns being either only in existing '
                     'file or only in new data. This is not possible.')
        if (self.file_scheme == 'simple'
            or (self.file_scheme == 'empty' and self.fn[-9:] != '_metadata')):
            # Case 'simple'.
            write_simple(self.fn, data, self.fmd,
                         row_group_offsets=row_group_offsets,
                         compression=compression, open_with=open_with,
                         has_nulls=None, append=True, stats=stats)
        else:
            # Case 'hive' or 'drill'.
            write_multi(self.basepath, data, self.fmd,
                        row_group_offsets=row_group_offsets,
                        compression=compression, file_scheme=self.file_scheme,
                        write_fmd=False, open_with=open_with, mkdirs=mkdirs,
                        partition_on=partition_on, append=True, stats=stats)
            if sort_key:
                # Not using 'sort()' because 'row_groups' is a ThriftObject,
                # not a list.
                self.fmd.row_groups = sorted(self.fmd.row_groups, key=sort_key)
            if sort_pnames:
                self._sort_part_names(False, open_with)
            if write_fmd:
                self._write_common_metadata(open_with)
        self._set_attrs()

    def _sort_part_names(self, write_fmd:bool=True, open_with=default_open):
        """Align parquet files id to that of the first row group they contain.

        This method only manages files which name follows pattern
        "part.{id}.parquet". Field `id` is then aligned to the index of the
        first row group it contains. The index of a row groups is its position
        in row group list.

        Parameters
        ----------
        write_fmd : bool, default True
            Write updated common metadata to disk.
        open_with : function
            When called with a f(path, mode), returns an open file-like object.
            Only needed if `write_fmd` is `True`.
        """
        from .writer import part_ids
        pids = part_ids(self.fmd.row_groups)
        if pids:
            # Keep only items for which row group position does not match part
            # name id.
            pids = dict(filter(lambda item: item[0] != item[1][0],
                               pids.items()))
            basepath = self.basepath
            # Give temporary names in a 1st pass to prevent overwritings.
            for pid in pids:
                item = pids[pid]
                rgid, fname = item[0], item[1]
                src = f'{basepath}/{fname}'
                parts = partitions(fname)
                dst = join_path(basepath, parts, f'part.{rgid}.parquet.tmp')
                self.fs.rename(src, dst)
            # Give definitive names in a 2nd pass.
            for pid in pids:
                item = pids[pid]
                rgid, fname = item[0], item[1]
                parts = partitions(fname)
                src = join_path(basepath, parts, f'part.{rgid}.parquet.tmp')
                dst_part = join_path(parts, f'part.{rgid}.parquet')
                dst = join_path(basepath, dst_part)
                self.fs.rename(src, dst)
                for col in self.fmd.row_groups[rgid].columns:
                    col.file_path = dst_part
            if write_fmd:
                self._write_common_metadata(open_with)

    def _write_common_metadata(self, open_with=default_open):
        """
        Write common metadata to disk.
        
        Parameter
        ---------
        open_with: function
            When called with a f(path, mode), returns an open file-like object.
        """
        from .writer import write_common_metadata
        if self.file_scheme == 'simple':
            raise ValueError("Not possible to write common metadata when file \
scheme is 'simple'.")
        fmd = self.fmd
        write_common_metadata(self.fn, fmd, open_with, no_row_groups=False)
        # replace '_metadata' with '_common_metadata'
        fn = f'{self.fn[:-9]}_common_metadata'
        write_common_metadata(fn, fmd, open_with)

    def _get_index(self, index=None):
        if index is None:
            index = [i if isinstance(i, str) else i["name"]
                     for i in self.pandas_metadata.get('index_columns', [])
                     if isinstance(i, str) or i.get("kind") != "range"
                     ]
        if isinstance(index, str):
            index = [index]
        return index

    def _columns_from_filters(self, filters):
        return [
            c for c in
            set(sum([[f[0]]
                     if isinstance(f[0], str)
                     else [g[0] for g in f] for f in filters], []))
            if c not in self.cats
        ]

    def _column_filter(self, df, filters):
        out = np.zeros(len(df), dtype=bool)
        for or_part in filters:
            if isinstance(or_part[0], str):
                name, op, val = or_part
                if name in self.cats:
                    continue
                if op == 'in':
                    out |= df[name].isin(val).values
                elif op == "not in":
                    out |= ~df[name].isin(val).values
                elif op in ops:
                    out |= ops[op](df[name], val).values
                elif op == "~":
                    out |= ~df[name].values
            else:
                and_part = np.ones(len(df), dtype=bool)
                for name, op, val in or_part:
                    if name in self.cats:
                        continue
                    if op == 'in':
                        and_part &= df[name].isin(val).values
                    elif op == "not in":
                        and_part &= ~df[name].isin(val).values
                    elif op in ops:
                        and_part &= ops[op](df[name].values, val)
                    elif op == "~":
                        and_part &= ~df[name].values
                out |= and_part
        return out

    def to_pandas(self, columns=None, categories=None, filters=[],
                  index=None, row_filter=False, dtypes=None):
        """
        Read data from parquet into a Pandas dataframe.

        Parameters
        ----------
        columns: list of names or `None`
            Column to load (see `ParquetFile.columns`). Any columns in the
            data not in this list will be ignored. If `None`, read all columns.
        categories: list, dict or `None`
            If a column is encoded using dictionary encoding in every row-group
            and its name is also in this list, it will generate a Pandas
            Category-type column, potentially saving memory and time. If a
            dict {col: int}, the value indicates the number of categories,
            so that the optimal data-dtype can be allocated. If ``None``,
            will automatically set *if* the data was written from pandas.
        filters: list of list of tuples or list of tuples
            To filter out data.
            Filter syntax: [[(column, op, val), ...],...]
            where op is [==, =, >, >=, <, <=, !=, in, not in]
            The innermost tuples are transposed into a set of filters applied
            through an `AND` operation.
            The outer list combines these sets of filters through an `OR`
            operation.
            A single list of tuples can also be used, meaning that no `OR`
            operation between set of filters is to be conducted.
        index: string or list of strings or False or None
            Column(s) to assign to the (multi-)index. If None, index is
            inferred from the metadata (if this was originally pandas data); if
            the metadata does not exist or index is False, index is simple
            sequential integers.
        row_filter: bool or boolean ndarray
            Whether filters are applied to whole row-groups (False, default)
            or row-wise (True, experimental). The latter requires two passes of
            any row group that may contain valid rows, but can be much more
            memory-efficient, especially if the filter columns are not required
            in the output.
            If boolean array, it is applied as custom row filter. In this case,
            'filter' parameter is ignored, and length of the array has to be
            equal to the total number of rows.

        Returns
        -------
        Pandas data-frame
        """
        rgs = filter_row_groups(self, filters) if filters else self.row_groups
        index = self._get_index(index)
        if columns is not None:
            columns = columns[:]
        else:
            columns = self.columns + list(self.cats)
        if index:
            columns += [i for i in index if i not in columns]
        check_column_names(self.columns + list(self.cats), columns, categories)
        if row_filter is not False:
            if filters and row_filter is True:
                # Rows are selected as per filters.
                # TODO: special case when filter columns are also in output
                cs = self._columns_from_filters(filters)
                df = self.to_pandas(columns=cs, filters=filters, row_filter=False,
                                    index=False)
                sel = self._column_filter(df, filters=filters)
            else:
                # Row are selected as per custom 'sel'.
                if sum(rg.num_rows for rg in rgs) != len(row_filter):
                    raise ValueError('Provided boolean array for custom row \
selection does not match number of rows in DataFrame.')
                sel = row_filter
            size = sel.sum()
            selected = []
            start = 0
            for rg in rgs[:]:
                selected.append(sel[start:start+rg.num_rows])
                start += rg.num_rows
        else:
            size = sum(rg.num_rows for rg in rgs)
            selected = [None] * len(rgs)  # just to fill zip, below
        df, views = self.pre_allocate(size, columns, categories, index, dtypes=dtypes)
        if "PANDAS_ATTRS" in self.key_value_metadata:
            import json
            df.attrs = json.loads(self.key_value_metadata["PANDAS_ATTRS"])

        start = 0
        if self.file_scheme == 'simple':
            infile = self.open(self.fn, 'rb')
        else:
            infile = None
        for rg, sel in zip(rgs, selected):
            thislen = sel.sum() if sel is not None else rg.num_rows
            if thislen == rg.num_rows:
                # all good; noop if no row filtering
                sel = None
            elif thislen == 0:
                # no valid rows
                continue
            parts = {name: (v if name.endswith('-catdef')
                            else v[start:start + thislen])
                     for (name, v) in views.items()}
            self.read_row_group_file(rg, columns, categories, index,
                                     assign=parts, partition_meta=self.partition_meta,
                                     row_filter=sel, infile=infile)
            start += thislen
        return df

    def pre_allocate(self, size, columns, categories, index, dtypes=None):
        if dtypes is not None:
            columns = list(dtypes)
        else:
            dtypes = self._dtypes(categories)
        categories = self.check_categories(categories)
        cats = {k: v for k, v in self.cats.items() if k in columns}
        df, arrs = _pre_allocate(size, columns, categories, index, cats,
                                 dtypes, self.tz, columns_dtype=self._columns_dtype)
        i_no_name = re.compile(r"__index_level_\d+__")
        if self.has_pandas_metadata:
            md = self.pandas_metadata
            if categories:
                for c in md['columns']:
                    if c['name'] in categories and c['name'] in df and c['metadata']:
                        df[c['name']].dtype._ordered = c['metadata']['ordered']
            if md.get('index_columns', False) and not (index or index is False):
                if len(md['index_columns']) == 1:
                    ic = md['index_columns'][0]
                    if isinstance(ic, dict) and ic.get('kind') == 'range':
                        from pandas import RangeIndex
                        df.index = RangeIndex(
                            start=ic['start'],
                            stop=ic['start'] + size * ic['step'] + 1,
                            step=ic['step']
                        )[:size]
                names = [(c['name'] if isinstance(c, dict) else c)
                         for c in md['index_columns']]
                names = [None if n is None or i_no_name.match(n) else n
                         for n in names]
                df.index.names = names
            if md.get('column_indexes', False):
                names = [(c['name'] if isinstance(c, dict) else c)
                         for c in md['column_indexes']]
                names = [None if n is None or i_no_name.match(n) else n
                         for n in names]
                if len(names) > 1:
                    df.columns = pd.MultiIndex.from_tuples(
                        [ast.literal_eval(c) for c in df.columns if c not in df.index.names],
                        names=names
                    )
                else:
                    df.columns.names = names
        return df, arrs

    def count(self, filters=None, row_filter=False):
        """ Total number of rows

        filters and row_filters have the same meaning as in to_pandas. Unless both are given,
        this method will not need to decode any data
        """
        if row_filter:
            cs = self._columns_from_filters(filters)
            df = self.to_pandas(columns=cs, filters=filters, row_filter=False,
                                index=False)
            return self._column_filter(df, filters=filters).sum()

        rgs = filter_row_groups(self, filters)
        return sum(rg.num_rows for rg in rgs)

    @property
    def info(self):
        """ Dataset summary """
        return {'name': self.fn, 'columns': self.columns,
                'partitions': list(self.cats), 'rows': self.count(),
                "row_groups": len(self.row_groups)}

    def check_categories(self, cats):
        categ = self.categories
        if not self.has_pandas_metadata:
            return cats or {}
        if cats is None:
            return categ or {}
        if set(cats) - set(categ) and len(self.row_groups) > 1:
            raise TypeError("Attempt to read as category a field that "
                            "was not stored as such")
        if isinstance(cats, dict):
            return cats
        out = {k: v for k, v in categ.items() if k in cats}
        out.update({c: pd.RangeIndex(0, 2**14) for c in cats if c not in categ})
        return out

    @property
    def has_pandas_metadata(self):
        if self._pdm:
            return True
        if self.fmd.key_value_metadata is None:
            return False
        return bool(self.key_value_metadata.get('pandas', False))

    @property
    def pandas_metadata(self):
        if self._pdm is None:
            if self.has_pandas_metadata:
                self._pdm = json_decoder()(self.key_value_metadata['pandas'])
            else:
                self._pdm = {}
        return self._pdm

    @property
    def categories(self):
        if self._categories is not None:
            return self._categories
        if self.has_pandas_metadata:
            metadata = self.pandas_metadata
            if "column_indexes" in metadata and len(metadata["column_indexes"]) > 0:
                self._columns_dtype = metadata["column_indexes"][0]["numpy_type"]
            else:
                self._columns_dtype = None
            cats = {}
            for m in metadata['columns']:
                if m['pandas_type'] != 'categorical':
                    continue
                out = False
                if b"fastparquet" in self.created_by:
                    # if pandas was categorical, we will have used dict encoding
                    cats[m['name']] = m['metadata']['num_categories']
                    continue
                for rg in self.row_groups:
                    # but others (pyarrow) may have used dict for only some pages
                    if out:
                        break
                    for col in rg.columns:
                        if ".".join(col.meta_data.path_in_schema) != m['name']:
                            continue
                        if col.meta_data.encoding_stats:
                            if any(s.encoding not in [parquet_thrift.Encoding.PLAIN_DICTIONARY,
                                                  parquet_thrift.Encoding.RLE_DICTIONARY]
                                   for s in col.meta_data.encoding_stats
                                   if s.page_type in [parquet_thrift.PageType.DATA_PAGE_V2,
                                                      parquet_thrift.PageType.DATA_PAGE]):
                                out = True
                                break
                if out is False:
                    cats[m['name']] = m['metadata']['num_categories']
            self._categories = cats
            return cats
        # old track
        vals = self.key_value_metadata.get("fastparquet.cats", None)
        if vals:
            self._categories = json_decoder()(vals)
            return self._categories
        else:
            return {}

    def _dtypes(self, categories=None):
        """ Implied types of the columns in the schema """
        import pandas as pd
        if self._base_dtype is None:
            if self.has_pandas_metadata:
                md = self.pandas_metadata['columns']
                md = {c['name']: c for c in md}
                tz = {k: v["metadata"]['timezone'] for k, v in md.items()
                      if v.get('metadata', {}) and v.get('metadata', {}).get('timezone', None)}
            else:
                tz = None
                md = None
            self.tz = tz

            dtype = OrderedDict((name, (converted_types.typemap(f, md=md)
                                if f.num_children in [None, 0] else np.dtype("O")))
                                for name, f in self.schema.root["children"].items()
                                if getattr(f, 'isflat', False) is False)
            for i, (col, dt) in enumerate(dtype.copy().items()):
                # int and bool columns produce masked pandas types, no need to
                # promote types here
                if dt.kind == "M":
                    if self.pandas_metadata and PANDAS_VERSION.major >= 2:
                        # get original resolution when pandas supports non-ns
                        dt = md[col]["numpy_type"]
                    if tz is not None and tz.get(col, False):
                        z = dataframe.tz_to_dt_tz(tz[col])
                        dt_series = pd.Series([], dtype=dt)
                        if PANDAS_VERSION.major >= 2 and dt_series.dt.tz is not None:
                            dt = dt_series.dt.tz_convert(z).dtype
                        else:
                            dt = dt_series.dt.tz_localize(z).dtype
                    dtype[col] = dt
                elif dt in converted_types.nullable:
                    if self.pandas_metadata:
                        tt = md.get(col, {}).get("numpy_type")
                        if tt and ("int" in tt or "bool" in tt):
                            continue
                    # uint/int/bool columns that may have nulls become nullable
                    # skip is pandas_metadata gives original types
                    num_nulls = 0
                    for rg in self.row_groups:
                        if rg[3] == 0:
                            continue
                        st = rg[1][i][3].get(12)
                        if st is None:
                            num_nulls = True
                            break
                        if st.get(3):
                            num_nulls = True
                            break
                    if num_nulls:
                        if self.pandas_nulls:
                            dtype[col] = converted_types.nullable[dt]
                        else:
                            dtype[col] = np.float64()
                elif dt == 'S12':
                    dtype[col] = 'M8[ns]'
            self._base_dtype = dtype
        dtype = self._base_dtype.copy()
        categories = self.check_categories(categories)
        for field in categories:
            dtype[field] = 'category'
        for cat in self.cats:
            dtype[cat] = "category"
        self.dtypes = dtype
        return dtype

    def __getstate__(self):
        if self.fmd.row_groups is None:
            self.fmd.row_groups = []
        return {"fn": self.fn, "open": self.open, "fmd": self.fmd,
                "pandas_nulls": self.pandas_nulls, "_base_dtype": self._base_dtype,
                "tz": self.tz}

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Decode 'file_path'.
        rgs = self.fmd[4] or []
        # 4th condition should not be necessary, depends on 'deepcopy' version.
        # https://github.com/dask/fastparquet/pull/731#issuecomment-1013507287
        if (rgs and rgs[0][1] and rgs[0][1][0] and rgs[0][1][0].get(1)
                and isinstance(rgs[0][1][0].get(1), bytes)):
            # for rg in fmd.row_groups:
            for rg in rgs:
                # chunk = rg.columns[0]
                chunk = rg[1][0]
                # chunk.file_path = chunk.file_path.decode()
                chunk[1] = chunk.get(1).decode()
        self._set_attrs()

    def __str__(self):
        return "<Parquet File: %s>" % self.info

    __repr__ = __str__


def _pre_allocate(size, columns, categories, index, cs, dt, tz=None, columns_dtype=None):
    index = [index] if isinstance(index, str) else (index or [])
    cols = [c for c in columns if c not in index]
    categories = categories or {}
    cats = cs.copy()
    if isinstance(categories, dict):
        cats.update(categories)

    def get_type(name, index=False):
        if name in categories:
            return 'category'
        t = dt[name]
        if index and isinstance(t, pd.core.arrays.masked.BaseMaskedDtype):
            return "int64"
        return t

    dtypes = [get_type(c) for c in cols]
    index_types = [get_type(i, index=True) for i in index]
    cols.extend(cs)
    dtypes.extend(['category'] * len(cs))
    df, views = dataframe.empty(dtypes, size, cols=cols, index_names=index,
                                index_types=index_types, cats=cats, timezones=tz,
                                columns_dtype=columns_dtype)
    return df, views


def paths_to_cats(paths, partition_meta=None):
    """
    Extract categorical fields and labels from hive- or drill-style paths.

    Parameters
    ----------
    paths (Iterable[str]): file paths relative to root
    file_scheme (str):
    partition_meta (Dict[str, dict]):

    Returns
    -------
    cats (OrderedDict[str, List[Any]]): a dict of field names and their values
    """
    if len(paths) == 0:
        return "empty", {}

    if all(p in [None, ""] for p in paths):
        return "simple", {}
    paths = _strip_path_tail(paths)
    parts = [path.split("/") for path in paths if path]
    lparts = [len(part) for part in parts]
    if not lparts or max(lparts) < 1:
        return "flat", {}
    if len(set(lparts)) > 1:
        return "other", {}

    try:
        return "hive", _path_to_cats(paths, parts, partition_meta=partition_meta)
    except ValueError:
        return "drill", _path_to_cats(paths, parts, "drill", partition_meta=partition_meta)


def _path_to_cats(paths, parts, file_scheme="hive", partition_meta=None):
    partition_meta = partition_meta or {}
    cats = OrderedDict()
    s = ex_from_sep('/')
    string_types = set()
    meta = {"pandas_type": "string", "numpy_type": "object"}
    seen = set()
    for path, path_parts in zip(paths, parts):

        if file_scheme == "hive":
            hivehits = [p.split("=") for p in path.split("/") if "=" in p]  # s.findall(path)
            if not hivehits:
                raise ValueError("Not a hive scheme")
        if file_scheme == "drill":
            hivehits = [(f"dir{i}", v) for i, v in enumerate(path_parts)]
        for key, val in hivehits:

            if (key, val) in seen:
                continue
            seen.add((key, val))
            tp = val_to_num(val, meta if key in string_types else partition_meta.get(key))
            if isinstance(tp, str):
                string_types.add(key)
            cats.setdefault(key, set()).add(tp)

    return OrderedDict([(key, list(v)) for key, v in cats.items()])


def filter_out_stats(rg, filters, schema):
    """
    According to the filters, should this row-group be excluded

    Considers the statistics included in the metadata of this row-group

    Parameters
    ----------
    rg: thrift RowGroup structure
    filters: list of 3-tuples
        Structure of each tuple: (column, op, value) where op is one of
        ['==', '=', '!=', '<', '<=', '>', '>=', 'in', 'not in'] and value is
        appropriate for the column in question

    Returns
    -------
    True or False
    """
    if rg.num_rows == 0:
        # always ignore empty row-groups, don't bother loading
        return True
    if len(filters) == 0:
        return False
    for column in rg.columns:
        vmax, vmin = None, None
        name = ".".join(column.meta_data.path_in_schema)
        app_filters = [f[1:] for f in filters if f[0] == name]
        for op, val in app_filters:
            se = schema.schema_element(name)
            if column.meta_data.statistics is not None:
                s = column.meta_data.statistics
                if s.null_count == column.meta_data.num_values:
                    # skip row groups with no valid values
                    return True
                # we cache the converted valuesin the stats object
                # TODO: keep this somewhere not in a thrift object?
                max = s.max or s.max_value
                if max is not None:
                    if not hasattr(s, "converted_max"):
                        b = ensure_bytes(max)
                        vmax = encoding.read_plain(
                            b, column.meta_data.type, 1, stat=True)
                        if se.converted_type is not None or se.logicalType is not None:
                            vmax = converted_types.convert(vmax, se)
                        s["converted_max"] = vmax
                    vmax = s["converted_max"]
                min = s.min or s.min_value
                if min is not None:
                    if not hasattr(s, "converted_min"):
                        b = ensure_bytes(min)
                        vmin = encoding.read_plain(
                            b, column.meta_data.type, 1, stat=True)
                        if se.converted_type is not None or se.logicalType is not None:
                            vmin = converted_types.convert(vmin, se)
                        s["converted_min"] = vmin
                    vmin = s["converted_min"]
                if filter_val(op, val, vmin, vmax):
                    return True
    return False


def statistics(obj):
    """
    Return per-column statistics for a ParquetFile

    Parameters
    ----------
    obj: ParquetFile

    Returns
    -------
    dictionary mapping stats (min, max, distinct_count, null_count) to column
    names to lists of values.  ``None``s used if no statistics found.

    Examples
    --------
    >>> statistics(my_parquet_file)
    {'min': {'x': [1, 4], 'y': [5, 3]},
     'max': {'x': [2, 6], 'y': [8, 6]},
     'distinct_count': {'x': [None, None], 'y': [None, None]},
     'null_count': {'x': [0, 3], 'y': [0, 0]}}
    """
    if isinstance(obj, ThriftObject) and obj.thrift_name == "ColumnChunk":
        md = obj.meta_data
        s = obj.meta_data.statistics
        rv = {}
        if not s:
            return rv
        if s.max is not None:
            try:
                if md.type == parquet_thrift.Type.BYTE_ARRAY:
                    rv['max'] = ensure_bytes(s.max)
                else:
                    rv['max'] = encoding.read_plain(ensure_bytes(s.max),
                                                    md.type, 1, stat=True)[0]
            except:
                rv['max'] = None
        elif s.max_value is not None:
            try:
                if md.type == parquet_thrift.Type.BYTE_ARRAY:
                    rv['max'] = ensure_bytes(s.max_value)
                else:
                    rv['max'] = encoding.read_plain(ensure_bytes(s.max_value),
                                                    md.type, 1, stat=True)[0]
            except:
                rv['max'] = None
        if s.min is not None:
            try:
                if md.type == parquet_thrift.Type.BYTE_ARRAY:
                    rv['min'] = ensure_bytes(s.min)
                else:
                    rv['min'] = encoding.read_plain(ensure_bytes(s.min),
                                                    md.type, 1, stat=True)[0]
            except:
                rv['min'] = None
        elif s.min_value is not None:
            try:
                if md.type == parquet_thrift.Type.BYTE_ARRAY:
                    rv['min'] = ensure_bytes(s.min_value)
                else:
                    rv['min'] = encoding.read_plain(ensure_bytes(s.min_value),
                                                    md.type, 1, stat=True)[0]
            except:
                rv['min'] = None
        if s.null_count is not None:
            rv['null_count'] = s.null_count
        if s.distinct_count is not None:
            rv['distinct_count'] = s.distinct_count
        return rv

    if isinstance(obj, ThriftObject) and obj.thrift_name == "RowGroup":
        return {'.'.join(c.meta_data.path_in_schema): statistics(c)
                for c in obj.columns}

    if isinstance(obj, ParquetFile):
        L = list(map(statistics, obj.row_groups))
        d = {n: {col: [item.get(col, {}).get(n, None) for item in L]
                 for col in obj.columns}
             for n in ['min', 'max', 'null_count', 'distinct_count']}
        if not L:
            return d
        schema = obj.schema
        for col in obj.row_groups[0].columns:
            column = '.'.join(col.meta_data.path_in_schema)
            se = schema.schema_element(col.meta_data.path_in_schema)
            if (se.converted_type is not None or se.logicalType is not None
                    or se.type == parquet_thrift.Type.INT96):
                dtype = 'S12' if se.type == parquet_thrift.Type.INT96 else None
                for name in ['min', 'max']:
                    try:
                        d[name][column] = (
                            [None] if d[name][column] is None
                            or None in d[name][column]
                            else list(converted_types.convert(
                                np.array(d[name][column], dtype), se))
                        )
                    except (KeyError, ValueError):
                        # catch no stat and bad conversions
                        d[name][column] = [None]
        return d


def sorted_partitioned_columns(pf, filters=None):
    """
    The columns that are known to be sorted partition-by-partition

    They may not be sorted within each partition, but all elements in one
    row group are strictly greater than all elements in previous row groups.

    Examples
    --------
    >>> sorted_partitioned_columns(pf)
    {'id': {'min': [1, 5, 10], 'max': [4, 9, 20]}}

    Returns
    -------
    A set of column names

    See Also
    --------
    statistics
    """
    s = statistics(pf)
    if filters:
        rg_idx_list = filter_row_groups(pf, filters, as_idx=True)
        for stat in s.keys():
            for col in s[stat].keys():
                s[stat][col] = [s[stat][col][i] for i in rg_idx_list]
    columns = pf.columns
    out = dict()
    for c in columns:
        min, max = s['min'][c], s['max'][c]
        if any(x is None for x in min + max):
            continue
        try:
            if (min and sorted(min) == min and
                    sorted(max) == max and
                    all(mx < mn for mx, mn in zip(max[:-1], min[1:]))):
                out[c] = {'min': min, 'max': max}
        except TypeError:
            # because some types, e.g., dicts cannot be sorted/compared
            continue
    return out


def filter_row_groups(pf, filters, as_idx: bool = False):
    """
    Select row groups using set of filters.

    Parameters
    ----------
    pf: ParquetFile
        `ParquetFile` object.
    filters: list of list of tuples
        To filter out some of the row-groups.
        Filter syntax: [[(column, op, val), ...],...] where op is
        [==, >, >=, <, <=, !=, in, not in]
        The innermost tuples are transposed into a set of filters applied
        through an `AND` operation.
        The outer list combines these sets of filters through an `OR`
        operation.
        A single list of tuples can also be used, meaning that no `OR`
        operation between set of filters is to be conducted.
    as_idx: bool, False
        If `False`, returns a row group list, if `True`, returns row group
        index list.

    Returns
    -------
    Filtered list of row groups (or row group indexes)
    """
    filters = filters or [[]]
    # If 2nd level is already a column name, then transform
    # `filters` into a list (OR condition) of list (AND condition)
    # of filters (tuple or list with 1st component being a column
    # name).
    if filters[0] and isinstance(filters[0][0], str):
        filters = [filters]
    # Retrieve all column names onto which are applied filters, and check they
    # are existing columns of the dataset.
    as_cols = pf.columns + list(pf.cats.keys())
    known = [ands[0] in as_cols for ors in filters for ands in ors]
    if not all(known):
        falses = [i for i, x in enumerate(known) if not x]
        cols_in_filter = [ands[0] for ors in filters for ands in ors]
        wrong_cols = {cols_in_filter[i] for i in falses}
        raise ValueError('No filter can be applied on nonexistent column(s) \
{!s}.'.format(wrong_cols))
    if as_idx:
        return [i for i, rg in enumerate(pf.row_groups) if any([
                   not(filter_out_stats(rg, and_filters, pf.schema)) and
                   not(filter_out_cats(rg, and_filters, pf.partition_meta))
                   for and_filters in filters])]
    else:
        return [rg for rg in pf.row_groups if any([
                   not(filter_out_stats(rg, and_filters, pf.schema)) and
                   not(filter_out_cats(rg, and_filters, pf.partition_meta))
                   for and_filters in filters])]            


def filter_out_cats(rg, filters, partition_meta={}):
    """
    According to the filters, should this row-group be excluded

    Considers the partitioning category applicable to this row-group

    Parameters
    ----------
    rg: thrift RowGroup structure
    filters: list of 3-tuples
        Structure of each tuple: (column, op, value) where op is one of
        ['==', '!=', '<', '<=', '>', '>=', 'in', 'not in'] and value is
        appropriate for the column in question

    Returns
    -------
    True or False
    """
    if len(filters) == 0 or rg.columns[0].file_path is None:
        return False
    s = ex_from_sep('/')
    partitions = s.findall(rg.columns[0].file_path)
    pairs = [(p[0], p[1]) for p in partitions]
    for cat, v in pairs:

        app_filters = [f[1:] for f in filters if f[0] == cat]
        for op, val in app_filters:
            if isinstance(val, str) or (isinstance(val, (tuple, list)) and
                                        all(isinstance(x, str) for x in val)):
                v0 = v
            else:
                v0 = val_to_num(v)
            if cat in partition_meta:
                val = val_to_num(val, meta=partition_meta.get(cat))
                v0 = val_to_num(v0, meta=partition_meta.get(cat))
            if filter_val(op, val, v0, v0):
                return True
    return False


def filter_val(op, val, vmin=None, vmax=None):
    """
    Perform value comparison for filtering

    op: ['==', '!=', '<', '<=', '>', '>=', 'in', 'not in']
    val: appropriate value
    vmin, vmax: the range to compare within

    Returns
    -------
    True or False
    """
    vmin = _handle_np_array(vmin)
    vmax = _handle_np_array(vmax)
    if op == 'in':
        return filter_in(val, vmin, vmax)
    if op == 'not in':
        return filter_not_in(val, vmin, vmax)
    if vmax is not None:
        if op in ['==', '>=', '='] and val > vmax:
            return True
        if op == '>' and val >= vmax:
            return True
    if vmin is not None:
        if op in ['==', '<=', '='] and val < vmin:
            return True
        if op == '<' and val <= vmin:
            return True
    if (op == '!=' and vmax is not None and vmin is not None and
            vmax == vmin and val == vmax):
        return True

    # keep this row_group
    return False


def _handle_np_array(v):
    if v is not None and isinstance(v, np.ndarray):
        v = v[0]
    return v


def filter_in(values, vmin=None, vmax=None):
    """
    Handles 'in' filters

    op: ['in', 'not in']
    values: iterable of values
    vmin, vmax: the range to compare within

    Returns
    -------
    True or False
    """
    if len(values) == 0:
        return True
    if vmax == vmin and vmax is not None and vmax not in values:
        return True
    if vmin is None and vmax is None:
        return False

    sorted_values = sorted(values)
    if vmin is None and vmax is not None:
        return sorted_values[0] > vmax
    elif vmax is None and vmin is not None:
        return sorted_values[-1] < vmin

    vmin_insert = np.searchsorted(sorted_values, vmin, side='left')
    vmax_insert = np.searchsorted(sorted_values, vmax, side='right')

    # if the indexes are equal, then there are no values within the range
    return vmin_insert == vmax_insert


def filter_not_in(values, vmin=None, vmax=None):
    """
    Handles 'not in' filters

    op: ['in', 'not in']
    values: iterable of values
    vmin, vmax: the range to compare within

    Returns
    -------
    True or False
    """
    if len(values) == 0:
        return False
    if vmax is not None and vmax in values:
        return True
    elif vmin is not None and vmin in values:
        return True
    else:
        return False


def row_groups_map(rgs: list) -> dict:
    """
    Returns row group lists sorted by parquet files.

    Parameters
    ----------
    rgs: list
        List of row groups.

    Returns
    -------
    dict
        Per parquet file, list of row group stored in said file.
    """
    files_rgs = defaultdict(lambda: [])
    for rg in rgs:
        file = rg.columns[0].file_path
        files_rgs[file].append(rg)
    return files_rgs


def partitions(row_group, only_values=False) -> str:
    """Returns partition values as string.
    
    Values of partitions are separated with '/'.

    Parameters
    ----------
    row_group : obj or str
        Row group object or row group `file_path` as given by
        `rg.columns[0].file_path`.
    only_values: bool, default False
        If False, only values of partitions are returned;
        If True, names and values of partitions are returned (faster).

    Returns
    -------
    str
        Partitions values.
    """
    f_path = (row_group if isinstance(row_group, str)
              else row_group.columns[0].file_path)
    if '/' in f_path:
        return ('/'.join(re.split('/|=', f_path)[1::2]) if only_values
                else f_path.rsplit('/',1)[0])


def part_ids(row_groups) -> dict:
    """Return ids of parquet part files.
    
    Find the integer matching "**part.*.parquet" in referenced paths and
    returns them as keys of a dict.
    Values of the dict are tuples `(row_group_id, part_name)`.
    In case of files with multiple row groups, the position (index in row group
    list) of the 1st group only is kept.
    """
    max_rgidx = len(row_groups)-1
    paths = [rg.columns[0].file_path for rg in row_groups]
    matches = [(PART_ID.match(path), path) for path in paths]
    return {int(pid_path[0]['i']): (max_rgidx-i, pid_path[1])
            for i, pid_path in enumerate(reversed(matches))}
