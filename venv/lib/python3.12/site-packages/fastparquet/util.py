from collections import defaultdict
import copy
from packaging.version import Version
from functools import lru_cache
import io
import struct
import os
import operator
import re
import numbers
import zoneinfo

import numpy as np
import pandas as pd

import fsspec

from fastparquet import parquet_thrift
from fastparquet.cencoding import ThriftObject
from fastparquet import __version__

PANDAS_VERSION = Version(pd.__version__)
created_by = f"fastparquet-python version {__version__} (build 0)"


class ParquetException(Exception):
    """Generic Exception related to unexpected data format when
     reading parquet file."""
    pass


def default_mkdirs(f):
    os.makedirs(f, exist_ok=True)


PATH_DATE_FMT = '%Y%m%d_%H%M%S.%f'


def path_string(o):
    if isinstance(o, pd.Timestamp):
        return o.isoformat()
    return str(o)


default_open = open


def default_remove(paths):
    for path in paths:
        try:
            os.unlink(path)
        except IOError:
            pass
    

def val_from_meta(x, meta):
    try:
        if meta['pandas_type'] == 'categorical':
            return x
        t = np.dtype(meta['numpy_type'])
        if t == "bool":
            return x in [True, "true", "True", 't', "T", 1, "1"]
        return np.dtype(t).type(x)
    except ValueError:
        if meta['numpy_type'] == 'datetime64[ns]':
            return pd.to_datetime(x, format=PATH_DATE_FMT)
        else:
            raise


def val_to_num(x, meta=None):
    """Parse a string as a number, date or timedelta if possible"""
    if meta:
        return val_from_meta(x, meta)
    return _val_to_num(x)


@lru_cache(1000)
def _val_to_num(x):
    if isinstance(x, numbers.Real):
        return x
    if x in ['now', 'NOW', 'TODAY', '']:
        return x
    if type(x) == str and x.lower() == 'nan':
        return x
    if x == "True":
        return True
    if x == "False":
        return False
    try:
        return int(x, base=10)
    except:
        pass
    try:
        return float(x)
    except:
        pass
    try:
        return pd.Timestamp(x)
    except:
        pass
    try:
        # TODO: determine the valid usecases for this, then try to limit the set
        #  ofstrings which may get inadvertently converted to timedeltas
        return pd.Timedelta(x)
    except:
        return x


def ensure_bytes(s):
    return s.encode('utf-8') if isinstance(s, str) else s


def ensure_str(b, *, ignore_error=False):
    if isinstance(b, str):
        return b
    else:
        try:
            return b.decode('utf-8')
        except (UnicodeDecodeError, AttributeError):
            if not ignore_error:
                raise
            return b


def check_column_names(columns, *args):
    """Ensure that parameters listing column names have corresponding columns"""
    for arg in args:
        if isinstance(arg, (tuple, list)):
            missing = set(arg) - set(columns)
            if missing:
                raise ValueError("Following columns were requested but are "
                                 "not available: %s.\n"
                                 "All requested columns: %s\n"
                                 "Available columns: %s"
                                 "" % (missing, arg, columns))


def reset_row_idx(data: pd.DataFrame) -> pd.DataFrame:
    """Reset row (multi-)index as column(s) of the DataFrame.

    Multi-index are stored in columns, one per index level.

    Parameters
    ----------
    data : dataframe

    Returns
    -------
    dataframe
    """
    if isinstance(data.index, pd.MultiIndex):
        for name, cats, codes in zip(data.index.names, data.index.levels,
                                     data.index.codes):
            data = data.assign(**{name: pd.Categorical.from_codes(codes,
                                                                  cats)})
        data.reset_index(drop=True)
    else:
        data = data.reset_index()
    return data


def metadata_from_many(file_list, verify_schema=False, open_with=default_open,
                       root=False, fs=None):
    """
    Given list of parquet files, make a FileMetaData that points to them

    Parameters
    ----------
    file_list: list of paths of parquet files
    verify_schema: bool (False)
        Whether to assert that the schemas in each file are identical
    open_with: function
        Use this to open each path.
    root: str
        Top of the dataset's directory tree, for cases where it can't be
        automatically inferred.
    fs: fsspsec.AbstractFileSystem
        Used in preference to open_with, if given

    Returns
    -------
    basepath: the root path that other paths are relative to
    fmd: metadata thrift structure
    """
    from fastparquet import api

    legacy = True
    if all(isinstance(pf, api.ParquetFile) for pf in file_list):
        pfs = file_list
        file_list = [pf.fn for pf in pfs]
    elif all(not isinstance(pf, api.ParquetFile) for pf in file_list):

        if verify_schema or fs is None or len(file_list) < 3:
            pfs = [api.ParquetFile(fn, open_with=open_with) for fn in file_list]
        else:
            # activate new code path here
            f0 = file_list[0]
            pf0 = api.ParquetFile(f0, open_with=open_with)
            if pf0.file_scheme not in ['empty', 'simple']:
                # set of directories, revert
                pfs = [pf0] + [api.ParquetFile(fn, open_with=open_with) for fn in file_list[1:]]
            else:
                # permits concurrent fetch of footers; needs fsspec >= 2021.6
                size = int(1.4 * pf0._head_size)
                pieces = fs.cat(file_list[1:], start=-size)
                sizes = {path: int.from_bytes(piece[-8:-4], "little") + 8 for
                         path, piece in pieces.items()}
                not_bigenough = [path for path, s in sizes.items() if s > size]
                if not_bigenough:
                    new_pieces = fs.cat(not_bigenough, start=-max(sizes.values()))
                    pieces.update(new_pieces)
                pieces = {k: _get_fmd(v) for k, v in pieces.items()}
                pieces = [(fn, pieces[fn]) for fn in file_list[1:]]  # recover ordering
                legacy = False
    else:
        raise ValueError("Merge requires all ParquetFile instances or none")
    basepath, file_list = analyse_paths(file_list, root=root)

    if legacy:
        # legacy code path
        if verify_schema:
            for pf in pfs[1:]:
                if pf._schema != pfs[0]._schema:
                    raise ValueError('Incompatible schemas')

        fmd = copy.copy(pfs[0].fmd)  # we inherit "created by" field
        rgs = []

        for pf, fn in zip(pfs, file_list):
            if pf.file_scheme not in ['simple', 'empty']:
                for rg in pf.row_groups:
                    rg = copy.copy(rg)
                    rg.columns = [copy.copy(c) for c in rg.columns]
                    for chunk in rg.columns:
                        chunk.file_path = '/'.join(
                            [fn, chunk.file_path if isinstance(chunk.file_path, str) else chunk.file_path.decode()]
                        )
                    rgs.append(rg)

            else:
                for rg in pf.row_groups:
                    rg = copy.copy(rg)
                    rg.columns = [copy.copy(c) for c in rg.columns]
                    for chunk in rg.columns:
                        chunk.file_path = fn
                    rgs.append(rg)

        fmd.row_groups = rgs
        fmd.num_rows = sum(rg.num_rows for rg in fmd.row_groups)
        return basepath, fmd

    for rg in pf0.fmd.row_groups:
        # chunks of first file, which would have file_path=None
        rg.columns[0].file_path = f0[len(basepath):].lstrip("/")

    rgs0 = pf0.fmd.row_groups
    for k, v in pieces:
        # Set file paths on other files
        if len(v.schema) > len(pf0.fmd.schema):
            # or was UPDATED with supercast
            pf0.fmd.schema = v.schema
        rgs = v.row_groups or []
        for rg in rgs:
            rg.columns[0].file_path = k[len(basepath):].lstrip("/")
        rgs0.extend(rgs)
    pf0.fmd.row_groups = rgs0
    pf0.fmd.num_rows = sum(rg.num_rows for rg in pf0.fmd.row_groups)
    return basepath, pf0.fmd


def _get_fmd(inbytes):
    from .cencoding import from_buffer

    f = io.BytesIO(inbytes)
    f.seek(-8, 2)
    head_size = struct.unpack('<i', f.read(4))[0]
    f.seek(-(head_size + 8), 2)
    data = f.read(head_size)
    return from_buffer(data, "FileMetaData")


def update_custom_metadata(obj, custom_metadata : dict):
    """Update custom metadata stored in thrift object or parquet file.

    Update strategy depends if key found in new custom metadata is also found
    in already existing custom metadata within thrift object, as well as its
    value.
        
      - If not found in existing, it is added.
      - If found in existing, it is updated.
      - If its value is `None`, it is not added, and if found in existing,
        it is removed from existing.

    Parameters
    ----------
    obj : metadata ThriftObject or parquet file
        Thrift object or parquet file which metadata is to update.
    custom_metadata : dict
        Key-value metadata to update in thrift object.
        The values must be strings or binary. To pass a dictionary, serialize it as json string then encode it in binary.
    Notes
    -----
    Key-value metadata are expected binary encoded. This function ensures it
    is.
    """
    kvm = (obj.key_value_metadata if isinstance(obj, ThriftObject)
           else obj.fmd.key_value_metadata)
    
    if kvm is None:
        kvm = []

    # Spare list of keys.
    kvm_keys = [item.key for item in kvm]
    for key, value in custom_metadata.items():
        key_b = ensure_bytes(key)
        if key_b in kvm_keys:
            idx = kvm_keys.index(key_b)
            if value is None:
                # Remove item.
                del kvm[idx]
                # Update 'kvm_keys' as well, for keeping indexing
                # up-to-date.
                del kvm_keys[idx]
            else:
                # Replace item.
                kvm[idx] = parquet_thrift.KeyValue(key=key_b,
                                                   value=ensure_bytes(value))
        elif value is not None:
            kvm.append(parquet_thrift.KeyValue(key=key_b,
                                               value=ensure_bytes(value)))
    if isinstance(obj, ThriftObject):
        obj.key_value_metadata = kvm
    else:
        obj.fmd.key_value_metadata = kvm
        # Reset '_kvm' to refresh 'key_value_metadata' cached property.
        obj._kvm = None


# simple cache to avoid re compile every time
seps = {}


def ex_from_sep(sep):
    """Generate regex for category folder matching"""
    if sep not in seps:
        if sep in r'\^$.|?*+()[]':
            s = re.compile(r"([a-zA-Z_0-9]+)=([^\\{}]+)".format(sep))
        else:
            s = re.compile("([a-zA-Z_0-9]+)=([^{}]+)".format(sep))
        seps[sep] = s
    return seps[sep]


def analyse_paths(file_list, root=False):
    """Consolidate list of file-paths into  parquet relative paths"""
    path_parts_list = [join_path(fn).split('/') for fn in file_list]
    if root is False:
        basepath = path_parts_list[0][:-1]
        for i, path_parts in enumerate(path_parts_list):
            j = len(path_parts) - 1
            for k, (base_part, path_part) in enumerate(
                    zip(basepath, path_parts)):
                if base_part != path_part:
                    j = k
                    break
            basepath = basepath[:j]
        l = len(basepath)

    else:
        basepath = join_path(root).split('/')
        l = len(basepath)
        assert all(p[:l] == basepath for p in path_parts_list
                   ), "All paths must begin with the given root"
    out_list = []
    for path_parts in path_parts_list:
        out_list.append('/'.join(path_parts[l:]))  # use '/'.join() instead of join_path to be consistent with split('/')

    return '/'.join(basepath), out_list  # use '/'.join() instead of join_path to be consistent with split('/')


def infer_dtype(column):
    try:
        return pd.api.types.infer_dtype(column, skipna=False)
    except AttributeError:
        return pd.lib.infer_dtype(column)


def groupby_types(iterable):
    groups = defaultdict(list)
    for x in iterable:
        groups[type(x)].append(x)
    return groups


def get_column_metadata(column, name, object_dtype=None):
    """Produce pandas column metadata block"""
    inferred_dtypes = {
        "utf8": "unicode",
        "bytes": "bytes",
        "bool": "bool",
        "int": "int",
        "json": "object",
        "bson": "object"
    }
    dtype = column.dtype
    if object_dtype in inferred_dtypes and dtype == "object":
        inferred_dtype = inferred_dtypes.get(object_dtype, "mixed")
    else:
        inferred_dtype = infer_dtype(column)
    if str(dtype) == "bool":
        # pandas accidentally calls this "boolean"
        inferred_dtype = "bool"

    if isinstance(dtype, pd.CategoricalDtype):
        extra_metadata = {
            'num_categories': len(column.cat.categories),
            'ordered': column.cat.ordered,
        }
        dtype = column.cat.codes.dtype
    elif isinstance(dtype, pd.DatetimeTZDtype):
        if isinstance(dtype.tz, zoneinfo.ZoneInfo):
            extra_metadata = {'timezone': dtype.tz.key}
        else:
            try:
                stz = str(dtype.tz)
                if "UTC" in stz and ":" in stz:
                    extra_metadata = {'timezone': stz.strip("UTC")}
                elif len(str(stz)) == 3:  # like "UTC", "CET", ...
                    extra_metadata = {'timezone': str(stz)}
                elif getattr(dtype.tz, "zone", False):
                    extra_metadata = {'timezone': dtype.tz.zone}
                elif "pytz" not in stz:
                    pd.Series([pd.to_datetime('now', utc=True)]).dt.tz_localize(stz)
                    extra_metadata = {'timezone': stz}
                elif "Offset" in stz:
                    extra_metadata = {'timezone': f"{dtype.tz._minutes // 60:+03}:00"}
                else:
                    raise KeyError
            except Exception as e:
                raise ValueError("Time-zone information could not be serialised: "
                                "%s, please use another" % str(dtype.tz)) from e
    else:
        extra_metadata = None

    if isinstance(name, tuple):
        name = str(name)
    elif not isinstance(name, str):
        raise TypeError(
            'Column name must be a string. Got column {} of type {}'.format(
                name, type(name).__name__
            )
        )

    return {
        'name': name,
        'field_name': name,
        'pandas_type': {
            'string': 'unicode',
            'datetime64': (
                'datetimetz' if hasattr(dtype, 'tz')
                else 'datetime'
            ),
            'integer': str(dtype),
            'floating': str(dtype),
        }.get(inferred_dtype, inferred_dtype),
        'numpy_type': get_numpy_type(dtype),
        'metadata': extra_metadata,
    }


def get_numpy_type(dtype):
    if isinstance(dtype, pd.CategoricalDtype):
        return 'category'
    elif "Int" in str(dtype):
        return str(dtype).lower()
    elif str(dtype) == "boolean":
        return "bool"
    elif str(dtype) == "string":
        return "object"
    else:
        return str(dtype)


def get_file_scheme(paths):
    """For the given row groups, figure out if the partitioning scheme

    Parameters
    ----------
    paths: list of str
        normally from row_group.columns[0].file_path

    Returns
    -------
    'empty': no rgs at all
    'simple': all rgs in a single file
    'flat': multiple files in one directory
    'hive': directories are all `key=value`; all files are at the same
        directory depth
    'drill': assume directory names are labels, and field names are of the
        form dir0, dir1; all files are at the same directory depth
    'other': none of the above, assume no partitioning
    """
    if not paths:
        return 'empty'
    if set(paths) == {None}:
        return 'simple'
    if None in paths:
        return 'other'
    parts = [p.split('/') for p in paths]
    lens = [len(p) for p in parts]
    if len(set(lens)) > 1:
        return 'other'
    if set(lens) == {1}:
        return 'flat'
    matches = all(all("=" in p[1:-1] for p in part[:-1]) for part in parts)
    return "hive" if matches else "drill"


def join_path(*path):
    return "/".join([str(p).replace("\\", "/").rstrip("/") for p in path if p])


def _strip_path_tail(paths) -> set:
    return {path.rsplit("/", 1)[0] if "/" in path else "" for path in paths}


ops = {
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le
}


def norm_col_name(name, is_index:bool=None):
    if isinstance(name, tuple):
        if is_index:
            return name[0]
        else:
            return str(name)
    return name


def get_fs(fn, open_with, mkdirs):
    fs = None
    if "FastParquetImpl.write.<locals>.<lambda>" in str(open_with):
        import inspect
        so = inspect.getclosurevars(open_with).nonlocals["storage_options"] or {}
        fs, fn = fsspec.core.url_to_fs(fn, **so)
        open_with = fs.open
        mkdirs = mkdirs or (lambda d: fs.mkdirs(d, exist_ok=True))
    return fs, fn, open_with, mkdirs

