import re
from collections import OrderedDict
from packaging.version import Version
import numpy as np
from pandas import (
    Categorical, DataFrame, Series,
    CategoricalIndex, RangeIndex, Index, MultiIndex,
    DatetimeIndex, CategoricalDtype,
    DatetimeTZDtype
)
from pandas.core.arrays.masked import BaseMaskedDtype
import warnings

from fastparquet.util import PANDAS_VERSION


class Dummy(object):
    pass


def empty(types, size, cats=None, cols=None, index_types=None, index_names=None,
          timezones=None, columns_dtype=None):
    """
    Create empty DataFrame to assign into

    In the simplest case, will return a Pandas dataframe of the given size,
    with columns of the given names and types. The second return value `views`
    is a dictionary of numpy arrays into which you can assign values that
    show up in the dataframe.

    For categorical columns, you get two views to assign into: if the
    column name is "col", you get both "col" (the category codes) and
    "col-catdef" (the category labels).

    For a single categorical index, you should use the `.set_categories`
    method of the appropriate "-catdef" columns, passing an Index of values

    ``views['index-catdef'].set_categories(pd.Index(newvalues), fastpath=True)``

    Multi-indexes work a lot like categoricals, even if the types of each
    index are not themselves categories, and will also have "-catdef" entries
    in the views. However, these will be Dummy instances, providing only a
    ``.set_categories`` method, to be used as above.

    Parameters
    ----------
    types: like np record structure, 'i4,u2,f4,f2,f4,M8,m8', or using tuples
        applies to non-categorical columns. If there are only categorical
        columns, an empty string of None will do.
    size: int
        Number of rows to allocate
    cats: dict {col: labels}
        Location and labels for categorical columns, e.g., {1: ['mary', 'mo]}
        will create column index 1 (inserted amongst the numerical columns)
        with two possible values. If labels is an integers, `{'col': 5}`,
        will generate temporary labels using range. If None, or column name
        is missing, will assume 16-bit integers (a reasonable default).
    cols: list of labels
        assigned column names, including categorical ones.
    index_types: list of str
        For one of more index columns, make them have this type. See general
        description, above, for caveats about multi-indexing. If None, the
        index will be the default RangeIndex.
    index_names: list of str
        Names of the index column(s), if using
    timezones: dict {col: timezone_str}
        for timestamp type columns, apply this timezone to the pandas series;
        the numpy view will be UTC.
    file_has_columns: bool, default False
        for files that are filtered but had columns before

    Returns
    -------
    - dataframe with correct shape and data-types
    - list of numpy views, in order, of the columns of the dataframe. Assign
        to this.
    """
    views = {}
    timezones = timezones or {}

    if isinstance(types, str):
        types = types.split(',')
    cols = cols if cols is not None else range(len(types))

    def cat(col):
        if cats is None or col not in cats:
            return RangeIndex(0, 2**14)
        elif isinstance(cats[col], int):
            return RangeIndex(0, cats[col])
        else:  # explicit labels list
            return cats[col]

    df = OrderedDict()
    for t, col in zip(types, cols):
        if str(t) == 'category':
            df[str(col)] = Categorical.from_codes([], categories=cat(col))
        elif isinstance(t, BaseMaskedDtype):
            # pandas masked types
            arr_type = t.construct_array_type()
            df[str(col)] = arr_type(
                values=np.empty(0, dtype=t.numpy_dtype),
                mask=np.empty(0, dtype=np.bool_),
                copy=False
            )
        else:
            if hasattr(t, 'base'):
                # funky pandas not-dtype
                t = t.base
            if ("M" in str(t) or "time" in str(t)) and "[" not in str(t):
                t = str(t) + "[ns]"
            d = np.empty(0, dtype=t)
            if d.dtype.kind == "M" and str(col) in timezones:
                try:
                    z = tz_to_dt_tz(timezones[str(col)])
                    d = Series(d).dt.tz_localize(z)
                except:
                    warnings.warn("Inferring time-zone from %s in column %s "
                                  "failed, using time-zone-agnostic"
                                  "" % (timezones[str(col)], col))
            df[str(col)] = d

    columns = Index(df.keys(), dtype=columns_dtype) if columns_dtype is not None else None
    df = DataFrame(df, columns=columns)
    if not index_types:
        index = RangeIndex(size)
    elif len(index_types) == 1:
        t, col = index_types[0], index_names[0]
        if col is None:
            raise ValueError('If using an index, must give an index name')
        if str(t) == 'category':
            # https://github.com/dask/fastparquet/issues/576#issuecomment-805579337
            temp = Categorical.from_codes([], categories=cat(col))
            vals = np.zeros(size, dtype=temp.codes.dtype)
            c = Categorical.from_codes(vals, dtype=temp.dtype)
            index = CategoricalIndex(c)

            views[col] = vals
            views[col+'-catdef'] = index._data
        else:
            if hasattr(t, 'base'):
                # funky pandas not-dtype
                 t = t.base
            # Initialize datetime index to zero: uninitialized data might fail
            # validation due to being an out-of-bounds datetime. xref
            # https://github.com/dask/fastparquet/issues/778
            dtype = np.dtype(t)
            d = np.zeros(size, dtype=dtype) if dtype.kind == "M" else np.empty(size, dtype=dtype)
            if d.dtype.kind == "M" and str(col) in timezones:
                # 1) create the DatetimeIndex in UTC as no datetime conversion is needed and
                # it works with d uninitialised data (no NonExistentTimeError or AmbiguousTimeError)
                # 2) convert to timezone (if UTC=noop, if None=remove tz, if other=change tz)
                index = DatetimeIndex(d, tz="UTC").tz_convert(
                    tz_to_dt_tz(timezones[str(col)]))
            else:
                index = Index(d)
            views[col] = d
    else:
        index = MultiIndex([[]], [[]])
        # index = MultiIndex.from_arrays(indexes)
        index._levels = list()
        index._labels = list()
        index._codes = list()
        index._names = list(index_names)
        for i, col in enumerate(index_names):
            index._levels.append(Index([None]))

            def set_cats(values, i=i, col=col, **kwargs):
                values.name = col
                if index._levels[i][0] is None:
                    index._levels[i] = values
                elif not index._levels[i].equals(values):
                    raise RuntimeError("Different dictionaries encountered"
                                       " while building categorical")

            x = Dummy()
            x._set_categories = set_cats
            x._multiindex = True

            d = np.zeros(size, dtype=int)
            if PANDAS_VERSION >= Version("0.24.0"):
                index._codes = list(index._codes) + [d]
            else:
                index._labels.append(d)
            views[col] = d
            views[col+'-catdef'] = x

    # Patch our blocks with desired-length arrays.  Kids: don't try this at home.
    mgr = df._mgr
    for block in mgr.blocks:
        bvalues = block.values
        shape = list(bvalues.shape)
        shape[-1] = size

        if isinstance(bvalues, Categorical):
            code = np.full(fill_value=-1, shape=shape, dtype=bvalues.codes.dtype)

            values = Categorical.from_codes(codes=code, dtype=bvalues.dtype)

        elif isinstance(bvalues.dtype, DatetimeTZDtype):
            dt = "M8[ns]" if PANDAS_VERSION.major < 2 else f'M8[{bvalues.dtype.unit}]'
            values = np.zeros(shape=shape, dtype=dt)
            values = type(bvalues)._from_sequence(values.view("int64"), copy=False, dtype=bvalues.dtype)
        else:
            if not isinstance(bvalues, np.ndarray):
                # e.g. DatetimeLikeBlock backed by DatetimeArray/TimedeltaArray
                if bvalues.dtype.kind == "m":
                    dt = "m8[ns]" if PANDAS_VERSION.major < 2 else bvalues.dtype
                    values = np.zeros(shape=shape, dtype=dt)
                    values = type(bvalues)._from_sequence(values.view("int64"), copy=False, dtype=bvalues.dtype)
                elif bvalues.dtype.kind == "M":
                    dt = "M8[ns]" if PANDAS_VERSION.major < 2 else bvalues.dtype
                    values = np.zeros(shape=shape, dtype=dt)
                    values = type(bvalues)._from_sequence(values.view("int64"), copy=False, dtype=bvalues.dtype)
                elif str(bvalues.dtype)[0] in {"I", "U"} or str(bvalues.dtype) == "boolean":
                    arr_type = bvalues.dtype.construct_array_type()
                    values = arr_type(
                        values=np.empty(size, dtype=bvalues.dtype.numpy_dtype),
                        mask=np.zeros(size, dtype=np.bool_)
                    )
                else:
                    raise NotImplementedError
            else:
                values = np.empty(shape=shape, dtype=bvalues.dtype)

        block.values = values

    mgr.axes[-1] = index

    # create views
    for block in df._mgr.blocks:
        dtype = block.dtype
        inds = block.mgr_locs.indexer
        if isinstance(inds, slice):
            inds = list(range(inds.start, inds.stop, inds.step))
        for i, ind in enumerate(inds):
            col = df.columns[ind]
            if isinstance(dtype, CategoricalDtype):
                views[col] = block.values._codes
                views[col+'-catdef'] = block.values
            elif getattr(block.dtype, 'tz', None):
                arr = np.asarray(block.values, dtype='M8[ns]')
                if len(arr.shape) > 1:
                    # pandas >= 1.3 does this for some reason
                    arr = arr.squeeze(axis=0)
                views[col] = arr
            elif str(dtype)[0] in {"I", "U"} or str(dtype) == "boolean":
                views[col] = block.values
            else:
                views[col] = block.values[i]

    if index_names:
        df.index.names = [
            None if re.match(r'__index_level_\d+__', n) else n
            for n in index_names
        ]
    return df, views


def tz_to_dt_tz(z):
    if ":" in z:
        import datetime
        hours, mins = z.split(":", 1)
        sign = z.startswith("-")
        z = int(hours) * 3600
        z += (1, -1)[sign] * int(mins) * 60
        z = datetime.timezone(datetime.timedelta(seconds=z))
    return z
