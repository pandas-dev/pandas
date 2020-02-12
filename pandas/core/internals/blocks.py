from datetime import datetime, timedelta
import functools
import inspect
import re
from typing import Any, List
import warnings

import numpy as np

from pandas._libs import NaT, Timestamp, algos as libalgos, lib, tslib, writers
import pandas._libs.internals as libinternals
from pandas._libs.tslibs import Timedelta, conversion
from pandas._libs.tslibs.timezones import tz_compare
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import (
    astype_nansafe,
    convert_scalar_for_putitemlike,
    find_common_type,
    infer_dtype_from,
    infer_dtype_from_scalar,
    maybe_downcast_numeric,
    maybe_downcast_to_dtype,
    maybe_infer_dtype_type,
    maybe_promote,
    maybe_upcast,
    soft_convert_objects,
)
from pandas.core.dtypes.common import (
    _NS_DTYPE,
    _TD_DTYPE,
    ensure_platform_int,
    is_bool_dtype,
    is_categorical,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_interval_dtype,
    is_list_like,
    is_object_dtype,
    is_period_dtype,
    is_re,
    is_re_compilable,
    is_sparse,
    is_timedelta64_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.concat import concat_categorical, concat_datetime
from pandas.core.dtypes.dtypes import CategoricalDtype, ExtensionDtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCExtensionArray,
    ABCPandasArray,
    ABCSeries,
)
from pandas.core.dtypes.missing import (
    _isna_compat,
    array_equivalent,
    is_valid_nat_for_dtype,
    isna,
)

import pandas.core.algorithms as algos
from pandas.core.arrays import (
    Categorical,
    DatetimeArray,
    ExtensionArray,
    PandasArray,
    PandasDtype,
    TimedeltaArray,
)
from pandas.core.base import PandasObject
import pandas.core.common as com
from pandas.core.construction import extract_array
from pandas.core.indexers import (
    check_setitem_lengths,
    is_empty_indexer,
    is_scalar_indexer,
)
import pandas.core.missing as missing
from pandas.core.nanops import nanpercentile


class Block(PandasObject):
    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas
    data structure

    Index-ignorant; let the container take care of that
    """

    __slots__ = ["_mgr_locs", "values", "ndim"]
    is_numeric = False
    is_float = False
    is_integer = False
    is_complex = False
    is_datetime = False
    is_datetimetz = False
    is_timedelta = False
    is_bool = False
    is_object = False
    is_categorical = False
    is_extension = False
    _can_hold_na = False
    _can_consolidate = True
    _verify_integrity = True
    _validate_ndim = True
    _ftype = "dense"
    _concatenator = staticmethod(np.concatenate)

    def __init__(self, values, placement, ndim=None):
        self.ndim = self._check_ndim(values, ndim)
        self.mgr_locs = placement
        self.values = values

        if self._validate_ndim and self.ndim and len(self.mgr_locs) != len(self.values):
            raise ValueError(
                f"Wrong number of items passed {len(self.values)}, "
                f"placement implies {len(self.mgr_locs)}"
            )

    def _check_ndim(self, values, ndim):
        """
        ndim inference and validation.

        Infers ndim from 'values' if not provided to __init__.
        Validates that values.ndim and ndim are consistent if and only if
        the class variable '_validate_ndim' is True.

        Parameters
        ----------
        values : array-like
        ndim : int or None

        Returns
        -------
        ndim : int

        Raises
        ------
        ValueError : the number of dimensions do not match
        """
        if ndim is None:
            ndim = values.ndim

        if self._validate_ndim and values.ndim != ndim:
            raise ValueError(
                "Wrong number of dimensions. "
                f"values.ndim != ndim [{values.ndim} != {ndim}]"
            )
        return ndim

    @property
    def _holder(self):
        """
        The array-like that can hold the underlying values.

        None for 'Block', overridden by subclasses that don't
        use an ndarray.
        """
        return None

    @property
    def _consolidate_key(self):
        return (self._can_consolidate, self.dtype.name)

    @property
    def _is_single_block(self):
        return self.ndim == 1

    @property
    def is_view(self):
        """ return a boolean if I am possibly a view """
        return self.values.base is not None

    @property
    def is_datelike(self):
        """ return True if I am a non-datelike """
        return self.is_datetime or self.is_timedelta

    def is_categorical_astype(self, dtype):
        """
        validate that we have a astypeable to categorical,
        returns a boolean if we are a categorical
        """
        if dtype is Categorical or dtype is CategoricalDtype:
            # this is a pd.Categorical, but is not
            # a valid type for astypeing
            raise TypeError(f"invalid type {dtype} for astype")

        elif is_categorical_dtype(dtype):
            return True

        return False

    def external_values(self):
        """
        The array that Series.values returns (public attribute).

        This has some historical constraints, and is overridden in block
        subclasses to return the correct array (e.g. period returns
        object ndarray and datetimetz a datetime64[ns] ndarray instead of
        proper extension array).
        """
        return self.values

    def internal_values(self):
        """
        The array that Series._values returns (internal values).
        """
        return self.values

    def array_values(self) -> ExtensionArray:
        """
        The array that Series.array returns. Always an ExtensionArray.
        """
        return PandasArray(self.values)

    def get_values(self, dtype=None):
        """
        return an internal format, currently just the ndarray
        this is often overridden to handle to_dense like operations
        """
        if is_object_dtype(dtype):
            return self.values.astype(object)
        return self.values

    def get_block_values(self, dtype=None):
        """
        This is used in the JSON C code
        """
        return self.get_values(dtype=dtype)

    def to_dense(self):
        return self.values.view()

    @property
    def fill_value(self):
        return np.nan

    @property
    def mgr_locs(self):
        return self._mgr_locs

    @mgr_locs.setter
    def mgr_locs(self, new_mgr_locs):
        if not isinstance(new_mgr_locs, libinternals.BlockPlacement):
            new_mgr_locs = libinternals.BlockPlacement(new_mgr_locs)

        self._mgr_locs = new_mgr_locs

    @property
    def array_dtype(self):
        """ the dtype to return if I want to construct this block as an
        array
        """
        return self.dtype

    def make_block(self, values, placement=None) -> "Block":
        """
        Create a new block, with type inference propagate any values that are
        not specified
        """
        if placement is None:
            placement = self.mgr_locs

        return make_block(values, placement=placement, ndim=self.ndim)

    def make_block_same_class(self, values, placement=None, ndim=None):
        """ Wrap given values in a block of same type as self. """
        if placement is None:
            placement = self.mgr_locs
        if ndim is None:
            ndim = self.ndim
        return make_block(values, placement=placement, ndim=ndim, klass=type(self))

    def __repr__(self) -> str:
        # don't want to print out all of the items here
        name = type(self).__name__
        if self._is_single_block:
            result = f"{name}: {len(self)} dtype: {self.dtype}"
        else:

            shape = " x ".join(str(s) for s in self.shape)
            result = f"{name}: {self.mgr_locs.indexer}, {shape}, dtype: {self.dtype}"

        return result

    def __len__(self) -> int:
        return len(self.values)

    def __getstate__(self):
        return self.mgr_locs.indexer, self.values

    def __setstate__(self, state):
        self.mgr_locs = libinternals.BlockPlacement(state[0])
        self.values = state[1]
        self.ndim = self.values.ndim

    def _slice(self, slicer):
        """ return a slice of my values """
        return self.values[slicer]

    def getitem_block(self, slicer, new_mgr_locs=None):
        """
        Perform __getitem__-like, return result as block.

        As of now, only supports slices that preserve dimensionality.
        """
        if new_mgr_locs is None:
            axis0_slicer = slicer[0] if isinstance(slicer, tuple) else slicer
            new_mgr_locs = self.mgr_locs[axis0_slicer]

        new_values = self._slice(slicer)

        if self._validate_ndim and new_values.ndim != self.ndim:
            raise ValueError("Only same dim slicing is allowed")

        return self.make_block_same_class(new_values, new_mgr_locs)

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    @property
    def ftype(self):
        if getattr(self.values, "_pandas_ftype", False):
            dtype = self.dtype.subtype
        else:
            dtype = self.dtype
        return f"{dtype}:{self._ftype}"

    def merge(self, other):
        return _merge_blocks([self, other])

    def concat_same_type(self, to_concat, placement=None):
        """
        Concatenate list of single blocks of the same type.
        """
        values = self._concatenator(
            [blk.values for blk in to_concat], axis=self.ndim - 1
        )
        return self.make_block_same_class(
            values, placement=placement or slice(0, len(values), 1)
        )

    def iget(self, i):
        return self.values[i]

    def set(self, locs, values):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        self.values[locs] = values

    def delete(self, loc):
        """
        Delete given loc(-s) from block in-place.
        """
        self.values = np.delete(self.values, loc, 0)
        self.mgr_locs = self.mgr_locs.delete(loc)

    def apply(self, func, **kwargs) -> List["Block"]:
        """ apply the function to my values; return a block if we are not
        one
        """
        with np.errstate(all="ignore"):
            result = func(self.values, **kwargs)

        return self._split_op_result(result)

    def _split_op_result(self, result) -> List["Block"]:
        # See also: split_and_operate
        if is_extension_array_dtype(result) and result.ndim > 1:
            # if we get a 2D ExtensionArray, we need to split it into 1D pieces
            nbs = []
            for i, loc in enumerate(self.mgr_locs):
                vals = result[i]
                nv = _block_shape(vals, ndim=self.ndim)
                block = self.make_block(values=nv, placement=[loc])
                nbs.append(block)
            return nbs

        if not isinstance(result, Block):
            result = self.make_block(values=_block_shape(result, ndim=self.ndim))

        return [result]

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        """ fillna on the block with the value. If we fail, then convert to
        ObjectBlock and try again
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        mask = isna(self.values)
        if limit is not None:
            limit = libalgos._validate_limit(None, limit=limit)
            mask[mask.cumsum(self.ndim - 1) > limit] = False

        if not self._can_hold_na:
            if inplace:
                return self
            else:
                return self.copy()

        if self._can_hold_element(value):
            # equivalent: _try_coerce_args(value) would not raise
            blocks = self.putmask(mask, value, inplace=inplace)
            return self._maybe_downcast(blocks, downcast)

        # we can't process the value, but nothing to do
        if not mask.any():
            return self if inplace else self.copy()

        # operate column-by-column
        def f(mask, val, idx):
            block = self.coerce_to_target_dtype(value)

            # slice out our block
            if idx is not None:
                # i.e. self.ndim == 2
                block = block.getitem_block(slice(idx, idx + 1))
            return block.fillna(value, limit=limit, inplace=inplace, downcast=None)

        return self.split_and_operate(None, f, inplace)

    def split_and_operate(self, mask, f, inplace: bool):
        """
        split the block per-column, and apply the callable f
        per-column, return a new block for each. Handle
        masking which will not change a block unless needed.

        Parameters
        ----------
        mask : 2-d boolean mask
        f : callable accepting (1d-mask, 1d values, indexer)
        inplace : boolean

        Returns
        -------
        list of blocks
        """
        if mask is None:
            mask = np.broadcast_to(True, shape=self.shape)

        new_values = self.values

        def make_a_block(nv, ref_loc):
            if isinstance(nv, list):
                assert len(nv) == 1, nv
                assert isinstance(nv[0], Block)
                block = nv[0]
            else:
                # Put back the dimension that was taken from it and make
                # a block out of the result.
                nv = _block_shape(nv, ndim=self.ndim)
                block = self.make_block(values=nv, placement=ref_loc)
            return block

        # ndim == 1
        if self.ndim == 1:
            if mask.any():
                nv = f(mask, new_values, None)
            else:
                nv = new_values if inplace else new_values.copy()
            block = make_a_block(nv, self.mgr_locs)
            return [block]

        # ndim > 1
        new_blocks = []
        for i, ref_loc in enumerate(self.mgr_locs):
            m = mask[i]
            v = new_values[i]

            # need a new block
            if m.any():
                nv = f(m, v, i)
            else:
                nv = v if inplace else v.copy()

            block = make_a_block(nv, [ref_loc])
            new_blocks.append(block)

        return new_blocks

    def _maybe_downcast(self, blocks: List["Block"], downcast=None) -> List["Block"]:

        # no need to downcast our float
        # unless indicated
        if downcast is None and (
            self.is_float or self.is_timedelta or self.is_datetime
        ):
            return blocks

        return _extend_blocks([b.downcast(downcast) for b in blocks])

    def downcast(self, dtypes=None):
        """ try to downcast each item to the dict of dtypes if present """
        # turn it off completely
        if dtypes is False:
            return self

        values = self.values

        # single block handling
        if self._is_single_block:

            # try to cast all non-floats here
            if dtypes is None:
                dtypes = "infer"

            nv = maybe_downcast_to_dtype(values, dtypes)
            return self.make_block(nv)

        # ndim > 1
        if dtypes is None:
            return self

        if not (dtypes == "infer" or isinstance(dtypes, dict)):
            raise ValueError(
                "downcast must have a dictionary or 'infer' as its argument"
            )
        elif dtypes != "infer":
            raise AssertionError("dtypes as dict is not supported yet")

        # operate column-by-column
        # this is expensive as it splits the blocks items-by-item
        def f(mask, val, idx):
            val = maybe_downcast_to_dtype(val, dtype="infer")
            return val

        return self.split_and_operate(None, f, False)

    def astype(self, dtype, copy: bool = False, errors: str = "raise"):
        """
        Coerce to the new dtype.

        Parameters
        ----------
        dtype : str, dtype convertible
        copy : bool, default False
            copy if indicated
        errors : str, {'raise', 'ignore'}, default 'ignore'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object

        Returns
        -------
        Block
        """
        errors_legal_values = ("raise", "ignore")

        if errors not in errors_legal_values:
            invalid_arg = (
                "Expected value of kwarg 'errors' to be one of "
                f"{list(errors_legal_values)}. Supplied value is '{errors}'"
            )
            raise ValueError(invalid_arg)

        if inspect.isclass(dtype) and issubclass(dtype, ExtensionDtype):
            msg = (
                f"Expected an instance of {dtype.__name__}, "
                "but got the class instead. Try instantiating 'dtype'."
            )
            raise TypeError(msg)

        # may need to convert to categorical
        if self.is_categorical_astype(dtype):

            if is_categorical_dtype(self.values):
                # GH 10696/18593: update an existing categorical efficiently
                return self.make_block(self.values.astype(dtype, copy=copy))

            return self.make_block(Categorical(self.values, dtype=dtype))

        dtype = pandas_dtype(dtype)

        # astype processing
        if is_dtype_equal(self.dtype, dtype):
            if copy:
                return self.copy()
            return self

        # force the copy here
        if self.is_extension:
            # TODO: Should we try/except this astype?
            values = self.values.astype(dtype)
        else:
            if issubclass(dtype.type, str):

                # use native type formatting for datetime/tz/timedelta
                if self.is_datelike:
                    values = self.to_native_types()

                # astype formatting
                else:
                    values = self.get_values()

            else:
                values = self.get_values(dtype=dtype)

            # _astype_nansafe works fine with 1-d only
            vals1d = values.ravel()
            try:
                values = astype_nansafe(vals1d, dtype, copy=True)
            except (ValueError, TypeError):
                # e.g. astype_nansafe can fail on object-dtype of strings
                #  trying to convert to float
                if errors == "raise":
                    raise
                newb = self.copy() if copy else self
                return newb

        # TODO(extension)
        # should we make this attribute?
        if isinstance(values, np.ndarray):
            values = values.reshape(self.shape)

        newb = make_block(values, placement=self.mgr_locs, ndim=self.ndim)

        if newb.is_numeric and self.is_numeric:
            if newb.shape != self.shape:
                raise TypeError(
                    f"cannot set astype for copy = [{copy}] for dtype "
                    f"({self.dtype.name} [{self.shape}]) to different shape "
                    f"({newb.dtype.name} [{newb.shape}])"
                )
        return newb

    def convert(
        self,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
        coerce: bool = False,
    ):
        """ attempt to coerce any object types to better types return a copy
        of the block (if copy = True) by definition we are not an ObjectBlock
        here!
        """
        return self.copy() if copy else self

    def _can_hold_element(self, element: Any) -> bool:
        """ require the same dtype as ourselves """
        dtype = self.values.dtype.type
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, dtype)
        return isinstance(element, dtype)

    def to_native_types(self, slicer=None, na_rep="nan", quoting=None, **kwargs):
        """ convert to our native types format, slicing if desired """
        values = self.get_values()

        if slicer is not None:
            values = values[:, slicer]
        mask = isna(values)
        itemsize = writers.word_len(na_rep)

        if not self.is_object and not quoting and itemsize:
            values = values.astype(str)
            if values.dtype.itemsize / np.dtype("U1").itemsize < itemsize:
                # enlarge for the na_rep
                values = values.astype(f"<U{itemsize}")
        else:
            values = np.array(values, dtype="object")

        values[mask] = na_rep
        return values

    # block actions #
    def copy(self, deep=True):
        """ copy constructor """
        values = self.values
        if deep:
            values = values.copy()
        return self.make_block_same_class(values, ndim=self.ndim)

    def replace(
        self, to_replace, value, inplace=False, filter=None, regex=False, convert=True
    ):
        """replace the to_replace value with value, possible to create new
        blocks here this is just a call to putmask. regex is not used here.
        It is used in ObjectBlocks.  It is here for API compatibility.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        original_to_replace = to_replace

        # If we cannot replace with own dtype, convert to ObjectBlock and
        # retry
        if not self._can_hold_element(to_replace):
            if not isinstance(to_replace, list):
                if inplace:
                    return [self]
                return [self.copy()]

            to_replace = [x for x in to_replace if self._can_hold_element(x)]
            if not len(to_replace):
                # GH#28084 avoid costly checks since we can infer
                #  that there is nothing to replace in this block
                if inplace:
                    return [self]
                return [self.copy()]

            if len(to_replace) == 1:
                # _can_hold_element checks have reduced this back to the
                #  scalar case and we can avoid a costly object cast
                return self.replace(
                    to_replace[0],
                    value,
                    inplace=inplace,
                    filter=filter,
                    regex=regex,
                    convert=convert,
                )

            # GH 22083, TypeError or ValueError occurred within error handling
            # causes infinite loop. Cast and retry only if not objectblock.
            if is_object_dtype(self):
                raise AssertionError

            # try again with a compatible block
            block = self.astype(object)
            return block.replace(
                to_replace=to_replace,
                value=value,
                inplace=inplace,
                filter=filter,
                regex=regex,
                convert=convert,
            )

        values = self.values
        if lib.is_scalar(to_replace) and isinstance(values, np.ndarray):
            # The only non-DatetimeLike class that also has a non-trivial
            #  try_coerce_args is ObjectBlock, but that overrides replace,
            #  so does not get here.
            to_replace = convert_scalar_for_putitemlike(to_replace, values.dtype)

        mask = missing.mask_missing(values, to_replace)
        if filter is not None:
            filtered_out = ~self.mgr_locs.isin(filter)
            mask[filtered_out.nonzero()[0]] = False

        if not mask.any():
            if inplace:
                return [self]
            return [self.copy()]

        try:
            blocks = self.putmask(mask, value, inplace=inplace)
            # Note: it is _not_ the case that self._can_hold_element(value)
            #  is always true at this point.  In particular, that can fail
            #  for:
            #   "2u" with bool-dtype, float-dtype
            #   0.5 with int64-dtype
            #   np.nan with int64-dtype
        except (TypeError, ValueError):
            # GH 22083, TypeError or ValueError occurred within error handling
            # causes infinite loop. Cast and retry only if not objectblock.
            if is_object_dtype(self):
                raise

            assert not self._can_hold_element(value), value

            # try again with a compatible block
            block = self.astype(object)
            return block.replace(
                to_replace=original_to_replace,
                value=value,
                inplace=inplace,
                filter=filter,
                regex=regex,
                convert=convert,
            )
        if convert:
            blocks = [b.convert(numeric=False, copy=not inplace) for b in blocks]
        return blocks

    def _replace_single(self, *args, **kwargs):
        """ no-op on a non-ObjectBlock """
        return self if kwargs["inplace"] else self.copy()

    def setitem(self, indexer, value):
        """
        Set the value inplace, returning a a maybe different typed block.

        Parameters
        ----------
        indexer : tuple, list-like, array-like, slice
            The subset of self.values to set
        value : object
            The value being set

        Returns
        -------
        Block

        Notes
        -----
        `indexer` is a direct slice/positional indexer. `value` must
        be a compatible shape.
        """
        transpose = self.ndim == 2

        if isinstance(indexer, np.ndarray) and indexer.ndim > self.ndim:
            raise ValueError(f"Cannot set values with ndim > {self.ndim}")

        # coerce None values, if appropriate
        if value is None:
            if self.is_numeric:
                value = np.nan

        # coerce if block dtype can store value
        values = self.values
        if self._can_hold_element(value):
            # We only get here for non-Extension Blocks, so _try_coerce_args
            #  is only relevant for DatetimeBlock and TimedeltaBlock
            if lib.is_scalar(value):
                value = convert_scalar_for_putitemlike(value, values.dtype)

        else:
            # current dtype cannot store value, coerce to common dtype
            find_dtype = False

            if hasattr(value, "dtype"):
                dtype = value.dtype
                find_dtype = True

            elif lib.is_scalar(value) and not isna(value):
                dtype, _ = infer_dtype_from_scalar(value, pandas_dtype=True)
                find_dtype = True

            if find_dtype:
                dtype = find_common_type([values.dtype, dtype])
                if not is_dtype_equal(self.dtype, dtype):
                    b = self.astype(dtype)
                    return b.setitem(indexer, value)

        # value must be storeable at this moment
        if is_extension_array_dtype(getattr(value, "dtype", None)):
            # We need to be careful not to allow through strings that
            #  can be parsed to EADtypes
            arr_value = value
        else:
            arr_value = np.array(value)

        # cast the values to a type that can hold nan (if necessary)
        if not self._can_hold_element(value):
            dtype, _ = maybe_promote(arr_value.dtype)
            values = values.astype(dtype)

        if transpose:
            values = values.T

        # length checking
        check_setitem_lengths(indexer, value, values)
        exact_match = (
            len(arr_value.shape)
            and arr_value.shape[0] == values.shape[0]
            and arr_value.size == values.size
        )
        if is_empty_indexer(indexer, arr_value):
            # GH#8669 empty indexers
            pass

        elif is_scalar_indexer(indexer, arr_value):
            # setting a single element for each dim and with a rhs that could
            #  be e.g. a list; see GH#6043
            values[indexer] = value

        elif (
            exact_match
            and is_categorical_dtype(arr_value.dtype)
            and not is_categorical_dtype(values)
        ):
            # GH25495 - If the current dtype is not categorical,
            # we need to create a new categorical block
            values[indexer] = value
            return self.make_block(Categorical(self.values, dtype=arr_value.dtype))

        # if we are an exact match (ex-broadcasting),
        # then use the resultant dtype
        elif exact_match:
            values[indexer] = value

            try:
                values = values.astype(arr_value.dtype)
            except ValueError:
                pass

        # set
        else:
            values[indexer] = value

        if transpose:
            values = values.T
        block = self.make_block(values)
        return block

    def putmask(self, mask, new, align=True, inplace=False, axis=0, transpose=False):
        """ putmask the data to the block; it is possible that we may create a
        new dtype of block

        return the resulting block(s)

        Parameters
        ----------
        mask  : the condition to respect
        new : a ndarray/object
        align : boolean, perform alignment on other/cond, default is True
        inplace : perform inplace modification, default is False
        axis : int
        transpose : boolean
            Set to True if self is stored with axes reversed

        Returns
        -------
        a list of new blocks, the result of the putmask
        """
        new_values = self.values if inplace else self.values.copy()

        new = getattr(new, "values", new)
        mask = getattr(mask, "values", mask)

        # if we are passed a scalar None, convert it here
        if not is_list_like(new) and isna(new) and not self.is_object:
            # FIXME: make sure we have compatible NA
            new = self.fill_value

        if self._can_hold_element(new):
            # We only get here for non-Extension Blocks, so _try_coerce_args
            #  is only relevant for DatetimeBlock and TimedeltaBlock
            if lib.is_scalar(new):
                new = convert_scalar_for_putitemlike(new, new_values.dtype)

            if transpose:
                new_values = new_values.T

            # If the default repeat behavior in np.putmask would go in the
            # wrong direction, then explicitly repeat and reshape new instead
            if getattr(new, "ndim", 0) >= 1:
                if self.ndim - 1 == new.ndim and axis == 1:
                    new = np.repeat(new, new_values.shape[-1]).reshape(self.shape)
                new = new.astype(new_values.dtype)

            # we require exact matches between the len of the
            # values we are setting (or is compat). np.putmask
            # doesn't check this and will simply truncate / pad
            # the output, but we want sane error messages
            #
            # TODO: this prob needs some better checking
            # for 2D cases
            if (
                is_list_like(new)
                and np.any(mask[mask])
                and getattr(new, "ndim", 1) == 1
            ):
                if mask[mask].shape[-1] == len(new):
                    # GH 30567
                    # If length of ``new`` is less than the length of ``new_values``,
                    # `np.putmask` would first repeat the ``new`` array and then
                    # assign the masked values hence produces incorrect result.
                    # `np.place` on the other hand uses the ``new`` values at it is
                    # to place in the masked locations of ``new_values``
                    np.place(new_values, mask, new)
                elif mask.shape[-1] == len(new) or len(new) == 1:
                    np.putmask(new_values, mask, new)
                else:
                    raise ValueError("cannot assign mismatch length to masked array")
            else:
                np.putmask(new_values, mask, new)

        # maybe upcast me
        elif mask.any():
            if transpose:
                mask = mask.T
                if isinstance(new, np.ndarray):
                    new = new.T
                axis = new_values.ndim - axis - 1

            # Pseudo-broadcast
            if getattr(new, "ndim", 0) >= 1:
                if self.ndim - 1 == new.ndim:
                    new_shape = list(new.shape)
                    new_shape.insert(axis, 1)
                    new = new.reshape(tuple(new_shape))

            # operate column-by-column
            def f(mask, val, idx):

                if idx is None:
                    # ndim==1 case.
                    n = new
                else:

                    if isinstance(new, np.ndarray):
                        n = np.squeeze(new[idx % new.shape[0]])
                    else:
                        n = np.array(new)

                    # type of the new block
                    dtype, _ = maybe_promote(n.dtype)

                    # we need to explicitly astype here to make a copy
                    n = n.astype(dtype)

                nv = _putmask_smart(val, mask, n)
                return nv

            new_blocks = self.split_and_operate(mask, f, inplace)
            return new_blocks

        if inplace:
            return [self]

        if transpose:
            new_values = new_values.T

        return [self.make_block(new_values)]

    def coerce_to_target_dtype(self, other):
        """
        coerce the current block to a dtype compat for other
        we will return a block, possibly object, and not raise

        we can also safely try to coerce to the same dtype
        and will receive the same block
        """
        # if we cannot then coerce to object
        dtype, _ = infer_dtype_from(other, pandas_dtype=True)

        if is_dtype_equal(self.dtype, dtype):
            return self

        if self.is_bool or is_object_dtype(dtype) or is_bool_dtype(dtype):
            # we don't upcast to bool
            return self.astype(object)

        elif (self.is_float or self.is_complex) and (
            is_integer_dtype(dtype) or is_float_dtype(dtype)
        ):
            # don't coerce float/complex to int
            return self

        elif (
            self.is_datetime
            or is_datetime64_dtype(dtype)
            or is_datetime64tz_dtype(dtype)
        ):

            # not a datetime
            if not (
                (is_datetime64_dtype(dtype) or is_datetime64tz_dtype(dtype))
                and self.is_datetime
            ):
                return self.astype(object)

            # don't upcast timezone with different timezone or no timezone
            mytz = getattr(self.dtype, "tz", None)
            othertz = getattr(dtype, "tz", None)

            if not tz_compare(mytz, othertz):
                return self.astype(object)

            raise AssertionError(
                f"possible recursion in coerce_to_target_dtype: {self} {other}"
            )

        elif self.is_timedelta or is_timedelta64_dtype(dtype):

            # not a timedelta
            if not (is_timedelta64_dtype(dtype) and self.is_timedelta):
                return self.astype(object)

            raise AssertionError(
                f"possible recursion in coerce_to_target_dtype: {self} {other}"
            )

        try:
            return self.astype(dtype)
        except (ValueError, TypeError, OverflowError):
            return self.astype(object)

    def interpolate(
        self,
        method="pad",
        axis=0,
        index=None,
        values=None,
        inplace=False,
        limit=None,
        limit_direction="forward",
        limit_area=None,
        fill_value=None,
        coerce=False,
        downcast=None,
        **kwargs,
    ):

        inplace = validate_bool_kwarg(inplace, "inplace")

        def check_int_bool(self, inplace):
            # Only FloatBlocks will contain NaNs.
            # timedelta subclasses IntBlock
            if (self.is_bool or self.is_integer) and not self.is_timedelta:
                if inplace:
                    return self
                else:
                    return self.copy()

        # a fill na type method
        try:
            m = missing.clean_fill_method(method)
        except ValueError:
            m = None

        if m is not None:
            r = check_int_bool(self, inplace)
            if r is not None:
                return r
            return self._interpolate_with_fill(
                method=m,
                axis=axis,
                inplace=inplace,
                limit=limit,
                fill_value=fill_value,
                coerce=coerce,
                downcast=downcast,
            )
        # validate the interp method
        m = missing.clean_interp_method(method, **kwargs)

        r = check_int_bool(self, inplace)
        if r is not None:
            return r
        return self._interpolate(
            method=m,
            index=index,
            values=values,
            axis=axis,
            limit=limit,
            limit_direction=limit_direction,
            limit_area=limit_area,
            fill_value=fill_value,
            inplace=inplace,
            downcast=downcast,
            **kwargs,
        )

    def _interpolate_with_fill(
        self,
        method="pad",
        axis=0,
        inplace=False,
        limit=None,
        fill_value=None,
        coerce=False,
        downcast=None,
    ):
        """ fillna but using the interpolate machinery """
        inplace = validate_bool_kwarg(inplace, "inplace")

        # if we are coercing, then don't force the conversion
        # if the block can't hold the type
        if coerce:
            if not self._can_hold_na:
                if inplace:
                    return [self]
                else:
                    return [self.copy()]

        values = self.values if inplace else self.values.copy()

        # We only get here for non-ExtensionBlock
        fill_value = convert_scalar_for_putitemlike(fill_value, self.values.dtype)

        values = missing.interpolate_2d(
            values,
            method=method,
            axis=axis,
            limit=limit,
            fill_value=fill_value,
            dtype=self.dtype,
        )

        blocks = [self.make_block_same_class(values, ndim=self.ndim)]
        return self._maybe_downcast(blocks, downcast)

    def _interpolate(
        self,
        method=None,
        index=None,
        values=None,
        fill_value=None,
        axis=0,
        limit=None,
        limit_direction="forward",
        limit_area=None,
        inplace=False,
        downcast=None,
        **kwargs,
    ):
        """ interpolate using scipy wrappers """
        inplace = validate_bool_kwarg(inplace, "inplace")
        data = self.values if inplace else self.values.copy()

        # only deal with floats
        if not self.is_float:
            if not self.is_integer:
                return self
            data = data.astype(np.float64)

        if fill_value is None:
            fill_value = self.fill_value

        if method in ("krogh", "piecewise_polynomial", "pchip"):
            if not index.is_monotonic:
                raise ValueError(
                    f"{method} interpolation requires that the index be monotonic."
                )
        # process 1-d slices in the axis direction

        def func(x):

            # process a 1-d slice, returning it
            # should the axis argument be handled below in apply_along_axis?
            # i.e. not an arg to missing.interpolate_1d
            return missing.interpolate_1d(
                index,
                x,
                method=method,
                limit=limit,
                limit_direction=limit_direction,
                limit_area=limit_area,
                fill_value=fill_value,
                bounds_error=False,
                **kwargs,
            )

        # interp each column independently
        interp_values = np.apply_along_axis(func, axis, data)

        blocks = [self.make_block_same_class(interp_values)]
        return self._maybe_downcast(blocks, downcast)

    def take_nd(self, indexer, axis, new_mgr_locs=None, fill_tuple=None):
        """
        Take values according to indexer and return them as a block.bb

        """
        # algos.take_nd dispatches for DatetimeTZBlock, CategoricalBlock
        # so need to preserve types
        # sparse is treated like an ndarray, but needs .get_values() shaping

        values = self.values

        if fill_tuple is None:
            fill_value = self.fill_value
            allow_fill = False
        else:
            fill_value = fill_tuple[0]
            allow_fill = True

        new_values = algos.take_nd(
            values, indexer, axis=axis, allow_fill=allow_fill, fill_value=fill_value
        )

        # Called from three places in managers, all of which satisfy
        #  this assertion
        assert not (axis == 0 and new_mgr_locs is None)
        if new_mgr_locs is None:
            new_mgr_locs = self.mgr_locs

        if not is_dtype_equal(new_values.dtype, self.dtype):
            return self.make_block(new_values, new_mgr_locs)
        else:
            return self.make_block_same_class(new_values, new_mgr_locs)

    def diff(self, n: int, axis: int = 1) -> List["Block"]:
        """ return block for the diff of the values """
        new_values = algos.diff(self.values, n, axis=axis, stacklevel=7)
        # We use block_shape for ExtensionBlock subclasses, which may call here
        # via a super.
        new_values = _block_shape(new_values, ndim=self.ndim)
        return [self.make_block(values=new_values)]

    def shift(self, periods, axis=0, fill_value=None):
        """ shift the block by periods, possibly upcast """
        # convert integer to float if necessary. need to do a lot more than
        # that, handle boolean etc also
        new_values, fill_value = maybe_upcast(self.values, fill_value)

        # make sure array sent to np.roll is c_contiguous
        f_ordered = new_values.flags.f_contiguous
        if f_ordered:
            new_values = new_values.T
            axis = new_values.ndim - axis - 1

        if np.prod(new_values.shape):
            new_values = np.roll(new_values, ensure_platform_int(periods), axis=axis)

        axis_indexer = [slice(None)] * self.ndim
        if periods > 0:
            axis_indexer[axis] = slice(None, periods)
        else:
            axis_indexer[axis] = slice(periods, None)
        new_values[tuple(axis_indexer)] = fill_value

        # restore original order
        if f_ordered:
            new_values = new_values.T

        return [self.make_block(new_values)]

    def where(
        self,
        other,
        cond,
        align=True,
        errors="raise",
        try_cast: bool = False,
        axis: int = 0,
    ) -> List["Block"]:
        """
        evaluate the block; return result block(s) from the result

        Parameters
        ----------
        other : a ndarray/object
        cond  : the condition to respect
        align : boolean, perform alignment on other/cond
        errors : str, {'raise', 'ignore'}, default 'raise'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object
        axis : int

        Returns
        -------
        a new block(s), the result of the func
        """
        import pandas.core.computation.expressions as expressions

        assert errors in ["raise", "ignore"]
        transpose = self.ndim == 2

        values = self.values
        orig_other = other
        if transpose:
            values = values.T

        other = getattr(other, "_values", getattr(other, "values", other))
        cond = getattr(cond, "values", cond)

        # If the default broadcasting would go in the wrong direction, then
        # explicitly reshape other instead
        if getattr(other, "ndim", 0) >= 1:
            if values.ndim - 1 == other.ndim and axis == 1:
                other = other.reshape(tuple(other.shape + (1,)))
            elif transpose and values.ndim == self.ndim - 1:
                cond = cond.T

        if not hasattr(cond, "shape"):
            raise ValueError("where must have a condition that is ndarray like")

        def where_func(cond, values, other):

            if not (
                (self.is_integer or self.is_bool)
                and lib.is_float(other)
                and np.isnan(other)
            ):
                # np.where will cast integer array to floats in this case
                if not self._can_hold_element(other):
                    raise TypeError
                if lib.is_scalar(other) and isinstance(values, np.ndarray):
                    # convert datetime to datetime64, timedelta to timedelta64
                    other = convert_scalar_for_putitemlike(other, values.dtype)

            # By the time we get here, we should have all Series/Index
            #  args extracted to  ndarray
            fastres = expressions.where(cond, values, other)
            return fastres

        if cond.ravel().all():
            result = values
        else:
            # see if we can operate on the entire block, or need item-by-item
            # or if we are a single block (ndim == 1)
            try:
                result = where_func(cond, values, other)
            except TypeError:

                # we cannot coerce, return a compat dtype
                # we are explicitly ignoring errors
                block = self.coerce_to_target_dtype(other)
                blocks = block.where(
                    orig_other,
                    cond,
                    align=align,
                    errors=errors,
                    try_cast=try_cast,
                    axis=axis,
                )
                return self._maybe_downcast(blocks, "infer")

        if self._can_hold_na or self.ndim == 1:

            if transpose:
                result = result.T

            return [self.make_block(result)]

        # might need to separate out blocks
        axis = cond.ndim - 1
        cond = cond.swapaxes(axis, 0)
        mask = np.array([cond[i].all() for i in range(cond.shape[0])], dtype=bool)

        result_blocks = []
        for m in [mask, ~mask]:
            if m.any():
                taken = result.take(m.nonzero()[0], axis=axis)
                r = maybe_downcast_numeric(taken, self.dtype)
                nb = self.make_block(r.T, placement=self.mgr_locs[m])
                result_blocks.append(nb)

        return result_blocks

    def equals(self, other) -> bool:
        if self.dtype != other.dtype or self.shape != other.shape:
            return False
        return array_equivalent(self.values, other.values)

    def _unstack(self, unstacker_func, new_columns, n_rows, fill_value):
        """Return a list of unstacked blocks of self

        Parameters
        ----------
        unstacker_func : callable
            Partially applied unstacker.
        new_columns : Index
            All columns of the unstacked BlockManager.
        n_rows : int
            Only used in ExtensionBlock._unstack
        fill_value : int
            Only used in ExtensionBlock._unstack

        Returns
        -------
        blocks : list of Block
            New blocks of unstacked values.
        mask : array_like of bool
            The mask of columns of `blocks` we should keep.
        """
        unstacker = unstacker_func(self.values.T)
        new_items = unstacker.get_new_columns()
        new_placement = new_columns.get_indexer(new_items)
        new_values, mask = unstacker.get_new_values()

        mask = mask.any(0)
        new_values = new_values.T[mask]
        new_placement = new_placement[mask]

        blocks = [make_block(new_values, placement=new_placement)]
        return blocks, mask

    def quantile(self, qs, interpolation="linear", axis=0):
        """
        compute the quantiles of the

        Parameters
        ----------
        qs: a scalar or list of the quantiles to be computed
        interpolation: type of interpolation, default 'linear'
        axis: axis to compute, default 0

        Returns
        -------
        Block
        """
        # We should always have ndim == 2 because Series dispatches to DataFrame
        assert self.ndim == 2

        values = self.get_values()

        is_empty = values.shape[axis] == 0
        orig_scalar = not is_list_like(qs)
        if orig_scalar:
            # make list-like, unpack later
            qs = [qs]

        if is_empty:
            # create the array of na_values
            # 2d len(values) * len(qs)
            result = np.repeat(
                np.array([self.fill_value] * len(qs)), len(values)
            ).reshape(len(values), len(qs))
        else:
            # asarray needed for Sparse, see GH#24600
            mask = np.asarray(isna(values))
            result = nanpercentile(
                values,
                np.array(qs) * 100,
                axis=axis,
                na_value=self.fill_value,
                mask=mask,
                ndim=values.ndim,
                interpolation=interpolation,
            )

            result = np.array(result, copy=False)
            result = result.T

        if orig_scalar and not lib.is_scalar(result):
            # result could be scalar in case with is_empty and self.ndim == 1
            assert result.shape[-1] == 1, result.shape
            result = result[..., 0]
            result = lib.item_from_zerodim(result)

        ndim = np.ndim(result)
        return make_block(result, placement=np.arange(len(result)), ndim=ndim)

    def _replace_coerce(
        self, to_replace, value, inplace=True, regex=False, convert=False, mask=None
    ):
        """
        Replace value corresponding to the given boolean array with another
        value.

        Parameters
        ----------
        to_replace : object or pattern
            Scalar to replace or regular expression to match.
        value : object
            Replacement object.
        inplace : bool, default False
            Perform inplace modification.
        regex : bool, default False
            If true, perform regular expression substitution.
        convert : bool, default True
            If true, try to coerce any object types to better types.
        mask : array-like of bool, optional
            True indicate corresponding element is ignored.

        Returns
        -------
        A new block if there is anything to replace or the original block.
        """
        if mask.any():
            if not regex:
                self = self.coerce_to_target_dtype(value)
                return self.putmask(mask, value, inplace=inplace)
            else:
                return self._replace_single(
                    to_replace,
                    value,
                    inplace=inplace,
                    regex=regex,
                    convert=convert,
                    mask=mask,
                )
        return self


class NonConsolidatableMixIn:
    """ hold methods for the nonconsolidatable blocks """

    _can_consolidate = False
    _verify_integrity = False
    _validate_ndim = False

    def __init__(self, values, placement, ndim=None):
        """Initialize a non-consolidatable block.

        'ndim' may be inferred from 'placement'.

        This will call continue to call __init__ for the other base
        classes mixed in with this Mixin.
        """
        # Placement must be converted to BlockPlacement so that we can check
        # its length
        if not isinstance(placement, libinternals.BlockPlacement):
            placement = libinternals.BlockPlacement(placement)

        # Maybe infer ndim from placement
        if ndim is None:
            if len(placement) != 1:
                ndim = 1
            else:
                ndim = 2
        super().__init__(values, placement, ndim=ndim)

    @property
    def shape(self):
        if self.ndim == 1:
            return ((len(self.values)),)
        return (len(self.mgr_locs), len(self.values))

    def iget(self, col):

        if self.ndim == 2 and isinstance(col, tuple):
            col, loc = col
            if not com.is_null_slice(col) and col != 0:
                raise IndexError(f"{self} only contains one item")
            elif isinstance(col, slice):
                if col != slice(None):
                    raise NotImplementedError(col)
                return self.values[[loc]]
            return self.values[loc]
        else:
            if col != 0:
                raise IndexError(f"{self} only contains one item")
            return self.values

    def should_store(self, value):
        return isinstance(value, self._holder)

    def set(self, locs, values, check=False):
        assert locs.tolist() == [0]
        self.values = values

    def putmask(self, mask, new, align=True, inplace=False, axis=0, transpose=False):
        """
        putmask the data to the block; we must be a single block and not
        generate other blocks

        return the resulting block

        Parameters
        ----------
        mask  : the condition to respect
        new : a ndarray/object
        align : boolean, perform alignment on other/cond, default is True
        inplace : perform inplace modification, default is False

        Returns
        -------
        a new block, the result of the putmask
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        # use block's copy logic.
        # .values may be an Index which does shallow copy by default
        new_values = self.values if inplace else self.copy().values

        if isinstance(new, np.ndarray) and len(new) == len(mask):
            new = new[mask]

        mask = _safe_reshape(mask, new_values.shape)

        new_values[mask] = new
        return [self.make_block(values=new_values)]

    def _get_unstack_items(self, unstacker, new_columns):
        """
        Get the placement, values, and mask for a Block unstack.

        This is shared between ObjectBlock and ExtensionBlock. They
        differ in that ObjectBlock passes the values, while ExtensionBlock
        passes the dummy ndarray of positions to be used by a take
        later.

        Parameters
        ----------
        unstacker : pandas.core.reshape.reshape._Unstacker
        new_columns : Index
            All columns of the unstacked BlockManager.

        Returns
        -------
        new_placement : ndarray[int]
            The placement of the new columns in `new_columns`.
        new_values : Union[ndarray, ExtensionArray]
            The first return value from _Unstacker.get_new_values.
        mask : ndarray[bool]
            The second return value from _Unstacker.get_new_values.
        """
        # shared with ExtensionBlock
        new_items = unstacker.get_new_columns()
        new_placement = new_columns.get_indexer(new_items)
        new_values, mask = unstacker.get_new_values()

        mask = mask.any(0)
        return new_placement, new_values, mask


class ExtensionBlock(NonConsolidatableMixIn, Block):
    """Block for holding extension types.

    Notes
    -----
    This holds all 3rd-party extension array types. It's also the immediate
    parent class for our internal extension types' blocks, CategoricalBlock.

    ExtensionArrays are limited to 1-D.
    """

    is_extension = True

    def __init__(self, values, placement, ndim=None):
        values = self._maybe_coerce_values(values)
        super().__init__(values, placement, ndim)

    def _maybe_coerce_values(self, values):
        """
        Unbox to an extension array.

        This will unbox an ExtensionArray stored in an Index or Series.
        ExtensionArrays pass through. No dtype coercion is done.

        Parameters
        ----------
        values : Index, Series, ExtensionArray

        Returns
        -------
        ExtensionArray
        """
        return extract_array(values)

    @property
    def _holder(self):
        # For extension blocks, the holder is values-dependent.
        return type(self.values)

    @property
    def fill_value(self):
        # Used in reindex_indexer
        return self.values.dtype.na_value

    @property
    def _can_hold_na(self):
        # The default ExtensionArray._can_hold_na is True
        return self._holder._can_hold_na

    @property
    def is_view(self):
        """Extension arrays are never treated as views."""
        return False

    @property
    def is_numeric(self):
        return self.values.dtype._is_numeric

    def setitem(self, indexer, value):
        """Set the value inplace, returning a same-typed block.

        This differs from Block.setitem by not allowing setitem to change
        the dtype of the Block.

        Parameters
        ----------
        indexer : tuple, list-like, array-like, slice
            The subset of self.values to set
        value : object
            The value being set

        Returns
        -------
        Block

        Notes
        -----
        `indexer` is a direct slice/positional indexer. `value` must
        be a compatible shape.
        """
        if isinstance(indexer, tuple):
            # we are always 1-D
            indexer = indexer[0]

        check_setitem_lengths(indexer, value, self.values)
        self.values[indexer] = value
        return self

    def get_values(self, dtype=None):
        # ExtensionArrays must be iterable, so this works.
        values = np.asarray(self.values)
        if values.ndim == self.ndim - 1:
            values = values.reshape((1,) + values.shape)
        return values

    def array_values(self) -> ExtensionArray:
        return self.values

    def to_dense(self):
        return np.asarray(self.values)

    def to_native_types(self, slicer=None, na_rep="nan", quoting=None, **kwargs):
        """override to use ExtensionArray astype for the conversion"""
        values = self.values
        if slicer is not None:
            values = values[slicer]
        mask = isna(values)

        values = np.asarray(values.astype(object))
        values[mask] = na_rep

        # we are expected to return a 2-d ndarray
        return values.reshape(1, len(values))

    def take_nd(self, indexer, axis=0, new_mgr_locs=None, fill_tuple=None):
        """
        Take values according to indexer and return them as a block.
        """
        if fill_tuple is None:
            fill_value = None
        else:
            fill_value = fill_tuple[0]

        # axis doesn't matter; we are really a single-dim object
        # but are passed the axis depending on the calling routing
        # if its REALLY axis 0, then this will be a reindex and not a take
        new_values = self.values.take(indexer, fill_value=fill_value, allow_fill=True)

        # Called from three places in managers, all of which satisfy
        #  this assertion
        assert not (self.ndim == 1 and new_mgr_locs is None)
        if new_mgr_locs is None:
            new_mgr_locs = self.mgr_locs

        return self.make_block_same_class(new_values, new_mgr_locs)

    def _can_hold_element(self, element: Any) -> bool:
        # XXX: We may need to think about pushing this onto the array.
        # We're doing the same as CategoricalBlock here.
        return True

    def _slice(self, slicer):
        """ return a slice of my values """
        # slice the category
        # return same dims as we currently have

        if isinstance(slicer, tuple) and len(slicer) == 2:
            if not com.is_null_slice(slicer[0]):
                raise AssertionError("invalid slicing for a 1-ndim categorical")
            slicer = slicer[1]

        return self.values[slicer]

    def concat_same_type(self, to_concat, placement=None):
        """
        Concatenate list of single blocks of the same type.
        """
        values = self._holder._concat_same_type([blk.values for blk in to_concat])
        placement = placement or slice(0, len(values), 1)
        return self.make_block_same_class(values, ndim=self.ndim, placement=placement)

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        values = self.values if inplace else self.values.copy()
        values = values.fillna(value=value, limit=limit)
        return [
            self.make_block_same_class(
                values=values, placement=self.mgr_locs, ndim=self.ndim
            )
        ]

    def interpolate(
        self, method="pad", axis=0, inplace=False, limit=None, fill_value=None, **kwargs
    ):

        values = self.values if inplace else self.values.copy()
        return self.make_block_same_class(
            values=values.fillna(value=fill_value, method=method, limit=limit),
            placement=self.mgr_locs,
        )

    def diff(self, n: int, axis: int = 1) -> List["Block"]:
        if axis == 1:
            # we are by definition 1D.
            axis = 0
        return super().diff(n, axis)

    def shift(
        self,
        periods: int,
        axis: libinternals.BlockPlacement = 0,
        fill_value: Any = None,
    ) -> List["ExtensionBlock"]:
        """
        Shift the block by `periods`.

        Dispatches to underlying ExtensionArray and re-boxes in an
        ExtensionBlock.
        """
        return [
            self.make_block_same_class(
                self.values.shift(periods=periods, fill_value=fill_value),
                placement=self.mgr_locs,
                ndim=self.ndim,
            )
        ]

    def where(
        self,
        other,
        cond,
        align=True,
        errors="raise",
        try_cast: bool = False,
        axis: int = 0,
    ) -> List["Block"]:
        if isinstance(other, ABCDataFrame):
            # ExtensionArrays are 1-D, so if we get here then
            # `other` should be a DataFrame with a single column.
            assert other.shape[1] == 1
            other = other.iloc[:, 0]

        other = extract_array(other, extract_numpy=True)

        if isinstance(cond, ABCDataFrame):
            assert cond.shape[1] == 1
            cond = cond.iloc[:, 0]

        cond = extract_array(cond, extract_numpy=True)

        if lib.is_scalar(other) and isna(other):
            # The default `other` for Series / Frame is np.nan
            # we want to replace that with the correct NA value
            # for the type
            other = self.dtype.na_value

        if is_sparse(self.values):
            # TODO(SparseArray.__setitem__): remove this if condition
            # We need to re-infer the type of the data after doing the
            # where, for cases where the subtypes don't match
            dtype = None
        else:
            dtype = self.dtype

        result = self.values.copy()
        icond = ~cond
        if lib.is_scalar(other):
            set_other = other
        else:
            set_other = other[icond]
        try:
            result[icond] = set_other
        except (NotImplementedError, TypeError):
            # NotImplementedError for class not implementing `__setitem__`
            # TypeError for SparseArray, which implements just to raise
            # a TypeError
            result = self._holder._from_sequence(
                np.where(cond, self.values, other), dtype=dtype
            )

        return [self.make_block_same_class(result, placement=self.mgr_locs)]

    @property
    def _ftype(self):
        return getattr(self.values, "_pandas_ftype", Block._ftype)

    def _unstack(self, unstacker_func, new_columns, n_rows, fill_value):
        # ExtensionArray-safe unstack.
        # We override ObjectBlock._unstack, which unstacks directly on the
        # values of the array. For EA-backed blocks, this would require
        # converting to a 2-D ndarray of objects.
        # Instead, we unstack an ndarray of integer positions, followed by
        # a `take` on the actual values.
        dummy_arr = np.arange(n_rows)
        dummy_unstacker = functools.partial(unstacker_func, fill_value=-1)
        unstacker = dummy_unstacker(dummy_arr)

        new_placement, new_values, mask = self._get_unstack_items(
            unstacker, new_columns
        )

        blocks = [
            self.make_block_same_class(
                self.values.take(indices, allow_fill=True, fill_value=fill_value),
                [place],
            )
            for indices, place in zip(new_values.T, new_placement)
        ]
        return blocks, mask


class ObjectValuesExtensionBlock(ExtensionBlock):
    """
    Block providing backwards-compatibility for `.values`.

    Used by PeriodArray and IntervalArray to ensure that
    Series[T].values is an ndarray of objects.
    """

    def external_values(self):
        return self.values.astype(object)


class NumericBlock(Block):
    __slots__ = ()
    is_numeric = True
    _can_hold_na = True


class FloatOrComplexBlock(NumericBlock):
    __slots__ = ()

    def equals(self, other) -> bool:
        if self.dtype != other.dtype or self.shape != other.shape:
            return False
        left, right = self.values, other.values
        return ((left == right) | (np.isnan(left) & np.isnan(right))).all()


class FloatBlock(FloatOrComplexBlock):
    __slots__ = ()
    is_float = True

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, (np.floating, np.integer)) and not issubclass(
                tipo.type, (np.datetime64, np.timedelta64)
            )
        return isinstance(
            element, (float, int, np.floating, np.int_)
        ) and not isinstance(
            element,
            (bool, np.bool_, datetime, timedelta, np.datetime64, np.timedelta64),
        )

    def to_native_types(
        self,
        slicer=None,
        na_rep="",
        float_format=None,
        decimal=".",
        quoting=None,
        **kwargs,
    ):
        """ convert to our native types format, slicing if desired """
        values = self.values
        if slicer is not None:
            values = values[:, slicer]

        # see gh-13418: no special formatting is desired at the
        # output (important for appropriate 'quoting' behaviour),
        # so do not pass it through the FloatArrayFormatter
        if float_format is None and decimal == ".":
            mask = isna(values)

            if not quoting:
                values = values.astype(str)
            else:
                values = np.array(values, dtype="object")

            values[mask] = na_rep
            return values

        from pandas.io.formats.format import FloatArrayFormatter

        formatter = FloatArrayFormatter(
            values,
            na_rep=na_rep,
            float_format=float_format,
            decimal=decimal,
            quoting=quoting,
            fixed_width=False,
        )
        return formatter.get_result_as_array()

    def should_store(self, value):
        # when inserting a column should not coerce integers to floats
        # unnecessarily
        return issubclass(value.dtype.type, np.floating) and value.dtype == self.dtype


class ComplexBlock(FloatOrComplexBlock):
    __slots__ = ()
    is_complex = True

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, (np.floating, np.integer, np.complexfloating))
        return isinstance(
            element, (float, int, complex, np.float_, np.int_)
        ) and not isinstance(element, (bool, np.bool_))

    def should_store(self, value):
        return issubclass(value.dtype.type, np.complexfloating)


class IntBlock(NumericBlock):
    __slots__ = ()
    is_integer = True
    _can_hold_na = False

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return (
                issubclass(tipo.type, np.integer)
                and not issubclass(tipo.type, (np.datetime64, np.timedelta64))
                and self.dtype.itemsize >= tipo.itemsize
            )
        return is_integer(element)

    def should_store(self, value):
        return is_integer_dtype(value) and value.dtype == self.dtype


class DatetimeLikeBlockMixin:
    """Mixin class for DatetimeBlock, DatetimeTZBlock, and TimedeltaBlock."""

    @property
    def _holder(self):
        return DatetimeArray

    @property
    def fill_value(self):
        return np.datetime64("NaT", "ns")

    def get_values(self, dtype=None):
        """
        return object dtype as boxed values, such as Timestamps/Timedelta
        """
        if is_object_dtype(dtype):
            values = self.values.ravel()
            result = self._holder(values).astype(object)
            return result.reshape(self.values.shape)
        return self.values

    def internal_values(self):
        # Override to return DatetimeArray and TimedeltaArray
        return self.array_values()

    def iget(self, key):
        # GH#31649 we need to wrap scalars in Timestamp/Timedelta
        # TODO: this can be removed if we ever have 2D EA
        result = super().iget(key)
        if isinstance(result, np.datetime64):
            result = Timestamp(result)
        elif isinstance(result, np.timedelta64):
            result = Timedelta(result)
        return result


class DatetimeBlock(DatetimeLikeBlockMixin, Block):
    __slots__ = ()
    is_datetime = True

    def __init__(self, values, placement, ndim=None):
        values = self._maybe_coerce_values(values)
        super().__init__(values, placement=placement, ndim=ndim)

    @property
    def _can_hold_na(self):
        return True

    def _maybe_coerce_values(self, values):
        """
        Input validation for values passed to __init__. Ensure that
        we have datetime64ns, coercing if necessary.

        Parameters
        ----------
        values : array-like
            Must be convertible to datetime64

        Returns
        -------
        values : ndarray[datetime64ns]

        Overridden by DatetimeTZBlock.
        """
        if values.dtype != _NS_DTYPE:
            values = conversion.ensure_datetime64ns(values)

        if isinstance(values, DatetimeArray):
            values = values._data

        assert isinstance(values, np.ndarray), type(values)
        return values

    def astype(self, dtype, copy: bool = False, errors: str = "raise"):
        """
        these automatically copy, so copy=True has no effect
        raise on an except if raise == True
        """
        dtype = pandas_dtype(dtype)

        # if we are passed a datetime64[ns, tz]
        if is_datetime64tz_dtype(dtype):
            values = self.values
            if getattr(values, "tz", None) is None:
                values = DatetimeArray(values).tz_localize("UTC")
            values = values.tz_convert(dtype.tz)
            return self.make_block(values)

        # delegate
        return super().astype(dtype=dtype, copy=copy, errors=errors)

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            if self.is_datetimetz:
                # require exact match, since non-nano does not exist
                return is_dtype_equal(tipo, self.dtype) or is_valid_nat_for_dtype(
                    element, self.dtype
                )

            # GH#27419 if we get a non-nano datetime64 object
            return is_datetime64_dtype(tipo)
        elif element is NaT:
            return True
        elif isinstance(element, datetime):
            if self.is_datetimetz:
                return tz_compare(element.tzinfo, self.dtype.tz)
            return element.tzinfo is None

        return is_valid_nat_for_dtype(element, self.dtype)

    def to_native_types(
        self, slicer=None, na_rep=None, date_format=None, quoting=None, **kwargs
    ):
        """ convert to our native types format, slicing if desired """
        values = self.values
        i8values = self.values.view("i8")

        if slicer is not None:
            values = values[..., slicer]
            i8values = i8values[..., slicer]

        from pandas.io.formats.format import _get_format_datetime64_from_values

        fmt = _get_format_datetime64_from_values(values, date_format)

        result = tslib.format_array_from_datetime(
            i8values.ravel(),
            tz=getattr(self.values, "tz", None),
            format=fmt,
            na_rep=na_rep,
        ).reshape(i8values.shape)
        return np.atleast_2d(result)

    def should_store(self, value):
        return (
            issubclass(value.dtype.type, np.datetime64)
            and not is_datetime64tz_dtype(value)
            and not is_extension_array_dtype(value)
        )

    def set(self, locs, values):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        values = conversion.ensure_datetime64ns(values, copy=False)

        self.values[locs] = values

    def external_values(self):
        return np.asarray(self.values.astype("datetime64[ns]", copy=False))

    def array_values(self) -> ExtensionArray:
        return DatetimeArray._simple_new(self.values)


class DatetimeTZBlock(ExtensionBlock, DatetimeBlock):
    """ implement a datetime64 block with a tz attribute """

    __slots__ = ()
    is_datetimetz = True
    is_extension = True

    internal_values = Block.internal_values
    _can_hold_element = DatetimeBlock._can_hold_element
    to_native_types = DatetimeBlock.to_native_types
    fill_value = np.datetime64("NaT", "ns")

    @property
    def _holder(self):
        return DatetimeArray

    def _maybe_coerce_values(self, values):
        """Input validation for values passed to __init__. Ensure that
        we have datetime64TZ, coercing if necessary.

        Parameters
        ----------
        values : array-like
            Must be convertible to datetime64

        Returns
        -------
        values : DatetimeArray
        """
        if not isinstance(values, self._holder):
            values = self._holder(values)

        if values.tz is None:
            raise ValueError("cannot create a DatetimeTZBlock without a tz")

        return values

    @property
    def is_view(self):
        """ return a boolean if I am possibly a view """
        # check the ndarray values of the DatetimeIndex values
        return self.values._data.base is not None

    def get_values(self, dtype=None):
        """
        Returns an ndarray of values.

        Parameters
        ----------
        dtype : np.dtype
            Only `object`-like dtypes are respected here (not sure
            why).

        Returns
        -------
        values : ndarray
            When ``dtype=object``, then and object-dtype ndarray of
            boxed values is returned. Otherwise, an M8[ns] ndarray
            is returned.

            DatetimeArray is always 1-d. ``get_values`` will reshape
            the return value to be the same dimensionality as the
            block.
        """
        values = self.values
        if is_object_dtype(dtype):
            values = values.astype(object)

        values = np.asarray(values)

        if self.ndim == 2:
            # Ensure that our shape is correct for DataFrame.
            # ExtensionArrays are always 1-D, even in a DataFrame when
            # the analogous NumPy-backed column would be a 2-D ndarray.
            values = values.reshape(1, -1)
        return values

    def to_dense(self):
        # we request M8[ns] dtype here, even though it discards tzinfo,
        # as lots of code (e.g. anything using values_from_object)
        # expects that behavior.
        return np.asarray(self.values, dtype=_NS_DTYPE)

    def _slice(self, slicer):
        """ return a slice of my values """
        if isinstance(slicer, tuple):
            col, loc = slicer
            if not com.is_null_slice(col) and col != 0:
                raise IndexError(f"{self} only contains one item")
            return self.values[loc]
        return self.values[slicer]

    def diff(self, n: int, axis: int = 0) -> List["Block"]:
        """
        1st discrete difference.

        Parameters
        ----------
        n : int
            Number of periods to diff.
        axis : int, default 0
            Axis to diff upon.

        Returns
        -------
        A list with a new TimeDeltaBlock.

        Notes
        -----
        The arguments here are mimicking shift so they are called correctly
        by apply.
        """
        if axis == 0:
            # Cannot currently calculate diff across multiple blocks since this
            # function is invoked via apply
            raise NotImplementedError
        new_values = (self.values - self.shift(n, axis=axis)[0].values).asi8

        # Reshape the new_values like how algos.diff does for timedelta data
        new_values = new_values.reshape(1, len(new_values))
        new_values = new_values.astype("timedelta64[ns]")
        return [TimeDeltaBlock(new_values, placement=self.mgr_locs.indexer)]

    def concat_same_type(self, to_concat, placement=None):
        # need to handle concat([tz1, tz2]) here, since DatetimeArray
        # only handles cases where all the tzs are the same.
        # Instead of placing the condition here, it could also go into the
        # is_uniform_join_units check, but I'm not sure what is better.
        if len({x.dtype for x in to_concat}) > 1:
            values = concat_datetime([x.values for x in to_concat])
            placement = placement or slice(0, len(values), 1)

            if self.ndim > 1:
                values = np.atleast_2d(values)
            return ObjectBlock(values, ndim=self.ndim, placement=placement)
        return super().concat_same_type(to_concat, placement)

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        # We support filling a DatetimeTZ with a `value` whose timezone
        # is different by coercing to object.
        if self._can_hold_element(value):
            return super().fillna(value, limit, inplace, downcast)

        # different timezones, or a non-tz
        return self.astype(object).fillna(
            value, limit=limit, inplace=inplace, downcast=downcast
        )

    def setitem(self, indexer, value):
        # https://github.com/pandas-dev/pandas/issues/24020
        # Need a dedicated setitem until #24020 (type promotion in setitem
        # for extension arrays) is designed and implemented.
        if self._can_hold_element(value) or (
            isinstance(indexer, np.ndarray) and indexer.size == 0
        ):
            return super().setitem(indexer, value)

        obj_vals = self.values.astype(object)
        newb = make_block(
            obj_vals, placement=self.mgr_locs, klass=ObjectBlock, ndim=self.ndim
        )
        return newb.setitem(indexer, value)

    def equals(self, other) -> bool:
        # override for significant performance improvement
        if self.dtype != other.dtype or self.shape != other.shape:
            return False
        return (self.values.view("i8") == other.values.view("i8")).all()

    def quantile(self, qs, interpolation="linear", axis=0):
        naive = self.values.view("M8[ns]")

        # kludge for 2D block with 1D values
        naive = naive.reshape(self.shape)

        blk = self.make_block(naive)
        res_blk = blk.quantile(qs, interpolation=interpolation, axis=axis)

        # ravel is kludge for 2D block with 1D values, assumes column-like
        aware = self._holder(res_blk.values.ravel(), dtype=self.dtype)
        return self.make_block_same_class(aware, ndim=res_blk.ndim)


class TimeDeltaBlock(DatetimeLikeBlockMixin, IntBlock):
    __slots__ = ()
    is_timedelta = True
    _can_hold_na = True
    is_numeric = False
    fill_value = np.timedelta64("NaT", "ns")

    def __init__(self, values, placement, ndim=None):
        if values.dtype != _TD_DTYPE:
            values = conversion.ensure_timedelta64ns(values)
        if isinstance(values, TimedeltaArray):
            values = values._data
        assert isinstance(values, np.ndarray), type(values)
        super().__init__(values, placement=placement, ndim=ndim)

    @property
    def _holder(self):
        return TimedeltaArray

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, np.timedelta64)
        elif element is NaT:
            return True
        elif isinstance(element, (timedelta, np.timedelta64)):
            return True
        return is_valid_nat_for_dtype(element, self.dtype)

    def fillna(self, value, **kwargs):

        # allow filling with integers to be
        # interpreted as nanoseconds
        if is_integer(value):
            # Deprecation GH#24694, GH#19233
            raise TypeError(
                "Passing integers to fillna for timedelta64[ns] dtype is no "
                "longer supported.  To obtain the old behavior, pass "
                "`pd.Timedelta(seconds=n)` instead."
            )
        return super().fillna(value, **kwargs)

    def should_store(self, value):
        return issubclass(
            value.dtype.type, np.timedelta64
        ) and not is_extension_array_dtype(value)

    def to_native_types(self, slicer=None, na_rep=None, quoting=None, **kwargs):
        """ convert to our native types format, slicing if desired """
        values = self.values
        if slicer is not None:
            values = values[:, slicer]
        mask = isna(values)

        rvalues = np.empty(values.shape, dtype=object)
        if na_rep is None:
            na_rep = "NaT"
        rvalues[mask] = na_rep
        imask = (~mask).ravel()

        # FIXME:
        # should use the formats.format.Timedelta64Formatter here
        # to figure what format to pass to the Timedelta
        # e.g. to not show the decimals say
        rvalues.flat[imask] = np.array(
            [Timedelta(val)._repr_base(format="all") for val in values.ravel()[imask]],
            dtype=object,
        )
        return rvalues

    def external_values(self):
        return np.asarray(self.values.astype("timedelta64[ns]", copy=False))

    def array_values(self) -> ExtensionArray:
        return TimedeltaArray._simple_new(self.values)


class BoolBlock(NumericBlock):
    __slots__ = ()
    is_bool = True
    _can_hold_na = False

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, np.bool_)
        return isinstance(element, (bool, np.bool_))

    def should_store(self, value):
        return issubclass(value.dtype.type, np.bool_) and not is_extension_array_dtype(
            value
        )

    def replace(
        self, to_replace, value, inplace=False, filter=None, regex=False, convert=True
    ):
        inplace = validate_bool_kwarg(inplace, "inplace")
        to_replace_values = np.atleast_1d(to_replace)
        if not np.can_cast(to_replace_values, bool):
            return self
        return super().replace(
            to_replace,
            value,
            inplace=inplace,
            filter=filter,
            regex=regex,
            convert=convert,
        )


class ObjectBlock(Block):
    __slots__ = ()
    is_object = True
    _can_hold_na = True

    def __init__(self, values, placement=None, ndim=2):
        if issubclass(values.dtype.type, str):
            values = np.array(values, dtype=object)

        super().__init__(values, ndim=ndim, placement=placement)

    @property
    def is_bool(self):
        """ we can be a bool if we have only bool values but are of type
        object
        """
        return lib.is_bool_array(self.values.ravel())

    def convert(
        self,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
        coerce: bool = False,
    ):
        """ attempt to coerce any object types to better types return a copy of
        the block (if copy = True) by definition we ARE an ObjectBlock!!!!!

        can return multiple blocks!
        """
        # operate column-by-column
        def f(mask, val, idx):
            shape = val.shape
            values = soft_convert_objects(
                val.ravel(),
                datetime=datetime,
                numeric=numeric,
                timedelta=timedelta,
                coerce=coerce,
                copy=copy,
            )
            if isinstance(values, np.ndarray):
                # TODO: allow EA once reshape is supported
                values = values.reshape(shape)

            values = _block_shape(values, ndim=self.ndim)
            return values

        if self.ndim == 2:
            blocks = self.split_and_operate(None, f, False)
        else:
            values = f(None, self.values.ravel(), None)
            blocks = [make_block(values, ndim=self.ndim, placement=self.mgr_locs)]

        return blocks

    def _maybe_downcast(self, blocks: List["Block"], downcast=None) -> List["Block"]:

        if downcast is not None:
            return blocks

        # split and convert the blocks
        return _extend_blocks([b.convert(datetime=True, numeric=False) for b in blocks])

    def _can_hold_element(self, element: Any) -> bool:
        return True

    def should_store(self, value):
        return not (
            issubclass(
                value.dtype.type,
                (np.integer, np.floating, np.complexfloating, np.datetime64, np.bool_),
            )
            or is_extension_array_dtype(value)
        )

    def replace(
        self, to_replace, value, inplace=False, filter=None, regex=False, convert=True
    ):
        to_rep_is_list = is_list_like(to_replace)
        value_is_list = is_list_like(value)
        both_lists = to_rep_is_list and value_is_list
        either_list = to_rep_is_list or value_is_list

        result_blocks = []
        blocks = [self]

        if not either_list and is_re(to_replace):
            return self._replace_single(
                to_replace,
                value,
                inplace=inplace,
                filter=filter,
                regex=True,
                convert=convert,
            )
        elif not (either_list or regex):
            return super().replace(
                to_replace,
                value,
                inplace=inplace,
                filter=filter,
                regex=regex,
                convert=convert,
            )
        elif both_lists:
            for to_rep, v in zip(to_replace, value):
                result_blocks = []
                for b in blocks:
                    result = b._replace_single(
                        to_rep,
                        v,
                        inplace=inplace,
                        filter=filter,
                        regex=regex,
                        convert=convert,
                    )
                    result_blocks = _extend_blocks(result, result_blocks)
                blocks = result_blocks
            return result_blocks

        elif to_rep_is_list and regex:
            for to_rep in to_replace:
                result_blocks = []
                for b in blocks:
                    result = b._replace_single(
                        to_rep,
                        value,
                        inplace=inplace,
                        filter=filter,
                        regex=regex,
                        convert=convert,
                    )
                    result_blocks = _extend_blocks(result, result_blocks)
                blocks = result_blocks
            return result_blocks

        return self._replace_single(
            to_replace,
            value,
            inplace=inplace,
            filter=filter,
            convert=convert,
            regex=regex,
        )

    def _replace_single(
        self,
        to_replace,
        value,
        inplace=False,
        filter=None,
        regex=False,
        convert=True,
        mask=None,
    ):
        """
        Replace elements by the given value.

        Parameters
        ----------
        to_replace : object or pattern
            Scalar to replace or regular expression to match.
        value : object
            Replacement object.
        inplace : bool, default False
            Perform inplace modification.
        filter : list, optional
        regex : bool, default False
            If true, perform regular expression substitution.
        convert : bool, default True
            If true, try to coerce any object types to better types.
        mask : array-like of bool, optional
            True indicate corresponding element is ignored.

        Returns
        -------
        a new block, the result after replacing
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        # to_replace is regex compilable
        to_rep_re = regex and is_re_compilable(to_replace)

        # regex is regex compilable
        regex_re = is_re_compilable(regex)

        # only one will survive
        if to_rep_re and regex_re:
            raise AssertionError(
                "only one of to_replace and regex can be regex compilable"
            )

        # if regex was passed as something that can be a regex (rather than a
        # boolean)
        if regex_re:
            to_replace = regex

        regex = regex_re or to_rep_re

        # try to get the pattern attribute (compiled re) or it's a string
        if is_re(to_replace):
            pattern = to_replace.pattern
        else:
            pattern = to_replace

        # if the pattern is not empty and to_replace is either a string or a
        # regex
        if regex and pattern:
            rx = re.compile(to_replace)
        else:
            # if the thing to replace is not a string or compiled regex call
            # the superclass method -> to_replace is some kind of object
            return super().replace(
                to_replace, value, inplace=inplace, filter=filter, regex=regex
            )

        new_values = self.values if inplace else self.values.copy()

        # deal with replacing values with objects (strings) that match but
        # whose replacement is not a string (numeric, nan, object)
        if isna(value) or not isinstance(value, str):

            def re_replacer(s):
                if is_re(rx) and isinstance(s, str):
                    return value if rx.search(s) is not None else s
                else:
                    return s

        else:
            # value is guaranteed to be a string here, s can be either a string
            # or null if it's null it gets returned
            def re_replacer(s):
                if is_re(rx) and isinstance(s, str):
                    return rx.sub(value, s)
                else:
                    return s

        f = np.vectorize(re_replacer, otypes=[self.dtype])

        if filter is None:
            filt = slice(None)
        else:
            filt = self.mgr_locs.isin(filter).nonzero()[0]

        if mask is None:
            new_values[filt] = f(new_values[filt])
        else:
            new_values[filt][mask] = f(new_values[filt][mask])

        # convert
        block = self.make_block(new_values)
        if convert:
            block = block.convert(numeric=False)
        return block

    def _replace_coerce(
        self, to_replace, value, inplace=True, regex=False, convert=False, mask=None
    ):
        """
        Replace value corresponding to the given boolean array with another
        value.

        Parameters
        ----------
        to_replace : object or pattern
            Scalar to replace or regular expression to match.
        value : object
            Replacement object.
        inplace : bool, default False
            Perform inplace modification.
        regex : bool, default False
            If true, perform regular expression substitution.
        convert : bool, default True
            If true, try to coerce any object types to better types.
        mask : array-like of bool, optional
            True indicate corresponding element is ignored.

        Returns
        -------
        A new block if there is anything to replace or the original block.
        """
        if mask.any():
            block = super()._replace_coerce(
                to_replace=to_replace,
                value=value,
                inplace=inplace,
                regex=regex,
                convert=convert,
                mask=mask,
            )
            if convert:
                block = [b.convert(numeric=False, copy=True) for b in block]
            return block
        if convert:
            return [self.convert(numeric=False, copy=True)]
        return self


class CategoricalBlock(ExtensionBlock):
    __slots__ = ()
    is_categorical = True
    _verify_integrity = True
    _can_hold_na = True
    _concatenator = staticmethod(concat_categorical)

    def __init__(self, values, placement, ndim=None):
        # coerce to categorical if we can
        values = extract_array(values)
        assert isinstance(values, Categorical), type(values)
        super().__init__(values, placement=placement, ndim=ndim)

    @property
    def _holder(self):
        return Categorical

    @property
    def array_dtype(self):
        """ the dtype to return if I want to construct this block as an
        array
        """
        return np.object_

    def to_dense(self):
        # Categorical.get_values returns a DatetimeIndex for datetime
        # categories, so we can't simply use `np.asarray(self.values)` like
        # other types.
        return self.values._internal_get_values()

    def to_native_types(self, slicer=None, na_rep="", quoting=None, **kwargs):
        """ convert to our native types format, slicing if desired """
        values = self.values
        if slicer is not None:
            # Categorical is always one dimension
            values = values[slicer]
        mask = isna(values)
        values = np.array(values, dtype="object")
        values[mask] = na_rep

        # we are expected to return a 2-d ndarray
        return values.reshape(1, len(values))

    def concat_same_type(self, to_concat, placement=None):
        """
        Concatenate list of single blocks of the same type.

        Note that this CategoricalBlock._concat_same_type *may* not
        return a CategoricalBlock. When the categories in `to_concat`
        differ, this will return an object ndarray.

        If / when we decide we don't like that behavior:

        1. Change Categorical._concat_same_type to use union_categoricals
        2. Delete this method.
        """
        values = self._concatenator(
            [blk.values for blk in to_concat], axis=self.ndim - 1
        )
        # not using self.make_block_same_class as values can be object dtype
        return make_block(
            values, placement=placement or slice(0, len(values), 1), ndim=self.ndim
        )

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        filter=None,
        regex: bool = False,
        convert: bool = True,
    ):
        inplace = validate_bool_kwarg(inplace, "inplace")
        result = self if inplace else self.copy()
        if filter is None:  # replace was called on a series
            result.values.replace(to_replace, value, inplace=True)
            if convert:
                return result.convert(numeric=False, copy=not inplace)
            else:
                return result
        else:  # replace was called on a DataFrame
            if not isna(value):
                result.values.add_categories(value, inplace=True)
            return super(CategoricalBlock, result).replace(
                to_replace, value, inplace, filter, regex, convert
            )


# -----------------------------------------------------------------
# Constructor Helpers


def get_block_type(values, dtype=None):
    """
    Find the appropriate Block subclass to use for the given values and dtype.

    Parameters
    ----------
    values : ndarray-like
    dtype : numpy or pandas dtype

    Returns
    -------
    cls : class, subclass of Block
    """
    dtype = dtype or values.dtype
    vtype = dtype.type

    if is_sparse(dtype):
        # Need this first(ish) so that Sparse[datetime] is sparse
        cls = ExtensionBlock
    elif is_categorical(values):
        cls = CategoricalBlock
    elif issubclass(vtype, np.datetime64):
        assert not is_datetime64tz_dtype(values)
        cls = DatetimeBlock
    elif is_datetime64tz_dtype(values):
        cls = DatetimeTZBlock
    elif is_interval_dtype(dtype) or is_period_dtype(dtype):
        cls = ObjectValuesExtensionBlock
    elif is_extension_array_dtype(values):
        cls = ExtensionBlock
    elif issubclass(vtype, np.floating):
        cls = FloatBlock
    elif issubclass(vtype, np.timedelta64):
        assert issubclass(vtype, np.integer)
        cls = TimeDeltaBlock
    elif issubclass(vtype, np.complexfloating):
        cls = ComplexBlock
    elif issubclass(vtype, np.integer):
        cls = IntBlock
    elif dtype == np.bool_:
        cls = BoolBlock
    else:
        cls = ObjectBlock
    return cls


def make_block(values, placement, klass=None, ndim=None, dtype=None):
    # Ensure that we don't allow PandasArray / PandasDtype in internals.
    # For now, blocks should be backed by ndarrays when possible.
    if isinstance(values, ABCPandasArray):
        values = values.to_numpy()
        if ndim and ndim > 1:
            values = np.atleast_2d(values)

    if isinstance(dtype, PandasDtype):
        dtype = dtype.numpy_dtype

    if klass is None:
        dtype = dtype or values.dtype
        klass = get_block_type(values, dtype)

    elif klass is DatetimeTZBlock and not is_datetime64tz_dtype(values):
        # TODO: This is no longer hit internally; does it need to be retained
        #  for e.g. pyarrow?
        values = DatetimeArray._simple_new(values, dtype=dtype)

    return klass(values, ndim=ndim, placement=placement)


# -----------------------------------------------------------------


def _extend_blocks(result, blocks=None):
    """ return a new extended blocks, given the result """
    if blocks is None:
        blocks = []
    if isinstance(result, list):
        for r in result:
            if isinstance(r, list):
                blocks.extend(r)
            else:
                blocks.append(r)
    else:
        assert isinstance(result, Block), type(result)
        blocks.append(result)
    return blocks


def _block_shape(values, ndim=1, shape=None):
    """ guarantee the shape of the values to be at least 1 d """
    if values.ndim < ndim:
        if shape is None:
            shape = values.shape
        if not is_extension_array_dtype(values):
            # TODO: https://github.com/pandas-dev/pandas/issues/23023
            # block.shape is incorrect for "2D" ExtensionArrays
            # We can't, and don't need to, reshape.
            values = values.reshape(tuple((1,) + shape))
    return values


def _merge_blocks(blocks, dtype=None, _can_consolidate=True):

    if len(blocks) == 1:
        return blocks[0]

    if _can_consolidate:

        if dtype is None:
            if len({b.dtype for b in blocks}) != 1:
                raise AssertionError("_merge_blocks are invalid!")

        # FIXME: optimization potential in case all mgrs contain slices and
        # combination of those slices is a slice, too.
        new_mgr_locs = np.concatenate([b.mgr_locs.as_array for b in blocks])
        new_values = np.vstack([b.values for b in blocks])

        argsort = np.argsort(new_mgr_locs)
        new_values = new_values[argsort]
        new_mgr_locs = new_mgr_locs[argsort]

        return make_block(new_values, placement=new_mgr_locs)

    # no merge
    return blocks


def _safe_reshape(arr, new_shape):
    """
    If possible, reshape `arr` to have shape `new_shape`,
    with a couple of exceptions (see gh-13012):

    1) If `arr` is a ExtensionArray or Index, `arr` will be
       returned as is.
    2) If `arr` is a Series, the `_values` attribute will
       be reshaped and returned.

    Parameters
    ----------
    arr : array-like, object to be reshaped
    new_shape : int or tuple of ints, the new shape
    """
    if isinstance(arr, ABCSeries):
        arr = arr._values
    if not isinstance(arr, ABCExtensionArray):
        arr = arr.reshape(new_shape)
    return arr


def _putmask_smart(v, mask, n):
    """
    Return a new ndarray, try to preserve dtype if possible.

    Parameters
    ----------
    v : `values`, updated in-place (array like)
    mask : np.ndarray
        Applies to both sides (array like).
    n : `new values` either scalar or an array like aligned with `values`

    Returns
    -------
    values : ndarray with updated values
        this *may* be a copy of the original

    See Also
    --------
    ndarray.putmask
    """
    # we cannot use np.asarray() here as we cannot have conversions
    # that numpy does when numeric are mixed with strings

    # n should be the length of the mask or a scalar here
    if not is_list_like(n):
        n = np.repeat(n, len(mask))

    # see if we are only masking values that if putted
    # will work in the current dtype
    try:
        nn = n[mask]
    except TypeError:
        # TypeError: only integer scalar arrays can be converted to a scalar index
        pass
    else:
        # make sure that we have a nullable type
        # if we have nulls
        if not _isna_compat(v, nn[0]):
            pass
        elif not (is_float_dtype(nn.dtype) or is_integer_dtype(nn.dtype)):
            # only compare integers/floats
            pass
        elif not (is_float_dtype(v.dtype) or is_integer_dtype(v.dtype)):
            # only compare integers/floats
            pass
        else:

            # we ignore ComplexWarning here
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("ignore", np.ComplexWarning)
                nn_at = nn.astype(v.dtype)

            comp = nn == nn_at
            if is_list_like(comp) and comp.all():
                nv = v.copy()
                nv[mask] = nn_at
                return nv

    n = np.asarray(n)

    def _putmask_preserve(nv, n):
        try:
            nv[mask] = n[mask]
        except (IndexError, ValueError):
            nv[mask] = n
        return nv

    # preserves dtype if possible
    if v.dtype.kind == n.dtype.kind:
        return _putmask_preserve(v, n)

    # change the dtype if needed
    dtype, _ = maybe_promote(n.dtype)

    if is_extension_array_dtype(v.dtype) and is_object_dtype(dtype):
        v = v._internal_get_values(dtype)
    else:
        v = v.astype(dtype)

    return _putmask_preserve(v, n)
