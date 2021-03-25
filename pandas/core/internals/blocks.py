from datetime import datetime, timedelta
import inspect
import re
from typing import TYPE_CHECKING, Any, List, Optional, Type, Union, cast
import warnings

import numpy as np

from pandas._libs import NaT, algos as libalgos, internals as libinternals, lib, writers
from pandas._libs.internals import BlockPlacement
from pandas._libs.tslibs import conversion
from pandas._libs.tslibs.timezones import tz_compare
from pandas._typing import ArrayLike, Scalar, Shape
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import (
    astype_nansafe,
    convert_scalar_for_putitemlike,
    find_common_type,
    infer_dtype_from,
    infer_dtype_from_scalar,
    maybe_box_datetimelike,
    maybe_downcast_numeric,
    maybe_downcast_to_dtype,
    maybe_infer_dtype_type,
    maybe_promote,
    maybe_upcast,
    soft_convert_objects,
)
from pandas.core.dtypes.common import (
    DT64NS_DTYPE,
    TD64NS_DTYPE,
    is_bool_dtype,
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float,
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
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndexClass,
    ABCPandasArray,
    ABCSeries,
)
from pandas.core.dtypes.missing import is_valid_nat_for_dtype, isna, isna_compat

import pandas.core.algorithms as algos
from pandas.core.array_algos.replace import compare_or_regex_search, replace_regex
from pandas.core.array_algos.transforms import shift
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

if TYPE_CHECKING:
    from pandas import Index


class Block(PandasObject):
    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas
    data structure

    Index-ignorant; let the container take care of that
    """

    values: Union[np.ndarray, ExtensionArray]

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
    is_extension = False
    _can_hold_na = False
    _can_consolidate = True
    _validate_ndim = True

    @classmethod
    def _simple_new(
        cls, values: ArrayLike, placement: BlockPlacement, ndim: int
    ) -> "Block":
        """
        Fastpath constructor, does *no* validation
        """
        obj = object.__new__(cls)
        obj.ndim = ndim
        obj.values = values
        obj._mgr_locs = placement
        return obj

    def __init__(self, values, placement, ndim: int):
        """
        Parameters
        ----------
        values : np.ndarray or ExtensionArray
        placement : BlockPlacement (or castable)
        ndim : int
            1 for SingleBlockManager/Series, 2 for BlockManager/DataFrame
        """
        # TODO(EA2D): ndim will be unnecessary with 2D EAs
        self.ndim = self._check_ndim(values, ndim)
        self.mgr_locs = placement
        self.values = self._maybe_coerce_values(values)

        if self._validate_ndim and self.ndim and len(self.mgr_locs) != len(self.values):
            raise ValueError(
                f"Wrong number of items passed {len(self.values)}, "
                f"placement implies {len(self.mgr_locs)}"
            )

    def _maybe_coerce_values(self, values):
        """
        Ensure we have correctly-typed values.

        Parameters
        ----------
        values : np.ndarray, ExtensionArray, Index

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        return values

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
        return self._can_consolidate, self.dtype.name

    @property
    def is_view(self) -> bool:
        """ return a boolean if I am possibly a view """
        values = self.values
        values = cast(np.ndarray, values)
        return values.base is not None

    @property
    def is_categorical(self) -> bool:
        return self._holder is Categorical

    @property
    def is_datelike(self) -> bool:
        """ return True if I am a non-datelike """
        return self.is_datetime or self.is_timedelta

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

    def get_block_values_for_json(self) -> np.ndarray:
        """
        This is used in the JSON C code.
        """
        # TODO(EA2D): reshape will be unnecessary with 2D EAs
        return np.asarray(self.values).reshape(self.shape)

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

    def make_block(self, values, placement=None) -> "Block":
        """
        Create a new block, with type inference propagate any values that are
        not specified
        """
        if placement is None:
            placement = self.mgr_locs
        if self.is_extension:
            values = _block_shape(values, ndim=self.ndim)

        return make_block(values, placement=placement, ndim=self.ndim)

    def make_block_same_class(self, values, placement=None, ndim=None):
        """ Wrap given values in a block of same type as self. """
        if placement is None:
            placement = self.mgr_locs
        if ndim is None:
            ndim = self.ndim
        return type(self)(values, placement=placement, ndim=ndim)

    def __repr__(self) -> str:
        # don't want to print out all of the items here
        name = type(self).__name__
        if self.ndim == 1:
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
        elif not isinstance(new_mgr_locs, BlockPlacement):
            new_mgr_locs = BlockPlacement(new_mgr_locs)

        new_values = self._slice(slicer)

        if self._validate_ndim and new_values.ndim != self.ndim:
            raise ValueError("Only same dim slicing is allowed")

        return type(self)._simple_new(new_values, new_mgr_locs, self.ndim)

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def iget(self, i):
        return self.values[i]

    def set_inplace(self, locs, values):
        """
        Modify block values in-place with new item value.

        Notes
        -----
        `set` never creates a new array or new Block, whereas `setitem` _may_
        create a new array and always creates a new Block.
        """
        self.values[locs] = values

    def delete(self, loc) -> None:
        """
        Delete given loc(-s) from block in-place.
        """
        self.values = np.delete(self.values, loc, 0)
        self.mgr_locs = self.mgr_locs.delete(loc)

    def apply(self, func, **kwargs) -> List["Block"]:
        """
        apply the function to my values; return a block if we are not
        one
        """
        with np.errstate(all="ignore"):
            result = func(self.values, **kwargs)

        return self._split_op_result(result)

    def reduce(self, func, ignore_failures: bool = False) -> List["Block"]:
        # We will apply the function and reshape the result into a single-row
        #  Block with the same mgr_locs; squeezing will be done at a higher level
        assert self.ndim == 2

        try:
            result = func(self.values)
        except (TypeError, NotImplementedError):
            if ignore_failures:
                return []
            raise

        if np.ndim(result) == 0:
            # TODO(EA2D): special case not needed with 2D EAs
            res_values = np.array([[result]])
        else:
            res_values = result.reshape(-1, 1)

        nb = self.make_block(res_values)
        return [nb]

    def _split_op_result(self, result) -> List["Block"]:
        # See also: split_and_operate
        if is_extension_array_dtype(result) and result.ndim > 1:
            # TODO(EA2D): unnecessary with 2D EAs
            # if we get a 2D ExtensionArray, we need to split it into 1D pieces
            nbs = []
            for i, loc in enumerate(self.mgr_locs):
                vals = result[i]
                block = self.make_block(values=vals, placement=[loc])
                nbs.append(block)
            return nbs

        if not isinstance(result, Block):
            result = self.make_block(result)

        return [result]

    def fillna(
        self, value, limit=None, inplace: bool = False, downcast=None
    ) -> List["Block"]:
        """
        fillna on the block with the value. If we fail, then convert to
        ObjectBlock and try again
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        mask = isna(self.values)
        mask = _extract_bool_array(mask)
        if limit is not None:
            limit = libalgos.validate_limit(None, limit=limit)
            mask[mask.cumsum(self.ndim - 1) > limit] = False

        if not self._can_hold_na:
            if inplace:
                return [self]
            else:
                return [self.copy()]

        if self._can_hold_element(value):
            nb = self if inplace else self.copy()
            nb._putmask_simple(mask, value)
            # TODO: should be nb._maybe_downcast?
            return self._maybe_downcast([nb], downcast)

        # we can't process the value, but nothing to do
        if not mask.any():
            return [self] if inplace else [self.copy()]

        # operate column-by-column
        def f(mask, val, idx):
            block = self.coerce_to_target_dtype(value)

            # slice out our block
            if idx is not None:
                # i.e. self.ndim == 2
                block = block.getitem_block(slice(idx, idx + 1))
            return block.fillna(value, limit=limit, inplace=inplace, downcast=None)

        return self.split_and_operate(None, f, inplace)

    def _split(self) -> List["Block"]:
        """
        Split a block into a list of single-column blocks.
        """
        assert self.ndim == 2

        new_blocks = []
        for i, ref_loc in enumerate(self.mgr_locs):
            vals = self.values[slice(i, i + 1)]

            nb = self.make_block(vals, [ref_loc])
            new_blocks.append(nb)
        return new_blocks

    def split_and_operate(
        self, mask, f, inplace: bool, ignore_failures: bool = False
    ) -> List["Block"]:
        """
        split the block per-column, and apply the callable f
        per-column, return a new block for each. Handle
        masking which will not change a block unless needed.

        Parameters
        ----------
        mask : 2-d boolean mask
        f : callable accepting (1d-mask, 1d values, indexer)
        inplace : bool
        ignore_failures : bool, default False

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
            if m.any() or m.size == 0:
                # Apply our function; we may ignore_failures if this is a
                #  reduction that is dropping nuisance columns GH#37827
                try:
                    nv = f(m, v, i)
                except TypeError:
                    if ignore_failures:
                        continue
                    else:
                        raise
            else:
                nv = v if inplace else v.copy()

            block = make_a_block(nv, [ref_loc])
            new_blocks.append(block)

        return new_blocks

    def _maybe_downcast(self, blocks: List["Block"], downcast=None) -> List["Block"]:

        # no need to downcast our float
        # unless indicated
        if downcast is None and (self.is_float or self.is_datelike):
            return blocks

        return extend_blocks([b.downcast(downcast) for b in blocks])

    def downcast(self, dtypes=None) -> List["Block"]:
        """ try to downcast each item to the dict of dtypes if present """
        # turn it off completely
        if dtypes is False:
            return [self]

        values = self.values

        if self.ndim == 1:

            # try to cast all non-floats here
            if dtypes is None:
                dtypes = "infer"

            nv = maybe_downcast_to_dtype(values, dtypes)
            return [self.make_block(nv)]

        # ndim > 1
        if dtypes is None:
            return [self]

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

        if dtype is not None:
            dtype = pandas_dtype(dtype)

        # may need to convert to categorical
        if is_categorical_dtype(dtype):

            if is_categorical_dtype(self.values.dtype):
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
            try:
                values = self.values.astype(dtype)
            except (ValueError, TypeError):
                if errors == "ignore":
                    values = self.values
                else:
                    raise
        else:
            if issubclass(dtype.type, str):

                # use native type formatting for datetime/tz/timedelta
                if self.is_datelike:
                    values = self.to_native_types().values

                # astype formatting
                else:
                    # Because we have neither is_extension nor is_datelike,
                    #  self.values already has the correct shape
                    values = self.values

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

        # TODO(EA2D): special case not needed with 2D EAs
        if isinstance(values, np.ndarray):
            values = values.reshape(self.shape)

        newb = self.make_block(values)

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
    ) -> List["Block"]:
        """
        attempt to coerce any object types to better types return a copy
        of the block (if copy = True) by definition we are not an ObjectBlock
        here!
        """
        return [self.copy()] if copy else [self]

    def _can_hold_element(self, element: Any) -> bool:
        """ require the same dtype as ourselves """
        dtype = self.values.dtype.type
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, dtype)
        return isinstance(element, dtype)

    def should_store(self, value: ArrayLike) -> bool:
        """
        Should we set self.values[indexer] = value inplace or do we need to cast?

        Parameters
        ----------
        value : np.ndarray or ExtensionArray

        Returns
        -------
        bool
        """
        return is_dtype_equal(value.dtype, self.dtype)

    def to_native_types(self, na_rep="nan", quoting=None, **kwargs):
        """ convert to our native types format """
        values = self.values

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
        return self.make_block(values)

    # block actions #
    def copy(self, deep: bool = True):
        """ copy constructor """
        values = self.values
        if deep:
            values = values.copy()
        return self.make_block_same_class(values, ndim=self.ndim)

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        regex: bool = False,
    ) -> List["Block"]:
        """
        replace the to_replace value with value, possible to create new
        blocks here this is just a call to putmask. regex is not used here.
        It is used in ObjectBlocks.  It is here for API compatibility.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        original_to_replace = to_replace

        if not self._can_hold_element(to_replace):
            # We cannot hold `to_replace`, so we know immediately that
            #  replacing it is a no-op.
            # Note: If to_replace were a list, NDFrame.replace would call
            #  replace_list instead of replace.
            return [self] if inplace else [self.copy()]

        values = self.values
        if lib.is_scalar(to_replace) and isinstance(values, np.ndarray):
            # The only non-DatetimeLike class that also has a non-trivial
            #  try_coerce_args is ObjectBlock, but that overrides replace,
            #  so does not get here.
            to_replace = convert_scalar_for_putitemlike(to_replace, values.dtype)

        mask = missing.mask_missing(values, to_replace)
        if not mask.any():
            # Note: we get here with test_replace_extension_other incorrectly
            #  bc _can_hold_element is incorrect.
            return [self] if inplace else [self.copy()]

        if not self._can_hold_element(value):
            blk = self.astype(object)
            return blk.replace(
                to_replace=original_to_replace,
                value=value,
                inplace=True,
                regex=regex,
            )

        blk = self if inplace else self.copy()
        blk._putmask_simple(mask, value)
        blocks = blk.convert(numeric=False, copy=not inplace)
        return blocks

    def _replace_regex(
        self,
        to_replace,
        value,
        inplace: bool = False,
        convert: bool = True,
        mask=None,
    ) -> List["Block"]:
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
        convert : bool, default True
            If true, try to coerce any object types to better types.
        mask : array-like of bool, optional
            True indicate corresponding element is ignored.

        Returns
        -------
        List[Block]
        """
        if not self._can_hold_element(to_replace):
            # i.e. only ObjectBlock, but could in principle include a
            #  String ExtensionBlock
            return [self] if inplace else [self.copy()]

        rx = re.compile(to_replace)

        new_values = self.values if inplace else self.values.copy()
        replace_regex(new_values, rx, value, mask)

        block = self.make_block(new_values)
        if convert:
            nbs = block.convert(numeric=False)
        else:
            nbs = [block]
        return nbs

    def _replace_list(
        self,
        src_list: List[Any],
        dest_list: List[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> List["Block"]:
        """
        See BlockManager._replace_list docstring.
        """
        # Exclude anything that we know we won't contain
        pairs = [
            (x, y) for x, y in zip(src_list, dest_list) if self._can_hold_element(x)
        ]
        if not len(pairs):
            # shortcut, nothing to replace
            return [self] if inplace else [self.copy()]

        src_len = len(pairs) - 1

        def comp(s: Scalar, mask: np.ndarray, regex: bool = False) -> np.ndarray:
            """
            Generate a bool array by perform an equality check, or perform
            an element-wise regular expression matching
            """
            if isna(s):
                return ~mask

            s = maybe_box_datetimelike(s)
            return compare_or_regex_search(self.values, s, regex, mask)

        if self.is_object:
            # Calculate the mask once, prior to the call of comp
            # in order to avoid repeating the same computations
            mask = ~isna(self.values)
            masks = [comp(s[0], mask, regex) for s in pairs]
        else:
            # GH#38086 faster if we know we dont need to check for regex
            masks = [missing.mask_missing(self.values, s[0]) for s in pairs]

        masks = [_extract_bool_array(x) for x in masks]

        rb = [self if inplace else self.copy()]
        for i, (src, dest) in enumerate(pairs):
            convert = i == src_len  # only convert once at the end
            new_rb: List["Block"] = []

            # GH-39338: _replace_coerce can split a block into
            # single-column blocks, so track the index so we know
            # where to index into the mask
            for blk_num, blk in enumerate(rb):
                if len(rb) == 1:
                    m = masks[i]
                else:
                    mib = masks[i]
                    assert not isinstance(mib, bool)
                    m = mib[blk_num : blk_num + 1]

                result = blk._replace_coerce(
                    to_replace=src,
                    value=dest,
                    mask=m,
                    inplace=inplace,
                    regex=regex,
                )
                if convert and blk.is_object:
                    result = extend_blocks(
                        [b.convert(numeric=False, copy=True) for b in result]
                    )
                new_rb.extend(result)
            rb = new_rb
        return rb

    def setitem(self, indexer, value):
        """
        Attempt self.values[indexer] = value, possibly creating a new array.

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

            if hasattr(value, "dtype"):
                dtype = value.dtype

            elif lib.is_scalar(value) and not isna(value):
                dtype, _ = infer_dtype_from_scalar(value, pandas_dtype=True)

            else:
                # e.g. we are bool dtype and value is nan
                # TODO: watch out for case with listlike value and scalar/empty indexer
                dtype, _ = maybe_promote(np.array(value).dtype)
                return self.astype(dtype).setitem(indexer, value)

            dtype = find_common_type([values.dtype, dtype])
            assert not is_dtype_equal(self.dtype, dtype)
            # otherwise should have _can_hold_element

            return self.astype(dtype).setitem(indexer, value)

        # value must be storable at this moment
        if is_extension_array_dtype(getattr(value, "dtype", None)):
            # We need to be careful not to allow through strings that
            #  can be parsed to EADtypes
            is_ea_value = True
            arr_value = value
        else:
            is_ea_value = False
            arr_value = np.array(value)

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

        elif is_scalar_indexer(indexer, self.ndim):
            # setting a single element for each dim and with a rhs that could
            #  be e.g. a list; see GH#6043
            values[indexer] = value

        elif exact_match and is_categorical_dtype(arr_value.dtype):
            # GH25495 - If the current dtype is not categorical,
            # we need to create a new categorical block
            values[indexer] = value
            return self.make_block(Categorical(self.values, dtype=arr_value.dtype))

        elif exact_match and is_ea_value:
            # GH#32395 if we're going to replace the values entirely, just
            #  substitute in the new array
            return self.make_block(arr_value)

        # if we are an exact match (ex-broadcasting),
        # then use the resultant dtype
        elif exact_match:
            # We are setting _all_ of the array's values, so can cast to new dtype
            values[indexer] = value

            values = values.astype(arr_value.dtype, copy=False)

        # set
        else:
            values[indexer] = value

        if transpose:
            values = values.T
        block = self.make_block(values)
        return block

    def _putmask_simple(self, mask: np.ndarray, value: Any):
        """
        Like putmask but

        a) we do not cast on failure
        b) we do not handle repeating or truncating like numpy.

        Parameters
        ----------
        mask : np.ndarray[bool]
            We assume _extract_bool_array has already been called.
        value : Any
            We assume self._can_hold_element(value)
        """
        values = self.values

        if lib.is_scalar(value) and isinstance(values, np.ndarray):
            value = convert_scalar_for_putitemlike(value, values.dtype)

        if self.is_extension or (self.is_object and not lib.is_scalar(value)):
            # GH#19266 using np.putmask gives unexpected results with listlike value
            if is_list_like(value) and len(value) == len(values):
                values[mask] = value[mask]
            else:
                values[mask] = value
        else:
            # GH#37833 np.putmask is more performant than __setitem__
            np.putmask(values, mask, value)

    def putmask(
        self, mask, new, inplace: bool = False, axis: int = 0, transpose: bool = False
    ) -> List["Block"]:
        """
        putmask the data to the block; it is possible that we may create a
        new dtype of block

        Return the resulting block(s).

        Parameters
        ----------
        mask : np.ndarray[bool], SparseArray[bool], or BooleanArray
        new : a ndarray/object
        inplace : bool, default False
            Perform inplace modification.
        axis : int
        transpose : bool, default False
            Set to True if self is stored with axes reversed.

        Returns
        -------
        List[Block]
        """
        mask = _extract_bool_array(mask)
        assert not isinstance(new, (ABCIndexClass, ABCSeries, ABCDataFrame))

        new_values = self.values  # delay copy if possible.
        # if we are passed a scalar None, convert it here
        if not is_list_like(new) and isna(new) and not self.is_object:
            # FIXME: make sure we have compatible NA
            new = self.fill_value

        if self._can_hold_element(new):
            # We only get here for non-Extension Blocks, so _try_coerce_args
            #  is only relevant for DatetimeBlock and TimedeltaBlock
            if lib.is_scalar(new):
                new = convert_scalar_for_putitemlike(new, self.values.dtype)

            if transpose:
                new_values = new_values.T

            # If the default repeat behavior in np.putmask would go in the
            # wrong direction, then explicitly repeat and reshape new instead
            if getattr(new, "ndim", 0) >= 1:
                if self.ndim - 1 == new.ndim and axis == 1:
                    new = np.repeat(new, new_values.shape[-1]).reshape(self.shape)
                new = new.astype(new_values.dtype)

            if new_values is self.values and not inplace:
                new_values = new_values.copy()
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
            if new_values is None:
                new_values = self.values if inplace else self.values.copy()
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

        elif self.is_datetime or is_datetime64_any_dtype(dtype):
            # The is_dtype_equal check above ensures that at most one of
            #  these two conditions hold, so we must cast to object.
            return self.astype(object)

        elif self.is_timedelta or is_timedelta64_dtype(dtype):
            # The is_dtype_equal check above ensures that at most one of
            #  these two conditions hold, so we must cast to object.
            return self.astype(object)

        try:
            return self.astype(dtype)
        except (ValueError, TypeError, OverflowError):
            return self.astype(object)

    def interpolate(
        self,
        method: str = "pad",
        axis: int = 0,
        index: Optional["Index"] = None,
        inplace: bool = False,
        limit: Optional[int] = None,
        limit_direction: str = "forward",
        limit_area: Optional[str] = None,
        fill_value: Optional[Any] = None,
        coerce: bool = False,
        downcast: Optional[str] = None,
        **kwargs,
    ):

        inplace = validate_bool_kwarg(inplace, "inplace")

        if not self._can_hold_na:
            # If there are no NAs, then interpolate is a no-op
            return self if inplace else self.copy()

        # a fill na type method
        try:
            m = missing.clean_fill_method(method)
        except ValueError:
            m = None

        if m is not None:
            if fill_value is not None:
                # similar to validate_fillna_kwargs
                raise ValueError("Cannot pass both fill_value and method")

            return self._interpolate_with_fill(
                method=m,
                axis=axis,
                inplace=inplace,
                limit=limit,
                limit_area=limit_area,
                downcast=downcast,
            )
        # validate the interp method
        m = missing.clean_interp_method(method, **kwargs)

        assert index is not None  # for mypy

        return self._interpolate(
            method=m,
            index=index,
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
        method: str = "pad",
        axis: int = 0,
        inplace: bool = False,
        limit: Optional[int] = None,
        limit_area: Optional[str] = None,
        downcast: Optional[str] = None,
    ) -> List["Block"]:
        """ fillna but using the interpolate machinery """
        inplace = validate_bool_kwarg(inplace, "inplace")

        assert self._can_hold_na  # checked by caller

        values = self.values if inplace else self.values.copy()

        values = missing.interpolate_2d(
            values,
            method=method,
            axis=axis,
            limit=limit,
            limit_area=limit_area,
        )

        blocks = [self.make_block_same_class(values, ndim=self.ndim)]
        return self._maybe_downcast(blocks, downcast)

    def _interpolate(
        self,
        method: str,
        index: "Index",
        fill_value: Optional[Any] = None,
        axis: int = 0,
        limit: Optional[int] = None,
        limit_direction: str = "forward",
        limit_area: Optional[str] = None,
        inplace: bool = False,
        downcast: Optional[str] = None,
        **kwargs,
    ) -> List["Block"]:
        """ interpolate using scipy wrappers """
        inplace = validate_bool_kwarg(inplace, "inplace")
        data = self.values if inplace else self.values.copy()

        # only deal with floats
        if not self.is_float:
            if not self.is_integer:
                return [self]
            data = data.astype(np.float64)

        if fill_value is None:
            fill_value = self.fill_value

        if method in ("krogh", "piecewise_polynomial", "pchip"):
            if not index.is_monotonic:
                raise ValueError(
                    f"{method} interpolation requires that the index be monotonic."
                )
        # process 1-d slices in the axis direction

        def func(yvalues: np.ndarray) -> np.ndarray:

            # process a 1-d slice, returning it
            # should the axis argument be handled below in apply_along_axis?
            # i.e. not an arg to missing.interpolate_1d
            return missing.interpolate_1d(
                xvalues=index,
                yvalues=yvalues,
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

    def take_nd(self, indexer, axis: int, new_mgr_locs=None, fill_value=lib.no_default):
        """
        Take values according to indexer and return them as a block.bb

        """
        # algos.take_nd dispatches for DatetimeTZBlock, CategoricalBlock
        # so need to preserve types
        # sparse is treated like an ndarray, but needs .get_values() shaping

        values = self.values

        if fill_value is lib.no_default:
            fill_value = self.fill_value
            allow_fill = False
        else:
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
        return [self.make_block(values=new_values)]

    def shift(self, periods: int, axis: int = 0, fill_value=None):
        """ shift the block by periods, possibly upcast """
        # convert integer to float if necessary. need to do a lot more than
        # that, handle boolean etc also
        new_values, fill_value = maybe_upcast(self.values, fill_value)

        new_values = shift(new_values, periods, axis, fill_value)

        return [self.make_block(new_values)]

    def where(
        self, other, cond, errors="raise", try_cast: bool = False, axis: int = 0
    ) -> List["Block"]:
        """
        evaluate the block; return result block(s) from the result

        Parameters
        ----------
        other : a ndarray/object
        cond : np.ndarray[bool], SparseArray[bool], or BooleanArray
        errors : str, {'raise', 'ignore'}, default 'raise'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object
        try_cast: bool, default False
        axis : int, default 0

        Returns
        -------
        List[Block]
        """
        import pandas.core.computation.expressions as expressions

        cond = _extract_bool_array(cond)
        assert not isinstance(other, (ABCIndexClass, ABCSeries, ABCDataFrame))

        assert errors in ["raise", "ignore"]
        transpose = self.ndim == 2

        values = self.values
        orig_other = other
        if transpose:
            values = values.T

        # If the default broadcasting would go in the wrong direction, then
        # explicitly reshape other instead
        if getattr(other, "ndim", 0) >= 1:
            if values.ndim - 1 == other.ndim and axis == 1:
                other = other.reshape(tuple(other.shape + (1,)))
            elif transpose and values.ndim == self.ndim - 1:
                # TODO(EA2D): not neceesssary with 2D EAs
                cond = cond.T

        if not hasattr(cond, "shape"):
            raise ValueError("where must have a condition that is ndarray like")

        if cond.ravel("K").all():
            result = values.copy()
        else:
            # see if we can operate on the entire block, or need item-by-item
            # or if we are a single block (ndim == 1)
            if (
                (self.is_integer or self.is_bool)
                and lib.is_float(other)
                and np.isnan(other)
            ):
                # GH#3733 special case to avoid object-dtype casting
                #  and go through numexpr path instead.
                # In integer case, np.where will cast to floats
                pass
            elif not self._can_hold_element(other):
                # we cannot coerce, return a compat dtype
                # we are explicitly ignoring errors
                block = self.coerce_to_target_dtype(other)
                blocks = block.where(
                    orig_other, cond, errors=errors, try_cast=try_cast, axis=axis
                )
                return self._maybe_downcast(blocks, "infer")

            if not (
                (self.is_integer or self.is_bool)
                and lib.is_float(other)
                and np.isnan(other)
            ):
                # convert datetime to datetime64, timedelta to timedelta64
                other = convert_scalar_for_putitemlike(other, values.dtype)

            # By the time we get here, we should have all Series/Index
            #  args extracted to ndarray
            result = expressions.where(cond, values, other)

        if self._can_hold_na or self.ndim == 1:

            if transpose:
                result = result.T

            return [self.make_block(result)]

        # might need to separate out blocks
        axis = cond.ndim - 1
        cond = cond.swapaxes(axis, 0)
        mask = np.array([cond[i].all() for i in range(cond.shape[0])], dtype=bool)

        result_blocks: List["Block"] = []
        for m in [mask, ~mask]:
            if m.any():
                result = cast(np.ndarray, result)  # EABlock overrides where
                taken = result.take(m.nonzero()[0], axis=axis)
                r = maybe_downcast_numeric(taken, self.dtype)
                nb = self.make_block(r.T, placement=self.mgr_locs[m])
                result_blocks.append(nb)

        return result_blocks

    def _unstack(self, unstacker, fill_value, new_placement):
        """
        Return a list of unstacked blocks of self

        Parameters
        ----------
        unstacker : reshape._Unstacker
        fill_value : int
            Only used in ExtensionBlock._unstack

        Returns
        -------
        blocks : list of Block
            New blocks of unstacked values.
        mask : array_like of bool
            The mask of columns of `blocks` we should keep.
        """
        new_values, mask = unstacker.get_new_values(
            self.values.T, fill_value=fill_value
        )

        mask = mask.any(0)
        # TODO: in all tests we have mask.all(); can we rely on that?

        new_values = new_values.T[mask]
        new_placement = new_placement[mask]

        blocks = [make_block(new_values, placement=new_placement)]
        return blocks, mask

    def quantile(self, qs, interpolation="linear", axis: int = 0):
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
        self,
        to_replace,
        value,
        mask: np.ndarray,
        inplace: bool = True,
        regex: bool = False,
    ) -> List["Block"]:
        """
        Replace value corresponding to the given boolean array with another
        value.

        Parameters
        ----------
        to_replace : object or pattern
            Scalar to replace or regular expression to match.
        value : object
            Replacement object.
        mask : np.ndarray[bool]
            True indicate corresponding element is ignored.
        inplace : bool, default True
            Perform inplace modification.
        regex : bool, default False
            If true, perform regular expression substitution.

        Returns
        -------
        List[Block]
        """
        if mask.any():
            if not regex:
                nb = self.coerce_to_target_dtype(value)
                if nb is self and not inplace:
                    nb = nb.copy()
                nb._putmask_simple(mask, value)
                return [nb]
            else:
                regex = _should_use_regex(regex, to_replace)
                if regex:
                    return self._replace_regex(
                        to_replace,
                        value,
                        inplace=inplace,
                        convert=False,
                        mask=mask,
                    )
                return self.replace(to_replace, value, inplace=inplace, regex=False)
        return [self]


class ExtensionBlock(Block):
    """
    Block for holding extension types.

    Notes
    -----
    This holds all 3rd-party extension array types. It's also the immediate
    parent class for our internal extension types' blocks, CategoricalBlock.

    ExtensionArrays are limited to 1-D.
    """

    _can_consolidate = False
    _validate_ndim = False
    is_extension = True

    values: ExtensionArray

    def __init__(self, values, placement, ndim: int):
        """
        Initialize a non-consolidatable block.

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

        if self.ndim == 2 and len(self.mgr_locs) != 1:
            # TODO(EA2D): check unnecessary with 2D EAs
            raise AssertionError("block.size != values.size")

    @property
    def shape(self):
        # TODO(EA2D): override unnecessary with 2D EAs
        if self.ndim == 1:
            return (len(self.values),)
        return len(self.mgr_locs), len(self.values)

    def iget(self, col):

        if self.ndim == 2 and isinstance(col, tuple):
            # TODO(EA2D): unnecessary with 2D EAs
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

    def set_inplace(self, locs, values):
        # NB: This is a misnomer, is supposed to be inplace but is not,
        #  see GH#33457
        assert locs.tolist() == [0]
        self.values = values

    def putmask(
        self, mask, new, inplace: bool = False, axis: int = 0, transpose: bool = False
    ) -> List["Block"]:
        """
        See Block.putmask.__doc__
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        mask = _extract_bool_array(mask)

        new_values = self.values if inplace else self.values.copy()

        if isinstance(new, (np.ndarray, ExtensionArray)) and len(new) == len(mask):
            new = new[mask]

        mask = safe_reshape(mask, new_values.shape)

        new_values[mask] = new
        return [self.make_block(values=new_values)]

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
    def is_view(self) -> bool:
        """Extension arrays are never treated as views."""
        return False

    @property
    def is_numeric(self):
        return self.values.dtype._is_numeric

    def setitem(self, indexer, value):
        """
        Attempt self.values[indexer] = value, possibly creating a new array.

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
        if not self._can_hold_element(value):
            # This is only relevant for DatetimeTZBlock, which has a
            #  non-trivial `_can_hold_element`.
            # https://github.com/pandas-dev/pandas/issues/24020
            # Need a dedicated setitem until GH#24020 (type promotion in setitem
            #  for extension arrays) is designed and implemented.
            return self.astype(object).setitem(indexer, value)

        if isinstance(indexer, tuple):
            # TODO(EA2D): not needed with 2D EAs
            # we are always 1-D
            indexer = indexer[0]

        check_setitem_lengths(indexer, value, self.values)
        self.values[indexer] = value
        return self

    def get_values(self, dtype=None):
        # ExtensionArrays must be iterable, so this works.
        # TODO(EA2D): reshape not needed with 2D EAs
        return np.asarray(self.values).reshape(self.shape)

    def array_values(self) -> ExtensionArray:
        return self.values

    def to_native_types(self, na_rep="nan", quoting=None, **kwargs):
        """override to use ExtensionArray astype for the conversion"""
        values = self.values
        mask = isna(values)

        values = np.asarray(values.astype(object))
        values[mask] = na_rep

        # TODO(EA2D): reshape not needed with 2D EAs
        # we are expected to return a 2-d ndarray
        return self.make_block(values)

    def take_nd(
        self, indexer, axis: int = 0, new_mgr_locs=None, fill_value=lib.no_default
    ):
        """
        Take values according to indexer and return them as a block.
        """
        if fill_value is lib.no_default:
            fill_value = None

        # TODO(EA2D): special case not needed with 2D EAs
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
        # TODO: We may need to think about pushing this onto the array.
        # We're doing the same as CategoricalBlock here.
        return True

    def _slice(self, slicer):
        """
        Return a slice of my values.

        Parameters
        ----------
        slicer : slice, ndarray[int], or a tuple of these
            Valid (non-reducing) indexer for self.values.

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        # return same dims as we currently have
        if not isinstance(slicer, tuple) and self.ndim == 2:
            # reached via getitem_block via _slice_take_blocks_ax0
            # TODO(EA2D): wont be necessary with 2D EAs
            slicer = (slicer, slice(None))

        if isinstance(slicer, tuple) and len(slicer) == 2:
            first = slicer[0]
            if not isinstance(first, slice):
                raise AssertionError(
                    "invalid slicing for a 1-ndim ExtensionArray", first
                )
            # GH#32959 only full-slicers along fake-dim0 are valid
            # TODO(EA2D): wont be necessary with 2D EAs
            new_locs = self.mgr_locs[first]
            if len(new_locs):
                # effectively slice(None)
                slicer = slicer[1]
            else:
                raise AssertionError(
                    "invalid slicing for a 1-ndim ExtensionArray", slicer
                )

        return self.values[slicer]

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
        if axis == 0 and n != 0:
            # n==0 case will be a no-op so let is fall through
            # Since we only have one column, the result will be all-NA.
            #  Create this result by shifting along axis=0 past the length of
            #  our values.
            return super().diff(len(self.values), axis=0)
        if axis == 1:
            # TODO(EA2D): unnecessary with 2D EAs
            # we are by definition 1D.
            axis = 0
        return super().diff(n, axis)

    def shift(
        self, periods: int, axis: int = 0, fill_value: Any = None
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
        self, other, cond, errors="raise", try_cast: bool = False, axis: int = 0
    ) -> List["Block"]:

        cond = _extract_bool_array(cond)
        assert not isinstance(other, (ABCIndexClass, ABCSeries, ABCDataFrame))

        if isinstance(other, np.ndarray) and other.ndim == 2:
            # TODO(EA2D): unnecessary with 2D EAs
            assert other.shape[1] == 1
            other = other[:, 0]

        if isinstance(cond, np.ndarray) and cond.ndim == 2:
            # TODO(EA2D): unnecessary with 2D EAs
            assert cond.shape[1] == 1
            cond = cond[:, 0]

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

    def _unstack(self, unstacker, fill_value, new_placement):
        # ExtensionArray-safe unstack.
        # We override ObjectBlock._unstack, which unstacks directly on the
        # values of the array. For EA-backed blocks, this would require
        # converting to a 2-D ndarray of objects.
        # Instead, we unstack an ndarray of integer positions, followed by
        # a `take` on the actual values.
        n_rows = self.shape[-1]
        dummy_arr = np.arange(n_rows)

        new_values, mask = unstacker.get_new_values(dummy_arr, fill_value=-1)
        mask = mask.any(0)
        # TODO: in all tests we have mask.all(); can we rely on that?

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

    def _can_hold_element(self, element: Any) -> bool:
        if is_valid_nat_for_dtype(element, self.dtype):
            return True
        if isinstance(element, list) and len(element) == 0:
            return True
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, self.dtype.type)
        return isinstance(element, self.dtype.type)


class NumericBlock(Block):
    __slots__ = ()
    is_numeric = True
    _can_hold_na = True


class FloatBlock(NumericBlock):
    __slots__ = ()
    is_float = True

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, (np.floating, np.integer)) and not issubclass(
                tipo.type, np.timedelta64
            )
        return isinstance(
            element, (float, int, np.floating, np.int_)
        ) and not isinstance(
            element,
            (bool, np.bool_, np.timedelta64),
        )

    def to_native_types(
        self, na_rep="", float_format=None, decimal=".", quoting=None, **kwargs
    ):
        """ convert to our native types format """
        values = self.values

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
            return self.make_block(values)

        from pandas.io.formats.format import FloatArrayFormatter

        formatter = FloatArrayFormatter(
            values,
            na_rep=na_rep,
            float_format=float_format,
            decimal=decimal,
            quoting=quoting,
            fixed_width=False,
        )
        res = formatter.get_result_as_array()
        return self.make_block(res)


class ComplexBlock(NumericBlock):
    __slots__ = ()
    is_complex = True

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, (np.floating, np.integer, np.complexfloating))
        return isinstance(
            element, (float, int, complex, np.float_, np.int_)
        ) and not isinstance(element, (bool, np.bool_))


class IntBlock(NumericBlock):
    __slots__ = ()
    is_integer = True
    _can_hold_na = False

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return (
                issubclass(tipo.type, np.integer)
                and not issubclass(tipo.type, np.timedelta64)
                and self.dtype.itemsize >= tipo.itemsize
            )
        # We have not inferred an integer from the dtype
        # check if we have a builtin int or a float equal to an int
        return is_integer(element) or (is_float(element) and element.is_integer())


class DatetimeLikeBlockMixin(Block):
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
            # DTA/TDA constructor and astype can handle 2D
            return self._holder(self.values).astype(object)
        return self.values

    def internal_values(self):
        # Override to return DatetimeArray and TimedeltaArray
        return self.array_values()

    def array_values(self):
        return self._holder._simple_new(self.values)

    def iget(self, key):
        # GH#31649 we need to wrap scalars in Timestamp/Timedelta
        # TODO(EA2D): this can be removed if we ever have 2D EA
        return self.array_values().reshape(self.shape)[key]

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
        # TODO(EA2D): reshape not necessary with 2D EAs
        values = self.array_values().reshape(self.shape)

        new_values = values - values.shift(n, axis=axis)
        return [
            TimeDeltaBlock(new_values, placement=self.mgr_locs.indexer, ndim=self.ndim)
        ]

    def shift(self, periods, axis=0, fill_value=None):
        # TODO(EA2D) this is unnecessary if these blocks are backed by 2D EAs
        values = self.array_values()
        new_values = values.shift(periods, fill_value=fill_value, axis=axis)
        return self.make_block_same_class(new_values)

    def to_native_types(self, na_rep="NaT", **kwargs):
        """ convert to our native types format """
        arr = self.array_values()

        result = arr._format_native_types(na_rep=na_rep, **kwargs)
        return self.make_block(result)


class DatetimeBlock(DatetimeLikeBlockMixin):
    __slots__ = ()
    is_datetime = True

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
        if values.dtype != DT64NS_DTYPE:
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
            if copy:
                # this should be the only copy
                values = values.copy()
            values = DatetimeArray._simple_new(values.view("i8"), dtype=dtype)
            return self.make_block(values)

        # delegate
        return super().astype(dtype=dtype, copy=copy, errors=errors)

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            if isinstance(element, list) and len(element) == 0:
                # Following DatetimeArray._validate_setitem_value
                #  convention, we treat this as object-dtype
                #  (even though tipo is float64)
                return True

            elif self.is_datetimetz:
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

    def set_inplace(self, locs, values):
        """
        See Block.set.__doc__
        """
        values = conversion.ensure_datetime64ns(values, copy=False)

        self.values[locs] = values


class DatetimeTZBlock(ExtensionBlock, DatetimeBlock):
    """ implement a datetime64 block with a tz attribute """

    values: DatetimeArray

    __slots__ = ()
    is_datetimetz = True
    is_extension = True

    internal_values = Block.internal_values
    _can_hold_element = DatetimeBlock._can_hold_element
    to_native_types = DatetimeBlock.to_native_types
    diff = DatetimeBlock.diff
    fill_value = np.datetime64("NaT", "ns")

    array_values = ExtensionBlock.array_values

    @property
    def _holder(self):
        return DatetimeArray

    def _maybe_coerce_values(self, values):
        """
        Input validation for values passed to __init__. Ensure that
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
    def is_view(self) -> bool:
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

        # TODO(EA2D): reshape unnecessary with 2D EAs
        # Ensure that our shape is correct for DataFrame.
        # ExtensionArrays are always 1-D, even in a DataFrame when
        # the analogous NumPy-backed column would be a 2-D ndarray.
        return np.asarray(values).reshape(self.shape)

    def external_values(self):
        # NB: this is different from np.asarray(self.values), since that
        #  return an object-dtype ndarray of Timestamps.
        return np.asarray(self.values.astype("datetime64[ns]", copy=False))

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        # We support filling a DatetimeTZ with a `value` whose timezone
        # is different by coercing to object.
        if self._can_hold_element(value):
            return super().fillna(value, limit, inplace, downcast)

        # different timezones, or a non-tz
        return self.astype(object).fillna(
            value, limit=limit, inplace=inplace, downcast=downcast
        )

    def quantile(self, qs, interpolation="linear", axis=0):
        naive = self.values.view("M8[ns]")

        # TODO(EA2D): kludge for 2D block with 1D values
        naive = naive.reshape(self.shape)

        blk = self.make_block(naive)
        res_blk = blk.quantile(qs, interpolation=interpolation, axis=axis)

        # TODO(EA2D): ravel is kludge for 2D block with 1D values, assumes column-like
        aware = self._holder(res_blk.values.ravel(), dtype=self.dtype)
        return self.make_block_same_class(aware, ndim=res_blk.ndim)

    def _check_ndim(self, values, ndim):
        """
        ndim inference and validation.

        This is overriden by the DatetimeTZBlock to check the case of 2D
        data (values.ndim == 2), which should only be allowed if ndim is
        also 2.
        The case of 1D array is still allowed with both ndim of 1 or 2, as
        if the case for other EAs. Therefore, we are only checking
        `values.ndim > ndim` instead of `values.ndim != ndim` as for
        consolidated blocks.
        """
        if ndim is None:
            ndim = values.ndim

        if values.ndim > ndim:
            raise ValueError(
                "Wrong number of dimensions. "
                f"values.ndim != ndim [{values.ndim} != {ndim}]"
            )
        return ndim


class TimeDeltaBlock(DatetimeLikeBlockMixin, IntBlock):
    __slots__ = ()
    is_timedelta = True
    _can_hold_na = True
    is_numeric = False
    fill_value = np.timedelta64("NaT", "ns")

    def _maybe_coerce_values(self, values):
        if values.dtype != TD64NS_DTYPE:
            # non-nano we will convert to nano
            if values.dtype.kind != "m":
                # caller is responsible for ensuring timedelta64 dtype
                raise TypeError(values.dtype)  # pragma: no cover

            values = TimedeltaArray._from_sequence(values)._data
        if isinstance(values, TimedeltaArray):
            values = values._data
        assert isinstance(values, np.ndarray), type(values)
        return values

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
        # TODO(EA2D): if we operated on array_values, TDA.fillna would handle
        #  raising here.
        if is_integer(value):
            # Deprecation GH#24694, GH#19233
            raise TypeError(
                "Passing integers to fillna for timedelta64[ns] dtype is no "
                "longer supported.  To obtain the old behavior, pass "
                "`pd.Timedelta(seconds=n)` instead."
            )
        return super().fillna(value, **kwargs)


class BoolBlock(NumericBlock):
    __slots__ = ()
    is_bool = True
    _can_hold_na = False

    def _can_hold_element(self, element: Any) -> bool:
        tipo = maybe_infer_dtype_type(element)
        if tipo is not None:
            return issubclass(tipo.type, np.bool_)
        return isinstance(element, (bool, np.bool_))


class ObjectBlock(Block):
    __slots__ = ()
    is_object = True
    _can_hold_na = True

    def _maybe_coerce_values(self, values):
        if issubclass(values.dtype.type, str):
            values = np.array(values, dtype=object)
        return values

    @property
    def is_bool(self):
        """
        we can be a bool if we have only bool values but are of type
        object
        """
        return lib.is_bool_array(self.values.ravel("K"))

    def reduce(self, func, ignore_failures: bool = False) -> List[Block]:
        """
        For object-dtype, we operate column-wise.
        """
        assert self.ndim == 2

        values = self.values
        if len(values) > 1:
            # split_and_operate expects func with signature (mask, values, inplace)
            def mask_func(mask, values, inplace):
                if values.ndim == 1:
                    values = values.reshape(1, -1)
                return func(values)

            return self.split_and_operate(
                None, mask_func, False, ignore_failures=ignore_failures
            )

        try:
            res = func(values)
        except TypeError:
            if not ignore_failures:
                raise
            return []

        assert isinstance(res, np.ndarray)
        assert res.ndim == 1
        res = res.reshape(1, -1)
        return [self.make_block_same_class(res)]

    def convert(
        self,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
    ) -> List["Block"]:
        """
        attempt to cast any object types to better types return a copy of
        the block (if copy = True) by definition we ARE an ObjectBlock!!!!!
        """

        # operate column-by-column
        def f(mask, val, idx):
            shape = val.shape
            values = soft_convert_objects(
                val.ravel(),
                datetime=datetime,
                numeric=numeric,
                timedelta=timedelta,
                copy=copy,
            )
            if isinstance(values, np.ndarray):
                # TODO(EA2D): allow EA once reshape is supported
                values = values.reshape(shape)

            return values

        if self.ndim == 2:
            blocks = self.split_and_operate(None, f, False)
        else:
            values = f(None, self.values.ravel(), None)
            blocks = [self.make_block(values)]

        return blocks

    def _maybe_downcast(self, blocks: List["Block"], downcast=None) -> List["Block"]:

        if downcast is not None:
            return blocks

        # split and convert the blocks
        return extend_blocks([b.convert(datetime=True, numeric=False) for b in blocks])

    def _can_hold_element(self, element: Any) -> bool:
        return True

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        regex: bool = False,
    ) -> List["Block"]:
        # Note: the checks we do in NDFrame.replace ensure we never get
        #  here with listlike to_replace or value, as those cases
        #  go through _replace_list

        regex = _should_use_regex(regex, to_replace)

        if regex:
            return self._replace_regex(to_replace, value, inplace=inplace)
        else:
            return super().replace(to_replace, value, inplace=inplace, regex=False)


def _should_use_regex(regex: bool, to_replace: Any) -> bool:
    """
    Decide whether to treat `to_replace` as a regular expression.
    """
    if is_re(to_replace):
        regex = True

    regex = regex and is_re_compilable(to_replace)

    # Don't use regex if the pattern is empty.
    regex = regex and re.compile(to_replace).pattern != ""
    return regex


class CategoricalBlock(ExtensionBlock):
    __slots__ = ()

    def _replace_list(
        self,
        src_list: List[Any],
        dest_list: List[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> List["Block"]:
        if len(algos.unique(dest_list)) == 1:
            # We likely got here by tiling value inside NDFrame.replace,
            #  so un-tile here
            return self.replace(src_list, dest_list[0], inplace, regex)
        return super()._replace_list(src_list, dest_list, inplace, regex)

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        regex: bool = False,
    ) -> List["Block"]:
        inplace = validate_bool_kwarg(inplace, "inplace")
        result = self if inplace else self.copy()

        result.values.replace(to_replace, value, inplace=True)
        return [result]


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

    cls: Type[Block]

    if is_sparse(dtype):
        # Need this first(ish) so that Sparse[datetime] is sparse
        cls = ExtensionBlock
    elif is_categorical_dtype(values.dtype):
        cls = CategoricalBlock
    elif issubclass(vtype, np.datetime64):
        assert not is_datetime64tz_dtype(values.dtype)
        cls = DatetimeBlock
    elif is_datetime64tz_dtype(values.dtype):
        cls = DatetimeTZBlock
    elif is_interval_dtype(dtype) or is_period_dtype(dtype):
        cls = ObjectValuesExtensionBlock
    elif is_extension_array_dtype(values.dtype):
        # Note: need to be sure PandasArray is unwrapped before we get here
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
            # TODO(EA2D): special case not needed with 2D EAs
            values = np.atleast_2d(values)

    if isinstance(dtype, PandasDtype):
        dtype = dtype.numpy_dtype

    if klass is None:
        dtype = dtype or values.dtype
        klass = get_block_type(values, dtype)

    elif klass is DatetimeTZBlock and not is_datetime64tz_dtype(values.dtype):
        # TODO: This is no longer hit internally; does it need to be retained
        #  for e.g. pyarrow?
        values = DatetimeArray._simple_new(values, dtype=dtype)

    return klass(values, ndim=ndim, placement=placement)


# -----------------------------------------------------------------


def extend_blocks(result, blocks=None):
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


def _block_shape(values: ArrayLike, ndim: int = 1) -> ArrayLike:
    """ guarantee the shape of the values to be at least 1 d """
    if values.ndim < ndim:
        shape = values.shape
        if not is_extension_array_dtype(values.dtype):
            # TODO(EA2D): https://github.com/pandas-dev/pandas/issues/23023
            # block.shape is incorrect for "2D" ExtensionArrays
            # We can't, and don't need to, reshape.
            # error: "ExtensionArray" has no attribute "reshape"
            values = values.reshape(tuple((1,) + shape))  # type: ignore[attr-defined]
    return values


def safe_reshape(arr, new_shape: Shape):
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
    if not is_extension_array_dtype(arr.dtype):
        # Note: this will include TimedeltaArray and tz-naive DatetimeArray
        # TODO(EA2D): special case will be unnecessary with 2D EAs
        arr = np.asarray(arr).reshape(new_shape)
    return arr


def _putmask_smart(v: np.ndarray, mask: np.ndarray, n) -> np.ndarray:
    """
    Return a new ndarray, try to preserve dtype if possible.

    Parameters
    ----------
    v : np.ndarray
        `values`, updated in-place.
    mask : np.ndarray[bool]
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
        if not isna_compat(v, nn[0]):
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

    v = v.astype(dtype)

    return _putmask_preserve(v, n)


def _extract_bool_array(mask: ArrayLike) -> np.ndarray:
    """
    If we have a SparseArray or BooleanArray, convert it to ndarray[bool].
    """
    if isinstance(mask, ExtensionArray):
        # We could have BooleanArray, Sparse[bool], ...
        #  Except for BooleanArray, this is equivalent to just
        #  np.asarray(mask, dtype=bool)
        mask = mask.to_numpy(dtype=bool, na_value=False)

    assert isinstance(mask, np.ndarray), type(mask)
    assert mask.dtype == bool, mask.dtype
    return mask
