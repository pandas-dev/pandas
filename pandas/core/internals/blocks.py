from __future__ import annotations

import re
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    List,
    Optional,
    Tuple,
    Type,
    Union,
    cast,
)

import numpy as np

from pandas._libs import (
    Interval,
    Period,
    Timestamp,
    algos as libalgos,
    internals as libinternals,
    lib,
    writers,
)
from pandas._libs.internals import BlockPlacement
from pandas._libs.tslibs import conversion
from pandas._typing import (
    ArrayLike,
    Dtype,
    DtypeObj,
    Shape,
    final,
)
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import (
    astype_array_safe,
    can_hold_element,
    find_common_type,
    infer_dtype_from,
    maybe_downcast_numeric,
    maybe_downcast_to_dtype,
    maybe_upcast,
    sanitize_to_nanoseconds,
    soft_convert_objects,
)
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_list_like,
    is_object_dtype,
    is_sparse,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    ExtensionDtype,
    PandasDtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndex,
    ABCPandasArray,
    ABCSeries,
)
from pandas.core.dtypes.missing import (
    is_valid_na_for_dtype,
    isna,
    na_value_for_dtype,
)

import pandas.core.algorithms as algos
from pandas.core.array_algos.putmask import (
    extract_bool_array,
    putmask_inplace,
    putmask_smart,
    putmask_without_repeat,
    setitem_datetimelike_compat,
    validate_putmask,
)
from pandas.core.array_algos.quantile import quantile_compat
from pandas.core.array_algos.replace import (
    compare_or_regex_search,
    replace_regex,
    should_use_regex,
)
from pandas.core.array_algos.transforms import shift
from pandas.core.arrays import (
    Categorical,
    DatetimeArray,
    ExtensionArray,
    FloatingArray,
    IntegerArray,
    PandasArray,
)
from pandas.core.base import PandasObject
import pandas.core.common as com
from pandas.core.construction import (
    ensure_wrapped_if_datetimelike,
    extract_array,
)
from pandas.core.indexers import (
    check_setitem_lengths,
    is_empty_indexer,
    is_exact_shape_match,
    is_scalar_indexer,
)
import pandas.core.missing as missing

if TYPE_CHECKING:
    from pandas import (
        Float64Index,
        Index,
    )
    from pandas.core.arrays._mixins import NDArrayBackedExtensionArray

# comparison is faster than is_object_dtype

# error: Value of type variable "_DTypeScalar" of "dtype" cannot be "object"
_dtype_obj = np.dtype(object)  # type: ignore[type-var]


class Block(PandasObject):
    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas
    data structure

    Index-ignorant; let the container take care of that
    """

    values: Union[np.ndarray, ExtensionArray]

    __slots__ = ["_mgr_locs", "values", "ndim"]
    is_numeric = False
    is_bool = False
    is_object = False
    is_extension = False
    _can_hold_na = False
    _can_consolidate = True
    _validate_ndim = True

    @classmethod
    def _simple_new(
        cls, values: ArrayLike, placement: BlockPlacement, ndim: int
    ) -> Block:
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
        self.ndim = ndim
        self.mgr_locs = placement
        self.values = self._maybe_coerce_values(values)

    @classmethod
    def _maybe_coerce_values(cls, values):
        """
        Ensure we have correctly-typed values.

        Parameters
        ----------
        values : np.ndarray or ExtensionArray

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        return values

    @property
    def _holder(self):
        """
        The array-like that can hold the underlying values.

        None for 'Block', overridden by subclasses that don't
        use an ndarray.
        """
        return None

    @final
    @property
    def _consolidate_key(self):
        return self._can_consolidate, self.dtype.name

    @property
    def is_view(self) -> bool:
        """ return a boolean if I am possibly a view """
        values = self.values
        values = cast(np.ndarray, values)
        return values.base is not None

    @final
    @property
    def is_categorical(self) -> bool:
        return self._holder is Categorical

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
        # error: Argument 1 to "PandasArray" has incompatible type "Union[ndarray,
        # ExtensionArray]"; expected "Union[ndarray, PandasArray]"
        return PandasArray(self.values)  # type: ignore[arg-type]

    def get_values(self, dtype: Optional[DtypeObj] = None) -> np.ndarray:
        """
        return an internal format, currently just the ndarray
        this is often overridden to handle to_dense like operations
        """
        if dtype == _dtype_obj:
            return self.values.astype(_dtype_obj)
        # error: Incompatible return value type (got "Union[ndarray, ExtensionArray]",
        # expected "ndarray")
        return self.values  # type: ignore[return-value]

    @final
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

    @final
    def make_block(self, values, placement=None) -> Block:
        """
        Create a new block, with type inference propagate any values that are
        not specified
        """
        if placement is None:
            placement = self.mgr_locs
        if self.is_extension:
            values = ensure_block_shape(values, ndim=self.ndim)

        return new_block(values, placement=placement, ndim=self.ndim)

    @final
    def make_block_same_class(self, values, placement=None) -> Block:
        """ Wrap given values in a block of same type as self. """
        if placement is None:
            placement = self.mgr_locs
        return type(self)(values, placement=placement, ndim=self.ndim)

    @final
    def __repr__(self) -> str:
        # don't want to print out all of the items here
        name = type(self).__name__
        if self.ndim == 1:
            result = f"{name}: {len(self)} dtype: {self.dtype}"
        else:

            shape = " x ".join(str(s) for s in self.shape)
            result = f"{name}: {self.mgr_locs.indexer}, {shape}, dtype: {self.dtype}"

        return result

    @final
    def __len__(self) -> int:
        return len(self.values)

    @final
    def __getstate__(self):
        return self.mgr_locs.indexer, self.values

    @final
    def __setstate__(self, state):
        self.mgr_locs = libinternals.BlockPlacement(state[0])
        self.values = extract_array(state[1], extract_numpy=True)
        self.ndim = self.values.ndim

    def _slice(self, slicer):
        """ return a slice of my values """

        return self.values[slicer]

    @final
    def getitem_block(self, slicer, new_mgr_locs=None) -> Block:
        """
        Perform __getitem__-like, return result as block.

        Only supports slices that preserve dimensionality.
        """
        if new_mgr_locs is None:
            axis0_slicer = slicer[0] if isinstance(slicer, tuple) else slicer
            new_mgr_locs = self.mgr_locs[axis0_slicer]
        elif not isinstance(new_mgr_locs, BlockPlacement):
            new_mgr_locs = BlockPlacement(new_mgr_locs)

        new_values = self._slice(slicer)

        if new_values.ndim != self.values.ndim:
            raise ValueError("Only same dim slicing is allowed")

        return type(self)._simple_new(new_values, new_mgr_locs, self.ndim)

    @property
    def shape(self) -> Shape:
        return self.values.shape

    @final
    @property
    def dtype(self) -> DtypeObj:
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

    @final
    def delete(self, loc) -> None:
        """
        Delete given loc(-s) from block in-place.
        """
        self.values = np.delete(self.values, loc, 0)
        self.mgr_locs = self.mgr_locs.delete(loc)

    @final
    def apply(self, func, **kwargs) -> List[Block]:
        """
        apply the function to my values; return a block if we are not
        one
        """
        with np.errstate(all="ignore"):
            result = func(self.values, **kwargs)

        return self._split_op_result(result)

    def reduce(self, func, ignore_failures: bool = False) -> List[Block]:
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

    @final
    def _split_op_result(self, result) -> List[Block]:
        # See also: split_and_operate
        if is_extension_array_dtype(result) and result.ndim > 1:
            # TODO(EA2D): unnecessary with 2D EAs
            # if we get a 2D ExtensionArray, we need to split it into 1D pieces
            nbs = []
            for i, loc in enumerate(self.mgr_locs):
                vals = result[i]
                block = self.make_block(values=vals, placement=loc)
                nbs.append(block)
            return nbs

        if not isinstance(result, Block):
            result = self.make_block(result)

        return [result]

    def fillna(
        self, value, limit=None, inplace: bool = False, downcast=None
    ) -> List[Block]:
        """
        fillna on the block with the value. If we fail, then convert to
        ObjectBlock and try again
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        mask = isna(self.values)
        mask, noop = validate_putmask(self.values, mask)

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
            putmask_inplace(nb.values, mask, value)
            return nb._maybe_downcast([nb], downcast)

        if noop:
            # we can't process the value, but nothing to do
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

    @final
    def _split(self) -> List[Block]:
        """
        Split a block into a list of single-column blocks.
        """
        assert self.ndim == 2

        new_blocks = []
        for i, ref_loc in enumerate(self.mgr_locs):
            vals = self.values[slice(i, i + 1)]

            nb = self.make_block(vals, BlockPlacement(ref_loc))
            new_blocks.append(nb)
        return new_blocks

    @final
    def split_and_operate(
        self, mask, f, inplace: bool, ignore_failures: bool = False
    ) -> List[Block]:
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
                nv = ensure_block_shape(nv, ndim=self.ndim)
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

    def _maybe_downcast(self, blocks: List[Block], downcast=None) -> List[Block]:

        # no need to downcast our float
        # unless indicated
        if downcast is None and self.dtype.kind in ["f", "m", "M"]:
            # TODO: complex?  more generally, self._can_hold_na?
            return blocks

        return extend_blocks([b.downcast(downcast) for b in blocks])

    def downcast(self, dtypes=None) -> List[Block]:
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

    @final
    def astype(self, dtype, copy: bool = False, errors: str = "raise"):
        """
        Coerce to the new dtype.

        Parameters
        ----------
        dtype : str, dtype convertible
        copy : bool, default False
            copy if indicated
        errors : str, {'raise', 'ignore'}, default 'raise'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object

        Returns
        -------
        Block
        """
        values = self.values
        if values.dtype.kind in ["m", "M"]:
            values = self.array_values()

        new_values = astype_array_safe(values, dtype, copy=copy, errors=errors)

        newb = self.make_block(new_values)
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
    ) -> List[Block]:
        """
        attempt to coerce any object types to better types return a copy
        of the block (if copy = True) by definition we are not an ObjectBlock
        here!
        """
        return [self.copy()] if copy else [self]

    def _can_hold_element(self, element: Any) -> bool:
        """ require the same dtype as ourselves """
        raise NotImplementedError("Implemented on subclasses")

    @final
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
    @final
    def copy(self, deep: bool = True):
        """ copy constructor """
        values = self.values
        if deep:
            values = values.copy()
        return self.make_block_same_class(values)

    # ---------------------------------------------------------------------
    # Replace

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        regex: bool = False,
    ) -> List[Block]:
        """
        replace the to_replace value with value, possible to create new
        blocks here this is just a call to putmask. regex is not used here.
        It is used in ObjectBlocks.  It is here for API compatibility.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")

        if not self._can_hold_element(to_replace):
            # We cannot hold `to_replace`, so we know immediately that
            #  replacing it is a no-op.
            # Note: If to_replace were a list, NDFrame.replace would call
            #  replace_list instead of replace.
            return [self] if inplace else [self.copy()]

        values = self.values

        mask = missing.mask_missing(values, to_replace)
        if not mask.any():
            # Note: we get here with test_replace_extension_other incorrectly
            #  bc _can_hold_element is incorrect.
            return [self] if inplace else [self.copy()]

        if not self._can_hold_element(value):
            if self.ndim == 2 and self.shape[0] > 1:
                # split so that we only upcast where necessary
                nbs = self._split()
                res_blocks = extend_blocks(
                    [
                        blk.replace(to_replace, value, inplace=inplace, regex=regex)
                        for blk in nbs
                    ]
                )
                return res_blocks

            blk = self.coerce_to_target_dtype(value)
            return blk.replace(
                to_replace=to_replace,
                value=value,
                inplace=True,
                regex=regex,
            )

        blk = self if inplace else self.copy()
        putmask_inplace(blk.values, mask, value)
        blocks = blk.convert(numeric=False, copy=False)
        return blocks

    @final
    def _replace_regex(
        self,
        to_replace,
        value,
        inplace: bool = False,
        convert: bool = True,
        mask=None,
    ) -> List[Block]:
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
        return [block]

    @final
    def _replace_list(
        self,
        src_list: List[Any],
        dest_list: List[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> List[Block]:
        """
        See BlockManager._replace_list docstring.
        """
        # TODO: dont special-case Categorical
        if self.is_categorical and len(algos.unique(dest_list)) == 1:
            # We likely got here by tiling value inside NDFrame.replace,
            #  so un-tile here
            return self.replace(src_list, dest_list[0], inplace, regex)

        # Exclude anything that we know we won't contain
        pairs = [
            (x, y) for x, y in zip(src_list, dest_list) if self._can_hold_element(x)
        ]
        if not len(pairs):
            # shortcut, nothing to replace
            return [self] if inplace else [self.copy()]

        src_len = len(pairs) - 1

        if self.is_object:
            # Calculate the mask once, prior to the call of comp
            # in order to avoid repeating the same computations
            mask = ~isna(self.values)
            masks = [
                compare_or_regex_search(self.values, s[0], regex=regex, mask=mask)
                for s in pairs
            ]
        else:
            # GH#38086 faster if we know we dont need to check for regex
            masks = [missing.mask_missing(self.values, s[0]) for s in pairs]

        # error: Argument 1 to "extract_bool_array" has incompatible type
        # "Union[ExtensionArray, ndarray, bool]"; expected "Union[ExtensionArray,
        # ndarray]"
        masks = [extract_bool_array(x) for x in masks]  # type: ignore[arg-type]

        rb = [self if inplace else self.copy()]
        for i, (src, dest) in enumerate(pairs):
            new_rb: List["Block"] = []
            for blk in rb:
                m = masks[i]
                convert = i == src_len  # only convert once at the end
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

    @final
    def _replace_coerce(
        self,
        to_replace,
        value,
        mask: np.ndarray,
        inplace: bool = True,
        regex: bool = False,
    ) -> List[Block]:
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
                putmask_inplace(nb.values, mask, value)
                return [nb]
            else:
                regex = should_use_regex(regex, to_replace)
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

    # ---------------------------------------------------------------------

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
        if not self._can_hold_element(value):
            # current dtype cannot store value, coerce to common dtype
            return self.coerce_to_target_dtype(value).setitem(indexer, value)

        if self.dtype.kind in ["m", "M"]:
            arr = self.array_values().T
            arr[indexer] = value
            return self

        # value must be storable at this moment
        if is_extension_array_dtype(getattr(value, "dtype", None)):
            # We need to be careful not to allow through strings that
            #  can be parsed to EADtypes
            is_ea_value = True
            arr_value = value
        else:
            is_ea_value = False
            arr_value = np.asarray(value)

        if transpose:
            values = values.T

        # length checking
        check_setitem_lengths(indexer, value, values)
        exact_match = is_exact_shape_match(values, arr_value)

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

        elif exact_match and is_ea_value:
            # GH#32395 if we're going to replace the values entirely, just
            #  substitute in the new array
            if not self.is_object and isinstance(value, (IntegerArray, FloatingArray)):
                values[indexer] = value.to_numpy(value.dtype.numpy_dtype)
            else:
                values[indexer] = np.asarray(value)

        # if we are an exact match (ex-broadcasting),
        # then use the resultant dtype
        elif exact_match:
            # We are setting _all_ of the array's values, so can cast to new dtype
            values[indexer] = value

        elif is_ea_value:
            # GH#38952
            if values.ndim == 1:
                values[indexer] = value
            else:
                # TODO(EA2D): special case not needed with 2D EA
                values[indexer] = value.to_numpy(values.dtype).reshape(-1, 1)

        else:
            # error: Argument 1 to "setitem_datetimelike_compat" has incompatible type
            # "Union[ndarray, ExtensionArray]"; expected "ndarray"
            value = setitem_datetimelike_compat(
                values, len(values[indexer]), value  # type: ignore[arg-type]
            )
            values[indexer] = value

        if transpose:
            values = values.T
        block = self.make_block(values)
        return block

    def putmask(self, mask, new) -> List[Block]:
        """
        putmask the data to the block; it is possible that we may create a
        new dtype of block

        Return the resulting block(s).

        Parameters
        ----------
        mask : np.ndarray[bool], SparseArray[bool], or BooleanArray
        new : a ndarray/object

        Returns
        -------
        List[Block]
        """
        orig_mask = mask
        mask, noop = validate_putmask(self.values.T, mask)
        assert not isinstance(new, (ABCIndex, ABCSeries, ABCDataFrame))

        # if we are passed a scalar None, convert it here
        if not self.is_object and is_valid_na_for_dtype(new, self.dtype):
            new = self.fill_value

        if self._can_hold_element(new):

            # error: Argument 1 to "putmask_without_repeat" has incompatible type
            # "Union[ndarray, ExtensionArray]"; expected "ndarray"
            putmask_without_repeat(self.values.T, mask, new)  # type: ignore[arg-type]
            return [self]

        elif noop:
            return [self]

        dtype, _ = infer_dtype_from(new)
        if dtype.kind in ["m", "M"]:
            # using putmask with object dtype will incorrectly cast to object
            # Having excluded self._can_hold_element, we know we cannot operate
            #  in-place, so we are safe using `where`
            return self.where(new, ~mask)

        elif self.ndim == 1 or self.shape[0] == 1:
            # no need to split columns

            # error: Argument 1 to "putmask_smart" has incompatible type "Union[ndarray,
            # ExtensionArray]"; expected "ndarray"
            nv = putmask_smart(self.values.T, mask, new).T  # type: ignore[arg-type]
            return [self.make_block(nv)]

        else:
            is_array = isinstance(new, np.ndarray)

            res_blocks = []
            nbs = self._split()
            for i, nb in enumerate(nbs):
                n = new
                if is_array:
                    # we have a different value per-column
                    n = new[:, i : i + 1]

                submask = orig_mask[:, i : i + 1]
                rbs = nb.putmask(submask, n)
                res_blocks.extend(rbs)
            return res_blocks

    @final
    def coerce_to_target_dtype(self, other) -> Block:
        """
        coerce the current block to a dtype compat for other
        we will return a block, possibly object, and not raise

        we can also safely try to coerce to the same dtype
        and will receive the same block
        """
        # if we cannot then coerce to object
        dtype, _ = infer_dtype_from(other, pandas_dtype=True)

        new_dtype = find_common_type([self.dtype, dtype])

        return self.astype(new_dtype, copy=False)

    def interpolate(
        self,
        method: str = "pad",
        axis: int = 0,
        index: Optional[Index] = None,
        inplace: bool = False,
        limit: Optional[int] = None,
        limit_direction: str = "forward",
        limit_area: Optional[str] = None,
        fill_value: Optional[Any] = None,
        coerce: bool = False,
        downcast: Optional[str] = None,
        **kwargs,
    ) -> List[Block]:

        inplace = validate_bool_kwarg(inplace, "inplace")

        if not self._can_hold_na:
            # If there are no NAs, then interpolate is a no-op
            return [self] if inplace else [self.copy()]

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

    @final
    def _interpolate_with_fill(
        self,
        method: str = "pad",
        axis: int = 0,
        inplace: bool = False,
        limit: Optional[int] = None,
        limit_area: Optional[str] = None,
        downcast: Optional[str] = None,
    ) -> List[Block]:
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

        blocks = [self.make_block_same_class(values)]
        return self._maybe_downcast(blocks, downcast)

    @final
    def _interpolate(
        self,
        method: str,
        index: Index,
        fill_value: Optional[Any] = None,
        axis: int = 0,
        limit: Optional[int] = None,
        limit_direction: str = "forward",
        limit_area: Optional[str] = None,
        inplace: bool = False,
        downcast: Optional[str] = None,
        **kwargs,
    ) -> List[Block]:
        """ interpolate using scipy wrappers """
        inplace = validate_bool_kwarg(inplace, "inplace")
        data = self.values if inplace else self.values.copy()

        # only deal with floats
        if self.dtype.kind != "f":
            # bc we already checked that can_hold_na, we dont have int dtype here
            return [self]

        if is_valid_na_for_dtype(fill_value, self.dtype):
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

    def take_nd(
        self, indexer, axis: int, new_mgr_locs=None, fill_value=lib.no_default
    ) -> Block:
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

    def diff(self, n: int, axis: int = 1) -> List[Block]:
        """ return block for the diff of the values """
        new_values = algos.diff(self.values, n, axis=axis, stacklevel=7)
        return [self.make_block(values=new_values)]

    def shift(self, periods: int, axis: int = 0, fill_value: Any = None) -> List[Block]:
        """ shift the block by periods, possibly upcast """
        # convert integer to float if necessary. need to do a lot more than
        # that, handle boolean etc also

        # error: Argument 1 to "maybe_upcast" has incompatible type "Union[ndarray,
        # ExtensionArray]"; expected "ndarray"
        new_values, fill_value = maybe_upcast(
            self.values, fill_value  # type: ignore[arg-type]
        )

        new_values = shift(new_values, periods, axis, fill_value)

        return [self.make_block(new_values)]

    def where(self, other, cond, errors="raise", axis: int = 0) -> List[Block]:
        """
        evaluate the block; return result block(s) from the result

        Parameters
        ----------
        other : a ndarray/object
        cond : np.ndarray[bool], SparseArray[bool], or BooleanArray
        errors : str, {'raise', 'ignore'}, default 'raise'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object
        axis : int, default 0

        Returns
        -------
        List[Block]
        """
        import pandas.core.computation.expressions as expressions

        assert not isinstance(other, (ABCIndex, ABCSeries, ABCDataFrame))

        assert errors in ["raise", "ignore"]
        transpose = self.ndim == 2

        values = self.values
        orig_other = other
        if transpose:
            values = values.T

        icond, noop = validate_putmask(values, ~cond)

        if is_valid_na_for_dtype(other, self.dtype) and not self.is_object:
            other = self.fill_value

        if noop:
            # TODO: avoid the downcasting at the end in this case?
            result = values
        else:
            # see if we can operate on the entire block, or need item-by-item
            # or if we are a single block (ndim == 1)
            if not self._can_hold_element(other):
                # we cannot coerce, return a compat dtype
                # we are explicitly ignoring errors
                block = self.coerce_to_target_dtype(other)
                blocks = block.where(orig_other, cond, errors=errors, axis=axis)
                return self._maybe_downcast(blocks, "infer")

            # error: Argument 1 to "setitem_datetimelike_compat" has incompatible type
            # "Union[ndarray, ExtensionArray]"; expected "ndarray"
            # error: Argument 2 to "setitem_datetimelike_compat" has incompatible type
            # "number[Any]"; expected "int"
            alt = setitem_datetimelike_compat(
                values, icond.sum(), other  # type: ignore[arg-type]
            )
            if alt is not other:
                result = values.copy()
                np.putmask(result, icond, alt)
            else:
                # By the time we get here, we should have all Series/Index
                #  args extracted to ndarray
                result = expressions.where(~icond, values, other)

        if self._can_hold_na or self.ndim == 1:

            if transpose:
                result = result.T

            return [self.make_block(result)]

        # might need to separate out blocks
        cond = ~icond
        axis = cond.ndim - 1
        cond = cond.swapaxes(axis, 0)
        mask = np.array([cond[i].all() for i in range(cond.shape[0])], dtype=bool)

        result_blocks: List[Block] = []
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

        blocks = [new_block(new_values, placement=new_placement, ndim=2)]
        return blocks, mask

    def quantile(
        self, qs: Float64Index, interpolation="linear", axis: int = 0
    ) -> Block:
        """
        compute the quantiles of the

        Parameters
        ----------
        qs : Float64Index
            List of the quantiles to be computed.
        interpolation : str, default 'linear'
            Type of interpolation.
        axis : int, default 0
            Axis to compute.

        Returns
        -------
        Block
        """
        # We should always have ndim == 2 because Series dispatches to DataFrame
        assert self.ndim == 2
        assert axis == 1  # only ever called this way
        assert is_list_like(qs)  # caller is responsible for this

        result = quantile_compat(self.values, qs, interpolation, axis)

        return new_block(result, placement=self.mgr_locs, ndim=2)


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

    @property
    def shape(self) -> Shape:
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
                # error: Invalid index type "List[Any]" for "ExtensionArray"; expected
                # type "Union[int, slice, ndarray]"
                return self.values[[loc]]  # type: ignore[index]
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

    def putmask(self, mask, new) -> List[Block]:
        """
        See Block.putmask.__doc__
        """
        mask = extract_bool_array(mask)

        new_values = self.values

        if isinstance(new, (np.ndarray, ExtensionArray)) and len(new) == len(mask):
            new = new[mask]

        if mask.ndim == new_values.ndim + 1:
            # TODO(EA2D): unnecessary with 2D EAs
            mask = mask.reshape(new_values.shape)

        new_values[mask] = new
        return [self.make_block(values=new_values)]

    @classmethod
    def _maybe_coerce_values(cls, values):
        """
        Unbox to an extension array.

        This will unbox an ExtensionArray stored in an Index or Series.
        ExtensionArrays pass through. No dtype coercion is done.

        Parameters
        ----------
        values : np.ndarray or ExtensionArray

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
            # This is only relevant for DatetimeTZBlock, ObjectValuesExtensionBlock,
            #  which has a non-trivial `_can_hold_element`.
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

    def get_values(self, dtype: Optional[DtypeObj] = None) -> np.ndarray:
        # ExtensionArrays must be iterable, so this works.
        # TODO(EA2D): reshape not needed with 2D EAs
        return np.asarray(self.values).reshape(self.shape)

    def array_values(self) -> ExtensionArray:
        return self.values

    def to_native_types(self, na_rep="nan", quoting=None, **kwargs):
        """override to use ExtensionArray astype for the conversion"""
        values = self.values
        mask = isna(values)

        # error: Incompatible types in assignment (expression has type "ndarray",
        # variable has type "ExtensionArray")
        values = np.asarray(values.astype(object))  # type: ignore[assignment]
        values[mask] = na_rep

        # TODO(EA2D): reshape not needed with 2D EAs
        # we are expected to return a 2-d ndarray
        return self.make_block(values)

    def take_nd(
        self, indexer, axis: int = 0, new_mgr_locs=None, fill_value=lib.no_default
    ) -> Block:
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
            # TODO(EA2D): won't be necessary with 2D EAs
            slicer = (slicer, slice(None))

        if isinstance(slicer, tuple) and len(slicer) == 2:
            first = slicer[0]
            if not isinstance(first, slice):
                raise AssertionError(
                    "invalid slicing for a 1-ndim ExtensionArray", first
                )
            # GH#32959 only full-slicers along fake-dim0 are valid
            # TODO(EA2D): won't be necessary with 2D EAs
            new_locs = self.mgr_locs[first]
            if len(new_locs):
                # effectively slice(None)
                slicer = slicer[1]
            else:
                raise AssertionError(
                    "invalid slicing for a 1-ndim ExtensionArray", slicer
                )

        return self.values[slicer]

    def fillna(
        self, value, limit=None, inplace: bool = False, downcast=None
    ) -> List[Block]:
        values = self.values.fillna(value=value, limit=limit)
        return [self.make_block_same_class(values=values)]

    def interpolate(
        self, method="pad", axis=0, inplace=False, limit=None, fill_value=None, **kwargs
    ):
        new_values = self.values.fillna(value=fill_value, method=method, limit=limit)
        return self.make_block_same_class(new_values)

    def diff(self, n: int, axis: int = 1) -> List[Block]:
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

    def shift(self, periods: int, axis: int = 0, fill_value: Any = None) -> List[Block]:
        """
        Shift the block by `periods`.

        Dispatches to underlying ExtensionArray and re-boxes in an
        ExtensionBlock.
        """
        new_values = self.values.shift(periods=periods, fill_value=fill_value)
        return [self.make_block_same_class(new_values)]

    def where(self, other, cond, errors="raise", axis: int = 0) -> List[Block]:

        cond = extract_bool_array(cond)
        assert not isinstance(other, (ABCIndex, ABCSeries, ABCDataFrame))

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

            # error: Item "dtype[Any]" of "Union[dtype[Any], ExtensionDtype]" has no
            # attribute "na_value"
            other = self.dtype.na_value  # type: ignore[union-attr]

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

        return [self.make_block_same_class(result)]

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


class HybridMixin:
    """
    Mixin for Blocks backed (maybe indirectly) by ExtensionArrays.
    """

    array_values: Callable

    def _can_hold_element(self, element: Any) -> bool:
        values = self.array_values()

        try:
            values._validate_setitem_value(element)
            return True
        except (ValueError, TypeError):
            return False


class ObjectValuesExtensionBlock(HybridMixin, ExtensionBlock):
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

    def _can_hold_element(self, element: Any) -> bool:
        element = extract_array(element, extract_numpy=True)
        if isinstance(element, (IntegerArray, FloatingArray)):
            if element._mask.any():
                return False
        # error: Argument 1 to "can_hold_element" has incompatible type
        # "Union[dtype[Any], ExtensionDtype]"; expected "dtype[Any]"
        return can_hold_element(self.dtype, element)  # type: ignore[arg-type]

    @property
    def _can_hold_na(self):
        return self.dtype.kind not in ["b", "i", "u"]

    @property
    def is_bool(self):
        return self.dtype.kind == "b"


class FloatBlock(NumericBlock):
    __slots__ = ()

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


class NDArrayBackedExtensionBlock(HybridMixin, Block):
    """
    Block backed by an NDArrayBackedExtensionArray
    """

    def internal_values(self):
        # Override to return DatetimeArray and TimedeltaArray
        return self.array_values()

    def get_values(self, dtype: Optional[DtypeObj] = None) -> np.ndarray:
        """
        return object dtype as boxed values, such as Timestamps/Timedelta
        """
        values = self.array_values()
        if is_object_dtype(dtype):
            # DTA/TDA constructor and astype can handle 2D
            values = values.astype(object)
        # TODO(EA2D): reshape not needed with 2D EAs
        return np.asarray(values).reshape(self.shape)

    def iget(self, key):
        # GH#31649 we need to wrap scalars in Timestamp/Timedelta
        # TODO(EA2D): this can be removed if we ever have 2D EA
        return self.array_values().reshape(self.shape)[key]

    def putmask(self, mask, new) -> List[Block]:
        mask = extract_bool_array(mask)

        if not self._can_hold_element(new):
            return self.astype(object).putmask(mask, new)

        # TODO(EA2D): reshape unnecessary with 2D EAs
        arr = self.array_values().reshape(self.shape)
        arr = cast("NDArrayBackedExtensionArray", arr)
        arr.T.putmask(mask, new)
        return [self]

    def where(self, other, cond, errors="raise", axis: int = 0) -> List[Block]:
        # TODO(EA2D): reshape unnecessary with 2D EAs
        arr = self.array_values().reshape(self.shape)

        cond = extract_bool_array(cond)

        try:
            res_values = arr.T.where(cond, other).T
        except (ValueError, TypeError):
            return super().where(other, cond, errors=errors, axis=axis)

        # TODO(EA2D): reshape not needed with 2D EAs
        res_values = res_values.reshape(self.values.shape)
        nb = self.make_block_same_class(res_values)
        return [nb]

    def diff(self, n: int, axis: int = 0) -> List[Block]:
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
        return [self.make_block(new_values)]

    def shift(self, periods: int, axis: int = 0, fill_value: Any = None) -> List[Block]:
        # TODO(EA2D) this is unnecessary if these blocks are backed by 2D EAs
        values = self.array_values().reshape(self.shape)
        new_values = values.shift(periods, fill_value=fill_value, axis=axis)
        return [self.make_block_same_class(new_values)]

    def fillna(
        self, value, limit=None, inplace: bool = False, downcast=None
    ) -> List[Block]:

        if not self._can_hold_element(value) and self.dtype.kind != "m":
            # We support filling a DatetimeTZ with a `value` whose timezone
            #  is different by coercing to object.
            # TODO: don't special-case td64
            return self.astype(object).fillna(value, limit, inplace, downcast)

        values = self.array_values()
        values = values if inplace else values.copy()
        new_values = values.fillna(value=value, limit=limit)
        return [self.make_block_same_class(values=new_values)]


class DatetimeLikeBlockMixin(NDArrayBackedExtensionBlock):
    """Mixin class for DatetimeBlock, DatetimeTZBlock, and TimedeltaBlock."""

    is_numeric = False
    _can_hold_na = True

    @classmethod
    def _maybe_coerce_values(cls, values):
        """
        Input validation for values passed to __init__. Ensure that
        we have nanosecond datetime64/timedelta64, coercing if necessary.

        Parameters
        ----------
        values : np.ndarray or ExtensionArray
            Must be convertible to datetime64/timedelta64

        Returns
        -------
        values : ndarray[datetime64ns/timedelta64ns]
        """
        values = extract_array(values, extract_numpy=True)
        if isinstance(values, np.ndarray):
            values = sanitize_to_nanoseconds(values)
        elif isinstance(values.dtype, np.dtype):
            # i.e. not datetime64tz
            values = values._data

        return values

    def array_values(self):
        return ensure_wrapped_if_datetimelike(self.values)

    @property
    def _holder(self):
        return type(self.array_values())

    @property
    def fill_value(self):
        return na_value_for_dtype(self.dtype)

    def to_native_types(self, na_rep="NaT", **kwargs):
        """ convert to our native types format """
        arr = self.array_values()

        result = arr._format_native_types(na_rep=na_rep, **kwargs)
        return self.make_block(result)


class DatetimeBlock(DatetimeLikeBlockMixin):
    __slots__ = ()

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
    is_extension = True
    _can_hold_na = True
    is_numeric = False

    internal_values = Block.internal_values
    _can_hold_element = DatetimeBlock._can_hold_element
    to_native_types = DatetimeBlock.to_native_types
    diff = DatetimeBlock.diff
    where = DatetimeBlock.where
    putmask = DatetimeLikeBlockMixin.putmask
    fillna = DatetimeLikeBlockMixin.fillna

    array_values = ExtensionBlock.array_values

    @property
    def is_view(self) -> bool:
        """ return a boolean if I am possibly a view """
        # check the ndarray values of the DatetimeIndex values
        return self.values._data.base is not None

    def external_values(self):
        # NB: this is different from np.asarray(self.values), since that
        #  return an object-dtype ndarray of Timestamps.
        # Avoid FutureWarning in .astype in casting from dt64tz to dt64
        return self.values._data


class TimeDeltaBlock(DatetimeLikeBlockMixin):
    __slots__ = ()


class ObjectBlock(Block):
    __slots__ = ()
    is_object = True
    _can_hold_na = True

    @classmethod
    def _maybe_coerce_values(cls, values):
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
    ) -> List[Block]:
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

    def _maybe_downcast(self, blocks: List[Block], downcast=None) -> List[Block]:

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
    ) -> List[Block]:
        # Note: the checks we do in NDFrame.replace ensure we never get
        #  here with listlike to_replace or value, as those cases
        #  go through _replace_list

        regex = should_use_regex(regex, to_replace)

        if regex:
            return self._replace_regex(to_replace, value, inplace=inplace)
        else:
            return super().replace(to_replace, value, inplace=inplace, regex=False)


class CategoricalBlock(ExtensionBlock):
    __slots__ = ()

    def replace(
        self,
        to_replace,
        value,
        inplace: bool = False,
        regex: bool = False,
    ) -> List[Block]:
        inplace = validate_bool_kwarg(inplace, "inplace")
        result = self if inplace else self.copy()

        result.values.replace(to_replace, value, inplace=True)
        return [result]


# -----------------------------------------------------------------
# Constructor Helpers


def get_block_type(values, dtype: Optional[Dtype] = None):
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
    # We use vtype and kind checks because they are much more performant
    #  than is_foo_dtype
    dtype = cast(np.dtype, pandas_dtype(dtype) if dtype else values.dtype)
    vtype = dtype.type
    kind = dtype.kind

    cls: Type[Block]

    if is_sparse(dtype):
        # Need this first(ish) so that Sparse[datetime] is sparse
        cls = ExtensionBlock
    elif isinstance(dtype, CategoricalDtype):
        cls = CategoricalBlock
    elif vtype is Timestamp:
        cls = DatetimeTZBlock
    elif vtype is Interval or vtype is Period:
        cls = ObjectValuesExtensionBlock
    elif isinstance(dtype, ExtensionDtype):
        # Note: need to be sure PandasArray is unwrapped before we get here
        cls = ExtensionBlock

    elif kind == "M":
        cls = DatetimeBlock
    elif kind == "m":
        cls = TimeDeltaBlock
    elif kind == "f":
        cls = FloatBlock
    elif kind in ["c", "i", "u", "b"]:
        cls = NumericBlock
    else:
        cls = ObjectBlock
    return cls


def new_block(values, placement, *, ndim: int, klass=None) -> Block:

    if not isinstance(placement, BlockPlacement):
        placement = BlockPlacement(placement)

    values, _ = extract_pandas_array(values, None, ndim)
    check_ndim(values, placement, ndim)

    if klass is None:
        klass = get_block_type(values, values.dtype)

    return klass(values, ndim=ndim, placement=placement)


def check_ndim(values, placement: BlockPlacement, ndim: int):
    """
    ndim inference and validation.

    Validates that values.ndim and ndim are consistent.
    Validates that len(values) and len(placement) are consistent.

    Parameters
    ----------
    values : array-like
    placement : BlockPlacement
    ndim : int

    Raises
    ------
    ValueError : the number of dimensions do not match
    """

    if values.ndim > ndim:
        # Check for both np.ndarray and ExtensionArray
        raise ValueError(
            "Wrong number of dimensions. "
            f"values.ndim > ndim [{values.ndim} > {ndim}]"
        )

    elif isinstance(values.dtype, np.dtype):
        # TODO(EA2D): special case not needed with 2D EAs
        if values.ndim != ndim:
            raise ValueError(
                "Wrong number of dimensions. "
                f"values.ndim != ndim [{values.ndim} != {ndim}]"
            )
        if len(placement) != len(values):
            raise ValueError(
                f"Wrong number of items passed {len(values)}, "
                f"placement implies {len(placement)}"
            )
    elif ndim == 2 and len(placement) != 1:
        # TODO(EA2D): special case unnecessary with 2D EAs
        raise AssertionError("block.size != values.size")


def extract_pandas_array(
    values: ArrayLike, dtype: Optional[DtypeObj], ndim: int
) -> Tuple[ArrayLike, Optional[DtypeObj]]:
    """
    Ensure that we don't allow PandasArray / PandasDtype in internals.
    """
    # For now, blocks should be backed by ndarrays when possible.
    if isinstance(values, ABCPandasArray):
        values = values.to_numpy()
        if ndim and ndim > 1:
            # TODO(EA2D): special case not needed with 2D EAs
            values = np.atleast_2d(values)

    if isinstance(dtype, PandasDtype):
        dtype = dtype.numpy_dtype

    return values, dtype


# -----------------------------------------------------------------


def extend_blocks(result, blocks=None) -> List[Block]:
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


def ensure_block_shape(values: ArrayLike, ndim: int = 1) -> ArrayLike:
    """
    Reshape if possible to have values.ndim == ndim.
    """
    if values.ndim < ndim:
        if not is_extension_array_dtype(values.dtype):
            # TODO(EA2D): https://github.com/pandas-dev/pandas/issues/23023
            # block.shape is incorrect for "2D" ExtensionArrays
            # We can't, and don't need to, reshape.
            values = np.asarray(values).reshape(1, -1)
    return values
