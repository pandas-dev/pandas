from collections import defaultdict
from functools import partial
import itertools
import operator
import re
from typing import Dict, List, Optional, Sequence, Tuple, TypeVar, Union

import numpy as np

from pandas._libs import Timedelta, Timestamp, internals as libinternals, lib
from pandas._typing import ArrayLike, DtypeObj, Label
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import (
    find_common_type,
    infer_dtype_from_scalar,
    maybe_convert_objects,
    maybe_promote,
)
from pandas.core.dtypes.common import (
    _NS_DTYPE,
    is_datetimelike_v_numeric,
    is_extension_array_dtype,
    is_list_like,
    is_numeric_v_string_like,
    is_scalar,
)
from pandas.core.dtypes.concat import concat_compat
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCExtensionArray, ABCSeries
from pandas.core.dtypes.missing import isna

import pandas.core.algorithms as algos
from pandas.core.arrays.sparse import SparseDtype
from pandas.core.base import PandasObject
from pandas.core.indexers import maybe_convert_indices
from pandas.core.indexes.api import Index, ensure_index
from pandas.core.internals.blocks import (
    Block,
    CategoricalBlock,
    DatetimeTZBlock,
    ExtensionBlock,
    ObjectValuesExtensionBlock,
    _extend_blocks,
    _merge_blocks,
    _safe_reshape,
    get_block_type,
    make_block,
)
from pandas.core.internals.concat import (  # all for concatenate_block_managers
    combine_concat_plans,
    concatenate_join_units,
    get_mgr_concatenation_plan,
    is_uniform_join_units,
)

from pandas.io.formats.printing import pprint_thing

# TODO: flexible with index=None and/or items=None

T = TypeVar("T", bound="BlockManager")


class BlockManager(PandasObject):
    """
    Core internal data structure to implement DataFrame, Series, etc.

    Manage a bunch of labeled 2D mixed-type ndarrays. Essentially it's a
    lightweight blocked set of labeled data to be manipulated by the DataFrame
    public API class

    Attributes
    ----------
    shape
    ndim
    axes
    values
    items

    Methods
    -------
    set_axis(axis, new_labels)
    copy(deep=True)

    get_dtype_counts
    get_dtypes

    apply(func, axes, block_filter_fn)

    get_bool_data
    get_numeric_data

    get_slice(slice_like, axis)
    get(label)
    iget(loc)

    take(indexer, axis)
    reindex_axis(new_labels, axis)
    reindex_indexer(new_labels, indexer, axis)

    delete(label)
    insert(loc, label, value)
    set(label, value)

    Parameters
    ----------
    blocks: Sequence of Block
    axes: Sequence of Index
    do_integrity_check: bool, default True

    Notes
    -----
    This is *not* a public API class
    """

    __slots__ = [
        "axes",
        "blocks",
        "_ndim",
        "_shape",
        "_known_consolidated",
        "_is_consolidated",
        "_blknos",
        "_blklocs",
    ]

    _blknos: np.ndarray
    _blklocs: np.ndarray

    def __init__(
        self,
        blocks: Sequence[Block],
        axes: Sequence[Index],
        do_integrity_check: bool = True,
    ):
        self.axes = [ensure_index(ax) for ax in axes]
        self.blocks: Tuple[Block, ...] = tuple(blocks)

        for block in blocks:
            if self.ndim != block.ndim:
                raise AssertionError(
                    f"Number of Block dimensions ({block.ndim}) must equal "
                    f"number of axes ({self.ndim})"
                )

        if do_integrity_check:
            self._verify_integrity()

        # Populate known_consolidate, blknos, and blklocs lazily
        self._known_consolidated = False
        self._blknos = None
        self._blklocs = None

    @classmethod
    def from_blocks(cls, blocks: List[Block], axes: List[Index]):
        """
        Constructor for BlockManager and SingleBlockManager with same signature.
        """
        return cls(blocks, axes, do_integrity_check=False)

    @property
    def blknos(self):
        """
        Suppose we want to find the array corresponding to our i'th column.

        blknos[i] identifies the block from self.blocks that contains this column.

        blklocs[i] identifies the column of interest within
        self.blocks[self.blknos[i]]
        """
        if self._blknos is None:
            # Note: these can be altered by other BlockManager methods.
            self._rebuild_blknos_and_blklocs()

        return self._blknos

    @property
    def blklocs(self):
        """
        See blknos.__doc__
        """
        if self._blklocs is None:
            # Note: these can be altered by other BlockManager methods.
            self._rebuild_blknos_and_blklocs()

        return self._blklocs

    def make_empty(self: T, axes=None) -> T:
        """ return an empty BlockManager with the items axis of len 0 """
        if axes is None:
            axes = [Index([])] + self.axes[1:]

        # preserve dtype if possible
        if self.ndim == 1:
            assert isinstance(self, SingleBlockManager)  # for mypy
            blk = self.blocks[0]
            arr = blk.values[:0]
            nb = blk.make_block_same_class(arr, placement=slice(0, 0), ndim=1)
            blocks = [nb]
        else:
            blocks = []
        return type(self).from_blocks(blocks, axes)

    def __nonzero__(self) -> bool:
        return True

    # Python3 compat
    __bool__ = __nonzero__

    @property
    def shape(self) -> Tuple[int, ...]:
        return tuple(len(ax) for ax in self.axes)

    @property
    def ndim(self) -> int:
        return len(self.axes)

    def set_axis(self, axis: int, new_labels: Index) -> None:
        # Caller is responsible for ensuring we have an Index object.
        old_len = len(self.axes[axis])
        new_len = len(new_labels)

        if new_len != old_len:
            raise ValueError(
                f"Length mismatch: Expected axis has {old_len} elements, new "
                f"values have {new_len} elements"
            )

        self.axes[axis] = new_labels

    @property
    def _is_single_block(self) -> bool:
        if self.ndim == 1:
            return True

        if len(self.blocks) != 1:
            return False

        blk = self.blocks[0]
        return blk.mgr_locs.is_slice_like and blk.mgr_locs.as_slice == slice(
            0, len(self), 1
        )

    def _rebuild_blknos_and_blklocs(self) -> None:
        """
        Update mgr._blknos / mgr._blklocs.
        """
        new_blknos = np.empty(self.shape[0], dtype=np.int64)
        new_blklocs = np.empty(self.shape[0], dtype=np.int64)
        new_blknos.fill(-1)
        new_blklocs.fill(-1)

        for blkno, blk in enumerate(self.blocks):
            rl = blk.mgr_locs
            new_blknos[rl.indexer] = blkno
            new_blklocs[rl.indexer] = np.arange(len(rl))

        if (new_blknos == -1).any():
            # TODO: can we avoid this?  it isn't cheap
            raise AssertionError("Gaps in blk ref_locs")

        self._blknos = new_blknos
        self._blklocs = new_blklocs

    @property
    def items(self) -> Index:
        return self.axes[0]

    def _get_counts(self, f):
        """ return a dict of the counts of the function in BlockManager """
        self._consolidate_inplace()
        counts = dict()
        for b in self.blocks:
            v = f(b)
            counts[v] = counts.get(v, 0) + b.shape[0]
        return counts

    def get_dtype_counts(self):
        return self._get_counts(lambda b: b.dtype.name)

    def get_dtypes(self):
        dtypes = np.array([blk.dtype for blk in self.blocks])
        return algos.take_1d(dtypes, self.blknos, allow_fill=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [self.items[b.mgr_locs.indexer] for b in self.blocks]
        axes_array = list(self.axes)

        extra_state = {
            "0.14.1": {
                "axes": axes_array,
                "blocks": [
                    dict(values=b.values, mgr_locs=b.mgr_locs.indexer)
                    for b in self.blocks
                ],
            }
        }

        # First three elements of the state are to maintain forward
        # compatibility with 0.13.1.
        return axes_array, block_values, block_items, extra_state

    def __setstate__(self, state):
        def unpickle_block(values, mgr_locs):
            return make_block(values, placement=mgr_locs)

        if isinstance(state, tuple) and len(state) >= 4 and "0.14.1" in state[3]:
            state = state[3]["0.14.1"]
            self.axes = [ensure_index(ax) for ax in state["axes"]]
            self.blocks = tuple(
                unpickle_block(b["values"], b["mgr_locs"]) for b in state["blocks"]
            )
        else:
            raise NotImplementedError("pre-0.14.1 pickles are no longer supported")

        self._post_setstate()

    def _post_setstate(self) -> None:
        self._is_consolidated = False
        self._known_consolidated = False
        self._rebuild_blknos_and_blklocs()

    def __len__(self) -> int:
        return len(self.items)

    def __repr__(self) -> str:
        output = type(self).__name__
        for i, ax in enumerate(self.axes):
            if i == 0:
                output += f"\nItems: {ax}"
            else:
                output += f"\nAxis {i}: {ax}"

        for block in self.blocks:
            output += f"\n{pprint_thing(block)}"
        return output

    def _verify_integrity(self) -> None:
        mgr_shape = self.shape
        tot_items = sum(len(x.mgr_locs) for x in self.blocks)
        for block in self.blocks:
            if block._verify_integrity and block.shape[1:] != mgr_shape[1:]:
                construction_error(tot_items, block.shape[1:], self.axes)
        if len(self.items) != tot_items:
            raise AssertionError(
                "Number of manager items must equal union of "
                f"block items\n# manager items: {len(self.items)}, # "
                f"tot_items: {tot_items}"
            )

    def reduce(self, func, *args, **kwargs):
        # If 2D, we assume that we're operating column-wise
        if self.ndim == 1:
            # we'll be returning a scalar
            blk = self.blocks[0]
            return func(blk.values, *args, **kwargs)

        res = {}
        for blk in self.blocks:
            bres = func(blk.values, *args, **kwargs)

            if np.ndim(bres) == 0:
                # EA
                assert blk.shape[0] == 1
                new_res = zip(blk.mgr_locs.as_array, [bres])
            else:
                assert bres.ndim == 1, bres.shape
                assert blk.shape[0] == len(bres), (blk.shape, bres.shape, args, kwargs)
                new_res = zip(blk.mgr_locs.as_array, bres)

            nr = dict(new_res)
            assert not any(key in res for key in nr)
            res.update(nr)

        return res

    def apply(self: T, f, filter=None, **kwargs) -> T:
        """
        Iterate over the blocks, collect and create a new BlockManager.

        Parameters
        ----------
        f : str or callable
            Name of the Block method to apply.
        filter : list, if supplied, only call the block if the filter is in
                 the block

        Returns
        -------
        BlockManager
        """
        result_blocks = []
        # fillna: Series/DataFrame is responsible for making sure value is aligned

        # filter kwarg is used in replace-* family of methods
        if filter is not None:
            filter_locs = set(self.items.get_indexer_for(filter))
            if len(filter_locs) == len(self.items):
                # All items are included, as if there were no filtering
                filter = None
            else:
                kwargs["filter"] = filter_locs

        self._consolidate_inplace()

        if f == "where":
            align_copy = True
            if kwargs.get("align", True):
                align_keys = ["other", "cond"]
            else:
                align_keys = ["cond"]
        elif f == "putmask":
            align_copy = False
            if kwargs.get("align", True):
                align_keys = ["new", "mask"]
            else:
                align_keys = ["mask"]
        else:
            align_keys = []

        # TODO(EA): may interfere with ExtensionBlock.setitem for blocks
        # with a .values attribute.
        aligned_args = {
            k: kwargs[k]
            for k in align_keys
            if not isinstance(kwargs[k], ABCExtensionArray)
            and hasattr(kwargs[k], "values")
        }

        for b in self.blocks:
            if filter is not None:
                if not b.mgr_locs.isin(filter_locs).any():
                    result_blocks.append(b)
                    continue

            if aligned_args:
                b_items = self.items[b.mgr_locs.indexer]

                for k, obj in aligned_args.items():
                    axis = obj._info_axis_number
                    kwargs[k] = obj.reindex(b_items, axis=axis, copy=align_copy)

            if callable(f):
                applied = b.apply(f, **kwargs)
            else:
                applied = getattr(b, f)(**kwargs)
            result_blocks = _extend_blocks(applied, result_blocks)

        if len(result_blocks) == 0:
            return self.make_empty(self.axes)

        return type(self).from_blocks(result_blocks, self.axes)

    def quantile(
        self,
        axis: int = 0,
        consolidate: bool = True,
        transposed: bool = False,
        interpolation="linear",
        qs=None,
        numeric_only=None,
    ) -> "BlockManager":
        """
        Iterate over blocks applying quantile reduction.
        This routine is intended for reduction type operations and
        will do inference on the generated blocks.

        Parameters
        ----------
        axis: reduction axis, default 0
        consolidate: bool, default True. Join together blocks having same
            dtype
        transposed: bool, default False
            we are holding transposed data
        interpolation : type of interpolation, default 'linear'
        qs : a scalar or list of the quantiles to be computed
        numeric_only : ignored

        Returns
        -------
        BlockManager
        """
        # Series dispatches to DataFrame for quantile, which allows us to
        #  simplify some of the code here and in the blocks
        assert self.ndim >= 2

        if consolidate:
            self._consolidate_inplace()

        def get_axe(block, qs, axes):
            # Because Series dispatches to DataFrame, we will always have
            #  block.ndim == 2
            from pandas import Float64Index

            if is_list_like(qs):
                ax = Float64Index(qs)
            else:
                ax = axes[0]
            return ax

        axes, blocks = [], []
        for b in self.blocks:
            block = b.quantile(axis=axis, qs=qs, interpolation=interpolation)

            axe = get_axe(b, qs, axes=self.axes)

            axes.append(axe)
            blocks.append(block)

        # note that some DatetimeTZ, Categorical are always ndim==1
        ndim = {b.ndim for b in blocks}
        assert 0 not in ndim, ndim

        if 2 in ndim:

            new_axes = list(self.axes)

            # multiple blocks that are reduced
            if len(blocks) > 1:
                new_axes[1] = axes[0]

                # reset the placement to the original
                for b, sb in zip(blocks, self.blocks):
                    b.mgr_locs = sb.mgr_locs

            else:
                new_axes[axis] = Index(np.concatenate([ax.values for ax in axes]))

            if transposed:
                new_axes = new_axes[::-1]
                blocks = [
                    b.make_block(b.values.T, placement=np.arange(b.shape[1]))
                    for b in blocks
                ]

            return type(self)(blocks, new_axes)

        # single block, i.e. ndim == {1}
        values = concat_compat([b.values for b in blocks])

        # compute the orderings of our original data
        if len(self.blocks) > 1:

            indexer = np.empty(len(self.axes[0]), dtype=np.intp)
            i = 0
            for b in self.blocks:
                for j in b.mgr_locs:
                    indexer[j] = i
                    i = i + 1

            values = values.take(indexer)

        return SingleBlockManager(
            make_block(values, ndim=1, placement=np.arange(len(values))),
            axes[0],
            fastpath=True,
        )

    def isna(self, func) -> "BlockManager":
        return self.apply("apply", func=func)

    def where(self, **kwargs) -> "BlockManager":
        return self.apply("where", **kwargs)

    def setitem(self, indexer, value) -> "BlockManager":
        return self.apply("setitem", indexer=indexer, value=value)

    def putmask(self, **kwargs):
        return self.apply("putmask", **kwargs)

    def diff(self, n: int, axis: int) -> "BlockManager":
        return self.apply("diff", n=n, axis=axis)

    def interpolate(self, **kwargs) -> "BlockManager":
        return self.apply("interpolate", **kwargs)

    def shift(self, **kwargs) -> "BlockManager":
        return self.apply("shift", **kwargs)

    def fillna(self, **kwargs) -> "BlockManager":
        return self.apply("fillna", **kwargs)

    def downcast(self) -> "BlockManager":
        return self.apply("downcast")

    def astype(
        self, dtype, copy: bool = False, errors: str = "raise"
    ) -> "BlockManager":
        return self.apply("astype", dtype=dtype, copy=copy, errors=errors)

    def convert(
        self,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
        coerce: bool = False,
    ) -> "BlockManager":
        return self.apply(
            "convert",
            copy=copy,
            datetime=datetime,
            numeric=numeric,
            timedelta=timedelta,
            coerce=coerce,
        )

    def replace(self, value, **kwargs) -> "BlockManager":
        assert np.ndim(value) == 0, value
        return self.apply("replace", value=value, **kwargs)

    def replace_list(
        self, src_list, dest_list, inplace: bool = False, regex: bool = False
    ) -> "BlockManager":
        """ do a list replace """
        inplace = validate_bool_kwarg(inplace, "inplace")

        # figure out our mask a-priori to avoid repeated replacements
        values = self.as_array()

        def comp(s, regex=False):
            """
            Generate a bool array by perform an equality check, or perform
            an element-wise regular expression matching
            """
            if isna(s):
                return isna(values)
            if isinstance(s, (Timedelta, Timestamp)) and getattr(s, "tz", None) is None:

                return _compare_or_regex_search(
                    maybe_convert_objects(values), s.asm8, regex
                )
            return _compare_or_regex_search(values, s, regex)

        masks = [comp(s, regex) for s in src_list]

        result_blocks = []
        src_len = len(src_list) - 1
        for blk in self.blocks:

            # its possible to get multiple result blocks here
            # replace ALWAYS will return a list
            rb = [blk if inplace else blk.copy()]
            for i, (s, d) in enumerate(zip(src_list, dest_list)):
                # TODO: assert/validate that `d` is always a scalar?
                new_rb: List[Block] = []
                for b in rb:
                    m = masks[i][b.mgr_locs.indexer]
                    convert = i == src_len
                    result = b._replace_coerce(
                        mask=m,
                        to_replace=s,
                        value=d,
                        inplace=inplace,
                        convert=convert,
                        regex=regex,
                    )
                    if m.any() or convert:
                        new_rb = _extend_blocks(result, new_rb)
                    else:
                        new_rb.append(b)
                rb = new_rb
            result_blocks.extend(rb)

        bm = type(self).from_blocks(result_blocks, self.axes)
        bm._consolidate_inplace()
        return bm

    def is_consolidated(self) -> bool:
        """
        Return True if more than one block with the same dtype
        """
        if not self._known_consolidated:
            self._consolidate_check()
        return self._is_consolidated

    def _consolidate_check(self) -> None:
        dtypes = [blk.dtype for blk in self.blocks if blk._can_consolidate]
        self._is_consolidated = len(dtypes) == len(set(dtypes))
        self._known_consolidated = True

    @property
    def is_mixed_type(self) -> bool:
        # Warning, consolidation needs to get checked upstairs
        self._consolidate_inplace()
        return len(self.blocks) > 1

    @property
    def is_numeric_mixed_type(self) -> bool:
        # Warning, consolidation needs to get checked upstairs
        self._consolidate_inplace()
        return all(block.is_numeric for block in self.blocks)

    @property
    def any_extension_types(self) -> bool:
        """Whether any of the blocks in this manager are extension blocks"""
        return any(block.is_extension for block in self.blocks)

    @property
    def is_view(self) -> bool:
        """ return a boolean if we are a single block and are a view """
        if len(self.blocks) == 1:
            return self.blocks[0].is_view

        # It is technically possible to figure out which blocks are views
        # e.g. [ b.values.base is not None for b in self.blocks ]
        # but then we have the case of possibly some blocks being a view
        # and some blocks not. setting in theory is possible on the non-view
        # blocks w/o causing a SettingWithCopy raise/warn. But this is a bit
        # complicated

        return False

    def get_bool_data(self, copy: bool = False) -> "BlockManager":
        """
        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """
        self._consolidate_inplace()
        return self.combine([b for b in self.blocks if b.is_bool], copy)

    def get_numeric_data(self, copy: bool = False) -> "BlockManager":
        """
        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """
        self._consolidate_inplace()
        return self.combine([b for b in self.blocks if b.is_numeric], copy)

    def combine(self, blocks: List[Block], copy: bool = True) -> "BlockManager":
        """ return a new manager with the blocks """
        if len(blocks) == 0:
            return self.make_empty()

        # FIXME: optimization potential
        indexer = np.sort(np.concatenate([b.mgr_locs.as_array for b in blocks]))
        inv_indexer = lib.get_reverse_indexer(indexer, self.shape[0])

        new_blocks = []
        for b in blocks:
            b = b.copy(deep=copy)
            b.mgr_locs = algos.take_1d(
                inv_indexer, b.mgr_locs.as_array, axis=0, allow_fill=False
            )
            new_blocks.append(b)

        axes = list(self.axes)
        axes[0] = self.items.take(indexer)

        return type(self).from_blocks(new_blocks, axes)

    def get_slice(self, slobj: slice, axis: int = 0) -> "BlockManager":

        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(slobj)
        elif axis == 1:
            _slicer = [slice(None)] * (axis + 1)
            _slicer[axis] = slobj
            slicer = tuple(_slicer)
            new_blocks = [blk.getitem_block(slicer) for blk in self.blocks]
        else:
            raise IndexError("Requested axis not found in manager")

        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis][slobj]

        bm = type(self)(new_blocks, new_axes, do_integrity_check=False)
        return bm

    def __contains__(self, item) -> bool:
        return item in self.items

    @property
    def nblocks(self) -> int:
        return len(self.blocks)

    def copy(self: T, deep=True) -> T:
        """
        Make deep or shallow copy of BlockManager

        Parameters
        ----------
        deep : bool or string, default True
            If False, return shallow copy (do not copy data)
            If 'all', copy data and a deep copy of the index

        Returns
        -------
        BlockManager
        """
        # this preserves the notion of view copying of axes
        if deep:
            # hit in e.g. tests.io.json.test_pandas

            def copy_func(ax):
                return ax.copy(deep=True) if deep == "all" else ax.view()

            new_axes = [copy_func(ax) for ax in self.axes]
        else:
            new_axes = list(self.axes)

        res = self.apply("copy", deep=deep)
        res.axes = new_axes
        return res

    def as_array(self, transpose: bool = False) -> np.ndarray:
        """
        Convert the blockmanager data into an numpy array.

        Parameters
        ----------
        transpose : bool, default False
            If True, transpose the return array,

        Returns
        -------
        arr : ndarray
        """
        if len(self.blocks) == 0:
            arr = np.empty(self.shape, dtype=float)
            return arr.transpose() if transpose else arr

        if self._is_single_block and self.blocks[0].is_datetimetz:
            # TODO(Block.get_values): Make DatetimeTZBlock.get_values
            # always be object dtype. Some callers seem to want the
            # DatetimeArray (previously DTI)
            arr = self.blocks[0].get_values(dtype=object)
        elif self._is_single_block or not self.is_mixed_type:
            arr = np.asarray(self.blocks[0].get_values())
        else:
            arr = self._interleave()

        return arr.transpose() if transpose else arr

    def _interleave(self) -> np.ndarray:
        """
        Return ndarray from blocks with specified item order
        Items must be contained in the blocks
        """
        dtype = _interleaved_dtype(self.blocks)

        # TODO: https://github.com/pandas-dev/pandas/issues/22791
        # Give EAs some input on what happens here. Sparse needs this.
        if isinstance(dtype, SparseDtype):
            dtype = dtype.subtype
        elif is_extension_array_dtype(dtype):
            dtype = "object"

        result = np.empty(self.shape, dtype=dtype)

        itemmask = np.zeros(self.shape[0])

        for blk in self.blocks:
            rl = blk.mgr_locs
            result[rl.indexer] = blk.get_values(dtype)
            itemmask[rl.indexer] = 1

        if not itemmask.all():
            raise AssertionError("Some items were not contained in blocks")

        return result

    def to_dict(self, copy: bool = True):
        """
        Return a dict of str(dtype) -> BlockManager

        Parameters
        ----------
        copy : bool, default True

        Returns
        -------
        values : a dict of dtype -> BlockManager

        Notes
        -----
        This consolidates based on str(dtype)
        """
        self._consolidate_inplace()

        bd: Dict[str, List[Block]] = {}
        for b in self.blocks:
            bd.setdefault(str(b.dtype), []).append(b)

        return {dtype: self.combine(blocks, copy=copy) for dtype, blocks in bd.items()}

    def fast_xs(self, loc: int):
        """
        get a cross sectional for a given location in the
        items ; handle dups

        return the result, is *could* be a view in the case of a
        single block
        """
        if len(self.blocks) == 1:
            return self.blocks[0].iget((slice(None), loc))

        items = self.items

        # non-unique (GH4726)
        if not items.is_unique:
            result = self._interleave()
            if self.ndim == 2:
                result = result.T
            return result[loc]

        # unique
        dtype = _interleaved_dtype(self.blocks)

        n = len(items)
        if is_extension_array_dtype(dtype):
            # we'll eventually construct an ExtensionArray.
            result = np.empty(n, dtype=object)
        else:
            result = np.empty(n, dtype=dtype)

        for blk in self.blocks:
            # Such assignment may incorrectly coerce NaT to None
            # result[blk.mgr_locs] = blk._slice((slice(None), loc))
            for i, rl in enumerate(blk.mgr_locs):
                result[rl] = blk.iget((i, loc))

        if isinstance(dtype, ExtensionDtype):
            result = dtype.construct_array_type()._from_sequence(result, dtype=dtype)

        return result

    def consolidate(self) -> "BlockManager":
        """
        Join together blocks having same dtype

        Returns
        -------
        y : BlockManager
        """
        if self.is_consolidated():
            return self

        bm = type(self)(self.blocks, self.axes)
        bm._is_consolidated = False
        bm._consolidate_inplace()
        return bm

    def _consolidate_inplace(self) -> None:
        if not self.is_consolidated():
            self.blocks = tuple(_consolidate(self.blocks))
            self._is_consolidated = True
            self._known_consolidated = True
            self._rebuild_blknos_and_blklocs()

    def get(self, item):
        """
        Return values for selected item (ndarray or BlockManager).
        """
        if self.items.is_unique:

            if not isna(item):
                loc = self.items.get_loc(item)
            else:
                indexer = np.arange(len(self.items))[isna(self.items)]

                # allow a single nan location indexer
                if not is_scalar(indexer):
                    if len(indexer) == 1:
                        loc = indexer.item()
                    else:
                        raise ValueError("cannot label index with a null key")

            return self.iget(loc)
        else:

            if isna(item):
                raise TypeError("cannot label index with a null key")

            indexer = self.items.get_indexer_for([item])
            return self.reindex_indexer(
                new_axis=self.items[indexer], indexer=indexer, axis=0, allow_dups=True
            )

    def iget(self, i: int) -> "SingleBlockManager":
        """
        Return the data as a SingleBlockManager.
        """
        block = self.blocks[self.blknos[i]]
        values = block.iget(self.blklocs[i])

        # shortcut for select a single-dim from a 2-dim BM
        return SingleBlockManager(
            block.make_block_same_class(
                values, placement=slice(0, len(values)), ndim=1
            ),
            self.axes[1],
            fastpath=True,
        )

    def delete(self, item):
        """
        Delete selected item (items if non-unique) in-place.
        """
        indexer = self.items.get_loc(item)

        is_deleted = np.zeros(self.shape[0], dtype=np.bool_)
        is_deleted[indexer] = True
        ref_loc_offset = -is_deleted.cumsum()

        is_blk_deleted = [False] * len(self.blocks)

        if isinstance(indexer, int):
            affected_start = indexer
        else:
            affected_start = is_deleted.nonzero()[0][0]

        for blkno, _ in _fast_count_smallints(self.blknos[affected_start:]):
            blk = self.blocks[blkno]
            bml = blk.mgr_locs
            blk_del = is_deleted[bml.indexer].nonzero()[0]

            if len(blk_del) == len(bml):
                is_blk_deleted[blkno] = True
                continue
            elif len(blk_del) != 0:
                blk.delete(blk_del)
                bml = blk.mgr_locs

            blk.mgr_locs = bml.add(ref_loc_offset[bml.indexer])

        # FIXME: use Index.delete as soon as it uses fastpath=True
        self.axes[0] = self.items[~is_deleted]
        self.blocks = tuple(
            b for blkno, b in enumerate(self.blocks) if not is_blk_deleted[blkno]
        )
        self._rebuild_blknos_and_blklocs()

    def set(self, item: Label, value):
        """
        Set new item in-place.

        Notes
        -----
        Does not consolidate.
        Adds new Block if not contained in the current items Index.
        """
        try:
            loc = self.items.get_loc(item)
        except KeyError:
            # This item wasn't present, just insert at end
            self.insert(len(self.items), item, value)
            return

        self.iset(loc, value)

    def iset(self, loc: Union[int, slice, np.ndarray], value):
        """
        Set new item in-place. Does not consolidate. Adds new Block if not
        contained in the current set of items
        """
        # FIXME: refactor, clearly separate broadcasting & zip-like assignment
        #        can prob also fix the various if tests for sparse/categorical
        if self._blklocs is None and self.ndim > 1:
            self._rebuild_blknos_and_blklocs()

        value_is_extension_type = is_extension_array_dtype(value)

        # categorical/sparse/datetimetz
        if value_is_extension_type:

            def value_getitem(placement):
                return value

        else:
            if value.ndim == self.ndim - 1:
                value = _safe_reshape(value, (1,) + value.shape)

                def value_getitem(placement):
                    return value

            else:

                def value_getitem(placement):
                    return value[placement.indexer]

            if value.shape[1:] != self.shape[1:]:
                raise AssertionError(
                    "Shape of new values must be compatible with manager shape"
                )

        if lib.is_integer(loc):
            # We have 6 tests where loc is _not_ an int.
            # In this case, get_blkno_placements will yield only one tuple,
            #  containing (self._blknos[loc], BlockPlacement(slice(0, 1, 1)))
            loc = [loc]

        # Accessing public blknos ensures the public versions are initialized
        blknos = self.blknos[loc]
        blklocs = self.blklocs[loc].copy()

        unfit_mgr_locs = []
        unfit_val_locs = []
        removed_blknos = []
        for blkno, val_locs in libinternals.get_blkno_placements(blknos, group=True):
            blk = self.blocks[blkno]
            blk_locs = blklocs[val_locs.indexer]
            if blk.should_store(value):
                blk.set(blk_locs, value_getitem(val_locs))
            else:
                unfit_mgr_locs.append(blk.mgr_locs.as_array[blk_locs])
                unfit_val_locs.append(val_locs)

                # If all block items are unfit, schedule the block for removal.
                if len(val_locs) == len(blk.mgr_locs):
                    removed_blknos.append(blkno)
                else:
                    self._blklocs[blk.mgr_locs.indexer] = -1
                    blk.delete(blk_locs)
                    self._blklocs[blk.mgr_locs.indexer] = np.arange(len(blk))

        if len(removed_blknos):
            # Remove blocks & update blknos accordingly
            is_deleted = np.zeros(self.nblocks, dtype=np.bool_)
            is_deleted[removed_blknos] = True

            new_blknos = np.empty(self.nblocks, dtype=np.int64)
            new_blknos.fill(-1)
            new_blknos[~is_deleted] = np.arange(self.nblocks - len(removed_blknos))
            self._blknos = algos.take_1d(
                new_blknos, self._blknos, axis=0, allow_fill=False
            )
            self.blocks = tuple(
                blk for i, blk in enumerate(self.blocks) if i not in set(removed_blknos)
            )

        if unfit_val_locs:
            unfit_mgr_locs = np.concatenate(unfit_mgr_locs)
            unfit_count = len(unfit_mgr_locs)

            new_blocks: List[Block] = []
            if value_is_extension_type:
                # This code (ab-)uses the fact that sparse blocks contain only
                # one item.
                new_blocks.extend(
                    make_block(
                        values=value,
                        ndim=self.ndim,
                        placement=slice(mgr_loc, mgr_loc + 1),
                    )
                    for mgr_loc in unfit_mgr_locs
                )

                self._blknos[unfit_mgr_locs] = np.arange(unfit_count) + len(self.blocks)
                self._blklocs[unfit_mgr_locs] = 0

            else:
                # unfit_val_locs contains BlockPlacement objects
                unfit_val_items = unfit_val_locs[0].append(unfit_val_locs[1:])

                new_blocks.append(
                    make_block(
                        values=value_getitem(unfit_val_items),
                        ndim=self.ndim,
                        placement=unfit_mgr_locs,
                    )
                )

                self._blknos[unfit_mgr_locs] = len(self.blocks)
                self._blklocs[unfit_mgr_locs] = np.arange(unfit_count)

            self.blocks += tuple(new_blocks)

            # Newly created block's dtype may already be present.
            self._known_consolidated = False

    def insert(self, loc: int, item, value, allow_duplicates: bool = False):
        """
        Insert item at selected position.

        Parameters
        ----------
        loc : int
        item : hashable
        value : array_like
        allow_duplicates: bool
            If False, trying to insert non-unique item will raise

        """
        if not allow_duplicates and item in self.items:
            # Should this be a different kind of error??
            raise ValueError(f"cannot insert {item}, already exists")

        if not isinstance(loc, int):
            raise TypeError("loc must be int")

        # insert to the axis; this could possibly raise a TypeError
        new_axis = self.items.insert(loc, item)

        if value.ndim == self.ndim - 1 and not is_extension_array_dtype(value):
            value = _safe_reshape(value, (1,) + value.shape)

        block = make_block(values=value, ndim=self.ndim, placement=slice(loc, loc + 1))

        for blkno, count in _fast_count_smallints(self.blknos[loc:]):
            blk = self.blocks[blkno]
            if count == len(blk.mgr_locs):
                blk.mgr_locs = blk.mgr_locs.add(1)
            else:
                new_mgr_locs = blk.mgr_locs.as_array.copy()
                new_mgr_locs[new_mgr_locs >= loc] += 1
                blk.mgr_locs = new_mgr_locs

        # Accessing public blklocs ensures the public versions are initialized
        if loc == self.blklocs.shape[0]:
            # np.append is a lot faster, let's use it if we can.
            self._blklocs = np.append(self._blklocs, 0)
            self._blknos = np.append(self._blknos, len(self.blocks))
        else:
            self._blklocs = np.insert(self._blklocs, loc, 0)
            self._blknos = np.insert(self._blknos, loc, len(self.blocks))

        self.axes[0] = new_axis
        self.blocks += (block,)

        self._known_consolidated = False

        if len(self.blocks) > 100:
            self._consolidate_inplace()

    def reindex_axis(
        self,
        new_index,
        axis: int,
        method=None,
        limit=None,
        fill_value=None,
        copy: bool = True,
    ):
        """
        Conform block manager to new index.
        """
        new_index = ensure_index(new_index)
        new_index, indexer = self.axes[axis].reindex(
            new_index, method=method, limit=limit
        )

        return self.reindex_indexer(
            new_index, indexer, axis=axis, fill_value=fill_value, copy=copy
        )

    def reindex_indexer(
        self: T,
        new_axis,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
    ) -> T:
        """
        Parameters
        ----------
        new_axis : Index
        indexer : ndarray of int64 or None
        axis : int
        fill_value : object, default None
        allow_dups : bool, default False
        copy : bool, default True


        pandas-indexer with -1's only.
        """
        if indexer is None:
            if new_axis is self.axes[axis] and not copy:
                return self

            result = self.copy(deep=copy)
            result.axes = list(self.axes)
            result.axes[axis] = new_axis
            return result

        self._consolidate_inplace()

        # some axes don't allow reindexing with dups
        if not allow_dups:
            self.axes[axis]._can_reindex(indexer)

        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(indexer, fill_tuple=(fill_value,))
        else:
            new_blocks = [
                blk.take_nd(
                    indexer,
                    axis=axis,
                    fill_tuple=(
                        fill_value if fill_value is not None else blk.fill_value,
                    ),
                )
                for blk in self.blocks
            ]

        new_axes = list(self.axes)
        new_axes[axis] = new_axis

        return type(self).from_blocks(new_blocks, new_axes)

    def _slice_take_blocks_ax0(self, slice_or_indexer, fill_tuple=None):
        """
        Slice/take blocks along axis=0.

        Overloaded for SingleBlock

        Returns
        -------
        new_blocks : list of Block
        """
        allow_fill = fill_tuple is not None

        sl_type, slobj, sllen = _preprocess_slice_or_indexer(
            slice_or_indexer, self.shape[0], allow_fill=allow_fill
        )

        if self._is_single_block:
            blk = self.blocks[0]

            if sl_type in ("slice", "mask"):
                return [blk.getitem_block(slobj, new_mgr_locs=slice(0, sllen))]
            elif not allow_fill or self.ndim == 1:
                if allow_fill and fill_tuple[0] is None:
                    _, fill_value = maybe_promote(blk.dtype)
                    fill_tuple = (fill_value,)

                return [
                    blk.take_nd(
                        slobj,
                        axis=0,
                        new_mgr_locs=slice(0, sllen),
                        fill_tuple=fill_tuple,
                    )
                ]

        if sl_type in ("slice", "mask"):
            blknos = self.blknos[slobj]
            blklocs = self.blklocs[slobj]
        else:
            blknos = algos.take_1d(
                self.blknos, slobj, fill_value=-1, allow_fill=allow_fill
            )
            blklocs = algos.take_1d(
                self.blklocs, slobj, fill_value=-1, allow_fill=allow_fill
            )

        # When filling blknos, make sure blknos is updated before appending to
        # blocks list, that way new blkno is exactly len(blocks).
        #
        # FIXME: mgr_groupby_blknos must return mgr_locs in ascending order,
        # pytables serialization will break otherwise.
        blocks = []
        for blkno, mgr_locs in libinternals.get_blkno_placements(blknos, group=True):
            if blkno == -1:
                # If we've got here, fill_tuple was not None.
                fill_value = fill_tuple[0]

                blocks.append(
                    self._make_na_block(placement=mgr_locs, fill_value=fill_value)
                )
            else:
                blk = self.blocks[blkno]

                # Otherwise, slicing along items axis is necessary.
                if not blk._can_consolidate:
                    # A non-consolidatable block, it's easy, because there's
                    # only one item and each mgr loc is a copy of that single
                    # item.
                    for mgr_loc in mgr_locs:
                        newblk = blk.copy(deep=False)
                        newblk.mgr_locs = slice(mgr_loc, mgr_loc + 1)
                        blocks.append(newblk)

                else:
                    blocks.append(
                        blk.take_nd(
                            blklocs[mgr_locs.indexer],
                            axis=0,
                            new_mgr_locs=mgr_locs,
                            fill_tuple=None,
                        )
                    )

        return blocks

    def _make_na_block(self, placement, fill_value=None):

        if fill_value is None:
            fill_value = np.nan
        block_shape = list(self.shape)
        block_shape[0] = len(placement)

        dtype, fill_value = infer_dtype_from_scalar(fill_value)
        block_values = np.empty(block_shape, dtype=dtype)
        block_values.fill(fill_value)
        return make_block(block_values, placement=placement)

    def take(self, indexer, axis: int = 1, verify: bool = True, convert: bool = True):
        """
        Take items along any axis.
        """
        self._consolidate_inplace()
        indexer = (
            np.arange(indexer.start, indexer.stop, indexer.step, dtype="int64")
            if isinstance(indexer, slice)
            else np.asanyarray(indexer, dtype="int64")
        )

        n = self.shape[axis]
        if convert:
            indexer = maybe_convert_indices(indexer, n)

        if verify:
            if ((indexer == -1) | (indexer >= n)).any():
                raise Exception("Indices must be nonzero and less than the axis length")

        new_labels = self.axes[axis].take(indexer)
        return self.reindex_indexer(
            new_axis=new_labels, indexer=indexer, axis=axis, allow_dups=True
        )

    def equals(self, other) -> bool:
        self_axes, other_axes = self.axes, other.axes
        if len(self_axes) != len(other_axes):
            return False
        if not all(ax1.equals(ax2) for ax1, ax2 in zip(self_axes, other_axes)):
            return False
        self._consolidate_inplace()
        other._consolidate_inplace()
        if len(self.blocks) != len(other.blocks):
            return False

        # canonicalize block order, using a tuple combining the mgr_locs
        # then type name because there might be unconsolidated
        # blocks (say, Categorical) which can only be distinguished by
        # the iteration order
        def canonicalize(block):
            return (block.mgr_locs.as_array.tolist(), block.dtype.name)

        self_blocks = sorted(self.blocks, key=canonicalize)
        other_blocks = sorted(other.blocks, key=canonicalize)
        return all(
            block.equals(oblock) for block, oblock in zip(self_blocks, other_blocks)
        )

    def unstack(self, unstacker_func, fill_value) -> "BlockManager":
        """
        Return a BlockManager with all blocks unstacked..

        Parameters
        ----------
        unstacker_func : callable
            A (partially-applied) ``pd.core.reshape._Unstacker`` class.
        fill_value : Any
            fill_value for newly introduced missing values.

        Returns
        -------
        unstacked : BlockManager
        """
        n_rows = self.shape[-1]
        dummy = unstacker_func(np.empty((0, 0)), value_columns=self.items)
        new_columns = dummy.get_new_columns()
        new_index = dummy.get_new_index()
        new_blocks: List[Block] = []
        columns_mask: List[np.ndarray] = []

        for blk in self.blocks:
            blocks, mask = blk._unstack(
                partial(unstacker_func, value_columns=self.items[blk.mgr_locs.indexer]),
                new_columns,
                n_rows,
                fill_value,
            )

            new_blocks.extend(blocks)
            columns_mask.extend(mask)

        new_columns = new_columns[columns_mask]

        bm = BlockManager(new_blocks, [new_columns, new_index])
        return bm


class SingleBlockManager(BlockManager):
    """ manage a single block with """

    ndim = 1
    _is_consolidated = True
    _known_consolidated = True
    __slots__ = ()

    def __init__(
        self,
        block: Block,
        axis: Union[Index, List[Index]],
        do_integrity_check: bool = False,
        fastpath: bool = False,
    ):
        assert isinstance(block, Block), type(block)

        if isinstance(axis, list):
            if len(axis) != 1:
                raise ValueError(
                    "cannot create SingleBlockManager with more than 1 axis"
                )
            axis = axis[0]

        # passed from constructor, single block, single axis
        if fastpath:
            self.axes = [axis]
        else:
            self.axes = [ensure_index(axis)]

        self.blocks = tuple([block])

    @classmethod
    def from_blocks(
        cls, blocks: List[Block], axes: List[Index]
    ) -> "SingleBlockManager":
        """
        Constructor for BlockManager and SingleBlockManager with same signature.
        """
        assert len(blocks) == 1
        assert len(axes) == 1
        return cls(blocks[0], axes[0], do_integrity_check=False, fastpath=True)

    @classmethod
    def from_array(cls, array: ArrayLike, index: Index) -> "SingleBlockManager":
        """
        Constructor for if we have an array that is not yet a Block.
        """
        block = make_block(array, placement=slice(0, len(index)), ndim=1)
        return cls(block, index, fastpath=True)

    def _post_setstate(self):
        pass

    @property
    def _block(self) -> Block:
        return self.blocks[0]

    @property
    def _blknos(self):
        """ compat with BlockManager """
        return None

    @property
    def _blklocs(self):
        """ compat with BlockManager """
        return None

    def get_slice(self, slobj: slice, axis: int = 0) -> "SingleBlockManager":
        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        blk = self._block
        array = blk._slice(slobj)
        block = blk.make_block_same_class(array, placement=range(len(array)))
        return type(self)(block, self.index[slobj], fastpath=True)

    @property
    def index(self) -> Index:
        return self.axes[0]

    @property
    def dtype(self) -> DtypeObj:
        return self._block.dtype

    def get_dtype_counts(self):
        return {self.dtype.name: 1}

    def get_dtypes(self) -> np.ndarray:
        return np.array([self._block.dtype])

    def external_values(self):
        """The array that Series.values returns"""
        return self._block.external_values()

    def internal_values(self):
        """The array that Series._values returns"""
        return self._block.internal_values()

    @property
    def _can_hold_na(self) -> bool:
        return self._block._can_hold_na

    def is_consolidated(self) -> bool:
        return True

    def _consolidate_check(self):
        pass

    def _consolidate_inplace(self):
        pass

    def delete(self, item):
        """
        Delete single item from SingleBlockManager.

        Ensures that self.blocks doesn't become empty.
        """
        loc = self.items.get_loc(item)
        self._block.delete(loc)
        self.axes[0] = self.axes[0].delete(loc)

    def fast_xs(self, loc):
        """
        fast path for getting a cross-section
        return a view of the data
        """
        raise NotImplementedError("Use series._values[loc] instead")

    def concat(self, to_concat, new_axis: Index) -> "SingleBlockManager":
        """
        Concatenate a list of SingleBlockManagers into a single
        SingleBlockManager.

        Used for pd.concat of Series objects with axis=0.

        Parameters
        ----------
        to_concat : list of SingleBlockManagers
        new_axis : Index of the result

        Returns
        -------
        SingleBlockManager
        """
        non_empties = [x for x in to_concat if len(x) > 0]

        # check if all series are of the same block type:
        if len(non_empties) > 0:
            blocks = [obj.blocks[0] for obj in non_empties]
            if len({b.dtype for b in blocks}) == 1:
                new_block = blocks[0].concat_same_type(blocks)
            else:
                values = [x.values for x in blocks]
                values = concat_compat(values)
                new_block = make_block(values, placement=slice(0, len(values), 1))
        else:
            values = [x._block.values for x in to_concat]
            values = concat_compat(values)
            new_block = make_block(values, placement=slice(0, len(values), 1))

        mgr = SingleBlockManager(new_block, new_axis)
        return mgr


# --------------------------------------------------------------------
# Constructor Helpers


def create_block_manager_from_blocks(blocks, axes):
    try:
        if len(blocks) == 1 and not isinstance(blocks[0], Block):
            # if blocks[0] is of length 0, return empty blocks
            if not len(blocks[0]):
                blocks = []
            else:
                # It's OK if a single block is passed as values, its placement
                # is basically "all items", but if there're many, don't bother
                # converting, it's an error anyway.
                blocks = [
                    make_block(values=blocks[0], placement=slice(0, len(axes[0])))
                ]

        mgr = BlockManager(blocks, axes)
        mgr._consolidate_inplace()
        return mgr

    except ValueError as e:
        blocks = [getattr(b, "values", b) for b in blocks]
        tot_items = sum(b.shape[0] for b in blocks)
        construction_error(tot_items, blocks[0].shape[1:], axes, e)


def create_block_manager_from_arrays(arrays, names, axes):

    try:
        blocks = form_blocks(arrays, names, axes)
        mgr = BlockManager(blocks, axes)
        mgr._consolidate_inplace()
        return mgr
    except ValueError as e:
        construction_error(len(arrays), arrays[0].shape, axes, e)


def construction_error(tot_items, block_shape, axes, e=None):
    """ raise a helpful message about our construction """
    passed = tuple(map(int, [tot_items] + list(block_shape)))
    # Correcting the user facing error message during dataframe construction
    if len(passed) <= 2:
        passed = passed[::-1]

    implied = tuple(len(ax) for ax in axes)
    # Correcting the user facing error message during dataframe construction
    if len(implied) <= 2:
        implied = implied[::-1]

    if passed == implied and e is not None:
        raise e
    if block_shape[0] == 0:
        raise ValueError("Empty data passed with indices specified.")
    raise ValueError(f"Shape of passed values is {passed}, indices imply {implied}")


# -----------------------------------------------------------------------


def form_blocks(arrays, names, axes):
    # put "leftover" items in float bucket, where else?
    # generalize?
    items_dict = defaultdict(list)
    extra_locs = []

    names_idx = ensure_index(names)
    if names_idx.equals(axes[0]):
        names_indexer = np.arange(len(names_idx))
    else:
        assert names_idx.intersection(axes[0]).is_unique
        names_indexer = names_idx.get_indexer_for(axes[0])

    for i, name_idx in enumerate(names_indexer):
        if name_idx == -1:
            extra_locs.append(i)
            continue

        k = names[name_idx]
        v = arrays[name_idx]

        block_type = get_block_type(v)
        items_dict[block_type.__name__].append((i, k, v))

    blocks = []
    if len(items_dict["FloatBlock"]):
        float_blocks = _multi_blockify(items_dict["FloatBlock"])
        blocks.extend(float_blocks)

    if len(items_dict["ComplexBlock"]):
        complex_blocks = _multi_blockify(items_dict["ComplexBlock"])
        blocks.extend(complex_blocks)

    if len(items_dict["TimeDeltaBlock"]):
        timedelta_blocks = _multi_blockify(items_dict["TimeDeltaBlock"])
        blocks.extend(timedelta_blocks)

    if len(items_dict["IntBlock"]):
        int_blocks = _multi_blockify(items_dict["IntBlock"])
        blocks.extend(int_blocks)

    if len(items_dict["DatetimeBlock"]):
        datetime_blocks = _simple_blockify(items_dict["DatetimeBlock"], _NS_DTYPE)
        blocks.extend(datetime_blocks)

    if len(items_dict["DatetimeTZBlock"]):
        dttz_blocks = [
            make_block(array, klass=DatetimeTZBlock, placement=i)
            for i, _, array in items_dict["DatetimeTZBlock"]
        ]
        blocks.extend(dttz_blocks)

    if len(items_dict["BoolBlock"]):
        bool_blocks = _simple_blockify(items_dict["BoolBlock"], np.bool_)
        blocks.extend(bool_blocks)

    if len(items_dict["ObjectBlock"]) > 0:
        object_blocks = _simple_blockify(items_dict["ObjectBlock"], np.object_)
        blocks.extend(object_blocks)

    if len(items_dict["CategoricalBlock"]) > 0:
        cat_blocks = [
            make_block(array, klass=CategoricalBlock, placement=i)
            for i, _, array in items_dict["CategoricalBlock"]
        ]
        blocks.extend(cat_blocks)

    if len(items_dict["ExtensionBlock"]):

        external_blocks = [
            make_block(array, klass=ExtensionBlock, placement=i)
            for i, _, array in items_dict["ExtensionBlock"]
        ]

        blocks.extend(external_blocks)

    if len(items_dict["ObjectValuesExtensionBlock"]):
        external_blocks = [
            make_block(array, klass=ObjectValuesExtensionBlock, placement=i)
            for i, _, array in items_dict["ObjectValuesExtensionBlock"]
        ]

        blocks.extend(external_blocks)

    if len(extra_locs):
        shape = (len(extra_locs),) + tuple(len(x) for x in axes[1:])

        # empty items -> dtype object
        block_values = np.empty(shape, dtype=object)
        block_values.fill(np.nan)

        na_block = make_block(block_values, placement=extra_locs)
        blocks.append(na_block)

    return blocks


def _simple_blockify(tuples, dtype):
    """
    return a single array of a block that has a single dtype; if dtype is
    not None, coerce to this dtype
    """
    values, placement = _stack_arrays(tuples, dtype)

    # TODO: CHECK DTYPE?
    if dtype is not None and values.dtype != dtype:  # pragma: no cover
        values = values.astype(dtype)

    block = make_block(values, placement=placement)
    return [block]


def _multi_blockify(tuples, dtype=None):
    """ return an array of blocks that potentially have different dtypes """
    # group by dtype
    grouper = itertools.groupby(tuples, lambda x: x[2].dtype)

    new_blocks = []
    for dtype, tup_block in grouper:

        values, placement = _stack_arrays(list(tup_block), dtype)

        block = make_block(values, placement=placement)
        new_blocks.append(block)

    return new_blocks


def _stack_arrays(tuples, dtype):

    # fml
    def _asarray_compat(x):
        if isinstance(x, ABCSeries):
            return x._values
        else:
            return np.asarray(x)

    def _shape_compat(x):
        if isinstance(x, ABCSeries):
            return (len(x),)
        else:
            return x.shape

    placement, names, arrays = zip(*tuples)

    first = arrays[0]
    shape = (len(arrays),) + _shape_compat(first)

    stacked = np.empty(shape, dtype=dtype)
    for i, arr in enumerate(arrays):
        stacked[i] = _asarray_compat(arr)

    return stacked, placement


def _interleaved_dtype(blocks: Sequence[Block]) -> Optional[DtypeObj]:
    """
    Find the common dtype for `blocks`.

    Parameters
    ----------
    blocks : List[Block]

    Returns
    -------
    dtype : np.dtype, ExtensionDtype, or None
        None is returned when `blocks` is empty.
    """
    if not len(blocks):
        return None

    return find_common_type([b.dtype for b in blocks])


def _consolidate(blocks):
    """
    Merge blocks having same dtype, exclude non-consolidating blocks
    """
    # sort by _can_consolidate, dtype
    gkey = lambda x: x._consolidate_key
    grouper = itertools.groupby(sorted(blocks, key=gkey), gkey)

    new_blocks = []
    for (_can_consolidate, dtype), group_blocks in grouper:
        merged_blocks = _merge_blocks(
            list(group_blocks), dtype=dtype, _can_consolidate=_can_consolidate
        )
        new_blocks = _extend_blocks(merged_blocks, new_blocks)
    return new_blocks


def _compare_or_regex_search(a, b, regex=False):
    """
    Compare two array_like inputs of the same shape or two scalar values

    Calls operator.eq or re.search, depending on regex argument. If regex is
    True, perform an element-wise regex matching.

    Parameters
    ----------
    a : array_like or scalar
    b : array_like or scalar
    regex : bool, default False

    Returns
    -------
    mask : array_like of bool
    """
    if not regex:
        op = lambda x: operator.eq(x, b)
    else:
        op = np.vectorize(
            lambda x: bool(re.search(b, x)) if isinstance(x, str) else False
        )

    is_a_array = isinstance(a, np.ndarray)
    is_b_array = isinstance(b, np.ndarray)

    if is_datetimelike_v_numeric(a, b) or is_numeric_v_string_like(a, b):
        # GH#29553 avoid deprecation warnings from numpy
        result = False
    else:
        result = op(a)

    if is_scalar(result) and (is_a_array or is_b_array):
        type_names = [type(a).__name__, type(b).__name__]

        if is_a_array:
            type_names[0] = f"ndarray(dtype={a.dtype})"

        if is_b_array:
            type_names[1] = f"ndarray(dtype={b.dtype})"

        raise TypeError(
            f"Cannot compare types {repr(type_names[0])} and {repr(type_names[1])}"
        )
    return result


def _fast_count_smallints(arr):
    """Faster version of set(arr) for sequences of small numbers."""
    counts = np.bincount(arr.astype(np.int_))
    nz = counts.nonzero()[0]
    return np.c_[nz, counts[nz]]


def _preprocess_slice_or_indexer(slice_or_indexer, length, allow_fill):
    if isinstance(slice_or_indexer, slice):
        return (
            "slice",
            slice_or_indexer,
            libinternals.slice_len(slice_or_indexer, length),
        )
    elif (
        isinstance(slice_or_indexer, np.ndarray) and slice_or_indexer.dtype == np.bool_
    ):
        return "mask", slice_or_indexer, slice_or_indexer.sum()
    else:
        indexer = np.asanyarray(slice_or_indexer, dtype=np.int64)
        if not allow_fill:
            indexer = maybe_convert_indices(indexer, length)
        return "fancy", indexer, len(indexer)


def concatenate_block_managers(mgrs_indexers, axes, concat_axis, copy):
    """
    Concatenate block managers into one.

    Parameters
    ----------
    mgrs_indexers : list of (BlockManager, {axis: indexer,...}) tuples
    axes : list of Index
    concat_axis : int
    copy : bool

    """
    concat_plans = [
        get_mgr_concatenation_plan(mgr, indexers) for mgr, indexers in mgrs_indexers
    ]
    concat_plan = combine_concat_plans(concat_plans, concat_axis)
    blocks = []

    for placement, join_units in concat_plan:

        if len(join_units) == 1 and not join_units[0].indexers:
            b = join_units[0].block
            values = b.values
            if copy:
                values = values.copy()
            else:
                values = values.view()
            b = b.make_block_same_class(values, placement=placement)
        elif is_uniform_join_units(join_units):
            b = join_units[0].block.concat_same_type(
                [ju.block for ju in join_units], placement=placement
            )
        else:
            b = make_block(
                concatenate_join_units(join_units, concat_axis, copy=copy),
                placement=placement,
            )
        blocks.append(b)

    return BlockManager(blocks, axes)
