import copy
import itertools
import re
import operator
from datetime import datetime, timedelta
from collections import defaultdict

import numpy as np
from pandas.core.base import PandasObject

from pandas.hashtable import Factorizer
from pandas.core.common import (_possibly_downcast_to_dtype, isnull, notnull,
                                _NS_DTYPE, _TD_DTYPE, ABCSeries, is_list_like,
                                ABCSparseSeries, _infer_dtype_from_scalar,
                                _is_null_datelike_scalar,
                                is_timedelta64_dtype, is_datetime64_dtype,)
from pandas.core.index import Index, MultiIndex, _ensure_index, _all_indexes_same
from pandas.core.indexing import (_maybe_convert_indices, _length_of_indexer)
import pandas.core.common as com
from pandas.sparse.array import _maybe_to_sparse, SparseArray
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.computation.expressions as expressions
from pandas.util.decorators import cache_readonly

from pandas.tslib import Timestamp
from pandas import compat
from pandas.compat import (range, lrange, lmap, callable, map, zip, u,
                           OrderedDict)
from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type


class Block(PandasObject):

    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas
    data structure

    Index-ignorant; let the container take care of that
    """
    __slots__ = ['_ref_locs', 'values', 'ndim']
    is_numeric = False
    is_float = False
    is_integer = False
    is_complex = False
    is_datetime = False
    is_timedelta = False
    is_bool = False
    is_object = False
    is_sparse = False
    _can_hold_na = False
    _downcast_dtype = None
    _can_consolidate = True
    _verify_integrity = True
    _ftype = 'dense'

    def __init__(self, values, placement, ndim=None, fastpath=False):

        if ndim is None:
            ndim = values.ndim

        if values.ndim != ndim:
            raise ValueError('Wrong number of dimensions')

        if len(placement) != len(values):
            raise ValueError('Wrong number of items passed %d, placement implies '
                             '%d' % (len(values), len(placement)))

        self._ref_locs = np.array(placement, dtype=np.int_, copy=True)
        self.values = values
        self.ndim = ndim

    @property
    def _consolidate_key(self):
        return (self._can_consolidate, self.dtype.name)

    @property
    def _is_single_block(self):
        return self.ndim == 1

    @property
    def is_datelike(self):
        """ return True if I am a non-datelike """
        return self.is_datetime or self.is_timedelta

    @property
    def fill_value(self):
        return np.nan

    @property
    def ref_locs(self):
        return self._ref_locs

    def __unicode__(self):

        # don't want to print out all of the items here
        name = com.pprint_thing(self.__class__.__name__)
        if self._is_single_block:

            result = '%s: %s dtype: %s' % (
                name, len(self), self.dtype)

        else:

            shape = ' x '.join([com.pprint_thing(s) for s in self.shape])
            result = '%s: %s, %s, dtype: %s' % (
                name, com.pprint_thing(self.ref_locs), shape, self.dtype)

        return result

    def __len__(self):
        return len(self.values)

    def __getstate__(self):
        return self.ref_locs, self.values

    def __setstate__(self, state):
        self._ref_locs, self.values = state
        self.ndim = self.values.ndim

    def _slice(self, slicer):
        """ return a slice of my values """
        return self.values[slicer]

    def _getitem_block(self, slicer):
        """
        Perform __getitem__-like, return result as block.
        """
        if isinstance(slicer, tuple):
            axis0_slicer = slicer[0]
        else:
            axis0_slicer = slicer

        return self.__class__(values=self.values[slicer],
                              ndim=self.ndim,
                              fastpath=True,
                              placement=self.ref_locs[axis0_slicer])

    @property
    def shape(self):
        return self.values.shape

    @property
    def itemsize(self):
        return self.values.itemsize

    @property
    def dtype(self):
        return self.values.dtype

    @property
    def ftype(self):
        return "%s:%s" % (self.dtype, self._ftype)

    def as_block(self, result):
        """ if we are not a block, then wrap as a block, must have compatible shape """
        if not isinstance(result, Block):
            result = make_block(values=result, placement=self.ref_locs,)
        return result

    def merge(self, other):
        return _merge_blocks([self, other])

    def reindex_axis(self, indexer, method=None, axis=1, fill_value=None,
                     limit=None, mask_info=None):
        """
        Reindex using pre-computed indexer information
        """
        if axis < 1:
            raise AssertionError('axis must be at least 1, got %d' % axis)
        if fill_value is None:
            fill_value = self.fill_value

        new_values = com.take_nd(self.values, indexer, axis,
                                 fill_value=fill_value, mask_info=mask_info)
        return make_block(new_values,
                          ndim=self.ndim, fastpath=True,
                          placement=self.ref_locs)

    def reindex_items_from(self, indexer, method=None,
                           fill_value=None, limit=None, copy=True):
        """
        Reindex to only those items contained in the input set of items

        E.g. if you have ['a', 'b'], and the input items is ['b', 'c', 'd'],
        then the resulting items will be ['b']

        Returns
        -------
        reindexed : Block
        """
        if fill_value is None:
            fill_value = self.fill_value

        # single block only
        assert self.ndim == 1
        new_values = com.take_1d(self.values, indexer,
                                 fill_value=fill_value)
        block = make_block(new_values,
                           ndim=self.ndim, fastpath=True,
                           placement=np.arange(len(new_values)))
        return block

    def get(self, item):
        loc = self.items.get_loc(item)
        return self.values[loc]

    def iget(self, i):
        return self.values[i]

    def set(self, locs, values, check=False):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        self.values[locs] = values

    def delete(self, loc):
        """
        Returns
        -------
        y : Block (new object)
        """
        new_values = np.delete(self.values, loc, 0)
        return make_block(new_values,
                          ndim=self.ndim, klass=self.__class__, fastpath=True,
                          placement=np.delete(self.ref_locs, loc))

    def split_block_at(self, item):
        """
        Split block into zero or more blocks around columns with given label,
        for "deleting" a column without having to copy data by returning views
        on the original array.

        Returns
        -------
        generator of Block
        """
        loc = self.items.get_loc(item)

        if type(loc) == slice or type(loc) == int:
            mask = [True] * len(self)
            mask[loc] = False
        else:  # already a mask, inverted
            mask = -loc

        for s, e in com.split_ranges(mask):
            # FIXME: drop this function
            yield make_block(self.values[s:e],
                             ndim=self.ndim,
                             klass=self.__class__,
                             fastpath=True)

    def apply(self, func, **kwargs):
        """ apply the function to my values; return a block if we are not one """
        return self.as_block(func(self.values))

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        if not self._can_hold_na:
            if inplace:
                return [self]
            else:
                return [self.copy()]

        mask = isnull(self.values)
        if limit is not None:
            if self.ndim > 2:
                raise NotImplementedError
            mask[mask.cumsum(self.ndim-1)>limit]=False

        value = self._try_fill(value)
        blocks = self.putmask(mask, value, inplace=inplace)
        return self._maybe_downcast(blocks, downcast)

    def _maybe_downcast(self, blocks, downcast=None):

        # no need to downcast our float
        # unless indicated
        if downcast is None and self.is_float:
            return blocks
        elif downcast is None and (self.is_timedelta or self.is_datetime):
            return blocks

        result_blocks = []
        for b in blocks:
            result_blocks.extend(b.downcast(downcast))

        return result_blocks

    def downcast(self, dtypes=None):
        """ try to downcast each item to the dict of dtypes if present """

        # turn it off completely
        if dtypes is False:
            return [self]

        values = self.values

        # single block handling
        if self._is_single_block:

            # try to cast all non-floats here
            if dtypes is None:
                dtypes = 'infer'

            nv = _possibly_downcast_to_dtype(values, dtypes)
            return [make_block(nv, ndim=self.ndim,
                               fastpath=True, placement=self.ref_locs)]

        # ndim > 1
        if dtypes is None:
            return [self]

        if not (dtypes == 'infer' or isinstance(dtypes, dict)):
            raise ValueError("downcast must have a dictionary or 'infer' as "
                             "its argument")

        # item-by-item
        # this is expensive as it splits the blocks items-by-item
        blocks = []
        for i, rl in enumerate(self.ref_locs):

            if dtypes == 'infer':
                dtype = 'infer'
            else:
                raise AssertionError("dtypes as dict is not supported yet")
                dtype = dtypes.get(item, self._downcast_dtype)

            if dtype is None:
                nv = _block_shape(values[i], ndim=self.ndim)
            else:
                nv = _possibly_downcast_to_dtype(values[i], dtype)
                nv = _block_shape(nv, ndim=self.ndim)

            blocks.append(make_block(nv,
                                     ndim=self.ndim, fastpath=True,
                                     placement=[rl]))

        return blocks

    def astype(self, dtype, copy=False, raise_on_error=True, values=None):
        return self._astype(dtype, copy=copy, raise_on_error=raise_on_error,
                            values=values)

    def _astype(self, dtype, copy=False, raise_on_error=True, values=None,
                klass=None):
        """
        Coerce to the new type (if copy=True, return a new copy)
        raise on an except if raise == True
        """
        dtype = np.dtype(dtype)
        if self.dtype == dtype:
            if copy:
                return self.copy()
            return self

        try:
            # force the copy here
            if values is None:
                # _astype_nansafe works fine with 1-d only
                values = com._astype_nansafe(self.values.ravel(), dtype, copy=True)
                values = values.reshape(self.values.shape)
            newb = make_block(values,
                              ndim=self.ndim, placement=self.ref_locs,
                              fastpath=True, dtype=dtype, klass=klass)
        except:
            if raise_on_error is True:
                raise
            newb = self.copy() if copy else self

        if newb.is_numeric and self.is_numeric:
            if newb.shape != self.shape:
                raise TypeError("cannot set astype for copy = [%s] for dtype "
                                "(%s [%s]) with smaller itemsize that current "
                                "(%s [%s])" % (copy, self.dtype.name,
                                               self.itemsize, newb.dtype.name,
                                               newb.itemsize))
        return newb

    def convert(self, copy=True, **kwargs):
        """ attempt to coerce any object types to better types
            return a copy of the block (if copy = True)
            by definition we are not an ObjectBlock here!  """

        return [self.copy()] if copy else [self]

    def prepare_for_merge(self, **kwargs):
        """ a regular block is ok to merge as is """
        return self

    def post_merge(self, items, **kwargs):
        """ we are non-sparse block, try to convert to a sparse block(s) """
        sparsified_mask = self.items.isin(items.keys())

        if not sparsified_mask.any():
            return self

        new_blocks = []
        for i in sparsified_mask.nonzero()[0]:
            item = self.items[i]
            ref_loc = self.ref_locs[i]

            dtypes = set(items[item])
            # this is a safe bet with multiple dtypes
            dtype = list(dtypes)[0] if len(dtypes) == 1 else np.float64

            new_blocks.append(make_block(
                values=SparseArray(self.iget(i), dtype=dtype),
                placement=[ref_loc]))

        nonsparsified_locs = (~sparsified_mask).nonzero()[0]
        if len(nonsparsified_locs):
            new_blocks.append(make_block(
                values=self.values[nonsparsified_locs],
                placement=self.ref_locs[nonsparsified_locs]))

        return new_blocks

    def _can_hold_element(self, value):
        raise NotImplementedError()

    def _try_cast(self, value):
        raise NotImplementedError()

    def _try_cast_result(self, result, dtype=None):
        """ try to cast the result to our original type,
        we may have roundtripped thru object in the mean-time """
        if dtype is None:
            dtype = self.dtype

        if self.is_integer or self.is_bool or self.is_datetime:
            pass
        elif self.is_float and result.dtype == self.dtype:

            # protect against a bool/object showing up here
            if isinstance(dtype, compat.string_types) and dtype == 'infer':
                return result
            if not isinstance(dtype, type):
                dtype = dtype.type
            if issubclass(dtype, (np.bool_, np.object_)):
                if issubclass(dtype, np.bool_):
                    if isnull(result).all():
                        return result.astype(np.bool_)
                    else:
                        result = result.astype(np.object_)
                        result[result == 1] = True
                        result[result == 0] = False
                        return result
                else:
                    return result.astype(np.object_)

            return result

        # may need to change the dtype here
        return _possibly_downcast_to_dtype(result, dtype)

    def _try_operate(self, values):
        """ return a version to operate on as the input """
        return values

    def _try_coerce_args(self, values, other):
        """ provide coercion to our input arguments """
        return values, other

    def _try_coerce_result(self, result):
        """ reverse of try_coerce_args """
        return result

    def _try_fill(self, value):
        return value

    def to_native_types(self, slicer=None, na_rep='', **kwargs):
        """ convert to our native types format, slicing if desired """

        values = self.values
        if slicer is not None:
            values = values[:, slicer]
        values = np.array(values, dtype=object)
        mask = isnull(values)
        values[mask] = na_rep
        return values.tolist()

    # block actions ####
    def copy(self, deep=True):
        values = self.values
        if deep:
            values = values.copy()
        return make_block(values, ndim=self.ndim,
                          klass=self.__class__, fastpath=True,
                          placement=self.ref_locs)

    def replace(self, to_replace, value, inplace=False, filter=None,
                regex=False):
        """ replace the to_replace value with value, possible to create new
        blocks here this is just a call to putmask. regex is not used here.
        It is used in ObjectBlocks.  It is here for API
        compatibility."""
        mask = com.mask_missing(self.values, to_replace)
        if filter is not None:
            filtered_out = ~Index(self.ref_locs, copy=False).isin(filter)
            mask[filtered_out.nonzero()[0]] = False

        if not mask.any():
            if inplace:
                return [self]
            return [self.copy()]
        return self.putmask(mask, value, inplace=inplace)

    def setitem(self, indexer, value):
        """ set the value inplace; return a new block (of a possibly different
        dtype)

        indexer is a direct slice/positional indexer; value must be a
        compatible shape
        """

        # coerce args
        values, value = self._try_coerce_args(self.values, value)
        arr_value = np.array(value)

        # cast the values to a type that can hold nan (if necessary)
        if not self._can_hold_element(value):
            dtype, _ = com._maybe_promote(arr_value.dtype)
            values = values.astype(dtype)

        transf = (lambda x: x.T) if self.ndim == 2 else (lambda x: x)
        values = transf(values)
        l = len(values)

        # length checking
        # boolean with truth values == len of the value is ok too
        if isinstance(indexer, (np.ndarray, list)):
            if is_list_like(value) and len(indexer) != len(value):
                if not (isinstance(indexer, np.ndarray) and
                        indexer.dtype == np.bool_ and
                        len(indexer[indexer]) == len(value)):
                    raise ValueError("cannot set using a list-like indexer "
                                     "with a different length than the value")

        # slice
        elif isinstance(indexer, slice):

            if is_list_like(value) and l:
                if len(value) != _length_of_indexer(indexer, values):
                    raise ValueError("cannot set using a slice indexer with a "
                                     "different length than the value")

        try:
            # setting a single element for each dim and with a rhs that could be say a list
            # GH 6043
            if arr_value.ndim == 1 and (
                np.isscalar(indexer) or (isinstance(indexer, tuple) and all([ np.isscalar(idx) for idx in indexer ]))):
                values[indexer] = value

            # if we are an exact match (ex-broadcasting),
            # then use the resultant dtype
            elif len(arr_value.shape) and arr_value.shape[0] == values.shape[0] and np.prod(arr_value.shape) == np.prod(values.shape):
                values[indexer] = value
                values = values.astype(arr_value.dtype)

            # set
            else:
                values[indexer] = value

            # coerce and try to infer the dtypes of the result
            if np.isscalar(value):
                dtype, _ = _infer_dtype_from_scalar(value)
            else:
                dtype = 'infer'
            values = self._try_coerce_result(values)
            values = self._try_cast_result(values, dtype)
            return [make_block(transf(values),
                               ndim=self.ndim, placement=self._ref_locs,
                               fastpath=True)]
        except (ValueError, TypeError) as detail:
            raise
        except Exception as detail:
            pass

        return [self]

    def putmask(self, mask, new, align=True, inplace=False):
        """ putmask the data to the block; it is possible that we may create a
        new dtype of block

        return the resulting block(s)

        Parameters
        ----------
        mask  : the condition to respect
        new : a ndarray/object
        align : boolean, perform alignment on other/cond, default is True
        inplace : perform inplace modification, default is False

        Returns
        -------
        a new block(s), the result of the putmask
        """

        new_values = self.values if inplace else self.values.copy()

        # may need to align the new
        if hasattr(new, 'reindex_axis'):
            new = new.values.T

        # may need to align the mask
        if hasattr(mask, 'reindex_axis'):
            mask = mask.values.T

        # if we are passed a scalar None, convert it here
        if not is_list_like(new) and isnull(new):
            new = self.fill_value

        if self._can_hold_element(new):
            new = self._try_cast(new)

            # pseudo-broadcast
            if isinstance(new, np.ndarray) and new.ndim == self.ndim - 1:
                new = np.repeat(new, self.shape[-1]).reshape(self.shape)

            np.putmask(new_values, mask, new)

        # maybe upcast me
        elif mask.any():

            # need to go column by column
            new_blocks = []
            if self.ndim > 1:
                for i, ref_loc in enumerate(self.ref_locs):
                    m = mask[i]
                    v = new_values[i]

                    # need a new block
                    if m.any():

                        n = new[i] if isinstance(
                            new, np.ndarray) else np.array(new)

                        # type of the new block
                        dtype, _ = com._maybe_promote(n.dtype)

                        # we need to exiplicty astype here to make a copy
                        n = n.astype(dtype)

                        nv = _putmask_smart(v, m, n)
                    else:
                        nv = v if inplace else v.copy()

                    # Put back the dimension that was taken from it and make
                    # a block out of the result.
                    block = make_block(values=nv[np.newaxis],
                                       placement=[ref_loc],
                                       fastpath=True)

                    new_blocks.append(block)

            else:
                nv = _putmask_smart(new_values, mask, new)
                new_blocks.append(make_block(values=nv,
                                             placement=self.ref_locs,
                                             fastpath=True))

            return new_blocks

        if inplace:
            return [self]

        return [make_block(new_values,
                           placement=self.ref_locs, fastpath=True)]

    def interpolate(self, method='pad', axis=0, index=None,
                    values=None, inplace=False, limit=None,
                    fill_value=None, coerce=False, downcast=None, **kwargs):

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
            m = com._clean_fill_method(method)
        except:
            m = None

        if m is not None:
            r = check_int_bool(self, inplace)
            if r is not None:
                return r
            return self._interpolate_with_fill(method=m,
                                               axis=axis,
                                               inplace=inplace,
                                               limit=limit,
                                               fill_value=fill_value,
                                               coerce=coerce,
                                               downcast=downcast)
        # try an interp method
        try:
            m = com._clean_interp_method(method, **kwargs)
        except:
            m = None

        if m is not None:
            r = check_int_bool(self, inplace)
            if r is not None:
                return r
            return self._interpolate(method=m,
                                     index=index,
                                     values=values,
                                     axis=axis,
                                     limit=limit,
                                     fill_value=fill_value,
                                     inplace=inplace,
                                     downcast=downcast,
                                     **kwargs)

        raise ValueError("invalid method '{0}' to interpolate.".format(method))

    def _interpolate_with_fill(self, method='pad', axis=0, inplace=False,
                               limit=None, fill_value=None, coerce=False,
                               downcast=None):
        """ fillna but using the interpolate machinery """

        # if we are coercing, then don't force the conversion
        # if the block can't hold the type
        if coerce:
            if not self._can_hold_na:
                if inplace:
                    return [self]
                else:
                    return [self.copy()]

        fill_value = self._try_fill(fill_value)
        values = self.values if inplace else self.values.copy()
        values = self._try_operate(values)
        values = com.interpolate_2d(values,
                                    method=method,
                                    axis=axis,
                                    limit=limit,
                                    fill_value=fill_value,
                                    dtype=self.dtype)
        values = self._try_coerce_result(values)

        blocks = [make_block(values,
                             ndim=self.ndim, klass=self.__class__,
                             fastpath=True, placement=self.ref_locs)]
        return self._maybe_downcast(blocks, downcast)

    def _interpolate(self, method=None, index=None, values=None,
                     fill_value=None, axis=0, limit=None,
                     inplace=False, downcast=None, **kwargs):
        """ interpolate using scipy wrappers """

        data = self.values if inplace else self.values.copy()

        # only deal with floats
        if not self.is_float:
            if not self.is_integer:
                return self
            data = data.astype(np.float64)

        if fill_value is None:
            fill_value = self.fill_value

        if method in ('krogh', 'piecewise_polynomial', 'pchip'):
            if not index.is_monotonic:
                raise ValueError("{0} interpolation requires that the "
                                 "index be monotonic.".format(method))
        # process 1-d slices in the axis direction

        def func(x):

            # process a 1-d slice, returning it
            # should the axis argument be handled below in apply_along_axis?
            # i.e. not an arg to com.interpolate_1d
            return com.interpolate_1d(index, x, method=method, limit=limit,
                                      fill_value=fill_value,
                                      bounds_error=False, **kwargs)

        # interp each column independently
        interp_values = np.apply_along_axis(func, axis, data)

        blocks = [make_block(interp_values,
                             ndim=self.ndim, klass=self.__class__,
                             fastpath=True, placement=self.ref_locs)]
        return self._maybe_downcast(blocks, downcast)

    def take(self, indexer, new_axis, axis=1):
        if axis < 1:
            raise AssertionError('axis must be at least 1, got %d' % axis)
        new_values = com.take_nd(self.values, indexer, axis=axis,
                                 allow_fill=False)

        # need to preserve the ref_locs and just shift them
        # GH6121
        ref_locs = None
        if not new_axis.is_unique:
            ref_locs = self._ref_locs

        return [make_block(new_values, ndim=self.ndim,
                           klass=self.__class__, placement=ref_locs, fastpath=True)]

    def get_values(self, dtype=None):
        return self.values

    def diff(self, n):
        """ return block for the diff of the values """
        new_values = com.diff(self.values, n, axis=1)
        return [make_block(values=new_values,
                           ndim=self.ndim, fastpath=True,
                           placement=self.ref_locs)]

    def shift(self, periods, axis=0):
        """ shift the block by periods, possibly upcast """
        # convert integer to float if necessary. need to do a lot more than
        # that, handle boolean etc also
        new_values, fill_value = com._maybe_upcast(self.values)
        # make sure array sent to np.roll is c_contiguous
        f_ordered = new_values.flags.f_contiguous
        if f_ordered:
            new_values = new_values.T
            axis = new_values.ndim - axis - 1
        new_values = np.roll(new_values, periods, axis=axis)
        axis_indexer = [ slice(None) ] * self.ndim
        if periods > 0:
            axis_indexer[axis] = slice(None,periods)
        else:
            axis_indexer[axis] = slice(periods,None)
        new_values[tuple(axis_indexer)] = fill_value

        # restore original order
        if f_ordered:
            new_values = new_values.T

        return [make_block(new_values,
                           ndim=self.ndim, fastpath=True,
                           placement=self.ref_locs)]

    def eval(self, func, other, raise_on_error=True, try_cast=False):
        """
        evaluate the block; return result block from the result

        Parameters
        ----------
        func  : how to combine self, other
        other : a ndarray/object
        raise_on_error : if True, raise when I can't perform the function,
            False by default (and just return the data that we had coming in)

        Returns
        -------
        a new block, the result of the func
        """
        values = self.values

        if hasattr(other, 'reindex_axis'):
            other = other.values

        # make sure that we can broadcast
        is_transposed = False
        if hasattr(other, 'ndim') and hasattr(values, 'ndim'):
            if values.ndim != other.ndim:
                    is_transposed = True
            else:
                if values.shape == other.shape[::-1]:
                    is_transposed = True
                elif values.shape[0] == other.shape[-1]:
                    is_transposed = True
                else:
                    # this is a broadcast error heree
                    raise ValueError("cannot broadcast shape [%s] with block "
                                     "values [%s]" % (values.T.shape,
                                                      other.shape))

        transf = (lambda x: x.T) if is_transposed else (lambda x: x)

        # coerce/transpose the args if needed
        values, other = self._try_coerce_args(transf(values), other)

        # get the result, may need to transpose the other
        def get_result(other):
            return self._try_coerce_result(func(values, other))

        # error handler if we have an issue operating with the function
        def handle_error():

            if raise_on_error:
                raise TypeError('Could not operate %s with block values %s'
                                % (repr(other), str(detail)))
            else:
                # return the values
                result = np.empty(values.shape, dtype='O')
                result.fill(np.nan)
                return result

        # get the result
        try:
            result = get_result(other)

        # if we have an invalid shape/broadcast error
        # GH4576, so raise instead of allowing to pass through
        except ValueError as detail:
            raise
        except Exception as detail:
            result = handle_error()

        # technically a broadcast error in numpy can 'work' by returning a
        # boolean False
        if not isinstance(result, np.ndarray):
            if not isinstance(result, np.ndarray):

                # differentiate between an invalid ndarray-ndarray comparison
                # and an invalid type comparison
                if isinstance(values, np.ndarray) and is_list_like(other):
                    raise ValueError('Invalid broadcasting comparison [%s] '
                                     'with block values' % repr(other))

                raise TypeError('Could not compare [%s] with block values'
                                % repr(other))

        # transpose if needed
        result = transf(result)

        # try to cast if requested
        if try_cast:
            result = self._try_cast_result(result)

        return [make_block(result, ndim=self.ndim,
                           fastpath=True, placement=self.ref_locs)]

    def where(self, other, cond, align=True, raise_on_error=True,
              try_cast=False):
        """
        evaluate the block; return result block(s) from the result

        Parameters
        ----------
        other : a ndarray/object
        cond  : the condition to respect
        align : boolean, perform alignment on other/cond
        raise_on_error : if True, raise when I can't perform the function,
            False by default (and just return the data that we had coming in)

        Returns
        -------
        a new block(s), the result of the func
        """

        values = self.values

        # see if we can align other
        if hasattr(other, 'reindex_axis'):
            other = other.values

        # make sure that we can broadcast
        is_transposed = False
        if hasattr(other, 'ndim') and hasattr(values, 'ndim'):
            if values.ndim != other.ndim or values.shape == other.shape[::-1]:

                # pseodo broadcast (its a 2d vs 1d say and where needs it in a
                # specific direction)
                if (other.ndim >= 1 and values.ndim - 1 == other.ndim and
                        values.shape[0] != other.shape[0]):
                    other = _block_shape(other).T
                else:
                    values = values.T
                    is_transposed = True

        # see if we can align cond
        if not hasattr(cond, 'shape'):
            raise ValueError(
                "where must have a condition that is ndarray like")

        if hasattr(cond, 'reindex_axis'):
            cond = cond.values

        # may need to undo transpose of values
        if hasattr(values, 'ndim'):
            if values.ndim != cond.ndim or values.shape == cond.shape[::-1]:
                values = values.T
                is_transposed = not is_transposed

        # our where function
        def func(c, v, o):
            if c.ravel().all():
                return v

            v, o = self._try_coerce_args(v, o)
            try:
                return self._try_coerce_result(
                    expressions.where(c, v, o, raise_on_error=True)
                )
            except Exception as detail:
                if raise_on_error:
                    raise TypeError('Could not operate [%s] with block values '
                                    '[%s]' % (repr(o), str(detail)))
                else:
                    # return the values
                    result = np.empty(v.shape, dtype='float64')
                    result.fill(np.nan)
                    return result

        # see if we can operate on the entire block, or need item-by-item
        # or if we are a single block (ndim == 1)
        result = func(cond, values, other)
        if self._can_hold_na or self.ndim == 1:

            if not isinstance(result, np.ndarray):
                raise TypeError('Could not compare [%s] with block values'
                                % repr(other))

            if is_transposed:
                result = result.T

            # try to cast if requested
            if try_cast:
                result = self._try_cast_result(result)

            return make_block(result,
                              ndim=self.ndim, placement=self.ref_locs)

        # might need to separate out blocks
        axis = cond.ndim - 1
        cond = cond.swapaxes(axis, 0)
        mask = np.array([cond[i].all() for i in range(cond.shape[0])],
                        dtype=bool)

        result_blocks = []
        for m in [mask, ~mask]:
            if m.any():
                r = self._try_cast_result(
                    result.take(m.nonzero()[0], axis=axis))
                result_blocks.append(make_block(r.T,
                                                placement=self.ref_locs[m]))

        return result_blocks

    def equals(self, other):
        if self.dtype != other.dtype or self.shape != other.shape: return False
        return np.array_equal(self.values, other.values)


class NumericBlock(Block):
    is_numeric = True
    _can_hold_na = True


class FloatOrComplexBlock(NumericBlock):

    def equals(self, other):
        if self.dtype != other.dtype or self.shape != other.shape: return False
        left, right = self.values, other.values
        return ((left == right) | (np.isnan(left) & np.isnan(right))).all()

class FloatBlock(FloatOrComplexBlock):
    is_float = True
    _downcast_dtype = 'int64'

    def _can_hold_element(self, element):
        if is_list_like(element):
            element = np.array(element)
            return issubclass(element.dtype.type, (np.floating, np.integer))
        return (isinstance(element, (float, int, np.float_, np.int_)) and
                not isinstance(bool, np.bool_))

    def _try_cast(self, element):
        try:
            return float(element)
        except:  # pragma: no cover
            return element

    def to_native_types(self, slicer=None, na_rep='', float_format=None,
                        **kwargs):
        """ convert to our native types format, slicing if desired """

        values = self.values
        if slicer is not None:
            values = values[:, slicer]
        values = np.array(values, dtype=object)
        mask = isnull(values)
        values[mask] = na_rep
        if float_format:
            imask = (-mask).ravel()
            values.flat[imask] = np.array(
                [float_format % val for val in values.ravel()[imask]])
        return values.tolist()

    def should_store(self, value):
        # when inserting a column should not coerce integers to floats
        # unnecessarily
        return (issubclass(value.dtype.type, np.floating) and
                value.dtype == self.dtype)


class ComplexBlock(FloatOrComplexBlock):
    is_complex = True

    def _can_hold_element(self, element):
        if is_list_like(element):
            element = np.array(element)
            return issubclass(element.dtype.type, (np.floating, np.integer, np.complexfloating))
        return (isinstance(element, (float, int, complex, np.float_, np.int_)) and
                not isinstance(bool, np.bool_))

    def _try_cast(self, element):
        try:
            return complex(element)
        except:  # pragma: no cover
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.complexfloating)


class IntBlock(NumericBlock):
    is_integer = True
    _can_hold_na = False

    def _can_hold_element(self, element):
        if is_list_like(element):
            element = np.array(element)
            return issubclass(element.dtype.type, np.integer)
        return com.is_integer(element)

    def _try_cast(self, element):
        try:
            return int(element)
        except:  # pragma: no cover
            return element

    def should_store(self, value):
        return com.is_integer_dtype(value) and value.dtype == self.dtype


class TimeDeltaBlock(IntBlock):
    is_timedelta = True
    _can_hold_na = True
    is_numeric = False

    @property
    def fill_value(self):
        return tslib.iNaT

    def _try_fill(self, value):
        """ if we are a NaT, return the actual fill value """
        if isinstance(value, type(tslib.NaT)) or np.array(isnull(value)).all():
            value = tslib.iNaT
        elif isinstance(value, np.timedelta64):
            pass
        elif com.is_integer(value):
            # coerce to seconds of timedelta
            value = np.timedelta64(int(value * 1e9))
        elif isinstance(value, timedelta):
            value = np.timedelta64(value)

        return value

    def _try_coerce_args(self, values, other):
        """ provide coercion to our input arguments
            we are going to compare vs i8, so coerce to floats
            repring NaT with np.nan so nans propagate
            values is always ndarray like, other may not be """
        def masker(v):
            mask = isnull(v)
            v = v.view('i8').astype('float64')
            v[mask] = np.nan
            return v

        values = masker(values)

        if _is_null_datelike_scalar(other):
            other = np.nan
        elif isinstance(other, np.timedelta64):
            other = _coerce_scalar_to_timedelta_type(other, unit='s').item()
            if other == tslib.iNaT:
                other = np.nan
        else:
            other = masker(other)

        return values, other

    def _try_operate(self, values):
        """ return a version to operate on """
        return values.view('i8')

    def _try_coerce_result(self, result):
        """ reverse of try_coerce_args / try_operate """
        if isinstance(result, np.ndarray):
            mask = isnull(result)
            if result.dtype.kind in ['i', 'f', 'O']:
                result = result.astype('m8[ns]')
            result[mask] = tslib.iNaT
        elif isinstance(result, np.integer):
            result = np.timedelta64(result)
        return result

    def should_store(self, value):
        return issubclass(value.dtype.type, np.timedelta64)

    def to_native_types(self, slicer=None, na_rep=None, **kwargs):
        """ convert to our native types format, slicing if desired """

        values = self.values
        if slicer is not None:
            values = values[:, slicer]
        mask = isnull(values)

        rvalues = np.empty(values.shape, dtype=object)
        if na_rep is None:
            na_rep = 'NaT'
        rvalues[mask] = na_rep
        imask = (-mask).ravel()
        rvalues.flat[imask] = np.array([lib.repr_timedelta64(val)
                                        for val in values.ravel()[imask]],
                                       dtype=object)
        return rvalues.tolist()


class BoolBlock(NumericBlock):
    is_bool = True
    _can_hold_na = False

    def _can_hold_element(self, element):
        if is_list_like(element):
            element = np.array(element)
            return issubclass(element.dtype.type, np.integer)
        return isinstance(element, (int, bool))

    def _try_cast(self, element):
        try:
            return bool(element)
        except:  # pragma: no cover
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.bool_)

    def replace(self, to_replace, value, inplace=False, filter=None,
                regex=False):
        to_replace_values = np.atleast_1d(to_replace)
        if not np.can_cast(to_replace_values, bool):
            to_replace = to_replace_values
        return super(BoolBlock, self).replace(to_replace, value,
                                              inplace=inplace, filter=filter,
                                              regex=regex)

class ObjectBlock(Block):
    is_object = True
    _can_hold_na = True

    def __init__(self, values, ndim=2, fastpath=False,
                 placement=None):
        if issubclass(values.dtype.type, compat.string_types):
            values = np.array(values, dtype=object)

        super(ObjectBlock, self).__init__(values, ndim=ndim,
                                          fastpath=fastpath,
                                          placement=placement)

    @property
    def is_bool(self):
        """ we can be a bool if we have only bool values but are of type
        object
        """
        return lib.is_bool_array(self.values.ravel())

    def convert(self, convert_dates=True, convert_numeric=True, convert_timedeltas=True,
                copy=True, by_item=True):
        """ attempt to coerce any object types to better types
            return a copy of the block (if copy = True)
            by definition we ARE an ObjectBlock!!!!!

            can return multiple blocks!
            """

        # attempt to create new type blocks
        blocks = []
        if by_item and not self._is_single_block:

            for i, rl in enumerate(self.ref_locs):
                values = self.iget(i)

                values = com._possibly_convert_objects(
                    values.ravel(), convert_dates=convert_dates,
                    convert_numeric=convert_numeric,
                    convert_timedeltas=convert_timedeltas,
                ).reshape(values.shape)
                values = _block_shape(values, ndim=self.ndim)
                newb = make_block(values,
                                  ndim=self.ndim, placement=[rl])
                blocks.append(newb)

        else:

            values = com._possibly_convert_objects(
                self.values.ravel(), convert_dates=convert_dates,
                convert_numeric=convert_numeric
            ).reshape(self.values.shape)
            blocks.append(make_block(values,
                                     ndim=self.ndim, placement=self.ref_locs))

        return blocks

    def set(self, locs, values, check=False):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """

        # GH6026
        if check:
            try:
                if (self.values[locs] == values).all():
                    return
            except:
                pass
        try:
            self.values[locs] = values
        except (ValueError):

            # broadcasting error
            # see GH6171
            new_shape = list(values.shape)
            new_shape[0] = len(self.items)
            self.values = np.empty(tuple(new_shape),dtype=self.dtype)
            self.values.fill(np.nan)
            self.values[locs] = values


    def _maybe_downcast(self, blocks, downcast=None):

        if downcast is not None:
            return blocks

        # split and convert the blocks
        result_blocks = []
        for blk in blocks:
            result_blocks.extend(blk.convert(convert_dates=True,
                                             convert_numeric=False))
        return result_blocks

    def _can_hold_element(self, element):
        return True

    def _try_cast(self, element):
        return element

    def should_store(self, value):
        return not issubclass(value.dtype.type,
                              (np.integer, np.floating, np.complexfloating,
                               np.datetime64, np.bool_))

    def replace(self, to_replace, value, inplace=False, filter=None,
                regex=False):
        blk = [self]
        to_rep_is_list = com.is_list_like(to_replace)
        value_is_list = com.is_list_like(value)
        both_lists = to_rep_is_list and value_is_list
        either_list = to_rep_is_list or value_is_list

        if not either_list and com.is_re(to_replace):
            blk[0], = blk[0]._replace_single(to_replace, value,
                                             inplace=inplace, filter=filter,
                                             regex=True)
        elif not (either_list or regex):
            blk = super(ObjectBlock, self).replace(to_replace, value,
                                                   inplace=inplace,
                                                   filter=filter, regex=regex)
        elif both_lists:
            for to_rep, v in zip(to_replace, value):
                blk[0], = blk[0]._replace_single(to_rep, v, inplace=inplace,
                                                 filter=filter, regex=regex)
        elif to_rep_is_list and regex:
            for to_rep in to_replace:
                blk[0], = blk[0]._replace_single(to_rep, value,
                                                 inplace=inplace,
                                                 filter=filter, regex=regex)
        else:
            blk[0], = blk[0]._replace_single(to_replace, value,
                                             inplace=inplace, filter=filter,
                                             regex=regex)
        return blk

    def _replace_single(self, to_replace, value, inplace=False, filter=None,
                        regex=False):
        # to_replace is regex compilable
        to_rep_re = regex and com.is_re_compilable(to_replace)

        # regex is regex compilable
        regex_re = com.is_re_compilable(regex)

        # only one will survive
        if to_rep_re and regex_re:
            raise AssertionError('only one of to_replace and regex can be '
                                 'regex compilable')

        # if regex was passed as something that can be a regex (rather than a
        # boolean)
        if regex_re:
            to_replace = regex

        regex = regex_re or to_rep_re

        # try to get the pattern attribute (compiled re) or it's a string
        try:
            pattern = to_replace.pattern
        except AttributeError:
            pattern = to_replace

        # if the pattern is not empty and to_replace is either a string or a
        # regex
        if regex and pattern:
            rx = re.compile(to_replace)
        else:
            # if the thing to replace is not a string or compiled regex call
            # the superclass method -> to_replace is some kind of object
            result = super(ObjectBlock, self).replace(to_replace, value,
                                                      inplace=inplace,
                                                      filter=filter,
                                                      regex=regex)
            if not isinstance(result, list):
                result = [result]
            return result

        new_values = self.values if inplace else self.values.copy()

        # deal with replacing values with objects (strings) that match but
        # whose replacement is not a string (numeric, nan, object)
        if isnull(value) or not isinstance(value, compat.string_types):
            def re_replacer(s):
                try:
                    return value if rx.search(s) is not None else s
                except TypeError:
                    return s
        else:
            # value is guaranteed to be a string here, s can be either a string
            # or null if it's null it gets returned
            def re_replacer(s):
                try:
                    return rx.sub(value, s)
                except TypeError:
                    return s

        f = np.vectorize(re_replacer, otypes=[self.dtype])

        if filter is None:
            filt = slice(None)
        else:
            filt = (Index(self.ref_locs, copy=False)
                    .isin(filter).nonzero()[0])

        new_values[filt] = f(new_values[filt])

        return [self if inplace else
                make_block(new_values,
                           fastpath=True, placement=self.ref_locs)]


class DatetimeBlock(Block):
    is_datetime = True
    _can_hold_na = True

    def __init__(self, values, placement,
                 fastpath=False, **kwargs):
        if values.dtype != _NS_DTYPE:
            values = tslib.cast_to_nanoseconds(values)

        super(DatetimeBlock, self).__init__(values,
                                            fastpath=True, placement=placement,
                                            **kwargs)

    def _can_hold_element(self, element):
        if is_list_like(element):
            element = np.array(element)
            return element.dtype == _NS_DTYPE or element.dtype == np.int64
        return (com.is_integer(element) or
                isinstance(element, datetime) or
                isnull(element))

    def _try_cast(self, element):
        try:
            return int(element)
        except:
            return element

    def _try_operate(self, values):
        """ return a version to operate on """
        return values.view('i8')

    def _try_coerce_args(self, values, other):
        """ provide coercion to our input arguments
            we are going to compare vs i8, so coerce to integer
            values is always ndarra like, other may not be """
        values = values.view('i8')
        if _is_null_datelike_scalar(other):
            other = tslib.iNaT
        elif isinstance(other, datetime):
            other = lib.Timestamp(other).asm8.view('i8')
        else:
            other = other.view('i8')

        return values, other

    def _try_coerce_result(self, result):
        """ reverse of try_coerce_args """
        if isinstance(result, np.ndarray):
            if result.dtype == 'i8':
                result = tslib.array_to_datetime(
                    result.astype(object).ravel()).reshape(result.shape)
            elif result.dtype.kind in ['i', 'f', 'O']:
                result = result.astype('M8[ns]')
        elif isinstance(result, (np.integer, np.datetime64)):
            result = lib.Timestamp(result)
        return result

    @property
    def fill_value(self):
        return tslib.iNaT

    def _try_fill(self, value):
        """ if we are a NaT, return the actual fill value """
        if isinstance(value, type(tslib.NaT)) or np.array(isnull(value)).all():
            value = tslib.iNaT
        return value

    def fillna(self, value, limit=None,
               inplace=False, downcast=None):

        # straight putmask here
        values = self.values if inplace else self.values.copy()
        mask = isnull(self.values)
        value = self._try_fill(value)
        if limit is not None:
            if self.ndim > 2:
                raise NotImplementedError
            mask[mask.cumsum(self.ndim-1)>limit]=False

        np.putmask(values, mask, value)
        return [self if inplace else
                make_block(values,
                           fastpath=True, placement=self.ref_locs)]

    def to_native_types(self, slicer=None, na_rep=None, date_format=None,
                        **kwargs):
        """ convert to our native types format, slicing if desired """

        values = self.values
        if slicer is not None:
            values = values[:, slicer]
        mask = isnull(values)

        rvalues = np.empty(values.shape, dtype=object)
        if na_rep is None:
            na_rep = 'NaT'
        rvalues[mask] = na_rep
        imask = (-mask).ravel()

        if date_format is None:
            date_formatter = lambda x: Timestamp(x)._repr_base
        else:
            date_formatter = lambda x: Timestamp(x).strftime(date_format)

        rvalues.flat[imask] = np.array([date_formatter(val) for val in
                                        values.ravel()[imask]], dtype=object)

        return rvalues.tolist()

    def should_store(self, value):
        return issubclass(value.dtype.type, np.datetime64)

    def astype(self, dtype, copy=False, raise_on_error=True):
        """
        handle convert to object as a special case
        """
        klass = None
        if np.dtype(dtype).type == np.object_:
            klass = ObjectBlock
        return self._astype(dtype, copy=copy, raise_on_error=raise_on_error,
                            klass=klass)

    def set(self, locs, values, check=False):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        if values.dtype != _NS_DTYPE:
            # Workaround for numpy 1.6 bug
            values = tslib.cast_to_nanoseconds(values)

        self.values[locs] = values

    def get_values(self, dtype=None):
        # return object dtype as Timestamps
        if dtype == object:
            return lib.map_infer(self.values.ravel(), lib.Timestamp)\
                      .reshape(self.values.shape)
        return self.values


class SparseBlock(Block):

    """ implement as a list of sparse arrays of the same dtype """
    __slots__ = ['_ref_locs', 'ndim', 'values']
    is_sparse = True
    is_numeric = True
    _can_hold_na = True
    _can_consolidate = False
    _verify_integrity = False
    _ftype = 'sparse'

    def __init__(self, values, placement,
                 ndim=None, fastpath=False,):

        # kludgetastic
        if ndim is not None:
            if ndim == 1:
                ndim = 1
            elif ndim > 2:
                ndim = ndim
        else:
            if len(placement) != 1:
                ndim = 1
            else:
                ndim = 2
        self.ndim = ndim

        self._ref_locs = np.array(placement, dtype=np.int_, copy=True)

        self.values = values

    @property
    def shape(self):
        return (len(self.ref_locs), self.sp_index.length)

    @property
    def itemsize(self):
        return self.dtype.itemsize

    @property
    def fill_value(self):
        return self.values.fill_value

    @fill_value.setter
    def fill_value(self, v):
        # we may need to upcast our fill to match our dtype
        if issubclass(self.dtype.type, np.floating):
            v = float(v)
        self.values.fill_value = v

    @property
    def sp_values(self):
        return self.values.sp_values

    @sp_values.setter
    def sp_values(self, v):
        # reset the sparse values
        self.values = SparseArray(v, sparse_index=self.sp_index,
                                  kind=self.kind, dtype=v.dtype,
                                  fill_value=self.fill_value, copy=False)

    def iget(self, col):
        if col != 0:
            raise IndexError("SparseBlock only contains one item")
        return self.values

    @property
    def sp_index(self):
        return self.values.sp_index

    @property
    def kind(self):
        return self.values.kind

    def __len__(self):
        try:
            return self.sp_index.length
        except:
            return 0

    def should_store(self, value):
        return isinstance(value, SparseArray)

    def set(self, locs, values, check=False):
        assert locs.tolist() == [0]
        self.values = values

    def get(self, item):
        if self.ndim == 1:
            loc = self.items.get_loc(item)
            return self.values[loc]
        else:
            return self.values

    def _slice(self, slicer):
        """ return a slice of my values (but densify first) """
        return self.get_values()[slicer]

    def get_values(self, dtype=None):
        """ need to to_dense myself (and always return a ndim sized object) """
        values = self.values.to_dense()
        if values.ndim == self.ndim - 1:
            values = values.reshape((1,) + values.shape)
        return values

    def copy(self, deep=True):
        return self.make_block(values=self.values,
                               sparse_index=self.sp_index,
                               kind=self.kind, copy=deep,
                               placement=self.ref_locs)

    def make_block(self, values, placement,
                   sparse_index=None, kind=None, dtype=None, fill_value=None,
                   copy=False, fastpath=True):
        """ return a new block """
        if dtype is None:
            dtype = self.dtype
        if fill_value is None:
            fill_value = self.fill_value
        new_values = SparseArray(values, sparse_index=sparse_index,
                                 kind=kind or self.kind, dtype=dtype,
                                 fill_value=fill_value, copy=copy)
        return make_block(new_values, ndim=self.ndim,
                          fastpath=fastpath, placement=placement)

    def interpolate(self, method='pad', axis=0, inplace=False,
                    limit=None, fill_value=None, **kwargs):

        values = com.interpolate_2d(
            self.values.to_dense(), method, axis, limit, fill_value)
        return self.make_block(values=values,
                               placement=self.ref_locs)

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        # we may need to upcast our fill to match our dtype
        if limit is not None:
            raise NotImplementedError
        if issubclass(self.dtype.type, np.floating):
            value = float(value)
        values = self.values if inplace else self.values.copy()
        return [self.make_block(values=values.get_values(value),
                                fill_value=value, placement=self.ref_locs)]


    def shift(self, periods, axis=0):
        """ shift the block by periods """
        N = len(self.values.T)
        indexer = np.zeros(N, dtype=int)
        if periods > 0:
            indexer[periods:] = np.arange(N - periods)
        else:
            indexer[:periods] = np.arange(-periods, N)
        new_values = self.values.to_dense().take(indexer)
        # convert integer to float if necessary. need to do a lot more than
        # that, handle boolean etc also
        new_values, fill_value = com._maybe_upcast(new_values)
        if periods > 0:
            new_values[:periods] = fill_value
        else:
            new_values[periods:] = fill_value
        return [self.make_block(new_values, placement=self.ref_locs)]

    def take(self, indexer, new_axis, axis=1):
        """ going to take our items
            along the long dimension"""
        if axis < 1:
            raise AssertionError('axis must be at least 1, got %d' % axis)

        return [self.make_block(self.values.take(indexer))]

    def reindex_axis(self, indexer, method=None, axis=1, fill_value=None,
                     limit=None, mask_info=None):
        """
        Reindex using pre-computed indexer information
        """
        if axis < 1:
            raise AssertionError('axis must be at least 1, got %d' % axis)

        # taking on the 0th axis always here
        if fill_value is None:
            fill_value = self.fill_value
        return self.make_block(self.values.take(indexer),
                               fill_value=fill_value,
                               placement=self.ref_locs)

    def reindex_items_from(self, indexer, method=None,
                           fill_value=None, limit=None, copy=True):
        """
        Reindex to only those items contained in the input set of items

        E.g. if you have ['a', 'b'], and the input items is ['b', 'c', 'd'],
        then the resulting items will be ['b']

        Returns
        -------
        reindexed : Block
        """

        # 1-d always
        if indexer is None:
            indexer = np.arange(len(self.ref_locs))

        # single block only
        assert self.ndim == 1
        new_values = com.take_1d(self.values.values, indexer)

        # fill if needed
        if method is not None or limit is not None:
            if fill_value is None:
                fill_value = self.fill_value
            new_values = com.interpolate_2d(new_values, method=method,
                                            limit=limit, fill_value=fill_value)

        return self.make_block(new_values,
                               copy=copy,
                               placement=np.arange(len(indexer)))

    def sparse_reindex(self, new_index):
        """ sparse reindex and return a new block
            current reindex only works for float64 dtype! """
        values = self.values
        values = values.sp_index.to_int_index().reindex(
            values.sp_values.astype('float64'), values.fill_value, new_index)
        return self.make_block(values, sparse_index=new_index,
                               placement=self.ref_locs)

    def split_block_at(self, item):
        if len(self.items) == 1 and item == self.items[0]:
            return []
        return super(SparseBlock, self).split_block_at(self, item)

    def _try_cast_result(self, result, dtype=None):
        return result


def make_block(values, placement, klass=None, ndim=None,
               dtype=None, fastpath=False):
    if klass is None:
        dtype = dtype or values.dtype
        vtype = dtype.type

        if isinstance(values, SparseArray):
            klass = SparseBlock
        elif issubclass(vtype, np.floating):
            klass = FloatBlock
        elif (issubclass(vtype, np.integer) and
                issubclass(vtype, np.timedelta64)):
            klass = TimeDeltaBlock
        elif (issubclass(vtype, np.integer) and
                not issubclass(vtype, np.datetime64)):
            klass = IntBlock
        elif dtype == np.bool_:
            klass = BoolBlock
        elif issubclass(vtype, np.datetime64):
            klass = DatetimeBlock
        elif issubclass(vtype, np.complexfloating):
            klass = ComplexBlock

        # try to infer a DatetimeBlock, or set to an ObjectBlock
        else:

            if np.prod(values.shape):
                flat = values.ravel()

                # try with just the first element; we just need to see if
                # this is a datetime or not
                inferred_type = lib.infer_dtype(flat[0:1])
                if inferred_type in ['datetime', 'datetime64']:

                    # we have an object array that has been inferred as
                    # datetime, so convert it
                    try:
                        values = tslib.array_to_datetime(
                            flat).reshape(values.shape)
                        if issubclass(values.dtype.type, np.datetime64):
                            klass = DatetimeBlock
                    except:  # it already object, so leave it
                        pass

            if klass is None:
                klass = ObjectBlock

    return klass(values, ndim=ndim, fastpath=fastpath,
                 placement=placement)


# TODO: flexible with index=None and/or items=None


class BlockManager(PandasObject):

    """
    Core internal data structure to implement DataFrame

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
    get_ftype_counts
    get_dtypes
    get_ftypes

    apply(func, axes, block_filter_fn)

    get_bool_data
    get_numeric_data

    get_slice(slice_like, axis)
    get(label)
    iget(loc)
    get_scalar(label_tup)

    take(indexer, axis)
    reindex_axis(new_labels, axis)
    reindex_indexer(new_labels, indexer, axis)

    delete(label)
    insert(loc, label, value)
    set(label, value)

    Parameters
    ----------


    Notes
    -----
    This is *not* a public API class
    """
    __slots__ = ['axes', 'blocks', '_ndim', '_shape', '_known_consolidated',
                 '_is_consolidated', '_has_sparse', '_ref_locs']

    def __init__(self, blocks, axes, do_integrity_check=True, fastpath=True):
        self.axes = [_ensure_index(ax) for ax in axes]
        self.blocks = blocks

        for block in blocks:
            if block.is_sparse:
                if len(block.ref_locs) != 1:
                    raise AssertionError("Sparse block refers to multiple items")
            else:
                if self.ndim != block.ndim:
                    raise AssertionError(('Number of Block dimensions (%d) must '
                                          'equal number of axes (%d)')
                                         % (block.ndim, self.ndim))

        if do_integrity_check:
            self._verify_integrity()

        self._has_sparse = False
        self._consolidate_check()

        self._rebuild_ref_locs()

    def make_empty(self, axes=None):
        """ return an empty BlockManager with the items axis of len 0 """
        if axes is None:
            axes = [_ensure_index([])] + [
                _ensure_index(a) for a in self.axes[1:]
            ]

        # preserve dtype if possible
        if self.ndim == 1:
            blocks = np.array([], dtype=self.dtype)
        else:
            blocks = []
        return self.__class__(blocks, axes)

    def __nonzero__(self):
        return True

    # Python3 compat
    __bool__ = __nonzero__

    @property
    def shape(self):
        return tuple(len(ax) for ax in self.axes)

    @property
    def ndim(self):
        return len(self.axes)

    def set_axis(self, axis, new_labels):
        new_labels = _ensure_index(new_labels)
        old_len = len(self.axes[axis])
        new_len = len(new_labels)

        if new_len != old_len:
            raise ValueError('Length mismatch: Expected axis has %d elements, '
                             'new values have %d elements' % (old_len, new_len))

        self.axes[axis] = new_labels

    def _rebuild_ref_locs(self):
        """
        Update mgr._ref_locs according to blk.ref_locs.
        """
        blocks = np.empty(self.shape[0], dtype=np.object_)
        blk_locs = np.empty(self.shape[0], dtype=np.int_)
        blk_locs.fill(-1)

        for blk in self.blocks:
            rl = blk.ref_locs
            blocks[rl] = blk
            blk_locs[rl] = np.arange(len(rl))

        if (blk_locs == -1).any():
            raise AssertionError("Gaps in blk ref_locs")

        self._ref_locs = lib.fast_zip([blocks, blk_locs])

    # make items read only for now
    def _get_items(self):
        return self.axes[0]
    items = property(fget=_get_items)

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

    def get_ftype_counts(self):
        return self._get_counts(lambda b: b.ftype)

    def get_dtypes(self):
        return [rl[0].dtype for rl in self._ref_locs]

    def get_ftypes(self):
        return [rl[0].ftype for rl in self._ref_locs]

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [self.items.take(b.ref_locs) for b in self.blocks]
        axes_array = [ax for ax in self.axes]
        return axes_array, block_values, block_items

    def __setstate__(self, state):
        # discard anything after 3rd, support beta pickling format for a little
        # while longer
        ax_arrays, bvalues, bitems = state[:3]

        self.axes = [_ensure_index(ax) for ax in ax_arrays]

        blocks = []
        for values, items in zip(bvalues, bitems):

            # numpy < 1.7 pickle compat
            if values.dtype == 'M8[us]':
                values = values.astype('M8[ns]')

            blk = make_block(values,
                             placement=self.axes[0].get_indexer(items))
            blocks.append(blk)
        self.blocks = blocks

        self._post_setstate()

    def _post_setstate(self):
        self._is_consolidated = False
        self._known_consolidated = False
        self._rebuild_ref_locs()
        self._set_has_sparse()

    def __len__(self):
        return len(self.items)

    def __unicode__(self):
        output = com.pprint_thing(self.__class__.__name__)
        for i, ax in enumerate(self.axes):
            if i == 0:
                output += u('\nItems: %s') % ax
            else:
                output += u('\nAxis %d: %s') % (i, ax)

        for block in self.blocks:
            output += u('\n%s') % com.pprint_thing(block)
        return output

    def _verify_integrity(self):
        mgr_shape = self.shape
        tot_items = sum(len(x.ref_locs) for x in self.blocks)
        for block in self.blocks:
            if not block.is_sparse and block.shape[1:] != mgr_shape[1:]:
                construction_error(tot_items, block.shape[1:], self.axes)
        if len(self.items) != tot_items:
            raise AssertionError('Number of manager items must equal union of '
                                 'block items\n# manager items: {0}, # '
                                 'tot_items: {1}'.format(len(self.items),
                                                         tot_items))

    def apply(self, f, axes=None, filter=None, do_integrity_check=False, **kwargs):
        """
        iterate over the blocks, collect and create a new block manager

        Parameters
        ----------
        f : the callable or function name to operate on at the block level
        axes : optional (if not supplied, use self.axes)
        filter : list, if supplied, only call the block if the filter is in
                 the block
        do_integrity_check : boolean, default False. Do the block manager integrity check

        Returns
        -------
        Block Manager (new object)

        """

        result_blocks = []

        if filter is not None:
            # filter kwarg is used in replace-* family of methods
            filter_locs = set(self.items.get_indexer_for(filter))
            kwargs['filter'] = filter_locs

        if f == 'where' and kwargs.get('align', True):
            align_copy = True
            align_keys = ['other', 'cond']
        elif f == 'putmask' and kwargs.get('align', True):
            align_copy = False
            align_keys = ['new', 'mask']
        elif f == 'eval':
            align_copy = False
            align_keys = ['other']
        elif f == 'fillna':
            # fillna internally does putmask, maybe it's better to do this
            # at mgr, not block level?
            align_copy = False
            align_keys = ['value']
        else:
            align_keys = []

        aligned_args = dict((k, kwargs[k]) for k in align_keys
                            if hasattr(kwargs[k], 'reindex_axis'))

        for b in self.blocks:
            if filter is not None:
                valid_locs = filter_locs.intersection(b.ref_locs)
                if not valid_locs:
                    result_blocks.append(b)
                    continue

            if aligned_args:
                b_items = self.items.take(b.ref_locs)

                for k, obj in aligned_args.items():
                    axis = getattr(obj, '_info_axis_number', 0)
                    kwargs[k] = obj.reindex_axis(b_items, axis=axis,
                                                 copy=align_copy)

            applied = getattr(b, f)(**kwargs)

            if isinstance(applied, list):
                result_blocks.extend(applied)
            else:
                result_blocks.append(applied)

        if len(result_blocks) == 0:
            return self.make_empty(axes or self.axes)
        bm = self.__class__(result_blocks, axes or self.axes,
                            do_integrity_check=do_integrity_check)
        bm._consolidate_inplace()
        return bm

    def isnull(self, **kwargs):
        return self.apply('apply', **kwargs)

    def where(self, **kwargs):
        return self.apply('where', **kwargs)

    def eval(self, **kwargs):
        return self.apply('eval', **kwargs)

    def setitem(self, **kwargs):
        return self.apply('setitem', **kwargs)

    def putmask(self, **kwargs):
        return self.apply('putmask', **kwargs)

    def diff(self, **kwargs):
        return self.apply('diff', **kwargs)

    def interpolate(self, **kwargs):
        return self.apply('interpolate', **kwargs)

    def shift(self, **kwargs):
        return self.apply('shift', **kwargs)

    def fillna(self, **kwargs):
        return self.apply('fillna', **kwargs)

    def downcast(self, **kwargs):
        return self.apply('downcast', **kwargs)

    def astype(self, dtype, **kwargs):
        return self.apply('astype', dtype=dtype, **kwargs)

    def convert(self, **kwargs):
        return self.apply('convert', **kwargs)

    def replace(self, **kwargs):
        return self.apply('replace', **kwargs)

    def replace_list(self, src_list, dest_list, inplace=False, regex=False):
        """ do a list replace """

        # figure out our mask a-priori to avoid repeated replacements
        values = self.as_matrix()

        def comp(s):
            if isnull(s):
                return isnull(values)
            return _possibly_compare(values, getattr(s, 'asm8', s),
                                     operator.eq)
        masks = [comp(s) for i, s in enumerate(src_list)]

        result_blocks = []
        for blk in self.blocks:

            # its possible to get multiple result blocks here
            # replace ALWAYS will return a list
            rb = [blk if inplace else blk.copy()]
            for i, (s, d) in enumerate(zip(src_list, dest_list)):
                new_rb = []
                for b in rb:
                    if b.dtype == np.object_:
                        result = b.replace(s, d, inplace=inplace,
                                           regex=regex)
                        if isinstance(result, list):
                            new_rb.extend(result)
                        else:
                            new_rb.append(result)
                    else:
                        # get our mask for this element, sized to this
                        # particular block
                        m = masks[i][b.ref_locs]
                        if m.any():
                            new_rb.extend(b.putmask(m, d, inplace=True))
                        else:
                            new_rb.append(b)
                rb = new_rb
            result_blocks.extend(rb)

        bm = self.__class__(result_blocks, self.axes)
        bm._consolidate_inplace()
        return bm

    def is_consolidated(self):
        """
        Return True if more than one block with the same dtype
        """
        if not self._known_consolidated:
            self._consolidate_check()
        return self._is_consolidated

    def _consolidate_check(self):
        ftypes = [blk.ftype for blk in self.blocks]
        self._is_consolidated = len(ftypes) == len(set(ftypes))
        self._known_consolidated = True
        self._set_has_sparse()

    def _set_has_sparse(self):
        self._has_sparse = any((blk.is_sparse for blk in self.blocks))

    @property
    def is_mixed_type(self):
        # Warning, consolidation needs to get checked upstairs
        self._consolidate_inplace()
        return len(self.blocks) > 1

    @property
    def is_numeric_mixed_type(self):
        # Warning, consolidation needs to get checked upstairs
        self._consolidate_inplace()
        return all([block.is_numeric for block in self.blocks])

    @property
    def is_datelike_mixed_type(self):
        # Warning, consolidation needs to get checked upstairs
        self._consolidate_inplace()
        return any([block.is_datelike for block in self.blocks])

    def get_bool_data(self, copy=False):
        """
        Parameters
        ----------
        copy : boolean, default False
            Whether to copy the blocks
        """
        self._consolidate_inplace()
        return self.combine([b for b in self.blocks if b.is_bool], copy)

    def get_numeric_data(self, copy=False):
        """
        Parameters
        ----------
        copy : boolean, default False
            Whether to copy the blocks
        """
        self._consolidate_inplace()
        return self.combine([b for b in self.blocks if b.is_numeric], copy)

    def combine(self, blocks, copy=True):
        """ return a new manager with the blocks """
        if len(blocks) == 0:
            return self.make_empty()

        indexer = np.sort(np.concatenate([b.ref_locs for b in blocks]))
        inv_indexer = _invert_reordering(indexer)
        new_items = self.items.take(indexer)

        new_blocks = []
        for b in blocks:
            b = b.copy(deep=copy)
            b._ref_locs = inv_indexer.take(b.ref_locs)
            new_blocks.append(b)

        new_axes = list(self.axes)
        new_axes[0] = new_items
        return self.__class__(new_blocks, new_axes, do_integrity_check=False)

    def get_slice(self, slobj, axis=0):
        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis][slobj]

        if axis == 0:
            new_items = new_axes[0]

            # we want to preserver the view of a single-block
            if (len(self.blocks) == 1 and
                (self.blocks[0]._ref_locs == np.arange(self.shape[0])).all()):
                blk = self.blocks[0]
                newb = make_block(blk._slice(slobj),
                                  klass=blk.__class__, fastpath=True,
                                  placement=np.arange(len(new_items)))

                new_blocks = [newb]
            else:
                return self.reindex_indexer(
                    new_items, indexer=np.arange(len(self.items))[slobj],
                    axis=0, allow_dups=True)
        else:
            slicer = [slice(None)] * self.ndim
            slicer[axis] = slobj

            new_blocks = [make_block(block._slice(slicer),
                                     klass=block.__class__,
                                     fastpath=True,
                                     placement=block.ref_locs)
                          for block in self.blocks]

        bm = self.__class__(new_blocks, new_axes, do_integrity_check=False)
        bm._consolidate_inplace()
        return bm

    def __contains__(self, item):
        return item in self.items

    @property
    def nblocks(self):
        return len(self.blocks)

    def copy(self, deep=True):
        """
        Make deep or shallow copy of BlockManager

        Parameters
        ----------
        deep : boolean, default True
            If False, return shallow copy (do not copy data)

        Returns
        -------
        copy : BlockManager
        """
        if deep:
            new_axes = [ax.view() for ax in self.axes]
        else:
            new_axes = list(self.axes)
        return self.apply('copy', axes=new_axes, deep=deep,
                          do_integrity_check=False)

    def as_matrix(self, items=None):
        if len(self.blocks) == 0:
            return np.empty(self.shape, dtype=float)

        if items is not None:
            mgr = self.reindex_axis(items, axis=0)
        else:
            mgr = self

        if (len(mgr.blocks) == 1 and
            (mgr.blocks[0]._ref_locs is None or
             (mgr.blocks[0]._ref_locs == np.arange(mgr.shape[0])).all())):
            return mgr.blocks[0].get_values()
        else:
            return mgr._interleave()

    def _interleave(self):
        """
        Return ndarray from blocks with specified item order
        Items must be contained in the blocks
        """
        dtype = _interleaved_dtype(self.blocks)

        result = np.empty(self.shape, dtype=dtype)
        itemmask = np.zeros(self.shape[0])

        for blk in self.blocks:
            rl = blk.ref_locs
            result[rl] = blk.get_values(dtype)
            itemmask[rl] = 1

        if not itemmask.all():
            raise AssertionError('Some items were not contained in blocks')

        return result

    def xs(self, key, axis=1, copy=True, takeable=False):
        if axis < 1:
            raise AssertionError('Can only take xs across axis >= 1, got %d'
                                 % axis)

        # take by position
        if takeable:
            loc = key
        else:
            loc = self.axes[axis].get_loc(key)

        slicer = [slice(None, None) for _ in range(self.ndim)]
        slicer[axis] = loc
        slicer = tuple(slicer)

        new_axes = list(self.axes)

        # could be an array indexer!
        if isinstance(loc, (slice, np.ndarray)):
            new_axes[axis] = new_axes[axis][loc]
        else:
            new_axes.pop(axis)

        new_blocks = []
        if len(self.blocks) > 1:
            # we must copy here as we are mixed type
            for blk in self.blocks:
                newb = make_block(values=blk.values[slicer],
                                  klass=blk.__class__, fastpath=True,
                                  placement=blk.ref_locs)
                new_blocks.append(newb)
        elif len(self.blocks) == 1:
            block = self.blocks[0]
            vals = block.values[slicer]
            if copy:
                vals = vals.copy()
            new_blocks = [make_block(values=vals, placement=block.ref_locs,
                                     klass=block.__class__, fastpath=True,)]

        return self.__class__(new_blocks, new_axes)

    def fast_xs(self, loc):
        """
        get a cross sectional for a given location in the
        items ; handle dups

        return the result, is *could* be a view in the case of a
        single block
        """
        if len(self.blocks) == 1:
            return self.blocks[0].values[:, loc]

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
        result = np.empty(n, dtype=dtype)
        for blk in self.blocks:
            # Such assignment may incorrectly coerce NaT to None
            # result[blk.ref_locs] = blk._slice((slice(None), loc))
            for i, rl in enumerate(blk.ref_locs):
                result[rl] = blk._try_coerce_result(blk.iget((i, loc)))

        return result

    def consolidate(self):
        """
        Join together blocks having same dtype

        Returns
        -------
        y : BlockManager
        """
        if self.is_consolidated():
            return self

        bm = self.__class__(self.blocks, self.axes)
        bm._consolidate_inplace()
        return bm

    def _consolidate_inplace(self):
        if not self.is_consolidated():
            self.blocks = _consolidate(self.blocks)

            self._is_consolidated = True
            self._known_consolidated = True
            self._set_has_sparse()
            self._rebuild_ref_locs()

    def get(self, item):
        """
        Return values for selected item (ndarray or BlockManager).
        """
        if self.items.is_unique:

            if not isnull(item):
                loc = self.items.get_loc(item)
            else:
                indexer = np.arange(len(self.items))[isnull(self.items)]

                # allow a single nan location indexer
                if not np.isscalar(indexer):
                    if len(indexer) == 1:
                        loc = indexer.item()
                    else:
                        raise ValueError("cannot label index with a null key")

            return self.iget(loc)
        else:

            if isnull(item):
                raise ValueError("cannot label index with a null key")

            indexer = self.items.get_indexer_for([item])
            return self.reindex_indexer(new_axis=self.items[indexer],
                                        indexer=indexer, axis=0, allow_dups=True)

    def iget(self, i):
        b, loc = self._ref_locs[i]
        return b.iget(loc)

    def get_scalar(self, tup):
        """
        Retrieve single item
        """
        full_loc = list(ax.get_loc(x)
                        for ax, x in zip(self.axes, tup))
        blk, blk_loc = self._ref_locs[full_loc[0]]
        full_loc[0] = blk_loc
        return blk.values[tuple(full_loc)]

    def delete(self, item):
        """
        Delete selected item (items if non-unique) in-place.
        """
        indexer = self.items.get_loc(item)

        is_deleted = np.zeros(self.shape[0], dtype=np.bool_)
        is_deleted[indexer] = True
        ref_loc_offset = is_deleted.cumsum()

        new_items = self.items[~is_deleted]
        new_blocks = []

        for blk in self.blocks:
            brl = blk.ref_locs
            blk_del = is_deleted[brl]
            blk_del_count = np.count_nonzero(blk_del)

            if blk_del_count == len(brl):
                continue

            blk._ref_locs -= ref_loc_offset[brl]
            if blk_del_count != 0:
                blk = blk._getitem_block(~blk_del)

            new_blocks.append(blk)

        self.axes[0] = new_items
        self.blocks = new_blocks
        self._shape = None
        self._rebuild_ref_locs()

    def set(self, item, value, check=False):
        """
        Set new item in-place. Does not consolidate. Adds new Block if not
        contained in the current set of items
        if check, then validate that we are not setting the same data in-place
        """
        # FIXME: refactor, clearly separate broadcasting & zip-like assignment
        is_sparse = isinstance(value, SparseArray)

        if is_sparse:
            assert self.ndim == 2

            def value_getitem(locs):
                return value
        else:
            if value.ndim == self.ndim - 1:
                value = value.reshape((1,) + value.shape)

                def value_getitem(locs):
                    return value
            else:
                def value_getitem(locs):
                    return value[locs]
            if value.shape[1:] != self.shape[1:]:
                raise AssertionError('Shape of new values must be compatible '
                                     'with manager shape')

        try:
            loc = self.items.get_loc(item)
        except KeyError:
            # This item wasn't present, just insert at end
            self.insert(len(self.items), item, value)
            return

        if isinstance(loc, int):
            loc = [loc]

        ref_locs = self._ref_locs[loc]

        unfit_mgr_locs = []
        unfit_val_locs = []
        for blk, blk_locs, val_locs in ref_loc_groupby_block(ref_locs):
            if blk.should_store(value):
                blk.set(blk_locs, value_getitem(val_locs), check=check)
            else:
                unfit_mgr_locs.append(blk.ref_locs[blk_locs])
                unfit_val_locs.append(val_locs)

                new_blk_ref_locs = np.delete(blk.ref_locs, blk_locs, axis=0)
                new_blk_len = len(new_blk_ref_locs)
                if not new_blk_len:
                    self.blocks.remove(blk)
                else:
                    blk.values = np.delete(blk.values, blk_locs, axis=0)
                    blk._ref_locs = new_blk_ref_locs
                    self._ref_locs[new_blk_ref_locs] = \
                        lib.fast_zip([np.array([blk] * new_blk_len),
                                      np.arange(new_blk_len)])

        if unfit_val_locs:
            unfit_val_locs = np.concatenate(unfit_val_locs)
            unfit_mgr_locs = np.concatenate(unfit_mgr_locs)
            unfit_count = len(unfit_val_locs)

            if is_sparse:
                for mgr_loc in unfit_mgr_locs:
                    new_block = make_block(values=value.copy(),
                                           ndim=self.ndim,
                                           placement=[mgr_loc])
                    self.blocks.append(new_block)
                    self._ref_locs[mgr_loc] = (new_block, 0)
            else:
                new_block = make_block(values=value[unfit_val_locs],
                                       ndim=self.ndim,
                                       placement=unfit_mgr_locs)

                self.blocks.append(new_block)
                self._ref_locs[unfit_mgr_locs] = lib.fast_zip([
                    np.array([new_block] * unfit_count, dtype=np.object_),
                    np.arange(unfit_count)])

            # Newly created block's dtype may already be present.
            self._known_consolidated = False

    def insert(self, loc, item, value, allow_duplicates=False):
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
            raise ValueError('cannot insert %s, already exists' % item)

        if not isinstance(loc, int):
            raise TypeError("loc must be int")

        new_items = self.items.insert(loc, item)
        block = make_block(values=value,
                           ndim=self.ndim,
                           placement=[loc])
        new_ref_locs = np.insert(self._ref_locs, loc, None, axis=0)
        new_ref_locs[loc] = (block, 0)

        for blk in self.blocks:
            blk._ref_locs[blk._ref_locs >= loc] += 1

        self.blocks.append(block)
        self.axes[0] = new_items
        self._shape = None
        self._ref_locs = new_ref_locs

        self._known_consolidated = False

        if len(self.blocks) > 100:
            self._consolidate_inplace()

    def reindex_axis(self, new_axis, axis, method=None, limit=None,
                     fill_value=None, copy=True):
        mgr = self if not copy else self.copy(deep=True)

        new_axis = _ensure_index(new_axis)
        new_axis, indexer = mgr.axes[axis].reindex(
            new_axis, method=method, limit=limit, copy_if_needed=True)

        return mgr.reindex_indexer(new_axis, indexer, axis=axis,
                                   fill_value=fill_value)

    def reindex_indexer(self, new_axis, indexer, axis, fill_value=None,
                        allow_dups=False):
        """
        pandas-indexer with -1's only.
        """
        # trying to reindex on an axis with duplicates
        if (not allow_dups and not self.axes[axis].is_unique
            and indexer is not None and len(indexer)):
            raise ValueError("cannot reindex from a duplicate axis")

        if axis >= self.ndim:
            raise AssertionError("Requested axis not found in manager")

        # FIXME: this code comes from generic.py, see if any of that is needed
        # elif (baxis == 0 and
        #         index is not new_data.axes[baxis]):
        #     new_data = new_data.reindex_items(index, copy=copy,
        #                                       fill_value=fill_value)

        # elif (baxis > 0 and index is not None and
        #         index is not new_data.axes[baxis]):
        #     new_data = new_data.copy(deep=copy)
        #     new_data.set_axis(baxis, index)

        if axis == 0:
            new_blocks = self._get_blocks_for_items_indexer(indexer,
                                                            fill_value)
        else:
            # TODO: is this faster than blk.reindex_axis?
            # return self.apply('take',
            #                   axes=new_axes,
            #                   indexer=indexer,
            #                   ref_items=new_axes[0],
            #                   new_axis=new_axes[axis],
            #                   axis=axis)
            new_blocks = [blk.reindex_axis(indexer, axis=axis,
                                           fill_value=fill_value)
                          for blk in self.blocks]

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return self.__class__(new_blocks, new_axes)

    def _get_blocks_for_items_indexer(self, indexer, fill_value):
        """
        Reindex blocks at axis=0 (overloaded for SingleBlockManager).

        Returns
        -------
        new_blocks : list of Block

        """
        # fill_value[0] == None will group soon-to-be-added items under None
        # fill_value[1] is an arbitrary integer (it's ignored)
        new_ref_locs = com.take_1d(self._ref_locs, indexer,
                                   fill_value=(None, 0))
        new_blocks = []
        for blk, blk_locs, mgr_locs in ref_loc_groupby_block(new_ref_locs):
            if blk is None:
                new_blocks.append(self._make_na_block(
                    placement=mgr_locs, fill_value=fill_value))
            else:
                # Otherwise, slicing along items axis is necessary.
                if blk.is_sparse:
                    # If it's a sparse block, it's easy:
                    #
                    # - it can only contain 1 item
                    # - if blk is here, the item wasn't deleted
                    # - if blk wasn't handled above, the item is multiplied
                    #
                    # Hence the block is replicated.
                    for mgr_loc in mgr_locs:
                        newblk = blk.copy(deep=True)
                        newblk._ref_locs = np.array([mgr_loc])
                        new_blocks.append(newblk)

                else:
                    # FIXME: this hack makes sure post-reindex blocks enumerate
                    # manager locs in ascending order.  It was implemented to
                    # make pytables serialization test happy and should be
                    # removed once the codebase successfully switches to
                    # axis-oblivious blocks & blockmanagers.
                    order = np.argsort(mgr_locs)
                    blk_locs = blk_locs.take(order)
                    mgr_locs = mgr_locs.take(order)

                    new_values = com.take_1d(blk.values, blk_locs,
                                             axis=0, allow_fill=False)
                    newblk = blk.__class__(values=new_values,
                                           ndim=blk.ndim,
                                           fastpath=True,
                                           placement=mgr_locs,)
                    new_blocks.append(newblk)

        return new_blocks

    def _make_na_block(self, placement, fill_value=None):
        # TODO: infer dtypes other than float64 from fill_value

        if fill_value is None:
            fill_value = np.nan
        block_shape = list(self.shape)
        block_shape[0] = len(placement)

        dtype, fill_value = com._infer_dtype_from_scalar(fill_value)
        block_values = np.empty(block_shape, dtype=dtype)
        block_values.fill(fill_value)
        return make_block(block_values, placement=placement)

    def take(self, indexer, axis=1, verify=True, convert=True):
        """
        Take items along any axis.
        """
        self._consolidate_inplace()
        indexer = np.asanyarray(indexer, dtype=np.int_)

        n = self.shape[axis]
        if convert:
            indexer = _maybe_convert_indices(indexer, n)

        if verify:
            if ((indexer == -1) | (indexer >= n)).any():
                raise Exception('Indices must be nonzero and less than '
                                'the axis length')

        new_labels = self.axes[axis].take(indexer)
        return self.reindex_indexer(new_axis=new_labels, indexer=indexer,
                                    axis=axis, allow_dups=True)

    def merge(self, other, lsuffix='', rsuffix=''):
        if not self._is_indexed_like(other):
            raise AssertionError('Must have same axes to merge managers')

        l, r = items_overlap_with_suffix(left=self.items, lsuffix=lsuffix,
                                         right=other.items, rsuffix=rsuffix)
        new_items = _concat_indexes([l, r])

        new_blocks = []
        for blocks, offset in [(self.blocks, 0),
                               (other.blocks, self.shape[0])]:
            for blk in blocks:
                blk = blk.copy(deep=False)
                blk._ref_locs += offset
                new_blocks.append(blk)

        new_axes = list(self.axes)
        new_axes[0] = new_items

        return self.__class__(_consolidate(new_blocks), new_axes)

    def _is_indexed_like(self, other):
        """
        Check all axes except items
        """
        if self.ndim != other.ndim:
            raise AssertionError(('Number of dimensions must agree '
                                  'got %d and %d') % (self.ndim, other.ndim))
        for ax, oax in zip(self.axes[1:], other.axes[1:]):
            if not ax.equals(oax):
                return False
        return True

    def rename_axis(self, mapper, axis, copy=True):
        """
        Rename one of axes.

        Parameters
        ----------
        mapper : unary callable
        axis : int
        copy : boolean, default True

        """
        new_axis = _transform_index(self.axes[axis], mapper)

        if axis != 0:
            new_blocks = self.blocks
        else:
            new_blocks = []
            for block in self.blocks:
                newb = block.copy(deep=copy)
                new_blocks.append(newb)

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return self.__class__(new_blocks, new_axes)

    def add_prefix(self, prefix):
        f = (('%s' % prefix) + '%s').__mod__
        return self.rename_axis(f, axis=0)

    def add_suffix(self, suffix):
        f = ('%s' + ('%s' % suffix)).__mod__
        return self.rename_axis(f, axis=0)

    def equals(self, other):
        self_axes, other_axes = self.axes, other.axes
        if len(self_axes) != len(other_axes):
           return False
        if not all (ax1.equals(ax2) for ax1, ax2 in zip(self_axes, other_axes)):
            return False
        self._consolidate_inplace()
        other._consolidate_inplace()
        return all(block.equals(oblock) for block, oblock in
                   zip(self.blocks, other.blocks))

    def group_blocks_by_ftype(self):
        """
        Combine blocks into map: ftype -> [blk0, blk1, ...].

        """
        bm = defaultdict(list)
        for b in self.blocks:
            bm[str(b.ftype)].append(b)
        return bm


class SingleBlockManager(BlockManager):

    """ manage a single block with """
    ndim = 1
    _is_consolidated = True
    _known_consolidated = True
    __slots__ = ['axes', 'blocks']

    def __init__(self, block, axis, do_integrity_check=False, fastpath=True):

        if isinstance(axis, list):
            if len(axis) != 1:
                raise ValueError(
                    "cannot create SingleBlockManager with more than 1 axis")
            axis = axis[0]

        # passed from constructor, single block, single axis
        if fastpath:
            self.axes = [axis]
            if isinstance(block, list):

                # empty block
                if len(block) == 0:
                    block = [np.array([])]
                elif len(block) != 1:
                    raise ValueError('Cannot create SingleBlockManager with '
                                     'more than 1 block')
                block = block[0]
            if not isinstance(block, Block):
                block = make_block(block, ndim=1, fastpath=True,
                                   placement=np.arange(len(axis)))

        else:

            self.axes = [_ensure_index(axis)]

            # create the block here
            if isinstance(block, list):

                # provide consolidation to the interleaved_dtype
                if len(block) > 1:
                    dtype = _interleaved_dtype(block)
                    block = [b.astype(dtype) for b in block]
                    block = _consolidate(block)

                if len(block) != 1:
                    raise ValueError('Cannot create SingleBlockManager with '
                                     'more than 1 block')
                block = block[0]

            if not isinstance(block, Block):
                block = make_block(block, axis, ndim=1,
                                   fastpath=True, placement=None)

        self.blocks = [block]

    def _post_setstate(self):
        pass

    @property
    def _block(self):
        return self.blocks[0]

    @property
    def _values(self):
        return self._block.values

    @property
    def _has_sparse(self):
        return self._block.is_sparse

    def _set_has_sparse(self):
        # _has_sparse is a property, nothing to set here
        pass

    # def apply(self, f, axes=None, do_integrity_check=False, **kwargs):
    #     """
    #     fast path for SingleBlock Manager

    #     ssee also BlockManager.apply
    #     """
    #     applied = getattr(self._block, f)(**kwargs)
    #     bm = self.__class__(applied, axes or self.axes,
    #                         do_integrity_check=do_integrity_check)
    #     bm._consolidate_inplace()
    #     return bm

    def reindex(self, new_axis, indexer=None, method=None, fill_value=None,
                limit=None, copy=True):
        # if we are the same and don't copy, just return
        if self.index.equals(new_axis):
            if copy:
                return self.copy(deep=True)
            else:
                return self

        values = self._block.get_values()

        if indexer is None:
            indexer = self.items.get_indexer_for(new_axis)

        if fill_value is None:
            # FIXME: is fill_value used correctly in sparse blocks?
            if not self._block.is_sparse:
                fill_value = self._block.fill_value
            else:
                fill_value = np.nan

        new_values = com.take_1d(values, indexer,
                                 fill_value=fill_value)

        # fill if needed
        if method is not None or limit is not None:
            new_values = com.interpolate_2d(new_values, method=method,
                                            limit=limit, fill_value=fill_value)

        if self._block.is_sparse:
            make_block = self._block.make_block

        block = make_block(new_values, copy=copy,
                           placement=np.arange(len(new_axis)))

        # block = self._block.reindex_items_from(new_axis, indexer=indexer,
        #                                        method=method,
        #                                        fill_value=fill_value,
        #                                        limit=limit, copy=copy)
        mgr = SingleBlockManager(block, new_axis)
        mgr._consolidate_inplace()
        return mgr

    def _reindex_indexer_items(self, new_items, indexer, fill_value):
        # equiv to a reindex
        return self.reindex(new_items, indexer=indexer, fill_value=fill_value,
                            copy=False)

    def _delete_from_block(self, i, item):
        super(SingleBlockManager, self)._delete_from_block(i, item)

        # possibly need to merge split blocks
        if len(self.blocks) > 1:
            new_values = np.concatenate([b.values for b in self.blocks])
            new_items = Index(np.concatenate([b.items for b in self.blocks]))

            block = make_block(values=new_values, placement=None,
                               dtype=self._block.dtype,)

        elif len(self.blocks):
            block = self.blocks[0]
        else:
            block = make_block(values=np.array([], dtype=self._block.dtype),
                               placement=None)

        self.blocks = [block]

    def get_slice(self, slobj):
        return self.__class__(self._block._slice(slobj),
                              self.index[slobj], fastpath=True)

    @property
    def index(self):
        return self.axes[0]

    def convert(self, **kwargs):
        """ convert the whole block as one """
        kwargs['by_item'] = False
        return self.apply('convert', **kwargs)

    @property
    def dtype(self):
        return self._values.dtype

    @property
    def ftype(self):
        return self._block.ftype

    def get_dtype_counts(self):
        return {self.dtype.name: 1}

    def get_ftype_counts(self):
        return {self.ftype: 1}

    def get_dtypes(self):
        return [self._block.dtype]

    def get_ftypes(self):
        return [self._block.ftype]

    @property
    def values(self):
        return self._values.view()

    @property
    def itemsize(self):
        return self._values.itemsize

    @property
    def _can_hold_na(self):
        return self._block._can_hold_na

    def is_consolidated(self):
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
        # Also, make sure dtype is preserved.
        dtype = self._block.dtype

        super(SingleBlockManager, self).delete(item)

        if not self.blocks:
            self.blocks = [make_block(values=np.empty(0, dtype=dtype),
                                      placement=np.arange(len(self.items)),
                                      ndim=1, dtype=dtype, fastpath=True)]

    def fast_xs(self, loc):
        """
        fast path for getting a cross-section
        return a view of the data
        """
        return self._block.values[loc]

    def _get_blocks_for_items_indexer(self, indexer, fill_value):
        """
        Reindex blocks at axis=0 (overloaded for SingleBlockManager).

        Returns
        -------
        new_blocks : list of Block

        """
        if indexer is None:
            new_values = self._values.copy()
        else:
            new_values = com.take_1d(self._values, indexer,
                                     fill_value=fill_value)

        return [make_block(values=new_values,
                           placement=np.arange(len(new_values)),
                           ndim=self.ndim, fastpath=True)]


def construction_error(tot_items, block_shape, axes, e=None):
    """ raise a helpful message about our construction """
    passed = tuple(map(int, [tot_items] + list(block_shape)))
    implied = tuple(map(int, [len(ax) for ax in axes]))
    if passed == implied and e is not None:
        raise e
    raise ValueError("Shape of passed values is {0}, indices imply {1}".format(
        passed,implied))


def create_block_manager_from_blocks(blocks, axes):
    try:
        if len(blocks) == 1 and not isinstance(blocks[0], Block):
            # It's OK if a single block is passed as values, its placement is
            # basically "all items", but if there're many, don't bother
            # converting, it's an error anyway.
            blocks = [make_block(values=blocks[0],
                                 placement=np.arange(len(axes[0])),)]

        mgr = BlockManager(blocks, axes)
        mgr._consolidate_inplace()
        return mgr

    except (ValueError) as e:
        blocks = [getattr(b, 'values', b) for b in blocks]
        tot_items = sum(b.shape[0] for b in blocks)
        construction_error(tot_items, blocks[0].shape[1:], axes, e)


def create_block_manager_from_arrays(arrays, names, axes):
    try:
        blocks = form_blocks(arrays, names, axes)
        mgr = BlockManager(blocks, axes)
        mgr._consolidate_inplace()
        return mgr
    except (ValueError) as e:
        construction_error(len(arrays), arrays[0].shape[1:], axes, e)


def form_blocks(arrays, names, axes):
    # put "leftover" items in float bucket, where else?
    # generalize?
    float_items = []
    complex_items = []
    int_items = []
    bool_items = []
    object_items = []
    sparse_items = []
    datetime_items = []
    extra_locs = []

    names_idx = Index(names)
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

        if isinstance(v, (SparseArray, ABCSparseSeries)):
            sparse_items.append((i, k, v))
        elif issubclass(v.dtype.type, np.floating):
            float_items.append((i, k, v))
        elif issubclass(v.dtype.type, np.complexfloating):
            complex_items.append((i, k, v))
        elif issubclass(v.dtype.type, np.datetime64):
            if v.dtype != _NS_DTYPE:
                v = tslib.cast_to_nanoseconds(v)

            if hasattr(v, 'tz') and v.tz is not None:
                object_items.append((i, k, v))
            else:
                datetime_items.append((i, k, v))
        elif issubclass(v.dtype.type, np.integer):
            if v.dtype == np.uint64:
                # HACK #2355 definite overflow
                if (v > 2 ** 63 - 1).any():
                    object_items.append((i, k, v))
                    continue
            int_items.append((i, k, v))
        elif v.dtype == np.bool_:
            bool_items.append((i, k, v))
        else:
            object_items.append((i, k, v))

    blocks = []
    if len(float_items):
        float_blocks = _multi_blockify(float_items)
        blocks.extend(float_blocks)

    if len(complex_items):
        complex_blocks = _simple_blockify(
            complex_items, np.complex128)
        blocks.extend(complex_blocks)

    if len(int_items):
        int_blocks = _multi_blockify(int_items)
        blocks.extend(int_blocks)

    if len(datetime_items):
        datetime_blocks = _simple_blockify(
            datetime_items, _NS_DTYPE)
        blocks.extend(datetime_blocks)

    if len(bool_items):
        bool_blocks = _simple_blockify(
            bool_items, np.bool_)
        blocks.extend(bool_blocks)

    if len(object_items) > 0:
        object_blocks = _simple_blockify(
            object_items, np.object_)
        blocks.extend(object_blocks)

    if len(sparse_items) > 0:
        sparse_blocks = _sparse_blockify(sparse_items)
        blocks.extend(sparse_blocks)

    if len(extra_locs):
        shape = (len(extra_locs),) + tuple(len(x) for x in axes[1:])

        # empty items -> dtype object
        block_values = np.empty(shape, dtype=object)
        block_values.fill(np.nan)

        na_block = make_block(block_values, placement=extra_locs)
        blocks.append(na_block)

    return blocks


def _simple_blockify(tuples, dtype):
    """ return a single array of a block that has a single dtype; if dtype is
    not None, coerce to this dtype
    """
    values, placement = _stack_arrays(tuples, dtype)

    # CHECK DTYPE?
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

        values, placement = _stack_arrays(
            list(tup_block), dtype)

        block = make_block(values, placement=placement)
        new_blocks.append(block)

    return new_blocks


def _sparse_blockify(tuples, dtype=None):
    """ return an array of blocks that potentially have different dtypes (and
    are sparse)
    """

    new_blocks = []
    for i, names, array in tuples:
        array = _maybe_to_sparse(array)
        block = make_block(
            array, klass=SparseBlock, fastpath=True,
            placement=[i])
        new_blocks.append(block)

    return new_blocks


def _stack_arrays(tuples, dtype):

    # fml
    def _asarray_compat(x):
        if isinstance(x, ABCSeries):
            return x.values
        else:
            return np.asarray(x)

    def _shape_compat(x):
        if isinstance(x, ABCSeries):
            return len(x),
        else:
            return x.shape

    placement, names, arrays = zip(*tuples)

    first = arrays[0]
    shape = (len(arrays),) + _shape_compat(first)

    stacked = np.empty(shape, dtype=dtype)
    for i, arr in enumerate(arrays):
        stacked[i] = _asarray_compat(arr)

    return stacked, placement


def _interleaved_dtype(blocks):
    if not len(blocks):
        return None

    counts = defaultdict(lambda: [])
    for x in blocks:
        counts[type(x)].append(x)

    def _lcd_dtype(l):
        """ find the lowest dtype that can accomodate the given types """
        m = l[0].dtype
        for x in l[1:]:
            if x.dtype.itemsize > m.itemsize:
                m = x.dtype
        return m

    have_int = len(counts[IntBlock]) > 0
    have_bool = len(counts[BoolBlock]) > 0
    have_object = len(counts[ObjectBlock]) > 0
    have_float = len(counts[FloatBlock]) > 0
    have_complex = len(counts[ComplexBlock]) > 0
    have_dt64 = len(counts[DatetimeBlock]) > 0
    have_td64 = len(counts[TimeDeltaBlock]) > 0
    have_sparse = len(counts[SparseBlock]) > 0
    have_numeric = have_float or have_complex or have_int

    if (have_object or
        (have_bool and have_numeric) or
            (have_numeric and (have_dt64 or have_td64))):
        return np.dtype(object)
    elif have_bool:
        return np.dtype(bool)
    elif have_int and not have_float and not have_complex:

        # if we are mixing unsigned and signed, then return
        # the next biggest int type (if we can)
        lcd = _lcd_dtype(counts[IntBlock])
        kinds = set([i.dtype.kind for i in counts[IntBlock]])
        if len(kinds) == 1:
            return lcd

        if lcd == 'uint64' or lcd == 'int64':
            return np.dtype('int64')

        # return 1 bigger on the itemsize if unsinged
        if lcd.kind == 'u':
            return np.dtype('int%s' % (lcd.itemsize * 8 * 2))
        return lcd

    elif have_dt64 and not have_float and not have_complex:
        return np.dtype('M8[ns]')
    elif have_td64 and not have_float and not have_complex:
        return np.dtype('m8[ns]')
    elif have_complex:
        return np.dtype('c16')
    else:
        return _lcd_dtype(counts[FloatBlock] + counts[SparseBlock])


def _consolidate(blocks):
    """
    Merge blocks having same dtype, exclude non-consolidating blocks
    """

    # sort by _can_consolidate, dtype
    gkey = lambda x: x._consolidate_key
    grouper = itertools.groupby(sorted(blocks, key=gkey), gkey)

    new_blocks = []
    for (_can_consolidate, dtype), group_blocks in grouper:
        merged_blocks = _merge_blocks(list(group_blocks), dtype=dtype,
                                      _can_consolidate=_can_consolidate)
        if isinstance(merged_blocks, list):
            new_blocks.extend(merged_blocks)
        else:
            new_blocks.append(merged_blocks)

    return new_blocks


def _merge_blocks(blocks, dtype=None, _can_consolidate=True):
    if len(blocks) == 1:
        return blocks[0]

    if _can_consolidate:

        if dtype is None:
            if len(set([b.dtype for b in blocks])) != 1:
                raise AssertionError("_merge_blocks are invalid!")
            dtype = blocks[0].dtype

        new_ref_locs = np.concatenate([b.ref_locs for b in blocks])
        new_values = _vstack([b.values for b in blocks], dtype)

        argsort = np.argsort(new_ref_locs)
        new_values = new_values[argsort]
        new_ref_locs = new_ref_locs[argsort]

        return make_block(new_values,
                          fastpath=True, placement=new_ref_locs)

    # no merge
    return blocks


def _block_shape(values, ndim=1, shape=None):
    """ guarantee the shape of the values to be at least 1 d """
    if values.ndim <= ndim:
        if shape is None:
            shape = values.shape
        values = values.reshape(tuple((1,) + shape))
    return values


def _vstack(to_stack, dtype):

    # work around NumPy 1.6 bug
    if dtype == _NS_DTYPE or dtype == _TD_DTYPE:
        new_values = np.vstack([x.view('i8') for x in to_stack])
        return new_values.view(dtype)

    else:
        return np.vstack(to_stack)


def _possibly_convert_to_indexer(loc):
    if com._is_bool_indexer(loc):
        loc = [i for i, v in enumerate(loc) if v]
    elif isinstance(loc, slice):
        loc = lrange(loc.start, loc.stop)
    return loc


def _possibly_compare(a, b, op):
    res = op(a, b)
    is_a_array = isinstance(a, np.ndarray)
    is_b_array = isinstance(b, np.ndarray)
    if np.isscalar(res) and (is_a_array or is_b_array):
        type_names = [type(a).__name__, type(b).__name__]

        if is_a_array:
            type_names[0] = 'ndarray(dtype=%s)' % a.dtype

        if is_b_array:
            type_names[1] = 'ndarray(dtype=%s)' % b.dtype

        raise TypeError("Cannot compare types %r and %r" % tuple(type_names))
    return res




def _concat_indexes(indexes):
    return indexes[0].append(indexes[1:])


def _invert_reordering(reordering, minlength=None):
    """
    Invert reordering operation.

    Given array `reordering`, make `reordering_inv` of it, such that::

        reordering_inv[reordering[x]] = x

    There are two types of indexers:

    source
        is when element *s* at position *i* means that values to fill *i-th*
        item of reindex operation should be taken from *s-th* item of the
        original (this is what is returned by `pandas.Index.reindex`).
    destination
        is when element *d* at position *i* means that values from *i-th* item
        of source should be used to fill *d-th* item of reindexing operation.

    This function will convert from *source* to *destination* and vice-versa.

    .. note:: trailing ``-1`` may be lost upon conversion (this is what
              `minlength` is there for).

    .. note:: if *source* indexer is not unique, corresponding *destination*
              indexer will have ``dtype=object`` and will contain lists.

    Examples:

    >>> _invert_reordering([3, -1, 2, 4, -1])
    array([-1, -1,  2,  0,  3])
    >>> _invert_reordering([-1, -1, 0, 2, 3])
    array([3, -1,  2,  4])
    >>> _invert_reordering([1,3,5])
    array([-1,  0, -1,  1, -1,  2])

    """
    reordering = np.asanyarray(reordering)
    if not com.is_integer_dtype(reordering):
        raise ValueError("Only integer indexers are supported")

    nonneg_indices = reordering[reordering >= 0]
    counts = np.bincount(nonneg_indices, minlength=minlength)
    has_non_unique = (counts > 1).any()

    dtype = np.dtype(np.object_) if has_non_unique else np.dtype(np.int_)
    inverted = np.empty_like(counts, dtype=dtype)
    inverted.fill(-1)

    nonneg_positions = np.arange(len(reordering), dtype=np.int_)[reordering >= 0]
    np.put(inverted, nonneg_indices, nonneg_positions)

    if has_non_unique:
        nonunique_elements = np.arange(len(counts))[counts > 1]
        for elt in nonunique_elements:
            inverted[elt] = nonneg_positions[nonneg_indices == elt].tolist()

    return inverted


def ref_loc_groupby_block(ref_locs):
    """
    Group given ref_locs by block.

    Returns
    -------
    iterator
        Yield (block, block_locs, original_locs)

    """
    if len(ref_locs) == 0:
        return

    blocks = com._ensure_object(lib.map_infer(ref_locs,
                                              operator.itemgetter(0)))
    indices = lib.map_infer(ref_locs, operator.itemgetter(1))

    factorizer = Factorizer(len(blocks))
    block_ids = factorizer.factorize(blocks, na_sentinel=-1)

    for i in range(factorizer.get_count()):
        locs = (block_ids == i).nonzero()[0]
        yield blocks[locs[0]], indices[locs], locs

    na_locs = (block_ids == -1).nonzero()[0]
    if len(na_locs):
        yield None, indices[na_locs], na_locs


def items_overlap_with_suffix(left, lsuffix, right, rsuffix):
    """
    If two indices overlap, add suffixes to overlapping entries.

    If corresponding suffix is empty, the entry is simply converted to string.

    """
    to_rename = left.intersection(right)
    if len(to_rename) == 0:
        return left, right
    else:
        if not lsuffix and not rsuffix:
            raise ValueError('columns overlap but no suffix specified: %s' %
                             to_rename)

        def lrenamer(x):
            if x in to_rename:
                return '%s%s' % (x, lsuffix)
            return x

        def rrenamer(x):
            if x in to_rename:
                return '%s%s' % (x, rsuffix)
            return x

        return (_transform_index(left, lrenamer),
                _transform_index(right, rrenamer))


def _transform_index(index, func):
    """
    Apply function to all values found in index.

    This includes transforming multiindex entries separately.

    """
    if isinstance(index, MultiIndex):
        items = [tuple(func(y) for y in x) for x in index]
        return MultiIndex.from_tuples(items, names=index.names)
    else:
        items = [func(x) for x in index]
        return Index(items, name=index.name)


def _putmask_smart(v, m, n):
    """
    Return a new block, try to preserve dtype if possible.

    Parameters
    ----------
    v : array_like
    m : array_like
    n : array_like
    """

    # n should be the length of the mask or a scalar here
    if not is_list_like(n):
        n = np.array([n] * len(m))

    # see if we are only masking values that if putted
    # will work in the current dtype
    try:
        nn = n[m]
        nn_at = nn.astype(v.dtype)
        if (nn == nn_at).all():
            nv = v.copy()
            nv[m] = nn_at
            return nv
    except (ValueError, IndexError, TypeError):
        pass

    # change the dtype
    dtype, _ = com._maybe_promote(n.dtype)
    nv = v.astype(dtype)
    try:
        nv[m] = n
    except ValueError:
        idx, = np.where(np.squeeze(m))
        for mask_index, new_val in zip(idx, n):
            nv[mask_index] = new_val
    return nv


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
    concat_plans = []

    for mgr, indexers in mgrs_indexers:
        plan = get_mgr_concatenation_plan(mgr, indexers)
        concat_plans = combine_concat_plans(concat_plans, plan, concat_axis)

    blocks = [concatenate_by_plan(plan, concat_axis, copy=copy)
              for plan in concat_plans]

    return BlockManager(blocks, axes)


def get_empty_dtype_and_na(join_units):
    """
    Return dtype and N/A values to use when concatenating specified units.

    Returned N/A value may be None which means there was no casting involved.

    Returns
    -------
    dtype
    na
    """

    has_none_blocks = False
    dtypes = set()
    upcast_classes = set()
    null_upcast_classes = set()
    for unit in join_units:
        if unit.block is None:
            # This value is not supposed to be used anywhere, it's here to make
            # sure "monotype" check (len(dtypes) == 1) fails and to indicate
            # that upcasting is required.
            has_none_blocks = True
            continue

        dtype = unit.dtype
        dtypes.add(unit.dtype)

        if issubclass(dtype.type, (np.object_, np.bool_)):
            upcast_cls = 'object'
        elif is_datetime64_dtype(dtype):
            upcast_cls = 'datetime'
        elif is_timedelta64_dtype(dtype):
            upcast_cls = 'timedelta'
        else:
            upcast_cls = 'float'

        # Null blocks should not influence upcast class selection, unless there
        # are only null blocks, when same upcasting rules must be applied to
        # null upcast classes.
        if unit.is_null:
            null_upcast_classes.add(upcast_cls)
        else:
            upcast_classes.add(upcast_cls)

    if not has_none_blocks and len(dtypes) == 1:
        # Unanimous decision, nothing to upcast.
        return next(iter(dtypes)), None

    if not upcast_classes:
        upcast_classes = null_upcast_classes

    # create the result
    if 'object' in upcast_classes:
        return np.dtype(np.object_), np.nan
    elif 'float' in upcast_classes:
        return np.dtype(np.float64), np.nan
    elif 'datetime' in upcast_classes:
        return np.dtype('M8[ns]'), tslib.iNaT
    elif 'timedelta' in upcast_classes:
        return np.dtype('m8[ns]'), tslib.iNaT
    else:  # pragma
        raise AssertionError("invalid dtype determination in get_concat_dtype")


def concatenate_by_plan(plan, concat_axis, copy):
    """
    Make block from concatenation plan.
    """
    concat_start, join_units = plan

    empty_dtype, upcasted_na = get_empty_dtype_and_na(join_units)

    to_concat = [ju.get_reindexed_values(empty_dtype=empty_dtype,
                                         upcasted_na=upcasted_na)
                 for ju in join_units]

    if len(to_concat) == 1:
        # Only one block, nothing to concatenate.
        if copy:
            concat_values = to_concat[0].copy()
        else:
            concat_values = to_concat[0]
    else:
        concat_values = com._concat_compat(to_concat, axis=concat_axis)

    rng = np.arange(concat_values.shape[0])

    if any(unit.is_sparse for unit in join_units):
        concat_values = SparseArray(concat_values[0])

    return make_block(concat_values,
                      placement=rng + concat_start)


def get_mgr_concatenation_plan(mgr, indexers):
    """
    Construct concatenation plan for given block manager and indexers.

    Parameters
    ----------
    mgr : BlockManager
    indexers : dict of {axis: indexer}

    Returns
    -------
    plan : list of (start_loc, [JoinUnit]) tuples

    """
    # Calculate post-reindex shape , save for item axis which will be separate
    # for each block anyway.
    mgr_shape = list(mgr.shape)
    for ax, indexer in indexers.items():
        mgr_shape[ax] = len(indexer)

    if 0 in indexers:
        indexer = indexers.pop(0)
        ref_locs = com.take_1d(mgr._ref_locs, indexer, fill_value=(None, 0))
    else:
        ref_locs = mgr._ref_locs

    plan = []
    for blk, blk_locs, concat_locs in ref_loc_groupby_block(ref_locs):
        # result_locs are assumed to be sorted
        slices = locs_to_contiguous_sequences(concat_locs)

        for slc in slices:
            join_unit_indexers = indexers.copy()
            axis0_blk_indexer = blk_locs[slc]

            # Omit indexer if no item reindexing is required.
            if (blk is None or
                np.array_equal(axis0_blk_indexer, np.arange(blk.shape[0]))):
                join_unit_indexers.pop(0, None)
            else:
                join_unit_indexers[0] = axis0_blk_indexer

            blk_shape = copy.copy(mgr_shape)
            blk_shape[0] = len(axis0_blk_indexer)
            unit = JoinUnit(blk, join_unit_indexers, shape=blk_shape)

            plan.append((concat_locs[slc.start], [unit]))

    plan.sort()
    return plan


def combine_concat_plans(existing_plan, new_plan, concat_axis):
    """
    Combine multiple concatenation plans into one.

    existing_plan is updated in-place.
    """
    if not existing_plan:
        # Shortcut: nothing to combine with
        return new_plan

    if concat_axis == 0:
        # Another shortcut: when concatenating along item axis, plans can be
        # simply appended.
        last_offset, last_units = existing_plan[-1]
        plan_offset = last_offset + last_units[0].shape[0]
        return existing_plan + [(off_i + plan_offset, units_i)
                                for off_i, units_i in new_plan]

    from collections import deque
    old_items = deque(existing_plan)
    new_items = deque(new_plan)
    result = []

    while new_items:
        old_start, old_units = old_items.popleft()
        new_start, new_units = new_items.popleft()

        assert old_start == new_start

        old_len = old_units[0].shape[0]
        new_len = new_units[0].shape[0]

        # Trim either old or new part as necessary
        common_len = min(old_len, new_len)
        if new_len > common_len:
            new_items.appendleft((new_start + common_len,
                                  [trim_join_unit(unit, common_len)
                                   for unit in new_units]))
        elif old_len > common_len:
            old_items.appendleft((old_start + common_len,
                                  [trim_join_unit(unit, common_len)
                                   for unit in old_units]))

        result.append((old_start, old_units + new_units))

    # The loop terminates when there's no new items, make sure that all old
    # items are processed.
    assert not old_items

    return result


def locs_to_contiguous_sequences(locs):
    """
    Return contiguous sequences found in locs as slices.
    """
    # FIXME: the code looks vaguely familiar, maybe there another version that
    # can be reused instead
    assert locs.ndim == 1
    length = len(locs)

    diff = np.diff(locs, axis=0)
    break_locs = (diff != 1).nonzero()[0] + 1

    if len(break_locs) == 0:
        return [slice(0, length)]
    else:
        return [slice(b, e)
                for b, e in lib.fast_zip([np.r_[0, break_locs],
                                          np.r_[break_locs, length]])]


def trim_join_unit(join_unit, length):
    """
    Reduce join_unit's shape along item axis to length.

    Extra items that didn't fit are returned as a separate block.
    """

    if 0 not in join_unit.indexers:
        join_unit.indexers[0] = np.arange(join_unit.shape[0])

    extra_indexers = copy.copy(join_unit.indexers)
    extra_shape = copy.copy(join_unit.shape)

    extra_shape[0] = join_unit.shape[0] - length
    extra_indexers[0] = extra_indexers[0][length:]

    join_unit.shape[0] = length
    join_unit.indexers[0] = join_unit.indexers[0][:length]

    return JoinUnit(block=join_unit.block, indexers=extra_indexers,
                    shape=extra_shape)


class JoinUnit(object):
    def __init__(self, block, indexers, shape):
        # Passing shape explicitly is required for cases when block is None.
        self.block = block
        self.indexers = indexers
        self.shape = shape

    def __repr__(self):
        return '%s(%r, %s)' % (self.__class__.__name__,
                               self.block, self.indexers)

    @cache_readonly
    def needs_filling(self):
        for indexer in self.indexers.values():
            # FIXME: cache results of indexer == -1 checks.
            if (indexer == -1).any():
                return True

        return False

    @cache_readonly
    def dtype(self):
        if self.block is None:
            raise AssertionError("Block is None, no dtype")

        if not self.needs_filling:
            return self.block.dtype
        else:
            return np.dtype(com._maybe_promote(self.block.dtype,
                                               self.block.fill_value)[0])
        return self._dtype

    @cache_readonly
    def is_null(self):
        return self.block is None or isnull(self.block.values).all()

    @cache_readonly
    def is_sparse(self):
        return self.block is not None and self.block.is_sparse

    def get_reindexed_values(self, empty_dtype, upcasted_na):
        if upcasted_na is not None:
            fill_value = upcasted_na
        else:
            # If upcasted_na is None, self.block should always exist.  If it
            # doesn't (i.e. is None), then it's a bug in get_empty_dtype_and_na
            # function.
            fill_value = self.block.fill_value

        if self.is_null:
            missing_arr = np.empty(self.shape, dtype=empty_dtype)
            if np.prod(self.shape):
                # NumPy 1.6 workaround: this statement gets strange if all
                # blocks are of same dtype and some of them are empty: empty
                # one are considered "null" so they must be filled, but no
                # dtype upcasting happens and the dtype may not allow NaNs.
                #
                # In general, no one should get hurt when one tries to put
                # incorrect values into empty array, but numpy 1.6 is strict
                # about that.
                missing_arr.fill(fill_value)
            return missing_arr
        else:
            if upcasted_na is not None and self.block.is_bool:
                # External code requested filling/upcasting, bool values must
                # be upcasted to object to avoid being upcasted to numeric.
                values = self.block.astype(np.object_).values
            else:
                values = self.block.get_values()

            for ax, indexer in self.indexers.items():
                values = com.take_nd(values, indexer, axis=ax,
                                     fill_value=fill_value)

            return values


# def _align_kwargs(blocks, items, kwargs, align_keys, copy):
#     aligned_objs = dict((k, kwargs[k]) for k in align_keys.items()
#                         if hasattr(kwargs[k], 'reindex_axis'))

#     if aligned_objs:
#         kwargs = kwargs.copy()

#     for b in blocks:
#         if aligned_objs:
#             b_items = items.take(b.ref_locs)

#             for k, obj in aligned_objs.items():
#                 axis = getattr(obj, '_info_axis_number', 0)
#                 kwargs[k] = obj.reindex_axis(b_items, axis=axis,
#                                              copy=copy)

#         yield b, kwargs
