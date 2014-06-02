import copy
import itertools
import re
import operator
from datetime import datetime, timedelta
from collections import defaultdict

import numpy as np
from pandas.core.base import PandasObject

from pandas.core.common import (_possibly_downcast_to_dtype, isnull,
                                _NS_DTYPE, _TD_DTYPE, ABCSeries, is_list_like,
                                ABCSparseSeries, _infer_dtype_from_scalar,
                                _is_null_datelike_scalar,
                                is_timedelta64_dtype, is_datetime64_dtype,)
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.indexing import (_maybe_convert_indices, _length_of_indexer)
import pandas.core.common as com
from pandas.sparse.array import _maybe_to_sparse, SparseArray
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.computation.expressions as expressions
from pandas.util.decorators import cache_readonly

from pandas.tslib import Timestamp
from pandas import compat
from pandas.compat import range, map, zip, u
from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type


from pandas.lib import BlockPlacement


class Block(PandasObject):

    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas
    data structure

    Index-ignorant; let the container take care of that
    """
    __slots__ = ['_mgr_locs', 'values', 'ndim']
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
        elif values.ndim != ndim:
            raise ValueError('Wrong number of dimensions')
        self.ndim = ndim

        self.mgr_locs = placement
        self.values = values

        if len(self.mgr_locs) != len(self.values):
            raise ValueError('Wrong number of items passed %d,'
                             ' placement implies %d' % (
                                 len(self.values), len(self.mgr_locs)))

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
    def mgr_locs(self):
        return self._mgr_locs

    def make_block_same_class(self, values, placement, copy=False,
                              **kwargs):
        """
        Wrap given values in a block of same type as self.

        `kwargs` are used in SparseBlock override.

        """
        if copy:
            values = values.copy()
        return make_block(values, placement, klass=self.__class__,
                          fastpath=True)

    @mgr_locs.setter
    def mgr_locs(self, new_mgr_locs):
        if not isinstance(new_mgr_locs, BlockPlacement):
            new_mgr_locs = BlockPlacement(new_mgr_locs)

        self._mgr_locs = new_mgr_locs

    def __unicode__(self):

        # don't want to print out all of the items here
        name = com.pprint_thing(self.__class__.__name__)
        if self._is_single_block:

            result = '%s: %s dtype: %s' % (
                name, len(self), self.dtype)

        else:

            shape = ' x '.join([com.pprint_thing(s) for s in self.shape])
            result = '%s: %s, %s, dtype: %s' % (
                name, com.pprint_thing(self.mgr_locs.indexer), shape,
                self.dtype)

        return result

    def __len__(self):
        return len(self.values)

    def __getstate__(self):
        return self.mgr_locs.indexer, self.values

    def __setstate__(self, state):
        self.mgr_locs = BlockPlacement(state[0])
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
            if isinstance(slicer, tuple):
                axis0_slicer = slicer[0]
            else:
                axis0_slicer = slicer
            new_mgr_locs = self.mgr_locs[axis0_slicer]

        new_values = self._slice(slicer)

        if new_values.ndim != self.ndim:
            raise ValueError("Only same dim slicing is allowed")

        return self.make_block_same_class(new_values, new_mgr_locs)

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
                          placement=self.mgr_locs)

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
        Delete given loc(-s) from block in-place.
        """
        self.values = np.delete(self.values, loc, 0)
        self.mgr_locs = self.mgr_locs.delete(loc)

    def apply(self, func, **kwargs):
        """ apply the function to my values; return a block if we are not one """
        result = func(self.values)
        if not isinstance(result, Block):
            result = make_block(values=result, placement=self.mgr_locs,)

        return result

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
                               fastpath=True, placement=self.mgr_locs)]

        # ndim > 1
        if dtypes is None:
            return [self]

        if not (dtypes == 'infer' or isinstance(dtypes, dict)):
            raise ValueError("downcast must have a dictionary or 'infer' as "
                             "its argument")

        # item-by-item
        # this is expensive as it splits the blocks items-by-item
        blocks = []
        for i, rl in enumerate(self.mgr_locs):

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
                              ndim=self.ndim, placement=self.mgr_locs,
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
                          placement=self.mgr_locs)

    def replace(self, to_replace, value, inplace=False, filter=None,
                regex=False):
        """ replace the to_replace value with value, possible to create new
        blocks here this is just a call to putmask. regex is not used here.
        It is used in ObjectBlocks.  It is here for API
        compatibility."""
        mask = com.mask_missing(self.values, to_replace)
        if filter is not None:
            filtered_out = ~self.mgr_locs.isin(filter)
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
                               ndim=self.ndim, placement=self.mgr_locs,
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
                for i, ref_loc in enumerate(self.mgr_locs):
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
                                             placement=self.mgr_locs,
                                             fastpath=True))

            return new_blocks

        if inplace:
            return [self]

        return [make_block(new_values,
                           placement=self.mgr_locs, fastpath=True)]

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
                             fastpath=True, placement=self.mgr_locs)]
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
                             fastpath=True, placement=self.mgr_locs)]
        return self._maybe_downcast(blocks, downcast)

    def take_nd(self, indexer, axis, new_mgr_locs=None, fill_tuple=None):
        """
        Take values according to indexer and return them as a block.bb

        """
        if fill_tuple is None:
            fill_value = self.fill_value
            new_values = com.take_nd(self.get_values(), indexer, axis=axis,
                                     allow_fill=False)
        else:
            fill_value = fill_tuple[0]
            new_values = com.take_nd(self.get_values(), indexer, axis=axis,
                                     allow_fill=True, fill_value=fill_value)

        if new_mgr_locs is None:
            if axis == 0:
                slc = lib.indexer_as_slice(indexer)
                if slc is not None:
                    new_mgr_locs = self.mgr_locs[slc]
                else:
                    new_mgr_locs = self.mgr_locs[indexer]
            else:
                new_mgr_locs = self.mgr_locs

        if new_values.dtype != self.dtype:
            return make_block(new_values, new_mgr_locs)
        else:
            return self.make_block_same_class(new_values, new_mgr_locs)

    def get_values(self, dtype=None):
        return self.values

    def diff(self, n):
        """ return block for the diff of the values """
        new_values = com.diff(self.values, n, axis=1)
        return [make_block(values=new_values,
                           ndim=self.ndim, fastpath=True,
                           placement=self.mgr_locs)]

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
                           placement=self.mgr_locs)]

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
                           fastpath=True, placement=self.mgr_locs)]

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
                              ndim=self.ndim, placement=self.mgr_locs)

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
                                                placement=self.mgr_locs[m]))

        return result_blocks

    def equals(self, other):
        if self.dtype != other.dtype or self.shape != other.shape: return False
        return np.array_equal(self.values, other.values)


class NumericBlock(Block):
    __slots__ = ()
    is_numeric = True
    _can_hold_na = True


class FloatOrComplexBlock(NumericBlock):
    __slots__ = ()

    def equals(self, other):
        if self.dtype != other.dtype or self.shape != other.shape: return False
        left, right = self.values, other.values
        return ((left == right) | (np.isnan(left) & np.isnan(right))).all()


class FloatBlock(FloatOrComplexBlock):
    __slots__ = ()
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
            imask = (~mask).ravel()
            values.flat[imask] = np.array(
                [float_format % val for val in values.ravel()[imask]])
        return values.tolist()

    def should_store(self, value):
        # when inserting a column should not coerce integers to floats
        # unnecessarily
        return (issubclass(value.dtype.type, np.floating) and
                value.dtype == self.dtype)


class ComplexBlock(FloatOrComplexBlock):
    __slots__ = ()
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
    __slots__ = ()
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
    __slots__ = ()
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
        imask = (~mask).ravel()
        rvalues.flat[imask] = np.array([lib.repr_timedelta64(val)
                                        for val in values.ravel()[imask]],
                                       dtype=object)
        return rvalues.tolist()


class BoolBlock(NumericBlock):
    __slots__ = ()
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
            return self
        return super(BoolBlock, self).replace(to_replace, value,
                                              inplace=inplace, filter=filter,
                                              regex=regex)


class ObjectBlock(Block):
    __slots__ = ()
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

            for i, rl in enumerate(self.mgr_locs):
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
                                     ndim=self.ndim, placement=self.mgr_locs))

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
            filt = self.mgr_locs.isin(filter).nonzero()[0]

        new_values[filt] = f(new_values[filt])

        return [self if inplace else
                make_block(new_values,
                           fastpath=True, placement=self.mgr_locs)]


class DatetimeBlock(Block):
    __slots__ = ()
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
                           fastpath=True, placement=self.mgr_locs)]

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
        imask = (~mask).ravel()

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
    __slots__ = ()
    is_sparse = True
    is_numeric = True
    _can_hold_na = True
    _can_consolidate = False
    _verify_integrity = False
    _ftype = 'sparse'

    def __init__(self, values, placement,
                 ndim=None, fastpath=False,):

        # kludgetastic
        if ndim is None:
            if len(placement) != 1:
                ndim = 1
            else:
                ndim = 2
        self.ndim = ndim

        self.mgr_locs = placement

        if not isinstance(values, SparseArray):
            raise TypeError("values must be SparseArray")

        self.values = values

    @property
    def shape(self):
        return (len(self.mgr_locs), self.sp_index.length)

    @property
    def itemsize(self):
        return self.dtype.itemsize

    @property
    def fill_value(self):
        #return np.nan
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
                                  fill_value=self.values.fill_value,
                                  copy=False)

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
        return self.make_block_same_class(values=self.values,
                                          sparse_index=self.sp_index,
                                          kind=self.kind, copy=deep,
                                          placement=self.mgr_locs)

    def make_block_same_class(self, values, placement,
                              sparse_index=None, kind=None, dtype=None,
                              fill_value=None, copy=False, fastpath=True):
        """ return a new block """
        if dtype is None:
            dtype = self.dtype
        if fill_value is None:
            fill_value = self.values.fill_value

        # if not isinstance(values, SparseArray) and values.ndim != self.ndim:
        #     raise ValueError("ndim mismatch")

        if values.ndim == 2:
            nitems = values.shape[0]

            if nitems == 0:
                # kludgy, but SparseBlocks cannot handle slices, where the
                # output is 0-item, so let's convert it to a dense block: it
                # won't take space since there's 0 items, plus it will preserve
                # the dtype.
                return make_block(np.empty(values.shape, dtype=dtype),
                                  placement, fastpath=True,)
            elif nitems > 1:
                raise ValueError("Only 1-item 2d sparse blocks are supported")
            else:
                values = values.reshape(values.shape[1])

        new_values = SparseArray(values, sparse_index=sparse_index,
                                 kind=kind or self.kind, dtype=dtype,
                                 fill_value=fill_value, copy=copy)
        return make_block(new_values, ndim=self.ndim,
                          fastpath=fastpath, placement=placement)

    def interpolate(self, method='pad', axis=0, inplace=False,
                    limit=None, fill_value=None, **kwargs):

        values = com.interpolate_2d(
            self.values.to_dense(), method, axis, limit, fill_value)
        return self.make_block_same_class(values=values,
                                          placement=self.mgr_locs)

    def fillna(self, value, limit=None, inplace=False, downcast=None):
        # we may need to upcast our fill to match our dtype
        if limit is not None:
            raise NotImplementedError
        if issubclass(self.dtype.type, np.floating):
            value = float(value)
        values = self.values if inplace else self.values.copy()
        return [self.make_block_same_class(values=values.get_values(value),
                                           fill_value=value,
                                           placement=self.mgr_locs)]

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
        return [self.make_block_same_class(new_values, placement=self.mgr_locs)]

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
        return self.make_block_same_class(self.values.take(indexer),
                                          fill_value=fill_value,
                                          placement=self.mgr_locs)

    def sparse_reindex(self, new_index):
        """ sparse reindex and return a new block
            current reindex only works for float64 dtype! """
        values = self.values
        values = values.sp_index.to_int_index().reindex(
            values.sp_values.astype('float64'), values.fill_value, new_index)
        return self.make_block_same_class(values, sparse_index=new_index,
                               placement=self.mgr_locs)

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
                 '_is_consolidated', '_blknos', '_blklocs']

    def __init__(self, blocks, axes, do_integrity_check=True, fastpath=True):
        self.axes = [_ensure_index(ax) for ax in axes]
        self.blocks = tuple(blocks)

        for block in blocks:
            if block.is_sparse:
                if len(block.mgr_locs) != 1:
                    raise AssertionError("Sparse block refers to multiple items")
            else:
                if self.ndim != block.ndim:
                    raise AssertionError(('Number of Block dimensions (%d) must '
                                          'equal number of axes (%d)')
                                         % (block.ndim, self.ndim))

        if do_integrity_check:
            self._verify_integrity()

        self._consolidate_check()

        self._rebuild_blknos_and_blklocs()

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

    def rename_axis(self, mapper, axis, copy=True):
        """
        Rename one of axes.

        Parameters
        ----------
        mapper : unary callable
        axis : int
        copy : boolean, default True

        """
        obj = self.copy(deep=copy)
        obj.set_axis(axis, _transform_index(self.axes[axis], mapper))
        return obj

    def add_prefix(self, prefix):
        f = (str(prefix) + '%s').__mod__
        return self.rename_axis(f, axis=0)

    def add_suffix(self, suffix):
        f = ('%s' + str(suffix)).__mod__
        return self.rename_axis(f, axis=0)

    @property
    def _is_single_block(self):
        if self.ndim == 1:
            return True

        if len(self.blocks) != 1:
            return False

        blk = self.blocks[0]
        return (blk.mgr_locs.is_slice_like and
                blk.mgr_locs.as_slice == slice(0, len(self), 1))

    def _rebuild_blknos_and_blklocs(self):
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
            raise AssertionError("Gaps in blk ref_locs")

        self._blknos = new_blknos
        self._blklocs = new_blklocs

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
        dtypes = np.array([blk.dtype for blk in self.blocks])
        return com.take_1d(dtypes, self._blknos, allow_fill=False)

    def get_ftypes(self):
        ftypes = np.array([blk.ftype for blk in self.blocks])
        return com.take_1d(ftypes, self._blknos, allow_fill=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [self.items[b.mgr_locs.indexer] for b in self.blocks]
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
        self.blocks = tuple(blocks)

        self._post_setstate()

    def _post_setstate(self):
        self._is_consolidated = False
        self._known_consolidated = False
        self._rebuild_blknos_and_blklocs()

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
        tot_items = sum(len(x.mgr_locs) for x in self.blocks)
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

        # filter kwarg is used in replace-* family of methods
        if filter is not None:
            filter_locs = set(self.items.get_indexer_for(filter))
            if len(filter_locs) == len(self.items):
                # All items are included, as if there were no filtering
                filter = None
            else:
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
                if not b.mgr_locs.isin(filter_locs).any():
                    result_blocks.append(b)
                    continue

            if aligned_args:
                b_items = self.items[b.mgr_locs.indexer]

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
                        m = masks[i][b.mgr_locs.indexer]
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

    @property
    def is_view(self):
        """ return a boolean if we are a single block and are a view """
        if len(self.blocks) == 1:
            return self.blocks[0].values.base is not None
        return False

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

        # FIXME: optimization potential
        indexer = np.sort(np.concatenate([b.mgr_locs.as_array for b in blocks]))
        inv_indexer = lib.get_reverse_indexer(indexer, self.shape[0])
        new_items = self.items.take(indexer)

        new_blocks = []
        for b in blocks:
            b = b.copy(deep=copy)
            b.mgr_locs = com.take_1d(inv_indexer, b.mgr_locs.as_array, axis=0,
                                     allow_fill=False)
            new_blocks.append(b)

        new_axes = list(self.axes)
        new_axes[0] = new_items
        return self.__class__(new_blocks, new_axes, do_integrity_check=False)

    def get_slice(self, slobj, axis=0):
        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(slobj)
        else:
            slicer = [slice(None)] * (axis + 1)
            slicer[axis] = slobj
            slicer = tuple(slicer)
            new_blocks = [blk.getitem_block(slicer) for blk in self.blocks]

        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis][slobj]

        bm = self.__class__(new_blocks, new_axes, do_integrity_check=False,
                            fastpath=True)
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

        if self._is_single_block:
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

        if result.shape[0] == 0:
            # Workaround for numpy 1.7 bug:
            #
            #     >>> a = np.empty((0,10))
            #     >>> a[slice(0,0)]
            #     array([], shape=(0, 10), dtype=float64)
            #     >>> a[[]]
            #     Traceback (most recent call last):
            #       File "<stdin>", line 1, in <module>
            #     IndexError: index 0 is out of bounds for axis 0 with size 0
            return result

        itemmask = np.zeros(self.shape[0])

        for blk in self.blocks:
            rl = blk.mgr_locs
            result[rl.indexer] = blk.get_values(dtype)
            itemmask[rl.indexer] = 1

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
                                  placement=blk.mgr_locs)
                new_blocks.append(newb)
        elif len(self.blocks) == 1:
            block = self.blocks[0]
            vals = block.values[slicer]
            if copy:
                vals = vals.copy()
            new_blocks = [make_block(values=vals, placement=block.mgr_locs,
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
            # result[blk.mgr_locs] = blk._slice((slice(None), loc))
            for i, rl in enumerate(blk.mgr_locs):
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
            self.blocks = tuple(_consolidate(self.blocks))

            self._is_consolidated = True
            self._known_consolidated = True
            self._rebuild_blknos_and_blklocs()

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
        return self.blocks[self._blknos[i]].iget(self._blklocs[i])

    def get_scalar(self, tup):
        """
        Retrieve single item
        """
        full_loc = list(ax.get_loc(x)
                        for ax, x in zip(self.axes, tup))
        blk = self.blocks[self._blknos[full_loc[0]]]
        full_loc[0] = self._blklocs[full_loc[0]]

        # FIXME: this may return non-upcasted types?
        return blk.values[tuple(full_loc)]

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

        for blkno, _ in _fast_count_smallints(self._blknos[affected_start:]):
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
        self.blocks = tuple(b for blkno, b in enumerate(self.blocks)
                            if not is_blk_deleted[blkno])
        self._shape = None
        self._rebuild_blknos_and_blklocs()

    def set(self, item, value, check=False):
        """
        Set new item in-place. Does not consolidate. Adds new Block if not
        contained in the current set of items
        if check, then validate that we are not setting the same data in-place
        """
        # FIXME: refactor, clearly separate broadcasting & zip-like assignment
        value_is_sparse = isinstance(value, SparseArray)

        if value_is_sparse:
            assert self.ndim == 2

            def value_getitem(placement):
                return value
        else:
            if value.ndim == self.ndim - 1:
                value = value.reshape((1,) + value.shape)

                def value_getitem(placement):
                    return value
            else:
                def value_getitem(placement):
                    return value[placement.indexer]
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

        blknos = self._blknos[loc]
        blklocs = self._blklocs[loc]

        unfit_mgr_locs = []
        unfit_val_locs = []
        removed_blknos = []
        for blkno, val_locs in _get_blkno_placements(blknos, len(self.blocks),
                                                     group=True):
            blk = self.blocks[blkno]
            blk_locs = blklocs[val_locs.indexer]
            if blk.should_store(value):
                blk.set(blk_locs, value_getitem(val_locs), check=check)
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
            new_blknos[~is_deleted] = np.arange(self.nblocks -
                                                len(removed_blknos))
            self._blknos = com.take_1d(new_blknos, self._blknos, axis=0,
                                       allow_fill=False)
            self.blocks = tuple(blk for i, blk in enumerate(self.blocks)
                                if i not in set(removed_blknos))

        if unfit_val_locs:
            unfit_mgr_locs = np.concatenate(unfit_mgr_locs)
            unfit_count = len(unfit_mgr_locs)

            new_blocks = []
            if value_is_sparse:
                # This code (ab-)uses the fact that sparse blocks contain only
                # one item.
                new_blocks.extend(
                    make_block(values=value.copy(), ndim=self.ndim,
                               placement=slice(mgr_loc, mgr_loc + 1))
                    for mgr_loc in unfit_mgr_locs)

                self._blknos[unfit_mgr_locs] = (np.arange(unfit_count) +
                                                len(self.blocks))
                self._blklocs[unfit_mgr_locs] = 0

            else:
                # unfit_val_locs contains BlockPlacement objects
                unfit_val_items = unfit_val_locs[0].append(unfit_val_locs[1:])

                new_blocks.append(
                    make_block(values=value_getitem(unfit_val_items),
                               ndim=self.ndim, placement=unfit_mgr_locs))

                self._blknos[unfit_mgr_locs] = len(self.blocks)
                self._blklocs[unfit_mgr_locs] = np.arange(unfit_count)

            self.blocks += tuple(new_blocks)

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

        block = make_block(values=value,
                           ndim=self.ndim,
                           placement=slice(loc, loc+1))

        for blkno, count in _fast_count_smallints(self._blknos[loc:]):
            blk = self.blocks[blkno]
            if count == len(blk.mgr_locs):
                blk.mgr_locs = blk.mgr_locs.add(1)
            else:
                new_mgr_locs = blk.mgr_locs.as_array.copy()
                new_mgr_locs[new_mgr_locs >= loc] += 1
                blk.mgr_locs = new_mgr_locs

        if loc == self._blklocs.shape[0]:
            # np.append is a lot faster (at least in numpy 1.7.1), let's use it
            # if we can.
            self._blklocs = np.append(self._blklocs, 0)
            self._blknos = np.append(self._blknos, len(self.blocks))
        else:
            self._blklocs = np.insert(self._blklocs, loc, 0)
            self._blknos = np.insert(self._blknos, loc, len(self.blocks))

        self.axes[0] = self.items.insert(loc, item)

        self.blocks += (block,)
        self._shape = None

        self._known_consolidated = False

        if len(self.blocks) > 100:
            self._consolidate_inplace()

    def reindex_axis(self, new_index, axis, method=None, limit=None,
                     fill_value=None, copy=True):
        """
        Conform block manager to new index.
        """
        new_index = _ensure_index(new_index)
        new_index, indexer = self.axes[axis].reindex(
            new_index, method=method, limit=limit, copy_if_needed=True)

        return self.reindex_indexer(new_index, indexer, axis=axis,
                                    fill_value=fill_value, copy=copy)

    def reindex_indexer(self, new_axis, indexer, axis, fill_value=None,
                        allow_dups=False, copy=True):
        """
        Parameters
        ----------
        new_axis : Index
        indexer : ndarray of int64 or None
        axis : int
        fill_value : object
        allow_dups : bool

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

        # trying to reindex on an axis with duplicates
        if (not allow_dups and not self.axes[axis].is_unique
            and len(indexer)):
            raise ValueError("cannot reindex from a duplicate axis")

        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(
                indexer, fill_tuple=(fill_value,))
        else:
            new_blocks = [blk.take_nd(indexer, axis=axis,
                                      fill_tuple=(fill_value if fill_value is not None else
                                                  blk.fill_value,))
                          for blk in self.blocks]

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return self.__class__(new_blocks, new_axes)

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
            slice_or_indexer, self.shape[0], allow_fill=allow_fill)

        if self._is_single_block:
            blk = self.blocks[0]

            if sl_type in ('slice', 'mask'):
                return [blk.getitem_block(slobj,
                                          new_mgr_locs=slice(0, sllen))]
            elif not allow_fill or self.ndim == 1:
                if allow_fill and fill_tuple[0] is None:
                    _, fill_value = com._maybe_promote(blk.dtype)
                    fill_tuple = (fill_value,)

                return [blk.take_nd(slobj, axis=0,
                                    new_mgr_locs=slice(0, sllen),
                                    fill_tuple=fill_tuple)]

        if sl_type in ('slice', 'mask'):
            blknos = self._blknos[slobj]
            blklocs = self._blklocs[slobj]
        else:
            blknos = com.take_1d(self._blknos, slobj, fill_value=-1,
                                 allow_fill=allow_fill)
            blklocs = com.take_1d(self._blklocs, slobj, fill_value=-1,
                                  allow_fill=allow_fill)

        # When filling blknos, make sure blknos is updated before appending to
        # blocks list, that way new blkno is exactly len(blocks).
        #
        # FIXME: mgr_groupby_blknos must return mgr_locs in ascending order,
        # pytables serialization will break otherwise.
        blocks = []
        for blkno, mgr_locs in _get_blkno_placements(blknos, len(self.blocks),
                                                     group=True):
            if blkno == -1:
                # If we've got here, fill_tuple was not None.
                fill_value = fill_tuple[0]

                blocks.append(self._make_na_block(
                    placement=mgr_locs, fill_value=fill_value))
            else:
                blk = self.blocks[blkno]

                # Otherwise, slicing along items axis is necessary.
                if blk.is_sparse:
                    # A sparse block, it's easy, because there's only one item
                    # and each mgr loc is a copy of that single item.
                    for mgr_loc in mgr_locs:
                        newblk = blk.copy(deep=True)
                        newblk.mgr_locs = slice(mgr_loc, mgr_loc + 1)
                        blocks.append(newblk)

                else:
                    blocks.append(blk.take_nd(
                        blklocs[mgr_locs.indexer], axis=0,
                        new_mgr_locs=mgr_locs, fill_tuple=None))

        return blocks

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

        new_blocks = [blk.copy(deep=False)
                      for blk in self.blocks]

        offset = self.shape[0]
        for blk in other.blocks:
            blk = blk.copy(deep=False)
            blk.mgr_locs = blk.mgr_locs.add(offset)
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


class SingleBlockManager(BlockManager):
    """ manage a single block with """

    ndim = 1
    _is_consolidated = True
    _known_consolidated = True
    __slots__ = ()

    def __init__(self, block, axis, do_integrity_check=False, fastpath=False):

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
            block = make_block(block,
                               placement=slice(0, len(axis)),
                               ndim=1, fastpath=True)

        self.blocks = [block]

    def _post_setstate(self):
        pass

    @property
    def _block(self):
        return self.blocks[0]

    @property
    def _values(self):
        return self._block.values

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
            make_block = self._block.make_block_same_class

        block = make_block(new_values, copy=copy,
                           placement=slice(0, len(new_axis)))

        mgr = SingleBlockManager(block, new_axis)
        mgr._consolidate_inplace()
        return mgr

    def get_slice(self, slobj, axis=0):
        if axis >= self.ndim:
            raise IndexError("Requested axis not found in manager")

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
        return np.array([self._block.dtype])

    def get_ftypes(self):
        return np.array([self._block.ftype])

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
        loc = self.items.get_loc(item)
        self._block.delete(loc)
        self.axes[0] = self.axes[0].delete(loc)

    def fast_xs(self, loc):
        """
        fast path for getting a cross-section
        return a view of the data
        """
        return self._block.values[loc]


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
                                 placement=slice(0, len(axes[0])))]

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

        # FIXME: optimization potential in case all mgrs contain slices and
        # combination of those slices is a slice, too.
        new_mgr_locs = np.concatenate([b.mgr_locs.as_array for b in blocks])
        new_values = _vstack([b.values for b in blocks], dtype)

        argsort = np.argsort(new_mgr_locs)
        new_values = new_values[argsort]
        new_mgr_locs = new_mgr_locs[argsort]

        return make_block(new_values,
                          fastpath=True, placement=new_mgr_locs)

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


def _get_blkno_placements(blknos, blk_count, group=True):
    """

    Parameters
    ----------
    blknos : array of int64
    blk_count : int
    group : bool

    Returns
    -------
    iterator
        yield (BlockPlacement, blkno)

    """

    blknos = com._ensure_int64(blknos)

    # FIXME: blk_count is unused, but it may avoid the use of dicts in cython
    for blkno, indexer in lib.get_blkno_indexers(blknos, group):
        yield blkno, BlockPlacement(indexer)


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
    concat_plan = combine_concat_plans([get_mgr_concatenation_plan(mgr, indexers)
                                        for mgr, indexers in mgrs_indexers],
                                       concat_axis)

    blocks = [make_block(concatenate_join_units(join_units, concat_axis,
                                                copy=copy),
                         placement=placement)
              for placement, join_units in concat_plan]

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

    if len(join_units) == 1:
        blk = join_units[0].block
        if blk is None:
            return np.float64, np.nan
        else:
            return blk.dtype, None

    has_none_blocks = False
    dtypes = [None] * len(join_units)

    for i, unit in enumerate(join_units):
        if unit.block is None:
            has_none_blocks = True
        else:
            dtypes[i] = unit.dtype

    if not has_none_blocks and len(set(dtypes)) == 1:
        # Unanimous decision, nothing to upcast.
        return dtypes[0], None

    # dtypes = set()
    upcast_classes = set()
    null_upcast_classes = set()
    for dtype, unit in zip(dtypes, join_units):
        if dtype is None:
            continue

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


def concatenate_join_units(join_units, concat_axis, copy):
    """
    Concatenate values from several join units along selected axis.
    """
    if concat_axis == 0 and len(join_units) > 1:
        # Concatenating join units along ax0 is handled in _merge_blocks.
        raise AssertionError("Concatenating join units along axis0")

    empty_dtype, upcasted_na = get_empty_dtype_and_na(join_units)

    to_concat = [ju.get_reindexed_values(empty_dtype=empty_dtype,
                                         upcasted_na=upcasted_na)
                 for ju in join_units]

    if len(to_concat) == 1:
        # Only one block, nothing to concatenate.
        concat_values = to_concat[0]
        if copy and concat_values.base is not None:
            concat_values = concat_values.copy()
    else:
        concat_values = com._concat_compat(to_concat, axis=concat_axis)

    # FIXME: optimization potential: if len(join_units) == 1, single join unit
    # is densified and sparsified back.
    if any(unit.is_sparse for unit in join_units):
        # If one of the units was sparse, concat_values are 2d and there's only
        # one item.
        return SparseArray(concat_values[0])
    else:
        return concat_values


def get_mgr_concatenation_plan(mgr, indexers):
    """
    Construct concatenation plan for given block manager and indexers.

    Parameters
    ----------
    mgr : BlockManager
    indexers : dict of {axis: indexer}

    Returns
    -------
    plan : list of (BlockPlacement, JoinUnit) tuples

    """
    # Calculate post-reindex shape , save for item axis which will be separate
    # for each block anyway.
    mgr_shape = list(mgr.shape)
    for ax, indexer in indexers.items():
        mgr_shape[ax] = len(indexer)
    mgr_shape = tuple(mgr_shape)

    if 0 in indexers:
        ax0_indexer = indexers.pop(0)
        blknos = com.take_1d(mgr._blknos, ax0_indexer, fill_value=-1)
        blklocs = com.take_1d(mgr._blklocs, ax0_indexer, fill_value=-1)
    else:

        if mgr._is_single_block:
            blk = mgr.blocks[0]
            return [(blk.mgr_locs, JoinUnit(blk, mgr_shape, indexers))]

        ax0_indexer = None
        blknos = mgr._blknos
        blklocs = mgr._blklocs

    plan = []
    for blkno, placements in _get_blkno_placements(blknos, len(mgr.blocks),
                                                   group=False):
        assert placements.is_slice_like

        join_unit_indexers = indexers.copy()

        shape = list(mgr_shape)
        shape[0] = len(placements)
        shape = tuple(shape)

        if blkno == -1:
            unit = JoinUnit(None, shape)
        else:
            blk = mgr.blocks[blkno]
            ax0_blk_indexer = blklocs[placements.indexer]

            unit_no_ax0_reindexing = (
                len(placements) == len(blk.mgr_locs) and
                # Fastpath detection of join unit not needing to reindex its
                # block: no ax0 reindexing took place and block placement was
                # sequential before.
                ((ax0_indexer is None
                  and blk.mgr_locs.is_slice_like
                  and blk.mgr_locs.as_slice.step == 1) or
                 # Slow-ish detection: all indexer locs are sequential (and
                 # length match is checked above).
                 (np.diff(ax0_blk_indexer) == 1).all()))

            # Omit indexer if no item reindexing is required.
            if unit_no_ax0_reindexing:
                join_unit_indexers.pop(0, None)
            else:
                join_unit_indexers[0] = ax0_blk_indexer

            unit = JoinUnit(blk, shape, join_unit_indexers)

        plan.append((placements, unit))

    return plan


def combine_concat_plans(plans, concat_axis):
    """
    Combine multiple concatenation plans into one.

    existing_plan is updated in-place.
    """
    if len(plans) == 1:
        for p in plans[0]:
            yield p[0], [p[1]]

    elif concat_axis == 0:
        offset = 0
        for plan in plans:
            last_plc = None

            for plc, unit in plan:
                yield plc.add(offset), [unit]
                last_plc = plc

            if last_plc is not None:
                offset += last_plc.as_slice.stop

    else:
        num_ended = [0]
        def _next_or_none(seq):
            retval = next(seq, None)
            if retval is None:
                num_ended[0] += 1
            return retval

        plans = list(map(iter, plans))
        next_items = list(map(_next_or_none, plans))

        while num_ended[0] != len(next_items):
            if num_ended[0] > 0:
                raise ValueError("Plan shapes are not aligned")

            placements, units = zip(*next_items)

            lengths = list(map(len, placements))
            min_len, max_len = min(lengths), max(lengths)

            if min_len == max_len:
                yield placements[0], units
                next_items[:] = map(_next_or_none, plans)
            else:
                yielded_placement = None
                yielded_units = [None] * len(next_items)
                for i, (plc, unit) in enumerate(next_items):
                    yielded_units[i] = unit
                    if len(plc) > min_len:
                        # trim_join_unit updates unit in place, so only
                        # placement needs to be sliced to skip min_len.
                        next_items[i] = (plc[min_len:],
                                         trim_join_unit(unit, min_len))
                    else:
                        yielded_placement = plc
                        next_items[i] = _next_or_none(plans[i])

                yield yielded_placement, yielded_units


def trim_join_unit(join_unit, length):
    """
    Reduce join_unit's shape along item axis to length.

    Extra items that didn't fit are returned as a separate block.
    """

    if 0 not in join_unit.indexers:
        extra_indexers = join_unit.indexers

        if join_unit.block is None:
            extra_block = None
        else:
            extra_block = join_unit.block.getitem_block(slice(length, None))
            join_unit.block = join_unit.block.getitem_block(slice(length))
    else:
        extra_block = join_unit.block

        extra_indexers = copy.copy(join_unit.indexers)
        extra_indexers[0] = extra_indexers[0][length:]
        join_unit.indexers[0] = join_unit.indexers[0][:length]

    extra_shape = (join_unit.shape[0] - length,) + join_unit.shape[1:]
    join_unit.shape = (length,) + join_unit.shape[1:]

    return JoinUnit(block=extra_block, indexers=extra_indexers,
                    shape=extra_shape)


class JoinUnit(object):
    def __init__(self, block, shape, indexers={}):
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
        if self.block is None:
            return True

        if not self.block._can_hold_na:
            return False

        # Usually it's enough to check but a small fraction of values to see if
        # a block is NOT null, chunks should help in such cases.  1000 value
        # was chosen rather arbitrarily.
        values_flat = self.block.values.ravel()
        total_len = values_flat.shape[0]
        chunk_len = max(total_len // 40, 1000)
        for i in range(0, total_len, chunk_len):
            if not isnull(values_flat[i: i + chunk_len]).all():
                return False

        return True

    @cache_readonly
    def is_sparse(self):
        return self.block is not None and self.block.is_sparse

    def get_reindexed_values(self, empty_dtype, upcasted_na):
        if upcasted_na is None:
            # No upcasting is necessary
            fill_value = self.block.fill_value
            values = self.block.get_values()
        else:
            fill_value = upcasted_na

            if self.is_null:
                missing_arr = np.empty(self.shape, dtype=empty_dtype)
                if np.prod(self.shape):
                    # NumPy 1.6 workaround: this statement gets strange if all
                    # blocks are of same dtype and some of them are empty:
                    # empty one are considered "null" so they must be filled,
                    # but no dtype upcasting happens and the dtype may not
                    # allow NaNs.
                    #
                    # In general, no one should get hurt when one tries to put
                    # incorrect values into empty array, but numpy 1.6 is
                    # strict about that.
                    missing_arr.fill(fill_value)
                return missing_arr

            if self.block.is_bool:
                # External code requested filling/upcasting, bool values must
                # be upcasted to object to avoid being upcasted to numeric.
                values = self.block.astype(np.object_).values
            else:
                # No dtype upcasting is done here, it will be performed during
                # concatenation itself.
                values = self.block.get_values()

        if not self.indexers:
            # If there's no indexing to be done, we want to signal outside
            # code that this array must be copied explicitly.  This is done
            # by returning a view and checking `retval.base`.
            return values.view()
        else:
            for ax, indexer in self.indexers.items():
                values = com.take_nd(values, indexer, axis=ax,
                                     fill_value=fill_value)

            return values


def _fast_count_smallints(arr):
    """Faster version of set(arr) for sequences of small numbers."""
    if len(arr) == 0:
        # Handle empty arr case separately: numpy 1.6 chokes on that.
        return np.empty((0, 2), dtype=arr.dtype)
    else:
        counts = np.bincount(arr.astype(np.int_))
        nz = counts.nonzero()[0]
        return np.c_[nz, counts[nz]]


def _preprocess_slice_or_indexer(slice_or_indexer, length, allow_fill):
    if isinstance(slice_or_indexer, slice):
        return 'slice', slice_or_indexer, lib.slice_len(slice_or_indexer,
                                                        length)
    elif (isinstance(slice_or_indexer, np.ndarray) and
          slice_or_indexer.dtype == np.bool_):
        return 'mask', slice_or_indexer, slice_or_indexer.sum()
    else:
        indexer = np.asanyarray(slice_or_indexer, dtype=np.int64)
        if not allow_fill:
            indexer = _maybe_convert_indices(indexer, length)
        return 'fancy', indexer, len(indexer)
