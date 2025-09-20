import numpy as np
import operator
from collections import namedtuple

from numba.core import types, utils
from numba.core.typing.templates import (AttributeTemplate, AbstractTemplate,
                                         infer, infer_global, infer_getattr,
                                         signature, bound_function)
# import time side effect: array operations requires typing support of sequence
# defined in collections: e.g. array.shape[i]
from numba.core.typing import collections
from numba.core.errors import (TypingError, RequireLiteralValue, NumbaTypeError,
                               NumbaNotImplementedError, NumbaAssertionError,
                               NumbaKeyError, NumbaIndexError, NumbaValueError)
from numba.core.cgutils import is_nonelike

numpy_version = tuple(map(int, np.__version__.split('.')[:2]))


Indexing = namedtuple("Indexing", ("index", "result", "advanced"))


def get_array_index_type(ary, idx):
    """
    Returns None or a tuple-3 for the types of the input array, index, and
    resulting type of ``array[index]``.

    Note: This is shared logic for ndarray getitem and setitem.
    """
    if not isinstance(ary, types.Buffer):
        return

    ndim = ary.ndim

    left_indices = []
    right_indices = []
    ellipsis_met = False
    advanced = False
    num_newaxis = 0

    if not isinstance(idx, types.BaseTuple):
        idx = [idx]

    # Here, a subspace is considered as a contiguous group of advanced indices.
    # num_subspaces keeps track of the number of such
    # contiguous groups.
    in_subspace = False
    num_subspaces = 0
    array_indices = 0

    # Walk indices
    for ty in idx:
        if ty is types.ellipsis:
            if ellipsis_met:
                raise NumbaTypeError(
                    "Only one ellipsis allowed in array indices "
                    "(got %s)" % (idx,))
            ellipsis_met = True
            in_subspace = False
        elif isinstance(ty, types.SliceType):
            # If we encounter a non-advanced index while in a
            # subspace then that subspace ends.
            in_subspace = False
        # In advanced indexing, any index broadcastable to an
        # array is considered an advanced index. Hence all the
        # branches below are considered as advanced indices.
        elif isinstance(ty, types.Integer):
            # Normalize integer index
            ty = types.intp if ty.signed else types.uintp
            # Integer indexing removes the given dimension
            ndim -= 1
            # If we're within a subspace/contiguous group of
            # advanced indices then no action is necessary
            # since we've already counted that subspace once.
            if not in_subspace:
                # If we're not within a subspace and we encounter
                # this branch then we have a new subspace/group.
                num_subspaces += 1
                in_subspace = True
        elif (isinstance(ty, types.Array) and ty.ndim == 0
              and isinstance(ty.dtype, types.Integer)):
            # 0-d array used as integer index
            ndim -= 1
            if not in_subspace:
                num_subspaces += 1
                in_subspace = True
        elif (isinstance(ty, types.Array)
              and isinstance(ty.dtype, (types.Integer, types.Boolean))):
            if ty.ndim > 1:
                # Advanced indexing limitation # 1
                raise NumbaTypeError(
                    "Multi-dimensional indices are not supported.")
            array_indices += 1
            # The condition for activating advanced indexing is simply
            # having at least one array with size > 1.
            advanced = True
            if not in_subspace:
                num_subspaces += 1
                in_subspace = True
        elif (is_nonelike(ty)):
            ndim += 1
            num_newaxis += 1
        else:
            raise NumbaTypeError("Unsupported array index type %s in %s"
                                 % (ty, idx))
        (right_indices if ellipsis_met else left_indices).append(ty)

    if advanced:
        if array_indices > 1:
            # Advanced indexing limitation # 2
            msg = "Using more than one non-scalar array index is unsupported."
            raise NumbaTypeError(msg)

        if num_subspaces > 1:
            # Advanced indexing limitation # 3
            msg = ("Using more than one indexing subspace is unsupported."
                   " An indexing subspace is a group of one or more"
                   " consecutive indices comprising integer or array types.")
            raise NumbaTypeError(msg)

    # Only Numpy arrays support advanced indexing
    if advanced and not isinstance(ary, types.Array):
        return

    # Check indices and result dimensionality
    all_indices = left_indices + right_indices
    if ellipsis_met:
        assert right_indices[0] is types.ellipsis
        del right_indices[0]

    n_indices = len(all_indices) - ellipsis_met - num_newaxis
    if n_indices > ary.ndim:
        raise NumbaTypeError("cannot index %s with %d indices: %s"
                             % (ary, n_indices, idx))
    if n_indices == ary.ndim and ndim == 0 and not ellipsis_met:
        # Full integer indexing => scalar result
        # (note if ellipsis is present, a 0-d view is returned instead)
        res = ary.dtype

    elif advanced:
        # Result is a copy
        res = ary.copy(ndim=ndim, layout='C', readonly=False)

    else:
        # Result is a view
        if ary.slice_is_copy:
            # Avoid view semantics when the original type creates a copy
            # when slicing.
            return

        # Infer layout
        layout = ary.layout

        def keeps_contiguity(ty, is_innermost):
            # A slice can only keep an array contiguous if it is the
            # innermost index and it is not strided
            return (ty is types.ellipsis or isinstance(ty, types.Integer)
                    or (is_innermost and isinstance(ty, types.SliceType)
                        and not ty.has_step))

        def check_contiguity(outer_indices):
            """
            Whether indexing with the given indices (from outer to inner in
            physical layout order) can keep an array contiguous.
            """
            for ty in outer_indices[:-1]:
                if not keeps_contiguity(ty, False):
                    return False
            if outer_indices and not keeps_contiguity(outer_indices[-1], True):
                return False
            return True

        if layout == 'C':
            # Integer indexing on the left keeps the array C-contiguous
            if n_indices == ary.ndim:
                # If all indices are there, ellipsis's place is indifferent
                left_indices = left_indices + right_indices
                right_indices = []
            if right_indices:
                layout = 'A'
            elif not check_contiguity(left_indices):
                layout = 'A'
        elif layout == 'F':
            # Integer indexing on the right keeps the array F-contiguous
            if n_indices == ary.ndim:
                # If all indices are there, ellipsis's place is indifferent
                right_indices = left_indices + right_indices
                left_indices = []
            if left_indices:
                layout = 'A'
            elif not check_contiguity(right_indices[::-1]):
                layout = 'A'

        if ndim == 0:
            # Implicitly convert to a scalar if the output ndim==0
            res = ary.dtype
        else:
            res = ary.copy(ndim=ndim, layout=layout)

    # Re-wrap indices
    if isinstance(idx, types.BaseTuple):
        idx = types.BaseTuple.from_types(all_indices)
    else:
        idx, = all_indices

    return Indexing(idx, res, advanced)


@infer_global(operator.getitem)
class GetItemBuffer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        [ary, idx] = args
        out = get_array_index_type(ary, idx)
        if out is not None:
            return signature(out.result, ary, out.index)


@infer_global(operator.setitem)
class SetItemBuffer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        ary, idx, val = args
        if not isinstance(ary, types.Buffer):
            return
        if not ary.mutable:
            msg = f"Cannot modify readonly array of type: {ary}"
            raise NumbaTypeError(msg)
        out = get_array_index_type(ary, idx)
        if out is None:
            return

        idx = out.index
        res = out.result  # res is the result type of the access ary[idx]
        if isinstance(res, types.Array):
            # Indexing produces an array
            if isinstance(val, types.Array):
                if not self.context.can_convert(val.dtype, res.dtype):
                    # DType conversion not possible
                    return
                else:
                    res = val
            elif isinstance(val, types.Sequence):
                if (res.ndim == 1 and
                    self.context.can_convert(val.dtype, res.dtype)):
                    # Allow assignment of sequence to 1d array
                    res = val
                else:
                    # NOTE: sequence-to-array broadcasting is unsupported
                    return
            else:
                # Allow scalar broadcasting
                if self.context.can_convert(val, res.dtype):
                    res = res.dtype
                else:
                    # Incompatible scalar type
                    return
        elif not isinstance(val, types.Array):
            # Single item assignment
            if not self.context.can_convert(val, res):
                # if the array dtype is not yet defined
                if not res.is_precise():
                    # set the array type to use the dtype of value (RHS)
                    newary = ary.copy(dtype=val)
                    return signature(types.none, newary, idx, res)
                else:
                    return
            res = val
        elif (isinstance(val, types.Array) and val.ndim == 0
              and self.context.can_convert(val.dtype, res)):
            # val is an array(T, 0d, O), where T is the type of res, O is order
            res = val
        else:
            return
        return signature(types.none, ary, idx, res)


def normalize_shape(shape):
    if isinstance(shape, types.UniTuple):
        if isinstance(shape.dtype, types.Integer):
            dimtype = types.intp if shape.dtype.signed else types.uintp
            return types.UniTuple(dimtype, len(shape))

    elif isinstance(shape, types.Tuple) and shape.count == 0:
        # Force (0 x intp) for consistency with other shapes
        return types.UniTuple(types.intp, 0)


@infer_getattr
class ArrayAttribute(AttributeTemplate):
    key = types.Array

    def resolve_dtype(self, ary):
        return types.DType(ary.dtype)

    def resolve_nbytes(self, ary):
        return types.intp

    def resolve_itemsize(self, ary):
        return types.intp

    def resolve_shape(self, ary):
        return types.UniTuple(types.intp, ary.ndim)

    def resolve_strides(self, ary):
        return types.UniTuple(types.intp, ary.ndim)

    def resolve_ndim(self, ary):
        return types.intp

    def resolve_size(self, ary):
        return types.intp

    def resolve_flat(self, ary):
        return types.NumpyFlatType(ary)

    def resolve_ctypes(self, ary):
        return types.ArrayCTypes(ary)

    def resolve_flags(self, ary):
        return types.ArrayFlags(ary)

    def resolve_T(self, ary):
        if ary.ndim <= 1:
            retty = ary
        else:
            layout = {"C": "F", "F": "C"}.get(ary.layout, "A")
            retty = ary.copy(layout=layout)
        return retty

    def resolve_real(self, ary):
        return self._resolve_real_imag(ary, attr='real')

    def resolve_imag(self, ary):
        return self._resolve_real_imag(ary, attr='imag')

    def _resolve_real_imag(self, ary, attr):
        if ary.dtype in types.complex_domain:
            return ary.copy(dtype=ary.dtype.underlying_float, layout='A')
        elif ary.dtype in types.number_domain:
            res = ary.copy(dtype=ary.dtype)
            if attr == 'imag':
                res = res.copy(readonly=True)
            return res
        else:
            msg = "cannot access .{} of array of {}"
            raise TypingError(msg.format(attr, ary.dtype))

    @bound_function("array.transpose")
    def resolve_transpose(self, ary, args, kws):
        def sentry_shape_scalar(ty):
            if ty in types.number_domain:
                # Guard against non integer type
                if not isinstance(ty, types.Integer):
                    msg = "transpose() arg cannot be {0}".format(ty)
                    raise TypingError(msg)
                return True
            else:
                return False

        assert not kws
        if len(args) == 0:
            return signature(self.resolve_T(ary))

        if len(args) == 1:
            shape, = args

            if sentry_shape_scalar(shape):
                assert ary.ndim == 1
                return signature(ary, *args)

            if isinstance(shape, types.NoneType):
                return signature(self.resolve_T(ary))

            shape = normalize_shape(shape)
            if shape is None:
                return

            assert ary.ndim == shape.count
            return signature(self.resolve_T(ary).copy(layout="A"), shape)

        else:
            if any(not sentry_shape_scalar(a) for a in args):
                msg = "transpose({0}) is not supported".format(
                    ', '.join(args))
                raise TypingError(msg)
            assert ary.ndim == len(args)
            return signature(self.resolve_T(ary).copy(layout="A"), *args)

    @bound_function("array.copy")
    def resolve_copy(self, ary, args, kws):
        assert not args
        assert not kws
        retty = ary.copy(layout="C", readonly=False)
        return signature(retty)

    @bound_function("array.item")
    def resolve_item(self, ary, args, kws):
        assert not kws
        # We don't support explicit arguments as that's exactly equivalent
        # to regular indexing.  The no-argument form is interesting to
        # allow some degree of genericity when writing functions.
        if not args:
            return signature(ary.dtype)

    if numpy_version < (2, 0):
        @bound_function("array.itemset")
        def resolve_itemset(self, ary, args, kws):
            assert not kws
            # We don't support explicit arguments as that's exactly equivalent
            # to regular indexing.  The no-argument form is interesting to
            # allow some degree of genericity when writing functions.
            if len(args) == 1:
                return signature(types.none, ary.dtype)

    @bound_function("array.nonzero")
    def resolve_nonzero(self, ary, args, kws):
        assert not args
        assert not kws
        if ary.ndim == 0 and numpy_version >= (2, 1):
            raise NumbaValueError(
                "Calling nonzero on 0d arrays is not allowed."
                " Use np.atleast_1d(scalar).nonzero() instead."
            )
        # 0-dim arrays return one result array
        ndim = max(ary.ndim, 1)
        retty = types.UniTuple(types.Array(types.intp, 1, 'C'), ndim)
        return signature(retty)

    @bound_function("array.reshape")
    def resolve_reshape(self, ary, args, kws):
        def sentry_shape_scalar(ty):
            if ty in types.number_domain:
                # Guard against non integer type
                if not isinstance(ty, types.Integer):
                    raise TypingError("reshape() arg cannot be {0}".format(ty))
                return True
            else:
                return False

        assert not kws
        if ary.layout not in 'CF':
            # only work for contiguous array
            raise TypingError("reshape() supports contiguous array only")

        if len(args) == 1:
            # single arg
            shape, = args

            if sentry_shape_scalar(shape):
                ndim = 1
            else:
                shape = normalize_shape(shape)
                if shape is None:
                    return
                ndim = shape.count
            retty = ary.copy(ndim=ndim)
            return signature(retty, shape)

        elif len(args) == 0:
            # no arg
            raise TypingError("reshape() take at least one arg")

        else:
            # vararg case
            if any(not sentry_shape_scalar(a) for a in args):
                raise TypingError("reshape({0}) is not supported".format(
                    ', '.join(map(str, args))))

            retty = ary.copy(ndim=len(args))
            return signature(retty, *args)

    @bound_function("array.sort")
    def resolve_sort(self, ary, args, kws):
        assert not args
        assert not kws
        return signature(types.none)

    @bound_function("array.argsort")
    def resolve_argsort(self, ary, args, kws):
        assert not args
        kwargs = dict(kws)
        kind = kwargs.pop('kind', types.StringLiteral('quicksort'))
        if not isinstance(kind, types.StringLiteral):
            raise TypingError('"kind" must be a string literal')
        if kwargs:
            msg = "Unsupported keywords: {!r}"
            raise TypingError(msg.format([k for k in kwargs.keys()]))
        if ary.ndim == 1:
            def argsort_stub(kind='quicksort'):
                pass
            pysig = utils.pysignature(argsort_stub)
            sig = signature(types.Array(types.intp, 1, 'C'), kind).replace(pysig=pysig)
            return sig

    @bound_function("array.view")
    def resolve_view(self, ary, args, kws):
        from .npydecl import parse_dtype
        assert not kws
        dtype, = args
        dtype = parse_dtype(dtype)
        if dtype is None:
            return
        retty = ary.copy(dtype=dtype)
        return signature(retty, *args)

    @bound_function("array.astype")
    def resolve_astype(self, ary, args, kws):
        from .npydecl import parse_dtype
        assert not kws
        dtype, = args
        if isinstance(dtype, types.UnicodeType):
            raise RequireLiteralValue(("array.astype if dtype is a string it "
                                       "must be constant"))
        dtype = parse_dtype(dtype)
        if dtype is None:
            return
        if not self.context.can_convert(ary.dtype, dtype):
            raise TypingError("astype(%s) not supported on %s: "
                              "cannot convert from %s to %s"
                              % (dtype, ary, ary.dtype, dtype))
        layout = ary.layout if ary.layout in 'CF' else 'C'
        # reset the write bit irrespective of whether the cast type is the same
        # as the current dtype, this replicates numpy
        retty = ary.copy(dtype=dtype, layout=layout, readonly=False)
        return signature(retty, *args)

    @bound_function("array.ravel")
    def resolve_ravel(self, ary, args, kws):
        # Only support no argument version (default order='C')
        assert not kws
        assert not args
        copy_will_be_made = ary.layout != 'C'
        readonly = not (copy_will_be_made or ary.mutable)
        return signature(ary.copy(ndim=1, layout='C', readonly=readonly))

    @bound_function("array.flatten")
    def resolve_flatten(self, ary, args, kws):
        # Only support no argument version (default order='C')
        assert not kws
        assert not args
        # To ensure that Numba behaves exactly like NumPy,
        # we also clear the read-only flag when doing a "flatten"
        # Why? Two reasons:
        # Because flatten always returns a copy. (see NumPy docs for "flatten")
        # And because a copy always returns a writeable array.
        # ref: https://numpy.org/doc/stable/reference/generated/numpy.copy.html
        return signature(ary.copy(ndim=1, layout='C', readonly=False))

    def generic_resolve(self, ary, attr):
        # Resolution of other attributes, for record arrays
        if isinstance(ary.dtype, types.Record):
            if attr in ary.dtype.fields:
                attr_dtype = ary.dtype.typeof(attr)
                if isinstance(attr_dtype, types.NestedArray):
                    return ary.copy(
                        dtype=attr_dtype.dtype,
                        ndim=ary.ndim + attr_dtype.ndim,
                        layout='A'
                    )
                else:
                    return ary.copy(dtype=attr_dtype, layout='A')


@infer_getattr
class DTypeAttr(AttributeTemplate):
    key = types.DType

    def resolve_type(self, ary):
        # Wrap the numeric type in NumberClass
        return types.NumberClass(ary.dtype)

    def resolve_kind(self, ary):
        if isinstance(ary.key, types.scalars.Float):
            val = 'f'
        elif isinstance(ary.key, types.scalars.Integer):
            val = 'i'
        else:
            return None  # other types not supported yet
        return types.StringLiteral(val)


@infer
class StaticGetItemArray(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        # Resolution of members for record and structured arrays
        ary, idx = args
        if (isinstance(ary, types.Array) and isinstance(idx, str) and
                isinstance(ary.dtype, types.Record)):
            if idx in ary.dtype.fields:
                attr_dtype = ary.dtype.typeof(idx)
                if isinstance(attr_dtype, types.NestedArray):
                    ret = ary.copy(
                        dtype=attr_dtype.dtype,
                        ndim=ary.ndim + attr_dtype.ndim,
                        layout='A'
                    )
                    return signature(ret, *args)
                else:
                    ret = ary.copy(dtype=attr_dtype, layout='A')
                    return signature(ret, *args)


@infer_getattr
class RecordAttribute(AttributeTemplate):
    key = types.Record

    def generic_resolve(self, record, attr):
        ret = record.typeof(attr)
        assert ret
        return ret


@infer
class StaticGetItemRecord(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        # Resolution of members for records
        record, idx = args
        if isinstance(record, types.Record) and isinstance(idx, str):
            if idx not in record.fields:
                raise NumbaKeyError(f"Field '{idx}' was not found in record "
                                    "with fields "
                                    f"{tuple(record.fields.keys())}")
            ret = record.typeof(idx)
            assert ret
            return signature(ret, *args)


@infer_global(operator.getitem)
class StaticGetItemLiteralRecord(AbstractTemplate):
    def generic(self, args, kws):
        # Resolution of members for records
        record, idx = args
        if isinstance(record, types.Record):
            if isinstance(idx, types.StringLiteral):
                if idx.literal_value not in record.fields:
                    msg = (f"Field '{idx.literal_value}' was not found in "
                           f"record with fields {tuple(record.fields.keys())}")
                    raise NumbaKeyError(msg)
                ret = record.typeof(idx.literal_value)
                assert ret
                return signature(ret, *args)
            elif isinstance(idx, types.IntegerLiteral):
                if idx.literal_value >= len(record.fields):
                    msg = f"Requested index {idx.literal_value} is out of range"
                    raise NumbaIndexError(msg)
                field_names = list(record.fields)
                ret = record.typeof(field_names[idx.literal_value])
                assert ret
                return signature(ret, *args)


@infer
class StaticSetItemRecord(AbstractTemplate):
    key = "static_setitem"

    def generic(self, args, kws):
        # Resolution of members for record and structured arrays
        record, idx, value = args
        if isinstance(record, types.Record):
            if isinstance(idx, str):
                expectedty = record.typeof(idx)
                if self.context.can_convert(value, expectedty) is not None:
                    return signature(types.void, record, types.literal(idx),
                                     value)
            elif isinstance(idx, int):
                if idx >= len(record.fields):
                    msg = f"Requested index {idx} is out of range"
                    raise NumbaIndexError(msg)
                str_field = list(record.fields)[idx]
                expectedty = record.typeof(str_field)
                if self.context.can_convert(value, expectedty) is not None:
                    return signature(types.void, record, types.literal(idx),
                                     value)


@infer_global(operator.setitem)
class StaticSetItemLiteralRecord(AbstractTemplate):
    def generic(self, args, kws):
        # Resolution of members for records
        target, idx, value = args
        if isinstance(target, types.Record) and isinstance(idx, types.StringLiteral):
            if idx.literal_value not in target.fields:
                msg = (f"Field '{idx.literal_value}' was not found in record "
                       f"with fields {tuple(target.fields.keys())}")
                raise NumbaKeyError(msg)
            expectedty = target.typeof(idx.literal_value)
            if self.context.can_convert(value, expectedty) is not None:
                return signature(types.void, target, idx, value)


@infer_getattr
class ArrayCTypesAttribute(AttributeTemplate):
    key = types.ArrayCTypes

    def resolve_data(self, ctinfo):
        return types.uintp


@infer_getattr
class ArrayFlagsAttribute(AttributeTemplate):
    key = types.ArrayFlags

    def resolve_contiguous(self, ctflags):
        return types.boolean

    def resolve_c_contiguous(self, ctflags):
        return types.boolean

    def resolve_f_contiguous(self, ctflags):
        return types.boolean


@infer_getattr
class NestedArrayAttribute(ArrayAttribute):
    key = types.NestedArray


def _expand_integer(ty):
    """
    If *ty* is an integer, expand it to a machine int (like Numpy).
    """
    if isinstance(ty, types.Integer):
        if ty.signed:
            return max(types.intp, ty)
        else:
            return max(types.uintp, ty)
    elif isinstance(ty, types.Boolean):
        return types.intp
    else:
        return ty


def generic_homog(self, args, kws):
    if args:
        raise NumbaAssertionError("args not supported")
    if kws:
        raise NumbaAssertionError("kws not supported")

    return signature(self.this.dtype, recvr=self.this)


def generic_expand(self, args, kws):
    assert not args
    assert not kws
    return signature(_expand_integer(self.this.dtype), recvr=self.this)


def sum_expand(self, args, kws):
    """
    sum can be called with or without an axis parameter, and with or without
    a dtype parameter
    """
    pysig = None
    if 'axis' in kws and 'dtype' not in kws:
        def sum_stub(axis):
            pass
        pysig = utils.pysignature(sum_stub)
        # rewrite args
        args = list(args) + [kws['axis']]
    elif 'dtype' in kws and 'axis' not in kws:
        def sum_stub(dtype):
            pass
        pysig = utils.pysignature(sum_stub)
        # rewrite args
        args = list(args) + [kws['dtype']]
    elif 'dtype' in kws and 'axis' in kws:
        def sum_stub(axis, dtype):
            pass
        pysig = utils.pysignature(sum_stub)
        # rewrite args
        args = list(args) + [kws['axis'], kws['dtype']]

    args_len = len(args)
    assert args_len <= 2
    if args_len == 0:
        # No axis or dtype parameter so the return type of the summation is a scalar
        # of the type of the array.
        out = signature(_expand_integer(self.this.dtype), *args,
                        recvr=self.this)
    elif args_len == 1 and 'dtype' not in kws:
        # There is an axis parameter, either arg or kwarg
        if self.this.ndim == 1:
            # 1d reduces to a scalar
            return_type = _expand_integer(self.this.dtype)
        else:
            # the return type of this summation is  an array of dimension one
            # less than the input array.
            return_type = types.Array(dtype=_expand_integer(self.this.dtype),
                                    ndim=self.this.ndim-1, layout='C')
        out = signature(return_type, *args, recvr=self.this)

    elif args_len == 1 and 'dtype' in kws:
        # No axis parameter so the return type of the summation is a scalar
        # of the dtype parameter.
        from .npydecl import parse_dtype
        dtype, = args
        dtype = parse_dtype(dtype)
        out = signature(dtype, *args, recvr=self.this)

    elif args_len == 2:
        # There is an axis and dtype parameter, either arg or kwarg
        from .npydecl import parse_dtype
        dtype = parse_dtype(args[1])
        return_type = dtype
        if self.this.ndim != 1:
            # 1d reduces to a scalar, 2d and above reduce dim by 1
            # the return type of this summation is  an array of dimension one
            # less than the input array.
            return_type = types.Array(dtype=return_type,
                                    ndim=self.this.ndim-1, layout='C')
        out = signature(return_type, *args, recvr=self.this)
    else:
        pass
    return out.replace(pysig=pysig)


def generic_expand_cumulative(self, args, kws):
    if args:
        raise NumbaAssertionError("args unsupported")
    if kws:
        raise NumbaAssertionError("kwargs unsupported")
    assert isinstance(self.this, types.Array)
    return_type = types.Array(dtype=_expand_integer(self.this.dtype),
                              ndim=1, layout='C')
    return signature(return_type, recvr=self.this)


def generic_hetero_real(self, args, kws):
    assert not args
    assert not kws
    if isinstance(self.this.dtype, (types.Integer, types.Boolean)):
        return signature(types.float64, recvr=self.this)
    return signature(self.this.dtype, recvr=self.this)


def generic_hetero_always_real(self, args, kws):
    assert not args
    assert not kws
    if isinstance(self.this.dtype, (types.Integer, types.Boolean)):
        return signature(types.float64, recvr=self.this)
    if isinstance(self.this.dtype, types.Complex):
        return signature(self.this.dtype.underlying_float, recvr=self.this)
    return signature(self.this.dtype, recvr=self.this)


def generic_index(self, args, kws):
    assert not args
    assert not kws
    return signature(types.intp, recvr=self.this)


def install_array_method(name, generic, prefer_literal=True):
    my_attr = {"key": "array." + name, "generic": generic,
               "prefer_literal": prefer_literal}
    temp_class = type("Array_" + name, (AbstractTemplate,), my_attr)

    def array_attribute_attachment(self, ary):
        return types.BoundFunction(temp_class, ary)

    setattr(ArrayAttribute, "resolve_" + name, array_attribute_attachment)


# Functions that return a machine-width type, to avoid overflows
install_array_method("sum", sum_expand, prefer_literal=True)


@infer_global(operator.eq)
class CmpOpEqArray(AbstractTemplate):
    #key = operator.eq

    def generic(self, args, kws):
        assert not kws
        [va, vb] = args
        if isinstance(va, types.Array) and va == vb:
            return signature(va.copy(dtype=types.boolean), va, vb)
