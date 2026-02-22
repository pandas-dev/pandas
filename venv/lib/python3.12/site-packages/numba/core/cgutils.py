"""
Generic helpers for LLVM code generation.
"""


import collections
from contextlib import contextmanager, ExitStack
import functools

from llvmlite import ir

from numba.core import utils, types, config, debuginfo
import numba.core.datamodel


bool_t = ir.IntType(1)
int8_t = ir.IntType(8)
int32_t = ir.IntType(32)
intp_t = ir.IntType(utils.MACHINE_BITS)
voidptr_t = int8_t.as_pointer()

true_bit = bool_t(1)
false_bit = bool_t(0)
true_byte = int8_t(1)
false_byte = int8_t(0)


def as_bool_bit(builder, value):
    return builder.icmp_unsigned('!=', value, value.type(0))


def make_anonymous_struct(builder, values, struct_type=None):
    """
    Create an anonymous struct containing the given LLVM *values*.
    """
    if struct_type is None:
        struct_type = ir.LiteralStructType([v.type for v in values])
    struct_val = struct_type(ir.Undefined)
    for i, v in enumerate(values):
        struct_val = builder.insert_value(struct_val, v, i)
    return struct_val


def make_bytearray(buf):
    """
    Make a byte array constant from *buf*.
    """
    b = bytearray(buf)
    n = len(b)
    return ir.Constant(ir.ArrayType(ir.IntType(8), n), b)


_struct_proxy_cache = {}


def create_struct_proxy(fe_type, kind='value'):
    """
    Returns a specialized StructProxy subclass for the given fe_type.
    """
    cache_key = (fe_type, kind)
    res = _struct_proxy_cache.get(cache_key)
    if res is None:
        base = {'value': ValueStructProxy,
                'data': DataStructProxy,
                }[kind]
        clsname = base.__name__ + '_' + str(fe_type)
        bases = (base,)
        clsmembers = dict(_fe_type=fe_type)
        res = type(clsname, bases, clsmembers)

        _struct_proxy_cache[cache_key] = res
    return res


def copy_struct(dst, src, repl=None):
    """
    Copy structure from *src* to *dst* with replacement from *repl*.
    """
    if repl is None:
        repl = {}
    repl = repl.copy()
    # copy data from src or use those in repl
    for k in src._datamodel._fields:
        v = repl.pop(k, getattr(src, k))
        setattr(dst, k, v)
    # use remaining key-values in repl
    for k, v in repl.items():
        setattr(dst, k, v)
    return dst


class _StructProxy(object):
    """
    Creates a `Structure` like interface that is constructed with information
    from DataModel instance.  FE type must have a data model that is a
    subclass of StructModel.
    """
    # The following class members must be overridden by subclass
    _fe_type = None

    def __init__(self, context, builder, value=None, ref=None):
        self._context = context
        self._datamodel = self._context.data_model_manager[self._fe_type]
        if not isinstance(self._datamodel, numba.core.datamodel.StructModel):
            raise TypeError(
                "Not a structure model: {0}".format(self._datamodel))
        self._builder = builder

        self._be_type = self._get_be_type(self._datamodel)
        assert not is_pointer(self._be_type)

        outer_ref, ref = self._make_refs(ref)
        if ref.type.pointee != self._be_type:
            raise AssertionError("bad ref type: expected %s, got %s"
                                 % (self._be_type.as_pointer(), ref.type))

        if value is not None:
            if value.type != outer_ref.type.pointee:
                raise AssertionError("bad value type: expected %s, got %s"
                                     % (outer_ref.type.pointee, value.type))
            self._builder.store(value, outer_ref)

        self._value = ref
        self._outer_ref = outer_ref

    def _make_refs(self, ref):
        """
        Return an (outer ref, value ref) pair.  By default, these are
        the same pointers, but a derived class may override this.
        """
        if ref is None:
            ref = alloca_once(self._builder, self._be_type, zfill=True)
        return ref, ref

    def _get_be_type(self, datamodel):
        raise NotImplementedError

    def _cast_member_to_value(self, index, val):
        raise NotImplementedError

    def _cast_member_from_value(self, index, val):
        raise NotImplementedError

    def _get_ptr_by_index(self, index):
        return gep_inbounds(self._builder, self._value, 0, index)

    def _get_ptr_by_name(self, attrname):
        index = self._datamodel.get_field_position(attrname)
        return self._get_ptr_by_index(index)

    def __getattr__(self, field):
        """
        Load the LLVM value of the named *field*.
        """
        if not field.startswith('_'):
            return self[self._datamodel.get_field_position(field)]
        else:
            raise AttributeError(field)

    def __setattr__(self, field, value):
        """
        Store the LLVM *value* into the named *field*.
        """
        if field.startswith('_'):
            return super(_StructProxy, self).__setattr__(field, value)
        self[self._datamodel.get_field_position(field)] = value

    def __getitem__(self, index):
        """
        Load the LLVM value of the field at *index*.
        """
        member_val = self._builder.load(self._get_ptr_by_index(index))
        return self._cast_member_to_value(index, member_val)

    def __setitem__(self, index, value):
        """
        Store the LLVM *value* into the field at *index*.
        """
        ptr = self._get_ptr_by_index(index)
        value = self._cast_member_from_value(index, value)
        if value.type != ptr.type.pointee:
            if (is_pointer(value.type) and is_pointer(ptr.type.pointee)
                    and value.type.pointee == ptr.type.pointee.pointee):
                # Differ by address-space only
                # Auto coerce it
                value = self._context.addrspacecast(self._builder,
                                                    value,
                                                    ptr.type.pointee.addrspace)
            else:
                raise TypeError("Invalid store of {value.type} to "
                                "{ptr.type.pointee} in "
                                "{self._datamodel} "
                                "(trying to write member #{index})"
                                .format(value=value, ptr=ptr, self=self,
                                        index=index))
        self._builder.store(value, ptr)

    def __len__(self):
        """
        Return the number of fields.
        """
        return self._datamodel.field_count

    def _getpointer(self):
        """
        Return the LLVM pointer to the underlying structure.
        """
        return self._outer_ref

    def _getvalue(self):
        """
        Load and return the value of the underlying LLVM structure.
        """
        return self._builder.load(self._outer_ref)

    def _setvalue(self, value):
        """
        Store the value in this structure.
        """
        assert not is_pointer(value.type)
        assert value.type == self._be_type, (value.type, self._be_type)
        self._builder.store(value, self._value)


class ValueStructProxy(_StructProxy):
    """
    Create a StructProxy suitable for accessing regular values
    (e.g. LLVM values or alloca slots).
    """
    def _get_be_type(self, datamodel):
        return datamodel.get_value_type()

    def _cast_member_to_value(self, index, val):
        return val

    def _cast_member_from_value(self, index, val):
        return val


class DataStructProxy(_StructProxy):
    """
    Create a StructProxy suitable for accessing data persisted in memory.
    """
    def _get_be_type(self, datamodel):
        return datamodel.get_data_type()

    def _cast_member_to_value(self, index, val):
        model = self._datamodel.get_model(index)
        return model.from_data(self._builder, val)

    def _cast_member_from_value(self, index, val):
        model = self._datamodel.get_model(index)
        return model.as_data(self._builder, val)


class Structure(object):
    """
    A high-level object wrapping a alloca'ed LLVM structure, including
    named fields and attribute access.
    """

    # XXX Should this warrant several separate constructors?
    def __init__(self, context, builder, value=None, ref=None, cast_ref=False):
        self._type = context.get_struct_type(self)
        self._context = context
        self._builder = builder
        if ref is None:
            self._value = alloca_once(builder, self._type, zfill=True)
            if value is not None:
                assert not is_pointer(value.type)
                assert value.type == self._type, (value.type, self._type)
                builder.store(value, self._value)
        else:
            assert value is None
            assert is_pointer(ref.type)
            if self._type != ref.type.pointee:
                if cast_ref:
                    ref = builder.bitcast(ref, self._type.as_pointer())
                else:
                    raise TypeError(
                        "mismatching pointer type: got %s, expected %s"
                        % (ref.type.pointee, self._type))
            self._value = ref

        self._namemap = {}
        self._fdmap = []
        self._typemap = []
        base = int32_t(0)
        for i, (k, tp) in enumerate(self._fields):
            self._namemap[k] = i
            self._fdmap.append((base, int32_t(i)))
            self._typemap.append(tp)

    def _get_ptr_by_index(self, index):
        ptr = self._builder.gep(self._value, self._fdmap[index], inbounds=True)
        return ptr

    def _get_ptr_by_name(self, attrname):
        return self._get_ptr_by_index(self._namemap[attrname])

    def __getattr__(self, field):
        """
        Load the LLVM value of the named *field*.
        """
        if not field.startswith('_'):
            return self[self._namemap[field]]
        else:
            raise AttributeError(field)

    def __setattr__(self, field, value):
        """
        Store the LLVM *value* into the named *field*.
        """
        if field.startswith('_'):
            return super(Structure, self).__setattr__(field, value)
        self[self._namemap[field]] = value

    def __getitem__(self, index):
        """
        Load the LLVM value of the field at *index*.
        """

        return self._builder.load(self._get_ptr_by_index(index))

    def __setitem__(self, index, value):
        """
        Store the LLVM *value* into the field at *index*.
        """
        ptr = self._get_ptr_by_index(index)
        if ptr.type.pointee != value.type:
            fmt = "Type mismatch: __setitem__(%d, ...) expected %r but got %r"
            raise AssertionError(fmt % (index,
                                        str(ptr.type.pointee),
                                        str(value.type)))
        self._builder.store(value, ptr)

    def __len__(self):
        """
        Return the number of fields.
        """
        return len(self._namemap)

    def _getpointer(self):
        """
        Return the LLVM pointer to the underlying structure.
        """
        return self._value

    def _getvalue(self):
        """
        Load and return the value of the underlying LLVM structure.
        """
        return self._builder.load(self._value)

    def _setvalue(self, value):
        """Store the value in this structure"""
        assert not is_pointer(value.type)
        assert value.type == self._type, (value.type, self._type)
        self._builder.store(value, self._value)

    # __iter__ is derived by Python from __len__ and __getitem__


def alloca_once(builder, ty, size=None, name='', zfill=False):
    """Allocate stack memory at the entry block of the current function
    pointed by ``builder`` with llvm type ``ty``.  The optional ``size`` arg
    set the number of element to allocate.  The default is 1.  The optional
    ``name`` arg set the symbol name inside the llvm IR for debugging.
    If ``zfill`` is set, fill the memory with zeros at the current
    use-site location.  Note that the memory is always zero-filled after the
    ``alloca`` at init-site (the entry block).
    """
    if isinstance(size, int):
        size = ir.Constant(intp_t, size)
    # suspend debug metadata emission else it links up python source lines with
    # alloca in the entry block as well as their actual location and it makes
    # the debug info "jump about".
    with debuginfo.suspend_emission(builder):
        with builder.goto_entry_block():
            ptr = builder.alloca(ty, size=size, name=name)
            # Always zero-fill at init-site.  This is safe.
            builder.store(ty(None), ptr)
        # Also zero-fill at the use-site
        if zfill:
            builder.store(ptr.type.pointee(None), ptr)
        return ptr


def sizeof(builder, ptr_type):
    """Compute sizeof using GEP
    """
    null = ptr_type(None)
    offset = null.gep([int32_t(1)])
    return builder.ptrtoint(offset, intp_t)


def alloca_once_value(builder, value, name='', zfill=False):
    """
    Like alloca_once(), but passing a *value* instead of a type.  The
    type is inferred and the allocated slot is also initialized with the
    given value.
    """
    storage = alloca_once(builder, value.type, zfill=zfill)
    builder.store(value, storage)
    return storage


def insert_pure_function(module, fnty, name):
    """
    Insert a pure function (in the functional programming sense) in the
    given module.
    """
    fn = get_or_insert_function(module, fnty, name)
    fn.attributes.add("readonly")
    fn.attributes.add("nounwind")
    return fn


def get_or_insert_function(module, fnty, name):
    """
    Get the function named *name* with type *fnty* from *module*, or insert it
    if it doesn't exist.
    """
    fn = module.globals.get(name, None)
    if fn is None:
        fn = ir.Function(module, fnty, name)
    return fn


def get_or_insert_named_metadata(module, name):
    try:
        return module.get_named_metadata(name)
    except KeyError:
        return module.add_named_metadata(name)


def add_global_variable(module, ty, name, addrspace=0):
    unique_name = module.get_unique_name(name)
    return ir.GlobalVariable(module, ty, unique_name, addrspace)


def terminate(builder, bbend):
    bb = builder.basic_block
    if bb.terminator is None:
        builder.branch(bbend)


def get_null_value(ltype):
    return ltype(None)


def is_null(builder, val):
    null = get_null_value(val.type)
    return builder.icmp_unsigned('==', null, val)


def is_not_null(builder, val):
    null = get_null_value(val.type)
    return builder.icmp_unsigned('!=', null, val)


def if_unlikely(builder, pred):
    return builder.if_then(pred, likely=False)


def if_likely(builder, pred):
    return builder.if_then(pred, likely=True)


def ifnot(builder, pred):
    return builder.if_then(builder.not_(pred))


def increment_index(builder, val):
    """
    Increment an index *val*.
    """
    one = val.type(1)
    # We pass the "nsw" flag in the hope that LLVM understands the index
    # never changes sign.  Unfortunately this doesn't always work
    # (e.g. ndindex()).
    return builder.add(val, one, flags=['nsw'])


Loop = collections.namedtuple('Loop', ('index', 'do_break'))


@contextmanager
def for_range(builder, count, start=None, intp=None):
    """
    Generate LLVM IR for a for-loop in [start, count).
    *start* is equal to 0 by default.

    Yields a Loop namedtuple with the following members:
    - `index` is the loop index's value
    - `do_break` is a no-argument callable to break out of the loop
    """
    if intp is None:
        intp = count.type
    if start is None:
        start = intp(0)
    stop = count

    bbcond = builder.append_basic_block("for.cond")
    bbbody = builder.append_basic_block("for.body")
    bbend = builder.append_basic_block("for.end")

    def do_break():
        builder.branch(bbend)

    bbstart = builder.basic_block
    builder.branch(bbcond)

    with builder.goto_block(bbcond):
        index = builder.phi(intp, name="loop.index")
        pred = builder.icmp_signed('<', index, stop)
        builder.cbranch(pred, bbbody, bbend)

    with builder.goto_block(bbbody):
        yield Loop(index, do_break)
        # Update bbbody as a new basic block may have been activated
        bbbody = builder.basic_block
        incr = increment_index(builder, index)
        terminate(builder, bbcond)

    index.add_incoming(start, bbstart)
    index.add_incoming(incr, bbbody)

    builder.position_at_end(bbend)


@contextmanager
def for_range_slice(builder, start, stop, step, intp=None, inc=True):
    """
    Generate LLVM IR for a for-loop based on a slice.  Yields a
    (index, count) tuple where `index` is the slice index's value
    inside the loop, and `count` the iteration count.

    Parameters
    -------------
    builder : object
        IRBuilder object
    start : int
        The beginning value of the slice
    stop : int
        The end value of the slice
    step : int
        The step value of the slice
    intp :
        The data type
    inc : boolean, optional
        Signals whether the step is positive (True) or negative (False).

    Returns
    -----------
        None
    """
    if intp is None:
        intp = start.type

    bbcond = builder.append_basic_block("for.cond")
    bbbody = builder.append_basic_block("for.body")
    bbend = builder.append_basic_block("for.end")
    bbstart = builder.basic_block
    builder.branch(bbcond)

    with builder.goto_block(bbcond):
        index = builder.phi(intp, name="loop.index")
        count = builder.phi(intp, name="loop.count")
        if (inc):
            pred = builder.icmp_signed('<', index, stop)
        else:
            pred = builder.icmp_signed('>', index, stop)
        builder.cbranch(pred, bbbody, bbend)

    with builder.goto_block(bbbody):
        yield index, count
        bbbody = builder.basic_block
        incr = builder.add(index, step)
        next_count = increment_index(builder, count)
        terminate(builder, bbcond)

    index.add_incoming(start, bbstart)
    index.add_incoming(incr, bbbody)
    count.add_incoming(ir.Constant(intp, 0), bbstart)
    count.add_incoming(next_count, bbbody)
    builder.position_at_end(bbend)


@contextmanager
def for_range_slice_generic(builder, start, stop, step):
    """
    A helper wrapper for for_range_slice().  This is a context manager which
    yields two for_range_slice()-alike context managers, the first for
    the positive step case, the second for the negative step case.

    Use:
        with for_range_slice_generic(...) as (pos_range, neg_range):
            with pos_range as (idx, count):
                ...
            with neg_range as (idx, count):
                ...
    """
    intp = start.type
    is_pos_step = builder.icmp_signed('>=', step, ir.Constant(intp, 0))

    pos_for_range = for_range_slice(builder, start, stop, step, intp, inc=True)
    neg_for_range = for_range_slice(builder, start, stop, step, intp, inc=False)

    @contextmanager
    def cm_cond(cond, inner_cm):
        with cond:
            with inner_cm as value:
                yield value

    with builder.if_else(is_pos_step, likely=True) as (then, otherwise):
        yield cm_cond(then, pos_for_range), cm_cond(otherwise, neg_for_range)


@contextmanager
def loop_nest(builder, shape, intp, order='C'):
    """
    Generate a loop nest walking a N-dimensional array.
    Yields a tuple of N indices for use in the inner loop body,
    iterating over the *shape* space.

    If *order* is 'C' (the default), indices are incremented inside-out
    (i.e. (0,0), (0,1), (0,2), (1,0) etc.).
    If *order* is 'F', they are incremented outside-in
    (i.e. (0,0), (1,0), (2,0), (0,1) etc.).
    This has performance implications when walking an array as it impacts
    the spatial locality of memory accesses.
    """
    assert order in 'CF'
    if not shape:
        # 0-d array
        yield ()
    else:
        if order == 'F':
            _swap = lambda x: x[::-1]
        else:
            _swap = lambda x: x
        with _loop_nest(builder, _swap(shape), intp) as indices:
            assert len(indices) == len(shape)
            yield _swap(indices)


@contextmanager
def _loop_nest(builder, shape, intp):
    with for_range(builder, shape[0], intp=intp) as loop:
        if len(shape) > 1:
            with _loop_nest(builder, shape[1:], intp) as indices:
                yield (loop.index,) + indices
        else:
            yield (loop.index,)


def pack_array(builder, values, ty=None):
    """
    Pack a sequence of values in a LLVM array.  *ty* should be given
    if the array may be empty, in which case the type can't be inferred
    from the values.
    """
    n = len(values)
    if ty is None:
        ty = values[0].type
    ary = ir.ArrayType(ty, n)(ir.Undefined)
    for i, v in enumerate(values):
        ary = builder.insert_value(ary, v, i)
    return ary


def pack_struct(builder, values):
    """
    Pack a sequence of values into a LLVM struct.
    """
    structty = ir.LiteralStructType([v.type for v in values])
    st = structty(ir.Undefined)
    for i, v in enumerate(values):
        st = builder.insert_value(st, v, i)
    return st


def unpack_tuple(builder, tup, count=None):
    """
    Unpack an array or structure of values, return a Python tuple.
    """
    if count is None:
        # Assuming *tup* is an aggregate
        count = len(tup.type.elements)
    vals = [builder.extract_value(tup, i)
            for i in range(count)]
    return vals


def get_item_pointer(context, builder, aryty, ary, inds, wraparound=False,
                     boundscheck=False):
    # Set boundscheck=True for any pointer access that should be
    # boundschecked. do_boundscheck() will handle enabling or disabling the
    # actual boundschecking based on the user config.
    shapes = unpack_tuple(builder, ary.shape, count=aryty.ndim)
    strides = unpack_tuple(builder, ary.strides, count=aryty.ndim)
    return get_item_pointer2(context, builder, data=ary.data, shape=shapes,
                             strides=strides, layout=aryty.layout, inds=inds,
                             wraparound=wraparound, boundscheck=boundscheck)


def do_boundscheck(context, builder, ind, dimlen, axis=None):
    def _dbg():
        # Remove this when we figure out how to include this information
        # in the error message.
        if axis is not None:
            if isinstance(axis, int):
                printf(builder, "debug: IndexError: index %d is out of bounds "
                       "for axis {} with size %d\n".format(axis), ind, dimlen)
            else:
                printf(builder, "debug: IndexError: index %d is out of bounds "
                       "for axis %d with size %d\n", ind, axis,
                       dimlen)
        else:
            printf(builder,
                   "debug: IndexError: index %d is out of bounds for size %d\n",
                   ind, dimlen)

    msg = "index is out of bounds"
    out_of_bounds_upper = builder.icmp_signed('>=', ind, dimlen)
    with if_unlikely(builder, out_of_bounds_upper):
        if config.FULL_TRACEBACKS:
            _dbg()
        context.call_conv.return_user_exc(builder, IndexError, (msg,))
    out_of_bounds_lower = builder.icmp_signed('<', ind, ind.type(0))
    with if_unlikely(builder, out_of_bounds_lower):
        if config.FULL_TRACEBACKS:
            _dbg()
        context.call_conv.return_user_exc(builder, IndexError, (msg,))


def get_item_pointer2(context, builder, data, shape, strides, layout, inds,
                      wraparound=False, boundscheck=False):
    # Set boundscheck=True for any pointer access that should be
    # boundschecked. do_boundscheck() will handle enabling or disabling the
    # actual boundschecking based on the user config.
    if wraparound:
        # Wraparound
        indices = []
        for ind, dimlen in zip(inds, shape):
            negative = builder.icmp_signed('<', ind, ind.type(0))
            wrapped = builder.add(dimlen, ind)
            selected = builder.select(negative, wrapped, ind)
            indices.append(selected)
    else:
        indices = inds
    if boundscheck:
        for axis, (ind, dimlen) in enumerate(zip(indices, shape)):
            do_boundscheck(context, builder, ind, dimlen, axis)

    if not indices:
        # Indexing with empty tuple
        return builder.gep(data, [int32_t(0)])
    intp = indices[0].type
    # Indexing code
    if layout in 'CF':
        steps = []
        # Compute steps for each dimension
        if layout == 'C':
            # C contiguous
            for i in range(len(shape)):
                last = intp(1)
                for j in shape[i + 1:]:
                    last = builder.mul(last, j)
                steps.append(last)
        elif layout == 'F':
            # F contiguous
            for i in range(len(shape)):
                last = intp(1)
                for j in shape[:i]:
                    last = builder.mul(last, j)
                steps.append(last)
        else:
            raise Exception("unreachable")

        # Compute index
        loc = intp(0)
        for i, s in zip(indices, steps):
            tmp = builder.mul(i, s)
            loc = builder.add(loc, tmp)
        ptr = builder.gep(data, [loc])
        return ptr
    else:
        # Any layout
        dimoffs = [builder.mul(s, i) for s, i in zip(strides, indices)]
        offset = functools.reduce(builder.add, dimoffs)
        return pointer_add(builder, data, offset)


def _scalar_pred_against_zero(builder, value, fpred, icond):
    nullval = value.type(0)
    if isinstance(value.type, (ir.FloatType, ir.DoubleType)):
        isnull = fpred(value, nullval)
    elif isinstance(value.type, ir.IntType):
        isnull = builder.icmp_signed(icond, value, nullval)
    else:
        raise TypeError("unexpected value type %s" % (value.type,))
    return isnull


def is_scalar_zero(builder, value):
    """
    Return a predicate representing whether *value* is equal to zero.
    """
    return _scalar_pred_against_zero(
        builder, value, functools.partial(builder.fcmp_ordered, '=='), '==')


def is_not_scalar_zero(builder, value):
    """
    Return a predicate representing whether a *value* is not equal to zero.
    (not exactly "not is_scalar_zero" because of nans)
    """
    return _scalar_pred_against_zero(
        builder, value, functools.partial(builder.fcmp_unordered, '!='), '!=')


def is_scalar_zero_or_nan(builder, value):
    """
    Return a predicate representing whether *value* is equal to either zero
    or NaN.
    """
    return _scalar_pred_against_zero(
        builder, value, functools.partial(builder.fcmp_unordered, '=='), '==')


is_true = is_not_scalar_zero
is_false = is_scalar_zero


def is_scalar_neg(builder, value):
    """
    Is *value* negative?  Assumes *value* is signed.
    """
    return _scalar_pred_against_zero(
        builder, value, functools.partial(builder.fcmp_ordered, '<'), '<')


@contextmanager
def early_exit_if(builder, stack: ExitStack, cond):
    """
    The Python code::

        with contextlib.ExitStack() as stack:
            with early_exit_if(builder, stack, cond):
                cleanup()
            body()

    emits the code::

        if (cond) {
            <cleanup>
        }
        else {
            <body>
        }

    This can be useful for generating code with lots of early exits, without
    having to increase the indentation each time.
    """
    then, otherwise = stack.enter_context(builder.if_else(cond, likely=False))
    with then:
        yield
    stack.enter_context(otherwise)


def early_exit_if_null(builder, stack, obj):
    """
    A convenience wrapper for :func:`early_exit_if`, for the common case where
    the CPython API indicates an error by returning ``NULL``.
    """
    return early_exit_if(builder, stack, is_null(builder, obj))


def guard_null(context, builder, value, exc_tuple):
    """
    Guard against *value* being null or zero.
    *exc_tuple* should be a (exception type, arguments...) tuple.
    """
    with builder.if_then(is_scalar_zero(builder, value), likely=False):
        exc = exc_tuple[0]
        exc_args = exc_tuple[1:] or None
        context.call_conv.return_user_exc(builder, exc, exc_args)


def guard_memory_error(context, builder, pointer, msg=None):
    """
    Guard against *pointer* being NULL (and raise a MemoryError).
    """
    assert isinstance(pointer.type, ir.PointerType), pointer.type
    exc_args = (msg,) if msg else ()
    with builder.if_then(is_null(builder, pointer), likely=False):
        context.call_conv.return_user_exc(builder, MemoryError, exc_args)


@contextmanager
def if_zero(builder, value, likely=False):
    """
    Execute the given block if the scalar value is zero.
    """
    with builder.if_then(is_scalar_zero(builder, value), likely=likely):
        yield


guard_zero = guard_null


def is_pointer(ltyp):
    """
    Whether the LLVM type *typ* is a struct type.
    """
    return isinstance(ltyp, ir.PointerType)


def get_record_member(builder, record, offset, typ):
    pval = gep_inbounds(builder, record, 0, offset)
    assert not is_pointer(pval.type.pointee)
    return builder.bitcast(pval, typ.as_pointer())


def is_neg_int(builder, val):
    return builder.icmp_signed('<', val, val.type(0))


def gep_inbounds(builder, ptr, *inds, **kws):
    """
    Same as *gep*, but add the `inbounds` keyword.
    """
    return gep(builder, ptr, *inds, inbounds=True, **kws)


def gep(builder, ptr, *inds, **kws):
    """
    Emit a getelementptr instruction for the given pointer and indices.
    The indices can be LLVM values or Python int constants.
    """
    name = kws.pop('name', '')
    inbounds = kws.pop('inbounds', False)
    assert not kws
    idx = []
    for i in inds:
        if isinstance(i, int):
            # NOTE: llvm only accepts int32 inside structs, not int64
            ind = int32_t(i)
        else:
            ind = i
        idx.append(ind)
    return builder.gep(ptr, idx, name=name, inbounds=inbounds)


def pointer_add(builder, ptr, offset, return_type=None):
    """
    Add an integral *offset* to pointer *ptr*, and return a pointer
    of *return_type* (or, if omitted, the same type as *ptr*).

    Note the computation is done in bytes, and ignores the width of
    the pointed item type.
    """
    intptr = builder.ptrtoint(ptr, intp_t)
    if isinstance(offset, int):
        offset = intp_t(offset)
    intptr = builder.add(intptr, offset)
    return builder.inttoptr(intptr, return_type or ptr.type)


def memset(builder, ptr, size, value):
    """
    Fill *size* bytes starting from *ptr* with *value*.
    """
    fn = builder.module.declare_intrinsic('llvm.memset', (voidptr_t, size.type))
    ptr = builder.bitcast(ptr, voidptr_t)
    if isinstance(value, int):
        value = int8_t(value)
    builder.call(fn, [ptr, value, size, bool_t(0)])


def memset_padding(builder, ptr):
    """
    Fill padding bytes of the pointee with zeros.
    """
    # Load existing value
    val = builder.load(ptr)
    # Fill pointee with zeros
    memset(builder, ptr, sizeof(builder, ptr.type), 0)
    # Store value back
    builder.store(val, ptr)


def global_constant(builder_or_module, name, value, linkage='internal'):
    """
    Get or create a (LLVM module-)global constant with *name* or *value*.
    """
    if isinstance(builder_or_module, ir.Module):
        module = builder_or_module
    else:
        module = builder_or_module.module
    data = add_global_variable(module, value.type, name)
    data.linkage = linkage
    data.global_constant = True
    data.initializer = value
    return data


def divmod_by_constant(builder, val, divisor):
    """
    Compute the (quotient, remainder) of *val* divided by the constant
    positive *divisor*.  The semantics reflects those of Python integer
    floor division, rather than C's / LLVM's signed division and modulo.
    The difference lies with a negative *val*.
    """
    assert divisor > 0
    divisor = val.type(divisor)
    one = val.type(1)

    quot = alloca_once(builder, val.type)

    with builder.if_else(is_neg_int(builder, val)) as (if_neg, if_pos):
        with if_pos:
            # quot = val / divisor
            quot_val = builder.sdiv(val, divisor)
            builder.store(quot_val, quot)
        with if_neg:
            # quot = -1 + (val + 1) / divisor
            val_plus_one = builder.add(val, one)
            quot_val = builder.sdiv(val_plus_one, divisor)
            builder.store(builder.sub(quot_val, one), quot)

    # rem = val - quot * divisor
    # (should be slightly faster than a separate modulo operation)
    quot_val = builder.load(quot)
    rem_val = builder.sub(val, builder.mul(quot_val, divisor))
    return quot_val, rem_val


def cbranch_or_continue(builder, cond, bbtrue):
    """
    Branch conditionally or continue.

    Note: a new block is created and builder is moved to the end of the new
          block.
    """
    bbcont = builder.append_basic_block('.continue')
    builder.cbranch(cond, bbtrue, bbcont)
    builder.position_at_end(bbcont)
    return bbcont


def memcpy(builder, dst, src, count):
    """
    Emit a memcpy to the builder.

    Copies each element of dst to src. Unlike the C equivalent, each element
    can be any LLVM type.

    Assumes
    -------
    * dst.type == src.type
    * count is positive
    """
    # Note this does seem to be optimized as a raw memcpy() by LLVM
    # whenever possible...
    assert dst.type == src.type
    with for_range(builder, count, intp=count.type) as loop:
        out_ptr = builder.gep(dst, [loop.index])
        in_ptr = builder.gep(src, [loop.index])
        builder.store(builder.load(in_ptr), out_ptr)


def _raw_memcpy(builder, func_name, dst, src, count, itemsize, align):
    size_t = count.type
    if isinstance(itemsize, int):
        itemsize = ir.Constant(size_t, itemsize)

    memcpy = builder.module.declare_intrinsic(func_name,
                                              [voidptr_t, voidptr_t, size_t])
    is_volatile = false_bit
    builder.call(memcpy, [builder.bitcast(dst, voidptr_t),
                          builder.bitcast(src, voidptr_t),
                          builder.mul(count, itemsize),
                          is_volatile])


def raw_memcpy(builder, dst, src, count, itemsize, align=1):
    """
    Emit a raw memcpy() call for `count` items of size `itemsize`
    from `src` to `dest`.
    """
    return _raw_memcpy(builder, 'llvm.memcpy', dst, src, count, itemsize, align)


def raw_memmove(builder, dst, src, count, itemsize, align=1):
    """
    Emit a raw memmove() call for `count` items of size `itemsize`
    from `src` to `dest`.
    """
    return _raw_memcpy(builder, 'llvm.memmove', dst, src, count,
                       itemsize, align)


def muladd_with_overflow(builder, a, b, c):
    """
    Compute (a * b + c) and return a (result, overflow bit) pair.
    The operands must be signed integers.
    """
    p = builder.smul_with_overflow(a, b)
    prod = builder.extract_value(p, 0)
    prod_ovf = builder.extract_value(p, 1)
    s = builder.sadd_with_overflow(prod, c)
    res = builder.extract_value(s, 0)
    ovf = builder.or_(prod_ovf, builder.extract_value(s, 1))
    return res, ovf


def printf(builder, format, *args):
    """
    Calls printf().
    Argument `format` is expected to be a Python string.
    Values to be printed are listed in `args`.

    Note: There is no checking to ensure there is correct number of values
    in `args` and there type matches the declaration in the format string.
    """
    assert isinstance(format, str)
    mod = builder.module
    # Make global constant for format string
    cstring = voidptr_t
    fmt_bytes = make_bytearray((format + '\00').encode('ascii'))
    global_fmt = global_constant(mod, "printf_format", fmt_bytes)
    fnty = ir.FunctionType(int32_t, [cstring], var_arg=True)
    # Insert printf()
    try:
        fn = mod.get_global('printf')
    except KeyError:
        fn = ir.Function(mod, fnty, name="printf")
    # Call
    ptr_fmt = builder.bitcast(global_fmt, cstring)
    return builder.call(fn, [ptr_fmt] + list(args))


def snprintf(builder, buffer, bufsz, format, *args):
    """Calls libc snprintf(buffer, bufsz, format, ...args)
    """
    assert isinstance(format, str)
    mod = builder.module
    # Make global constant for format string
    cstring = voidptr_t
    fmt_bytes = make_bytearray((format + '\00').encode('ascii'))
    global_fmt = global_constant(mod, "snprintf_format", fmt_bytes)
    fnty = ir.FunctionType(
        int32_t, [cstring, intp_t, cstring], var_arg=True,
    )
    # Actual symbol name of snprintf is different on win32.
    symbol = 'snprintf'
    if config.IS_WIN32:
        symbol = '_' + symbol
    # Insert snprintf()
    try:
        fn = mod.get_global(symbol)
    except KeyError:
        fn = ir.Function(mod, fnty, name=symbol)
    # Call
    ptr_fmt = builder.bitcast(global_fmt, cstring)
    return builder.call(fn, [buffer, bufsz, ptr_fmt] + list(args))


def snprintf_stackbuffer(builder, bufsz, format, *args):
    """Similar to `snprintf()` but the buffer is stack allocated to size
    *bufsz*.

    Returns the buffer pointer as i8*.
    """
    assert isinstance(bufsz, int)
    spacety = ir.ArrayType(ir.IntType(8), bufsz)
    space = alloca_once(builder, spacety, zfill=True)
    buffer = builder.bitcast(space, voidptr_t)
    snprintf(builder, buffer, intp_t(bufsz), format, *args)
    return buffer


def normalize_ir_text(text):
    """
    Normalize the given string to latin1 compatible encoding that is
    suitable for use in LLVM IR.
    """
    # Just re-encoding to latin1 is enough
    return text.encode('utf8').decode('latin1')


def hexdump(builder, ptr, nbytes):
    """Debug print the memory region in *ptr* to *ptr + nbytes*
    as hex.
    """
    bytes_per_line = 16
    nbytes = builder.zext(nbytes, intp_t)
    printf(builder, "hexdump p=%p n=%zu",
           ptr, nbytes)
    byte_t = ir.IntType(8)
    ptr = builder.bitcast(ptr, byte_t.as_pointer())
    # Loop to print the bytes in *ptr* as hex
    with for_range(builder, nbytes) as idx:
        div_by = builder.urem(idx.index, intp_t(bytes_per_line))
        do_new_line = builder.icmp_unsigned("==", div_by, intp_t(0))
        with builder.if_then(do_new_line):
            printf(builder, "\n")

        offset = builder.gep(ptr, [idx.index])
        val = builder.load(offset)
        printf(builder, " %02x", val)
    printf(builder, "\n")


def is_nonelike(ty):
    """ returns if 'ty' is none """
    return (
        ty is None or
        isinstance(ty, types.NoneType) or
        isinstance(ty, types.Omitted)
    )


def is_empty_tuple(ty):
    """ returns if 'ty' is an empty tuple """
    return (
        isinstance(ty, types.Tuple) and
        len(ty.types) == 0
    )


def create_constant_array(ty, val):
    """
    Create an LLVM-constant of a fixed-length array from Python values.

    The type provided is the type of the elements.
    """
    return ir.Constant(ir.ArrayType(ty, len(val)), val)
