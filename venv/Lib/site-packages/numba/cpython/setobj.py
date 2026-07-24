"""
Support for native homogeneous sets.
"""


import collections
import contextlib
import math
import operator
from functools import cached_property

from llvmlite import ir
from numba.core import types, typing, cgutils
from numba.core.imputils import (lower_builtin, lower_cast,
                                    iternext_impl, impl_ret_borrowed,
                                    impl_ret_new_ref, impl_ret_untracked,
                                    for_iter, call_len, RefType)
from numba.misc import quicksort
from numba.cpython import slicing
from numba.core.errors import NumbaValueError, TypingError
from numba.core.extending import overload, overload_method, intrinsic


def get_payload_struct(context, builder, set_type, ptr):
    """
    Given a set value and type, get its payload structure (as a
    reference, so that mutations are seen by all).
    """
    payload_type = types.SetPayload(set_type)
    ptrty = context.get_data_type(payload_type).as_pointer()
    payload = builder.bitcast(ptr, ptrty)
    return context.make_data_helper(builder, payload_type, ref=payload)


def get_entry_size(context, set_type):
    """
    Return the entry size for the given set type.
    """
    llty = context.get_data_type(types.SetEntry(set_type))
    return context.get_abi_sizeof(llty)


# Note these values are special:
# - EMPTY is obtained by issuing memset(..., 0xFF)
# - (unsigned) EMPTY > (unsigned) DELETED > any other hash value
EMPTY = -1
DELETED = -2
FALLBACK = -43

# Minimal size of entries table.  Must be a power of 2!
MINSIZE = 16

# Number of cache-friendly linear probes before switching to non-linear probing
LINEAR_PROBES = 3

DEBUG_ALLOCS = False


def get_hash_value(context, builder, typ, value):
    """
    Compute the hash of the given value.
    """
    typingctx = context.typing_context
    fnty = typingctx.resolve_value_type(hash)
    sig = fnty.get_call_type(typingctx, (typ,), {})
    fn = context.get_function(fnty, sig)
    h = fn(builder, (value,))
    # Fixup reserved values
    is_ok = is_hash_used(context, builder, h)
    fallback = ir.Constant(h.type, FALLBACK)
    return builder.select(is_ok, h, fallback)


@intrinsic
def _get_hash_value_intrinsic(typingctx, value):
    def impl(context, builder, typ, args):
        return get_hash_value(context, builder, value, args[0])
    fnty = typingctx.resolve_value_type(hash)
    sig = fnty.get_call_type(typingctx, (value,), {})
    return sig, impl


def is_hash_empty(context, builder, h):
    """
    Whether the hash value denotes an empty entry.
    """
    empty = ir.Constant(h.type, EMPTY)
    return builder.icmp_unsigned('==', h, empty)

def is_hash_deleted(context, builder, h):
    """
    Whether the hash value denotes a deleted entry.
    """
    deleted = ir.Constant(h.type, DELETED)
    return builder.icmp_unsigned('==', h, deleted)

def is_hash_used(context, builder, h):
    """
    Whether the hash value denotes an active entry.
    """
    # Everything below DELETED is an used entry
    deleted = ir.Constant(h.type, DELETED)
    return builder.icmp_unsigned('<', h, deleted)


def check_all_set(*args):
    if not all([isinstance(typ, types.Set) for typ in args]):
        raise TypingError(f"All arguments must be Sets, got {args}")

    if not all([args[0].dtype == s.dtype for s in args]):
        raise TypingError(f"All Sets must be of the same type, got {args}")


SetLoop = collections.namedtuple('SetLoop', ('index', 'entry', 'do_break'))


class _SetPayload(object):

    def __init__(self, context, builder, set_type, ptr):
        payload = get_payload_struct(context, builder, set_type, ptr)
        self._context = context
        self._builder = builder
        self._ty = set_type
        self._payload = payload
        self._entries = payload._get_ptr_by_name('entries')
        self._ptr = ptr

    @property
    def mask(self):
        return self._payload.mask

    @mask.setter
    def mask(self, value):
        # CAUTION: mask must be a power of 2 minus 1
        self._payload.mask = value

    @property
    def used(self):
        return self._payload.used

    @used.setter
    def used(self, value):
        self._payload.used = value

    @property
    def fill(self):
        return self._payload.fill

    @fill.setter
    def fill(self, value):
        self._payload.fill = value

    @property
    def finger(self):
        return self._payload.finger

    @finger.setter
    def finger(self, value):
        self._payload.finger = value

    @property
    def dirty(self):
        return self._payload.dirty

    @dirty.setter
    def dirty(self, value):
        self._payload.dirty = value

    @property
    def entries(self):
        """
        A pointer to the start of the entries array.
        """
        return self._entries

    @property
    def ptr(self):
        """
        A pointer to the start of the NRT-allocated area.
        """
        return self._ptr

    def get_entry(self, idx):
        """
        Get entry number *idx*.
        """
        entry_ptr = cgutils.gep(self._builder, self._entries, idx)
        entry = self._context.make_data_helper(self._builder,
                                               types.SetEntry(self._ty),
                                               ref=entry_ptr)
        return entry

    def _lookup(self, item, h, for_insert=False):
        """
        Lookup the *item* with the given hash values in the entries.

        Return a (found, entry index) tuple:
        - If found is true, <entry index> points to the entry containing
          the item.
        - If found is false, <entry index> points to the empty entry that
          the item can be written to (only if *for_insert* is true)
        """
        context = self._context
        builder = self._builder

        intp_t = h.type

        mask = self.mask
        dtype = self._ty.dtype
        tyctx = context.typing_context
        fnty = tyctx.resolve_value_type(operator.eq)
        sig = fnty.get_call_type(tyctx, (dtype, dtype), {})
        eqfn = context.get_function(fnty, sig)

        one = ir.Constant(intp_t, 1)
        five = ir.Constant(intp_t, 5)

        # The perturbation value for probing
        perturb = cgutils.alloca_once_value(builder, h)
        # The index of the entry being considered: start with (hash & mask)
        index = cgutils.alloca_once_value(builder,
                                          builder.and_(h, mask))
        if for_insert:
            # The index of the first deleted entry in the lookup chain
            free_index_sentinel = mask.type(-1)  # highest unsigned index
            free_index = cgutils.alloca_once_value(builder, free_index_sentinel)

        bb_body = builder.append_basic_block("lookup.body")
        bb_found = builder.append_basic_block("lookup.found")
        bb_not_found = builder.append_basic_block("lookup.not_found")
        bb_end = builder.append_basic_block("lookup.end")

        def check_entry(i):
            """
            Check entry *i* against the value being searched for.
            """
            entry = self.get_entry(i)
            entry_hash = entry.hash

            with builder.if_then(builder.icmp_unsigned('==', h, entry_hash)):
                # Hashes are equal, compare values
                # (note this also ensures the entry is used)
                eq = eqfn(builder, (item, entry.key))
                with builder.if_then(eq):
                    builder.branch(bb_found)

            with builder.if_then(is_hash_empty(context, builder, entry_hash)):
                builder.branch(bb_not_found)

            if for_insert:
                # Memorize the index of the first deleted entry
                with builder.if_then(is_hash_deleted(context, builder, entry_hash)):
                    j = builder.load(free_index)
                    j = builder.select(builder.icmp_unsigned('==', j, free_index_sentinel),
                                       i, j)
                    builder.store(j, free_index)

        # First linear probing.  When the number of collisions is small,
        # the lineary probing loop achieves better cache locality and
        # is also slightly cheaper computationally.
        with cgutils.for_range(builder, ir.Constant(intp_t, LINEAR_PROBES)):
            i = builder.load(index)
            check_entry(i)
            i = builder.add(i, one)
            i = builder.and_(i, mask)
            builder.store(i, index)

        # If not found after linear probing, switch to a non-linear
        # perturbation keyed on the unmasked hash value.
        # XXX how to tell LLVM this branch is unlikely?
        builder.branch(bb_body)
        with builder.goto_block(bb_body):
            i = builder.load(index)
            check_entry(i)

            # Perturb to go to next entry:
            #   perturb >>= 5
            #   i = (i * 5 + 1 + perturb) & mask
            p = builder.load(perturb)
            p = builder.lshr(p, five)
            i = builder.add(one, builder.mul(i, five))
            i = builder.and_(mask, builder.add(i, p))
            builder.store(i, index)
            builder.store(p, perturb)
            # Loop
            builder.branch(bb_body)

        with builder.goto_block(bb_not_found):
            if for_insert:
                # Not found => for insertion, return the index of the first
                # deleted entry (if any), to avoid creating an infinite
                # lookup chain (issue #1913).
                i = builder.load(index)
                j = builder.load(free_index)
                i = builder.select(builder.icmp_unsigned('==', j, free_index_sentinel),
                                   i, j)
                builder.store(i, index)
            builder.branch(bb_end)

        with builder.goto_block(bb_found):
            builder.branch(bb_end)

        builder.position_at_end(bb_end)

        found = builder.phi(ir.IntType(1), 'found')
        found.add_incoming(cgutils.true_bit, bb_found)
        found.add_incoming(cgutils.false_bit, bb_not_found)

        return found, builder.load(index)

    @contextlib.contextmanager
    def _iterate(self, start=None):
        """
        Iterate over the payload's entries.  Yield a SetLoop.
        """
        context = self._context
        builder = self._builder

        intp_t = context.get_value_type(types.intp)
        one = ir.Constant(intp_t, 1)
        size = builder.add(self.mask, one)

        with cgutils.for_range(builder, size, start=start) as range_loop:
            entry = self.get_entry(range_loop.index)
            is_used = is_hash_used(context, builder, entry.hash)
            with builder.if_then(is_used):
                loop = SetLoop(index=range_loop.index, entry=entry,
                               do_break=range_loop.do_break)
                yield loop

    @contextlib.contextmanager
    def _next_entry(self):
        """
        Yield a random entry from the payload.  Caller must ensure the
        set isn't empty, otherwise the function won't end.
        """
        context = self._context
        builder = self._builder

        intp_t = context.get_value_type(types.intp)
        zero = ir.Constant(intp_t, 0)
        one = ir.Constant(intp_t, 1)
        mask = self.mask

        # Start walking the entries from the stored "search finger" and
        # break as soon as we find a used entry.

        bb_body = builder.append_basic_block('next_entry_body')
        bb_end = builder.append_basic_block('next_entry_end')

        index = cgutils.alloca_once_value(builder, self.finger)
        builder.branch(bb_body)

        with builder.goto_block(bb_body):
            i = builder.load(index)
            # ANDing with mask ensures we stay inside the table boundaries
            i = builder.and_(mask, builder.add(i, one))
            builder.store(i, index)
            entry = self.get_entry(i)
            is_used = is_hash_used(context, builder, entry.hash)
            builder.cbranch(is_used, bb_end, bb_body)

        builder.position_at_end(bb_end)

        # Update the search finger with the next position.  This avoids
        # O(n**2) behaviour when pop() is called in a loop.
        i = builder.load(index)
        self.finger = i
        yield self.get_entry(i)


class SetInstance(object):

    def __init__(self, context, builder, set_type, set_val):
        self._context = context
        self._builder = builder
        self._ty = set_type
        self._entrysize = get_entry_size(context, set_type)
        self._set = context.make_helper(builder, set_type, set_val)

    @property
    def dtype(self):
        return self._ty.dtype

    @property
    def payload(self):
        """
        The _SetPayload for this set.
        """
        # This cannot be cached as the pointer can move around!
        context = self._context
        builder = self._builder

        ptr = self._context.nrt.meminfo_data(builder, self.meminfo)
        return _SetPayload(context, builder, self._ty, ptr)

    @property
    def value(self):
        return self._set._getvalue()

    @property
    def meminfo(self):
        return self._set.meminfo

    @property
    def parent(self):
        return self._set.parent

    @parent.setter
    def parent(self, value):
        self._set.parent = value

    def get_size(self):
        """
        Return the number of elements in the size.
        """
        return self.payload.used

    def set_dirty(self, val):
        if self._ty.reflected:
            self.payload.dirty = cgutils.true_bit if val else cgutils.false_bit

    def _add_entry(self, payload, entry, item, h, do_resize=True):
        context = self._context
        builder = self._builder

        old_hash = entry.hash
        entry.hash = h
        self.incref_value(item)
        entry.key = item
        # used++
        used = payload.used
        one = ir.Constant(used.type, 1)
        used = payload.used = builder.add(used, one)
        # fill++ if entry wasn't a deleted one
        with builder.if_then(is_hash_empty(context, builder, old_hash),
                             likely=True):
            payload.fill = builder.add(payload.fill, one)
        # Grow table if necessary
        if do_resize:
            self.upsize(used)
        self.set_dirty(True)

    def _add_key(self, payload, item, h, do_resize=True, do_incref=True):
        context = self._context
        builder = self._builder

        found, i = payload._lookup(item, h, for_insert=True)
        not_found = builder.not_(found)

        with builder.if_then(not_found):
            # Not found => add it
            entry = payload.get_entry(i)
            old_hash = entry.hash
            entry.hash = h
            if do_incref:
                self.incref_value(item)
            entry.key = item
            # used++
            used = payload.used
            one = ir.Constant(used.type, 1)
            used = payload.used = builder.add(used, one)
            # fill++ if entry wasn't a deleted one
            with builder.if_then(is_hash_empty(context, builder, old_hash),
                                 likely=True):
                payload.fill = builder.add(payload.fill, one)
            # Grow table if necessary
            if do_resize:
                self.upsize(used)
            self.set_dirty(True)

    def _remove_entry(self, payload, entry, do_resize=True, do_decref=True):
        # Mark entry deleted
        entry.hash = ir.Constant(entry.hash.type, DELETED)
        if do_decref:
            self.decref_value(entry.key)
        # used--
        used = payload.used
        one = ir.Constant(used.type, 1)
        used = payload.used = self._builder.sub(used, one)
        # Shrink table if necessary
        if do_resize:
            self.downsize(used)
        self.set_dirty(True)

    def _remove_key(self, payload, item, h, do_resize=True):
        context = self._context
        builder = self._builder

        found, i = payload._lookup(item, h)

        with builder.if_then(found):
            entry = payload.get_entry(i)
            self._remove_entry(payload, entry, do_resize)

        return found

    def add(self, item, do_resize=True):
        context = self._context
        builder = self._builder

        payload = self.payload
        h = get_hash_value(context, builder, self._ty.dtype, item)
        self._add_key(payload, item, h, do_resize)

    def add_pyapi(self, pyapi, item, do_resize=True):
        """A version of .add for use inside functions following Python calling
        convention.
        """
        context = self._context
        builder = self._builder

        payload = self.payload
        h = self._pyapi_get_hash_value(pyapi, context, builder, item)
        self._add_key(payload, item, h, do_resize)

    def _pyapi_get_hash_value(self, pyapi, context, builder, item):
        """Python API compatible version of `get_hash_value()`.
        """
        argtypes = [self._ty.dtype]
        resty = types.intp

        def wrapper(val):
            return _get_hash_value_intrinsic(val)

        args = [item]
        sig = typing.signature(resty, *argtypes)
        is_error, retval = pyapi.call_jit_code(wrapper, sig, args)
        # Handle return status
        with builder.if_then(is_error, likely=False):
            # Raise nopython exception as a Python exception
            builder.ret(pyapi.get_null_object())
        return retval

    def contains(self, item):
        context = self._context
        builder = self._builder

        payload = self.payload
        h = get_hash_value(context, builder, self._ty.dtype, item)
        found, i = payload._lookup(item, h)
        return found

    def discard(self, item):
        context = self._context
        builder = self._builder

        payload = self.payload
        h = get_hash_value(context, builder, self._ty.dtype, item)
        found = self._remove_key(payload, item, h)
        return found

    def pop(self):
        context = self._context
        builder = self._builder

        lty = context.get_value_type(self._ty.dtype)
        key = cgutils.alloca_once(builder, lty)

        payload = self.payload
        with payload._next_entry() as entry:
            builder.store(entry.key, key)
            # since the value is returned don't decref in _remove_entry()
            self._remove_entry(payload, entry, do_decref=False)

        return builder.load(key)

    def clear(self):
        context = self._context
        builder = self._builder

        intp_t = context.get_value_type(types.intp)
        minsize = ir.Constant(intp_t, MINSIZE)
        self._replace_payload(minsize)
        self.set_dirty(True)

    def copy(self):
        """
        Return a copy of this set.
        """
        context = self._context
        builder = self._builder

        payload = self.payload
        used = payload.used
        fill = payload.fill

        other = type(self)(context, builder, self._ty, None)

        no_deleted_entries = builder.icmp_unsigned('==', used, fill)
        with builder.if_else(no_deleted_entries, likely=True) \
            as (if_no_deleted, if_deleted):
            with if_no_deleted:
                # No deleted entries => raw copy the payload
                ok = other._copy_payload(payload)
                with builder.if_then(builder.not_(ok), likely=False):
                    context.call_conv.return_user_exc(builder, MemoryError,
                                                      ("cannot copy set",))

            with if_deleted:
                # Deleted entries => re-insert entries one by one
                nentries = self.choose_alloc_size(context, builder, used)
                ok = other._allocate_payload(nentries)
                with builder.if_then(builder.not_(ok), likely=False):
                    context.call_conv.return_user_exc(builder, MemoryError,
                                                      ("cannot copy set",))

                other_payload = other.payload
                with payload._iterate() as loop:
                    entry = loop.entry
                    other._add_key(other_payload, entry.key, entry.hash,
                                   do_resize=False)

        return other

    def intersect(self, other):
        """
        In-place intersection with *other* set.
        """
        context = self._context
        builder = self._builder
        payload = self.payload
        other_payload = other.payload

        with payload._iterate() as loop:
            entry = loop.entry
            found, _ = other_payload._lookup(entry.key, entry.hash)
            with builder.if_then(builder.not_(found)):
                self._remove_entry(payload, entry, do_resize=False)

        # Final downsize
        self.downsize(payload.used)

    def difference(self, other):
        """
        In-place difference with *other* set.
        """
        context = self._context
        builder = self._builder
        payload = self.payload
        other_payload = other.payload

        with other_payload._iterate() as loop:
            entry = loop.entry
            self._remove_key(payload, entry.key, entry.hash, do_resize=False)

        # Final downsize
        self.downsize(payload.used)

    def symmetric_difference(self, other):
        """
        In-place symmetric difference with *other* set.
        """
        context = self._context
        builder = self._builder
        other_payload = other.payload

        with other_payload._iterate() as loop:
            key = loop.entry.key
            h = loop.entry.hash
            # We must reload our payload as it may be resized during the loop
            payload = self.payload
            found, i = payload._lookup(key, h, for_insert=True)
            entry = payload.get_entry(i)
            with builder.if_else(found) as (if_common, if_not_common):
                with if_common:
                    self._remove_entry(payload, entry, do_resize=False)
                with if_not_common:
                    self._add_entry(payload, entry, key, h)

        # Final downsize
        self.downsize(self.payload.used)

    def issubset(self, other, strict=False):
        context = self._context
        builder = self._builder
        payload = self.payload
        other_payload = other.payload

        cmp_op = '<' if strict else '<='

        res = cgutils.alloca_once_value(builder, cgutils.true_bit)
        with builder.if_else(
            builder.icmp_unsigned(cmp_op, payload.used, other_payload.used)
            ) as (if_smaller, if_larger):
            with if_larger:
                # self larger than other => self cannot possibly a subset
                builder.store(cgutils.false_bit, res)
            with if_smaller:
                # check whether each key of self is in other
                with payload._iterate() as loop:
                    entry = loop.entry
                    found, _ = other_payload._lookup(entry.key, entry.hash)
                    with builder.if_then(builder.not_(found)):
                        builder.store(cgutils.false_bit, res)
                        loop.do_break()

        return builder.load(res)

    def isdisjoint(self, other):
        context = self._context
        builder = self._builder
        payload = self.payload
        other_payload = other.payload

        res = cgutils.alloca_once_value(builder, cgutils.true_bit)

        def check(smaller, larger):
            # Loop over the smaller of the two, and search in the larger
            with smaller._iterate() as loop:
                entry = loop.entry
                found, _ = larger._lookup(entry.key, entry.hash)
                with builder.if_then(found):
                    builder.store(cgutils.false_bit, res)
                    loop.do_break()

        with builder.if_else(
            builder.icmp_unsigned('>', payload.used, other_payload.used)
            ) as (if_larger, otherwise):

            with if_larger:
                # len(self) > len(other)
                check(other_payload, payload)

            with otherwise:
                # len(self) <= len(other)
                check(payload, other_payload)

        return builder.load(res)

    def equals(self, other):
        context = self._context
        builder = self._builder
        payload = self.payload
        other_payload = other.payload

        res = cgutils.alloca_once_value(builder, cgutils.true_bit)
        with builder.if_else(
            builder.icmp_unsigned('==', payload.used, other_payload.used)
            ) as (if_same_size, otherwise):
            with if_same_size:
                # same sizes => check whether each key of self is in other
                with payload._iterate() as loop:
                    entry = loop.entry
                    found, _ = other_payload._lookup(entry.key, entry.hash)
                    with builder.if_then(builder.not_(found)):
                        builder.store(cgutils.false_bit, res)
                        loop.do_break()
            with otherwise:
                # different sizes => cannot possibly be equal
                builder.store(cgutils.false_bit, res)

        return builder.load(res)

    @classmethod
    def allocate_ex(cls, context, builder, set_type, nitems=None):
        """
        Allocate a SetInstance with its storage.
        Return a (ok, instance) tuple where *ok* is a LLVM boolean and
        *instance* is a SetInstance object (the object's contents are
        only valid when *ok* is true).
        """
        intp_t = context.get_value_type(types.intp)

        if nitems is None:
            nentries = ir.Constant(intp_t, MINSIZE)
        else:
            if isinstance(nitems, int):
                nitems = ir.Constant(intp_t, nitems)
            nentries = cls.choose_alloc_size(context, builder, nitems)

        self = cls(context, builder, set_type, None)
        ok = self._allocate_payload(nentries)
        return ok, self

    @classmethod
    def allocate(cls, context, builder, set_type, nitems=None):
        """
        Allocate a SetInstance with its storage.  Same as allocate_ex(),
        but return an initialized *instance*.  If allocation failed,
        control is transferred to the caller using the target's current
        call convention.
        """
        ok, self = cls.allocate_ex(context, builder, set_type, nitems)
        with builder.if_then(builder.not_(ok), likely=False):
            context.call_conv.return_user_exc(builder, MemoryError,
                                              ("cannot allocate set",))
        return self

    @classmethod
    def from_meminfo(cls, context, builder, set_type, meminfo):
        """
        Allocate a new set instance pointing to an existing payload
        (a meminfo pointer).
        Note the parent field has to be filled by the caller.
        """
        self = cls(context, builder, set_type, None)
        self._set.meminfo = meminfo
        self._set.parent = context.get_constant_null(types.pyobject)
        context.nrt.incref(builder, set_type, self.value)
        # Payload is part of the meminfo, no need to touch it
        return self

    @classmethod
    def choose_alloc_size(cls, context, builder, nitems):
        """
        Choose a suitable number of entries for the given number of items.
        """
        intp_t = nitems.type
        one = ir.Constant(intp_t, 1)
        minsize = ir.Constant(intp_t, MINSIZE)

        # Ensure number of entries >= 2 * used
        min_entries = builder.shl(nitems, one)
        # Find out first suitable power of 2, starting from MINSIZE
        size_p = cgutils.alloca_once_value(builder, minsize)

        bb_body = builder.append_basic_block("calcsize.body")
        bb_end = builder.append_basic_block("calcsize.end")

        builder.branch(bb_body)

        with builder.goto_block(bb_body):
            size = builder.load(size_p)
            is_large_enough = builder.icmp_unsigned('>=', size, min_entries)
            with builder.if_then(is_large_enough, likely=False):
                builder.branch(bb_end)
            next_size = builder.shl(size, one)
            builder.store(next_size, size_p)
            builder.branch(bb_body)

        builder.position_at_end(bb_end)
        return builder.load(size_p)

    def upsize(self, nitems):
        """
        When adding to the set, ensure it is properly sized for the given
        number of used entries.
        """
        context = self._context
        builder = self._builder
        intp_t = nitems.type

        one = ir.Constant(intp_t, 1)
        two = ir.Constant(intp_t, 2)

        payload = self.payload

        # Ensure number of entries >= 2 * used
        min_entries = builder.shl(nitems, one)
        size = builder.add(payload.mask, one)
        need_resize = builder.icmp_unsigned('>=', min_entries, size)

        with builder.if_then(need_resize, likely=False):
            # Find out next suitable size
            new_size_p = cgutils.alloca_once_value(builder, size)

            bb_body = builder.append_basic_block("calcsize.body")
            bb_end = builder.append_basic_block("calcsize.end")

            builder.branch(bb_body)

            with builder.goto_block(bb_body):
                # Multiply by 4 (ensuring size remains a power of two)
                new_size = builder.load(new_size_p)
                new_size = builder.shl(new_size, two)
                builder.store(new_size, new_size_p)
                is_too_small = builder.icmp_unsigned('>=', min_entries, new_size)
                builder.cbranch(is_too_small, bb_body, bb_end)

            builder.position_at_end(bb_end)

            new_size = builder.load(new_size_p)
            if DEBUG_ALLOCS:
                context.printf(builder,
                               "upsize to %zd items: current size = %zd, "
                               "min entries = %zd, new size = %zd\n",
                               nitems, size, min_entries, new_size)
            self._resize(payload, new_size, "cannot grow set")

    def downsize(self, nitems):
        """
        When removing from the set, ensure it is properly sized for the given
        number of used entries.
        """
        context = self._context
        builder = self._builder
        intp_t = nitems.type

        one = ir.Constant(intp_t, 1)
        two = ir.Constant(intp_t, 2)
        minsize = ir.Constant(intp_t, MINSIZE)

        payload = self.payload

        # Ensure entries >= max(2 * used, MINSIZE)
        min_entries = builder.shl(nitems, one)
        min_entries = builder.select(builder.icmp_unsigned('>=', min_entries, minsize),
                                     min_entries, minsize)
        # Shrink only if size >= 4 * min_entries && size > MINSIZE
        max_size = builder.shl(min_entries, two)
        size = builder.add(payload.mask, one)
        need_resize = builder.and_(
            builder.icmp_unsigned('<=', max_size, size),
            builder.icmp_unsigned('<', minsize, size))

        with builder.if_then(need_resize, likely=False):
            # Find out next suitable size
            new_size_p = cgutils.alloca_once_value(builder, size)

            bb_body = builder.append_basic_block("calcsize.body")
            bb_end = builder.append_basic_block("calcsize.end")

            builder.branch(bb_body)

            with builder.goto_block(bb_body):
                # Divide by 2 (ensuring size remains a power of two)
                new_size = builder.load(new_size_p)
                new_size = builder.lshr(new_size, one)
                # Keep current size if new size would be < min_entries
                is_too_small = builder.icmp_unsigned('>', min_entries, new_size)
                with builder.if_then(is_too_small):
                    builder.branch(bb_end)
                builder.store(new_size, new_size_p)
                builder.branch(bb_body)

            builder.position_at_end(bb_end)

            # Ensure new_size >= MINSIZE
            new_size = builder.load(new_size_p)
            # At this point, new_size should be < size if the factors
            # above were chosen carefully!

            if DEBUG_ALLOCS:
                context.printf(builder,
                               "downsize to %zd items: current size = %zd, "
                               "min entries = %zd, new size = %zd\n",
                               nitems, size, min_entries, new_size)
            self._resize(payload, new_size, "cannot shrink set")

    def _resize(self, payload, nentries, errmsg):
        """
        Resize the payload to the given number of entries.

        CAUTION: *nentries* must be a power of 2!
        """
        context = self._context
        builder = self._builder

        # Allocate new entries
        old_payload = payload

        ok = self._allocate_payload(nentries, realloc=True)
        with builder.if_then(builder.not_(ok), likely=False):
            context.call_conv.return_user_exc(builder, MemoryError,
                                              (errmsg,))

        # Re-insert old entries
        # No incref since they already were the first time they were inserted
        payload = self.payload
        with old_payload._iterate() as loop:
            entry = loop.entry
            self._add_key(payload, entry.key, entry.hash,
                          do_resize=False, do_incref=False)

        self._free_payload(old_payload.ptr)

    def _replace_payload(self, nentries):
        """
        Replace the payload with a new empty payload with the given number
        of entries.

        CAUTION: *nentries* must be a power of 2!
        """
        context = self._context
        builder = self._builder

        # decref all of the previous entries
        with self.payload._iterate() as loop:
            entry = loop.entry
            self.decref_value(entry.key)

        # Free old payload
        self._free_payload(self.payload.ptr)

        ok = self._allocate_payload(nentries, realloc=True)
        with builder.if_then(builder.not_(ok), likely=False):
            context.call_conv.return_user_exc(builder, MemoryError,
                                              ("cannot reallocate set",))

    def _allocate_payload(self, nentries, realloc=False):
        """
        Allocate and initialize payload for the given number of entries.
        If *realloc* is True, the existing meminfo is reused.

        CAUTION: *nentries* must be a power of 2!
        """
        context = self._context
        builder = self._builder

        ok = cgutils.alloca_once_value(builder, cgutils.true_bit)

        intp_t = context.get_value_type(types.intp)
        zero = ir.Constant(intp_t, 0)
        one = ir.Constant(intp_t, 1)

        payload_type = context.get_data_type(types.SetPayload(self._ty))
        payload_size = context.get_abi_sizeof(payload_type)
        entry_size = self._entrysize
        # Account for the fact that the payload struct already contains an entry
        payload_size -= entry_size

        # Total allocation size = <payload header size> + nentries * entry_size
        allocsize, ovf = cgutils.muladd_with_overflow(builder, nentries,
                                                      ir.Constant(intp_t, entry_size),
                                                      ir.Constant(intp_t, payload_size))
        with builder.if_then(ovf, likely=False):
            builder.store(cgutils.false_bit, ok)

        with builder.if_then(builder.load(ok), likely=True):
            if realloc:
                meminfo = self._set.meminfo
                ptr = context.nrt.meminfo_varsize_alloc_unchecked(builder,
                                                                  meminfo,
                                                        size=allocsize)
                alloc_ok = cgutils.is_null(builder, ptr)
            else:
                # create destructor to be called upon set destruction
                dtor = self._imp_dtor(context, builder.module)
                meminfo = context.nrt.meminfo_new_varsize_dtor_unchecked(
                    builder, allocsize, builder.bitcast(dtor, cgutils.voidptr_t))
                alloc_ok = cgutils.is_null(builder, meminfo)

            with builder.if_else(alloc_ok,
                                 likely=False) as (if_error, if_ok):
                with if_error:
                    builder.store(cgutils.false_bit, ok)
                with if_ok:
                    if not realloc:
                        self._set.meminfo = meminfo
                        self._set.parent = context.get_constant_null(types.pyobject)
                    payload = self.payload
                    # Initialize entries to 0xff (EMPTY)
                    cgutils.memset(builder, payload.ptr, allocsize, 0xFF)
                    payload.used = zero
                    payload.fill = zero
                    payload.finger = zero
                    new_mask = builder.sub(nentries, one)
                    payload.mask = new_mask

                    if DEBUG_ALLOCS:
                        context.printf(builder,
                                       "allocated %zd bytes for set at %p: mask = %zd\n",
                                       allocsize, payload.ptr, new_mask)

        return builder.load(ok)

    def _free_payload(self, ptr):
        """
        Free an allocated old payload at *ptr*.
        """
        self._context.nrt.meminfo_varsize_free(self._builder, self.meminfo, ptr)

    def _copy_payload(self, src_payload):
        """
        Raw-copy the given payload into self.
        """
        context = self._context
        builder = self._builder

        ok = cgutils.alloca_once_value(builder, cgutils.true_bit)

        intp_t = context.get_value_type(types.intp)
        zero = ir.Constant(intp_t, 0)
        one = ir.Constant(intp_t, 1)

        payload_type = context.get_data_type(types.SetPayload(self._ty))
        payload_size = context.get_abi_sizeof(payload_type)
        entry_size = self._entrysize
        # Account for the fact that the payload struct already contains an entry
        payload_size -= entry_size

        mask = src_payload.mask
        nentries = builder.add(one, mask)

        # Total allocation size = <payload header size> + nentries * entry_size
        # (note there can't be any overflow since we're reusing an existing
        #  payload's parameters)
        allocsize = builder.add(ir.Constant(intp_t, payload_size),
                                builder.mul(ir.Constant(intp_t, entry_size),
                                            nentries))

        with builder.if_then(builder.load(ok), likely=True):
            # create destructor for new meminfo
            dtor = self._imp_dtor(context, builder.module)
            meminfo = context.nrt.meminfo_new_varsize_dtor_unchecked(
                builder, allocsize, builder.bitcast(dtor, cgutils.voidptr_t))
            alloc_ok = cgutils.is_null(builder, meminfo)

            with builder.if_else(alloc_ok, likely=False) as (if_error, if_ok):
                with if_error:
                    builder.store(cgutils.false_bit, ok)
                with if_ok:
                    self._set.meminfo = meminfo
                    payload = self.payload
                    payload.used = src_payload.used
                    payload.fill = src_payload.fill
                    payload.finger = zero
                    payload.mask = mask

                    # instead of using `_add_key` for every entry, since the
                    # size of the new set is the same, we can just copy the
                    # data directly without having to re-compute the hash
                    cgutils.raw_memcpy(builder, payload.entries,
                                       src_payload.entries, nentries,
                                       entry_size)
                    # increment the refcounts to simulate `_add_key` for each
                    # element
                    with src_payload._iterate() as loop:
                        self.incref_value(loop.entry.key)

                    if DEBUG_ALLOCS:
                        context.printf(builder,
                                       "allocated %zd bytes for set at %p: mask = %zd\n",
                                       allocsize, payload.ptr, mask)

        return builder.load(ok)

    def _imp_dtor(self, context, module):
        """Define the dtor for set
        """
        llvoidptr = cgutils.voidptr_t
        llsize_t= context.get_value_type(types.size_t)
        # create a dtor function that takes (void* set, size_t size, void* dtor_info)
        fnty = ir.FunctionType(
            ir.VoidType(),
            [llvoidptr, llsize_t, llvoidptr],
        )
        # create type-specific name
        fname = f".dtor.set.{self._ty.dtype}"

        fn = cgutils.get_or_insert_function(module, fnty, name=fname)

        if fn.is_declaration:
            # Set linkage
            fn.linkage = 'linkonce_odr'
            # Define
            builder = ir.IRBuilder(fn.append_basic_block())
            payload = _SetPayload(context, builder, self._ty, fn.args[0])
            with payload._iterate() as loop:
                entry = loop.entry
                context.nrt.decref(builder, self._ty.dtype, entry.key)
            builder.ret_void()

        return fn

    def incref_value(self, val):
        """Incref an element value
        """
        self._context.nrt.incref(self._builder, self._ty.dtype, val)

    def decref_value(self, val):
        """Decref an element value
        """
        self._context.nrt.decref(self._builder, self._ty.dtype, val)


class SetIterInstance(object):

    def __init__(self, context, builder, iter_type, iter_val):
        self._context = context
        self._builder = builder
        self._ty = iter_type
        self._iter = context.make_helper(builder, iter_type, iter_val)
        ptr = self._context.nrt.meminfo_data(builder, self.meminfo)
        self._payload = _SetPayload(context, builder, self._ty.container, ptr)

    @classmethod
    def from_set(cls, context, builder, iter_type, set_val):
        set_inst = SetInstance(context, builder, iter_type.container, set_val)
        self = cls(context, builder, iter_type, None)
        index = context.get_constant(types.intp, 0)
        self._iter.index = cgutils.alloca_once_value(builder, index)
        self._iter.meminfo = set_inst.meminfo
        return self

    @property
    def value(self):
        return self._iter._getvalue()

    @property
    def meminfo(self):
        return self._iter.meminfo

    @property
    def index(self):
        return self._builder.load(self._iter.index)

    @index.setter
    def index(self, value):
        self._builder.store(value, self._iter.index)

    def iternext(self, result):
        index = self.index
        payload = self._payload
        one = ir.Constant(index.type, 1)

        result.set_exhausted()

        with payload._iterate(start=index) as loop:
            # An entry was found
            entry = loop.entry
            result.set_valid()
            result.yield_(entry.key)
            self.index = self._builder.add(loop.index, one)
            loop.do_break()


#-------------------------------------------------------------------------------
# Constructors

def build_set(context, builder, set_type, items):
    """
    Build a set of the given type, containing the given items.
    """
    nitems = len(items)
    inst = SetInstance.allocate(context, builder, set_type, nitems)

    if nitems > 0:

        # Populate set.  Inlining the insertion code for each item would be very
        # costly, instead we create a LLVM array and iterate over it.
        array = cgutils.pack_array(builder, items)
        array_ptr = cgutils.alloca_once_value(builder, array)

        count = context.get_constant(types.intp, nitems)
        with cgutils.for_range(builder, count) as loop:
            item = builder.load(cgutils.gep(builder, array_ptr, 0, loop.index))
            inst.add(item)

    return impl_ret_new_ref(context, builder, set_type, inst.value)


@lower_builtin(set)
def set_empty_constructor(context, builder, sig, args):
    set_type = sig.return_type
    inst = SetInstance.allocate(context, builder, set_type)
    return impl_ret_new_ref(context, builder, set_type, inst.value)

@lower_builtin(set, types.IterableType)
def set_constructor(context, builder, sig, args):
    set_type = sig.return_type
    items_type, = sig.args
    items, = args

    # If the argument has a len(), preallocate the set so as to
    # avoid resizes.
    # `for_iter` increfs each item in the set, so a `decref` is required each
    # iteration to balance. Because the `incref` from `.add` is dependent on
    # the item not already existing in the set, just removing its incref is not
    # enough to guarantee all memory is freed
    n = call_len(context, builder, items_type, items)
    inst = SetInstance.allocate(context, builder, set_type, n)
    with for_iter(context, builder, items_type, items) as loop:
        inst.add(loop.value)
        context.nrt.decref(builder, set_type.dtype, loop.value)

    return impl_ret_new_ref(context, builder, set_type, inst.value)


#-------------------------------------------------------------------------------
# Various operations

@lower_builtin(len, types.Set)
def set_len(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    return inst.get_size()

@lower_builtin(operator.contains, types.Set, types.Any)
def in_set(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    return inst.contains(args[1])

@lower_builtin('getiter', types.Set)
def getiter_set(context, builder, sig, args):
    inst = SetIterInstance.from_set(context, builder, sig.return_type, args[0])
    return impl_ret_borrowed(context, builder, sig.return_type, inst.value)

@lower_builtin('iternext', types.SetIter)
@iternext_impl(RefType.BORROWED)
def iternext_listiter(context, builder, sig, args, result):
    inst = SetIterInstance(context, builder, sig.args[0], args[0])
    inst.iternext(result)


#-------------------------------------------------------------------------------
# Methods

# One-item-at-a-time operations

@lower_builtin("set.add", types.Set, types.Any)
def set_add(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    item = args[1]
    inst.add(item)

    return context.get_dummy_value()


@intrinsic
def _set_discard(typingctx, s, item):
    sig = types.none(s, item)

    def set_discard(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        item = args[1]
        inst.discard(item)

        return context.get_dummy_value()

    return sig, set_discard


@overload_method(types.Set, "discard")
def ol_set_discard(s, item):
    return lambda s, item: _set_discard(s, item)


@intrinsic
def _set_pop(typingctx, s):
    sig = s.dtype(s)

    def set_pop(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        used = inst.payload.used
        with builder.if_then(cgutils.is_null(builder, used), likely=False):
            context.call_conv.return_user_exc(builder, KeyError,
                                            ("set.pop(): empty set",))

        return inst.pop()

    return sig, set_pop


@overload_method(types.Set, "pop")
def ol_set_pop(s):
    return lambda s: _set_pop(s)


@intrinsic
def _set_remove(typingctx, s, item):
    sig = types.none(s, item)

    def set_remove(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        item = args[1]
        found = inst.discard(item)
        with builder.if_then(builder.not_(found), likely=False):
            context.call_conv.return_user_exc(builder, KeyError,
                                            ("set.remove(): key not in set",))

        return context.get_dummy_value()

    return sig, set_remove


@overload_method(types.Set, "remove")
def ol_set_remove(s, item):
    if s.dtype == item:
        return lambda s, item: _set_remove(s, item)


# Mutating set operations

@intrinsic
def _set_clear(typingctx, s):
    sig = types.none(s)

    def set_clear(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        inst.clear()
        return context.get_dummy_value()

    return sig, set_clear


@overload_method(types.Set, "clear")
def ol_set_clear(s):
    return lambda s: _set_clear(s)


@intrinsic
def _set_copy(typingctx, s):
    sig = s(s)

    def set_copy(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        other = inst.copy()
        return impl_ret_new_ref(context, builder, sig.return_type, other.value)

    return sig, set_copy


@overload_method(types.Set, "copy")
def ol_set_copy(s):
    return lambda s: _set_copy(s)


def set_difference_update(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    other = SetInstance(context, builder, sig.args[1], args[1])

    inst.difference(other)

    return context.get_dummy_value()


@intrinsic
def _set_difference_update(typingctx, a, b):
    sig = types.none(a, b)
    return sig, set_difference_update


@overload_method(types.Set, "difference_update")
def set_difference_update_impl(a, b):
    check_all_set(a, b)
    return lambda a, b: _set_difference_update(a, b)


def set_intersection_update(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    other = SetInstance(context, builder, sig.args[1], args[1])
    inst.intersect(other)
    return context.get_dummy_value()


@intrinsic
def _set_intersection_update(typingctx, a, b):
    sig = types.none(a, b)
    return sig, set_intersection_update


@overload_method(types.Set, "intersection_update")
def set_intersection_update_impl(a, b):
    check_all_set(a, b)
    return lambda a, b: _set_intersection_update(a, b)


def set_symmetric_difference_update(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    other = SetInstance(context, builder, sig.args[1], args[1])
    inst.symmetric_difference(other)
    return context.get_dummy_value()


@intrinsic
def _set_symmetric_difference_update(typingctx, a, b):
    sig = types.none(a, b)
    return sig, set_symmetric_difference_update


@overload_method(types.Set, "symmetric_difference_update")
def set_symmetric_difference_update_impl(a, b):
    check_all_set(a, b)
    return lambda a, b: _set_symmetric_difference_update(a, b)


@lower_builtin("set.update", types.Set, types.IterableType)
def set_update(context, builder, sig, args):
    inst = SetInstance(context, builder, sig.args[0], args[0])
    items_type = sig.args[1]
    items = args[1]

    # If the argument has a len(), assume there are few collisions and
    # presize to len(set) + len(items)
    n = call_len(context, builder, items_type, items)
    if n is not None:
        new_size = builder.add(inst.payload.used, n)
        inst.upsize(new_size)

    with for_iter(context, builder, items_type, items) as loop:
        # make sure that the items being added are of the same dtype as the
        # set instance
        casted = context.cast(builder, loop.value, items_type.dtype, inst.dtype)
        inst.add(casted)
        # decref each item to counter balance the incref from `for_iter`
        # `.add` will conditionally incref when the item does not already exist
        # in the set, therefore removing its incref is not enough to guarantee
        # all memory is freed
        context.nrt.decref(builder, items_type.dtype, loop.value)

    if n is not None:
        # If we pre-grew the set, downsize in case there were many collisions
        inst.downsize(inst.payload.used)

    return context.get_dummy_value()

def gen_operator_impl(op, impl):
    @intrinsic
    def _set_operator_intr(typingctx, a, b):
        sig = a(a, b)
        def codegen(context, builder, sig, args):
            assert sig.return_type == sig.args[0]
            impl(context, builder, sig, args)
            return impl_ret_borrowed(context, builder, sig.args[0], args[0])
        return sig, codegen

    @overload(op)
    def _ol_set_operator(a, b):
        check_all_set(a, b)
        return lambda a, b: _set_operator_intr(a, b)


for op_, op_impl in [
    (operator.iand, set_intersection_update),
    (operator.ior, set_update),
    (operator.isub, set_difference_update),
    (operator.ixor, set_symmetric_difference_update),
    ]:
    gen_operator_impl(op_, op_impl)


# Set operations creating a new set

@overload(operator.sub)
@overload_method(types.Set, "difference")
def impl_set_difference(a, b):
    check_all_set(a, b)

    def difference_impl(a, b):
        s = a.copy()
        s.difference_update(b)
        return s

    return difference_impl

@overload(operator.and_)
@overload_method(types.Set, "intersection")
def set_intersection(a, b):
    check_all_set(a, b)

    def intersection_impl(a, b):
        if len(a) < len(b):
            s = a.copy()
            s.intersection_update(b)
            return s
        else:
            s = b.copy()
            s.intersection_update(a)
            return s

    return intersection_impl

@overload(operator.xor)
@overload_method(types.Set, "symmetric_difference")
def set_symmetric_difference(a, b):
    check_all_set(a, b)

    def symmetric_difference_impl(a, b):
        if len(a) > len(b):
            s = a.copy()
            s.symmetric_difference_update(b)
            return s
        else:
            s = b.copy()
            s.symmetric_difference_update(a)
            return s

    return symmetric_difference_impl

@overload(operator.or_)
@overload_method(types.Set, "union")
def set_union(a, b):
    check_all_set(a, b)

    def union_impl(a, b):
        if len(a) > len(b):
            s = a.copy()
            s.update(b)
            return s
        else:
            s = b.copy()
            s.update(a)
            return s

    return union_impl


# Predicates

@intrinsic
def _set_isdisjoint(typingctx, a, b):
    sig = types.boolean(a, b)

    def codegen(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        other = SetInstance(context, builder, sig.args[1], args[1])

        return inst.isdisjoint(other)

    return sig, codegen


@overload_method(types.Set, "isdisjoint")
def set_isdisjoint(a, b):
    check_all_set(a, b)

    return lambda a, b: _set_isdisjoint(a, b)


@intrinsic
def _set_issubset(typingctx, a, b):
    sig = types.boolean(a, b)

    def codegen(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        other = SetInstance(context, builder, sig.args[1], args[1])

        return inst.issubset(other)

    return sig, codegen

@overload(operator.le)
@overload_method(types.Set, "issubset")
def set_issubset(a, b):
    check_all_set(a, b)

    return lambda a, b: _set_issubset(a, b)


@overload(operator.ge)
@overload_method(types.Set, "issuperset")
def set_issuperset(a, b):
    check_all_set(a, b)

    def superset_impl(a, b):
        return b.issubset(a)

    return superset_impl

@intrinsic
def _set_eq(typingctx, a, b):
    sig = types.boolean(a, b)

    def codegen(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        other = SetInstance(context, builder, sig.args[1], args[1])

        return inst.equals(other)

    return sig, codegen

@overload(operator.eq)
def set_eq(a, b):
    check_all_set(a, b)

    return lambda a, b: _set_eq(a, b)

@overload(operator.ne)
def set_ne(a, b):
    check_all_set(a, b)

    def ne_impl(a, b):
        return not a == b

    return ne_impl

@intrinsic
def _set_lt(typingctx, a, b):
    sig = types.boolean(a, b)

    def codegen(context, builder, sig, args):
        inst = SetInstance(context, builder, sig.args[0], args[0])
        other = SetInstance(context, builder, sig.args[1], args[1])

        return inst.issubset(other, strict=True)

    return sig, codegen

@overload(operator.lt)
def set_lt(a, b):
    check_all_set(a, b)

    return lambda a, b: _set_lt(a, b)

@overload(operator.gt)
def set_gt(a, b):
    check_all_set(a, b)

    def gt_impl(a, b):
        return b < a

    return gt_impl

@lower_builtin(operator.is_, types.Set, types.Set)
def set_is(context, builder, sig, args):
    a = SetInstance(context, builder, sig.args[0], args[0])
    b = SetInstance(context, builder, sig.args[1], args[1])
    ma = builder.ptrtoint(a.meminfo, cgutils.intp_t)
    mb = builder.ptrtoint(b.meminfo, cgutils.intp_t)
    return builder.icmp_signed('==', ma, mb)


# -----------------------------------------------------------------------------
# Implicit casting

@lower_cast(types.Set, types.Set)
def set_to_set(context, builder, fromty, toty, val):
    # Casting from non-reflected to reflected
    assert fromty.dtype == toty.dtype
    return val
