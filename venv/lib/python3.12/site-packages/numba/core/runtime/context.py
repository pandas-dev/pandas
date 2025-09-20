import functools
from collections import namedtuple

from llvmlite import ir
from numba.core import types, cgutils, errors, config
from numba.core.utils import PYVERSION


_NRT_Meminfo_Functions = namedtuple("_NRT_Meminfo_Functions",
                                    ("alloc",
                                     "alloc_dtor",
                                     "alloc_aligned"))


_NRT_MEMINFO_SAFE_API = _NRT_Meminfo_Functions("NRT_MemInfo_alloc_safe",
                                               "NRT_MemInfo_alloc_dtor_safe",
                                               "NRT_MemInfo_alloc_safe_aligned")


_NRT_MEMINFO_DEFAULT_API = _NRT_Meminfo_Functions("NRT_MemInfo_alloc",
                                                  "NRT_MemInfo_alloc_dtor",
                                                  "NRT_MemInfo_alloc_aligned")


class NRTContext(object):
    """
    An object providing access to NRT APIs in the lowering pass.
    """

    def __init__(self, context, enabled):
        self._context = context
        self._enabled = enabled
        # If DEBUG_NRT is set, use the safe function variants which use memset
        # to inject a few known bytes into the start of allocated regions.
        if config.DEBUG_NRT:
            self._meminfo_api = _NRT_MEMINFO_SAFE_API
        else:
            self._meminfo_api = _NRT_MEMINFO_DEFAULT_API

    def _require_nrt(self):
        if not self._enabled:
            raise errors.NumbaRuntimeError("NRT required but not enabled")

    def _check_null_result(func):
        @functools.wraps(func)
        def wrap(self, builder, *args, **kwargs):
            memptr = func(self, builder, *args, **kwargs)
            msg = "Allocation failed (probably too large)."
            cgutils.guard_memory_error(self._context, builder, memptr, msg=msg)
            return memptr
        return wrap

    @_check_null_result
    def allocate(self, builder, size):
        """
        Low-level allocate a new memory area of `size` bytes. The result of the
        call is checked and if it is NULL, i.e. allocation failed, then a
        MemoryError is raised.
        """
        return self.allocate_unchecked(builder, size)

    def allocate_unchecked(self, builder, size):
        """
        Low-level allocate a new memory area of `size` bytes. Returns NULL to
        indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t, [cgutils.intp_t])
        fn = cgutils.get_or_insert_function(mod, fnty, "NRT_Allocate")
        fn.return_value.add_attribute("noalias")
        return builder.call(fn, [size])

    def free(self, builder, ptr):
        """
        Low-level free a memory area allocated with allocate().
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(ir.VoidType(), [cgutils.voidptr_t])
        fn = cgutils.get_or_insert_function(mod, fnty, "NRT_Free")
        return builder.call(fn, [ptr])

    @_check_null_result
    def meminfo_alloc(self, builder, size):
        """
        Allocate a new MemInfo with a data payload of `size` bytes.

        A pointer to the MemInfo is returned.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_alloc_unchecked(builder, size)

    def meminfo_alloc_unchecked(self, builder, size):
        """
        Allocate a new MemInfo with a data payload of `size` bytes.

        A pointer to the MemInfo is returned.

        Returns NULL to indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t, [cgutils.intp_t])
        fn = cgutils.get_or_insert_function(mod, fnty,
                                            self._meminfo_api.alloc)
        fn.return_value.add_attribute("noalias")
        return builder.call(fn, [size])

    @_check_null_result
    def meminfo_alloc_dtor(self, builder, size, dtor):
        """
        Allocate a new MemInfo with a data payload of `size` bytes and a
        destructor `dtor`.

        A pointer to the MemInfo is returned.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_alloc_dtor_unchecked(builder, size, dtor)

    def meminfo_alloc_dtor_unchecked(self, builder, size, dtor):
        """
        Allocate a new MemInfo with a data payload of `size` bytes and a
        destructor `dtor`.

        A pointer to the MemInfo is returned.

        Returns NULL to indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t,
                               [cgutils.intp_t, cgutils.voidptr_t])
        fn = cgutils.get_or_insert_function(mod, fnty,
                                            self._meminfo_api.alloc_dtor)
        fn.return_value.add_attribute("noalias")
        return builder.call(fn, [size,
                                 builder.bitcast(dtor, cgutils.voidptr_t)])

    @_check_null_result
    def meminfo_alloc_aligned(self, builder, size, align):
        """
        Allocate a new MemInfo with an aligned data payload of `size` bytes.
        The data pointer is aligned to `align` bytes.  `align` can be either
        a Python int or a LLVM uint32 value.

        A pointer to the MemInfo is returned.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_alloc_aligned_unchecked(builder, size, align)

    def meminfo_alloc_aligned_unchecked(self, builder, size, align):
        """
        Allocate a new MemInfo with an aligned data payload of `size` bytes.
        The data pointer is aligned to `align` bytes.  `align` can be either
        a Python int or a LLVM uint32 value.

        A pointer to the MemInfo is returned.

        Returns NULL to indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        u32 = ir.IntType(32)
        fnty = ir.FunctionType(cgutils.voidptr_t, [cgutils.intp_t, u32])
        fn = cgutils.get_or_insert_function(mod, fnty,
                                            self._meminfo_api.alloc_aligned)
        fn.return_value.add_attribute("noalias")
        if isinstance(align, int):
            align = self._context.get_constant(types.uint32, align)
        else:
            assert align.type == u32, "align must be a uint32"
        return builder.call(fn, [size, align])

    @_check_null_result
    def meminfo_new_varsize(self, builder, size):
        """
        Allocate a MemInfo pointing to a variable-sized data area.  The area
        is separately allocated (i.e. two allocations are made) so that
        re-allocating it doesn't change the MemInfo's address.

        A pointer to the MemInfo is returned.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_new_varsize_unchecked(builder, size)

    def meminfo_new_varsize_unchecked(self, builder, size):
        """
        Allocate a MemInfo pointing to a variable-sized data area.  The area
        is separately allocated (i.e. two allocations are made) so that
        re-allocating it doesn't change the MemInfo's address.

        A pointer to the MemInfo is returned.

        Returns NULL to indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t, [cgutils.intp_t])
        fn = cgutils.get_or_insert_function(mod, fnty,
                                            "NRT_MemInfo_new_varsize")
        fn.return_value.add_attribute("noalias")
        return builder.call(fn, [size])

    @_check_null_result
    def meminfo_new_varsize_dtor(self, builder, size, dtor):
        """
        Like meminfo_new_varsize() but also set the destructor for
        cleaning up references to objects inside the allocation.

        A pointer to the MemInfo is returned.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_new_varsize_dtor_unchecked(builder, size, dtor)

    def meminfo_new_varsize_dtor_unchecked(self, builder, size, dtor):
        """
        Like meminfo_new_varsize() but also set the destructor for
        cleaning up references to objects inside the allocation.

        A pointer to the MemInfo is returned.

        Returns NULL to indicate error/failure to allocate.
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t,
                               [cgutils.intp_t, cgutils.voidptr_t])
        fn = cgutils.get_or_insert_function(
            mod, fnty, "NRT_MemInfo_new_varsize_dtor")
        return builder.call(fn, [size, dtor])

    @_check_null_result
    def meminfo_varsize_alloc(self, builder, meminfo, size):
        """
        Allocate a new data area for a MemInfo created by meminfo_new_varsize().
        The new data pointer is returned, for convenience.

        Contrary to realloc(), this always allocates a new area and doesn't
        copy the old data.  This is useful if resizing a container needs
        more than simply copying the data area (e.g. for hash tables).

        The old pointer will have to be freed with meminfo_varsize_free().

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_varsize_alloc_unchecked(builder, meminfo, size)

    def meminfo_varsize_alloc_unchecked(self, builder, meminfo, size):
        """
        Allocate a new data area for a MemInfo created by meminfo_new_varsize().
        The new data pointer is returned, for convenience.

        Contrary to realloc(), this always allocates a new area and doesn't
        copy the old data.  This is useful if resizing a container needs
        more than simply copying the data area (e.g. for hash tables).

        The old pointer will have to be freed with meminfo_varsize_free().

        Returns NULL to indicate error/failure to allocate.
        """
        return self._call_varsize_alloc(builder, meminfo, size,
                                        "NRT_MemInfo_varsize_alloc")

    @_check_null_result
    def meminfo_varsize_realloc(self, builder, meminfo, size):
        """
        Reallocate a data area allocated by meminfo_new_varsize().
        The new data pointer is returned, for convenience.

        The result of the call is checked and if it is NULL, i.e. allocation
        failed, then a MemoryError is raised.
        """
        return self.meminfo_varsize_realloc_unchecked(builder, meminfo, size)

    def meminfo_varsize_realloc_unchecked(self, builder, meminfo, size):
        """
        Reallocate a data area allocated by meminfo_new_varsize().
        The new data pointer is returned, for convenience.

        Returns NULL to indicate error/failure to allocate.
        """
        return self._call_varsize_alloc(builder, meminfo, size,
                                        "NRT_MemInfo_varsize_realloc")

    def meminfo_varsize_free(self, builder, meminfo, ptr):
        """
        Free a memory area allocated for a NRT varsize object.
        Note this does *not* free the NRT object itself!
        """
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(ir.VoidType(),
                               [cgutils.voidptr_t, cgutils.voidptr_t])
        fn = cgutils.get_or_insert_function(mod, fnty,
                                            "NRT_MemInfo_varsize_free")
        return builder.call(fn, (meminfo, ptr))

    def _call_varsize_alloc(self, builder, meminfo, size, funcname):
        self._require_nrt()

        mod = builder.module
        fnty = ir.FunctionType(cgutils.voidptr_t,
                               [cgutils.voidptr_t, cgutils.intp_t])
        fn = cgutils.get_or_insert_function(mod, fnty, funcname)
        fn.return_value.add_attribute("noalias")
        return builder.call(fn, [meminfo, size])

    def meminfo_data(self, builder, meminfo):
        """
        Given a MemInfo pointer, return a pointer to the allocated data
        managed by it.  This works for MemInfos allocated with all the
        above methods.
        """
        self._require_nrt()

        from numba.core.runtime.nrtdynmod import meminfo_data_ty

        mod = builder.module
        fn = cgutils.get_or_insert_function(mod, meminfo_data_ty,
                                            "NRT_MemInfo_data_fast")
        return builder.call(fn, [meminfo])

    def get_meminfos(self, builder, ty, val):
        """Return a list of *(type, meminfo)* inside the given value.
        """
        datamodel = self._context.data_model_manager[ty]
        members = datamodel.traverse(builder)

        meminfos = []
        if datamodel.has_nrt_meminfo():
            mi = datamodel.get_nrt_meminfo(builder, val)
            meminfos.append((ty, mi))

        for mtyp, getter in members:
            field = getter(val)
            inner_meminfos = self.get_meminfos(builder, mtyp, field)
            meminfos.extend(inner_meminfos)
        return meminfos

    def _call_incref_decref(self, builder, typ, value, funcname):
        """Call function of *funcname* on every meminfo found in *value*.
        """
        self._require_nrt()

        from numba.core.runtime.nrtdynmod import incref_decref_ty

        meminfos = self.get_meminfos(builder, typ, value)
        for _, mi in meminfos:
            mod = builder.module
            fn = cgutils.get_or_insert_function(mod, incref_decref_ty,
                                                funcname)
            # XXX "nonnull" causes a crash in test_dyn_array: can this
            # function be called with a NULL pointer?
            fn.args[0].add_attribute("noalias")
            fn.args[0].add_attribute("nocapture")
            builder.call(fn, [mi])

    def incref(self, builder, typ, value):
        """
        Recursively incref the given *value* and its members.
        """
        self._call_incref_decref(builder, typ, value, "NRT_incref")

    def decref(self, builder, typ, value):
        """
        Recursively decref the given *value* and its members.
        """
        self._call_incref_decref(builder, typ, value, "NRT_decref")

    def get_nrt_api(self, builder):
        """Calls NRT_get_api(), which returns the NRT API function table.
        """
        self._require_nrt()

        fnty = ir.FunctionType(cgutils.voidptr_t, ())
        mod = builder.module
        fn = cgutils.get_or_insert_function(mod, fnty, "NRT_get_api")
        return builder.call(fn, ())

    def eh_check(self, builder):
        """Check if an exception is raised
        """
        ctx = self._context
        cc = ctx.call_conv
        # Inspect the excinfo argument on the function
        trystatus = cc.check_try_status(builder)
        excinfo = trystatus.excinfo
        has_raised = builder.not_(cgutils.is_null(builder, excinfo))
        if PYVERSION < (3, 11):
            with builder.if_then(has_raised):
                self.eh_end_try(builder)
        return has_raised

    def eh_try(self, builder):
        """Begin a try-block.
        """
        ctx = self._context
        cc = ctx.call_conv
        cc.set_try_status(builder)

    def eh_end_try(self, builder):
        """End a try-block
        """
        ctx = self._context
        cc = ctx.call_conv
        cc.unset_try_status(builder)
