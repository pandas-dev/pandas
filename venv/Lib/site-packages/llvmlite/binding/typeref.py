from ctypes import c_int, c_bool, c_void_p, c_uint64, c_uint, POINTER
import enum

from llvmlite import ir
from llvmlite.binding import ffi

# FIXME: Remove `opaque_pointers_enabled' when TP's are removed.
from llvmlite import opaque_pointers_enabled


class TypeKind(enum.IntEnum):
    # The LLVMTypeKind enum from llvm-c/Core.h

    void = 0
    half = 1
    float = 2
    double = 3
    x86_fp80 = 4
    fp128 = 5
    ppc_fp128 = 6
    label = 7
    integer = 8
    function = 9
    struct = 10
    array = 11
    pointer = 12
    vector = 13
    metadata = 14
    x86_mmx = 15
    token = 16
    scalable_vector = 17
    bfloat = 18
    x86_amx = 19


_TypeKindToIRType = {
    # All TypeKind here must have a TypeRef.as_ir() implementation
    TypeKind.void: ir.VoidType,
    TypeKind.half: ir.HalfType,
    TypeKind.float: ir.FloatType,
    TypeKind.double: ir.DoubleType,
    TypeKind.integer: ir.IntType,
    TypeKind.function: ir.FunctionType,
    TypeKind.pointer: ir.PointerType,
    TypeKind.array: ir.ArrayType,
    TypeKind.vector: ir.VectorType,
    TypeKind.struct: ir.LiteralStructType,
}


class TypeRef(ffi.ObjectRef):
    """A weak reference to a LLVM type
    """
    @property
    def name(self):
        """
        Get type name
        """
        return ffi.ret_string(ffi.lib.LLVMPY_GetTypeName(self))

    @property
    def is_struct(self):
        """
        Returns true if the type is a struct type.
        """
        return ffi.lib.LLVMPY_TypeIsStruct(self)

    @property
    def is_pointer(self):
        """
        Returns true if the type is a pointer type.
        """
        return ffi.lib.LLVMPY_TypeIsPointer(self)

    @property
    def is_array(self):
        """
        Returns true if the type is an array type.
        """
        return ffi.lib.LLVMPY_TypeIsArray(self)

    @property
    def is_vector(self):
        """
        Returns true if the type is a vector type.
        """
        return ffi.lib.LLVMPY_TypeIsVector(self)

    @property
    def is_function(self):
        """
        Returns true if the type is a function type.
        """
        return ffi.lib.LLVMPY_TypeIsFunction(self)

    @property
    def is_function_vararg(self):
        """
        Returns true if a function type accepts a variable number of arguments.
        When the type is not a function, raises exception.
        """
        if self.type_kind != TypeKind.function:
            raise ValueError("Type {} is not a function".format(self))
        return ffi.lib.LLVMPY_IsFunctionVararg(self)

    @property
    def elements(self):
        """
        Returns iterator over enclosing types
        """
        if self.is_pointer and opaque_pointers_enabled:
            raise ValueError("Type {} doesn't contain elements.".format(self))
        return _TypeListIterator(ffi.lib.LLVMPY_ElementIter(self))

    # FIXME: Remove me once typed pointers support is removed.
    @property
    def element_type(self):
        """
        Returns the pointed-to type. When the type is not a pointer,
        raises exception.
        """
        if not self.is_pointer:
            raise ValueError("Type {} is not a pointer".format(self))
        return TypeRef(ffi.lib.LLVMPY_GetElementType(self))

    @property
    def element_count(self):
        """
        Returns the number of elements in an array or a vector. For scalable
        vectors, returns minimum number of elements. When the type is neither
        an array nor a vector, raises exception.
        """
        if not self.is_array and not self.is_vector:
            raise ValueError("Type {} is not an array nor vector".format(self))
        return ffi.lib.LLVMPY_GetTypeElementCount(self)

    @property
    def type_width(self):
        """
        Return the basic size of this type if it is a primitive type. These are
        fixed by LLVM and are not target-dependent.
        This will return zero if the type does not have a size or is not a
        primitive type.

        If this is a scalable vector type, the scalable property will be set and
        the runtime size will be a positive integer multiple of the base size.

        Note that this may not reflect the size of memory allocated for an
        instance of the type or the number of bytes that are written when an
        instance of the type is stored to memory.
        """
        return ffi.lib.LLVMPY_GetTypeBitWidth(self)

    @property
    def type_kind(self):
        """
        Returns the LLVMTypeKind enumeration of this type.
        """
        return TypeKind(ffi.lib.LLVMPY_GetTypeKind(self))

    @property
    def is_packed_struct(self):
        return ffi.lib.LLVMPY_IsPackedStruct(self)

    @property
    def is_literal_struct(self):
        return ffi.lib.LLVMPY_IsLiteralStruct(self)

    @property
    def is_opaque_struct(self):
        return ffi.lib.LLVMPY_IsOpaqueStruct(self)

    def get_function_parameters(self) -> tuple["TypeRef"]:
        nparams = ffi.lib.LLVMPY_CountParamTypes(self)
        if nparams > 0:
            out_buffer = (ffi.LLVMTypeRef * nparams)(None)
            ffi.lib.LLVMPY_GetParamTypes(self, out_buffer)
            return tuple(map(TypeRef, out_buffer))
        else:
            return ()

    def get_function_return(self) -> "TypeRef":
        return TypeRef(ffi.lib.LLVMPY_GetReturnType(self))

    def as_ir(self, ir_ctx: ir.Context) -> ir.Type:
        """Convert into a ``llvmlite.ir.Type``.
        """
        try:
            cls = _TypeKindToIRType[self.type_kind]
        except KeyError:
            msg = f"as_ir() unsupported for TypeRef of {self.type_kind}"
            raise TypeError(msg)
        else:
            return cls.from_llvm(self, ir_ctx)

    def __str__(self):
        return ffi.ret_string(ffi.lib.LLVMPY_PrintType(self))


class _TypeIterator(ffi.ObjectRef):

    def __next__(self):
        vp = self._next()
        if vp:
            return TypeRef(vp)
        else:
            raise StopIteration

    next = __next__

    def __iter__(self):
        return self


class _TypeListIterator(_TypeIterator):

    def _dispose(self):
        self._capi.LLVMPY_DisposeElementIter(self)

    def _next(self):
        return ffi.lib.LLVMPY_ElementIterNext(self)


# FFI

ffi.lib.LLVMPY_PrintType.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_PrintType.restype = c_void_p

# FIXME: Remove me once typed pointers support is removed.
ffi.lib.LLVMPY_GetElementType.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_GetElementType.restype = ffi.LLVMTypeRef

ffi.lib.LLVMPY_TypeIsPointer.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_TypeIsPointer.restype = c_bool

ffi.lib.LLVMPY_TypeIsArray.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_TypeIsArray.restype = c_bool

ffi.lib.LLVMPY_TypeIsVector.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_TypeIsVector.restype = c_bool

ffi.lib.LLVMPY_TypeIsStruct.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_TypeIsStruct.restype = c_bool

ffi.lib.LLVMPY_TypeIsFunction.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_TypeIsFunction.restype = c_bool

ffi.lib.LLVMPY_IsPackedStruct.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_IsPackedStruct.restype = c_bool

ffi.lib.LLVMPY_IsOpaqueStruct.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_IsOpaqueStruct.restype = c_bool

ffi.lib.LLVMPY_IsLiteralStruct.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_IsLiteralStruct.restype = c_bool

ffi.lib.LLVMPY_GetReturnType.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_GetReturnType.restype = ffi.LLVMTypeRef

ffi.lib.LLVMPY_CountParamTypes.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_CountParamTypes.restype = c_uint

ffi.lib.LLVMPY_GetParamTypes.argtypes = [ffi.LLVMTypeRef,
                                         POINTER(ffi.LLVMTypeRef)]
ffi.lib.LLVMPY_GetParamTypes.restype = None

ffi.lib.LLVMPY_IsFunctionVararg.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_IsFunctionVararg.restype = c_bool

ffi.lib.LLVMPY_GetTypeKind.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_GetTypeKind.restype = c_int

ffi.lib.LLVMPY_GetTypeElementCount.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_GetTypeElementCount.restype = c_int

ffi.lib.LLVMPY_GetTypeBitWidth.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_GetTypeBitWidth.restype = c_uint64

ffi.lib.LLVMPY_ElementIter.argtypes = [ffi.LLVMTypeRef]
ffi.lib.LLVMPY_ElementIter.restype = ffi.LLVMElementIterator

ffi.lib.LLVMPY_ElementIterNext.argtypes = [ffi.LLVMElementIterator]
ffi.lib.LLVMPY_ElementIterNext.restype = ffi.LLVMTypeRef

ffi.lib.LLVMPY_DisposeElementIter.argtypes = [ffi.LLVMElementIterator]
