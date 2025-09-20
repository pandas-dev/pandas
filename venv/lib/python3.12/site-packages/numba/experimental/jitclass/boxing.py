"""
Implement logic relating to wrapping (box) and unwrapping (unbox) instances
of jitclasses for use inside the python interpreter.
"""

from functools import wraps, partial

from llvmlite import ir

from numba.core import types, cgutils
from numba.core.decorators import njit
from numba.core.pythonapi import box, unbox, NativeValue
from numba.core.typing.typeof import typeof_impl
from numba.experimental.jitclass import _box


_getter_code_template = """
def accessor(__numba_self_):
    return __numba_self_.{0}
"""

_setter_code_template = """
def mutator(__numba_self_, __numba_val):
    __numba_self_.{0} = __numba_val
"""

_method_code_template = """
def method(__numba_self_, *args):
    return __numba_self_.{method}(*args)
"""


def _generate_property(field, template, fname):
    """
    Generate simple function that get/set a field of the instance
    """
    source = template.format(field)
    glbls = {}
    exec(source, glbls)
    return njit(glbls[fname])


_generate_getter = partial(_generate_property, template=_getter_code_template,
                           fname='accessor')
_generate_setter = partial(_generate_property, template=_setter_code_template,
                           fname='mutator')


def _generate_method(name, func):
    """
    Generate a wrapper for calling a method.  Note the wrapper will only
    accept positional arguments.
    """
    source = _method_code_template.format(method=name)
    glbls = {}
    exec(source, glbls)
    method = njit(glbls['method'])

    @wraps(func)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)

    return wrapper


_cache_specialized_box = {}


def _specialize_box(typ):
    """
    Create a subclass of Box that is specialized to the jitclass.

    This function caches the result to avoid code bloat.
    """
    # Check cache
    if typ in _cache_specialized_box:
        return _cache_specialized_box[typ]
    dct = {'__slots__': (),
           '_numba_type_': typ,
           '__doc__': typ.class_type.class_doc,
           }
    # Inject attributes as class properties
    for field in typ.struct:
        getter = _generate_getter(field)
        setter = _generate_setter(field)
        dct[field] = property(getter, setter)
    # Inject properties as class properties
    for field, impdct in typ.jit_props.items():
        getter = None
        setter = None
        if 'get' in impdct:
            getter = _generate_getter(field)
        if 'set' in impdct:
            setter = _generate_setter(field)
        # get docstring from either the fget or fset
        imp = impdct.get('get') or impdct.get('set') or None
        doc = getattr(imp, '__doc__', None)
        dct[field] = property(getter, setter, doc=doc)
    # Inject methods as class members
    supported_dunders = {
        "__abs__",
        "__bool__",
        "__complex__",
        "__contains__",
        "__float__",
        "__getitem__",
        "__hash__",
        "__index__",
        "__invert__",
        "__int__",
        "__len__",
        "__setitem__",
        "__str__",
        "__eq__",
        "__ne__",
        "__ge__",
        "__gt__",
        "__le__",
        "__lt__",
        "__add__",
        "__floordiv__",
        "__lshift__",
        "__matmul__",
        "__mod__",
        "__mul__",
        "__neg__",
        "__pos__",
        "__pow__",
        "__rshift__",
        "__sub__",
        "__truediv__",
        "__and__",
        "__or__",
        "__xor__",
        "__iadd__",
        "__ifloordiv__",
        "__ilshift__",
        "__imatmul__",
        "__imod__",
        "__imul__",
        "__ipow__",
        "__irshift__",
        "__isub__",
        "__itruediv__",
        "__iand__",
        "__ior__",
        "__ixor__",
        "__radd__",
        "__rfloordiv__",
        "__rlshift__",
        "__rmatmul__",
        "__rmod__",
        "__rmul__",
        "__rpow__",
        "__rrshift__",
        "__rsub__",
        "__rtruediv__",
        "__rand__",
        "__ror__",
        "__rxor__",
    }
    for name, func in typ.methods.items():
        if name == "__init__":
            continue
        if (
            name.startswith("__")
            and name.endswith("__")
            and name not in supported_dunders
        ):
            raise TypeError(f"Method '{name}' is not supported.")
        dct[name] = _generate_method(name, func)

    # Inject static methods as class members
    for name, func in typ.static_methods.items():
        dct[name] = _generate_method(name, func)

    # Create subclass
    subcls = type(typ.classname, (_box.Box,), dct)
    # Store to cache
    _cache_specialized_box[typ] = subcls

    # Pre-compile attribute getter.
    # Note: This must be done after the "box" class is created because
    #       compiling the getter requires the "box" class to be defined.
    for k, v in dct.items():
        if isinstance(v, property):
            prop = getattr(subcls, k)
            if prop.fget is not None:
                fget = prop.fget
                fast_fget = fget.compile((typ,))
                fget.disable_compile()
                setattr(subcls, k,
                        property(fast_fget, prop.fset, prop.fdel,
                                 doc=prop.__doc__))

    return subcls


###############################################################################
# Implement box/unbox for call wrapper

@box(types.ClassInstanceType)
def _box_class_instance(typ, val, c):
    meminfo, dataptr = cgutils.unpack_tuple(c.builder, val)

    # Create Box instance
    box_subclassed = _specialize_box(typ)
    # Note: the ``box_subclassed`` is kept alive by the cache
    voidptr_boxcls = c.context.add_dynamic_addr(
        c.builder,
        id(box_subclassed),
        info="box_class_instance",
    )
    box_cls = c.builder.bitcast(voidptr_boxcls, c.pyapi.pyobj)

    box = c.pyapi.call_function_objargs(box_cls, ())

    # Initialize Box instance
    llvoidptr = ir.IntType(8).as_pointer()
    addr_meminfo = c.builder.bitcast(meminfo, llvoidptr)
    addr_data = c.builder.bitcast(dataptr, llvoidptr)

    def set_member(member_offset, value):
        # Access member by byte offset
        offset = c.context.get_constant(types.uintp, member_offset)
        ptr = cgutils.pointer_add(c.builder, box, offset)
        casted = c.builder.bitcast(ptr, llvoidptr.as_pointer())
        c.builder.store(value, casted)

    set_member(_box.box_meminfoptr_offset, addr_meminfo)
    set_member(_box.box_dataptr_offset, addr_data)
    return box


@unbox(types.ClassInstanceType)
def _unbox_class_instance(typ, val, c):
    def access_member(member_offset):
        # Access member by byte offset
        offset = c.context.get_constant(types.uintp, member_offset)
        llvoidptr = ir.IntType(8).as_pointer()
        ptr = cgutils.pointer_add(c.builder, val, offset)
        casted = c.builder.bitcast(ptr, llvoidptr.as_pointer())
        return c.builder.load(casted)

    struct_cls = cgutils.create_struct_proxy(typ)
    inst = struct_cls(c.context, c.builder)

    # load from Python object
    ptr_meminfo = access_member(_box.box_meminfoptr_offset)
    ptr_dataptr = access_member(_box.box_dataptr_offset)

    # store to native structure
    inst.meminfo = c.builder.bitcast(ptr_meminfo, inst.meminfo.type)
    inst.data = c.builder.bitcast(ptr_dataptr, inst.data.type)

    ret = inst._getvalue()

    c.context.nrt.incref(c.builder, typ, ret)

    return NativeValue(ret, is_error=c.pyapi.c_api_error())


# Add a typeof_impl implementation for boxed jitclasses to short-circut the
# various tests in typeof. This is needed for jitclasses which implement a
# custom hash method. Without this, typeof_impl will return None, and one of the
# later attempts to determine the type of the jitclass (before checking for
# _numba_type_) will look up the object in a dictionary, triggering the hash
# method. This will cause the dispatcher to determine the call signature of the
# jit decorated obj.__hash__ method, which will call typeof(obj), and thus
# infinite loop.
# This implementation is here instead of in typeof.py to avoid circular imports.
@typeof_impl.register(_box.Box)
def _typeof_jitclass_box(val, c):
    return getattr(type(val), "_numba_type_")
