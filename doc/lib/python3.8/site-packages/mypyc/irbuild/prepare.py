"""Prepare for IR transform.

This needs to run after type checking and before generating IR.

For example, construct partially initialized FuncIR and ClassIR
objects for all functions and classes. This allows us to bind
references to functions and classes before we've generated full IR for
functions or classes.  The actual IR transform will then populate all
the missing bits, such as function bodies (basic blocks).

Also build a mapping from mypy TypeInfos to ClassIR objects.
"""

from typing import List, Dict, Iterable, Optional, Union

from mypy.nodes import (
    MypyFile, TypeInfo, FuncDef, ClassDef, Decorator, OverloadedFuncDef, MemberExpr, Var,
    Expression, SymbolNode, ARG_STAR, ARG_STAR2
)
from mypy.types import Type
from mypy.build import Graph

from mypyc.ir.ops import DeserMaps
from mypyc.ir.rtypes import RInstance, tuple_rprimitive, dict_rprimitive
from mypyc.ir.func_ir import (
    FuncDecl, FuncSignature, RuntimeArg, FUNC_NORMAL, FUNC_STATICMETHOD, FUNC_CLASSMETHOD
)
from mypyc.ir.class_ir import ClassIR
from mypyc.common import PROPSET_PREFIX
from mypyc.irbuild.mapper import Mapper
from mypyc.irbuild.util import (
    get_func_def, is_dataclass, is_trait, is_extension_class, get_mypyc_attrs
)
from mypyc.errors import Errors
from mypyc.options import CompilerOptions
from mypyc.crash import catch_errors


def build_type_map(mapper: Mapper,
                   modules: List[MypyFile],
                   graph: Graph,
                   types: Dict[Expression, Type],
                   options: CompilerOptions,
                   errors: Errors) -> None:
    # Collect all classes defined in everything we are compiling
    classes = []
    for module in modules:
        module_classes = [node for node in module.defs if isinstance(node, ClassDef)]
        classes.extend([(module, cdef) for cdef in module_classes])

    # Collect all class mappings so that we can bind arbitrary class name
    # references even if there are import cycles.
    for module, cdef in classes:
        class_ir = ClassIR(cdef.name, module.fullname, is_trait(cdef),
                           is_abstract=cdef.info.is_abstract)
        class_ir.is_ext_class = is_extension_class(cdef)
        # If global optimizations are disabled, turn of tracking of class children
        if not options.global_opts:
            class_ir.children = None
        mapper.type_to_ir[cdef.info] = class_ir

    # Populate structural information in class IR for extension classes.
    for module, cdef in classes:
        with catch_errors(module.path, cdef.line):
            if mapper.type_to_ir[cdef.info].is_ext_class:
                prepare_class_def(module.path, module.fullname, cdef, errors, mapper)
            else:
                prepare_non_ext_class_def(module.path, module.fullname, cdef, errors, mapper)

    # Collect all the functions also. We collect from the symbol table
    # so that we can easily pick out the right copy of a function that
    # is conditionally defined.
    for module in modules:
        for func in get_module_func_defs(module):
            prepare_func_def(module.fullname, None, func, mapper)
            # TODO: what else?


def is_from_module(node: SymbolNode, module: MypyFile) -> bool:
    return node.fullname == module.fullname + '.' + node.name


def load_type_map(mapper: 'Mapper',
                  modules: List[MypyFile],
                  deser_ctx: DeserMaps) -> None:
    """Populate a Mapper with deserialized IR from a list of modules."""
    for module in modules:
        for name, node in module.names.items():
            if isinstance(node.node, TypeInfo) and is_from_module(node.node, module):
                ir = deser_ctx.classes[node.node.fullname]
                mapper.type_to_ir[node.node] = ir
                mapper.func_to_decl[node.node] = ir.ctor

    for module in modules:
        for func in get_module_func_defs(module):
            mapper.func_to_decl[func] = deser_ctx.functions[func.fullname].decl


def get_module_func_defs(module: MypyFile) -> Iterable[FuncDef]:
    """Collect all of the (non-method) functions declared in a module."""
    for name, node in module.names.items():
        # We need to filter out functions that are imported or
        # aliases.  The best way to do this seems to be by
        # checking that the fullname matches.
        if (isinstance(node.node, (FuncDef, Decorator, OverloadedFuncDef))
                and is_from_module(node.node, module)):
            yield get_func_def(node.node)


def prepare_func_def(module_name: str, class_name: Optional[str],
                     fdef: FuncDef, mapper: Mapper) -> FuncDecl:
    kind = FUNC_STATICMETHOD if fdef.is_static else (
        FUNC_CLASSMETHOD if fdef.is_class else FUNC_NORMAL)
    decl = FuncDecl(fdef.name, class_name, module_name, mapper.fdef_to_sig(fdef), kind)
    mapper.func_to_decl[fdef] = decl
    return decl


def prepare_method_def(ir: ClassIR, module_name: str, cdef: ClassDef, mapper: Mapper,
                       node: Union[FuncDef, Decorator]) -> None:
    if isinstance(node, FuncDef):
        ir.method_decls[node.name] = prepare_func_def(module_name, cdef.name, node, mapper)
    elif isinstance(node, Decorator):
        # TODO: do something about abstract methods here. Currently, they are handled just like
        # normal methods.
        decl = prepare_func_def(module_name, cdef.name, node.func, mapper)
        if not node.decorators:
            ir.method_decls[node.name] = decl
        elif isinstance(node.decorators[0], MemberExpr) and node.decorators[0].name == 'setter':
            # Make property setter name different than getter name so there are no
            # name clashes when generating C code, and property lookup at the IR level
            # works correctly.
            decl.name = PROPSET_PREFIX + decl.name
            decl.is_prop_setter = True
            ir.method_decls[PROPSET_PREFIX + node.name] = decl

        if node.func.is_property:
            assert node.func.type
            decl.is_prop_getter = True
            ir.property_types[node.name] = decl.sig.ret_type


def is_valid_multipart_property_def(prop: OverloadedFuncDef) -> bool:
    # Checks to ensure supported property decorator semantics
    if len(prop.items) == 2:
        getter = prop.items[0]
        setter = prop.items[1]
        if isinstance(getter, Decorator) and isinstance(setter, Decorator):
            if getter.func.is_property and len(setter.decorators) == 1:
                if isinstance(setter.decorators[0], MemberExpr):
                    if setter.decorators[0].name == "setter":
                        return True
    return False


def can_subclass_builtin(builtin_base: str) -> bool:
    # BaseException and dict are special cased.
    return builtin_base in (
        ('builtins.Exception', 'builtins.LookupError', 'builtins.IndexError',
        'builtins.Warning', 'builtins.UserWarning', 'builtins.ValueError',
        'builtins.object', ))


def prepare_class_def(path: str, module_name: str, cdef: ClassDef,
                      errors: Errors, mapper: Mapper) -> None:

    ir = mapper.type_to_ir[cdef.info]
    info = cdef.info

    attrs = get_mypyc_attrs(cdef)
    if attrs.get("allow_interpreted_subclasses") is True:
        ir.allow_interpreted_subclasses = True

    # We sort the table for determinism here on Python 3.5
    for name, node in sorted(info.names.items()):
        # Currently all plugin generated methods are dummies and not included.
        if node.plugin_generated:
            continue

        if isinstance(node.node, Var):
            assert node.node.type, "Class member %s missing type" % name
            if not node.node.is_classvar and name != '__slots__':
                ir.attributes[name] = mapper.type_to_rtype(node.node.type)
        elif isinstance(node.node, (FuncDef, Decorator)):
            prepare_method_def(ir, module_name, cdef, mapper, node.node)
        elif isinstance(node.node, OverloadedFuncDef):
            # Handle case for property with both a getter and a setter
            if node.node.is_property:
                if is_valid_multipart_property_def(node.node):
                    for item in node.node.items:
                        prepare_method_def(ir, module_name, cdef, mapper, item)
                else:
                    errors.error("Unsupported property decorator semantics", path, cdef.line)

            # Handle case for regular function overload
            else:
                assert node.node.impl
                prepare_method_def(ir, module_name, cdef, mapper, node.node.impl)

    # Check for subclassing from builtin types
    for cls in info.mro:
        # Special case exceptions and dicts
        # XXX: How do we handle *other* things??
        if cls.fullname == 'builtins.BaseException':
            ir.builtin_base = 'PyBaseExceptionObject'
        elif cls.fullname == 'builtins.dict':
            ir.builtin_base = 'PyDictObject'
        elif cls.fullname.startswith('builtins.'):
            if not can_subclass_builtin(cls.fullname):
                # Note that if we try to subclass a C extension class that
                # isn't in builtins, bad things will happen and we won't
                # catch it here! But this should catch a lot of the most
                # common pitfalls.
                errors.error("Inheriting from most builtin types is unimplemented",
                             path, cdef.line)

    if ir.builtin_base:
        ir.attributes.clear()

    # Set up a constructor decl
    init_node = cdef.info['__init__'].node
    if not ir.is_trait and not ir.builtin_base and isinstance(init_node, FuncDef):
        init_sig = mapper.fdef_to_sig(init_node)

        defining_ir = mapper.type_to_ir.get(init_node.info)
        # If there is a nontrivial __init__ that wasn't defined in an
        # extension class, we need to make the constructor take *args,
        # **kwargs so it can call tp_init.
        if ((defining_ir is None or not defining_ir.is_ext_class
             or cdef.info['__init__'].plugin_generated)
                and init_node.info.fullname != 'builtins.object'):
            init_sig = FuncSignature(
                [init_sig.args[0],
                 RuntimeArg("args", tuple_rprimitive, ARG_STAR),
                 RuntimeArg("kwargs", dict_rprimitive, ARG_STAR2)],
                init_sig.ret_type)

        ctor_sig = FuncSignature(init_sig.args[1:], RInstance(ir))
        ir.ctor = FuncDecl(cdef.name, None, module_name, ctor_sig)
        mapper.func_to_decl[cdef.info] = ir.ctor

    # Set up the parent class
    bases = [mapper.type_to_ir[base.type] for base in info.bases
             if base.type in mapper.type_to_ir]
    if not all(c.is_trait for c in bases[1:]):
        errors.error("Non-trait bases must appear first in parent list", path, cdef.line)
    ir.traits = [c for c in bases if c.is_trait]

    mro = []
    base_mro = []
    for cls in info.mro:
        if cls not in mapper.type_to_ir:
            if cls.fullname != 'builtins.object':
                ir.inherits_python = True
            continue
        base_ir = mapper.type_to_ir[cls]
        if not base_ir.is_trait:
            base_mro.append(base_ir)
        mro.append(base_ir)

        if cls.defn.removed_base_type_exprs or not base_ir.is_ext_class:
            ir.inherits_python = True

    base_idx = 1 if not ir.is_trait else 0
    if len(base_mro) > base_idx:
        ir.base = base_mro[base_idx]
    ir.mro = mro
    ir.base_mro = base_mro

    for base in bases:
        if base.children is not None:
            base.children.append(ir)

    if is_dataclass(cdef):
        ir.is_augmented = True


def prepare_non_ext_class_def(path: str, module_name: str, cdef: ClassDef,
                              errors: Errors, mapper: Mapper) -> None:

    ir = mapper.type_to_ir[cdef.info]
    info = cdef.info

    for name, node in info.names.items():
        if isinstance(node.node, (FuncDef, Decorator)):
            prepare_method_def(ir, module_name, cdef, mapper, node.node)
        elif isinstance(node.node, OverloadedFuncDef):
            # Handle case for property with both a getter and a setter
            if node.node.is_property:
                if not is_valid_multipart_property_def(node.node):
                    errors.error("Unsupported property decorator semantics", path, cdef.line)
                for item in node.node.items:
                    prepare_method_def(ir, module_name, cdef, mapper, item)
            # Handle case for regular function overload
            else:
                prepare_method_def(ir, module_name, cdef, mapper, get_func_def(node.node))

    if any(
        cls in mapper.type_to_ir and mapper.type_to_ir[cls].is_ext_class for cls in info.mro
    ):
        errors.error(
            "Non-extension classes may not inherit from extension classes", path, cdef.line)
