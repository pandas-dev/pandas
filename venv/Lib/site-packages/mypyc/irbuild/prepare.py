"""Prepare for IR transform.

This needs to run after type checking and before generating IR.

For example, construct partially initialized FuncIR and ClassIR
objects for all functions and classes. This allows us to bind
references to functions and classes before we've generated full IR for
functions or classes.  The actual IR transform will then populate all
the missing bits, such as function bodies (basic blocks).

Also build a mapping from mypy TypeInfos to ClassIR objects.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Final, NamedTuple

from mypy.build import Graph
from mypy.nodes import (
    ARG_STAR,
    ARG_STAR2,
    CallExpr,
    ClassDef,
    Decorator,
    Expression,
    FuncDef,
    IndexExpr,
    MemberExpr,
    MypyFile,
    NameExpr,
    OverloadedFuncDef,
    RefExpr,
    SymbolNode,
    TypeInfo,
    Var,
)
from mypy.semanal import refers_to_fullname
from mypy.traverser import TraverserVisitor
from mypy.types import Instance, Type, get_proper_type
from mypyc.common import FAST_PREFIX, PROPSET_PREFIX, SELF_NAME, get_id_from_name
from mypyc.crash import catch_errors
from mypyc.errors import Errors
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.func_ir import (
    FUNC_CLASSMETHOD,
    FUNC_NORMAL,
    FUNC_STATICMETHOD,
    FuncDecl,
    FuncSignature,
    RuntimeArg,
)
from mypyc.ir.ops import DeserMaps
from mypyc.ir.rtypes import (
    RInstance,
    RType,
    dict_rprimitive,
    none_rprimitive,
    object_pointer_rprimitive,
    object_rprimitive,
    tuple_rprimitive,
)
from mypyc.irbuild.mapper import Mapper
from mypyc.irbuild.util import (
    get_func_def,
    get_mypyc_attrs,
    is_dataclass,
    is_extension_class,
    is_trait,
)
from mypyc.options import CompilerOptions
from mypyc.sametype import is_same_type

GENERATOR_HELPER_NAME: Final = "__mypyc_generator_helper__"


def build_type_map(
    mapper: Mapper,
    modules: list[MypyFile],
    graph: Graph,
    types: dict[Expression, Type],
    options: CompilerOptions,
    errors: Errors,
) -> None:
    # Collect all classes defined in everything we are compiling
    classes = []
    for module in modules:
        module_classes = [node for node in module.defs if isinstance(node, ClassDef)]
        classes.extend([(module, cdef) for cdef in module_classes])

    # Collect all class mappings so that we can bind arbitrary class name
    # references even if there are import cycles.
    for module, cdef in classes:
        class_ir = ClassIR(
            cdef.name,
            module.fullname,
            is_trait(cdef),
            is_abstract=cdef.info.is_abstract,
            is_final_class=cdef.info.is_final,
        )
        class_ir.is_ext_class = is_extension_class(module.path, cdef, errors)
        if class_ir.is_ext_class:
            class_ir.deletable = cdef.info.deletable_attributes.copy()
        # If global optimizations are disabled, turn of tracking of class children
        if not options.global_opts:
            class_ir.children = None
        mapper.type_to_ir[cdef.info] = class_ir
        mapper.symbol_fullnames.add(class_ir.fullname)
        class_ir.is_enum = cdef.info.is_enum and len(cdef.info.enum_members) > 0

    # Populate structural information in class IR for extension classes.
    for module, cdef in classes:
        with catch_errors(module.path, cdef.line):
            if mapper.type_to_ir[cdef.info].is_ext_class:
                prepare_class_def(module.path, module.fullname, cdef, errors, mapper, options)
            else:
                prepare_non_ext_class_def(
                    module.path, module.fullname, cdef, errors, mapper, options
                )

    # Validate cross-class properties after all ClassIR flags are populated.
    for module, cdef in classes:
        with catch_errors(module.path, cdef.line):
            if mapper.type_to_ir[cdef.info].is_ext_class:
                validate_acyclic_class_bases(module.path, cdef, errors, mapper)

    # Prepare implicit attribute accessors as needed if an attribute overrides a property.
    for module, cdef in classes:
        class_ir = mapper.type_to_ir[cdef.info]
        if class_ir.is_ext_class:
            prepare_implicit_property_accessors(cdef.info, class_ir, module.fullname, mapper)

    # Collect all the functions also. We collect from the symbol table
    # so that we can easily pick out the right copy of a function that
    # is conditionally defined. This doesn't include nested functions!
    for module in modules:
        for func in get_module_func_defs(module):
            prepare_func_def(module.fullname, None, func, mapper, options)
            # TODO: what else?

    # Check for incompatible attribute definitions that were not
    # flagged by mypy but can't be supported when compiling.
    for module, cdef in classes:
        class_ir = mapper.type_to_ir[cdef.info]
        for attr in class_ir.attributes:
            for base_ir in class_ir.mro[1:]:
                if attr in base_ir.attributes:
                    if not is_same_type(class_ir.attributes[attr], base_ir.attributes[attr]):
                        node = cdef.info.names[attr].node
                        assert node is not None
                        kind = "trait" if base_ir.is_trait else "class"
                        errors.error(
                            f'Type of "{attr}" is incompatible with '
                            f'definition in {kind} "{base_ir.name}"',
                            module.path,
                            node.line,
                        )


def is_from_module(node: SymbolNode, module: MypyFile) -> bool:
    return node.fullname == module.fullname + "." + node.name


def load_type_map(mapper: Mapper, modules: list[MypyFile], deser_ctx: DeserMaps) -> None:
    """Populate a Mapper with deserialized IR from a list of modules."""
    for module in modules:
        for node in module.names.values():
            if (
                isinstance(node.node, TypeInfo)
                and is_from_module(node.node, module)
                and not node.node.is_newtype
                and not node.node.is_named_tuple
                and node.node.typeddict_type is None
            ):
                # Some TypeInfo entries are mypy-synthetic (e.g. anonymous
                # intersection classes like "<subclass of X and Y>") and have
                # no corresponding mypyc ClassIR. Skip those rather than
                # aborting the whole cache load.
                ir = deser_ctx.classes.get(node.node.fullname)
                if ir is None:
                    continue
                mapper.type_to_ir[node.node] = ir
                mapper.symbol_fullnames.add(node.node.fullname)
                # Trait/builtin-base classes have an ir.ctor FuncDecl
                # but no emitted CPyDef_<ctor>, so a cross-group direct
                # call would hit an undefined symbol. Mirror the skip
                # in prepare_init_method.
                if not ir.is_trait and not ir.builtin_base:
                    mapper.func_to_decl[node.node] = ir.ctor

    for module in modules:
        for func in get_module_func_defs(module):
            func_id = get_id_from_name(func.name, func.fullname, func.line)
            mapper.func_to_decl[func] = deser_ctx.functions[func_id].decl


def get_module_func_defs(module: MypyFile) -> Iterable[FuncDef]:
    """Collect all of the (non-method) functions declared in a module."""
    for node in module.names.values():
        # We need to filter out functions that are imported or
        # aliases.  The best way to do this seems to be by
        # checking that the fullname matches.
        if isinstance(node.node, (FuncDef, Decorator, OverloadedFuncDef)) and is_from_module(
            node.node, module
        ):
            yield get_func_def(node.node)


def prepare_func_def(
    module_name: str,
    class_name: str | None,
    fdef: FuncDef,
    mapper: Mapper,
    options: CompilerOptions,
) -> FuncDecl:
    kind = (
        FUNC_CLASSMETHOD
        if fdef.is_class
        else (FUNC_STATICMETHOD if fdef.is_static else FUNC_NORMAL)
    )
    sig = mapper.fdef_to_sig(fdef, options.strict_dunders_typing)
    decl = FuncDecl(
        fdef.name,
        class_name,
        module_name,
        sig,
        kind,
        is_generator=fdef.is_generator,
        is_coroutine=fdef.is_coroutine,
    )
    mapper.func_to_decl[fdef] = decl
    return decl


def create_generator_class_for_func(
    module_name: str, class_name: str | None, fdef: FuncDef, mapper: Mapper, name_suffix: str = ""
) -> ClassIR:
    """For a generator/async function, declare a generator class.

    Each generator and async function gets a dedicated class that implements the
    generator protocol with generated methods.
    """
    assert fdef.is_coroutine or fdef.is_generator
    name = "_".join(x for x in [fdef.name, class_name] if x) + "_gen" + name_suffix
    cir = ClassIR(name, module_name, is_generated=True, is_final_class=class_name is None)
    cir.reuse_freed_instance = True
    mapper.fdef_to_generator[fdef] = cir

    helper_sig = FuncSignature(
        (
            RuntimeArg(SELF_NAME, object_rprimitive),
            RuntimeArg("type", object_rprimitive),
            RuntimeArg("value", object_rprimitive),
            RuntimeArg("traceback", object_rprimitive),
            RuntimeArg("arg", object_rprimitive),
            # If non-NULL, used to store return value instead of raising StopIteration(retv)
            RuntimeArg("stop_iter_ptr", object_pointer_rprimitive),
        ),
        object_rprimitive,
    )

    # The implementation of most generator functionality is behind this magic method.
    helper_fn_decl = FuncDecl(GENERATOR_HELPER_NAME, name, module_name, helper_sig, internal=True)
    cir.method_decls[helper_fn_decl.name] = helper_fn_decl
    return cir


def prepare_method_def(
    ir: ClassIR,
    module_name: str,
    cdef: ClassDef,
    mapper: Mapper,
    node: FuncDef | Decorator,
    options: CompilerOptions,
) -> None:
    if isinstance(node, FuncDef):
        ir.method_decls[node.name] = prepare_func_def(
            module_name, cdef.name, node, mapper, options
        )
    elif isinstance(node, Decorator):
        # TODO: do something about abstract methods here. Currently, they are handled just like
        # normal methods.
        decl = prepare_func_def(module_name, cdef.name, node.func, mapper, options)
        if not node.decorators:
            ir.method_decls[node.name] = decl
        elif isinstance(node.decorators[0], MemberExpr) and node.decorators[0].name == "setter":
            # Make property setter name different than getter name so there are no
            # name clashes when generating C code, and property lookup at the IR level
            # works correctly.
            decl.name = PROPSET_PREFIX + decl.name
            decl.is_prop_setter = True
            # Making the argument implicitly positional-only avoids unnecessary glue methods
            decl.sig.args[1].pos_only = True
            ir.method_decls[PROPSET_PREFIX + node.name] = decl

        if node.func.is_property:
            assert node.func.type, f"Expected return type annotation for property '{node.name}'"
            decl.is_prop_getter = True
            ir.property_types[node.name] = decl.sig.ret_type


def prepare_fast_path(
    ir: ClassIR,
    module_name: str,
    cdef: ClassDef,
    mapper: Mapper,
    node: SymbolNode | None,
    options: CompilerOptions,
) -> None:
    """Add fast (direct) variants of methods in non-extension classes."""
    if ir.is_enum:
        # We check that non-empty enums are implicitly final in mypy, so we
        # can generate direct calls to enum methods.
        if isinstance(node, OverloadedFuncDef):
            if node.is_property:
                return
            node = node.impl
        if not isinstance(node, FuncDef):
            # TODO: support decorated methods (at least @classmethod and @staticmethod).
            return
        # The simplest case is a regular or overloaded method without decorators. In this
        # case we can generate practically identical IR method body, but with a signature
        # suitable for direct calls (usual non-extension class methods are converted to
        # callable classes, and thus have an extra __mypyc_self__ argument).
        name = FAST_PREFIX + node.name
        sig = mapper.fdef_to_sig(node, options.strict_dunders_typing)
        decl = FuncDecl(name, cdef.name, module_name, sig, FUNC_NORMAL)
        ir.method_decls[name] = decl
    return


def is_valid_multipart_property_def(prop: OverloadedFuncDef) -> bool:
    # Checks to ensure supported property decorator semantics
    if len(prop.items) != 2:
        return False

    getter = prop.items[0]
    setter = prop.items[1]

    return (
        isinstance(getter, Decorator)
        and isinstance(setter, Decorator)
        and getter.func.is_property
        and len(setter.decorators) == 1
        and isinstance(setter.decorators[0], MemberExpr)
        and setter.decorators[0].name == "setter"
    )


def can_subclass_builtin(builtin_base: str) -> bool:
    # BaseException and dict are special cased.
    return builtin_base in (
        (
            "builtins.Exception",
            "builtins.LookupError",
            "builtins.IndexError",
            "builtins.Warning",
            "builtins.UserWarning",
            "builtins.ValueError",
            "builtins.object",
        )
    )


def get_removed_base_fullname(expr: Expression) -> str | None:
    if isinstance(expr, IndexExpr):
        expr = expr.base
    if isinstance(expr, RefExpr):
        return expr.fullname
    return None


def find_non_acyclic_base(cdef: ClassDef, mapper: Mapper) -> str | None:
    if cdef.type_args:
        return "typing.Generic"

    for expr in cdef.removed_base_type_exprs:
        if fullname := get_removed_base_fullname(expr):
            return fullname
        return "a removed base class"

    for base in cdef.info.mro[1:]:
        if base.fullname == "builtins.object":
            continue

        base_ir = mapper.type_to_ir.get(base)
        if base_ir is not None and base_ir.is_acyclic:
            continue

        return base.fullname

    return None


def validate_acyclic_class_bases(
    path: str, cdef: ClassDef, errors: Errors, mapper: Mapper
) -> None:
    ir = mapper.type_to_ir[cdef.info]
    if not ir.is_acyclic:
        return

    if fullname := find_non_acyclic_base(cdef, mapper):
        errors.error(
            f'"acyclic" can\'t be used in a class that inherits from non-acyclic type "{fullname}"',
            path,
            cdef.line,
        )


def prepare_class_def(
    path: str,
    module_name: str,
    cdef: ClassDef,
    errors: Errors,
    mapper: Mapper,
    options: CompilerOptions,
) -> None:
    """Populate the interface-level information in a class IR.

    This includes attribute and method declarations, and the MRO, among other things, but
    method bodies are generated in a later pass.
    """

    ir = mapper.type_to_ir[cdef.info]
    info = cdef.info

    attrs, attrs_lines = get_mypyc_attrs(cdef, path, errors)
    if attrs.get("allow_interpreted_subclasses") is True:
        ir.allow_interpreted_subclasses = True
    if attrs.get("serializable") is True:
        # Supports copy.copy and pickle (including subclasses)
        ir._serializable = True

    if attrs.get("acyclic") is True:
        ir.is_acyclic = True

    # Check for subclassing from builtin types
    for cls in info.mro:
        # Special case exceptions and dicts
        # XXX: How do we handle *other* things??
        if cls.fullname == "builtins.BaseException":
            ir.builtin_base = "PyBaseExceptionObject"
        elif cls.fullname == "builtins.dict":
            ir.builtin_base = "PyDictObject"
        elif cls.fullname.startswith("builtins."):
            if not can_subclass_builtin(cls.fullname):
                # Note that if we try to subclass a C extension class that
                # isn't in builtins, bad things will happen and we won't
                # catch it here! But this should catch a lot of the most
                # common pitfalls.
                errors.error(
                    "Inheriting from most builtin types is unimplemented", path, cdef.line
                )
                errors.note(
                    "Potential workaround: @mypy_extensions.mypyc_attr(native_class=False)",
                    path,
                    cdef.line,
                )
                errors.note(
                    "https://mypyc.readthedocs.io/en/stable/native_classes.html#defining-non-native-classes",
                    path,
                    cdef.line,
                )

    free_list_len = attrs.get("free_list_len")
    if free_list_len is not None:
        line = attrs_lines["free_list_len"]
        if ir.is_trait:
            errors.error('"free_list_len" can\'t be used with traits', path, line)
        if ir.allow_interpreted_subclasses:
            errors.error(
                '"free_list_len" can\'t be used in a class that allows interpreted subclasses',
                path,
                line,
            )
        if ir.builtin_base:
            errors.error(
                '"free_list_len" can\'t be used in a class that inherits from a built-in type',
                path,
                line,
            )
        if free_list_len == 1:
            ir.reuse_freed_instance = True
        else:
            errors.error(f'Unsupported value for "free_list_len": {free_list_len}', path, line)

    # Set up the parent class
    bases = [mapper.type_to_ir[base.type] for base in info.bases if base.type in mapper.type_to_ir]
    if len(bases) > 1 and any(not c.is_trait for c in bases) and bases[0].is_trait:
        # If the first base is a non-trait, don't ever error here. While it is correct
        # to error if a trait comes before the next non-trait base (e.g. non-trait, trait,
        # non-trait), it's pointless, confusing noise from the bigger issue: multiple
        # inheritance is *not* supported.
        errors.error("Non-trait base must appear first in parent list", path, cdef.line)
    ir.traits = [c for c in bases if c.is_trait]

    mro = []  # All mypyc base classes
    base_mro = []  # Non-trait mypyc base classes
    for cls in info.mro:
        if cls not in mapper.type_to_ir:
            if cls.fullname != "builtins.object":
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

    prepare_methods_and_attributes(cdef, ir, path, module_name, errors, mapper, options)
    prepare_init_method(cdef, ir, module_name, mapper)

    for base in bases:
        if base.children is not None:
            base.children.append(ir)

    if is_dataclass(cdef):
        ir.is_augmented = True


def prepare_methods_and_attributes(
    cdef: ClassDef,
    ir: ClassIR,
    path: str,
    module_name: str,
    errors: Errors,
    mapper: Mapper,
    options: CompilerOptions,
) -> None:
    """Populate attribute and method declarations."""
    info = cdef.info
    for name, node in info.names.items():
        # Currently all plugin generated methods are dummies and not included.
        if node.plugin_generated:
            continue

        if isinstance(node.node, Var):
            assert node.node.type, "Class member %s missing type" % name
            if not node.node.is_classvar and name not in ("__slots__", "__deletable__"):
                attr_rtype = mapper.type_to_rtype(node.node.type)
                if ir.is_trait and attr_rtype.error_overlap:
                    # Traits don't have attribute definedness bitmaps, so use
                    # property accessor methods to access attributes that need them.
                    # We will generate accessor implementations that use the class bitmap
                    # for any concrete subclasses.
                    add_getter_declaration(ir, name, attr_rtype, module_name)
                    add_setter_declaration(ir, name, attr_rtype, module_name)
                ir.attributes[name] = attr_rtype
        elif isinstance(node.node, (FuncDef, Decorator)):
            prepare_method_def(ir, module_name, cdef, mapper, node.node, options)
        elif isinstance(node.node, OverloadedFuncDef):
            # Handle case for property with both a getter and a setter
            if node.node.is_property:
                if is_valid_multipart_property_def(node.node):
                    for item in node.node.items:
                        prepare_method_def(ir, module_name, cdef, mapper, item, options)
                else:
                    errors.error("Unsupported property decorator semantics", path, cdef.line)

            # Handle case for regular function overload
            else:
                if not node.node.impl:
                    errors.error(
                        "Overloads without implementation are not supported", path, cdef.line
                    )
                else:
                    prepare_method_def(ir, module_name, cdef, mapper, node.node.impl, options)

    if ir.builtin_base:
        ir.attributes.clear()


def prepare_implicit_property_accessors(
    info: TypeInfo, ir: ClassIR, module_name: str, mapper: Mapper
) -> None:
    concrete_attributes = set()
    for base in ir.base_mro:
        for name, attr_rtype in base.attributes.items():
            concrete_attributes.add(name)
            add_property_methods_for_attribute_if_needed(
                info, ir, name, attr_rtype, module_name, mapper
            )
    for base in ir.mro[1:]:
        if base.is_trait:
            for name, attr_rtype in base.attributes.items():
                if name not in concrete_attributes:
                    add_property_methods_for_attribute_if_needed(
                        info, ir, name, attr_rtype, module_name, mapper
                    )


def add_property_methods_for_attribute_if_needed(
    info: TypeInfo,
    ir: ClassIR,
    attr_name: str,
    attr_rtype: RType,
    module_name: str,
    mapper: Mapper,
) -> None:
    """Add getter and/or setter for attribute if defined as property in a base class.

    Only add declarations. The body IR will be synthesized later during irbuild.
    """
    for base in info.mro[1:]:
        if base in mapper.type_to_ir:
            base_ir = mapper.type_to_ir[base]
            n = base.names.get(attr_name)
            if n is None:
                continue
            node = n.node
            if isinstance(node, Decorator) and node.name not in ir.method_decls:
                # Defined as a read-only property in base class/trait
                add_getter_declaration(ir, attr_name, attr_rtype, module_name)
            elif isinstance(node, OverloadedFuncDef) and is_valid_multipart_property_def(node):
                # Defined as a read-write property in base class/trait
                add_getter_declaration(ir, attr_name, attr_rtype, module_name)
                add_setter_declaration(ir, attr_name, attr_rtype, module_name)
            elif base_ir.is_trait and attr_rtype.error_overlap:
                add_getter_declaration(ir, attr_name, attr_rtype, module_name)
                add_setter_declaration(ir, attr_name, attr_rtype, module_name)


def add_getter_declaration(
    ir: ClassIR, attr_name: str, attr_rtype: RType, module_name: str
) -> None:
    self_arg = RuntimeArg("self", RInstance(ir), pos_only=True)
    sig = FuncSignature([self_arg], attr_rtype)
    decl = FuncDecl(attr_name, ir.name, module_name, sig, FUNC_NORMAL)
    decl.is_prop_getter = True
    decl.implicit = True  # Triggers synthesization
    ir.method_decls[attr_name] = decl
    ir.property_types[attr_name] = attr_rtype  # TODO: Needed??


def add_setter_declaration(
    ir: ClassIR, attr_name: str, attr_rtype: RType, module_name: str
) -> None:
    self_arg = RuntimeArg("self", RInstance(ir), pos_only=True)
    value_arg = RuntimeArg("value", attr_rtype, pos_only=True)
    sig = FuncSignature([self_arg, value_arg], none_rprimitive)
    setter_name = PROPSET_PREFIX + attr_name
    decl = FuncDecl(setter_name, ir.name, module_name, sig, FUNC_NORMAL)
    decl.is_prop_setter = True
    decl.implicit = True  # Triggers synthesization
    ir.method_decls[setter_name] = decl


def check_matching_args(init_sig: FuncSignature, new_sig: FuncSignature) -> bool:
    num_init_args = len(init_sig.args) - init_sig.num_bitmap_args
    num_new_args = len(new_sig.args) - new_sig.num_bitmap_args
    if num_init_args != num_new_args:
        return False

    for idx in range(1, num_init_args):
        init_arg = init_sig.args[idx]
        new_arg = new_sig.args[idx]
        if init_arg.type != new_arg.type:
            return False

        if init_arg.kind != new_arg.kind:
            return False

    return True


def prepare_init_method(cdef: ClassDef, ir: ClassIR, module_name: str, mapper: Mapper) -> None:
    # Set up a constructor decl
    init_node = cdef.info["__init__"].node

    new_node: SymbolNode | None = None
    new_symbol = cdef.info.get("__new__")
    # We are only interested in __new__ method defined in a user-defined class,
    # so we ignore it if it comes from a builtin type. It's usually builtins.object
    # but could also be builtins.type for metaclasses so we detect the prefix which
    # matches both.
    if new_symbol and new_symbol.fullname and not new_symbol.fullname.startswith("builtins."):
        new_node = new_symbol.node
    if isinstance(new_node, (Decorator, OverloadedFuncDef)):
        new_node = get_func_def(new_node)
    if not ir.is_trait and not ir.builtin_base and isinstance(init_node, FuncDef):
        init_sig = mapper.fdef_to_sig(init_node, True)
        args_match = True
        if isinstance(new_node, FuncDef):
            new_sig = mapper.fdef_to_sig(new_node, True)
            args_match = check_matching_args(init_sig, new_sig)

        defining_ir = mapper.type_to_ir.get(init_node.info)
        # If there is a nontrivial __init__ that wasn't defined in an
        # extension class, we need to make the constructor take *args,
        # **kwargs so it can call tp_init.
        if (
            (
                defining_ir is None
                or not defining_ir.is_ext_class
                or cdef.info["__init__"].plugin_generated
            )
            and init_node.info.fullname != "builtins.object"
        ) or not args_match:
            init_sig = FuncSignature(
                [
                    init_sig.args[0],
                    RuntimeArg("args", tuple_rprimitive, ARG_STAR),
                    RuntimeArg("kwargs", dict_rprimitive, ARG_STAR2),
                ],
                init_sig.ret_type,
            )

        last_arg = len(init_sig.args) - init_sig.num_bitmap_args
        ctor_sig = FuncSignature(init_sig.args[1:last_arg], RInstance(ir))
        ir.ctor = FuncDecl(cdef.name, None, module_name, ctor_sig)
        mapper.func_to_decl[cdef.info] = ir.ctor


def prepare_non_ext_class_def(
    path: str,
    module_name: str,
    cdef: ClassDef,
    errors: Errors,
    mapper: Mapper,
    options: CompilerOptions,
) -> None:
    ir = mapper.type_to_ir[cdef.info]
    info = cdef.info

    for node in info.names.values():
        if isinstance(node.node, (FuncDef, Decorator)):
            prepare_method_def(ir, module_name, cdef, mapper, node.node, options)
        elif isinstance(node.node, OverloadedFuncDef):
            # Handle case for property with both a getter and a setter
            if node.node.is_property:
                if not is_valid_multipart_property_def(node.node):
                    errors.error("Unsupported property decorator semantics", path, cdef.line)
                for item in node.node.items:
                    prepare_method_def(ir, module_name, cdef, mapper, item, options)
            # Handle case for regular function overload
            else:
                prepare_method_def(ir, module_name, cdef, mapper, get_func_def(node.node), options)

        prepare_fast_path(ir, module_name, cdef, mapper, node.node, options)

    if any(cls in mapper.type_to_ir and mapper.type_to_ir[cls].is_ext_class for cls in info.mro):
        errors.error(
            "Non-extension classes may not inherit from extension classes", path, cdef.line
        )


RegisterImplInfo = tuple[TypeInfo, FuncDef]


class SingledispatchInfo(NamedTuple):
    singledispatch_impls: dict[FuncDef, list[RegisterImplInfo]]
    decorators_to_remove: dict[FuncDef, list[int]]


def find_singledispatch_register_impls(
    modules: list[MypyFile], errors: Errors
) -> SingledispatchInfo:
    visitor = SingledispatchVisitor(errors)
    for module in modules:
        visitor.current_path = module.path
        module.accept(visitor)
    return SingledispatchInfo(visitor.singledispatch_impls, visitor.decorators_to_remove)


class SingledispatchVisitor(TraverserVisitor):
    current_path: str

    def __init__(self, errors: Errors) -> None:
        super().__init__()

        # Map of main singledispatch function to list of registered implementations
        self.singledispatch_impls: defaultdict[FuncDef, list[RegisterImplInfo]] = defaultdict(list)

        # Map of decorated function to the indices of any decorators to remove
        self.decorators_to_remove: dict[FuncDef, list[int]] = {}

        self.errors: Errors = errors
        self.func_stack_depth = 0

    def visit_func_def(self, o: FuncDef) -> None:
        self.func_stack_depth += 1
        super().visit_func_def(o)
        self.func_stack_depth -= 1

    def visit_decorator(self, dec: Decorator) -> None:
        if dec.decorators:
            decorators_to_store = dec.decorators.copy()
            decorators_to_remove: list[int] = []
            # the index of the last non-register decorator before finding a register decorator
            # when going through decorators from top to bottom
            last_non_register: int | None = None
            for i, d in enumerate(decorators_to_store):
                impl = get_singledispatch_register_call_info(d, dec.func)
                if impl is not None:
                    if self.func_stack_depth > 0:
                        self.errors.error(
                            "Registering nested functions not supported", self.current_path, d.line
                        )
                    self.singledispatch_impls[impl.singledispatch_func].append(
                        (impl.dispatch_type, dec.func)
                    )
                    decorators_to_remove.append(i)
                    if last_non_register is not None:
                        # found a register decorator after a non-register decorator, which we
                        # don't support because we'd have to make a copy of the function before
                        # calling the decorator so that we can call it later, which complicates
                        # the implementation for something that is probably not commonly used
                        self.errors.error(
                            "Calling decorator after registering function not supported",
                            self.current_path,
                            decorators_to_store[last_non_register].line,
                        )
                else:
                    if refers_to_fullname(d, "functools.singledispatch"):
                        if self.func_stack_depth > 0:
                            self.errors.error(
                                "Nested singledispatch functions not supported",
                                self.current_path,
                                d.line,
                            )
                        decorators_to_remove.append(i)
                        # make sure that we still treat the function as a singledispatch function
                        # even if we don't find any registered implementations (which might happen
                        # if all registered implementations are registered dynamically)
                        self.singledispatch_impls.setdefault(dec.func, [])
                    last_non_register = i

            if decorators_to_remove:
                # calling register on a function that tries to dispatch based on type annotations
                # raises a TypeError because compiled functions don't have an __annotations__
                # attribute
                self.decorators_to_remove[dec.func] = decorators_to_remove

        super().visit_decorator(dec)


class RegisteredImpl(NamedTuple):
    singledispatch_func: FuncDef
    dispatch_type: TypeInfo


def get_singledispatch_register_call_info(
    decorator: Expression, func: FuncDef
) -> RegisteredImpl | None:
    # @fun.register(complex)
    # def g(arg): ...
    if (
        isinstance(decorator, CallExpr)
        and len(decorator.args) == 1
        and isinstance(decorator.args[0], RefExpr)
    ):
        callee = decorator.callee
        dispatch_type = decorator.args[0].node
        if not isinstance(dispatch_type, TypeInfo):
            return None

        if isinstance(callee, MemberExpr):
            return registered_impl_from_possible_register_call(callee, dispatch_type)
    # @fun.register
    # def g(arg: int): ...
    elif isinstance(decorator, MemberExpr):
        # we don't know if this is a register call yet, so we can't be sure that the function
        # actually has arguments
        if not func.arguments:
            return None
        arg_type = get_proper_type(func.arguments[0].variable.type)
        if not isinstance(arg_type, Instance):
            return None
        info = arg_type.type
        return registered_impl_from_possible_register_call(decorator, info)
    return None


def registered_impl_from_possible_register_call(
    expr: MemberExpr, dispatch_type: TypeInfo
) -> RegisteredImpl | None:
    if expr.name == "register" and isinstance(expr.expr, NameExpr):
        node = expr.expr.node
        if isinstance(node, Decorator):
            return RegisteredImpl(node.func, dispatch_type)
    return None


def adjust_generator_classes_of_methods(mapper: Mapper) -> None:
    """Make optimizations and adjustments to generated generator classes of methods.

    This is a separate pass after type map has been built, since we need all classes
    to be processed to analyze class hierarchies.
    """

    generator_methods = []

    for fdef, fn_ir in mapper.func_to_decl.items():
        if isinstance(fdef, FuncDef) and (fdef.is_coroutine or fdef.is_generator):
            gen_ir = create_generator_class_for_func(
                fn_ir.module_name, fn_ir.class_name, fdef, mapper
            )
            # TODO: We could probably support decorators sometimes (static and class method?)
            if not fdef.is_decorated:
                name = fn_ir.name
                precise_ret_type = True
                if fn_ir.class_name is not None:
                    class_ir = mapper.type_to_ir[fdef.info]
                    subcls = class_ir.subclasses()
                    if subcls is None:
                        # Override could be of a different type, so we can't make assumptions.
                        precise_ret_type = False
                    elif class_ir.is_trait:
                        # Give up on traits. We could possibly have an abstract base class
                        # for generator return types to make this use precise types.
                        precise_ret_type = False
                    else:
                        for s in subcls:
                            if name in s.method_decls:
                                m = s.method_decls[name]
                                if (
                                    m.is_generator != fn_ir.is_generator
                                    or m.is_coroutine != fn_ir.is_coroutine
                                ):
                                    # Override is of a different kind, and the optimization
                                    # to use a precise generator return type doesn't work.
                                    precise_ret_type = False
                else:
                    class_ir = None

                if precise_ret_type:
                    # Give a more precise type for generators, so that we can optimize
                    # code that uses them. They return a generator object, which has a
                    # specific class. Without this, the type would have to be 'object'.
                    fn_ir.sig.ret_type = RInstance(gen_ir)
                    if fn_ir.bound_sig:
                        fn_ir.bound_sig.ret_type = RInstance(gen_ir)
                    if class_ir is not None:
                        if class_ir.is_method_final(name):
                            gen_ir.is_final_class = True
                        generator_methods.append((name, class_ir, gen_ir))

    new_bases = {}

    for name, class_ir, gen in generator_methods:
        # For generator methods, we need to have subclass generator classes inherit from
        # baseclass generator classes when there are overrides to maintain LSP.
        base = class_ir.real_base()
        if base is not None:
            if base.has_method(name):
                base_sig = base.method_sig(name)
                if isinstance(base_sig.ret_type, RInstance):
                    base_gen = base_sig.ret_type.class_ir
                    new_bases[gen] = base_gen

    # Add generator inheritance relationships by adjusting MROs.
    for deriv, base in new_bases.items():
        if base.children is not None:
            base.children.append(deriv)
        while True:
            deriv.mro.append(base)
            deriv.base_mro.append(base)
            if base not in new_bases:
                break
            base = new_bases[base]
