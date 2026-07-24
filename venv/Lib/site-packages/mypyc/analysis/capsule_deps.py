from __future__ import annotations

from mypyc.ir.class_ir import ClassIR
from mypyc.ir.deps import Dependency
from mypyc.ir.func_ir import FuncIR
from mypyc.ir.ops import Assign, CallC, PrimitiveOp
from mypyc.ir.rtypes import RStruct, RTuple, RType, RUnion, RVec


def find_implicit_op_dependencies(fn: FuncIR) -> set[Dependency] | None:
    """Find implicit dependencies that need to be imported.

    Using primitives or types defined in librt submodules such as "librt.base64"
    requires dependency imports (e.g., capsule imports).

    Note that a module can depend on a librt module even if it doesn't explicitly
    import it, for example via re-exported names or via return types of functions
    defined in other modules.
    """
    deps: set[Dependency] | None = None
    # Check function signature types for dependencies
    deps = find_type_dependencies(fn, deps)
    # Check ops for dependencies
    for block in fn.blocks:
        for op in block.ops:
            assert not isinstance(op, PrimitiveOp), "Lowered IR is expected"
            if isinstance(op, CallC) and op.dependencies is not None:
                for dep in op.dependencies:
                    if deps is None:
                        deps = set()
                    deps.add(dep)
            deps = collect_type_deps(op.type, deps)
            if isinstance(op, Assign):
                deps = collect_type_deps(op.dest.type, deps)
    return deps


def find_type_dependencies(fn: FuncIR, deps: set[Dependency] | None) -> set[Dependency] | None:
    """Find dependencies from RTypes in function signatures.

    Some RTypes (e.g., those for librt types) have associated dependencies
    that need to be imported when the type is used.
    """
    # Check parameter types
    for arg in fn.decl.sig.args:
        deps = collect_type_deps(arg.type, deps)
    # Check return type
    deps = collect_type_deps(fn.decl.sig.ret_type, deps)
    return deps


def find_class_dependencies(cl: ClassIR) -> set[Dependency] | None:
    """Find dependencies from class attribute types."""
    deps: set[Dependency] | None = None
    for base in cl.mro:
        for attr_type in base.attributes.values():
            deps = collect_type_deps(attr_type, deps)
    return deps


def collect_type_deps(typ: RType, deps: set[Dependency] | None) -> set[Dependency] | None:
    """Collect dependencies from an RType, recursively checking compound types."""
    if typ.dependencies is not None:
        for dep in typ.dependencies:
            if deps is None:
                deps = set()
            deps.add(dep)
    if isinstance(typ, RUnion):
        for item in typ.items:
            deps = collect_type_deps(item, deps)
    elif isinstance(typ, RTuple):
        for item in typ.types:
            deps = collect_type_deps(item, deps)
    elif isinstance(typ, RStruct):
        for item in typ.types:
            deps = collect_type_deps(item, deps)
    elif isinstance(typ, RVec):
        deps = collect_type_deps(typ.item_type, deps)
    return deps
