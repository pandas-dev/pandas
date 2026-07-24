"""Intermediate representation of modules."""

from __future__ import annotations

from mypyc.common import JsonDict
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.deps import Capsule, Dependency, HeaderDep, SourceDep
from mypyc.ir.func_ir import FuncDecl, FuncIR
from mypyc.ir.ops import DeserMaps
from mypyc.ir.rtypes import RType, deserialize_type


class ModuleIR:
    """Intermediate representation of a module."""

    def __init__(
        self,
        fullname: str,
        imports: list[str],
        functions: list[FuncIR],
        classes: list[ClassIR],
        final_names: list[tuple[str, RType]],
        type_var_names: list[str],
    ) -> None:
        self.fullname = fullname
        self.imports = imports.copy()
        self.functions = functions
        self.classes = classes
        self.final_names = final_names
        # Names of C statics used for Python 3.12 type variable objects.
        # These are only visible in the module that defined them, so no need
        # to serialize.
        self.type_var_names = type_var_names
        # Dependencies needed by the module (such as capsules or source files)
        self.dependencies: set[Dependency] = set()

    def serialize(self) -> JsonDict:
        # Serialize dependencies as a list of dicts with type information
        serialized_deps = []
        for dep in sorted(self.dependencies, key=lambda d: (type(d).__name__, str(d))):
            if isinstance(dep, Capsule):
                serialized_deps.append({"type": "Capsule", "name": dep.name})
            elif isinstance(dep, SourceDep):
                source_dep: JsonDict = {
                    "type": "SourceDep",
                    "path": dep.path,
                    "include_dirs": dep.include_dirs,
                    "internal": dep.internal,
                }
                serialized_deps.append(source_dep)
            elif isinstance(dep, HeaderDep):
                header_dep: JsonDict = {
                    "type": "HeaderDep",
                    "path": dep.path,
                    "include_dirs": dep.include_dirs,
                    "internal": dep.internal,
                }
                serialized_deps.append(header_dep)

        return {
            "fullname": self.fullname,
            "imports": self.imports,
            "functions": [f.serialize() for f in self.functions],
            "classes": [c.serialize() for c in self.classes],
            "final_names": [(k, t.serialize()) for k, t in self.final_names],
            "dependencies": serialized_deps,
        }

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> ModuleIR:
        module = ModuleIR(
            data["fullname"],
            data["imports"],
            [ctx.functions[FuncDecl.get_id_from_json(f)] for f in data["functions"]],
            [ClassIR.deserialize(c, ctx) for c in data["classes"]],
            [(k, deserialize_type(t, ctx)) for k, t in data["final_names"]],
            [],
        )

        # Deserialize dependencies
        deps: set[Dependency] = set()
        for dep_dict in data["dependencies"]:
            if dep_dict["type"] == "Capsule":
                deps.add(Capsule(dep_dict["name"]))
            elif dep_dict["type"] == "SourceDep":
                deps.add(
                    SourceDep(
                        dep_dict["path"],
                        include_dirs=dep_dict["include_dirs"],
                        internal=dep_dict["internal"],
                    )
                )
            elif dep_dict["type"] == "HeaderDep":
                deps.add(
                    HeaderDep(
                        dep_dict["path"],
                        include_dirs=dep_dict["include_dirs"],
                        internal=dep_dict["internal"],
                    )
                )
        module.dependencies = deps

        return module


def deserialize_modules(data: dict[str, JsonDict], ctx: DeserMaps) -> dict[str, ModuleIR]:
    """Deserialize a collection of modules.

    The modules can contain dependencies on each other.

    Arguments:
        data: A dict containing the modules to deserialize.
        ctx: The deserialization maps to use and to populate.
             They are populated with information from the deserialized
             modules and as a precondition must have been populated by
             deserializing any dependencies of the modules being deserialized
             (outside of dependencies between the modules themselves).

    Returns a map containing the deserialized modules.
    """
    for mod in data.values():
        # First create ClassIRs for every class so that we can construct types and whatnot
        for cls in mod["classes"]:
            ir = ClassIR(cls["name"], cls["module_name"])
            assert ir.fullname not in ctx.classes, "Class %s already in map" % ir.fullname
            ctx.classes[ir.fullname] = ir

    for mod in data.values():
        # Then deserialize all of the functions so that methods are available
        # to the class deserialization.
        for method in mod["functions"]:
            func = FuncIR.deserialize(method, ctx)
            assert func.decl.id not in ctx.functions, (
                "Method %s already in map" % func.decl.fullname
            )
            ctx.functions[func.decl.id] = func

    return {k: ModuleIR.deserialize(v, ctx) for k, v in data.items()}


# ModulesIRs should also always be an *OrderedDict*, but if we
# declared it that way we would need to put it in quotes everywhere...
ModuleIRs = dict[str, ModuleIR]
