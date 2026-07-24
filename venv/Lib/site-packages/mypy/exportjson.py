"""Tool to convert binary mypy cache files (.ff) to JSON (.ff.json).

Usage:
   python -m mypy.exportjson .mypy_cache/.../my_module.data.ff

The idea is to make caches introspectable once we've switched to a binary
cache format and removed support for the older JSON cache format.

This is primarily to support existing use cases that need to inspect
cache files, and to support debugging mypy caching issues. This means that
this doesn't necessarily need to be kept 1:1 up to date with changes in the
binary cache format (to simplify maintenance -- we don't want this to slow
down mypy development).
"""

import argparse
import json
import sys
from typing import Any, TypeAlias as _TypeAlias

from librt.internal import ReadBuffer, cache_version

from mypy.cache import CACHE_VERSION, CacheMeta
from mypy.nodes import (
    FUNCBASE_FLAGS,
    FUNCDEF_FLAGS,
    VAR_FLAGS,
    ClassDef,
    DataclassTransformSpec,
    Decorator,
    FuncDef,
    MypyFile,
    OverloadedFuncDef,
    OverloadPart,
    ParamSpecExpr,
    SymbolNode,
    SymbolTable,
    SymbolTableNode,
    TypeAlias,
    TypeInfo,
    TypeVarExpr,
    TypeVarTupleExpr,
    Var,
    get_flags,
    node_kinds,
)
from mypy.types import (
    NOT_READY,
    AnyType,
    CallableType,
    ExtraAttrs,
    Instance,
    LiteralType,
    NoneType,
    Overloaded,
    Parameters,
    ParamSpecType,
    TupleType,
    Type,
    TypeAliasType,
    TypedDictType,
    TypeType,
    TypeVarTupleType,
    TypeVarType,
    UnboundType,
    UninhabitedType,
    UnionType,
    UnpackType,
    get_proper_type,
)

Json: _TypeAlias = dict[str, Any] | str


class Config:
    def __init__(self, *, implicit_names: bool = True) -> None:
        self.implicit_names = implicit_names


def convert_binary_cache_to_json(data: bytes, *, implicit_names: bool = True) -> Json:
    tree = MypyFile.read(ReadBuffer(data))
    return convert_mypy_file_to_json(tree, Config(implicit_names=implicit_names))


def convert_mypy_file_to_json(self: MypyFile, cfg: Config) -> Json:
    return {
        ".class": "MypyFile",
        "_fullname": self._fullname,
        "names": convert_symbol_table(self.names, cfg),
        "is_stub": self.is_stub,
        "path": self.path,
        "is_partial_stub_package": self.is_partial_stub_package,
        "future_import_flags": sorted(self.future_import_flags),
    }


def convert_symbol_table(self: SymbolTable, cfg: Config) -> Json:
    data: dict[str, Any] = {".class": "SymbolTable"}
    for key, value in self.items():
        # Skip __builtins__: it's a reference to the builtins
        # module that gets added to every module by
        # SemanticAnalyzerPass2.visit_file(), but it shouldn't be
        # accessed by users of the module.
        if key == "__builtins__" or value.no_serialize:
            continue
        if not cfg.implicit_names and key in {
            "__spec__",
            "__package__",
            "__file__",
            "__doc__",
            "__annotations__",
            "__name__",
        }:
            continue
        data[key] = convert_symbol_table_node(value, cfg)
    return data


def convert_symbol_table_node(self: SymbolTableNode, cfg: Config) -> Json:
    data: dict[str, Any] = {".class": "SymbolTableNode", "kind": node_kinds[self.kind]}
    if self.module_hidden:
        data["module_hidden"] = True
    if not self.module_public:
        data["module_public"] = False
    if self.implicit:
        data["implicit"] = True
    if self.plugin_generated:
        data["plugin_generated"] = True
    if self.cross_ref:
        data["cross_ref"] = self.cross_ref
    elif self.node is not None:
        data["node"] = convert_symbol_node(self.node, cfg)
    return data


def convert_symbol_node(self: SymbolNode, cfg: Config) -> Json:
    if isinstance(self, FuncDef):
        return convert_func_def(self)
    elif isinstance(self, OverloadedFuncDef):
        return convert_overloaded_func_def(self)
    elif isinstance(self, Decorator):
        return convert_decorator(self)
    elif isinstance(self, Var):
        return convert_var(self)
    elif isinstance(self, TypeInfo):
        return convert_type_info(self, cfg)
    elif isinstance(self, TypeAlias):
        return convert_type_alias(self)
    elif isinstance(self, TypeVarExpr):
        return convert_type_var_expr(self)
    elif isinstance(self, ParamSpecExpr):
        return convert_param_spec_expr(self)
    elif isinstance(self, TypeVarTupleExpr):
        return convert_type_var_tuple_expr(self)
    return {"ERROR": f"{type(self)!r} unrecognized"}


def convert_func_def(self: FuncDef) -> Json:
    return {
        ".class": "FuncDef",
        "name": self._name,
        "fullname": self._fullname,
        "arg_names": self.arg_names,
        "arg_kinds": [int(x.value) for x in self.arg_kinds],
        "type": None if self.type is None else convert_type(self.type),
        "flags": get_flags(self, FUNCDEF_FLAGS),
        "abstract_status": self.abstract_status,
        # TODO: Do we need expanded, original_def?
        "dataclass_transform_spec": (
            None
            if self.dataclass_transform_spec is None
            else convert_dataclass_transform_spec(self.dataclass_transform_spec)
        ),
        "deprecated": self.deprecated,
        "original_first_arg": self.original_first_arg,
    }


def convert_dataclass_transform_spec(self: DataclassTransformSpec) -> Json:
    return {
        "eq_default": self.eq_default,
        "order_default": self.order_default,
        "kw_only_default": self.kw_only_default,
        "frozen_default": self.frozen_default,
        "field_specifiers": list(self.field_specifiers),
    }


def convert_overloaded_func_def(self: OverloadedFuncDef) -> Json:
    return {
        ".class": "OverloadedFuncDef",
        "items": [convert_overload_part(i) for i in self.items],
        "type": None if self.type is None else convert_type(self.type),
        "fullname": self._fullname,
        "impl": None if self.impl is None else convert_overload_part(self.impl),
        "flags": get_flags(self, FUNCBASE_FLAGS),
        "deprecated": self.deprecated,
        "setter_index": self.setter_index,
    }


def convert_overload_part(self: OverloadPart) -> Json:
    if isinstance(self, FuncDef):
        return convert_func_def(self)
    else:
        return convert_decorator(self)


def convert_decorator(self: Decorator) -> Json:
    return {
        ".class": "Decorator",
        "func": convert_func_def(self.func),
        "var": convert_var(self.var),
        "is_overload": self.is_overload,
    }


def convert_var(self: Var) -> Json:
    data: dict[str, Any] = {
        ".class": "Var",
        "name": self._name,
        "fullname": self._fullname,
        "type": None if self.type is None else convert_type(self.type),
        "setter_type": None if self.setter_type is None else convert_type(self.setter_type),
        "flags": get_flags(self, VAR_FLAGS),
    }
    if self.final_value is not None:
        data["final_value"] = self.final_value
    return data


def convert_type_info(self: TypeInfo, cfg: Config) -> Json:
    data = {
        ".class": "TypeInfo",
        "module_name": self.module_name,
        "fullname": self.fullname,
        "names": convert_symbol_table(self.names, cfg),
        "defn": convert_class_def(self.defn),
        "abstract_attributes": self.abstract_attributes,
        "type_vars": self.type_vars,
        "has_param_spec_type": self.has_param_spec_type,
        "bases": [convert_type(b) for b in self.bases],
        "mro": self._mro_refs,
        "_promote": [convert_type(p) for p in self._promote],
        "alt_promote": None if self.alt_promote is None else convert_type(self.alt_promote),
        "declared_metaclass": (
            None if self.declared_metaclass is None else convert_type(self.declared_metaclass)
        ),
        "metaclass_type": (
            None if self.metaclass_type is None else convert_type(self.metaclass_type)
        ),
        "tuple_type": None if self.tuple_type is None else convert_type(self.tuple_type),
        "typeddict_type": (
            None if self.typeddict_type is None else convert_typeddict_type(self.typeddict_type)
        ),
        "flags": get_flags(self, TypeInfo.FLAGS),
        "metadata": self.metadata,
        "slots": sorted(self.slots) if self.slots is not None else None,
        "deletable_attributes": self.deletable_attributes,
        "self_type": convert_type(self.self_type) if self.self_type is not None else None,
        "dataclass_transform_spec": (
            convert_dataclass_transform_spec(self.dataclass_transform_spec)
            if self.dataclass_transform_spec is not None
            else None
        ),
        "deprecated": self.deprecated,
    }
    return data


def convert_class_def(self: ClassDef) -> Json:
    return {
        ".class": "ClassDef",
        "name": self.name,
        "fullname": self.fullname,
        "type_vars": [convert_type(v) for v in self.type_vars],
    }


def convert_type_alias(self: TypeAlias) -> Json:
    data: Json = {
        ".class": "TypeAlias",
        "fullname": self._fullname,
        "module": self.module,
        "target": convert_type(self.target),
        "alias_tvars": [convert_type(v) for v in self.alias_tvars],
        "no_args": self.no_args,
        "normalized": self.normalized,
        "python_3_12_type_alias": self.python_3_12_type_alias,
    }
    return data


def convert_type_var_expr(self: TypeVarExpr) -> Json:
    return {
        ".class": "TypeVarExpr",
        "name": self._name,
        "fullname": self._fullname,
        "values": [convert_type(t) for t in self.values],
        "upper_bound": convert_type(self.upper_bound),
        "default": convert_type(self.default),
        "variance": self.variance,
    }


def convert_param_spec_expr(self: ParamSpecExpr) -> Json:
    return {
        ".class": "ParamSpecExpr",
        "name": self._name,
        "fullname": self._fullname,
        "upper_bound": convert_type(self.upper_bound),
        "default": convert_type(self.default),
        "variance": self.variance,
    }


def convert_type_var_tuple_expr(self: TypeVarTupleExpr) -> Json:
    return {
        ".class": "TypeVarTupleExpr",
        "name": self._name,
        "fullname": self._fullname,
        "upper_bound": convert_type(self.upper_bound),
        "tuple_fallback": convert_type(self.tuple_fallback),
        "default": convert_type(self.default),
        "variance": self.variance,
    }


def convert_type(typ: Type) -> Json:
    if type(typ) is TypeAliasType:
        return convert_type_alias_type(typ)
    typ = get_proper_type(typ)
    if isinstance(typ, Instance):
        return convert_instance(typ)
    elif isinstance(typ, AnyType):
        return convert_any_type(typ)
    elif isinstance(typ, NoneType):
        return convert_none_type(typ)
    elif isinstance(typ, UnionType):
        return convert_union_type(typ)
    elif isinstance(typ, TupleType):
        return convert_tuple_type(typ)
    elif isinstance(typ, CallableType):
        return convert_callable_type(typ)
    elif isinstance(typ, Overloaded):
        return convert_overloaded(typ)
    elif isinstance(typ, LiteralType):
        return convert_literal_type(typ)
    elif isinstance(typ, TypeVarType):
        return convert_type_var_type(typ)
    elif isinstance(typ, TypeType):
        return convert_type_type(typ)
    elif isinstance(typ, UninhabitedType):
        return convert_uninhabited_type(typ)
    elif isinstance(typ, UnpackType):
        return convert_unpack_type(typ)
    elif isinstance(typ, ParamSpecType):
        return convert_param_spec_type(typ)
    elif isinstance(typ, TypeVarTupleType):
        return convert_type_var_tuple_type(typ)
    elif isinstance(typ, Parameters):
        return convert_parameters(typ)
    elif isinstance(typ, TypedDictType):
        return convert_typeddict_type(typ)
    elif isinstance(typ, UnboundType):
        return convert_unbound_type(typ)
    return {"ERROR": f"{type(typ)!r} unrecognized"}


def convert_instance(self: Instance) -> Json:
    ready = self.type is not NOT_READY
    if not self.args and not self.last_known_value and not self.extra_attrs:
        if ready:
            return self.type.fullname
        elif self.type_ref:
            return self.type_ref

    data: dict[str, Any] = {
        ".class": "Instance",
        "type_ref": self.type.fullname if ready else self.type_ref,
        "args": [convert_type(arg) for arg in self.args],
    }
    if self.last_known_value is not None:
        data["last_known_value"] = convert_type(self.last_known_value)
    data["extra_attrs"] = convert_extra_attrs(self.extra_attrs) if self.extra_attrs else None
    return data


def convert_extra_attrs(self: ExtraAttrs) -> Json:
    return {
        ".class": "ExtraAttrs",
        "attrs": {k: convert_type(v) for k, v in self.attrs.items()},
        "immutable": sorted(self.immutable),
        "mod_name": self.mod_name,
    }


def convert_type_alias_type(self: TypeAliasType) -> Json:
    data: Json = {
        ".class": "TypeAliasType",
        "type_ref": self.type_ref,
        "args": [convert_type(arg) for arg in self.args],
    }
    return data


def convert_any_type(self: AnyType) -> Json:
    return {
        ".class": "AnyType",
        "type_of_any": self.type_of_any,
        "source_any": convert_type(self.source_any) if self.source_any is not None else None,
        "missing_import_name": self.missing_import_name,
    }


def convert_none_type(self: NoneType) -> Json:
    return {".class": "NoneType"}


def convert_union_type(self: UnionType) -> Json:
    return {
        ".class": "UnionType",
        "items": [convert_type(t) for t in self.items],
        "uses_pep604_syntax": self.uses_pep604_syntax,
    }


def convert_tuple_type(self: TupleType) -> Json:
    return {
        ".class": "TupleType",
        "items": [convert_type(t) for t in self.items],
        "partial_fallback": convert_type(self.partial_fallback),
        "implicit": self.implicit,
    }


def convert_literal_type(self: LiteralType) -> Json:
    return {".class": "LiteralType", "value": self.value, "fallback": convert_type(self.fallback)}


def convert_type_var_type(self: TypeVarType) -> Json:
    assert not self.id.is_meta_var()
    return {
        ".class": "TypeVarType",
        "name": self.name,
        "fullname": self.fullname,
        "id": self.id.raw_id,
        "namespace": self.id.namespace,
        "values": [convert_type(v) for v in self.values],
        "upper_bound": convert_type(self.upper_bound),
        "default": convert_type(self.default),
        "variance": self.variance,
    }


def convert_callable_type(self: CallableType) -> Json:
    return {
        ".class": "CallableType",
        "arg_types": [convert_type(t) for t in self.arg_types],
        "arg_kinds": [int(x.value) for x in self.arg_kinds],
        "arg_names": self.arg_names,
        "ret_type": convert_type(self.ret_type),
        "fallback": convert_type(self.fallback),
        "name": self.name,
        # We don't serialize the definition (only used for error messages).
        "variables": [convert_type(v) for v in self.variables],
        "is_ellipsis_args": self.is_ellipsis_args,
        "implicit": self.implicit,
        "is_bound": self.is_bound,
        "type_guard": convert_type(self.type_guard) if self.type_guard is not None else None,
        "type_is": convert_type(self.type_is) if self.type_is is not None else None,
        "from_concatenate": self.from_concatenate,
        "imprecise_arg_kinds": self.imprecise_arg_kinds,
        "unpack_kwargs": self.unpack_kwargs,
    }


def convert_overloaded(self: Overloaded) -> Json:
    return {".class": "Overloaded", "items": [convert_type(t) for t in self.items]}


def convert_type_type(self: TypeType) -> Json:
    return {".class": "TypeType", "item": convert_type(self.item)}


def convert_uninhabited_type(self: UninhabitedType) -> Json:
    return {".class": "UninhabitedType"}


def convert_unpack_type(self: UnpackType) -> Json:
    return {".class": "UnpackType", "type": convert_type(self.type)}


def convert_param_spec_type(self: ParamSpecType) -> Json:
    assert not self.id.is_meta_var()
    return {
        ".class": "ParamSpecType",
        "name": self.name,
        "fullname": self.fullname,
        "id": self.id.raw_id,
        "namespace": self.id.namespace,
        "flavor": self.flavor,
        "upper_bound": convert_type(self.upper_bound),
        "default": convert_type(self.default),
        "prefix": convert_type(self.prefix),
    }


def convert_type_var_tuple_type(self: TypeVarTupleType) -> Json:
    assert not self.id.is_meta_var()
    return {
        ".class": "TypeVarTupleType",
        "name": self.name,
        "fullname": self.fullname,
        "id": self.id.raw_id,
        "namespace": self.id.namespace,
        "upper_bound": convert_type(self.upper_bound),
        "tuple_fallback": convert_type(self.tuple_fallback),
        "default": convert_type(self.default),
        "min_len": self.min_len,
    }


def convert_parameters(self: Parameters) -> Json:
    return {
        ".class": "Parameters",
        "arg_types": [convert_type(t) for t in self.arg_types],
        "arg_kinds": [int(x.value) for x in self.arg_kinds],
        "arg_names": self.arg_names,
        "variables": [convert_type(tv) for tv in self.variables],
        "imprecise_arg_kinds": self.imprecise_arg_kinds,
    }


def convert_typeddict_type(self: TypedDictType) -> Json:
    return {
        ".class": "TypedDictType",
        "items": [[n, convert_type(t)] for (n, t) in self.items.items()],
        "required_keys": sorted(self.required_keys),
        "readonly_keys": sorted(self.readonly_keys),
        "fallback": convert_type(self.fallback),
    }


def convert_unbound_type(self: UnboundType) -> Json:
    return {
        ".class": "UnboundType",
        "name": self.name,
        "args": [convert_type(a) for a in self.args],
        "expr": self.original_str_expr,
        "expr_fallback": self.original_str_fallback,
    }


def convert_binary_cache_meta_to_json(data: bytes, data_file: str) -> Json:
    assert (
        data[0] == cache_version() and data[1] == CACHE_VERSION
    ), "Cache file created by an incompatible mypy version"
    meta = CacheMeta.read(ReadBuffer(data[2:]), data_file)
    assert meta is not None, f"Error reading meta cache file associated with {data_file}"
    return {
        "id": meta.id,
        "path": meta.path,
        "mtime": meta.mtime,
        "size": meta.size,
        "hash": meta.hash,
        "data_mtime": meta.data_mtime,
        "dependencies": meta.dependencies,
        "suppressed": meta.suppressed,
        "options": meta.options,
        "dep_prios": meta.dep_prios,
        "dep_lines": meta.dep_lines,
        "dep_hashes": [dep.hex() for dep in meta.dep_hashes],
        "interface_hash": meta.interface_hash.hex(),
        "version_id": meta.version_id,
        "ignore_all": meta.ignore_all,
        "plugin_data": meta.plugin_data,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert binary cache files to JSON. "
        "Create files in the same directory with extra .json extension."
    )
    parser.add_argument(
        "path", nargs="+", help="mypy cache data file to convert (.data.ff extension)"
    )
    args = parser.parse_args()
    fnams: list[str] = args.path
    for fnam in fnams:
        if fnam.endswith(".data.ff"):
            is_data = True
        elif fnam.endswith(".meta.ff"):
            is_data = False
        else:
            sys.exit(f"error: Expected .data.ff or .meta.ff extension, but got {fnam}")
        with open(fnam, "rb") as f:
            data = f.read()
        if is_data:
            json_data = convert_binary_cache_to_json(data)
        else:
            data_file = fnam.removesuffix(".meta.ff") + ".data.ff"
            json_data = convert_binary_cache_meta_to_json(data, data_file)
        new_fnam = fnam + ".json"
        with open(new_fnam, "w") as f:
            json.dump(json_data, f)
        print(f"{fnam} -> {new_fnam}")


if __name__ == "__main__":
    main()
