from __future__ import annotations

from collections.abc import Iterable

import mypy.types as types
from mypy.types import TypeVisitor


class TypeIndirectionVisitor(TypeVisitor[None]):
    """Returns all module references within a particular type."""

    def __init__(self) -> None:
        # Module references are collected here
        self.modules: set[str] = set()
        # User to avoid infinite recursion with recursive types
        self.seen_types: set[types.TypeAliasType | types.Instance] = set()

    def find_modules(self, typs: Iterable[types.Type]) -> set[str]:
        self.modules = set()
        self.seen_types = set()
        for typ in typs:
            self._visit(typ)
        return self.modules

    def _visit(self, typ: types.Type) -> None:
        # Note: instances are needed for `class str(Sequence[str]): ...`
        if (
            isinstance(typ, types.TypeAliasType)
            or isinstance(typ, types.ProperType)
            and isinstance(typ, types.Instance)
        ):
            # Avoid infinite recursion for recursive types.
            if typ in self.seen_types:
                return
            self.seen_types.add(typ)
        typ.accept(self)

    def _visit_type_tuple(self, typs: tuple[types.Type, ...]) -> None:
        # Micro-optimization: Specialized version of _visit for lists
        for typ in typs:
            if (
                isinstance(typ, types.TypeAliasType)
                or isinstance(typ, types.ProperType)
                and isinstance(typ, types.Instance)
            ):
                # Avoid infinite recursion for recursive types.
                if typ in self.seen_types:
                    continue
                self.seen_types.add(typ)
            typ.accept(self)

    def _visit_type_list(self, typs: list[types.Type]) -> None:
        # Micro-optimization: Specialized version of _visit for tuples
        for typ in typs:
            if (
                isinstance(typ, types.TypeAliasType)
                or isinstance(typ, types.ProperType)
                and isinstance(typ, types.Instance)
            ):
                # Avoid infinite recursion for recursive types.
                if typ in self.seen_types:
                    continue
                self.seen_types.add(typ)
            typ.accept(self)

    def visit_unbound_type(self, t: types.UnboundType) -> None:
        self._visit_type_tuple(t.args)

    def visit_any(self, t: types.AnyType) -> None:
        pass

    def visit_none_type(self, t: types.NoneType) -> None:
        pass

    def visit_uninhabited_type(self, t: types.UninhabitedType) -> None:
        pass

    def visit_erased_type(self, t: types.ErasedType) -> None:
        pass

    def visit_deleted_type(self, t: types.DeletedType) -> None:
        pass

    def visit_type_var(self, t: types.TypeVarType) -> None:
        self._visit_type_list(t.values)
        self._visit(t.upper_bound)
        self._visit(t.default)

    def visit_param_spec(self, t: types.ParamSpecType) -> None:
        self._visit(t.upper_bound)
        self._visit(t.default)
        self._visit(t.prefix)

    def visit_type_var_tuple(self, t: types.TypeVarTupleType) -> None:
        self._visit(t.upper_bound)
        self._visit(t.default)

    def visit_unpack_type(self, t: types.UnpackType) -> None:
        t.type.accept(self)

    def visit_parameters(self, t: types.Parameters) -> None:
        self._visit_type_list(t.arg_types)

    def visit_instance(self, t: types.Instance) -> None:
        # Instance is named, record its definition and continue digging into
        # components that constitute semantic meaning of this type: bases, metaclass,
        # tuple type, and typeddict type.
        # Note: we cannot simply record the MRO, in case an intermediate base contains
        # a reference to type alias, this affects meaning of map_instance_to_supertype(),
        # see e.g. testDoubleReexportGenericUpdated.
        self._visit_type_tuple(t.args)
        if t.type:
            # Important optimization: instead of simply recording the definition and
            # recursing into bases, record the MRO and only traverse generic bases.
            for s in t.type.mro:
                self.modules.add(s.module_name)
                for base in s.bases:
                    if base.args:
                        self._visit_type_tuple(base.args)
            if t.type.metaclass_type:
                self._visit(t.type.metaclass_type)
            if t.type.typeddict_type:
                self._visit(t.type.typeddict_type)
            if t.type.tuple_type:
                self._visit(t.type.tuple_type)
            if t.type.is_protocol:
                # For protocols, member types constitute the semantic meaning of the type.
                # TODO: this doesn't cover some edge cases, like setter types and exotic nodes.
                for m in t.type.protocol_members:
                    node = t.type.names.get(m)
                    if node and node.type:
                        self._visit(node.type)

    def visit_callable_type(self, t: types.CallableType) -> None:
        self._visit_type_list(t.arg_types)
        self._visit(t.ret_type)
        self._visit_type_tuple(t.variables)

    def visit_overloaded(self, t: types.Overloaded) -> None:
        for item in t.items:
            self._visit(item)
        self._visit(t.fallback)

    def visit_tuple_type(self, t: types.TupleType) -> None:
        self._visit_type_list(t.items)
        self._visit(t.partial_fallback)

    def visit_typeddict_type(self, t: types.TypedDictType) -> None:
        self._visit_type_list(list(t.items.values()))
        self._visit(t.fallback)

    def visit_literal_type(self, t: types.LiteralType) -> None:
        self._visit(t.fallback)

    def visit_union_type(self, t: types.UnionType) -> None:
        self._visit_type_list(t.items)

    def visit_partial_type(self, t: types.PartialType) -> None:
        pass

    def visit_type_type(self, t: types.TypeType) -> None:
        self._visit(t.item)

    def visit_type_alias_type(self, t: types.TypeAliasType) -> None:
        # Type alias is named, record its definition and continue digging into
        # components that constitute semantic meaning of this type: target and args.
        if t.alias:
            self.modules.add(t.alias.module)
            self._visit(t.alias.target)
        self._visit_type_list(t.args)
