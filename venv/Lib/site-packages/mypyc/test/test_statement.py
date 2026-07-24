from __future__ import annotations

import unittest
from typing import cast

from mypy.build import Graph
from mypy.nodes import Import, MypyFile
from mypy.options import Options
from mypyc.errors import Errors
from mypyc.irbuild.builder import IRBuilder
from mypyc.irbuild.mapper import Mapper
from mypyc.irbuild.prebuildvisitor import PreBuildVisitor
from mypyc.irbuild.statement import (
    IMPORT_NATIVE_ATTR,
    IMPORT_NATIVE_SUBMODULE,
    IMPORT_NON_NATIVE,
    classify_import_from,
    group_consecutive,
    import_globals_id_and_name,
    split_import_group_to_python_and_native,
)
from mypyc.irbuild.visitor import IRBuilderVisitor
from mypyc.options import CompilerOptions


def make_builder(
    *,
    module_name: str = "pkg.current",
    native_modules: set[str] | None = None,
    same_group_modules: set[str] | None = None,
    graph: set[str] | None = None,
) -> IRBuilder:
    native_modules = native_modules or set()
    same_group_modules = same_group_modules or set()
    group_map: dict[str, str | None] = {module_name: "current-group"}
    for module in native_modules:
        group_map[module] = "current-group" if module in same_group_modules else "other-group"

    errors = Errors(Options())
    current_file = MypyFile([], [])
    current_file._fullname = module_name
    pbv = PreBuildVisitor(errors, current_file, {}, {})
    builder = IRBuilder(
        module_name,
        {},
        cast(Graph, {name: object() for name in (graph or set())}),
        errors,
        Mapper(group_map),
        pbv,
        IRBuilderVisitor(),
        CompilerOptions(),
        {},
    )
    builder.set_module(module_name, module_name.replace(".", "/") + ".py")
    return builder


class TestStatementHelpers(unittest.TestCase):
    def test_import_globals_id_and_name_for_plain_import(self) -> None:
        assert import_globals_id_and_name("foo.bar", None) == ("foo", "foo")

    def test_import_globals_id_and_name_for_import_as(self) -> None:
        assert import_globals_id_and_name("foo.bar", "baz") == ("foo.bar", "baz")

    def test_split_import_group_to_python_and_native_preserves_runs(self) -> None:
        builder = make_builder(
            native_modules={"pkg.alpha", "pkg.beta", "pkg.gamma"},
            same_group_modules={"pkg.alpha", "pkg.beta", "pkg.gamma"},
        )
        group = [
            Import([("pkg.alpha", None), ("py_mod", None)]),
            Import([("pkg.beta", "beta_alias"), ("foreign.mod", None), ("pkg.gamma", None)]),
        ]
        group[0].line = 10
        group[1].line = 20

        result = split_import_group_to_python_and_native(builder, group)

        assert result == [
            ([("pkg.alpha", None, 10)], True),
            ([("py_mod", None, 10)], False),
            ([("pkg.beta", "beta_alias", 20)], True),
            ([("foreign.mod", None, 20)], False),
            ([("pkg.gamma", None, 20)], True),
        ]

    def test_group_consecutive_groups_by_kind_and_preserves_aliases(self) -> None:
        buckets = group_consecutive(
            [
                (IMPORT_NATIVE_SUBMODULE, "a", "a"),
                (IMPORT_NATIVE_SUBMODULE, "b", "b_alias"),
                (IMPORT_NON_NATIVE, "c", "c"),
                (IMPORT_NATIVE_ATTR, "d", "d_alias"),
                (IMPORT_NATIVE_ATTR, "e", "e"),
            ]
        )

        assert [(bucket.kind, bucket.names, bucket.as_names) for bucket in buckets] == [
            (IMPORT_NATIVE_SUBMODULE, ["a", "b"], ["a", "b_alias"]),
            (IMPORT_NON_NATIVE, ["c"], ["c"]),
            (IMPORT_NATIVE_ATTR, ["d", "e"], ["d_alias", "e"]),
        ]

    def test_classify_import_from_groups_consecutive_kinds(self) -> None:
        builder = make_builder(
            native_modules={"pkg.native_a", "pkg.native_b"},
            same_group_modules={"pkg.native_a", "pkg.native_b"},
            graph={"pkg.native_a", "pkg.native_b", "pkg.foreign_a", "pkg.foreign_b"},
        )

        buckets = classify_import_from(
            builder,
            "pkg",
            ["native_a", "native_b", "foreign_a", "foreign_b"],
            ["native_a", "native_b_alias", "foreign_a", "foreign_b_alias"],
            parent_is_native=True,
        )

        assert [(bucket.kind, bucket.names, bucket.as_names) for bucket in buckets] == [
            (IMPORT_NATIVE_SUBMODULE, ["native_a", "native_b"], ["native_a", "native_b_alias"]),
            (IMPORT_NON_NATIVE, ["foreign_a", "foreign_b"], ["foreign_a", "foreign_b_alias"]),
        ]

    def test_classify_import_from_treats_missing_name_under_native_parent_as_attr(self) -> None:
        builder = make_builder(graph={"pkg.foreign"})

        buckets = classify_import_from(
            builder,
            "pkg",
            ["attr_name", "foreign"],
            ["attr_alias", "foreign_alias"],
            parent_is_native=True,
        )

        assert [(bucket.kind, bucket.names, bucket.as_names) for bucket in buckets] == [
            (IMPORT_NATIVE_ATTR, ["attr_name"], ["attr_alias"]),
            (IMPORT_NON_NATIVE, ["foreign"], ["foreign_alias"]),
        ]

    def test_classify_import_from_without_native_parent_never_uses_native_attr(self) -> None:
        builder = make_builder()

        buckets = classify_import_from(
            builder, "pkg", ["attr_name"], ["attr_alias"], parent_is_native=False
        )

        assert [(bucket.kind, bucket.names, bucket.as_names) for bucket in buckets] == [
            (IMPORT_NON_NATIVE, ["attr_name"], ["attr_alias"])
        ]
