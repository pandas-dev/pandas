from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pytest

from mypy import build
from mypy.options import Options
from mypyc.build import construct_groups
from mypyc.codegen import emitmodule
from mypyc.errors import Errors
from mypyc.irbuild.mapper import Mapper
from mypyc.options import CompilerOptions


class FakeSCC:
    def __init__(self, mod_ids: list[str]) -> None:
        self.mod_ids = mod_ids


class TestEmitModule(unittest.TestCase):
    def test_compile_modules_to_ir_orders_scc_members_deterministically(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir, pytest.MonkeyPatch.context() as monkeypatch:
            tmp_path = Path(tmp_dir)
            a_py = tmp_path / "a.py"
            b_py = tmp_path / "b.py"
            a_py.write_text("import b\n\nclass A: pass\nclass C(A): pass\n", encoding="utf-8")
            b_py.write_text(
                "import a\n\nclass B(a.A): pass\nclass D(a.A): pass\n", encoding="utf-8"
            )

            sources = [
                build.BuildSource(str(a_py), "a", None),
                build.BuildSource(str(b_py), "b", None),
            ]
            options = Options()
            options.preserve_asts = True
            options.mypy_path = [str(tmp_path)]
            options.cache_dir = str(tmp_path / ".mypy_cache")
            for source in sources:
                options.per_module_options.setdefault(source.module, {})["mypyc"] = True

            compiler_options = CompilerOptions(strict_traceback_checks=True)
            groups = construct_groups(
                sources, False, use_shared_lib=True, group_name_override=None
            )
            result = emitmodule.parse_and_typecheck(sources, options, compiler_options, groups)
            try:
                group_map = {
                    source.module: lib_name for group, lib_name in groups for source in group
                }
                children_by_order = []
                for order in (["a", "b"], ["b", "a"]):
                    monkeypatch.setattr(
                        emitmodule,
                        "sorted_components",
                        lambda graph, order=order: [FakeSCC(order)],
                    )
                    mapper = Mapper(group_map)
                    errors = Errors(options)
                    modules = emitmodule.compile_modules_to_ir(
                        result, mapper, compiler_options, errors
                    )
                    assert errors.num_errors == 0, errors.new_messages()
                    classes = {
                        cl.fullname: cl for module in modules.values() for cl in module.classes
                    }
                    children = classes["a.A"].children
                    assert children is not None
                    children_by_order.append([child.fullname for child in children])

                assert children_by_order[1] == children_by_order[0]
            finally:
                result.manager.metastore.close()
