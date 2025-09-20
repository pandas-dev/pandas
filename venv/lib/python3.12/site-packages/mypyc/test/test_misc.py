from __future__ import annotations

import unittest

from mypyc.ir.ops import BasicBlock
from mypyc.ir.pprint import format_blocks, generate_names_for_ir
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.options import CompilerOptions


class TestMisc(unittest.TestCase):
    def test_debug_op(self) -> None:
        block = BasicBlock()
        builder = LowLevelIRBuilder(errors=None, options=CompilerOptions())
        builder.activate_block(block)
        builder.debug_print("foo")

        names = generate_names_for_ir([], [block])
        code = format_blocks([block], names, {})
        assert code[:-1] == ["L0:", "    r0 = 'foo'", "    CPyDebug_PrintObject(r0)"]
