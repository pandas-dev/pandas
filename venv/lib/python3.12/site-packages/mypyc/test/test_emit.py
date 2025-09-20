from __future__ import annotations

import unittest

from mypyc.codegen.emit import Emitter, EmitterContext
from mypyc.common import HAVE_IMMORTAL
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.ops import BasicBlock, Register, Value
from mypyc.ir.rtypes import (
    RInstance,
    RTuple,
    RUnion,
    bool_rprimitive,
    int_rprimitive,
    list_rprimitive,
    none_rprimitive,
    object_rprimitive,
    str_rprimitive,
)
from mypyc.irbuild.vtable import compute_vtable
from mypyc.namegen import NameGenerator


class TestEmitter(unittest.TestCase):
    def setUp(self) -> None:
        self.n = Register(int_rprimitive, "n")
        self.context = EmitterContext(NameGenerator([["mod"]]))
        self.emitter = Emitter(self.context, {})

        ir = ClassIR("A", "mod")
        compute_vtable(ir)
        ir.mro = [ir]
        self.instance_a = RInstance(ir)

    def test_label(self) -> None:
        assert self.emitter.label(BasicBlock(4)) == "CPyL4"

    def test_reg(self) -> None:
        names: dict[Value, str] = {self.n: "n"}
        emitter = Emitter(self.context, names)
        assert emitter.reg(self.n) == "cpy_r_n"

    def test_object_annotation(self) -> None:
        assert self.emitter.object_annotation("hello, world", "line;") == " /* 'hello, world' */"
        assert (
            self.emitter.object_annotation(list(range(30)), "line;")
            == """\
 /* [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
         23, 24, 25, 26, 27, 28, 29] */"""
        )

    def test_emit_line(self) -> None:
        emitter = self.emitter
        emitter.emit_line("line;")
        emitter.emit_line("a {")
        emitter.emit_line("f();")
        emitter.emit_line("}")
        assert emitter.fragments == ["line;\n", "a {\n", "    f();\n", "}\n"]
        emitter = Emitter(self.context, {})
        emitter.emit_line("CPyStatics[0];", ann="hello, world")
        emitter.emit_line("CPyStatics[1];", ann=list(range(30)))
        assert emitter.fragments[0] == "CPyStatics[0]; /* 'hello, world' */\n"
        assert (
            emitter.fragments[1]
            == """\
CPyStatics[1]; /* [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 22, 23, 24, 25, 26, 27, 28, 29] */\n"""
        )

    def test_emit_undefined_value_for_simple_type(self) -> None:
        emitter = self.emitter
        assert emitter.c_undefined_value(int_rprimitive) == "CPY_INT_TAG"
        assert emitter.c_undefined_value(str_rprimitive) == "NULL"
        assert emitter.c_undefined_value(bool_rprimitive) == "2"

    def test_emit_undefined_value_for_tuple(self) -> None:
        emitter = self.emitter
        assert (
            emitter.c_undefined_value(RTuple([str_rprimitive, int_rprimitive, bool_rprimitive]))
            == "(tuple_T3OIC) { NULL, CPY_INT_TAG, 2 }"
        )
        assert emitter.c_undefined_value(RTuple([str_rprimitive])) == "(tuple_T1O) { NULL }"
        assert (
            emitter.c_undefined_value(RTuple([RTuple([str_rprimitive]), bool_rprimitive]))
            == "(tuple_T2T1OC) { { NULL }, 2 }"
        )

    def test_emit_inc_ref_object(self) -> None:
        self.emitter.emit_inc_ref("x", object_rprimitive)
        self.assert_output("CPy_INCREF(x);\n")

    def test_emit_inc_ref_int(self) -> None:
        self.emitter.emit_inc_ref("x", int_rprimitive)
        self.assert_output("CPyTagged_INCREF(x);\n")

    def test_emit_inc_ref_rare(self) -> None:
        self.emitter.emit_inc_ref("x", object_rprimitive, rare=True)
        self.assert_output("CPy_INCREF(x);\n")
        self.emitter.emit_inc_ref("x", int_rprimitive, rare=True)
        self.assert_output("CPyTagged_IncRef(x);\n")

    def test_emit_inc_ref_list(self) -> None:
        self.emitter.emit_inc_ref("x", list_rprimitive)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_INCREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_INCREF(x);\n")

    def test_emit_inc_ref_instance(self) -> None:
        self.emitter.emit_inc_ref("x", self.instance_a)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_INCREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_INCREF(x);\n")

    def test_emit_inc_ref_optional(self) -> None:
        optional = RUnion([self.instance_a, none_rprimitive])
        self.emitter.emit_inc_ref("o", optional)
        self.assert_output("CPy_INCREF(o);\n")

    def test_emit_dec_ref_object(self) -> None:
        self.emitter.emit_dec_ref("x", object_rprimitive)
        self.assert_output("CPy_DECREF(x);\n")
        self.emitter.emit_dec_ref("x", object_rprimitive, is_xdec=True)
        self.assert_output("CPy_XDECREF(x);\n")

    def test_emit_dec_ref_int(self) -> None:
        self.emitter.emit_dec_ref("x", int_rprimitive)
        self.assert_output("CPyTagged_DECREF(x);\n")
        self.emitter.emit_dec_ref("x", int_rprimitive, is_xdec=True)
        self.assert_output("CPyTagged_XDECREF(x);\n")

    def test_emit_dec_ref_rare(self) -> None:
        self.emitter.emit_dec_ref("x", object_rprimitive, rare=True)
        self.assert_output("CPy_DecRef(x);\n")
        self.emitter.emit_dec_ref("x", int_rprimitive, rare=True)
        self.assert_output("CPyTagged_DecRef(x);\n")

    def test_emit_dec_ref_list(self) -> None:
        self.emitter.emit_dec_ref("x", list_rprimitive)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_DECREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_DECREF(x);\n")
        self.emitter.emit_dec_ref("x", list_rprimitive, is_xdec=True)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_XDECREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_XDECREF(x);\n")

    def test_emit_dec_ref_instance(self) -> None:
        self.emitter.emit_dec_ref("x", self.instance_a)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_DECREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_DECREF(x);\n")
        self.emitter.emit_dec_ref("x", self.instance_a, is_xdec=True)
        if HAVE_IMMORTAL:
            self.assert_output("CPy_XDECREF_NO_IMM(x);\n")
        else:
            self.assert_output("CPy_XDECREF(x);\n")

    def test_emit_dec_ref_optional(self) -> None:
        optional = RUnion([self.instance_a, none_rprimitive])
        self.emitter.emit_dec_ref("o", optional)
        self.assert_output("CPy_DECREF(o);\n")

    def assert_output(self, expected: str) -> None:
        assert "".join(self.emitter.fragments) == expected
        self.emitter.fragments = []
