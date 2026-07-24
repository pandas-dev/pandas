"""This optional pass adds logging of various executed operations.

Some subset of the executed operations are logged to the mypyc_trace.txt file.

This is useful for performance analysis. For example, it's possible
to identify how frequently various primitive functions are called,
and in which code locations they are called.
"""

from __future__ import annotations

from typing import Final

from mypyc.ir.func_ir import FuncIR
from mypyc.ir.ops import (
    Box,
    Call,
    CallC,
    Cast,
    CString,
    DecRef,
    GetAttr,
    IncRef,
    LoadLiteral,
    LoadStatic,
    Op,
    PrimitiveOp,
    SetAttr,
    Unbox,
    Value,
)
from mypyc.ir.rtypes import none_rprimitive
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.options import CompilerOptions
from mypyc.primitives.misc_ops import log_trace_event
from mypyc.transform.ir_transform import IRTransform


def insert_event_trace_logging(fn: FuncIR, options: CompilerOptions) -> None:
    builder = LowLevelIRBuilder(None, options)
    transform = LogTraceEventTransform(builder, fn.decl.fullname)
    transform.transform_blocks(fn.blocks)
    fn.blocks = builder.blocks


def get_load_global_name(op: CallC) -> str | None:
    name = op.function_name
    if name == "CPyDict_GetItem":
        arg = op.args[0]
        if (
            isinstance(arg, LoadStatic)
            and arg.namespace == "static"
            and arg.identifier == "globals"
            and isinstance(op.args[1], LoadLiteral)
        ):
            return str(op.args[1].value)
    return None


# These primitives perform an implicit IncRef for the return value. Only some of the most common ones
# are included, and mostly ops that could be switched to use borrowing in some contexts.
primitives_that_inc_ref: Final = {
    "list_get_item_unsafe",
    "CPyList_GetItemShort",
    "CPyDict_GetWithNone",
    "CPyList_GetItem",
    "CPyDict_GetItem",
    "CPyList_PopLast",
}


class LogTraceEventTransform(IRTransform):
    def __init__(self, builder: LowLevelIRBuilder, fullname: str) -> None:
        super().__init__(builder)
        self.fullname = fullname.encode("utf-8")

    def visit_call(self, op: Call) -> Value:
        # TODO: Use different op name when constructing an instance
        return self.log(op, "call", op.fn.fullname)

    def visit_primitive_op(self, op: PrimitiveOp) -> Value:
        value = self.log(op, "primitive_op", op.desc.name)
        if op.desc.name in primitives_that_inc_ref:
            self.log_inc_ref(value)
        return value

    def visit_call_c(self, op: CallC) -> Value:
        if global_name := get_load_global_name(op):
            return self.log(op, "globals_dict_get_item", global_name)

        func_name = op.function_name
        if func_name == "PyObject_Vectorcall" and isinstance(op.args[0], CallC):
            if global_name := get_load_global_name(op.args[0]):
                return self.log(op, "python_call_global", global_name)
        elif func_name == "CPyObject_GetAttr" and isinstance(op.args[1], LoadLiteral):
            return self.log(op, "python_get_attr", str(op.args[1].value))
        elif func_name == "PyObject_VectorcallMethod" and isinstance(op.args[0], LoadLiteral):
            return self.log(op, "python_call_method", str(op.args[0].value))

        value = self.log(op, "call_c", func_name)
        if func_name in primitives_that_inc_ref:
            self.log_inc_ref(value)
        return value

    def visit_get_attr(self, op: GetAttr) -> Value:
        value = self.log(op, "get_attr", f"{op.class_type.name}.{op.attr}")
        if not op.is_borrowed and op.type.is_refcounted:
            self.log_inc_ref(op)
        return value

    def visit_set_attr(self, op: SetAttr) -> Value:
        name = "set_attr" if not op.is_init else "set_attr_init"
        return self.log(op, name, f"{op.class_type.name}.{op.attr}")

    def visit_box(self, op: Box) -> Value:
        if op.src.type is none_rprimitive:
            # Boxing 'None' is a very quick operation, so we don't log it.
            return self.add(op)
        else:
            return self.log(op, "box", str(op.src.type))

    def visit_unbox(self, op: Unbox) -> Value:
        return self.log(op, "unbox", str(op.type))

    def visit_cast(self, op: Cast) -> Value | None:
        value = self.log(op, "cast", str(op.type))
        if not op.is_borrowed:
            self.log_inc_ref(value)
        return value

    def visit_inc_ref(self, op: IncRef) -> Value:
        return self.log(op, "inc_ref", str(op.src.type))

    def visit_dec_ref(self, op: DecRef) -> Value:
        return self.log(op, "dec_ref", str(op.src.type))

    def log_inc_ref(self, value: Value) -> None:
        self.log_event("inc_ref", str(value.type), value.line)

    def log(self, op: Op, name: str, details: str) -> Value:
        self.log_event(name, details, op.line)
        return self.add(op)

    def log_event(self, name: str, details: str, line: int) -> None:
        if line >= 0:
            line_str = str(line)
        else:
            line_str = ""
        self.builder.primitive_op(
            log_trace_event,
            [
                CString(self.fullname),
                CString(line_str.encode("ascii")),
                CString(name.encode("utf-8")),
                CString(details.encode("utf-8")),
            ],
            line,
        )
