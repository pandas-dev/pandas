"""Code generation for native function bodies."""

from __future__ import annotations

from typing import Final

from mypyc.analysis.blockfreq import frequently_executed_blocks
from mypyc.codegen.emit import (
    DEBUG_ERRORS,
    PREFIX_MAP,
    Emitter,
    TracebackAndGotoHandler,
    c_array_initializer,
)
from mypyc.common import GENERATOR_ATTRIBUTE_PREFIX, HAVE_IMMORTAL, NATIVE_PREFIX, REG_PREFIX
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.func_ir import FUNC_CLASSMETHOD, FUNC_STATICMETHOD, FuncDecl, FuncIR, all_values
from mypyc.ir.ops import (
    ERR_FALSE,
    NAMESPACE_TYPE,
    Assign,
    AssignMulti,
    BasicBlock,
    Box,
    Branch,
    Call,
    CallC,
    Cast,
    ComparisonOp,
    ControlOp,
    CString,
    DecRef,
    Extend,
    Float,
    FloatComparisonOp,
    FloatNeg,
    FloatOp,
    GetAttr,
    GetElement,
    GetElementPtr,
    Goto,
    IncRef,
    InitStatic,
    Integer,
    IntOp,
    KeepAlive,
    LoadAddress,
    LoadErrorValue,
    LoadGlobal,
    LoadLiteral,
    LoadMem,
    LoadStatic,
    MethodCall,
    Op,
    OpVisitor,
    PrimitiveOp,
    RaiseStandardError,
    Register,
    Return,
    SetAttr,
    SetElement,
    SetMem,
    Truncate,
    TupleGet,
    TupleSet,
    Unborrow,
    Unbox,
    Undef,
    Unreachable,
    Value,
)
from mypyc.ir.pprint import generate_names_for_ir
from mypyc.ir.rtypes import (
    RArray,
    RInstance,
    RStruct,
    RTuple,
    RType,
    RVec,
    is_bool_or_bit_rprimitive,
    is_int32_rprimitive,
    is_int64_rprimitive,
    is_int_rprimitive,
    is_none_rprimitive,
    is_pointer_rprimitive,
    is_tagged,
)

VEC_ITEMS_C_TYPE: Final = {
    "VecI64": "int64_t *",
    "VecI32": "int32_t *",
    "VecI16": "int16_t *",
    "VecU8": "uint8_t *",
    "VecFloat": "double *",
    "VecBool": "char *",
    "VecT": "PyObject **",
    "VecNested": "VecNestedBufItem *",
    "VecNestedBufItem": "void *",
}


def native_function_type(fn: FuncIR, emitter: Emitter) -> str:
    return native_function_type_from_decl(fn.decl, emitter)


def native_function_type_from_decl(decl: FuncDecl, emitter: Emitter) -> str:
    args = ", ".join(emitter.ctype(arg.type) for arg in decl.sig.args) or "void"
    ret = emitter.ctype(decl.sig.ret_type)
    return f"{ret} (*)({args})"


def native_function_header(fn: FuncDecl, emitter: Emitter) -> str:
    args = []
    for arg in fn.sig.args:
        args.append(f"{emitter.ctype_spaced(arg.type)}{REG_PREFIX}{arg.name}")

    return "{ret_type}{name}({args})".format(
        ret_type=emitter.ctype_spaced(fn.sig.ret_type),
        name=emitter.native_function_name(fn),
        args=", ".join(args) or "void",
    )


def generate_native_function(
    fn: FuncIR, emitter: Emitter, source_path: str, module_name: str
) -> None:
    declarations = Emitter(emitter.context)
    names = generate_names_for_ir(fn.arg_regs, fn.blocks)
    body = Emitter(emitter.context, names)
    visitor = FunctionEmitterVisitor(body, declarations, source_path, module_name)

    declarations.emit_line(f"{native_function_header(fn.decl, emitter)} {{")
    body.indent()

    for r in all_values(fn.arg_regs, fn.blocks):
        if isinstance(r.type, RTuple):
            emitter.declare_tuple_struct(r.type)
        if isinstance(r.type, RArray):
            continue  # Special: declared on first assignment

        if r in fn.arg_regs:
            continue  # Skip the arguments

        ctype = emitter.ctype_spaced(r.type)
        init = ""
        declarations.emit_line(
            "{ctype}{prefix}{name}{init};".format(
                ctype=ctype, prefix=REG_PREFIX, name=names[r], init=init
            )
        )

    # Before we emit the blocks, give them all labels
    blocks = fn.blocks
    for i, block in enumerate(blocks):
        block.label = i

    # Find blocks that are never jumped to or are only jumped to from the
    # block directly above it. This allows for more labels and gotos to be
    # eliminated during code generation.
    for block in fn.blocks:
        terminator = block.terminator
        assert isinstance(terminator, ControlOp), terminator

        for target in terminator.targets():
            is_next_block = target.label == block.label + 1

            # Always emit labels for GetAttr error checks since the emit code that
            # generates them will add instructions between the branch and the
            # next label, causing the label to be wrongly removed. A better
            # solution would be to change the IR so that it adds a basic block
            # in between the calls.
            is_problematic_op = isinstance(terminator, Branch) and any(
                isinstance(s, GetAttr) for s in terminator.sources()
            )

            if not is_next_block or is_problematic_op:
                fn.blocks[target.label].referenced = True

    common = frequently_executed_blocks(fn.blocks[0])

    for i in range(len(blocks)):
        block = blocks[i]
        visitor.rare = block not in common
        next_block = None
        if i + 1 < len(blocks):
            next_block = blocks[i + 1]
        body.emit_label(block)
        visitor.next_block = next_block

        ops = block.ops
        visitor.ops = ops
        visitor.op_index = 0
        while visitor.op_index < len(ops):
            ops[visitor.op_index].accept(visitor)
            visitor.op_index += 1

    body.emit_line("}")

    emitter.emit_from_emitter(declarations)
    emitter.emit_from_emitter(body)


class FunctionEmitterVisitor(OpVisitor[None]):
    def __init__(
        self, emitter: Emitter, declarations: Emitter, source_path: str, module_name: str
    ) -> None:
        self.emitter = emitter
        self.names = emitter.names
        self.declarations = declarations
        self.source_path = source_path
        self.module_name = module_name
        self.literals = emitter.context.literals
        self.rare = False
        # Next basic block to be processed after the current one (if any), set by caller
        self.next_block: BasicBlock | None = None
        # Ops in the basic block currently being processed, set by caller
        self.ops: list[Op] = []
        # Current index within ops; visit methods can increment this to skip/merge ops
        self.op_index = 0

    def temp_name(self) -> str:
        return self.emitter.temp_name()

    def visit_goto(self, op: Goto) -> None:
        if op.label is not self.next_block:
            self.emit_line("goto %s;" % self.label(op.label))

    def error_value_check(self, value: Value, compare: str) -> str:
        typ = value.type
        if isinstance(typ, RTuple):
            # TODO: What about empty tuple?
            return self.emitter.tuple_undefined_check_cond(
                typ, self.reg(value), self.c_error_value, compare
            )
        elif isinstance(typ, RVec):
            # Error values for vecs are represented by a negative length.
            vec_compare = ">=" if compare == "!=" else "<"
            return f"{self.reg(value)}.len {vec_compare} 0"
        else:
            return f"{self.reg(value)} {compare} {self.c_error_value(typ)}"

    def visit_branch(self, op: Branch) -> None:
        true, false = op.true, op.false
        negated = op.negated
        negated_rare = False
        if true is self.next_block and op.traceback_entry is None:
            # Switch true/false since it avoids an else block.
            true, false = false, true
            negated = not negated
            negated_rare = True

        neg = "!" if negated else ""
        cond = ""
        if op.op == Branch.BOOL:
            expr_result = self.reg(op.value)
            cond = f"{neg}{expr_result}"
        elif op.op == Branch.IS_ERROR:
            compare = "!=" if negated else "=="
            cond = self.error_value_check(op.value, compare)
        else:
            assert False, "Invalid branch"

        # For error checks, tell the compiler the branch is unlikely
        if op.traceback_entry is not None or op.rare:
            if not negated_rare:
                cond = f"unlikely({cond})"
            else:
                cond = f"likely({cond})"

        if false is self.next_block:
            if op.traceback_entry is None:
                if true is not self.next_block:
                    self.emit_line(f"if ({cond}) goto {self.label(true)};")
            else:
                self.emit_line(f"if ({cond}) {{")
                self.emit_traceback(op)
                self.emit_lines("goto %s;" % self.label(true), "}")
        else:
            self.emit_line(f"if ({cond}) {{")
            self.emit_traceback(op)

            if true is not self.next_block:
                self.emit_line("goto %s;" % self.label(true))

            self.emit_lines("} else", "    goto %s;" % self.label(false))

    def visit_return(self, op: Return) -> None:
        value_str = self.reg(op.value)
        self.emit_line("return %s;" % value_str)

    def visit_tuple_set(self, op: TupleSet) -> None:
        dest = self.reg(op)
        tuple_type = op.tuple_type
        self.emitter.declare_tuple_struct(tuple_type)
        if len(op.items) == 0:  # empty tuple
            self.emit_line(f"{dest}.empty_struct_error_flag = 0;")
        else:
            for i, item in enumerate(op.items):
                self.emit_line(f"{dest}.f{i} = {self.reg(item)};")

    def visit_assign(self, op: Assign) -> None:
        dest = self.reg(op.dest)
        src = self.reg(op.src)
        # clang whines about self assignment (which we might generate
        # for some casts), so don't emit it.
        if dest != src:
            src_type = op.src.type
            dest_type = op.dest.type
            if src_type.is_unboxed and not dest_type.is_unboxed:
                # We sometimes assign from an integer prepresentation of a pointer
                # to a real pointer, and C compilers insist on a cast.
                src = f"(void *){src}"
            elif not src_type.is_unboxed and dest_type.is_unboxed:
                # We sometimes assign a pointer to an integer type (e.g. to create
                # tagged pointers), and here we need an explicit cast.
                src = f"({self.emitter.ctype(dest_type)}){src}"
            self.emit_line(f"{dest} = {src};")

    def visit_assign_multi(self, op: AssignMulti) -> None:
        typ = op.dest.type
        assert isinstance(typ, RArray), typ
        dest = self.reg(op.dest)
        # RArray values can only be assigned to once, so we can always
        # declare them on initialization.
        self.emit_line(
            "%s%s[%d] = %s;"
            % (
                self.emitter.ctype_spaced(typ.item_type),
                dest,
                len(op.src),
                c_array_initializer([self.reg(s) for s in op.src], indented=True),
            )
        )

    def visit_load_error_value(self, op: LoadErrorValue) -> None:
        reg = self.reg(op)
        if isinstance(op.type, RTuple):
            values = [self.c_undefined_value(item) for item in op.type.types]
            tmp = self.temp_name()
            self.emit_line("{} {} = {{ {} }};".format(self.ctype(op.type), tmp, ", ".join(values)))
            self.emit_line(f"{reg} = {tmp};")
        elif isinstance(op.type, RVec):
            self.emitter.set_undefined_value(reg, op.type)
        else:
            self.emit_line(f"{self.reg(op)} = {self.c_error_value(op.type)};")

    def visit_load_literal(self, op: LoadLiteral) -> None:
        index = self.literals.literal_index(op.value)
        if not is_int_rprimitive(op.type):
            self.emit_line("%s = CPyStatics[%d];" % (self.reg(op), index), ann=op.value)
        else:
            self.emit_line(
                "%s = (CPyTagged)CPyStatics[%d] | 1;" % (self.reg(op), index), ann=op.value
            )

    def get_attr_expr(self, obj: str, op: GetAttr | SetAttr, decl_cl: ClassIR) -> str:
        """Generate attribute accessor for normal (non-property) access.

        This either has a form like obj->attr_name for attributes defined in non-trait
        classes, and *(obj + attr_offset) for attributes defined by traits. We also
        insert all necessary C casts here.
        """
        # The struct cast below needs the defining group's __native.h
        # included by the consuming .c file. Record both the receiver
        # and declaring classes as cross-group deps.
        self.emitter.register_group_dep(op.class_type.class_ir)
        self.emitter.register_group_dep(decl_cl)
        cast = f"({op.class_type.struct_name(self.emitter.names)} *)"
        if decl_cl.is_trait and op.class_type.class_ir.is_trait:
            # For pure trait access find the offset first, offsets
            # are ordered by attribute position in the cl.attributes dict.
            # TODO: pre-calculate the mapping to make this faster.
            trait_attr_index = list(decl_cl.attributes).index(op.attr)
            # TODO: reuse these names somehow?
            offset = self.emitter.temp_name()
            self.declarations.emit_line(f"size_t {offset};")
            self.emitter.emit_line(
                "{} = {};".format(
                    offset,
                    "CPy_FindAttrOffset({}, {}, {})".format(
                        self.emitter.type_struct_name(decl_cl),
                        f"({cast}{obj})->vtable",
                        trait_attr_index,
                    ),
                )
            )
            attr_cast = f"({self.ctype(op.class_type.attr_type(op.attr))} *)"
            return f"*{attr_cast}((char *){obj} + {offset})"
        else:
            # Cast to something non-trait. Note: for this to work, all struct
            # members for non-trait classes must obey monotonic linear growth.
            if op.class_type.class_ir.is_trait:
                assert not decl_cl.is_trait
                cast = f"({decl_cl.struct_name(self.emitter.names)} *)"
            return f"({cast}{obj})->{self.emitter.attr(op.attr)}"

    def visit_get_attr(self, op: GetAttr) -> None:
        if op.allow_error_value:
            self.get_attr_with_allow_error_value(op)
            return
        dest = self.reg(op)
        obj = self.reg(op.obj)
        rtype = op.class_type
        cl = rtype.class_ir
        attr_rtype, decl_cl = cl.attr_details(op.attr)
        prefer_method = cl.is_trait and attr_rtype.error_overlap
        if cl.get_method(op.attr, prefer_method=prefer_method):
            # Properties are essentially methods, so use vtable access for them
            # (except for non-ext classes, which have no vtable — see emit_method_call).
            if not cl.is_ext_class or cl.is_method_final(op.attr):
                self.emit_method_call(f"{dest} = ", op.obj, op.attr, [])
            else:
                version = "_TRAIT" if cl.is_trait else ""
                self.emit_line(
                    "%s = CPY_GET_ATTR%s(%s, %s, %d, %s, %s); /* %s */"
                    % (
                        dest,
                        version,
                        obj,
                        self.emitter.type_struct_name(rtype.class_ir),
                        rtype.getter_index(op.attr),
                        rtype.struct_name(self.names),
                        self.ctype(rtype.attr_type(op.attr)),
                        op.attr,
                    )
                )
        else:
            # Otherwise, use direct or offset struct access.
            attr_expr = self.get_attr_expr(obj, op, decl_cl)
            self.emitter.emit_line(f"{dest} = {attr_expr};")
            always_defined = cl.is_always_defined(op.attr)
            merged_branch = None
            if not always_defined:
                self.emitter.emit_undefined_attr_check(
                    attr_rtype, dest, "==", obj, op.attr, cl, unlikely=True
                )
                branch = self.next_branch()
                if branch is not None:
                    if (
                        branch.value is op
                        and branch.op == Branch.IS_ERROR
                        and branch.traceback_entry is not None
                        and not branch.negated
                    ):
                        # Generate code for the following branch here to avoid
                        # redundant branches in the generated code.
                        self.emit_attribute_error(branch, cl.name, op.attr)
                        self.emit_line("goto %s;" % self.label(branch.true))
                        merged_branch = branch
                        self.emitter.emit_line("}")
                if not merged_branch:
                    exc_class = "PyExc_AttributeError"
                    self.emitter.emit_line(
                        'PyErr_SetString({}, "attribute {} of {} undefined");'.format(
                            exc_class,
                            repr(op.attr.removeprefix(GENERATOR_ATTRIBUTE_PREFIX)),
                            repr(cl.name),
                        )
                    )

            if attr_rtype.is_refcounted and not op.is_borrowed:
                if not merged_branch and not always_defined:
                    self.emitter.emit_line("} else {")
                self.emitter.emit_inc_ref(dest, attr_rtype)
            if merged_branch:
                if merged_branch.false is not self.next_block:
                    self.emit_line("goto %s;" % self.label(merged_branch.false))
                self.op_index += 1
            elif not always_defined:
                self.emitter.emit_line("}")

    def get_attr_with_allow_error_value(self, op: GetAttr) -> None:
        """Handle GetAttr with allow_error_value=True.

        This allows NULL or other error value without raising AttributeError.
        """
        dest = self.reg(op)
        obj = self.reg(op.obj)
        rtype = op.class_type
        cl = rtype.class_ir
        attr_rtype, decl_cl = cl.attr_details(op.attr)

        # Direct struct access without NULL check
        attr_expr = self.get_attr_expr(obj, op, decl_cl)
        self.emitter.emit_line(f"{dest} = {attr_expr};")

        # Only emit inc_ref if not NULL
        if attr_rtype.is_refcounted and not op.is_borrowed:
            check = self.error_value_check(op, "!=")
            self.emitter.emit_line(f"if ({check}) {{")
            self.emitter.emit_inc_ref(dest, attr_rtype)
            self.emitter.emit_line("}")

    def next_branch(self) -> Branch | None:
        if self.op_index + 1 < len(self.ops):
            next_op = self.ops[self.op_index + 1]
            if isinstance(next_op, Branch):
                return next_op
        return None

    def visit_set_attr(self, op: SetAttr) -> None:
        if op.error_kind == ERR_FALSE:
            dest = self.reg(op)
        obj = self.reg(op.obj)
        src = self.reg(op.src)
        rtype = op.class_type
        cl = rtype.class_ir
        attr_rtype, decl_cl = cl.attr_details(op.attr)
        if op.is_propset:
            # Again, use vtable access for properties...
            assert not op.is_init and op.error_kind == ERR_FALSE, "%s %d %d %s" % (
                op.attr,
                op.is_init,
                op.error_kind,
                rtype,
            )
            version = "_TRAIT" if cl.is_trait else ""
            self.emit_line(
                "%s = CPY_SET_ATTR%s(%s, %s, %d, %s, %s, %s); /* %s */"
                % (
                    dest,
                    version,
                    obj,
                    self.emitter.type_struct_name(rtype.class_ir),
                    rtype.setter_index(op.attr),
                    src,
                    rtype.struct_name(self.names),
                    self.ctype(rtype.attr_type(op.attr)),
                    op.attr,
                )
            )
        else:
            # ...and struct access for normal attributes.
            attr_expr = self.get_attr_expr(obj, op, decl_cl)
            if not op.is_init and attr_rtype.is_refcounted:
                # This is not an initialization (where we know that the attribute was
                # previously undefined), so decref the old value.
                always_defined = cl.is_always_defined(op.attr)
                if not always_defined:
                    self.emitter.emit_undefined_attr_check(
                        attr_rtype, attr_expr, "!=", obj, op.attr, cl
                    )
                self.emitter.emit_dec_ref(attr_expr, attr_rtype)
                if not always_defined:
                    self.emitter.emit_line("}")
            elif attr_rtype.error_overlap and not cl.is_always_defined(op.attr):
                # If there is overlap with the error value, update bitmap to mark
                # attribute as defined.
                self.emitter.emit_attr_bitmap_set(src, obj, attr_rtype, cl, op.attr)

            # This steals the reference to src, so we don't need to increment the arg
            self.emitter.emit_line(f"{attr_expr} = {src};")
            if op.error_kind == ERR_FALSE:
                self.emitter.emit_line(f"{dest} = 1;")

    def visit_load_static(self, op: LoadStatic) -> None:
        dest = self.reg(op)
        prefix = PREFIX_MAP[op.namespace]
        name = self.emitter.static_name(op.identifier, op.module_name, prefix)
        if op.namespace == NAMESPACE_TYPE:
            name = "(PyObject *)%s" % name
        self.emit_line(f"{dest} = {name};", ann=op.ann)

    def visit_init_static(self, op: InitStatic) -> None:
        value = self.reg(op.value)
        prefix = PREFIX_MAP[op.namespace]
        name = self.emitter.static_name(op.identifier, op.module_name, prefix)
        if op.namespace == NAMESPACE_TYPE:
            value = "(PyTypeObject *)%s" % value
        self.emit_line(f"{name} = {value};")
        self.emit_inc_ref(name, op.value.type)

    def visit_tuple_get(self, op: TupleGet) -> None:
        dest = self.reg(op)
        src = self.reg(op.src)
        self.emit_line(f"{dest} = {src}.f{op.index};")
        if not op.is_borrowed and op.type.is_refcounted:
            self.emit_inc_ref(dest, op.type)

    def get_dest_assign(self, dest: Value) -> str:
        if not dest.is_void:
            return self.reg(dest) + " = "
        else:
            return ""

    def visit_call(self, op: Call) -> None:
        """Call native function."""
        dest = self.get_dest_assign(op)
        args = ", ".join(self.reg(arg) for arg in op.args)
        lib = self.emitter.get_group_prefix(op.fn)
        cname = op.fn.cname(self.names)
        self.emit_line(f"{dest}{lib}{NATIVE_PREFIX}{cname}({args});")

    def visit_method_call(self, op: MethodCall) -> None:
        """Call native method."""
        dest = self.get_dest_assign(op)
        self.emit_method_call(dest, op.obj, op.method, op.args)

    def emit_method_call(self, dest: str, op_obj: Value, name: str, op_args: list[Value]) -> None:
        obj = self.reg(op_obj)
        rtype = op_obj.type
        assert isinstance(rtype, RInstance), rtype
        class_ir = rtype.class_ir
        # Use method_decl (not get_method) because under separate compilation the
        # FuncIR body may live in a different group — only its declaration is
        # visible here, and a decl is all we need to emit a direct C call
        # (the symbol resolves through that group's exports table).
        method_decl = rtype.class_ir.method_decl(name)

        # Can we call the method directly, bypassing vtable? Non-extension classes
        # don't have a vtable (compute_vtable is skipped for them), so the only
        # way to dispatch is a direct C call.
        is_direct = not class_ir.is_ext_class or class_ir.is_method_final(name)

        # The first argument gets omitted for static methods and
        # turned into the class for class methods
        obj_args = (
            []
            if method_decl.kind == FUNC_STATICMETHOD
            else [f"(PyObject *)Py_TYPE({obj})"] if method_decl.kind == FUNC_CLASSMETHOD else [obj]
        )
        args = ", ".join(obj_args + [self.reg(arg) for arg in op_args])
        mtype = native_function_type_from_decl(method_decl, self.emitter)
        version = "_TRAIT" if rtype.class_ir.is_trait else ""
        if is_direct:
            # Directly call method, without going through the vtable.
            self.emit_line(f"{dest}{self.emitter.native_function_call(method_decl)}({args});")
        else:
            # Call using vtable.
            method_idx = rtype.method_index(name)
            self.emit_line(
                "{}CPY_GET_METHOD{}({}, {}, {}, {}, {})({}); /* {} */".format(
                    dest,
                    version,
                    obj,
                    self.emitter.type_struct_name(rtype.class_ir),
                    method_idx,
                    rtype.struct_name(self.names),
                    mtype,
                    args,
                    name,
                )
            )

    def visit_inc_ref(self, op: IncRef) -> None:
        if (
            isinstance(op.src, Box)
            and (is_none_rprimitive(op.src.src.type) or is_bool_or_bit_rprimitive(op.src.src.type))
            and HAVE_IMMORTAL
        ):
            # On Python 3.12+, None/True/False are immortal, and we can skip inc ref
            return

        if isinstance(op.src, LoadLiteral) and HAVE_IMMORTAL:
            value = op.src.value
            # We can skip inc ref for immortal literals on Python 3.12+
            if type(value) is int and -5 <= value <= 256:
                # Small integers are immortal
                return

        src = self.reg(op.src)
        self.emit_inc_ref(src, op.src.type)

    def visit_dec_ref(self, op: DecRef) -> None:
        src = self.reg(op.src)
        self.emit_dec_ref(src, op.src.type, is_xdec=op.is_xdec)

    def visit_box(self, op: Box) -> None:
        self.emitter.emit_box(self.reg(op.src), self.reg(op), op.src.type, can_borrow=True)

    def visit_cast(self, op: Cast) -> None:
        if op.is_unchecked and op.is_borrowed:
            self.emit_line(f"{self.reg(op)} = {self.reg(op.src)};")
            return
        branch = self.next_branch()
        handler = None
        if branch is not None:
            if (
                branch.value is op
                and branch.op == Branch.IS_ERROR
                and branch.traceback_entry is not None
                and not branch.negated
                and branch.false is self.next_block
            ):
                # Generate code also for the following branch here to avoid
                # redundant branches in the generated code.
                handler = TracebackAndGotoHandler(
                    self.label(branch.true),
                    self.source_path,
                    self.module_name,
                    branch.traceback_entry,
                )
                self.op_index += 1

        self.emitter.emit_cast(
            self.reg(op.src), self.reg(op), op.type, src_type=op.src.type, error=handler
        )

    def visit_unbox(self, op: Unbox) -> None:
        self.emitter.emit_unbox(self.reg(op.src), self.reg(op), op.type)

    def visit_unreachable(self, op: Unreachable) -> None:
        self.emitter.emit_line("CPy_Unreachable();")

    def visit_raise_standard_error(self, op: RaiseStandardError) -> None:
        # TODO: Better escaping of backspaces and such
        if op.value is not None:
            if isinstance(op.value, str):
                message = op.value.replace('"', '\\"')
                self.emitter.emit_line(f'PyErr_SetString(PyExc_{op.class_name}, "{message}");')
            elif isinstance(op.value, Value):
                self.emitter.emit_line(
                    "PyErr_SetObject(PyExc_{}, {});".format(
                        op.class_name, self.emitter.reg(op.value)
                    )
                )
            else:
                assert False, "op value type must be either str or Value"
        else:
            self.emitter.emit_line(f"PyErr_SetNone(PyExc_{op.class_name});")
        self.emitter.emit_line(f"{self.reg(op)} = 0;")

    def visit_call_c(self, op: CallC) -> None:
        if op.is_void:
            dest = ""
        else:
            dest = self.get_dest_assign(op)
        args = ", ".join(self.reg(arg) for arg in op.args)
        self.emitter.emit_line(f"{dest}{op.function_name}({args});")

    def visit_primitive_op(self, op: PrimitiveOp) -> None:
        raise RuntimeError(
            f"unexpected PrimitiveOp {op.desc.name}: they must be lowered before codegen"
        )

    def visit_truncate(self, op: Truncate) -> None:
        dest = self.reg(op)
        value = self.reg(op.src)
        # for C backend the generated code are straight assignments
        self.emit_line(f"{dest} = {value};")

    def visit_extend(self, op: Extend) -> None:
        dest = self.reg(op)
        value = self.reg(op.src)
        if op.signed:
            src_cast = self.emit_signed_int_cast(op.src.type)
        else:
            src_cast = self.emit_unsigned_int_cast(op.src.type)
        self.emit_line(f"{dest} = {src_cast}{value};")

    def visit_load_global(self, op: LoadGlobal) -> None:
        dest = self.reg(op)
        self.emit_line(f"{dest} = {op.identifier};", ann=op.ann)

    def visit_int_op(self, op: IntOp) -> None:
        dest = self.reg(op)
        lhs = self.reg(op.lhs)
        rhs = self.reg(op.rhs)
        if op.op == IntOp.RIGHT_SHIFT:
            # Signed right shift
            lhs = self.emit_signed_int_cast(op.lhs.type) + lhs
            rhs = self.emit_signed_int_cast(op.rhs.type) + rhs
        self.emit_line(f"{dest} = {lhs} {op.op_str[op.op]} {rhs};")

    def visit_comparison_op(self, op: ComparisonOp) -> None:
        dest = self.reg(op)
        lhs = self.reg(op.lhs)
        rhs = self.reg(op.rhs)
        lhs_cast = ""
        rhs_cast = ""
        if op.op in (ComparisonOp.SLT, ComparisonOp.SGT, ComparisonOp.SLE, ComparisonOp.SGE):
            # Always signed comparison op
            lhs_cast = self.emit_signed_int_cast(op.lhs.type)
            rhs_cast = self.emit_signed_int_cast(op.rhs.type)
        elif op.op in (ComparisonOp.ULT, ComparisonOp.UGT, ComparisonOp.ULE, ComparisonOp.UGE):
            # Always unsigned comparison op
            lhs_cast = self.emit_unsigned_int_cast(op.lhs.type)
            rhs_cast = self.emit_unsigned_int_cast(op.rhs.type)
        elif isinstance(op.lhs, Integer) and op.lhs.value < 0:
            # Force signed ==/!= with negative operand
            rhs_cast = self.emit_signed_int_cast(op.rhs.type)
        elif isinstance(op.rhs, Integer) and op.rhs.value < 0:
            # Force signed ==/!= with negative operand
            lhs_cast = self.emit_signed_int_cast(op.lhs.type)
        self.emit_line(f"{dest} = {lhs_cast}{lhs} {op.op_str[op.op]} {rhs_cast}{rhs};")

    def visit_float_op(self, op: FloatOp) -> None:
        dest = self.reg(op)
        lhs = self.reg(op.lhs)
        rhs = self.reg(op.rhs)
        if op.op != FloatOp.MOD:
            self.emit_line(f"{dest} = {lhs} {op.op_str[op.op]} {rhs};")
        else:
            # TODO: This may set errno as a side effect, that is a little sketchy.
            self.emit_line(f"{dest} = fmod({lhs}, {rhs});")

    def visit_float_neg(self, op: FloatNeg) -> None:
        dest = self.reg(op)
        src = self.reg(op.src)
        self.emit_line(f"{dest} = -{src};")

    def visit_float_comparison_op(self, op: FloatComparisonOp) -> None:
        dest = self.reg(op)
        lhs = self.reg(op.lhs)
        rhs = self.reg(op.rhs)
        self.emit_line(f"{dest} = {lhs} {op.op_str[op.op]} {rhs};")

    def visit_load_mem(self, op: LoadMem) -> None:
        dest = self.reg(op)
        src = self.reg(op.src)
        # TODO: we shouldn't dereference to type that are pointer type so far
        type = self.ctype(op.type)
        self.emit_line(f"{dest} = *({type} *){src};")
        if not op.is_borrowed and op.type.is_refcounted:
            self.emit_inc_ref(dest, op.type)

    def visit_set_mem(self, op: SetMem) -> None:
        dest = self.reg(op.dest)
        src = self.reg(op.src)
        dest_type = self.ctype(op.dest_type)
        # clang whines about self assignment (which we might generate
        # for some casts), so don't emit it.
        if dest != src:
            self.emit_line(f"*({dest_type} *){dest} = {src};")

    def visit_get_element(self, op: GetElement) -> None:
        dest = self.reg(op)
        src = self.reg(op.src)
        dest_type = self.ctype(op.type)
        self.emit_line(f"{dest} = ({dest_type}){src}.{op.field};")

    def visit_get_element_ptr(self, op: GetElementPtr) -> None:
        dest = self.reg(op)
        src = self.reg(op.src)
        # TODO: support tuple type
        assert isinstance(op.src_type, RStruct), op.src_type
        assert op.field in op.src_type.names, "Invalid field name."
        # Use offsetof to avoid undefined behavior when src is NULL
        # (e.g., vec buf pointer for empty vecs). The &((T*)p)->field
        # pattern is UB when p is NULL, which GCC -O3 can exploit.
        self.emit_line(
            "{} = ({})((CPyPtr){} + offsetof({}, {}));".format(
                dest, op.type._ctype, src, op.src_type.name, op.field
            )
        )

    def visit_set_element(self, op: SetElement) -> None:
        dest = self.reg(op)
        field = op.field
        item = self.set_element_item(op.src.type, field, self.reg(op.item))
        if isinstance(op.src, Undef):
            # First assignment to an undefined struct is trivial.
            self.emit_line(f"{dest}.{field} = {item};")
        else:
            # In the general case create a copy of the struct with a single
            # item modified.
            #
            # TODO: Can we do better if only a subset of fields are initialized?
            # TODO: Make this less verbose in the common case
            # TODO: Support tuples (or use RStruct for tuples)?
            src = self.reg(op.src)
            src_type = op.src.type
            assert isinstance(src_type, (RStruct, RVec)), src_type
            init_items = []
            for n in src_type.names:
                if n != field:
                    init_items.append(f"{src}.{n}")
                else:
                    init_items.append(item)
            self.emit_line(f"{dest} = ({self.ctype(src_type)}) {{ {', '.join(init_items)} }};")

    def set_element_item(self, src_type: RType, field: str, item: str) -> str:
        if field == "items" and src_type._ctype in VEC_ITEMS_C_TYPE:
            return f"({VEC_ITEMS_C_TYPE[src_type._ctype]}){item}"
        return item

    def visit_load_address(self, op: LoadAddress) -> None:
        typ = op.type
        dest = self.reg(op)
        if isinstance(op.src, Register):
            src = self.reg(op.src)
        elif isinstance(op.src, LoadStatic):
            prefix = PREFIX_MAP[op.src.namespace]
            src = self.emitter.static_name(op.src.identifier, op.src.module_name, prefix)
        else:
            src = op.src
        self.emit_line(f"{dest} = ({typ._ctype})&{src};")

    def visit_keep_alive(self, op: KeepAlive) -> None:
        # This is a no-op.
        pass

    def visit_unborrow(self, op: Unborrow) -> None:
        # This is a no-op that propagates the source value.
        dest = self.reg(op)
        src = self.reg(op.src)
        self.emit_line(f"{dest} = {src};")

    # Helpers

    def label(self, label: BasicBlock) -> str:
        return self.emitter.label(label)

    def reg(self, reg: Value) -> str:
        if isinstance(reg, Integer):
            val = reg.value
            if val == 0 and is_pointer_rprimitive(reg.type):
                return "NULL"
            s = str(val)
            if val >= (1 << 31):
                # Avoid overflowing signed 32-bit int
                if val >= (1 << 63):
                    s += "ULL"
                else:
                    s += "LL"
            elif val == -(1 << 63):
                # Avoid overflowing C integer literal
                s = "(-9223372036854775807LL - 1)"
            elif val <= -(1 << 31):
                s += "LL"
            return s
        elif isinstance(reg, Float):
            r = repr(reg.value)
            if r == "inf":
                return "INFINITY"
            elif r == "-inf":
                return "-INFINITY"
            elif r == "nan":
                return "NAN"
            return r
        elif isinstance(reg, CString):
            return '"' + encode_c_string_literal(reg.value) + '"'
        else:
            return self.emitter.reg(reg)

    def ctype(self, rtype: RType) -> str:
        return self.emitter.ctype(rtype)

    def c_error_value(self, rtype: RType) -> str:
        return self.emitter.c_error_value(rtype)

    def c_undefined_value(self, rtype: RType) -> str:
        return self.emitter.c_undefined_value(rtype)

    def emit_line(self, line: str, *, ann: object = None) -> None:
        self.emitter.emit_line(line, ann=ann)

    def emit_lines(self, *lines: str) -> None:
        self.emitter.emit_lines(*lines)

    def emit_inc_ref(self, dest: str, rtype: RType) -> None:
        self.emitter.emit_inc_ref(dest, rtype, rare=self.rare)

    def emit_dec_ref(self, dest: str, rtype: RType, is_xdec: bool) -> None:
        self.emitter.emit_dec_ref(dest, rtype, is_xdec=is_xdec, rare=self.rare)

    def emit_declaration(self, line: str) -> None:
        self.declarations.emit_line(line)

    def emit_traceback(self, op: Branch) -> None:
        if op.traceback_entry is not None:
            self.emitter.emit_traceback(self.source_path, self.module_name, op.traceback_entry)

    def emit_attribute_error(self, op: Branch, class_name: str, attr: str) -> None:
        assert op.traceback_entry is not None
        if self.emitter.context.strict_traceback_checks:
            assert (
                op.traceback_entry[1] >= 0
            ), "AttributeError traceback cannot have a negative line number"
        globals_static = self.emitter.static_name("globals", self.module_name)
        self.emit_line(
            'CPy_AttributeError("%s", "%s", "%s", "%s", %d, %s);'
            % (
                self.source_path.replace("\\", "\\\\"),
                op.traceback_entry[0],
                class_name,
                attr.removeprefix(GENERATOR_ATTRIBUTE_PREFIX),
                op.traceback_entry[1],
                globals_static,
            )
        )
        if DEBUG_ERRORS:
            self.emit_line('assert(PyErr_Occurred() != NULL && "failure w/o err!");')

    def emit_signed_int_cast(self, type: RType) -> str:
        if is_tagged(type):
            return "(Py_ssize_t)"
        else:
            return ""

    def emit_unsigned_int_cast(self, type: RType) -> str:
        if is_int32_rprimitive(type):
            return "(uint32_t)"
        elif is_int64_rprimitive(type):
            return "(uint64_t)"
        else:
            return ""


_translation_table: Final[dict[int, str]] = {}


def encode_c_string_literal(b: bytes) -> str:
    """Convert bytestring to the C string literal syntax (with necessary escaping).

    For example, b'foo\n' gets converted to 'foo\\n' (note that double quotes are not added).
    """
    if not _translation_table:
        # Initialize the translation table on the first call.
        d = {
            ord("\n"): "\\n",
            ord("\r"): "\\r",
            ord("\t"): "\\t",
            ord('"'): '\\"',
            ord("\\"): "\\\\",
        }
        for i in range(256):
            if i not in d:
                if i < 32 or i >= 127:
                    d[i] = "\\x%.2x" % i
                else:
                    d[i] = chr(i)
        _translation_table.update(str.maketrans(d))
    return b.decode("latin1").translate(_translation_table)
