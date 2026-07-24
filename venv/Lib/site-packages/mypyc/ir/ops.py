"""Low-level opcodes for compiler intermediate representation (IR).

Opcodes operate on abstract values (Value) in a register machine. Each
value has a type (RType). A value can hold various things, such as:

- local variables or temporaries (Register)
- intermediate values of expressions (RegisterOp subclasses)
- condition flags (true/false)
- literals (integer literals, True, False, etc.)

NOTE: As a convention, we don't create subclasses of concrete Value/Op
      subclasses (e.g. you shouldn't define a subclass of Integer, which
      is a concrete class).

      If you want to introduce a variant of an existing class, you'd
      typically add an attribute (e.g. a flag) to an existing concrete
      class to enable the new behavior. Sometimes adding a new abstract
      base class is also an option, or just creating a new subclass
      without any inheritance relationship (some duplication of code
      is preferred over introducing complex implementation inheritance).

      This makes it possible to use isinstance(x, <concrete Value
      subclass>) checks without worrying about potential subclasses.
"""

from __future__ import annotations

from abc import abstractmethod
from collections.abc import Sequence
from typing import TYPE_CHECKING, Final, Generic, NamedTuple, TypeVar, final

from mypy_extensions import trait

from mypyc.common import PROPSET_PREFIX
from mypyc.ir.deps import Dependency
from mypyc.ir.rtypes import (
    RArray,
    RInstance,
    RStruct,
    RTuple,
    RType,
    RUnion,
    RVec,
    RVoid,
    bit_rprimitive,
    bool_rprimitive,
    cstring_rprimitive,
    float_rprimitive,
    int_rprimitive,
    is_bool_or_bit_rprimitive,
    is_fixed_width_rtype,
    is_int_rprimitive,
    is_none_rprimitive,
    is_pointer_rprimitive,
    is_short_int_rprimitive,
    object_rprimitive,
    pointer_rprimitive,
    short_int_rprimitive,
    void_rtype,
)

if TYPE_CHECKING:
    from mypyc.codegen.literals import LiteralValue
    from mypyc.ir.class_ir import ClassIR
    from mypyc.ir.func_ir import FuncDecl, FuncIR

T = TypeVar("T")


@final
class BasicBlock:
    """IR basic block.

    Contains a sequence of Ops and ends with a ControlOp (Goto,
    Branch, Return or Unreachable). Only the last op can be a
    ControlOp.

    All generated Ops live in basic blocks. Basic blocks determine the
    order of evaluation and control flow within a function. A basic
    block is always associated with a single function/method (FuncIR).

    When building the IR, ops that raise exceptions can be included in
    the middle of a basic block, but the exceptions aren't checked.
    Afterwards we perform a transform that inserts explicit checks for
    all error conditions and splits basic blocks accordingly to preserve
    the invariant that a jump, branch or return can only ever appear
    as the final op in a block. Manually inserting error checking ops
    would be boring and error-prone.

    BasicBlocks have an error_handler attribute that determines where
    to jump if an error occurs. If none is specified, an error will
    propagate up out of the function. This is compiled away by the
    `exceptions` module.

    Block labels are used for pretty printing and emitting C code, and
    get filled in by those passes.

    Ops that may terminate the program aren't treated as exits.
    """

    def __init__(self, label: int = -1) -> None:
        self.label = label
        self.ops: list[Op] = []
        self.error_handler: BasicBlock | None = None
        self.referenced = False

    @property
    def terminated(self) -> bool:
        """Does the block end with a jump, branch or return?

        This should always be true after the basic block has been fully built, but
        this is false during construction.
        """
        return bool(self.ops) and isinstance(self.ops[-1], ControlOp)

    @property
    def terminator(self) -> ControlOp:
        """The terminator operation of the block."""
        assert bool(self.ops) and isinstance(self.ops[-1], ControlOp)
        return self.ops[-1]


# Never generates an exception
ERR_NEVER: Final = 0
# Generates magic value (c_error_value) based on target RType on exception
ERR_MAGIC: Final = 1
# Generates false (bool) on exception
ERR_FALSE: Final = 2
# Always fails
ERR_ALWAYS: Final = 3
# Like ERR_MAGIC, but the magic return overlaps with a possible return value, and
# an extra PyErr_Occurred() check is also required
ERR_MAGIC_OVERLAPPING: Final = 4

# Hack: using this line number for an op will suppress it in tracebacks
NO_TRACEBACK_LINE_NO = -10000


class Value:
    """Abstract base class for all IR values.

    These include references to registers, literals, and all
    operations (Ops), such as assignments, calls and branches.

    Values are often used as inputs of Ops. Register can be used as an
    assignment target.

    A Value is part of the IR being compiled if it's included in a BasicBlock
    that is reachable from a FuncIR (i.e., is part of a function).

    See also: Op is a subclass of Value that is the base class of all
    operations.
    """

    # Source line number (-1 for no/unknown line)
    line = -1
    # Type of the value or the result of the operation
    type: RType = void_rtype
    is_borrowed = False

    @property
    def is_void(self) -> bool:
        return isinstance(self.type, RVoid)


@final
class Register(Value):
    """A Register holds a value of a specific type, and it can be read and mutated.

    A Register is always local to a function. Each local variable maps
    to a Register, and they are also used for some (but not all)
    temporary values.

    Note that the term 'register' is overloaded and is sometimes used
    to refer to arbitrary Values (for example, in RegisterOp).
    """

    def __init__(self, type: RType, name: str = "", is_arg: bool = False, line: int = -1) -> None:
        self.type = type
        self.name = name
        self.is_arg = is_arg
        self.is_borrowed = is_arg
        self.line = line

    @property
    def is_void(self) -> bool:
        return False

    def __repr__(self) -> str:
        return f"<Register {self.name!r} at {hex(id(self))}>"


@final
class Integer(Value):
    """Short integer literal.

    Integer literals are treated as constant values and are generally
    not included in data flow analyses and such, unlike Register and
    Op subclasses.

    Integer can represent multiple types:

     * Short tagged integers (short_int_primitive type; the tag bit is clear)
     * Ordinary fixed-width integers (e.g., int32_rprimitive)
     * Values of other unboxed primitive types that are represented as integers
       (none_rprimitive, bool_rprimitive)
     * Null pointers (value 0) of various types, including object_rprimitive
    """

    def __init__(self, value: int, rtype: RType = short_int_rprimitive, line: int = -1) -> None:
        if is_short_int_rprimitive(rtype) or is_int_rprimitive(rtype):
            self.value = value * 2
        else:
            self.value = value
        self.type = rtype
        self.line = line

    def numeric_value(self) -> int:
        if is_short_int_rprimitive(self.type) or is_int_rprimitive(self.type):
            return self.value // 2
        return self.value


@final
class Float(Value):
    """Float literal.

    Floating point literals are treated as constant values and are generally
    not included in data flow analyses and such, unlike Register and
    Op subclasses.
    """

    def __init__(self, value: float, line: int = -1) -> None:
        self.value = value
        self.type = float_rprimitive
        self.line = line


@final
class CString(Value):
    """C string literal (zero-terminated).

    You can also include zero values in the value, but then you'll need to track
    the length of the string separately.
    """

    def __init__(self, value: bytes, line: int = -1) -> None:
        self.value = value
        self.type = cstring_rprimitive
        self.line = line


@final
class Undef(Value):
    """An undefined value.

    Use Undef() as the initial value followed by one or more SetElement
    ops to initialize a struct. Pseudocode example:

      r0 = set_element undef MyStruct, "field1", f1
      r1 = set_element r0, "field2", f2
      # r1 now has new struct value with two fields set

    Warning: Always initialize undefined values before using them,
    as otherwise the values are garbage. You shouldn't expect that
    undefined values are zeroed, in particular.
    """

    def __init__(self, rtype: RType) -> None:
        self.type = rtype


class Op(Value):
    """Abstract base class for all IR operations.

    Each operation must be stored in a BasicBlock (in 'ops') to be
    active in the IR. This is different from non-Op values, including
    Register and Integer, where a reference from an active Op is
    sufficient to be considered active.

    In well-formed IR an active Op has no references to inactive ops
    or ops used in another function.
    """

    def __init__(self, line: int) -> None:
        self.line = line

    def can_raise(self) -> bool:
        # Override this is if Op may raise an exception. Note that currently the fact that
        # only RegisterOps may raise an exception in hard coded in some places.
        return False

    @abstractmethod
    def sources(self) -> list[Value]:
        """All the values the op may read."""

    @abstractmethod
    def set_sources(self, new: list[Value]) -> None:
        """Rewrite the sources of an op"""

    def stolen(self) -> list[Value]:
        """Return arguments that have a reference count stolen by this op"""
        return []

    def unique_sources(self) -> list[Value]:
        result: list[Value] = []
        for reg in self.sources():
            if reg not in result:
                result.append(reg)
        return result

    @abstractmethod
    def accept(self, visitor: OpVisitor[T]) -> T:
        pass


class BaseAssign(Op):
    """Abstract base class for ops that assign to a register."""

    def __init__(self, dest: Register, line: int = -1) -> None:
        super().__init__(line)
        self.dest = dest


@final
class Assign(BaseAssign):
    """Assign a value to a Register (dest = src)."""

    error_kind = ERR_NEVER

    def __init__(self, dest: Register, src: Value, line: int = -1) -> None:
        super().__init__(dest, line)
        self.src = src

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        if not self.dest.type.is_refcounted:
            return []
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_assign(self)


@final
class AssignMulti(BaseAssign):
    """Assign multiple values to a Register (dest = src1, src2, ...).

    This is used to initialize RArray values. It's provided to avoid
    very verbose IR for common vectorcall operations.

    Note that this interacts atypically with reference counting. We
    assume that each RArray register is initialized exactly once
    with this op.
    """

    error_kind = ERR_NEVER

    def __init__(self, dest: Register, src: list[Value], line: int = -1) -> None:
        super().__init__(dest, line)
        assert src
        assert isinstance(dest.type, RArray)
        assert dest.type.length == len(src)
        self.src = src

    def sources(self) -> list[Value]:
        return self.src.copy()

    def set_sources(self, new: list[Value]) -> None:
        self.src = new[:]

    def stolen(self) -> list[Value]:
        return []

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_assign_multi(self)


class ControlOp(Op):
    """Abstract base class for control flow operations."""

    def targets(self) -> Sequence[BasicBlock]:
        """Get all basic block targets of the control operation."""
        return ()

    def set_target(self, i: int, new: BasicBlock) -> None:
        """Update a basic block target."""
        raise AssertionError(f"Invalid set_target({self}, {i})")


@final
class Goto(ControlOp):
    """Unconditional jump."""

    error_kind = ERR_NEVER

    def __init__(self, label: BasicBlock, line: int = -1) -> None:
        super().__init__(line)
        self.label = label

    def targets(self) -> Sequence[BasicBlock]:
        return (self.label,)

    def set_target(self, i: int, new: BasicBlock) -> None:
        assert i == 0
        self.label = new

    def __repr__(self) -> str:
        return "<Goto %s>" % self.label.label

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_goto(self)


@final
class Branch(ControlOp):
    """Branch based on a value.

    If op is BOOL, branch based on a bit/bool value:
       if [not] r1 goto L1 else goto L2

    If op is IS_ERROR, branch based on whether there is an error value:
       if [not] is_error(r1) goto L1 else goto L2
    """

    # Branch ops never raise an exception.
    error_kind = ERR_NEVER

    BOOL: Final = 100
    IS_ERROR: Final = 101

    def __init__(
        self,
        value: Value,
        true_label: BasicBlock,
        false_label: BasicBlock,
        op: int,
        line: int = -1,
        *,
        rare: bool = False,
    ) -> None:
        super().__init__(line)
        # Target value being checked
        self.value = value
        # Branch here if the condition is true
        self.true = true_label
        # Branch here if the condition is false
        self.false = false_label
        # Branch.BOOL (boolean check) or Branch.IS_ERROR (error value check)
        self.op = op
        # If True, the condition is negated
        self.negated = False
        # If not None, the true label should generate a traceback entry (func name, line number)
        self.traceback_entry: tuple[str, int] | None = None
        # If True, we expect to usually take the false branch (for optimization purposes);
        # this is implicitly treated as true if there is a traceback entry
        self.rare = rare

    def targets(self) -> Sequence[BasicBlock]:
        return (self.true, self.false)

    def set_target(self, i: int, new: BasicBlock) -> None:
        assert i == 0 or i == 1
        if i == 0:
            self.true = new
        else:
            self.false = new

    def sources(self) -> list[Value]:
        return [self.value]

    def set_sources(self, new: list[Value]) -> None:
        (self.value,) = new

    def invert(self) -> None:
        self.negated = not self.negated

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_branch(self)


@final
class Return(ControlOp):
    """Return a value from a function."""

    error_kind = ERR_NEVER

    def __init__(
        self, value: Value, line: int = -1, *, yield_target: BasicBlock | None = None
    ) -> None:
        super().__init__(line)
        self.value = value
        # If this return is created by a yield, keep track of the next
        # basic block. This doesn't affect the code we generate but
        # can feed into analysis that need to understand the
        # *original* CFG.
        self.yield_target = yield_target

    def sources(self) -> list[Value]:
        return [self.value]

    def set_sources(self, new: list[Value]) -> None:
        (self.value,) = new

    def stolen(self) -> list[Value]:
        return [self.value]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_return(self)


@final
class Unreachable(ControlOp):
    """Mark the end of basic block as unreachable.

    This is sometimes necessary when the end of a basic block is never
    reached. This can also be explicitly added to the end of non-None
    returning functions (in None-returning function we can just return
    None).

    Mypy statically guarantees that the end of the function is not
    unreachable if there is not a return statement.

    This prevents the block formatter from being confused due to lack
    of a leave and also leaves a nifty note in the IR. It is not
    generally processed by visitors.
    """

    error_kind = ERR_NEVER

    def __init__(self, line: int = -1) -> None:
        super().__init__(line)

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_unreachable(self)


class RegisterOp(Op):
    """Abstract base class for operations that can be written as r1 = f(r2, ..., rn).

    Takes some values, performs an operation, and generates an output
    (unless the 'type' attribute is void_rtype, which is the default).
    Other ops can refer to the result of the Op by referring to the Op
    instance. This doesn't do any explicit control flow, but can raise an
    error.

    Note that the operands can be arbitrary Values, not just Register
    instances, even though the naming may suggest otherwise.
    """

    error_kind = -1  # Can this raise exception and how is it signalled; one of ERR_*

    _type: RType | None = None

    def __init__(self, line: int) -> None:
        super().__init__(line)
        assert self.error_kind != -1, "error_kind not defined"

    def can_raise(self) -> bool:
        return self.error_kind != ERR_NEVER


@final
class IncRef(RegisterOp):
    """Increase reference count (inc_ref src)."""

    error_kind = ERR_NEVER

    def __init__(self, src: Value, line: int = -1) -> None:
        assert src.type.is_refcounted
        super().__init__(line)
        self.src = src

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_inc_ref(self)


@final
class DecRef(RegisterOp):
    """Decrease reference count and free object if zero (dec_ref src).

    The is_xdec flag says to use an XDECREF, which checks if the
    pointer is NULL first.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, is_xdec: bool = False, line: int = -1) -> None:
        assert src.type.is_refcounted
        super().__init__(line)
        self.src = src
        self.is_xdec = is_xdec

    def __repr__(self) -> str:
        return "<{}DecRef {!r}>".format("X" if self.is_xdec else "", self.src)

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_dec_ref(self)


@final
class Call(RegisterOp):
    """Native call f(arg, ...).

    The call target can be a module-level function or a class.
    """

    def __init__(self, fn: FuncDecl, args: Sequence[Value], line: int) -> None:
        self.fn = fn
        self.args = list(args)
        assert len(self.args) == len(fn.sig.args)
        self.type = fn.sig.ret_type
        ret_type = fn.sig.ret_type
        if not ret_type.error_overlap:
            self.error_kind = ERR_MAGIC
        else:
            self.error_kind = ERR_MAGIC_OVERLAPPING
        super().__init__(line)

    def sources(self) -> list[Value]:
        return list(self.args.copy())

    def set_sources(self, new: list[Value]) -> None:
        self.args = new[:]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_call(self)


@final
class MethodCall(RegisterOp):
    """Native method call obj.method(arg, ...)"""

    def __init__(self, obj: Value, method: str, args: list[Value], line: int = -1) -> None:
        self.obj = obj
        self.method = method
        self.args = args
        assert isinstance(obj.type, RInstance), "Methods can only be called on instances"
        self.receiver_type = obj.type
        method_ir = self.receiver_type.class_ir.method_sig(method)
        assert method_ir is not None, "{} doesn't have method {}".format(
            self.receiver_type.name, method
        )
        ret_type = method_ir.ret_type
        self.type = ret_type
        if not ret_type.error_overlap:
            self.error_kind = ERR_MAGIC
        else:
            self.error_kind = ERR_MAGIC_OVERLAPPING
        super().__init__(line)

    def sources(self) -> list[Value]:
        return self.args.copy() + [self.obj]

    def set_sources(self, new: list[Value]) -> None:
        *self.args, self.obj = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_method_call(self)


@final
class PrimitiveDescription:
    """Description of a primitive op.

    Primitives get lowered into lower-level ops before code generation.

    If c_function_name is provided, a primitive will be lowered into a CallC op.
    Otherwise, custom logic will need to be implemented to transform the
    primitive into lower-level ops.
    """

    def __init__(
        self,
        name: str,
        arg_types: list[RType],
        return_type: RType,  # TODO: What about generic?
        var_arg_type: RType | None,
        truncated_type: RType | None,
        c_function_name: str | None,
        error_kind: int,
        steals: StealsDescription,
        is_borrowed: bool,
        ordering: list[int] | None,
        extra_int_constants: list[tuple[int, RType]],
        priority: int,
        is_pure: bool,
        experimental: bool,
        dependencies: list[Dependency] | None,
    ) -> None:
        # Each primitive much have a distinct name, but otherwise they are arbitrary.
        self.name: Final = name
        self.arg_types: Final = arg_types
        self.return_type: Final = return_type
        self.var_arg_type: Final = var_arg_type
        self.truncated_type: Final = truncated_type
        # If non-None, this will map to a call of a C helper function; if None,
        # there must be a custom handler function that gets invoked during the lowering
        # pass to generate low-level IR for the primitive (in the mypyc.lower package)
        self.c_function_name: Final = c_function_name
        self.error_kind: Final = error_kind
        self.steals: Final = steals
        self.is_borrowed: Final = is_borrowed
        self.ordering: Final = ordering
        self.extra_int_constants: Final = extra_int_constants
        self.priority: Final = priority
        # Pure primitives have no side effects, take immutable arguments, and
        # never fail. They support additional optimizations.
        self.is_pure: Final = is_pure
        if is_pure:
            assert error_kind == ERR_NEVER
        # Experimental primitives are not used unless mypyc experimental features are
        # explicitly enabled
        self.experimental = experimental
        # Dependencies for the primitive, such as a capsule that needs to imported
        # and configured to call the primitive.
        self.dependencies = dependencies
        # Native integer types such as u8 can cause ambiguity in primitive
        # matching, since these are assignable to plain int *and* vice versa.
        # If this flag is set, the primitive has native integer types and must
        # be matched using more complex rules.
        self.is_ambiguous = any(has_fixed_width_int(t) for t in arg_types)

    def __repr__(self) -> str:
        return f"<PrimitiveDescription {self.name!r}: {self.arg_types}>"


def has_fixed_width_int(t: RType) -> bool:
    if isinstance(t, RTuple):
        return any(has_fixed_width_int(t) for t in t.types)
    elif isinstance(t, RUnion):
        return any(has_fixed_width_int(t) for t in t.items)
    return is_fixed_width_rtype(t)


@final
class PrimitiveOp(RegisterOp):
    """A higher-level primitive operation.

    Some of these have special compiler support. These will be lowered
    (transformed) into lower-level IR ops before code generation, and after
    reference counting op insertion. Others will be transformed into CallC
    ops.

    Tagged integer equality is a typical primitive op with non-trivial
    lowering. It gets transformed into a tag check, followed by different
    code paths for short and long representations.
    """

    def __init__(self, args: list[Value], desc: PrimitiveDescription, line: int = -1) -> None:
        self.error_kind = desc.error_kind
        super().__init__(line)
        self.args = args
        self.type = desc.return_type
        self.desc = desc

    def sources(self) -> list[Value]:
        return self.args

    def set_sources(self, new: list[Value]) -> None:
        self.args = new[:]

    def stolen(self) -> list[Value]:
        steals = self.desc.steals
        if isinstance(steals, list):
            assert len(steals) == len(self.args)
            return [arg for arg, steal in zip(self.args, steals) if steal]
        else:
            return [] if not steals else self.sources()

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_primitive_op(self)


@final
class LoadErrorValue(RegisterOp):
    """Load an error value.

    Each type has one reserved value that signals an error (exception). This
    loads the error value for a specific type.
    """

    error_kind = ERR_NEVER

    def __init__(
        self, rtype: RType, line: int = -1, is_borrowed: bool = False, undefines: bool = False
    ) -> None:
        super().__init__(line)
        self.type = rtype
        self.is_borrowed = is_borrowed
        # Undefines is true if this should viewed by the definedness
        # analysis pass as making the register it is assigned to
        # undefined (and thus checks should be added on uses).
        self.undefines = undefines

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_error_value(self)


@final
class LoadLiteral(RegisterOp):
    """Load a Python literal object (dest = 'foo' / b'foo' / ...).

    This is used to load a static PyObject * value corresponding to
    a literal of one of the supported types.

    Tuple / frozenset literals must contain only valid literal values as items.

    NOTE: You can use this to load boxed (Python) int objects. Use
          Integer to load unboxed, tagged integers or fixed-width,
          low-level integers.

          For int literals, both int_rprimitive (CPyTagged) and
          object_primitive (PyObject *) are supported as rtype. However,
          when using int_rprimitive, the value must *not* be small enough
          to fit in an unboxed integer.
    """

    error_kind = ERR_NEVER
    is_borrowed = True

    def __init__(self, value: LiteralValue, rtype: RType, line: int = -1) -> None:
        super().__init__(line)
        self.value = value
        self.type = rtype

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_literal(self)


@final
class GetAttr(RegisterOp):
    """obj.attr (for a native object)"""

    error_kind = ERR_MAGIC

    def __init__(
        self,
        obj: Value,
        attr: str,
        line: int,
        *,
        borrow: bool = False,
        allow_error_value: bool = False,
    ) -> None:
        super().__init__(line)
        self.obj = obj
        self.attr = attr
        self.allow_error_value = allow_error_value
        assert isinstance(obj.type, RInstance), "Attribute access not supported: %s" % obj.type
        self.class_type = obj.type
        attr_type = obj.type.attr_type(attr)
        self.type = attr_type
        if allow_error_value:
            self.error_kind = ERR_NEVER
        elif attr_type.error_overlap:
            self.error_kind = ERR_MAGIC_OVERLAPPING
        self.is_borrowed = borrow and attr_type.is_refcounted

    def sources(self) -> list[Value]:
        return [self.obj]

    def set_sources(self, new: list[Value]) -> None:
        (self.obj,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_get_attr(self)


@final
class SetAttr(RegisterOp):
    """obj.attr = src (for a native object)"""

    error_kind = ERR_FALSE

    def __init__(self, obj: Value, attr: str, src: Value, line: int) -> None:
        super().__init__(line)
        self.obj = obj
        self.attr = attr
        self.src = src
        assert isinstance(obj.type, RInstance), "Attribute access not supported: %s" % obj.type
        self.class_type = obj.type
        self.type = bool_rprimitive
        # If True, we can safely assume that the attribute is previously undefined
        # and we don't use a setter
        self.is_init = False

        cl = self.class_type.class_ir
        is_propset = False
        for ir in cl.mro:
            propset = ir.method_decls.get(PROPSET_PREFIX + attr)
            if propset is not None:
                is_propset = not propset.implicit
                break
        # If True, this op represents calling a property setter.
        self.is_propset = is_propset

    def mark_as_initializer(self) -> None:
        self.is_init = True
        self.error_kind = ERR_NEVER
        self.type = void_rtype

    def sources(self) -> list[Value]:
        return [self.obj, self.src]

    def set_sources(self, new: list[Value]) -> None:
        self.obj, self.src = new

    def stolen(self) -> list[Value]:
        # The property setter method increfs the passed value so don't treat it as a steal
        # to avoid leaking.
        if self.is_propset:
            return []
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_set_attr(self)


# Default name space for statics, variables
NAMESPACE_STATIC: Final = "static"

# Static namespace for pointers to native type objects
NAMESPACE_TYPE: Final = "type"

# Namespace for modules
NAMESPACE_MODULE: Final = "module"

# Namespace for Python 3.12 type variable objects (implicitly created TypeVar instances, etc.)
NAMESPACE_TYPE_VAR: Final = "typevar"


@final
class LoadStatic(RegisterOp):
    """Load a static name (name :: static).

    Load a C static variable/pointer. The namespace for statics is shared
    for the entire compilation group. You can optionally provide a module
    name and a sub-namespace identifier for additional namespacing to avoid
    name conflicts. The static namespace does not overlap with other C names,
    since the final C name will get a prefix, so conflicts only must be
    avoided with other statics.
    """

    error_kind = ERR_NEVER
    is_borrowed = True

    def __init__(
        self,
        type: RType,
        identifier: str,
        module_name: str | None = None,
        namespace: str = NAMESPACE_STATIC,
        line: int = -1,
        ann: object = None,
    ) -> None:
        super().__init__(line)
        self.identifier = identifier
        self.module_name = module_name
        self.namespace = namespace
        self.type = type
        self.ann = ann  # An object to pretty print with the load

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_static(self)


@final
class InitStatic(RegisterOp):
    """static = value :: static

    Initialize a C static variable/pointer. See everything in LoadStatic.
    """

    error_kind = ERR_NEVER

    def __init__(
        self,
        value: Value,
        identifier: str,
        module_name: str | None = None,
        namespace: str = NAMESPACE_STATIC,
        line: int = -1,
    ) -> None:
        super().__init__(line)
        self.identifier = identifier
        self.module_name = module_name
        self.namespace = namespace
        self.value = value

    def sources(self) -> list[Value]:
        return [self.value]

    def set_sources(self, new: list[Value]) -> None:
        (self.value,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_init_static(self)


@final
class TupleSet(RegisterOp):
    """dest = (reg, ...) (for fixed-length tuple)"""

    error_kind = ERR_NEVER

    def __init__(self, items: list[Value], line: int) -> None:
        super().__init__(line)
        self.items = items
        # Don't keep track of the fact that an int is short after it
        # is put into a tuple, since we don't properly implement
        # runtime subtyping for tuples.
        self.tuple_type = RTuple(
            [
                arg.type if not is_short_int_rprimitive(arg.type) else int_rprimitive
                for arg in items
            ]
        )
        self.type = self.tuple_type

    def sources(self) -> list[Value]:
        return self.items.copy()

    def stolen(self) -> list[Value]:
        return self.items.copy()

    def set_sources(self, new: list[Value]) -> None:
        self.items = new[:]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_tuple_set(self)


@final
class TupleGet(RegisterOp):
    """Get item of a fixed-length tuple (src[index])."""

    error_kind = ERR_NEVER

    def __init__(self, src: Value, index: int, line: int = -1, *, borrow: bool = False) -> None:
        super().__init__(line)
        assert isinstance(
            src.type, RTuple
        ), f"TupleGet only operates on tuples, not {type(src.type).__name__}"
        src_len = len(src.type.types)
        self.src = src
        self.index = index
        if index < 0:
            self.index += src_len
        assert (
            self.index <= src_len - 1
        ), f"Index out of range.\nsource type: {src.type}\nindex: {index}"
        self.type = src.type.types[index]
        self.is_borrowed = borrow

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_tuple_get(self)


@final
class Cast(RegisterOp):
    """cast(type, src)

    Perform a runtime type check (no representation or value conversion).

    DO NOT increment reference counts.
    """

    error_kind = ERR_MAGIC

    def __init__(
        self, src: Value, typ: RType, line: int, *, borrow: bool = False, unchecked: bool = False
    ) -> None:
        super().__init__(line)
        self.src = src
        self.type = typ
        # If true, don't incref the result.
        self.is_borrowed = borrow
        # If true, don't perform a runtime type check (only changes the static type of
        # the operand). Used when we know that the cast will always succeed.
        self.is_unchecked = unchecked
        if unchecked:
            self.error_kind = ERR_NEVER

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        if self.is_borrowed:
            return []
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_cast(self)


@final
class Box(RegisterOp):
    """box(type, src)

    This converts from a potentially unboxed representation to a straight Python object.
    Only supported for types with an unboxed representation.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, line: int = -1) -> None:
        super().__init__(line)
        self.src = src
        self.type = object_rprimitive
        # When we box None and bool values, we produce a borrowed result
        if is_none_rprimitive(self.src.type) or is_bool_or_bit_rprimitive(self.src.type):
            self.is_borrowed = True

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_box(self)


@final
class Unbox(RegisterOp):
    """unbox(type, src)

    This is similar to a cast, but it also changes to a (potentially) unboxed runtime
    representation. Only supported for types with an unboxed representation.
    """

    def __init__(self, src: Value, typ: RType, line: int) -> None:
        self.src = src
        self.type = typ
        if not typ.error_overlap:
            self.error_kind = ERR_MAGIC
        else:
            self.error_kind = ERR_MAGIC_OVERLAPPING
        super().__init__(line)

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_unbox(self)


@final
class RaiseStandardError(RegisterOp):
    """Raise built-in exception with an optional error string.

    We have a separate opcode for this for convenience and to
    generate smaller, more idiomatic C code.
    """

    # TODO: Make it more explicit at IR level that this always raises

    error_kind = ERR_FALSE

    VALUE_ERROR: Final = "ValueError"
    ASSERTION_ERROR: Final = "AssertionError"
    STOP_ITERATION: Final = "StopIteration"
    UNBOUND_LOCAL_ERROR: Final = "UnboundLocalError"
    RUNTIME_ERROR: Final = "RuntimeError"
    NAME_ERROR: Final = "NameError"
    ZERO_DIVISION_ERROR: Final = "ZeroDivisionError"
    INDEX_ERROR: Final = "IndexError"

    def __init__(self, class_name: str, value: str | Value | None, line: int) -> None:
        super().__init__(line)
        self.class_name = class_name
        self.value = value
        self.type = bool_rprimitive

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_raise_standard_error(self)


# True steals all arguments, False steals none, a list steals those in matching positions
StealsDescription = bool | list[bool]


@final
class CallC(RegisterOp):
    """result = function(arg0, arg1, ...)

    Call a C function that is not a compiled/native function (for
    example, a Python C API function). Use Call to call native
    functions.
    """

    def __init__(
        self,
        function_name: str,
        args: list[Value],
        ret_type: RType,
        steals: StealsDescription,
        is_borrowed: bool,
        error_kind: int,
        line: int,
        var_arg_idx: int = -1,
        *,
        is_pure: bool = False,
        returns_null: bool = False,
        dependencies: list[Dependency] | None = None,
    ) -> None:
        self.error_kind = error_kind
        super().__init__(line)
        self.function_name = function_name
        self.args = args
        self.type = ret_type
        self.steals = steals
        self.is_borrowed = is_borrowed
        # The position of the first variable argument in args (if >= 0)
        self.var_arg_idx = var_arg_idx
        # Is the function pure? Pure functions have no side effects
        # and all the arguments are immutable. Pure functions support
        # additional optimizations. Pure functions never fail.
        self.is_pure = is_pure
        # The function might return a null value that does not indicate
        # an error.
        self.returns_null = returns_null
        # Dependencies (such as capsules) that must be imported and initialized before
        # calling this function (used for C functions exported from librt).
        self.dependencies = dependencies
        if is_pure or returns_null:
            assert error_kind == ERR_NEVER

    def sources(self) -> list[Value]:
        return self.args[:]

    def set_sources(self, new: list[Value]) -> None:
        self.args = new[:]

    def stolen(self) -> list[Value]:
        if isinstance(self.steals, list):
            assert len(self.steals) == len(self.args)
            return [arg for arg, steal in zip(self.args, self.steals) if steal]
        else:
            return [] if not self.steals else self.sources()

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_call_c(self)


@final
class Truncate(RegisterOp):
    """result = truncate src from src_type to dst_type

    Truncate a value from type with more bits to type with less bits.

    dst_type and src_type can be native integer types, bools or tagged
    integers. Tagged integers should have the tag bit unset.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, dst_type: RType, line: int = -1) -> None:
        super().__init__(line)
        self.src = src
        self.type = dst_type
        self.src_type = src.type

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        return []

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_truncate(self)


@final
class Extend(RegisterOp):
    """result = extend src from src_type to dst_type

    Extend a value from a type with fewer bits to a type with more bits.

    dst_type and src_type can be native integer types, bools or tagged
    integers. Tagged integers should have the tag bit unset.

    If 'signed' is true, perform sign extension. Otherwise, the result will be
    zero extended.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, dst_type: RType, signed: bool, line: int = -1) -> None:
        super().__init__(line)
        self.src = src
        self.type = dst_type
        self.src_type = src.type
        self.signed = signed

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        return []

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_extend(self)


@final
class LoadGlobal(RegisterOp):
    """Load a low-level global variable/pointer.

    Note that can't be used to directly load Python module-level
    global variable, since they are stored in a globals dictionary
    and accessed using dictionary operations.
    """

    error_kind = ERR_NEVER
    is_borrowed = True

    def __init__(self, type: RType, identifier: str, line: int = -1, ann: object = None) -> None:
        super().__init__(line)
        self.identifier = identifier
        self.type = type
        self.ann = ann  # An object to pretty print with the load

    def sources(self) -> list[Value]:
        return []

    def set_sources(self, new: list[Value]) -> None:
        assert not new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_global(self)


@final
class IntOp(RegisterOp):
    """Binary arithmetic or bitwise op on integer operands (e.g., r1 = r2 + r3).

    These ops are low-level and are similar to the corresponding C
    operations.

    The left and right values must have low-level integer types with
    compatible representations. Fixed-width integers, short_int_rprimitive,
    bool_rprimitive and bit_rprimitive are supported.

    For tagged (arbitrary-precision) integer ops look at mypyc.primitives.int_ops.
    """

    error_kind = ERR_NEVER

    # Arithmetic ops
    ADD: Final = 0
    SUB: Final = 1
    MUL: Final = 2
    DIV: Final = 3
    MOD: Final = 4

    # Bitwise ops
    AND: Final = 200
    OR: Final = 201
    XOR: Final = 202
    LEFT_SHIFT: Final = 203
    RIGHT_SHIFT: Final = 204

    op_str: Final = {
        ADD: "+",
        SUB: "-",
        MUL: "*",
        DIV: "/",
        MOD: "%",
        AND: "&",
        OR: "|",
        XOR: "^",
        LEFT_SHIFT: "<<",
        RIGHT_SHIFT: ">>",
    }

    def __init__(self, type: RType, lhs: Value, rhs: Value, op: int, line: int = -1) -> None:
        super().__init__(line)
        self.type = type
        self.lhs = lhs
        self.rhs = rhs
        self.op = op

    def sources(self) -> list[Value]:
        return [self.lhs, self.rhs]

    def set_sources(self, new: list[Value]) -> None:
        self.lhs, self.rhs = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_int_op(self)


# We can't have this in the IntOp class body, because of
# https://github.com/mypyc/mypyc/issues/932.
int_op_to_id: Final = {op: op_id for op_id, op in IntOp.op_str.items()}


@final
class ComparisonOp(RegisterOp):
    """Low-level comparison op for integers and pointers.

    Both unsigned and signed comparisons are supported. Supports
    comparisons between fixed-width integer types and pointer types.
    The operands should have matching sizes.

    The result is always a bit (representing a boolean).

    Python semantics, such as calling __eq__, are not supported.
    """

    # Must be ERR_NEVER or ERR_FALSE. ERR_FALSE means that a false result
    # indicates that an exception has been raised and should be propagated.
    error_kind = ERR_NEVER

    # S for signed and U for unsigned
    EQ: Final = 100
    NEQ: Final = 101
    SLT: Final = 102
    SGT: Final = 103
    SLE: Final = 104
    SGE: Final = 105
    ULT: Final = 106
    UGT: Final = 107
    ULE: Final = 108
    UGE: Final = 109

    op_str: Final = {
        EQ: "==",
        NEQ: "!=",
        SLT: "<",
        SGT: ">",
        SLE: "<=",
        SGE: ">=",
        ULT: "<",
        UGT: ">",
        ULE: "<=",
        UGE: ">=",
    }

    signed_ops: Final = {"==": EQ, "!=": NEQ, "<": SLT, ">": SGT, "<=": SLE, ">=": SGE}
    unsigned_ops: Final = {"==": EQ, "!=": NEQ, "<": ULT, ">": UGT, "<=": ULE, ">=": UGE}

    def __init__(self, lhs: Value, rhs: Value, op: int, line: int = -1) -> None:
        super().__init__(line)
        self.type = bit_rprimitive
        self.lhs = lhs
        self.rhs = rhs
        self.op = op

    def sources(self) -> list[Value]:
        return [self.lhs, self.rhs]

    def set_sources(self, new: list[Value]) -> None:
        self.lhs, self.rhs = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_comparison_op(self)


@final
class FloatOp(RegisterOp):
    """Binary float arithmetic op (e.g., r1 = r2 + r3).

    These ops are low-level and are similar to the corresponding C
    operations (and somewhat different from Python operations).

    The left and right values must be floats.
    """

    error_kind = ERR_NEVER

    ADD: Final = 0
    SUB: Final = 1
    MUL: Final = 2
    DIV: Final = 3
    MOD: Final = 4

    op_str: Final = {ADD: "+", SUB: "-", MUL: "*", DIV: "/", MOD: "%"}

    def __init__(self, lhs: Value, rhs: Value, op: int, line: int = -1) -> None:
        super().__init__(line)
        self.type = float_rprimitive
        self.lhs = lhs
        self.rhs = rhs
        self.op = op

    def sources(self) -> list[Value]:
        return [self.lhs, self.rhs]

    def set_sources(self, new: list[Value]) -> None:
        self.lhs, self.rhs = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_float_op(self)


# We can't have this in the FloatOp class body, because of
# https://github.com/mypyc/mypyc/issues/932.
float_op_to_id: Final = {op: op_id for op_id, op in FloatOp.op_str.items()}


@final
class FloatNeg(RegisterOp):
    """Float negation op (r1 = -r2)."""

    error_kind = ERR_NEVER

    def __init__(self, src: Value, line: int = -1) -> None:
        super().__init__(line)
        self.type = float_rprimitive
        self.src = src

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_float_neg(self)


@final
class FloatComparisonOp(RegisterOp):
    """Low-level comparison op for floats."""

    error_kind = ERR_NEVER

    EQ: Final = 200
    NEQ: Final = 201
    LT: Final = 202
    GT: Final = 203
    LE: Final = 204
    GE: Final = 205

    op_str: Final = {EQ: "==", NEQ: "!=", LT: "<", GT: ">", LE: "<=", GE: ">="}

    def __init__(self, lhs: Value, rhs: Value, op: int, line: int = -1) -> None:
        super().__init__(line)
        self.type = bit_rprimitive
        self.lhs = lhs
        self.rhs = rhs
        self.op = op

    def sources(self) -> list[Value]:
        return [self.lhs, self.rhs]

    def set_sources(self, new: list[Value]) -> None:
        self.lhs, self.rhs = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_float_comparison_op(self)


# We can't have this in the FloatOp class body, because of
# https://github.com/mypyc/mypyc/issues/932.
float_comparison_op_to_id: Final = {op: op_id for op_id, op in FloatComparisonOp.op_str.items()}


@final
class LoadMem(RegisterOp):
    """Read a memory location: result = *(type *)src.

    Attributes:
      type: Type of the read value
      src: Pointer to memory to read
    """

    error_kind = ERR_NEVER

    def __init__(self, type: RType, src: Value, line: int = -1, *, borrow: bool = False) -> None:
        super().__init__(line)
        self.type = type
        # TODO: Support other native integer types
        assert is_pointer_rprimitive(src.type)
        self.src = src
        self.is_borrowed = borrow and type.is_refcounted

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_mem(self)


@final
class SetMem(Op):
    """Write to a memory location: *(type *)dest = src

    Attributes:
      type: Type of the written value
      dest: Pointer to memory to write
      src: Source value
    """

    error_kind = ERR_NEVER

    def __init__(self, type: RType, dest: Value, src: Value, line: int = -1) -> None:
        super().__init__(line)
        self.type = void_rtype
        self.dest_type = type
        self.src = src
        self.dest = dest

    def sources(self) -> list[Value]:
        return [self.src, self.dest]

    def set_sources(self, new: list[Value]) -> None:
        self.src, self.dest = new

    def stolen(self) -> list[Value]:
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_set_mem(self)


@final
class GetElement(RegisterOp):
    """Get the value of a struct element from a struct value."""

    error_kind = ERR_NEVER
    is_borrowed = True

    def __init__(self, src: Value, field: str, line: int = -1) -> None:
        super().__init__(line)
        assert isinstance(src.type, (RStruct, RVec))
        self.type = src.type.field_type(field)
        self.src = src
        self.src_type = src.type
        self.field = field

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_get_element(self)


@final
class GetElementPtr(RegisterOp):
    """Get the address of a struct element from a pointer to a struct.

    If you have a struct value, use GetElement instead.

    Note that you may need to use KeepAlive to avoid the struct
    being freed, if it's reference counted, such as PyObject *.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, src_type: RType, field: str, line: int = -1) -> None:
        super().__init__(line)
        assert not isinstance(src.type, (RStruct, RVec))
        self.type = pointer_rprimitive
        self.src = src
        self.src_type = src_type
        self.field = field

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_get_element_ptr(self)


@final
class SetElement(RegisterOp):
    """Set the value of a struct element.

    This evaluates to a new struct with the changed value.

    Use together with Undef to initialize a fresh struct value
    (see Undef for more details).
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, field: str, item: Value, line: int = -1) -> None:
        super().__init__(line)
        assert isinstance(src.type, (RStruct, RVec)), src.type
        self.type = src.type
        self.src = src
        self.item = item
        self.field = field

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        return [self.src]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_set_element(self)


@final
class LoadAddress(RegisterOp):
    """Get the address of a value: result = (type)&src

    Attributes:
      type: Type of the loaded address(e.g. ptr/object_ptr)
      src: Source value (str for globals like 'PyList_Type',
           Register for temporary values or locals, LoadStatic
           for statics.)
    """

    error_kind = ERR_NEVER
    is_borrowed = True

    def __init__(self, type: RType, src: str | Register | LoadStatic, line: int = -1) -> None:
        super().__init__(line)
        self.type = type
        self.src = src

    def sources(self) -> list[Value]:
        if isinstance(self.src, Register):
            return [self.src]
        else:
            return []

    def set_sources(self, new: list[Value]) -> None:
        if new:
            assert isinstance(new[0], Register)
            assert len(new) == 1
            self.src = new[0]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_load_address(self)


@final
class KeepAlive(RegisterOp):
    """A no-op operation that ensures source values aren't freed.

    This is sometimes useful to avoid decref when a reference is still
    being held but not seen by the compiler.

    A typical use case is like this (C-like pseudocode):

      ptr = &x.item
      r = *ptr
      keep_alive x  # x must not be freed here
      # x may be freed here

    If we didn't have "keep_alive x", x could be freed immediately
    after taking the address of 'item', resulting in a read after free
    on the second line.

    If 'steal' is true, the value is considered to be stolen at
    this op, i.e. it won't be decref'd. You need to ensure that
    the value is freed otherwise, perhaps by using borrowing
    followed by Unborrow.

    Be careful with steal=True -- this can cause memory leaks.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: list[Value], line: int = -1, *, steal: bool = False) -> None:
        super().__init__(line)
        assert src
        self.src = src
        self.steal = steal

    def sources(self) -> list[Value]:
        return self.src.copy()

    def stolen(self) -> list[Value]:
        if self.steal:
            return self.src.copy()
        return []

    def set_sources(self, new: list[Value]) -> None:
        self.src = new[:]

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_keep_alive(self)


@final
class Unborrow(RegisterOp):
    """A no-op op to create a regular reference from a borrowed one.

    Borrowed references can only be used temporarily and the reference
    counts won't be managed. This value will be refcounted normally.

    This is mainly useful if you split an aggregate value, such as
    a tuple, into components using borrowed values (to avoid increfs),
    and want to treat the components as sharing the original managed
    reference. You'll also need to use KeepAlive with steal=True to
    "consume" the original tuple reference:

      # t is a 2-tuple
      r0 = borrow t[0]
      r1 = borrow t[1]
      keep_alive steal t
      r2 = unborrow r0
      r3 = unborrow r1
      # now (r2, r3) represent the tuple as separate items, that are
      # managed again. (Note we need to steal before unborrow, to avoid
      # refcount briefly touching zero if r2 or r3 are unused.)

    Be careful with this -- this can easily cause double freeing.
    """

    error_kind = ERR_NEVER

    def __init__(self, src: Value, line: int = -1) -> None:
        super().__init__(line)
        assert src.is_borrowed
        self.src = src
        self.type = src.type

    def sources(self) -> list[Value]:
        return [self.src]

    def set_sources(self, new: list[Value]) -> None:
        (self.src,) = new

    def stolen(self) -> list[Value]:
        return []

    def accept(self, visitor: OpVisitor[T]) -> T:
        return visitor.visit_unborrow(self)


@trait
class OpVisitor(Generic[T]):
    """Generic visitor over ops (uses the visitor design pattern)."""

    @abstractmethod
    def visit_goto(self, op: Goto) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_branch(self, op: Branch) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_return(self, op: Return) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_unreachable(self, op: Unreachable) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_assign(self, op: Assign) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_assign_multi(self, op: AssignMulti) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_error_value(self, op: LoadErrorValue) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_literal(self, op: LoadLiteral) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_get_attr(self, op: GetAttr) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_set_attr(self, op: SetAttr) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_static(self, op: LoadStatic) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_init_static(self, op: InitStatic) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_tuple_get(self, op: TupleGet) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_tuple_set(self, op: TupleSet) -> T:
        raise NotImplementedError

    def visit_inc_ref(self, op: IncRef) -> T:
        raise NotImplementedError

    def visit_dec_ref(self, op: DecRef) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_call(self, op: Call) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_method_call(self, op: MethodCall) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_cast(self, op: Cast) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_box(self, op: Box) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_unbox(self, op: Unbox) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_raise_standard_error(self, op: RaiseStandardError) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_call_c(self, op: CallC) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_primitive_op(self, op: PrimitiveOp) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_truncate(self, op: Truncate) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_extend(self, op: Extend) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_global(self, op: LoadGlobal) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_int_op(self, op: IntOp) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_comparison_op(self, op: ComparisonOp) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_float_op(self, op: FloatOp) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_float_neg(self, op: FloatNeg) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_float_comparison_op(self, op: FloatComparisonOp) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_mem(self, op: LoadMem) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_set_mem(self, op: SetMem) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_get_element(self, op: GetElement) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_get_element_ptr(self, op: GetElementPtr) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_set_element(self, op: SetElement) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_load_address(self, op: LoadAddress) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_keep_alive(self, op: KeepAlive) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_unborrow(self, op: Unborrow) -> T:
        raise NotImplementedError


# TODO: Should the following definition live somewhere else?


# We do a three-pass deserialization scheme in order to resolve name
# references.
#  1. Create an empty ClassIR for each class in an SCC.
#  2. Deserialize all of the functions, which can contain references
#     to ClassIRs in their types
#  3. Deserialize all of the classes, which contain lots of references
#     to the functions they contain. (And to other classes.)
#
# Note that this approach differs from how we deserialize ASTs in mypy itself,
# where everything is deserialized in one pass then a second pass cleans up
# 'cross_refs'. We don't follow that approach here because it seems to be more
# code for not a lot of gain since it is easy in mypyc to identify all the objects
# we might need to reference.
#
# Because of these references, we need to maintain maps from class
# names to ClassIRs and func IDs to FuncIRs.
#
# These are tracked in a DeserMaps which is passed to every
# deserialization function.
#
# (Serialization and deserialization *will* be used for incremental
# compilation but so far it is not hooked up to anything.)
class DeserMaps(NamedTuple):
    classes: dict[str, ClassIR]
    functions: dict[str, FuncIR]
