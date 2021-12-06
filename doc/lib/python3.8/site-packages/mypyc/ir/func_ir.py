"""Intermediate representation of functions."""

from typing import List, Optional, Sequence
from typing_extensions import Final

from mypy.nodes import FuncDef, Block, ARG_POS, ARG_OPT, ARG_NAMED_OPT

from mypyc.common import JsonDict
from mypyc.ir.ops import (
    DeserMaps, BasicBlock, Value, Register, Assign, AssignMulti, ControlOp, LoadAddress
)
from mypyc.ir.rtypes import RType, deserialize_type
from mypyc.namegen import NameGenerator


class RuntimeArg:
    """Description of a function argument in IR.

    Argument kind is one of ARG_* constants defined in mypy.nodes.
    """

    def __init__(self, name: str, typ: RType, kind: int = ARG_POS) -> None:
        self.name = name
        self.type = typ
        self.kind = kind

    @property
    def optional(self) -> bool:
        return self.kind == ARG_OPT or self.kind == ARG_NAMED_OPT

    def __repr__(self) -> str:
        return 'RuntimeArg(name=%s, type=%s, optional=%r)' % (self.name, self.type, self.optional)

    def serialize(self) -> JsonDict:
        return {'name': self.name, 'type': self.type.serialize(), 'kind': self.kind}

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> 'RuntimeArg':
        return RuntimeArg(
            data['name'],
            deserialize_type(data['type'], ctx),
            data['kind'],
        )


class FuncSignature:
    """Signature of a function in IR."""

    # TODO: Track if method?

    def __init__(self, args: Sequence[RuntimeArg], ret_type: RType) -> None:
        self.args = tuple(args)
        self.ret_type = ret_type

    def __repr__(self) -> str:
        return 'FuncSignature(args=%r, ret=%r)' % (self.args, self.ret_type)

    def serialize(self) -> JsonDict:
        return {'args': [t.serialize() for t in self.args], 'ret_type': self.ret_type.serialize()}

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> 'FuncSignature':
        return FuncSignature(
            [RuntimeArg.deserialize(arg, ctx) for arg in data['args']],
            deserialize_type(data['ret_type'], ctx),
        )


FUNC_NORMAL = 0  # type: Final
FUNC_STATICMETHOD = 1  # type: Final
FUNC_CLASSMETHOD = 2  # type: Final


class FuncDecl:
    """Declaration of a function in IR (without body or implementation).

    A function can be a regular module-level function, a method, a
    static method, a class method, or a property getter/setter.
    """

    def __init__(self,
                 name: str,
                 class_name: Optional[str],
                 module_name: str,
                 sig: FuncSignature,
                 kind: int = FUNC_NORMAL,
                 is_prop_setter: bool = False,
                 is_prop_getter: bool = False) -> None:
        self.name = name
        self.class_name = class_name
        self.module_name = module_name
        self.sig = sig
        self.kind = kind
        self.is_prop_setter = is_prop_setter
        self.is_prop_getter = is_prop_getter
        if class_name is None:
            self.bound_sig = None  # type: Optional[FuncSignature]
        else:
            if kind == FUNC_STATICMETHOD:
                self.bound_sig = sig
            else:
                self.bound_sig = FuncSignature(sig.args[1:], sig.ret_type)

    @staticmethod
    def compute_shortname(class_name: Optional[str], name: str) -> str:
        return class_name + '.' + name if class_name else name

    @property
    def shortname(self) -> str:
        return FuncDecl.compute_shortname(self.class_name, self.name)

    @property
    def fullname(self) -> str:
        return self.module_name + '.' + self.shortname

    def cname(self, names: NameGenerator) -> str:
        return names.private_name(self.module_name, self.shortname)

    def serialize(self) -> JsonDict:
        return {
            'name': self.name,
            'class_name': self.class_name,
            'module_name': self.module_name,
            'sig': self.sig.serialize(),
            'kind': self.kind,
            'is_prop_setter': self.is_prop_setter,
            'is_prop_getter': self.is_prop_getter,
        }

    @staticmethod
    def get_name_from_json(f: JsonDict) -> str:
        return f['module_name'] + '.' + FuncDecl.compute_shortname(f['class_name'], f['name'])

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> 'FuncDecl':
        return FuncDecl(
            data['name'],
            data['class_name'],
            data['module_name'],
            FuncSignature.deserialize(data['sig'], ctx),
            data['kind'],
            data['is_prop_setter'],
            data['is_prop_getter'],
        )


class FuncIR:
    """Intermediate representation of a function with contextual information.

    Unlike FuncDecl, this includes the IR of the body (basic blocks).
    """

    def __init__(self,
                 decl: FuncDecl,
                 arg_regs: List[Register],
                 blocks: List[BasicBlock],
                 line: int = -1,
                 traceback_name: Optional[str] = None) -> None:
        # Declaration of the function, including the signature
        self.decl = decl
        # Registers for all the arguments to the function
        self.arg_regs = arg_regs
        # Body of the function
        self.blocks = blocks
        self.line = line
        # The name that should be displayed for tracebacks that
        # include this function. Function will be omitted from
        # tracebacks if None.
        self.traceback_name = traceback_name

    @property
    def args(self) -> Sequence[RuntimeArg]:
        return self.decl.sig.args

    @property
    def ret_type(self) -> RType:
        return self.decl.sig.ret_type

    @property
    def class_name(self) -> Optional[str]:
        return self.decl.class_name

    @property
    def sig(self) -> FuncSignature:
        return self.decl.sig

    @property
    def name(self) -> str:
        return self.decl.name

    @property
    def fullname(self) -> str:
        return self.decl.fullname

    def cname(self, names: NameGenerator) -> str:
        return self.decl.cname(names)

    def __repr__(self) -> str:
        if self.class_name:
            return '<FuncIR {}.{}>'.format(self.class_name, self.name)
        else:
            return '<FuncIR {}>'.format(self.name)

    def serialize(self) -> JsonDict:
        # We don't include blocks in the serialized version
        return {
            'decl': self.decl.serialize(),
            'line': self.line,
            'traceback_name': self.traceback_name,
        }

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> 'FuncIR':
        return FuncIR(
            FuncDecl.deserialize(data['decl'], ctx),
            [],
            [],
            data['line'],
            data['traceback_name'],
        )


INVALID_FUNC_DEF = FuncDef('<INVALID_FUNC_DEF>', [], Block([]))  # type: Final


def all_values(args: List[Register], blocks: List[BasicBlock]) -> List[Value]:
    """Return the set of all values that may be initialized in the blocks.

    This omits registers that are only read.
    """
    values = list(args)  # type: List[Value]
    seen_registers = set(args)

    for block in blocks:
        for op in block.ops:
            if not isinstance(op, ControlOp):
                if isinstance(op, (Assign, AssignMulti)):
                    if op.dest not in seen_registers:
                        values.append(op.dest)
                        seen_registers.add(op.dest)
                elif op.is_void:
                    continue
                else:
                    # If we take the address of a register, it might get initialized.
                    if (isinstance(op, LoadAddress)
                            and isinstance(op.src, Register)
                            and op.src not in seen_registers):
                        values.append(op.src)
                        seen_registers.add(op.src)
                    values.append(op)

    return values


def all_values_full(args: List[Register], blocks: List[BasicBlock]) -> List[Value]:
    """Return set of all values that are initialized or accessed."""
    values = list(args)  # type: List[Value]
    seen_registers = set(args)

    for block in blocks:
        for op in block.ops:
            for source in op.sources():
                # Look for unitialized registers that are accessed. Ignore
                # non-registers since we don't allow ops outside basic blocks.
                if isinstance(source, Register) and source not in seen_registers:
                    values.append(source)
                    seen_registers.add(source)
            if not isinstance(op, ControlOp):
                if isinstance(op, (Assign, AssignMulti)):
                    if op.dest not in seen_registers:
                        values.append(op.dest)
                        seen_registers.add(op.dest)
                elif op.is_void:
                    continue
                else:
                    values.append(op)

    return values
