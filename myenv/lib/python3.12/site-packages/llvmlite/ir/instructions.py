"""
Implementation of LLVM IR instructions.
"""

from llvmlite.ir import types
from llvmlite.ir.values import (Block, Function, Value, NamedValue, Constant,
                                MetaDataArgument, MetaDataString, AttributeSet,
                                Undefined, ArgumentAttributes)
from llvmlite.ir._utils import _HasMetadata


class Instruction(NamedValue, _HasMetadata):
    def __init__(self, parent, typ, opname, operands, name='', flags=()):
        super(Instruction, self).__init__(parent, typ, name=name)
        assert isinstance(parent, Block)
        assert isinstance(flags, (tuple, list))
        self.opname = opname
        self.operands = operands
        self.flags = list(flags)
        self.metadata = {}

    @property
    def function(self):
        return self.parent.function

    @property
    def module(self):
        return self.parent.function.module

    def descr(self, buf):
        opname = self.opname
        if self.flags:
            opname = ' '.join([opname] + self.flags)
        operands = ', '.join([op.get_reference() for op in self.operands])
        typ = self.type
        metadata = self._stringify_metadata(leading_comma=True)
        buf.append("{0} {1} {2}{3}\n"
                   .format(opname, typ, operands, metadata))

    def replace_usage(self, old, new):
        if old in self.operands:
            ops = []
            for op in self.operands:
                ops.append(new if op is old else op)
            self.operands = tuple(ops)
            self._clear_string_cache()

    def __repr__(self):
        return "<ir.%s %r of type '%s', opname %r, operands %r>" % (
            self.__class__.__name__, self.name, self.type,
            self.opname, self.operands)


class CallInstrAttributes(AttributeSet):
    _known = frozenset(['convergent', 'noreturn', 'nounwind', 'readonly',
                        'readnone', 'noinline', 'alwaysinline'])


TailMarkerOptions = frozenset(['tail', 'musttail', 'notail'])


class FastMathFlags(AttributeSet):
    _known = frozenset(['fast', 'nnan', 'ninf', 'nsz', 'arcp', 'contract',
                        'afn', 'reassoc'])


class CallInstr(Instruction):
    def __init__(self, parent, func, args, name='', cconv=None, tail=None,
                 fastmath=(), attrs=(), arg_attrs=None):
        self.cconv = (func.calling_convention
                      if cconv is None and isinstance(func, Function)
                      else cconv)

        # For backwards compatibility with previous API of accepting a "truthy"
        # value for a hint to the optimizer to potentially tail optimize.
        if isinstance(tail, str) and tail in TailMarkerOptions:
            pass
        elif tail:
            tail = "tail"
        else:
            tail = ""

        self.tail = tail
        self.fastmath = FastMathFlags(fastmath)
        self.attributes = CallInstrAttributes(attrs)
        self.arg_attributes = {}
        if arg_attrs:
            for idx, attrs in arg_attrs.items():
                if not (0 <= idx < len(args)):
                    raise ValueError("Invalid argument index {}"
                                     .format(idx))
                self.arg_attributes[idx] = ArgumentAttributes(attrs)

        # Fix and validate arguments
        args = list(args)
        for i in range(len(func.function_type.args)):
            arg = args[i]
            expected_type = func.function_type.args[i]
            if (isinstance(expected_type, types.MetaDataType) and
                    arg.type != expected_type):
                arg = MetaDataArgument(arg)
            if arg.type != expected_type:
                msg = ("Type of #{0} arg mismatch: {1} != {2}"
                       .format(1 + i, expected_type, arg.type))
                raise TypeError(msg)
            args[i] = arg

        super(CallInstr, self).__init__(parent, func.function_type.return_type,
                                        "call", [func] + list(args), name=name)

    @property
    def callee(self):
        return self.operands[0]

    @callee.setter
    def callee(self, newcallee):
        self.operands[0] = newcallee

    @property
    def args(self):
        return self.operands[1:]

    def replace_callee(self, newfunc):
        if newfunc.function_type != self.callee.function_type:
            raise TypeError("New function has incompatible type")
        self.callee = newfunc

    @property
    def called_function(self):
        """The callee function"""
        return self.callee

    def _descr(self, buf, add_metadata):
        def descr_arg(i, a):
            if i in self.arg_attributes:
                attrs = ' '.join(self.arg_attributes[i]._to_list(a.type)) + ' '
            else:
                attrs = ''
            return '{0} {1}{2}'.format(a.type, attrs, a.get_reference())
        args = ', '.join([descr_arg(i, a) for i, a in enumerate(self.args)])

        fnty = self.callee.function_type
        # Only print function type if variable-argument
        if fnty.var_arg:
            ty = fnty
        # Otherwise, just print the return type.
        else:
            # Fastmath flag work only in this case
            ty = fnty.return_type
        callee_ref = "{0} {1}".format(ty, self.callee.get_reference())
        if self.cconv:
            callee_ref = "{0} {1}".format(self.cconv, callee_ref)

        tail_marker = ""
        if self.tail:
            tail_marker = "{0} ".format(self.tail)

        fn_attrs = ' ' + ' '.join(self.attributes._to_list(fnty.return_type))\
            if self.attributes else ''

        fm_attrs = ' ' + ' '.join(self.fastmath._to_list(fnty.return_type))\
            if self.fastmath else ''

        buf.append("{tail}{op}{fastmath} {callee}({args}){attr}{meta}\n".format(
            tail=tail_marker,
            op=self.opname,
            callee=callee_ref,
            fastmath=fm_attrs,
            args=args,
            attr=fn_attrs,
            meta=(self._stringify_metadata(leading_comma=True)
                  if add_metadata else ""),
        ))

    def descr(self, buf):
        self._descr(buf, add_metadata=True)


class InvokeInstr(CallInstr):
    def __init__(self, parent, func, args, normal_to, unwind_to, name='',
                 cconv=None, fastmath=(), attrs=(), arg_attrs=None):
        assert isinstance(normal_to, Block)
        assert isinstance(unwind_to, Block)
        super(InvokeInstr, self).__init__(parent, func, args, name, cconv,
                                          tail=False, fastmath=fastmath,
                                          attrs=attrs, arg_attrs=arg_attrs)
        self.opname = "invoke"
        self.normal_to = normal_to
        self.unwind_to = unwind_to

    def descr(self, buf):
        super(InvokeInstr, self)._descr(buf, add_metadata=False)
        buf.append("      to label {0} unwind label {1}{metadata}\n".format(
            self.normal_to.get_reference(),
            self.unwind_to.get_reference(),
            metadata=self._stringify_metadata(leading_comma=True),
        ))


class Terminator(Instruction):
    def __init__(self, parent, opname, operands):
        super(Terminator, self).__init__(parent, types.VoidType(), opname,
                                         operands)

    def descr(self, buf):
        opname = self.opname
        operands = ', '.join(["{0} {1}".format(op.type, op.get_reference())
                              for op in self.operands])
        metadata = self._stringify_metadata(leading_comma=True)
        buf.append("{0} {1}{2}".format(opname, operands, metadata))


class PredictableInstr(Instruction):

    def set_weights(self, weights):
        operands = [MetaDataString(self.module, "branch_weights")]
        for w in weights:
            if w < 0:
                raise ValueError("branch weight must be a positive integer")
            operands.append(Constant(types.IntType(32), w))
        md = self.module.add_metadata(operands)
        self.set_metadata("prof", md)


class Ret(Terminator):
    def __init__(self, parent, opname, return_value=None):
        operands = [return_value] if return_value is not None else []
        super(Ret, self).__init__(parent, opname, operands)

    @property
    def return_value(self):
        if self.operands:
            return self.operands[0]
        else:
            return None

    def descr(self, buf):
        return_value = self.return_value
        metadata = self._stringify_metadata(leading_comma=True)
        if return_value is not None:
            buf.append("{0} {1} {2}{3}\n"
                       .format(self.opname, return_value.type,
                               return_value.get_reference(),
                               metadata))
        else:
            buf.append("{0}{1}\n".format(self.opname, metadata))


class Branch(Terminator):
    pass


class ConditionalBranch(PredictableInstr, Terminator):
    pass


class IndirectBranch(PredictableInstr, Terminator):
    def __init__(self, parent, opname, addr):
        super(IndirectBranch, self).__init__(parent, opname, [addr])
        self.destinations = []

    @property
    def address(self):
        return self.operands[0]

    def add_destination(self, block):
        assert isinstance(block, Block)
        self.destinations.append(block)

    def descr(self, buf):
        destinations = ["label {0}".format(blk.get_reference())
                        for blk in self.destinations]
        buf.append("indirectbr {0} {1}, [{2}]  {3}\n".format(
            self.address.type,
            self.address.get_reference(),
            ', '.join(destinations),
            self._stringify_metadata(leading_comma=True),
        ))


class SwitchInstr(PredictableInstr, Terminator):

    def __init__(self, parent, opname, val, default):
        super(SwitchInstr, self).__init__(parent, opname, [val])
        self.default = default
        self.cases = []

    @property
    def value(self):
        return self.operands[0]

    def add_case(self, val, block):
        assert isinstance(block, Block)
        if not isinstance(val, Value):
            val = Constant(self.value.type, val)
        self.cases.append((val, block))

    def descr(self, buf):
        cases = ["{0} {1}, label {2}".format(val.type, val.get_reference(),
                                             blk.get_reference())
                 for val, blk in self.cases]
        buf.append("switch {0} {1}, label {2} [{3}]  {4}\n".format(
            self.value.type,
            self.value.get_reference(),
            self.default.get_reference(),
            ' '.join(cases),
            self._stringify_metadata(leading_comma=True),
        ))


class Resume(Terminator):
    pass


class SelectInstr(Instruction):
    def __init__(self, parent, cond, lhs, rhs, name='', flags=()):
        assert lhs.type == rhs.type
        super(SelectInstr, self).__init__(parent, lhs.type, "select",
                                          [cond, lhs, rhs], name=name,
                                          flags=flags)

    @property
    def cond(self):
        return self.operands[0]

    @property
    def lhs(self):
        return self.operands[1]

    @property
    def rhs(self):
        return self.operands[2]

    def descr(self, buf):
        buf.append("select {0} {1} {2}, {3} {4}, {5} {6} {7}\n".format(
            ' '.join(self.flags),
            self.cond.type, self.cond.get_reference(),
            self.lhs.type, self.lhs.get_reference(),
            self.rhs.type, self.rhs.get_reference(),
            self._stringify_metadata(leading_comma=True),
        ))


class CompareInstr(Instruction):
    # Define the following in subclasses
    OPNAME = 'invalid-compare'
    VALID_OP = {}

    def __init__(self, parent, op, lhs, rhs, name='', flags=[]):
        if op not in self.VALID_OP:
            raise ValueError("invalid comparison %r for %s" % (op, self.OPNAME))
        for flag in flags:
            if flag not in self.VALID_FLAG:
                raise ValueError("invalid flag %r for %s" % (flag, self.OPNAME))
        opname = self.OPNAME
        if isinstance(lhs.type, types.VectorType):
            typ = types.VectorType(types.IntType(1), lhs.type.count)
        else:
            typ = types.IntType(1)
        super(CompareInstr, self).__init__(parent, typ,
                                           opname, [lhs, rhs], flags=flags,
                                           name=name)
        self.op = op

    def descr(self, buf):
        buf.append("{opname}{flags} {op} {ty} {lhs}, {rhs} {meta}\n".format(
            opname=self.opname,
            flags=''.join(' ' + it for it in self.flags),
            op=self.op,
            ty=self.operands[0].type,
            lhs=self.operands[0].get_reference(),
            rhs=self.operands[1].get_reference(),
            meta=self._stringify_metadata(leading_comma=True),
        ))


class ICMPInstr(CompareInstr):
    OPNAME = 'icmp'
    VALID_OP = {
        'eq': 'equal',
        'ne': 'not equal',
        'ugt': 'unsigned greater than',
        'uge': 'unsigned greater or equal',
        'ult': 'unsigned less than',
        'ule': 'unsigned less or equal',
        'sgt': 'signed greater than',
        'sge': 'signed greater or equal',
        'slt': 'signed less than',
        'sle': 'signed less or equal',
    }
    VALID_FLAG = set()


class FCMPInstr(CompareInstr):
    OPNAME = 'fcmp'
    VALID_OP = {
        'false': 'no comparison, always returns false',
        'oeq': 'ordered and equal',
        'ogt': 'ordered and greater than',
        'oge': 'ordered and greater than or equal',
        'olt': 'ordered and less than',
        'ole': 'ordered and less than or equal',
        'one': 'ordered and not equal',
        'ord': 'ordered (no nans)',
        'ueq': 'unordered or equal',
        'ugt': 'unordered or greater than',
        'uge': 'unordered or greater than or equal',
        'ult': 'unordered or less than',
        'ule': 'unordered or less than or equal',
        'une': 'unordered or not equal',
        'uno': 'unordered (either nans)',
        'true': 'no comparison, always returns true',
    }
    VALID_FLAG = {'nnan', 'ninf', 'nsz', 'arcp', 'contract', 'afn', 'reassoc',
                  'fast'}


class CastInstr(Instruction):
    def __init__(self, parent, op, val, typ, name=''):
        super(CastInstr, self).__init__(parent, typ, op, [val], name=name)

    def descr(self, buf):
        buf.append("{0} {1} {2} to {3} {4}\n".format(
            self.opname,
            self.operands[0].type,
            self.operands[0].get_reference(),
            self.type,
            self._stringify_metadata(leading_comma=True),
        ))


class LoadInstr(Instruction):

    def __init__(self, parent, ptr, name=''):
        super(LoadInstr, self).__init__(parent, ptr.type.pointee, "load",
                                        [ptr], name=name)
        self.align = None

    def descr(self, buf):
        [val] = self.operands
        if self.align is not None:
            align = ', align %d' % (self.align)
        else:
            align = ''
        buf.append("load {0}, {1} {2}{3}{4}\n".format(
            val.type.pointee,
            val.type,
            val.get_reference(),
            align,
            self._stringify_metadata(leading_comma=True),
        ))


class StoreInstr(Instruction):
    def __init__(self, parent, val, ptr):
        super(StoreInstr, self).__init__(parent, types.VoidType(), "store",
                                         [val, ptr])

    def descr(self, buf):
        val, ptr = self.operands
        if self.align is not None:
            align = ', align %d' % (self.align)
        else:
            align = ''
        buf.append("store {0} {1}, {2} {3}{4}{5}\n".format(
            val.type,
            val.get_reference(),
            ptr.type,
            ptr.get_reference(),
            align,
            self._stringify_metadata(leading_comma=True),
        ))


class LoadAtomicInstr(Instruction):
    def __init__(self, parent, ptr, ordering, align, name=''):
        super(LoadAtomicInstr, self).__init__(parent, ptr.type.pointee,
                                              "load atomic", [ptr], name=name)
        self.ordering = ordering
        self.align = align

    def descr(self, buf):
        [val] = self.operands
        buf.append("load atomic {0}, {1} {2} {3}, align {4}{5}\n".format(
            val.type.pointee,
            val.type,
            val.get_reference(),
            self.ordering,
            self.align,
            self._stringify_metadata(leading_comma=True),
        ))


class StoreAtomicInstr(Instruction):
    def __init__(self, parent, val, ptr, ordering, align):
        super(StoreAtomicInstr, self).__init__(parent, types.VoidType(),
                                               "store atomic", [val, ptr])
        self.ordering = ordering
        self.align = align

    def descr(self, buf):
        val, ptr = self.operands
        buf.append("store atomic {0} {1}, {2} {3} {4}, align {5}{6}\n".format(
            val.type,
            val.get_reference(),
            ptr.type,
            ptr.get_reference(),
            self.ordering,
            self.align,
            self._stringify_metadata(leading_comma=True),
        ))


class AllocaInstr(Instruction):
    def __init__(self, parent, typ, count, name):
        operands = [count] if count else ()
        super(AllocaInstr, self).__init__(parent, typ.as_pointer(), "alloca",
                                          operands, name)
        self.align = None

    def descr(self, buf):
        buf.append("{0} {1}".format(self.opname, self.type.pointee))
        if self.operands:
            op, = self.operands
            buf.append(", {0} {1}".format(op.type, op.get_reference()))
        if self.align is not None:
            buf.append(", align {0}".format(self.align))
        if self.metadata:
            buf.append(self._stringify_metadata(leading_comma=True))


class GEPInstr(Instruction):
    def __init__(self, parent, ptr, indices, inbounds, name):
        typ = ptr.type
        lasttyp = None
        lastaddrspace = 0
        for i in indices:
            lasttyp, typ = typ, typ.gep(i)
            # inherit the addrspace from the last seen pointer
            if isinstance(lasttyp, types.PointerType):
                lastaddrspace = lasttyp.addrspace

        if (not isinstance(typ, types.PointerType) and
                isinstance(lasttyp, types.PointerType)):
            typ = lasttyp
        else:
            typ = typ.as_pointer(lastaddrspace)

        super(GEPInstr, self).__init__(parent, typ, "getelementptr",
                                       [ptr] + list(indices), name=name)
        self.pointer = ptr
        self.indices = indices
        self.inbounds = inbounds

    def descr(self, buf):
        indices = ['{0} {1}'.format(i.type, i.get_reference())
                   for i in self.indices]
        op = "getelementptr inbounds" if self.inbounds else "getelementptr"
        buf.append("{0} {1}, {2} {3}, {4} {5}\n".format(
                   op,
                   self.pointer.type.pointee,
                   self.pointer.type,
                   self.pointer.get_reference(),
                   ', '.join(indices),
                   self._stringify_metadata(leading_comma=True),
                   ))


class PhiInstr(Instruction):
    def __init__(self, parent, typ, name, flags=()):
        super(PhiInstr, self).__init__(parent, typ, "phi", (), name=name,
                                       flags=flags)
        self.incomings = []

    def descr(self, buf):
        incs = ', '.join('[{0}, {1}]'.format(v.get_reference(),
                                             b.get_reference())
                         for v, b in self.incomings)
        buf.append("phi {0} {1} {2} {3}\n".format(
                   ' '.join(self.flags),
                   self.type,
                   incs,
                   self._stringify_metadata(leading_comma=True),
                   ))

    def add_incoming(self, value, block):
        assert isinstance(block, Block)
        self.incomings.append((value, block))

    def replace_usage(self, old, new):
        self.incomings = [((new if val is old else val), blk)
                          for (val, blk) in self.incomings]


class ExtractElement(Instruction):
    def __init__(self, parent, vector, index, name=''):
        if not isinstance(vector.type, types.VectorType):
            raise TypeError("vector needs to be of VectorType.")
        if not isinstance(index.type, types.IntType):
            raise TypeError("index needs to be of IntType.")
        typ = vector.type.element
        super(ExtractElement, self).__init__(parent, typ, "extractelement",
                                             [vector, index], name=name)

    def descr(self, buf):
        operands = ", ".join("{0} {1}".format(
                   op.type, op.get_reference()) for op in self.operands)
        buf.append("{opname} {operands}\n".format(
                   opname=self.opname, operands=operands))


class InsertElement(Instruction):
    def __init__(self, parent, vector, value, index, name=''):
        if not isinstance(vector.type, types.VectorType):
            raise TypeError("vector needs to be of VectorType.")
        if not value.type == vector.type.element:
            raise TypeError(
                "value needs to be of type {} not {}.".format(
                    vector.type.element, value.type))
        if not isinstance(index.type, types.IntType):
            raise TypeError("index needs to be of IntType.")
        typ = vector.type
        super(InsertElement, self).__init__(parent, typ, "insertelement",
                                            [vector, value, index], name=name)

    def descr(self, buf):
        operands = ", ".join("{0} {1}".format(
                   op.type, op.get_reference()) for op in self.operands)
        buf.append("{opname} {operands}\n".format(
                   opname=self.opname, operands=operands))


class ShuffleVector(Instruction):
    def __init__(self, parent, vector1, vector2, mask, name=''):
        if not isinstance(vector1.type, types.VectorType):
            raise TypeError("vector1 needs to be of VectorType.")
        if vector2 != Undefined:
            if vector2.type != vector1.type:
                raise TypeError("vector2 needs to be " +
                                "Undefined or of the same type as vector1.")
        if (not isinstance(mask, Constant) or
            not isinstance(mask.type, types.VectorType) or
            not (isinstance(mask.type.element, types.IntType) and
                 mask.type.element.width == 32)):
            raise TypeError("mask needs to be a constant i32 vector.")
        typ = types.VectorType(vector1.type.element, mask.type.count)
        index_range = range(vector1.type.count
                            if vector2 == Undefined
                            else 2 * vector1.type.count)
        if not all(ii.constant in index_range for ii in mask.constant):
            raise IndexError(
                "mask values need to be in {0}".format(index_range),
            )
        super(ShuffleVector, self).__init__(parent, typ, "shufflevector",
                                            [vector1, vector2, mask], name=name)

    def descr(self, buf):
        buf.append("shufflevector {0} {1}\n".format(
                   ", ".join("{0} {1}".format(op.type, op.get_reference())
                             for op in self.operands),
                   self._stringify_metadata(leading_comma=True),
                   ))


class ExtractValue(Instruction):
    def __init__(self, parent, agg, indices, name=''):
        typ = agg.type
        try:
            for i in indices:
                typ = typ.elements[i]
        except (AttributeError, IndexError):
            raise TypeError("Can't index at %r in %s"
                            % (list(indices), agg.type))

        super(ExtractValue, self).__init__(parent, typ, "extractvalue",
                                           [agg], name=name)

        self.aggregate = agg
        self.indices = indices

    def descr(self, buf):
        indices = [str(i) for i in self.indices]

        buf.append("extractvalue {0} {1}, {2} {3}\n".format(
                   self.aggregate.type,
                   self.aggregate.get_reference(),
                   ', '.join(indices),
                   self._stringify_metadata(leading_comma=True),
                   ))


class InsertValue(Instruction):
    def __init__(self, parent, agg, elem, indices, name=''):
        typ = agg.type
        try:
            for i in indices:
                typ = typ.elements[i]
        except (AttributeError, IndexError):
            raise TypeError("Can't index at %r in %s"
                            % (list(indices), agg.type))
        if elem.type != typ:
            raise TypeError("Can only insert %s at %r in %s: got %s"
                            % (typ, list(indices), agg.type, elem.type))
        super(InsertValue, self).__init__(parent, agg.type, "insertvalue",
                                          [agg, elem], name=name)

        self.aggregate = agg
        self.value = elem
        self.indices = indices

    def descr(self, buf):
        indices = [str(i) for i in self.indices]

        buf.append("insertvalue {0} {1}, {2} {3}, {4} {5}\n".format(
                   self.aggregate.type, self.aggregate.get_reference(),
                   self.value.type, self.value.get_reference(),
                   ', '.join(indices),
                   self._stringify_metadata(leading_comma=True),
                   ))


class Unreachable(Instruction):
    def __init__(self, parent):
        super(Unreachable, self).__init__(parent, types.VoidType(),
                                          "unreachable", (), name='')

    def descr(self, buf):
        buf += (self.opname, "\n")


class InlineAsm(object):
    def __init__(self, ftype, asm, constraint, side_effect=False):
        self.type = ftype.return_type
        self.function_type = ftype
        self.asm = asm
        self.constraint = constraint
        self.side_effect = side_effect

    def descr(self, buf):
        sideeffect = 'sideeffect' if self.side_effect else ''
        fmt = 'asm {sideeffect} "{asm}", "{constraint}"\n'
        buf.append(fmt.format(sideeffect=sideeffect, asm=self.asm,
                              constraint=self.constraint))

    def get_reference(self):
        buf = []
        self.descr(buf)
        return "".join(buf)

    def __str__(self):
        return "{0} {1}".format(self.type, self.get_reference())


class AtomicRMW(Instruction):
    def __init__(self, parent, op, ptr, val, ordering, name):
        super(AtomicRMW, self).__init__(parent, val.type, "atomicrmw",
                                        (ptr, val), name=name)
        self.operation = op
        self.ordering = ordering

    def descr(self, buf):
        ptr, val = self.operands
        fmt = ("atomicrmw {op} {ptrty} {ptr}, {valty} {val} {ordering} "
               "{metadata}\n")
        buf.append(fmt.format(op=self.operation,
                              ptrty=ptr.type,
                              ptr=ptr.get_reference(),
                              valty=val.type,
                              val=val.get_reference(),
                              ordering=self.ordering,
                              metadata=self._stringify_metadata(
                                  leading_comma=True),
                              ))


class CmpXchg(Instruction):
    """This instruction has changed since llvm3.5.  It is not compatible with
    older llvm versions.
    """

    def __init__(self, parent, ptr, cmp, val, ordering, failordering, name):
        outtype = types.LiteralStructType([val.type, types.IntType(1)])
        super(CmpXchg, self).__init__(parent, outtype, "cmpxchg",
                                      (ptr, cmp, val), name=name)
        self.ordering = ordering
        self.failordering = failordering

    def descr(self, buf):
        ptr, cmpval, val = self.operands
        fmt = "cmpxchg {ptrty} {ptr}, {ty} {cmp}, {ty} {val} {ordering} " \
              "{failordering} {metadata}\n"
        buf.append(fmt.format(ptrty=ptr.type,
                              ptr=ptr.get_reference(),
                              ty=cmpval.type,
                              cmp=cmpval.get_reference(),
                              val=val.get_reference(),
                              ordering=self.ordering,
                              failordering=self.failordering,
                              metadata=self._stringify_metadata(
                                  leading_comma=True),
                              ))


class _LandingPadClause(object):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return "{kind} {type} {value}".format(
            kind=self.kind,
            type=self.value.type,
            value=self.value.get_reference())


class CatchClause(_LandingPadClause):
    kind = 'catch'


class FilterClause(_LandingPadClause):
    kind = 'filter'

    def __init__(self, value):
        assert isinstance(value, Constant)
        assert isinstance(value.type, types.ArrayType)
        super(FilterClause, self).__init__(value)


class LandingPadInstr(Instruction):
    def __init__(self, parent, typ, name='', cleanup=False):
        super(LandingPadInstr, self).__init__(parent, typ, "landingpad", [],
                                              name=name)
        self.cleanup = cleanup
        self.clauses = []

    def add_clause(self, clause):
        assert isinstance(clause, _LandingPadClause)
        self.clauses.append(clause)

    def descr(self, buf):
        fmt = "landingpad {type}{cleanup}{clauses}\n"
        buf.append(fmt.format(type=self.type,
                              cleanup=' cleanup' if self.cleanup else '',
                              clauses=''.join(["\n      {0}".format(clause)
                                               for clause in self.clauses]),
                              ))


class Fence(Instruction):
    """
    The `fence` instruction.

    As of LLVM 5.0.1:

    fence [syncscope("<target-scope>")] <ordering>  ; yields void
    """

    VALID_FENCE_ORDERINGS = {"acquire", "release", "acq_rel", "seq_cst"}

    def __init__(self, parent, ordering, targetscope=None, name=''):
        super(Fence, self).__init__(parent, types.VoidType(), "fence", (),
                                    name=name)
        if ordering not in self.VALID_FENCE_ORDERINGS:
            msg = "Invalid fence ordering \"{0}\"! Should be one of {1}."
            raise ValueError(msg .format(ordering,
                                         ", ".join(self.VALID_FENCE_ORDERINGS)))
        self.ordering = ordering
        self.targetscope = targetscope

    def descr(self, buf):
        if self.targetscope is None:
            syncscope = ""
        else:
            syncscope = 'syncscope("{0}") '.format(self.targetscope)

        fmt = "fence {syncscope}{ordering}\n"
        buf.append(fmt.format(syncscope=syncscope,
                              ordering=self.ordering,
                              ))


class Comment(Instruction):
    """
    A line comment.
    """

    def __init__(self, parent, text):
        super(Comment, self).__init__(parent, types.VoidType(), ";", (),
                                      name='')
        assert "\n" not in text, "Comment cannot contain new line"
        self.text = text

    def descr(self, buf):
        buf.append(f"; {self.text}")
