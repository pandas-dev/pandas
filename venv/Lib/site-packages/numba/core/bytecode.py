from collections import namedtuple, OrderedDict
import dis
import inspect
import itertools

from types import CodeType, ModuleType

from numba.core import errors, utils, serialize
from numba.core.utils import PYVERSION


if PYVERSION in ((3, 12), (3, 13)):
    from opcode import _inline_cache_entries
    # Instruction/opcode length in bytes
    INSTR_LEN = 2
elif PYVERSION in ((3, 10), (3, 11)):
    pass
else:
    raise NotImplementedError(PYVERSION)


opcode_info = namedtuple('opcode_info', ['argsize'])
_ExceptionTableEntry = namedtuple("_ExceptionTableEntry",
                                  "start end target depth lasti")

# The following offset is used as a hack to inject a NOP at the start of the
# bytecode. So that function starting with `while True` will not have block-0
# as a jump target. The Lowerer puts argument initialization at block-0.
_FIXED_OFFSET = 2


def get_function_object(obj):
    """
    Objects that wraps function should provide a "__numba__" magic attribute
    that contains a name of an attribute that contains the actual python
    function object.
    """
    attr = getattr(obj, "__numba__", None)
    if attr:
        return getattr(obj, attr)
    return obj


def get_code_object(obj):
    "Shamelessly borrowed from llpython"
    return getattr(obj, '__code__', getattr(obj, 'func_code', None))


def _as_opcodes(seq):
    lst = []
    for s in seq:
        c = dis.opmap.get(s)
        if c is not None:
            lst.append(c)
    return lst


JREL_OPS = frozenset(dis.hasjrel)
JABS_OPS = frozenset(dis.hasjabs)
JUMP_OPS = JREL_OPS | JABS_OPS
TERM_OPS = frozenset(_as_opcodes(['RETURN_VALUE', 'RAISE_VARARGS']))
EXTENDED_ARG = dis.EXTENDED_ARG
HAVE_ARGUMENT = dis.HAVE_ARGUMENT


class ByteCodeInst(object):
    '''
    Attributes
    ----------
    - offset:
        byte offset of opcode
    - opcode:
        opcode integer value
    - arg:
        instruction arg
    - lineno:
        -1 means unknown
    '''
    __slots__ = 'offset', 'next', 'opcode', 'opname', 'arg', 'lineno'

    def __init__(self, offset, opcode, arg, nextoffset):
        self.offset = offset
        self.next = nextoffset
        self.opcode = opcode
        self.opname = dis.opname[opcode]
        self.arg = arg
        self.lineno = -1  # unknown line number

    @property
    def is_jump(self):
        return self.opcode in JUMP_OPS

    @property
    def is_terminator(self):
        return self.opcode in TERM_OPS

    def get_jump_target(self):
        # With Python 3.10 the addressing of "bytecode" instructions has
        # changed from using bytes to using 16-bit words instead. As a
        # consequence the code to determine where a jump will lead had to be
        # adapted.
        # See also:
        # https://bugs.python.org/issue26647
        # https://bugs.python.org/issue27129
        # https://github.com/python/cpython/pull/25069
        assert self.is_jump
        if PYVERSION in ((3, 13),):
            if self.opcode in (dis.opmap[k]
                               for k in ["JUMP_BACKWARD",
                                         "JUMP_BACKWARD_NO_INTERRUPT"]):
                return self.next - (self.arg * 2)
        elif PYVERSION in ((3, 12),):
            if self.opcode in (dis.opmap[k]
                               for k in ["JUMP_BACKWARD"]):
                return self.offset - (self.arg - 1) * 2
        elif PYVERSION in ((3, 11), ):
            if self.opcode in (dis.opmap[k]
                               for k in ("JUMP_BACKWARD",
                                         "POP_JUMP_BACKWARD_IF_TRUE",
                                         "POP_JUMP_BACKWARD_IF_FALSE",
                                         "POP_JUMP_BACKWARD_IF_NONE",
                                         "POP_JUMP_BACKWARD_IF_NOT_NONE",)):
                return self.offset - (self.arg - 1) * 2
        elif PYVERSION in ((3, 10),):
            pass
        else:
            raise NotImplementedError(PYVERSION)

        if PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            if self.opcode in JREL_OPS:
                return self.next + self.arg * 2
            else:
                assert self.opcode in JABS_OPS
                return self.arg * 2 - 2
        else:
            raise NotImplementedError(PYVERSION)

    def __repr__(self):
        return '%s(arg=%s, lineno=%d)' % (self.opname, self.arg, self.lineno)

    @property
    def block_effect(self):
        """Effect of the block stack
        Returns +1 (push), 0 (none) or -1 (pop)
        """
        if self.opname.startswith('SETUP_'):
            return 1
        elif self.opname == 'POP_BLOCK':
            return -1
        else:
            return 0


CODE_LEN = 1
ARG_LEN = 1
NO_ARG_LEN = 1

OPCODE_NOP = dis.opname.index('NOP')


if PYVERSION in ((3, 13),):

    def _unpack_opargs(code):
        buf = []
        for i, start_offset, op, arg in dis._unpack_opargs(code):
            buf.append((start_offset, op, arg))
        for i, (start_offset, op, arg) in enumerate(buf):
            if i + 1 < len(buf):
                next_offset = buf[i + 1][0]
            else:
                next_offset = len(code)
            yield (start_offset, op, arg, next_offset)

elif PYVERSION in ((3, 10), (3, 11), (3, 12)):

    # Adapted from Lib/dis.py
    def _unpack_opargs(code):
        """
        Returns a 4-int-tuple of
        (bytecode offset, opcode, argument, offset of next bytecode).
        """
        extended_arg = 0
        n = len(code)
        offset = i = 0
        while i < n:
            op = code[i]
            i += CODE_LEN
            if op >= HAVE_ARGUMENT:
                arg = code[i] | extended_arg
                for j in range(ARG_LEN):
                    arg |= code[i + j] << (8 * j)
                i += ARG_LEN
                if PYVERSION in ((3, 12),):
                    # Python 3.12 introduced cache slots. We need to account for
                    # cache slots when we determine the offset of the next
                    # opcode. The number of cache slots is specific to each
                    # opcode and can be looked up in the _inline_cache_entries
                    # dictionary.
                    i += _inline_cache_entries[op] * INSTR_LEN
                elif PYVERSION in ((3, 10), (3, 11)):
                    pass
                else:
                    raise NotImplementedError(PYVERSION)
                if op == EXTENDED_ARG:
                    # This is a deviation from what dis does...
                    # In python 3.11 it seems like EXTENDED_ARGs appear more
                    # often and are also used as jump targets. So as to not have
                    # to do "book keeping" for where EXTENDED_ARGs have been
                    # "skipped" they are replaced with NOPs so as to provide a
                    # legal jump target and also ensure that the bytecode
                    # offsets are correct.
                    yield (offset, OPCODE_NOP, arg, i)
                    extended_arg = arg << 8 * ARG_LEN
                    offset = i
                    continue
            else:
                arg = None
                i += NO_ARG_LEN
                if PYVERSION in ((3, 12),):
                    # Python 3.12 introduced cache slots. We need to account for
                    # cache slots when we determine the offset of the next
                    # opcode. The number of cache slots is specific to each
                    # opcode and can be looked up in the _inline_cache_entries
                    # dictionary.
                    i += _inline_cache_entries[op] * INSTR_LEN
                elif PYVERSION in ((3, 10), (3, 11)):
                    pass
                else:
                    raise NotImplementedError(PYVERSION)

            extended_arg = 0
            yield (offset, op, arg, i)
            offset = i  # Mark inst offset at first extended
else:
    raise NotImplementedError(PYVERSION)


def _patched_opargs(bc_stream):
    """Patch the bytecode stream.

    - Adds a NOP bytecode at the start to avoid jump target being at the entry.
    """
    # Injected NOP
    yield (0, OPCODE_NOP, None, _FIXED_OFFSET)
    # Adjust bytecode offset for the rest of the stream
    for offset, opcode, arg, nextoffset in bc_stream:
        # If the opcode has an absolute jump target, adjust it.
        if opcode in JABS_OPS:
            arg += _FIXED_OFFSET
        yield offset + _FIXED_OFFSET, opcode, arg, nextoffset + _FIXED_OFFSET


class ByteCodeIter(object):
    def __init__(self, code):
        self.code = code
        self.iter = iter(_patched_opargs(_unpack_opargs(self.code.co_code)))

    def __iter__(self):
        return self

    def _fetch_opcode(self):
        return next(self.iter)

    def next(self):
        offset, opcode, arg, nextoffset = self._fetch_opcode()
        return offset, ByteCodeInst(offset=offset, opcode=opcode, arg=arg,
                                    nextoffset=nextoffset)

    __next__ = next

    def read_arg(self, size):
        buf = 0
        for i in range(size):
            _offset, byte = next(self.iter)
            buf |= byte << (8 * i)
        return buf


class _ByteCode(object):
    """
    The decoded bytecode of a function, and related information.
    """
    __slots__ = ('func_id', 'co_names', 'co_varnames', 'co_consts',
                 'co_cellvars', 'co_freevars', 'exception_entries',
                 'table', 'labels')

    def __init__(self, func_id):
        code = func_id.code

        labels = set(x + _FIXED_OFFSET for x in dis.findlabels(code.co_code))
        labels.add(0)

        # A map of {offset: ByteCodeInst}
        table = OrderedDict(ByteCodeIter(code))
        self._compute_lineno(table, code)

        self.func_id = func_id
        self.co_names = code.co_names
        self.co_varnames = code.co_varnames
        self.co_consts = code.co_consts
        self.co_cellvars = code.co_cellvars
        self.co_freevars = code.co_freevars

        self.table = table
        self.labels = sorted(labels)

    @classmethod
    def _compute_lineno(cls, table, code):
        """
        Compute the line numbers for all bytecode instructions.
        """
        for offset, lineno in dis.findlinestarts(code):
            adj_offset = offset + _FIXED_OFFSET
            if adj_offset in table:
                table[adj_offset].lineno = lineno
        # Assign unfilled lineno
        # Start with first bytecode's lineno
        known = code.co_firstlineno
        for inst in table.values():
            if inst.lineno is not None and inst.lineno >= 0:
                known = inst.lineno
            else:
                inst.lineno = known
        return table

    def __iter__(self):
        return iter(self.table.values())

    def __getitem__(self, offset):
        return self.table[offset]

    def __contains__(self, offset):
        return offset in self.table

    def dump(self):
        def label_marker(i):
            if i[1].offset in self.labels:
                return '>'
            else:
                return ' '

        return '\n'.join('%s %10s\t%s' % ((label_marker(i),) + i)
                         for i in self.table.items()
                         if i[1].opname != "CACHE")

    @classmethod
    def _compute_used_globals(cls, func, table, co_consts, co_names):
        """
        Compute the globals used by the function with the given
        bytecode table.
        """
        d = {}
        globs = func.__globals__
        builtins = globs.get('__builtins__', utils.builtins)
        if isinstance(builtins, ModuleType):
            builtins = builtins.__dict__
        # Look for LOAD_GLOBALs in the bytecode
        for inst in table.values():
            if inst.opname == 'LOAD_GLOBAL':
                name = co_names[_fix_LOAD_GLOBAL_arg(inst.arg)]
                if name not in d:
                    try:
                        value = globs[name]
                    except KeyError:
                        value = builtins[name]
                    d[name] = value
        # Add globals used by any nested code object
        for co in co_consts:
            if isinstance(co, CodeType):
                subtable = OrderedDict(ByteCodeIter(co))
                d.update(cls._compute_used_globals(func, subtable,
                                                   co.co_consts, co.co_names))
        return d

    def get_used_globals(self):
        """
        Get a {name: value} map of the globals used by this code
        object and any nested code objects.
        """
        return self._compute_used_globals(self.func_id.func, self.table,
                                          self.co_consts, self.co_names)


def _fix_LOAD_GLOBAL_arg(arg):
    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        return arg >> 1
    elif PYVERSION in ((3, 10),):
        return arg
    else:
        raise NotImplementedError(PYVERSION)


class ByteCodePy311(_ByteCode):

    def __init__(self, func_id):
        super().__init__(func_id)
        entries = dis.Bytecode(func_id.code).exception_entries
        self.exception_entries = tuple(map(self.fixup_eh, entries))

    @staticmethod
    def fixup_eh(ent):
        # Patch up the exception table offset
        # because we add a NOP in _patched_opargs
        out = dis._ExceptionTableEntry(
            start=ent.start + _FIXED_OFFSET, end=ent.end + _FIXED_OFFSET,
            target=ent.target + _FIXED_OFFSET,
            depth=ent.depth, lasti=ent.lasti,
        )
        return out

    def find_exception_entry(self, offset):
        """
        Returns the exception entry for the given instruction offset
        """
        candidates = []
        for ent in self.exception_entries:
            if ent.start <= offset < ent.end:
                candidates.append((ent.depth, ent))
        if candidates:
            ent = max(candidates)[1]
            return ent


class ByteCodePy312(ByteCodePy311):

    def __init__(self, func_id):
        super().__init__(func_id)

        # initialize lazy property
        self._ordered_offsets = None

        # Fixup offsets for all exception entries.
        entries = [self.fixup_eh(e) for e in
                   dis.Bytecode(func_id.code).exception_entries
                   ]

        # Remove exceptions, innermost ones first
        # Can be done by using a stack
        entries = self.remove_build_list_swap_pattern(entries)

        # If this is a generator, we need to skip any exception table entries
        # that point to the exception handler with the highest offset.
        if func_id.is_generator:
            # Get the exception handler with the highest offset.
            max_exception_target = max([e.target for e in entries])
            # Remove any exception table entries that point to that exception
            # handler.
            entries = [e for e in entries if e.target != max_exception_target]

        self.exception_entries = tuple(entries)

    @property
    def ordered_offsets(self):
        if not self._ordered_offsets:
            # Get an ordered list of offsets.
            self._ordered_offsets = [o for o in self.table]
        return self._ordered_offsets

    def remove_build_list_swap_pattern(self, entries):
        """ Find the following bytecode pattern:

            BUILD_{LIST, MAP, SET}
            SWAP(2)
            FOR_ITER
            ...
            END_FOR
            SWAP(2)

            This pattern indicates that a list/dict/set comprehension has
            been inlined. In this case we can skip the exception blocks
            entirely along with the dead exceptions that it points to.
            A pair of exception that sandwiches these exception will
            also be merged into a single exception.

            Update for Python 3.13, the ending of the pattern has a extra
            POP_TOP:

            ...
            END_FOR
            POP_TOP
            SWAP(2)

            Update for Python 3.13.1, there's now a GET_ITER before FOR_ITER.
            This patch the GET_ITER to NOP to minimize changes downstream
            (e.g. array-comprehension).
        """
        def pop_and_merge_exceptions(entries: list,
                                     entry_to_remove: _ExceptionTableEntry):
            lower_entry_idx = entries.index(entry_to_remove) - 1
            upper_entry_idx = entries.index(entry_to_remove) + 1

            # Merge the upper and lower exceptions if possible.
            if lower_entry_idx >= 0 and upper_entry_idx < len(entries):
                lower_entry = entries[lower_entry_idx]
                upper_entry = entries[upper_entry_idx]
                if lower_entry.target == upper_entry.target:
                    entries[lower_entry_idx] = _ExceptionTableEntry(
                        lower_entry.start,
                        upper_entry.end,
                        lower_entry.target,
                        lower_entry.depth,
                        upper_entry.lasti)
                    entries.remove(upper_entry)

            # Remove the exception entry.
            entries.remove(entry_to_remove)
            # Remove dead exceptions, if any, that the entry above may point to.
            entries = [e for e in entries
                       if not e.start == entry_to_remove.target]
            return entries

        change_to_nop = set()
        work_remaining = True
        while work_remaining:
            # Temporarily set work_remaining to False, if we find a pattern
            # then work is not complete, hence we set it again to True.
            work_remaining = False
            current_nop_fixes = set()
            for entry in entries.copy():
                # Check start of pattern, three instructions.
                # Work out the index of the instruction.
                index = self.ordered_offsets.index(entry.start)
                # If there is a BUILD_{LIST, MAP, SET} instruction at this
                # location.
                curr_inst = self.table[self.ordered_offsets[index]]
                if curr_inst.opname not in ("BUILD_LIST",
                                            "BUILD_MAP",
                                            "BUILD_SET"):
                    continue
                # Check if the BUILD_{LIST, MAP, SET} instruction is followed
                # by a SWAP(2).
                next_inst = self.table[self.ordered_offsets[index + 1]]
                if not next_inst.opname == "SWAP" and next_inst.arg == 2:
                    continue
                next_inst = self.table[self.ordered_offsets[index + 2]]
                # Check if the SWAP is followed by a FOR_ITER
                # BUT Python3.13.1 introduced an extra GET_ITER.
                # If we see a GET_ITER here, check if the next thing is a
                # FOR_ITER.
                if next_inst.opname == "GET_ITER":
                    # Add the inst to potentially be replaced to NOP
                    current_nop_fixes.add(next_inst)
                    # Loop up next instruction.
                    next_inst = self.table[self.ordered_offsets[index + 3]]

                if not next_inst.opname == "FOR_ITER":
                    continue

                if PYVERSION in ((3, 13),):
                    # Check end of pattern, two instructions.
                    # Check for the corresponding END_FOR, exception table end
                    # is non-inclusive, so subtract one.
                    index = self.ordered_offsets.index(entry.end)
                    curr_inst = self.table[self.ordered_offsets[index - 2]]
                    if not curr_inst.opname == "END_FOR":
                        continue
                    next_inst = self.table[self.ordered_offsets[index - 1]]
                    if not next_inst.opname == "POP_TOP":
                        continue
                    # END_FOR must be followed by SWAP(2)
                    next_inst = self.table[self.ordered_offsets[index]]
                    if not next_inst.opname == "SWAP" and next_inst.arg == 2:
                        continue
                elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
                    # Check end of pattern, two instructions.
                    # Check for the corresponding END_FOR, exception table end
                    # is non-inclusive, so subtract one.
                    index = self.ordered_offsets.index(entry.end)
                    curr_inst = self.table[self.ordered_offsets[index - 1]]
                    if not curr_inst.opname == "END_FOR":
                        continue
                    # END_FOR must be followed by SWAP(2)
                    next_inst = self.table[self.ordered_offsets[index]]
                    if not next_inst.opname == "SWAP" and next_inst.arg == 2:
                        continue
                else:
                    raise NotImplementedError(PYVERSION)
                # If all conditions are met that means this exception entry
                # is for a list/dict/set comprehension and can be removed.
                # Also if there exist exception entries above and below this
                # entry pointing to the same target. those can be merged into
                # a single bigger exception block.
                entries = pop_and_merge_exceptions(entries, entry)
                work_remaining = True

                # Commit NOP fixes since we confirmed the suspects belong to
                # a comprehension code.
                change_to_nop |= current_nop_fixes

        # Complete fixes to NOPs
        for inst in change_to_nop:
            self.table[inst.offset] = ByteCodeInst(inst.offset,
                                                   dis.opmap["NOP"],
                                                   None,
                                                   inst.next)
        return entries


if PYVERSION == (3, 11):
    ByteCode = ByteCodePy311
elif PYVERSION in ((3, 12), (3, 13),):
    ByteCode = ByteCodePy312
elif PYVERSION < (3, 11):
    ByteCode = _ByteCode
else:
    raise NotImplementedError(PYVERSION)


class FunctionIdentity(serialize.ReduceMixin):
    """
    A function's identity and metadata.

    Note this typically represents a function whose bytecode is
    being compiled, not necessarily the top-level user function
    (the two might be distinct).
    """
    _unique_ids = itertools.count(1)

    @classmethod
    def from_function(cls, pyfunc):
        """
        Create the FunctionIdentity of the given function.
        """
        func = get_function_object(pyfunc)
        code = get_code_object(func)
        pysig = utils.pysignature(func)
        if not code:
            raise errors.ByteCodeSupportError(
                "%s does not provide its bytecode" % func)

        try:
            func_qualname = func.__qualname__
        except AttributeError:
            func_qualname = func.__name__

        self = cls()
        self.func = func
        self.func_qualname = func_qualname
        self.func_name = func_qualname.split('.')[-1]
        self.code = code
        self.module = inspect.getmodule(func)
        self.modname = (utils._dynamic_modname
                        if self.module is None
                        else self.module.__name__)
        self.is_generator = inspect.isgeneratorfunction(func)
        self.pysig = pysig
        self.filename = code.co_filename
        self.firstlineno = code.co_firstlineno
        self.arg_count = len(pysig.parameters)
        self.arg_names = list(pysig.parameters)

        # Even the same function definition can be compiled into
        # several different function objects with distinct closure
        # variables, so we make sure to disambiguate using an unique id.
        uid = next(cls._unique_ids)
        self.unique_name = '{}${}'.format(self.func_qualname, uid)
        self.unique_id = uid

        return self

    def derive(self):
        """Copy the object and increment the unique counter.
        """
        return self.from_function(self.func)

    def _reduce_states(self):
        """
        NOTE: part of ReduceMixin protocol
        """
        return dict(pyfunc=self.func)

    @classmethod
    def _rebuild(cls, pyfunc):
        """
        NOTE: part of ReduceMixin protocol
        """
        return cls.from_function(pyfunc)
