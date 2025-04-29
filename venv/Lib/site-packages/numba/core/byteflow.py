"""
Implement python 3.8+ bytecode analysis
"""
import dis
import logging
from collections import namedtuple, defaultdict, deque
from functools import total_ordering

from numba.core.utils import (UniqueDict, PYVERSION, ALL_BINOPS_TO_OPERATORS,
                              _lazy_pformat)
from numba.core.controlflow import NEW_BLOCKERS, CFGraph
from numba.core.ir import Loc
from numba.core.errors import UnsupportedBytecodeError


_logger = logging.getLogger(__name__)

_EXCEPT_STACK_OFFSET = 6
_FINALLY_POP = _EXCEPT_STACK_OFFSET
_NO_RAISE_OPS = frozenset({
    'LOAD_CONST',
    'NOP',
    'LOAD_DEREF',
    'PRECALL',
})

if PYVERSION in ((3, 12), (3, 13)):
    from enum import Enum

    # Operands for CALL_INTRINSIC_1
    class CALL_INTRINSIC_1_Operand(Enum):
        INTRINSIC_STOPITERATION_ERROR = 3
        UNARY_POSITIVE = 5
        INTRINSIC_LIST_TO_TUPLE = 6
    ci1op = CALL_INTRINSIC_1_Operand
elif PYVERSION in ((3, 10), (3, 11)):
    pass
else:
    raise NotImplementedError(PYVERSION)


@total_ordering
class BlockKind(object):
    """Kinds of block to make related code safer than just `str`.
    """
    _members = frozenset({
        'LOOP',
        'TRY', 'EXCEPT', 'FINALLY',
        'WITH', 'WITH_FINALLY',
    })

    def __init__(self, value):
        assert value in self._members
        self._value = value

    def __hash__(self):
        return hash((type(self), self._value))

    def __lt__(self, other):
        if isinstance(other, BlockKind):
            return self._value < other._value
        else:
            raise TypeError('cannot compare to {!r}'.format(type(other)))

    def __eq__(self, other):
        if isinstance(other, BlockKind):
            return self._value == other._value
        else:
            raise TypeError('cannot compare to {!r}'.format(type(other)))

    def __repr__(self):
        return "BlockKind({})".format(self._value)


class Flow(object):
    """Data+Control Flow analysis.

    Simulate execution to recover dataflow and controlflow information.
    """
    def __init__(self, bytecode):
        _logger.debug("bytecode dump:\n%s", bytecode.dump())
        self._bytecode = bytecode
        self.block_infos = UniqueDict()

    def run(self):
        """Run a trace over the bytecode over all reachable path.

        The trace starts at bytecode offset 0 and gathers stack and control-
        flow information by partially interpreting each bytecode.
        Each ``State`` instance in the trace corresponds to a basic-block.
        The State instances forks when a jump instruction is encountered.
        A newly forked state is then added to the list of pending states.
        The trace ends when there are no more pending states.
        """
        firststate = State(bytecode=self._bytecode, pc=0, nstack=0,
                           blockstack=())
        runner = TraceRunner(debug_filename=self._bytecode.func_id.filename)
        runner.pending.append(firststate)

        # Enforce unique-ness on initial PC to avoid re-entering the PC with
        # a different stack-depth. We don't know if such a case is ever
        # possible, but no such case has been encountered in our tests.
        first_encounter = UniqueDict()
        # Loop over each pending state at a initial PC.
        # Each state is tracing a basic block
        while runner.pending:
            _logger.debug("pending: %s", runner.pending)
            state = runner.pending.popleft()
            if state not in runner.finished:
                _logger.debug("stack: %s", state._stack)
                _logger.debug("state.pc_initial: %s", state)
                first_encounter[state.pc_initial] = state
                # Loop over the state until it is terminated.
                while True:
                    runner.dispatch(state)
                    # Terminated?
                    if state.has_terminated():
                        break
                    else:
                        if self._run_handle_exception(runner, state):
                            break

                        if self._is_implicit_new_block(state):
                            # check if this is a with...as, abort if so
                            self._guard_with_as(state)
                            # else split
                            state.split_new_block()
                            break
                _logger.debug("end state. edges=%s", state.outgoing_edges)
                runner.finished.add(state)
                out_states = state.get_outgoing_states()
                runner.pending.extend(out_states)

        # Complete controlflow
        self._build_cfg(runner.finished)
        # Prune redundant PHI-nodes
        self._prune_phis(runner)
        # Post process
        for state in sorted(runner.finished, key=lambda x: x.pc_initial):
            self.block_infos[state.pc_initial] = si = adapt_state_infos(state)
            _logger.debug("block_infos %s:\n%s", state, si)

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def _run_handle_exception(self, runner, state):
            if not state.in_with() and (
                    state.has_active_try() and
                    state.get_inst().opname not in _NO_RAISE_OPS):
                # Is in a *try* block
                state.fork(pc=state.get_inst().next)
                runner._adjust_except_stack(state)
                return True
            else:
                state.advance_pc()

                # Must the new PC be a new block?
                if not state.in_with() and state.is_in_exception():
                    _logger.debug("3.11 exception %s PC=%s",
                                  state.get_exception(), state._pc)
                    eh = state.get_exception()
                    eh_top = state.get_top_block('TRY')
                    if eh_top and eh_top['end'] == eh.target:
                        # Same exception
                        eh_block = None
                    else:
                        eh_block = state.make_block("TRY", end=eh.target)
                        eh_block['end_offset'] = eh.end
                        eh_block['stack_depth'] = eh.depth
                        eh_block['push_lasti'] = eh.lasti
                        state.fork(pc=state._pc, extra_block=eh_block)
                        return True
    elif PYVERSION in ((3, 10),):
        def _run_handle_exception(self, runner, state):
            if (state.has_active_try() and
                    state.get_inst().opname not in _NO_RAISE_OPS):
                # Is in a *try* block
                state.fork(pc=state.get_inst().next)
                tryblk = state.get_top_block('TRY')
                state.pop_block_and_above(tryblk)
                nstack = state.stack_depth
                kwargs = {}
                if nstack > tryblk['entry_stack']:
                    kwargs['npop'] = nstack - tryblk['entry_stack']
                handler = tryblk['handler']
                kwargs['npush'] = {
                    BlockKind('EXCEPT'): _EXCEPT_STACK_OFFSET,
                    BlockKind('FINALLY'): _FINALLY_POP
                }[handler['kind']]
                kwargs['extra_block'] = handler
                state.fork(pc=tryblk['end'], **kwargs)
                return True
            else:
                state.advance_pc()
    else:
        raise NotImplementedError(PYVERSION)

    def _build_cfg(self, all_states):
        graph = CFGraph()
        for state in all_states:
            b = state.pc_initial
            graph.add_node(b)
        for state in all_states:
            for edge in state.outgoing_edges:
                graph.add_edge(state.pc_initial, edge.pc, 0)
        graph.set_entry_point(0)
        graph.process()
        self.cfgraph = graph

    def _prune_phis(self, runner):
        # Find phis that are unused in the local block
        _logger.debug("Prune PHIs".center(60, '-'))

        # Compute dataflow for used phis and propagate

        # 1. Get used-phis for each block
        # Map block to used_phis
        def get_used_phis_per_state():
            used_phis = defaultdict(set)
            phi_set = set()
            for state in runner.finished:
                used = set(state._used_regs)
                phis = set(state._phis)
                used_phis[state] |= phis & used
                phi_set |= phis
            return used_phis, phi_set

        # Find use-defs
        def find_use_defs():
            defmap = {}
            phismap = defaultdict(set)
            for state in runner.finished:
                for phi, rhs in state._outgoing_phis.items():
                    if rhs not in phi_set:
                        # Is a definition
                        defmap[phi] = state
                    phismap[phi].add((rhs, state))
            _logger.debug("defmap: %s", _lazy_pformat(defmap))
            _logger.debug("phismap: %s", _lazy_pformat(phismap))
            return defmap, phismap

        def propagate_phi_map(phismap):
            """An iterative dataflow algorithm to find the definition
            (the source) of each PHI node.
            """
            blacklist = defaultdict(set)

            while True:
                changing = False
                for phi, defsites in sorted(list(phismap.items())):
                    for rhs, state in sorted(list(defsites)):
                        if rhs in phi_set:
                            defsites |= phismap[rhs]
                            blacklist[phi].add((rhs, state))
                    to_remove = blacklist[phi]
                    if to_remove & defsites:
                        defsites -= to_remove
                        changing = True

                _logger.debug("changing phismap: %s", _lazy_pformat(phismap))
                if not changing:
                    break

        def apply_changes(used_phis, phismap):
            keep = {}
            for state, used_set in used_phis.items():
                for phi in used_set:
                    keep[phi] = phismap[phi]
            _logger.debug("keep phismap: %s", _lazy_pformat(keep))
            new_out = defaultdict(dict)
            for phi in keep:
                for rhs, state in keep[phi]:
                    new_out[state][phi] = rhs

            _logger.debug("new_out: %s", _lazy_pformat(new_out))
            for state in runner.finished:
                state._outgoing_phis.clear()
                state._outgoing_phis.update(new_out[state])

        used_phis, phi_set = get_used_phis_per_state()
        _logger.debug("Used_phis: %s", _lazy_pformat(used_phis))
        defmap, phismap = find_use_defs()
        propagate_phi_map(phismap)
        apply_changes(used_phis, phismap)
        _logger.debug("DONE Prune PHIs".center(60, '-'))

    def _is_implicit_new_block(self, state):
        inst = state.get_inst()

        if inst.offset in self._bytecode.labels:
            return True
        elif inst.opname in NEW_BLOCKERS:
            return True
        else:
            return False

    def _guard_with_as(self, state):
        """Checks if the next instruction after a SETUP_WITH is something other
        than a POP_TOP, if it is something else it'll be some sort of store
        which is not supported (this corresponds to `with CTXMGR as VAR(S)`)."""
        current_inst = state.get_inst()
        if current_inst.opname in {"SETUP_WITH", "BEFORE_WITH"}:
            next_op = self._bytecode[current_inst.next].opname
            if next_op != "POP_TOP":
                msg = ("The 'with (context manager) as "
                       "(variable):' construct is not "
                       "supported.")
                raise UnsupportedBytecodeError(msg)


def _is_null_temp_reg(reg):
    return reg.startswith("$null$")


class TraceRunner(object):
    """Trace runner contains the states for the trace and the opcode dispatch.
    """
    def __init__(self, debug_filename):
        self.debug_filename = debug_filename
        self.pending = deque()
        self.finished = set()

    def get_debug_loc(self, lineno):
        return Loc(self.debug_filename, lineno)

    def dispatch(self, state):
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            if state._blockstack:
                state: State
                while state._blockstack:
                    topblk = state._blockstack[-1]
                    blk_end = topblk['end']
                    if blk_end is not None and blk_end <= state.pc_initial:
                        state._blockstack.pop()
                    else:
                        break
        elif PYVERSION in ((3, 10),):
            pass
        else:
            raise NotImplementedError(PYVERSION)
        inst = state.get_inst()
        if inst.opname != "CACHE":
            _logger.debug("dispatch pc=%s, inst=%s", state._pc, inst)
            _logger.debug("stack %s", state._stack)
        fn = getattr(self, "op_{}".format(inst.opname), None)
        if fn is not None:
            fn(state, inst)
        else:
            msg = "Use of unsupported opcode (%s) found" % inst.opname
            raise UnsupportedBytecodeError(msg,
                                           loc=self.get_debug_loc(inst.lineno))

    def _adjust_except_stack(self, state):
        """
        Adjust stack when entering an exception handler to match expectation
        by the bytecode.
        """
        tryblk = state.get_top_block('TRY')
        state.pop_block_and_above(tryblk)
        nstack = state.stack_depth
        kwargs = {}
        expected_depth = tryblk['stack_depth']
        if nstack > expected_depth:
            # Pop extra item in the stack
            kwargs['npop'] = nstack - expected_depth
        # Set extra stack itemcount due to the exception values.
        extra_stack = 1
        if tryblk['push_lasti']:
            extra_stack += 1
        kwargs['npush'] = extra_stack
        state.fork(pc=tryblk['end'], **kwargs)

    def op_NOP(self, state, inst):
        state.append(inst)

    def op_RESUME(self, state, inst):
        state.append(inst)

    def op_CACHE(self, state, inst):
        state.append(inst)

    def op_PRECALL(self, state, inst):
        state.append(inst)

    def op_PUSH_NULL(self, state, inst):
        state.push(state.make_null())
        state.append(inst)

    def op_RETURN_GENERATOR(self, state, inst):
        # This impl doesn't follow what CPython does. CPython is hacking
        # the frame stack in the interpreter. From usage, it always
        # has a POP_TOP after it so we push a dummy value to the stack.
        #
        # Example bytecode:
        # >          0	NOP(arg=None, lineno=80)
        #            2	RETURN_GENERATOR(arg=None, lineno=80)
        #            4	POP_TOP(arg=None, lineno=80)
        #            6	RESUME(arg=0, lineno=80)
        state.push(state.make_temp())
        state.append(inst)

    if PYVERSION in ((3, 13),):
        def op_FORMAT_SIMPLE(self, state, inst):
            assert PYVERSION == (3, 13)
            value = state.pop()
            strvar = state.make_temp()
            res = state.make_temp()
            state.append(inst, value=value, res=res, strvar=strvar)
            state.push(res)

    def op_FORMAT_VALUE(self, state, inst):
        """
        FORMAT_VALUE(flags): flags argument specifies format spec which is
        not supported yet. Currently, we just call str() on the value.
        Pops a value from stack and pushes results back.
        Required for supporting f-strings.
        https://docs.python.org/3/library/dis.html#opcode-FORMAT_VALUE
        """
        if inst.arg != 0:
            msg = "format spec in f-strings not supported yet"
            raise UnsupportedBytecodeError(msg,
                                           loc=self.get_debug_loc(inst.lineno))
        value = state.pop()
        strvar = state.make_temp()
        res = state.make_temp()
        state.append(inst, value=value, res=res, strvar=strvar)
        state.push(res)

    def op_BUILD_STRING(self, state, inst):
        """
        BUILD_STRING(count): Concatenates count strings from the stack and
        pushes the resulting string onto the stack.
        Required for supporting f-strings.
        https://docs.python.org/3/library/dis.html#opcode-BUILD_STRING
        """
        count = inst.arg
        strings = list(reversed([state.pop() for _ in range(count)]))
        # corner case: f""
        if count == 0:
            tmps = [state.make_temp()]
        else:
            tmps = [state.make_temp() for _ in range(count - 1)]
        state.append(inst, strings=strings, tmps=tmps)
        state.push(tmps[-1])

    def op_POP_TOP(self, state, inst):
        state.pop()

    if PYVERSION in ((3, 13),):
        def op_TO_BOOL(self, state, inst):
            res = state.make_temp()
            tos = state.pop()
            state.append(inst, val=tos, res=res)
            state.push(res)

    elif PYVERSION < (3, 13):
        pass

    if PYVERSION in ((3, 13),):
        def op_LOAD_GLOBAL(self, state, inst):
            # Ordering of the global value and NULL is swapped in Py3.13
            res = state.make_temp()
            idx = inst.arg >> 1
            state.append(inst, idx=idx, res=res)
            state.push(res)
            # ignoring the NULL
            if inst.arg & 1:
                state.push(state.make_null())
    elif PYVERSION in ((3, 11), (3, 12)):
        def op_LOAD_GLOBAL(self, state, inst):
            res = state.make_temp()
            idx = inst.arg >> 1
            state.append(inst, idx=idx, res=res)
            # ignoring the NULL
            if inst.arg & 1:
                state.push(state.make_null())
            state.push(res)
    elif PYVERSION in ((3, 10),):
        def op_LOAD_GLOBAL(self, state, inst):
            res = state.make_temp()
            state.append(inst, res=res)
            state.push(res)
    else:
        raise NotImplementedError(PYVERSION)

    def op_COPY_FREE_VARS(self, state, inst):
        state.append(inst)

    def op_MAKE_CELL(self, state, inst):
        state.append(inst)

    def op_LOAD_DEREF(self, state, inst):
        res = state.make_temp()
        state.append(inst, res=res)
        state.push(res)

    def op_LOAD_CONST(self, state, inst):
        # append const index for interpreter to read the const value
        res = state.make_temp("const") + f".{inst.arg}"
        state.push(res)
        state.append(inst, res=res)

    def op_LOAD_ATTR(self, state, inst):
        item = state.pop()
        res = state.make_temp()
        if PYVERSION in ((3, 13),):
            state.push(res)  # the attr
            if inst.arg & 1:
                state.push(state.make_null())
        elif PYVERSION in ((3, 12),):
            if inst.arg & 1:
                state.push(state.make_null())
            state.push(res)
        elif PYVERSION in ((3, 10), (3, 11)):
            state.push(res)
        else:
            raise NotImplementedError(PYVERSION)
        state.append(inst, item=item, res=res)

    def op_LOAD_FAST(self, state, inst):
        assert PYVERSION <= (3, 13)
        if PYVERSION in ((3, 13), ):
            try:
                name = state.get_varname(inst)
            except IndexError:   # oparg is out of range
                # Handle this like a LOAD_DEREF
                # Assume MAKE_CELL and COPY_FREE_VARS has correctly setup the
                # states.
                # According to https://github.com/python/cpython/blob/9ac606080a0074cdf7589d9b7c9413a73e0ddf37/Objects/codeobject.c#L730C9-L759 # noqa E501
                # localsplus is locals + cells + freevars
                bc = state._bytecode
                num_varnames = len(bc.co_varnames)
                num_freevars = len(bc.co_freevars)
                num_cellvars = len(bc.co_cellvars)
                max_fast_local = num_cellvars + num_freevars
                assert 0 <= inst.arg - num_varnames < max_fast_local
                res = state.make_temp()
                state.append(inst, res=res, as_load_deref=True)
                state.push(res)
                return
        else:
            name = state.get_varname(inst)
        res = state.make_temp(name)
        state.append(inst, res=res)
        state.push(res)

    if PYVERSION in ((3, 13),):
        def op_LOAD_FAST_LOAD_FAST(self, state, inst):
            oparg = inst.arg
            oparg1 = oparg >> 4
            oparg2 = oparg & 15
            name1 = state.get_varname_by_arg(oparg1)
            name2 = state.get_varname_by_arg(oparg2)
            res1 = state.make_temp(name1)
            res2 = state.make_temp(name2)
            state.append(inst, res1=res1, res2=res2)
            state.push(res1)
            state.push(res2)

        def op_STORE_FAST_LOAD_FAST(self, state, inst):
            oparg = inst.arg
            # oparg1 = oparg >> 4  # not needed
            oparg2 = oparg & 15
            store_value = state.pop()
            load_name = state.get_varname_by_arg(oparg2)
            load_res = state.make_temp(load_name)
            state.append(inst, store_value=store_value, load_res=load_res)
            state.push(load_res)

        def op_STORE_FAST_STORE_FAST(self, state, inst):
            value1 = state.pop()
            value2 = state.pop()
            state.append(inst, value1=value1, value2=value2)

    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 12), (3, 13)):
        op_LOAD_FAST_CHECK = op_LOAD_FAST
        op_LOAD_FAST_AND_CLEAR = op_LOAD_FAST
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_DELETE_FAST(self, state, inst):
        state.append(inst)

    def op_DELETE_ATTR(self, state, inst):
        target = state.pop()
        state.append(inst, target=target)

    def op_STORE_ATTR(self, state, inst):
        target = state.pop()
        value = state.pop()
        state.append(inst, target=target, value=value)

    def op_STORE_DEREF(self, state, inst):
        value = state.pop()
        state.append(inst, value=value)

    def op_STORE_FAST(self, state, inst):
        value = state.pop()
        state.append(inst, value=value)

    def op_SLICE_1(self, state, inst):
        """
        TOS = TOS1[TOS:]
        """
        tos = state.pop()
        tos1 = state.pop()
        res = state.make_temp()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            start=tos,
            res=res,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )
        state.push(res)

    def op_SLICE_2(self, state, inst):
        """
        TOS = TOS1[:TOS]
        """
        tos = state.pop()
        tos1 = state.pop()
        res = state.make_temp()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            stop=tos,
            res=res,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )
        state.push(res)

    def op_SLICE_3(self, state, inst):
        """
        TOS = TOS2[TOS1:TOS]
        """
        tos = state.pop()
        tos1 = state.pop()
        tos2 = state.pop()
        res = state.make_temp()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        state.append(
            inst,
            base=tos2,
            start=tos1,
            stop=tos,
            res=res,
            slicevar=slicevar,
            indexvar=indexvar,
        )
        state.push(res)

    def op_STORE_SLICE_0(self, state, inst):
        """
        TOS[:] = TOS1
        """
        tos = state.pop()
        value = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos,
            value=value,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_STORE_SLICE_1(self, state, inst):
        """
        TOS1[TOS:] = TOS2
        """
        tos = state.pop()
        tos1 = state.pop()
        value = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            start=tos,
            slicevar=slicevar,
            value=value,
            indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_STORE_SLICE_2(self, state, inst):
        """
        TOS1[:TOS] = TOS2
        """
        tos = state.pop()
        tos1 = state.pop()
        value = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            stop=tos,
            value=value,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_STORE_SLICE_3(self, state, inst):
        """
        TOS2[TOS1:TOS] = TOS3
        """
        tos = state.pop()
        tos1 = state.pop()
        tos2 = state.pop()
        value = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        state.append(
            inst,
            base=tos2,
            start=tos1,
            stop=tos,
            value=value,
            slicevar=slicevar,
            indexvar=indexvar,
        )

    def op_DELETE_SLICE_0(self, state, inst):
        """
        del TOS[:]
        """
        tos = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst, base=tos, slicevar=slicevar, indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_DELETE_SLICE_1(self, state, inst):
        """
        del TOS1[TOS:]
        """
        tos = state.pop()
        tos1 = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            start=tos,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_DELETE_SLICE_2(self, state, inst):
        """
        del TOS1[:TOS]
        """
        tos = state.pop()
        tos1 = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        nonevar = state.make_temp()
        state.append(
            inst,
            base=tos1,
            stop=tos,
            slicevar=slicevar,
            indexvar=indexvar,
            nonevar=nonevar,
        )

    def op_DELETE_SLICE_3(self, state, inst):
        """
        del TOS2[TOS1:TOS]
        """
        tos = state.pop()
        tos1 = state.pop()
        tos2 = state.pop()
        slicevar = state.make_temp()
        indexvar = state.make_temp()
        state.append(
            inst, base=tos2, start=tos1, stop=tos, slicevar=slicevar,
            indexvar=indexvar
        )

    def op_BUILD_SLICE(self, state, inst):
        """
        slice(TOS1, TOS) or slice(TOS2, TOS1, TOS)
        """
        argc = inst.arg
        if argc == 2:
            tos = state.pop()
            tos1 = state.pop()
            start = tos1
            stop = tos
            step = None
        elif argc == 3:
            tos = state.pop()
            tos1 = state.pop()
            tos2 = state.pop()
            start = tos2
            stop = tos1
            step = tos
        else:
            raise Exception("unreachable")
        slicevar = state.make_temp()
        res = state.make_temp()
        state.append(
            inst, start=start, stop=stop, step=step, res=res, slicevar=slicevar
        )
        state.push(res)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_BINARY_SLICE(self, state, inst):
            end = state.pop()
            start = state.pop()
            container = state.pop()
            temp_res = state.make_temp()
            res = state.make_temp()
            slicevar = state.make_temp()
            state.append(
                inst, start=start, end=end, container=container, res=res,
                slicevar=slicevar, temp_res=temp_res
            )
            state.push(res)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_STORE_SLICE(self, state, inst):
            end = state.pop()
            start = state.pop()
            container = state.pop()
            value = state.pop()

            slicevar = state.make_temp()
            res = state.make_temp()
            state.append(
                inst, start=start, end=end, container=container, value=value,
                res=res, slicevar=slicevar,
            )
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def _op_POP_JUMP_IF(self, state, inst):
        pred = state.pop()
        state.append(inst, pred=pred)

        target_inst = inst.get_jump_target()
        next_inst = inst.next
        # if the next inst and the jump target are the same location, issue one
        # fork else issue a fork for the next and the target.
        state.fork(pc=next_inst)
        if target_inst != next_inst:
            state.fork(pc=target_inst)

    op_POP_JUMP_IF_TRUE = _op_POP_JUMP_IF
    op_POP_JUMP_IF_FALSE = _op_POP_JUMP_IF

    if PYVERSION in ((3, 12), (3, 13)):
        op_POP_JUMP_IF_NONE = _op_POP_JUMP_IF
        op_POP_JUMP_IF_NOT_NONE = _op_POP_JUMP_IF
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def _op_JUMP_IF_OR_POP(self, state, inst):
        pred = state.get_tos()
        state.append(inst, pred=pred)
        state.fork(pc=inst.next, npop=1)
        state.fork(pc=inst.get_jump_target())

    op_JUMP_IF_FALSE_OR_POP = _op_JUMP_IF_OR_POP
    op_JUMP_IF_TRUE_OR_POP = _op_JUMP_IF_OR_POP

    def op_POP_JUMP_FORWARD_IF_NONE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_FORWARD_IF_NOT_NONE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_BACKWARD_IF_NONE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_BACKWARD_IF_NOT_NONE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_FORWARD_IF_FALSE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_FORWARD_IF_TRUE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_BACKWARD_IF_FALSE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_POP_JUMP_BACKWARD_IF_TRUE(self, state, inst):
        self._op_POP_JUMP_IF(state, inst)

    def op_JUMP_FORWARD(self, state, inst):
        state.append(inst)
        state.fork(pc=inst.get_jump_target())

    def op_JUMP_BACKWARD(self, state, inst):
        state.append(inst)
        state.fork(pc=inst.get_jump_target())

    op_JUMP_BACKWARD_NO_INTERRUPT = op_JUMP_BACKWARD

    def op_JUMP_ABSOLUTE(self, state, inst):
        state.append(inst)
        state.fork(pc=inst.get_jump_target())

    def op_BREAK_LOOP(self, state, inst):
        # NOTE: bytecode removed since py3.8
        end = state.get_top_block('LOOP')['end']
        state.append(inst, end=end)
        state.pop_block()
        state.fork(pc=end)

    def op_RETURN_VALUE(self, state, inst):
        state.append(inst, retval=state.pop(), castval=state.make_temp())
        state.terminate()

    if PYVERSION in ((3, 12), (3, 13)):
        def op_RETURN_CONST(self, state, inst):
            res = state.make_temp("const")
            state.append(inst, retval=res, castval=state.make_temp())
            state.terminate()
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_YIELD_VALUE(self, state, inst):
        val = state.pop()
        res = state.make_temp()
        state.append(inst, value=val, res=res)
        state.push(res)

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_RAISE_VARARGS(self, state, inst):
            if inst.arg == 0:
                exc = None
                # No re-raising within a try-except block.
                # But we allow bare reraise.
                if state.has_active_try():
                    raise UnsupportedBytecodeError(
                        "The re-raising of an exception is not yet supported.",
                        loc=self.get_debug_loc(inst.lineno),
                    )
            elif inst.arg == 1:
                exc = state.pop()
            else:
                raise ValueError("Multiple argument raise is not supported.")
            state.append(inst, exc=exc)

            if state.has_active_try():
                self._adjust_except_stack(state)
            else:
                state.terminate()

    elif PYVERSION in ((3, 10),):
        def op_RAISE_VARARGS(self, state, inst):
            in_exc_block = any([
                state.get_top_block("EXCEPT") is not None,
                state.get_top_block("FINALLY") is not None
            ])
            if inst.arg == 0:
                exc = None
                if in_exc_block:
                    raise UnsupportedBytecodeError(
                        "The re-raising of an exception is not yet supported.",
                        loc=self.get_debug_loc(inst.lineno),
                    )
            elif inst.arg == 1:
                exc = state.pop()
            else:
                raise ValueError("Multiple argument raise is not supported.")
            state.append(inst, exc=exc)
            state.terminate()
    else:
        raise NotImplementedError(PYVERSION)

    def op_BEGIN_FINALLY(self, state, inst):
        temps = []
        for i in range(_EXCEPT_STACK_OFFSET):
            tmp = state.make_temp()
            temps.append(tmp)
            state.push(tmp)
        state.append(inst, temps=temps)

    def op_END_FINALLY(self, state, inst):
        blk = state.pop_block()
        state.reset_stack(blk['entry_stack'])

    if PYVERSION in ((3, 13),):
        def op_END_FOR(self, state, inst):
            state.pop()
    elif PYVERSION in ((3, 12),):
        def op_END_FOR(self, state, inst):
            state.pop()
            state.pop()
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_POP_FINALLY(self, state, inst):
        # we don't emulate the exact stack behavior
        if inst.arg != 0:
            msg = ('Unsupported use of a bytecode related to try..finally'
                   ' or a with-context')
            raise UnsupportedBytecodeError(msg,
                                           loc=self.get_debug_loc(inst.lineno))

    def op_CALL_FINALLY(self, state, inst):
        pass

    def op_WITH_EXCEPT_START(self, state, inst):
        state.terminate()  # do not support

    def op_WITH_CLEANUP_START(self, state, inst):
        # we don't emulate the exact stack behavior
        state.append(inst)

    def op_WITH_CLEANUP_FINISH(self, state, inst):
        # we don't emulate the exact stack behavior
        state.append(inst)

    def op_SETUP_LOOP(self, state, inst):
        # NOTE: bytecode removed since py3.8
        state.push_block(
            state.make_block(
                kind='LOOP',
                end=inst.get_jump_target(),
            )
        )

    def op_BEFORE_WITH(self, state, inst):
        # Almost the same as py3.10 SETUP_WITH just lacking the finally block.
        cm = state.pop()    # the context-manager

        yielded = state.make_temp()
        exitfn = state.make_temp(prefix='setup_with_exitfn')

        state.push(exitfn)
        state.push(yielded)

        # Gather all exception entries for this WITH. There maybe multiple
        # entries; esp. for nested WITHs.
        bc = state._bytecode
        ehhead = bc.find_exception_entry(inst.next)
        ehrelated = [ehhead]
        for eh in bc.exception_entries:
            if eh.target == ehhead.target:
                ehrelated.append(eh)
        end = max(eh.end for eh in ehrelated)
        state.append(inst, contextmanager=cm, exitfn=exitfn, end=end)

        state.push_block(
            state.make_block(
                kind='WITH',
                end=end,
            )
        )
        # Forces a new block
        state.fork(pc=inst.next)

    def op_SETUP_WITH(self, state, inst):
        cm = state.pop()    # the context-manager

        yielded = state.make_temp()
        exitfn = state.make_temp(prefix='setup_with_exitfn')
        state.append(inst, contextmanager=cm, exitfn=exitfn)

        state.push(exitfn)
        state.push(yielded)

        state.push_block(
            state.make_block(
                kind='WITH',
                end=inst.get_jump_target(),
            )
        )
        # Forces a new block
        state.fork(pc=inst.next)

    def _setup_try(self, kind, state, next, end):
        # Forces a new block
        # Fork to the body of the finally
        handler_block = state.make_block(
            kind=kind,
            end=None,
            reset_stack=False,
        )
        # Forces a new block
        # Fork to the body of the finally
        state.fork(
            pc=next,
            extra_block=state.make_block(
                kind='TRY',
                end=end,
                reset_stack=False,
                handler=handler_block,
            )
        )

    def op_PUSH_EXC_INFO(self, state, inst):
        tos = state.pop()
        state.push(state.make_temp("exception"))
        state.push(tos)

    def op_SETUP_FINALLY(self, state, inst):
        state.append(inst)
        self._setup_try(
            'FINALLY', state, next=inst.next, end=inst.get_jump_target(),
        )

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_POP_EXCEPT(self, state, inst):
            state.pop()

    elif PYVERSION in ((3, 10),):
        def op_POP_EXCEPT(self, state, inst):
            blk = state.pop_block()
            if blk['kind'] not in {BlockKind('EXCEPT'), BlockKind('FINALLY')}:
                raise UnsupportedBytecodeError(
                    f"POP_EXCEPT got an unexpected block: {blk['kind']}",
                    loc=self.get_debug_loc(inst.lineno),
                )
            state.pop()
            state.pop()
            state.pop()
            # Forces a new block
            state.fork(pc=inst.next)
    else:
        raise NotImplementedError(PYVERSION)

    def op_POP_BLOCK(self, state, inst):
        blk = state.pop_block()
        if blk['kind'] == BlockKind('TRY'):
            state.append(inst, kind='try')
        elif blk['kind'] == BlockKind('WITH'):
            state.append(inst, kind='with')
        state.fork(pc=inst.next)

    def op_BINARY_SUBSCR(self, state, inst):
        index = state.pop()
        target = state.pop()
        res = state.make_temp()
        state.append(inst, index=index, target=target, res=res)
        state.push(res)

    def op_STORE_SUBSCR(self, state, inst):
        index = state.pop()
        target = state.pop()
        value = state.pop()
        state.append(inst, target=target, index=index, value=value)

    def op_DELETE_SUBSCR(self, state, inst):
        index = state.pop()
        target = state.pop()
        state.append(inst, target=target, index=index)

    def op_CALL(self, state, inst):
        narg = inst.arg
        args = list(reversed([state.pop() for _ in range(narg)]))
        if PYVERSION == (3, 13):
            null_or_self = state.pop()
            # position of the callable is fixed
            callable = state.pop()
            if not _is_null_temp_reg(null_or_self):
                args = [null_or_self, *args]
            kw_names = None
        elif PYVERSION < (3, 13):
            callable_or_firstarg = state.pop()
            null_or_callable = state.pop()
            if _is_null_temp_reg(null_or_callable):
                callable = callable_or_firstarg
            else:
                callable = null_or_callable
                args = [callable_or_firstarg, *args]
            kw_names = state.pop_kw_names()
        res = state.make_temp()

        state.append(inst, func=callable, args=args, kw_names=kw_names, res=res)
        state.push(res)

    def op_KW_NAMES(self, state, inst):
        state.set_kw_names(inst.arg)

    def op_CALL_FUNCTION(self, state, inst):
        narg = inst.arg
        args = list(reversed([state.pop() for _ in range(narg)]))
        func = state.pop()

        res = state.make_temp()
        state.append(inst, func=func, args=args, res=res)
        state.push(res)

    def op_CALL_FUNCTION_KW(self, state, inst):
        narg = inst.arg
        names = state.pop()  # tuple of names
        args = list(reversed([state.pop() for _ in range(narg)]))
        func = state.pop()

        res = state.make_temp()
        state.append(inst, func=func, args=args, names=names, res=res)
        state.push(res)

    if PYVERSION in ((3, 13),):
        def op_CALL_KW(self, state, inst):
            narg = inst.arg
            kw_names = state.pop()
            args = list(reversed([state.pop() for _ in range(narg)]))
            null_or_firstarg = state.pop()
            callable = state.pop()
            if not _is_null_temp_reg(null_or_firstarg):
                args = [null_or_firstarg, *args]

            res = state.make_temp()
            state.append(inst, func=callable, args=args, kw_names=kw_names,
                         res=res)
            state.push(res)

    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 13),):
        def op_CALL_FUNCTION_EX(self, state, inst):
            # (func, unused, callargs, kwargs if (oparg & 1) -- result))
            if inst.arg & 1:
                varkwarg = state.pop()
            else:
                varkwarg = None

            vararg = state.pop()
            state.pop()  # unused
            func = state.pop()

            res = state.make_temp()
            state.append(inst, func=func, vararg=vararg, varkwarg=varkwarg,
                         res=res)
            state.push(res)

    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):

        def op_CALL_FUNCTION_EX(self, state, inst):
            if inst.arg & 1:
                varkwarg = state.pop()
            else:
                varkwarg = None
            vararg = state.pop()
            func = state.pop()

            if PYVERSION in ((3, 11), (3, 12)):
                if _is_null_temp_reg(state.peek(1)):
                    state.pop() # pop NULL, it's not used
            elif PYVERSION in ((3, 10),):
                pass
            else:
                raise NotImplementedError(PYVERSION)

            res = state.make_temp()
            state.append(inst, func=func, vararg=vararg, varkwarg=varkwarg,
                         res=res)
            state.push(res)
    else:
        raise NotImplementedError(PYVERSION)

    def _dup_topx(self, state, inst, count):
        orig = [state.pop() for _ in range(count)]
        orig.reverse()
        # We need to actually create new temporaries if we want the
        # IR optimization pass to work correctly (see issue #580)
        duped = [state.make_temp() for _ in range(count)]
        state.append(inst, orig=orig, duped=duped)
        for val in orig:
            state.push(val)
        for val in duped:
            state.push(val)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_CALL_INTRINSIC_1(self, state, inst):
            # See https://github.com/python/cpython/blob/v3.12.0rc2/Include/
            # internal/pycore_intrinsics.h#L3-L17C36
            try:
                operand = CALL_INTRINSIC_1_Operand(inst.arg)
            except TypeError:
                msg = f"op_CALL_INTRINSIC_1({inst.arg})"
                loc = self.get_debug_loc(inst.lineno)
                raise UnsupportedBytecodeError(msg, loc=loc)
            if operand == ci1op.INTRINSIC_STOPITERATION_ERROR:
                state.append(inst, operand=operand)
                state.terminate()
                return
            elif operand == ci1op.UNARY_POSITIVE:
                val = state.pop()
                res = state.make_temp()
                state.append(inst, operand=operand,
                             value=val, res=res)
                state.push(res)
                return
            elif operand == ci1op.INTRINSIC_LIST_TO_TUPLE:
                tos = state.pop()
                res = state.make_temp()
                state.append(inst, operand=operand,
                             const_list=tos, res=res)
                state.push(res)
                return
            else:
                raise NotImplementedError(operand)

    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_DUP_TOPX(self, state, inst):
        count = inst.arg
        assert 1 <= count <= 5, "Invalid DUP_TOPX count"
        self._dup_topx(state, inst, count)

    def op_DUP_TOP(self, state, inst):
        self._dup_topx(state, inst, count=1)

    def op_DUP_TOP_TWO(self, state, inst):
        self._dup_topx(state, inst, count=2)

    def op_COPY(self, state, inst):
        state.push(state.peek(inst.arg))

    def op_SWAP(self, state, inst):
        state.swap(inst.arg)

    def op_ROT_TWO(self, state, inst):
        first = state.pop()
        second = state.pop()
        state.push(first)
        state.push(second)

    def op_ROT_THREE(self, state, inst):
        first = state.pop()
        second = state.pop()
        third = state.pop()
        state.push(first)
        state.push(third)
        state.push(second)

    def op_ROT_FOUR(self, state, inst):
        first = state.pop()
        second = state.pop()
        third = state.pop()
        forth = state.pop()
        state.push(first)
        state.push(forth)
        state.push(third)
        state.push(second)

    def op_UNPACK_SEQUENCE(self, state, inst):
        count = inst.arg
        iterable = state.pop()
        stores = [state.make_temp() for _ in range(count)]
        tupleobj = state.make_temp()
        state.append(inst, iterable=iterable, stores=stores, tupleobj=tupleobj)
        for st in reversed(stores):
            state.push(st)

    def op_BUILD_TUPLE(self, state, inst):
        count = inst.arg
        items = list(reversed([state.pop() for _ in range(count)]))
        tup = state.make_temp()
        state.append(inst, items=items, res=tup)
        state.push(tup)

    def _build_tuple_unpack(self, state, inst):
        # Builds tuple from other tuples on the stack
        tuples = list(reversed([state.pop() for _ in range(inst.arg)]))
        temps = [state.make_temp() for _ in range(len(tuples) - 1)]

        # if the unpack is assign-like, e.g. x = (*y,), it needs handling
        # differently.
        is_assign = len(tuples) == 1
        if is_assign:
            temps = [state.make_temp(),]

        state.append(inst, tuples=tuples, temps=temps, is_assign=is_assign)
        # The result is in the last temp var
        state.push(temps[-1])

    def op_BUILD_TUPLE_UNPACK_WITH_CALL(self, state, inst):
        # just unpack the input tuple, call inst will be handled afterwards
        self._build_tuple_unpack(state, inst)

    def op_BUILD_TUPLE_UNPACK(self, state, inst):
        self._build_tuple_unpack(state, inst)

    def op_LIST_TO_TUPLE(self, state, inst):
        # "Pops a list from the stack and pushes a tuple containing the same
        #  values."
        tos = state.pop()
        res = state.make_temp() # new tuple var
        state.append(inst, const_list=tos, res=res)
        state.push(res)

    def op_BUILD_CONST_KEY_MAP(self, state, inst):
        keys = state.pop()
        vals = list(reversed([state.pop() for _ in range(inst.arg)]))
        keytmps = [state.make_temp() for _ in range(inst.arg)]
        res = state.make_temp()
        state.append(inst, keys=keys, keytmps=keytmps, values=vals, res=res)
        state.push(res)

    def op_BUILD_LIST(self, state, inst):
        count = inst.arg
        items = list(reversed([state.pop() for _ in range(count)]))
        lst = state.make_temp()
        state.append(inst, items=items, res=lst)
        state.push(lst)

    def op_LIST_APPEND(self, state, inst):
        value = state.pop()
        index = inst.arg
        target = state.peek(index)
        appendvar = state.make_temp()
        res = state.make_temp()
        state.append(inst, target=target, value=value, appendvar=appendvar,
                     res=res)

    def op_LIST_EXTEND(self, state, inst):
        value = state.pop()
        index = inst.arg
        target = state.peek(index)
        extendvar = state.make_temp()
        res = state.make_temp()
        state.append(inst, target=target, value=value, extendvar=extendvar,
                     res=res)

    def op_BUILD_MAP(self, state, inst):
        dct = state.make_temp()
        count = inst.arg
        items = []
        # In 3.5+, BUILD_MAP takes <count> pairs from the stack
        for i in range(count):
            v, k = state.pop(), state.pop()
            items.append((k, v))
        state.append(inst, items=items[::-1], size=count, res=dct)
        state.push(dct)

    def op_MAP_ADD(self, state, inst):
        TOS = state.pop()
        TOS1 = state.pop()
        key, value = (TOS1, TOS)
        index = inst.arg
        target = state.peek(index)
        setitemvar = state.make_temp()
        res = state.make_temp()
        state.append(inst, target=target, key=key, value=value,
                     setitemvar=setitemvar, res=res)

    def op_BUILD_SET(self, state, inst):
        count = inst.arg
        # Note: related python bug http://bugs.python.org/issue26020
        items = list(reversed([state.pop() for _ in range(count)]))
        res = state.make_temp()
        state.append(inst, items=items, res=res)
        state.push(res)

    def op_SET_UPDATE(self, state, inst):
        value = state.pop()
        index = inst.arg
        target = state.peek(index)
        updatevar = state.make_temp()
        res = state.make_temp()
        state.append(inst, target=target, value=value, updatevar=updatevar,
                     res=res)

    def op_DICT_UPDATE(self, state, inst):
        value = state.pop()
        index = inst.arg
        target = state.peek(index)
        updatevar = state.make_temp()
        res = state.make_temp()
        state.append(inst, target=target, value=value, updatevar=updatevar,
                     res=res)

    def op_GET_ITER(self, state, inst):
        value = state.pop()
        res = state.make_temp()
        state.append(inst, value=value, res=res)
        state.push(res)

    def op_FOR_ITER(self, state, inst):
        iterator = state.get_tos()
        pair = state.make_temp()
        indval = state.make_temp()
        pred = state.make_temp()
        state.append(inst, iterator=iterator, pair=pair, indval=indval,
                     pred=pred)
        state.push(indval)
        end = inst.get_jump_target()
        if PYVERSION in ((3, 12), (3, 13)):
            # Changed in version 3.12: Up until 3.11 the iterator was
            # popped when it was exhausted. Now this is handled using END_FOR
            # op code.
            state.fork(pc=end)
        elif PYVERSION in ((3, 10), (3, 11)):
            state.fork(pc=end, npop=2)
        else:
            raise NotImplementedError(PYVERSION)
        state.fork(pc=inst.next)

    def op_GEN_START(self, state, inst):
        """Pops TOS. If TOS was not None, raises an exception. The kind
        operand corresponds to the type of generator or coroutine and
        determines the error message. The legal kinds are 0 for generator,
        1 for coroutine, and 2 for async generator.

        New in version 3.10.
        """
        # no-op in Numba
        pass

    def op_BINARY_OP(self, state, inst):
        op = dis._nb_ops[inst.arg][1]
        rhs = state.pop()
        lhs = state.pop()
        op_name = ALL_BINOPS_TO_OPERATORS[op].__name__
        res = state.make_temp(prefix=f"binop_{op_name}")
        state.append(inst, op=op, lhs=lhs, rhs=rhs, res=res)
        state.push(res)

    def _unaryop(self, state, inst):
        val = state.pop()
        res = state.make_temp()
        state.append(inst, value=val, res=res)
        state.push(res)

    op_UNARY_NEGATIVE = _unaryop
    op_UNARY_POSITIVE = _unaryop
    op_UNARY_NOT = _unaryop
    op_UNARY_INVERT = _unaryop

    def _binaryop(self, state, inst):
        rhs = state.pop()
        lhs = state.pop()
        res = state.make_temp()
        state.append(inst, lhs=lhs, rhs=rhs, res=res)
        state.push(res)

    op_COMPARE_OP = _binaryop
    op_IS_OP = _binaryop
    op_CONTAINS_OP = _binaryop

    op_INPLACE_ADD = _binaryop
    op_INPLACE_SUBTRACT = _binaryop
    op_INPLACE_MULTIPLY = _binaryop
    op_INPLACE_DIVIDE = _binaryop
    op_INPLACE_TRUE_DIVIDE = _binaryop
    op_INPLACE_FLOOR_DIVIDE = _binaryop
    op_INPLACE_MODULO = _binaryop
    op_INPLACE_POWER = _binaryop
    op_INPLACE_MATRIX_MULTIPLY = _binaryop

    op_INPLACE_LSHIFT = _binaryop
    op_INPLACE_RSHIFT = _binaryop
    op_INPLACE_AND = _binaryop
    op_INPLACE_OR = _binaryop
    op_INPLACE_XOR = _binaryop

    op_BINARY_ADD = _binaryop
    op_BINARY_SUBTRACT = _binaryop
    op_BINARY_MULTIPLY = _binaryop
    op_BINARY_DIVIDE = _binaryop
    op_BINARY_TRUE_DIVIDE = _binaryop
    op_BINARY_FLOOR_DIVIDE = _binaryop
    op_BINARY_MODULO = _binaryop
    op_BINARY_POWER = _binaryop
    op_BINARY_MATRIX_MULTIPLY = _binaryop

    op_BINARY_LSHIFT = _binaryop
    op_BINARY_RSHIFT = _binaryop
    op_BINARY_AND = _binaryop
    op_BINARY_OR = _binaryop
    op_BINARY_XOR = _binaryop

    def op_MAKE_FUNCTION(self, state, inst, MAKE_CLOSURE=False):
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            # https://github.com/python/cpython/commit/2f180ce
            # name set via co_qualname
            name = None
        elif PYVERSION in ((3, 10),):
            name = state.pop()
        else:
            raise NotImplementedError(PYVERSION)
        code = state.pop()
        closure = annotations = kwdefaults = defaults = None
        if PYVERSION in ((3, 13), ):
            assert inst.arg is None
            # SET_FUNCTION_ATTRIBUTE is responsible for setting
            # closure, annotations, kwdefaults and defaults.
        else:
            if inst.arg & 0x8:
                closure = state.pop()
            if inst.arg & 0x4:
                annotations = state.pop()
            if inst.arg & 0x2:
                kwdefaults = state.pop()
            if inst.arg & 0x1:
                defaults = state.pop()
        res = state.make_temp()
        state.append(
            inst,
            name=name,
            code=code,
            closure=closure,
            annotations=annotations,
            kwdefaults=kwdefaults,
            defaults=defaults,
            res=res,
        )
        state.push(res)

    def op_SET_FUNCTION_ATTRIBUTE(self, state, inst):
        assert PYVERSION in ((3, 13), )
        make_func_stack = state.pop()
        data = state.pop()
        if inst.arg == 0x1:
            # 0x01 a tuple of default values for positional-only and
            #      positional-or-keyword parameters in positional order
            state.set_function_attribute(make_func_stack, defaults=data)
        elif inst.arg & 0x2:
            # 0x02 a tuple of strings containing parameters’ annotations
            state.set_function_attribute(make_func_stack, kwdefaults=data)
        elif inst.arg & 0x4:
            # 0x04 a tuple of strings containing parameters’ annotations
            state.set_function_attribute(make_func_stack, annotations=data)
        elif inst.arg == 0x8:
            # 0x08 a tuple containing cells for free variables, making a closure
            state.set_function_attribute(make_func_stack, closure=data)
        else:
            raise AssertionError("unreachable")
        state.push(make_func_stack)

    def op_MAKE_CLOSURE(self, state, inst):
        self.op_MAKE_FUNCTION(state, inst, MAKE_CLOSURE=True)

    def op_LOAD_CLOSURE(self, state, inst):
        res = state.make_temp()
        state.append(inst, res=res)
        state.push(res)

    def op_LOAD_ASSERTION_ERROR(self, state, inst):
        res = state.make_temp("assertion_error")
        state.append(inst, res=res)
        state.push(res)

    def op_CHECK_EXC_MATCH(self, state, inst):
        pred = state.make_temp("predicate")
        tos = state.pop()
        tos1 = state.get_tos()
        state.append(inst, pred=pred, tos=tos, tos1=tos1)
        state.push(pred)

    def op_JUMP_IF_NOT_EXC_MATCH(self, state, inst):
        # Tests whether the second value on the stack is an exception matching
        # TOS, and jumps if it is not. Pops two values from the stack.
        pred = state.make_temp("predicate")
        tos = state.pop()
        tos1 = state.pop()
        state.append(inst, pred=pred, tos=tos, tos1=tos1)
        state.fork(pc=inst.next)
        state.fork(pc=inst.get_jump_target())

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_RERAISE(self, state, inst):
            # This isn't handled, but the state is set up anyway
            exc = state.pop()
            if inst.arg != 0:
                state.pop()     # lasti
            state.append(inst, exc=exc)

            if state.has_active_try():
                self._adjust_except_stack(state)
            else:
                state.terminate()

    elif PYVERSION in ((3, 10),):
        def op_RERAISE(self, state, inst):
            # This isn't handled, but the state is set up anyway
            exc = state.pop()
            state.append(inst, exc=exc)
            state.terminate()
    else:
        raise NotImplementedError(PYVERSION)

    # NOTE: Please see notes in `interpreter.py` surrounding the implementation
    # of LOAD_METHOD and CALL_METHOD.

    if PYVERSION in ((3, 12), (3, 13)):
        # LOAD_METHOD has become a pseudo-instruction in 3.12
        pass
    elif PYVERSION in ((3, 11), ):
        def op_LOAD_METHOD(self, state, inst):
            item = state.pop()
            extra = state.make_null()
            state.push(extra)
            res = state.make_temp()
            state.append(inst, item=item, res=res)
            state.push(res)
    elif PYVERSION in ((3, 10),):
        def op_LOAD_METHOD(self, state, inst):
            self.op_LOAD_ATTR(state, inst)
    else:
        raise NotImplementedError(PYVERSION)

    def op_CALL_METHOD(self, state, inst):
        self.op_CALL_FUNCTION(state, inst)


@total_ordering
class _State(object):
    """State of the trace
    """
    def __init__(self, bytecode, pc, nstack, blockstack, nullvals=()):
        """
        Parameters
        ----------
        bytecode : numba.bytecode.ByteCode
            function bytecode
        pc : int
            program counter
        nstack : int
            stackdepth at entry
        blockstack : Sequence[Dict]
            A sequence of dictionary denoting entries on the blockstack.
        """
        self._bytecode = bytecode
        self._pc_initial = pc
        self._pc = pc
        self._nstack_initial = nstack
        self._stack = []
        self._blockstack_initial = tuple(blockstack)
        self._blockstack = list(blockstack)
        self._temp_registers = []
        self._insts = []
        self._outedges = []
        self._terminated = False
        self._phis = {}
        self._outgoing_phis = UniqueDict()
        self._used_regs = set()
        for i in range(nstack):
            if i in nullvals:
                phi = self.make_temp("null$")
            else:
                phi = self.make_temp("phi")
            self._phis[phi] = i
            self.push(phi)

    def __repr__(self):
        return "State(pc_initial={} nstack_initial={})".format(
            self._pc_initial, self._nstack_initial
        )

    def get_identity(self):
        return (self._pc_initial, self._nstack_initial)

    def __hash__(self):
        return hash(self.get_identity())

    def __lt__(self, other):
        return self.get_identity() < other.get_identity()

    def __eq__(self, other):
        return self.get_identity() == other.get_identity()

    @property
    def pc_initial(self):
        """The starting bytecode offset of this State.
        The PC given to the constructor.
        """
        return self._pc_initial

    @property
    def instructions(self):
        """The list of instructions information as a 2-tuple of
        ``(pc : int, register_map : Dict)``
        """
        return self._insts

    @property
    def outgoing_edges(self):
        """The list of outgoing edges.

        Returns
        -------
        edges : List[State]
        """
        return self._outedges

    @property
    def outgoing_phis(self):
        """The dictionary of outgoing phi nodes.

        The keys are the name of the PHI nodes.
        The values are the outgoing states.
        """
        return self._outgoing_phis

    @property
    def blockstack_initial(self):
        """A copy of the initial state of the blockstack
        """
        return self._blockstack_initial

    @property
    def stack_depth(self):
        """The current size of the stack

        Returns
        -------
        res : int
        """
        return len(self._stack)

    def find_initial_try_block(self):
        """Find the initial *try* block.
        """
        for blk in reversed(self._blockstack_initial):
            if blk['kind'] == BlockKind('TRY'):
                return blk

    def has_terminated(self):
        return self._terminated

    def get_inst(self):
        return self._bytecode[self._pc]

    def advance_pc(self):
        inst = self.get_inst()
        self._pc = inst.next

    def make_temp(self, prefix=""):
        if not prefix:
            name = "${prefix}{offset}{opname}.{tempct}".format(
                prefix=prefix,
                offset=self._pc,
                opname=self.get_inst().opname.lower(),
                tempct=len(self._temp_registers),
            )
        else:
            name = "${prefix}{offset}.{tempct}".format(
                prefix=prefix,
                offset=self._pc,
                tempct=len(self._temp_registers),
            )

        self._temp_registers.append(name)
        return name

    def append(self, inst, **kwargs):
        """Append new inst"""
        self._insts.append((inst.offset, kwargs))
        self._used_regs |= set(_flatten_inst_regs(kwargs.values()))

    def get_tos(self):
        return self.peek(1)

    def peek(self, k):
        """Return the k'th element on the stack
        """
        return self._stack[-k]

    def push(self, item):
        """Push to stack"""
        self._stack.append(item)

    def pop(self):
        """Pop the stack"""
        return self._stack.pop()

    def swap(self, idx):
        """Swap stack[idx] with the tos"""
        s = self._stack
        s[-1], s[-idx] = s[-idx], s[-1]

    def push_block(self, synblk):
        """Push a block to blockstack
        """
        assert 'stack_depth' in synblk
        self._blockstack.append(synblk)

    def reset_stack(self, depth):
        """Reset the stack to the given stack depth.
        Returning the popped items.
        """
        self._stack, popped = self._stack[:depth], self._stack[depth:]
        return popped

    def make_block(self, kind, end, reset_stack=True, handler=None):
        """Make a new block
        """
        d = {
            'kind': BlockKind(kind),
            'end': end,
            'entry_stack': len(self._stack),
        }
        if reset_stack:
            d['stack_depth'] = len(self._stack)
        else:
            d['stack_depth'] = None
        d['handler'] = handler
        return d

    def pop_block(self):
        """Pop a block and unwind the stack
        """
        b = self._blockstack.pop()
        self.reset_stack(b['stack_depth'])
        return b

    def pop_block_and_above(self, blk):
        """Find *blk* in the blockstack and remove it and all blocks above it
        from the stack.
        """
        idx = self._blockstack.index(blk)
        assert 0 <= idx < len(self._blockstack)
        self._blockstack = self._blockstack[:idx]

    def get_top_block(self, kind):
        """Find the first block that matches *kind*
        """
        kind = BlockKind(kind)
        for bs in reversed(self._blockstack):
            if bs['kind'] == kind:
                return bs

    def get_top_block_either(self, *kinds):
        """Find the first block that matches *kind*
        """
        kinds = {BlockKind(kind) for kind in kinds}
        for bs in reversed(self._blockstack):
            if bs['kind'] in kinds:
                return bs

    def has_active_try(self):
        """Returns a boolean indicating if the top-block is a *try* block
        """
        return self.get_top_block('TRY') is not None

    def get_varname(self, inst):
        """Get referenced variable name from the instruction's oparg
        """
        return self.get_varname_by_arg(inst.arg)

    def get_varname_by_arg(self, oparg: int):
        """Get referenced variable name from the oparg
        """
        return self._bytecode.co_varnames[oparg]

    def terminate(self):
        """Mark block as terminated
        """
        self._terminated = True

    def fork(self, pc, npop=0, npush=0, extra_block=None):
        """Fork the state
        """
        # Handle changes on the stack
        stack = list(self._stack)
        if npop:
            assert 0 <= npop <= len(self._stack)
            nstack = len(self._stack) - npop
            stack = stack[:nstack]
        if npush:
            assert 0 <= npush
            for i in range(npush):
                stack.append(self.make_temp())
        # Handle changes on the blockstack
        blockstack = list(self._blockstack)
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            # pop expired block in destination pc
            while blockstack:
                top = blockstack[-1]
                end = top.get('end_offset') or top['end']
                if pc >= end:
                    blockstack.pop()
                else:
                    break
        elif PYVERSION in ((3, 10),):
            pass # intentionally bypass
        else:
            raise NotImplementedError(PYVERSION)

        if extra_block:
            blockstack.append(extra_block)
        self._outedges.append(Edge(
            pc=pc, stack=tuple(stack), npush=npush,
            blockstack=tuple(blockstack),
        ))
        self.terminate()

    def split_new_block(self):
        """Split the state
        """
        self.fork(pc=self._pc)

    def get_outgoing_states(self):
        """Get states for each outgoing edges
        """
        # Should only call once
        assert not self._outgoing_phis
        ret = []
        for edge in self._outedges:
            state = State(bytecode=self._bytecode, pc=edge.pc,
                          nstack=len(edge.stack), blockstack=edge.blockstack,
                          nullvals=[i for i, v in enumerate(edge.stack)
                                    if _is_null_temp_reg(v)])
            ret.append(state)
            # Map outgoing_phis
            for phi, i in state._phis.items():
                self._outgoing_phis[phi] = edge.stack[i]
        return ret

    def get_outgoing_edgepushed(self):
        """
        Returns
        -------
        Dict[int, int]
            where keys are the PC
            values are the edge-pushed stack values
        """

        return {edge.pc: tuple(edge.stack[-edge.npush:])
                for edge in self._outedges}


class StatePy311(_State):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._kw_names = None

    def pop_kw_names(self):
        out = self._kw_names
        self._kw_names = None
        return out

    def set_kw_names(self, val):
        assert self._kw_names is None
        self._kw_names = val

    def is_in_exception(self):
        bc = self._bytecode
        return bc.find_exception_entry(self._pc) is not None

    def get_exception(self):
        bc = self._bytecode
        return bc.find_exception_entry(self._pc)

    def in_with(self):
        for ent in self._blockstack_initial:
            if ent["kind"] == BlockKind("WITH"):
                return True

    def make_null(self):
        return self.make_temp(prefix="null$")


class StatePy313(StatePy311):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._make_func_attrs = defaultdict(dict)

    def set_function_attribute(self, make_func_res, **kwargs):
        self._make_func_attrs[make_func_res].update(kwargs)

    def get_function_attributes(self, make_func_res):
        return self._make_func_attrs[make_func_res]


if PYVERSION in ((3, 13), ):
    State = StatePy313
elif PYVERSION in ((3, 11), (3, 12)):
    State = StatePy311
elif PYVERSION < (3, 11):
    State = _State
else:
    raise NotImplementedError(PYVERSION)


Edge = namedtuple("Edge", ["pc", "stack", "blockstack", "npush"])


class AdaptDFA(object):
    """Adapt Flow to the old DFA class expected by Interpreter
    """
    def __init__(self, flow):
        self._flow = flow

    @property
    def infos(self):
        return self._flow.block_infos


AdaptBlockInfo = namedtuple(
    "AdaptBlockInfo",
    ["insts", "outgoing_phis", "blockstack", "active_try_block",
     "outgoing_edgepushed"],
)


def adapt_state_infos(state):
    def process_function_attributes(inst_pair):
        offset, data = inst_pair
        inst = state._bytecode[offset]
        if inst.opname == "MAKE_FUNCTION":
            data.update(state.get_function_attributes(data['res']))
        return offset, data
    if PYVERSION in ((3, 13), ):
        insts = tuple(map(process_function_attributes, state.instructions))
    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
        insts = tuple(state.instructions)
    else:
        raise NotImplementedError(PYVERSION)
    return AdaptBlockInfo(
        insts=insts,
        outgoing_phis=state.outgoing_phis,
        blockstack=state.blockstack_initial,
        active_try_block=state.find_initial_try_block(),
        outgoing_edgepushed=state.get_outgoing_edgepushed(),
    )


def _flatten_inst_regs(iterable):
    """Flatten an iterable of registers used in an instruction
    """
    for item in iterable:
        if isinstance(item, str):
            yield item
        elif isinstance(item, (tuple, list)):
            for x in _flatten_inst_regs(item):
                yield x


class AdaptCFA(object):
    """Adapt Flow to the old CFA class expected by Interpreter
    """
    def __init__(self, flow):
        self._flow = flow
        self._blocks = {}
        for offset, blockinfo in flow.block_infos.items():
            self._blocks[offset] = AdaptCFBlock(blockinfo, offset)
        backbone = self._flow.cfgraph.backbone()

        graph = flow.cfgraph
        # Find backbone
        backbone = graph.backbone()
        # Filter out in loop blocks (Assuming no other cyclic control blocks)
        # This is to unavoid variables defined in loops being considered as
        # function scope.
        inloopblocks = set()
        for b in self.blocks.keys():
            if graph.in_loops(b):
                inloopblocks.add(b)
        self._backbone = backbone - inloopblocks

    @property
    def graph(self):
        return self._flow.cfgraph

    @property
    def backbone(self):
        return self._backbone

    @property
    def blocks(self):
        return self._blocks

    def iterliveblocks(self):
        for b in sorted(self.blocks):
            yield self.blocks[b]

    def dump(self):
        self._flow.cfgraph.dump()


class AdaptCFBlock(object):
    def __init__(self, blockinfo, offset):
        self.offset = offset
        self.body = tuple(i for i, _ in blockinfo.insts)
