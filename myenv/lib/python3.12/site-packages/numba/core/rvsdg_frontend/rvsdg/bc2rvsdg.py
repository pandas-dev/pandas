import os
from dataclasses import dataclass, replace, field, fields
import dis
import operator
from functools import reduce
from typing import (
    Optional,
    Protocol,
    runtime_checkable,
    Union,
    NamedTuple,
    Mapping,
    TypeVar,
    TYPE_CHECKING,
)

from collections import ChainMap
from numba_rvsdg.core.datastructures.byte_flow import ByteFlow
from numba_rvsdg.core.datastructures.scfg import SCFG
from numba_rvsdg.core.datastructures.basic_block import (
    BasicBlock,
    PythonBytecodeBlock,
    RegionBlock,
    SyntheticBranch,
    SyntheticAssignment,
)
from numba_rvsdg.rendering.rendering import ByteFlowRenderer
from numba_rvsdg.core.datastructures import block_names


from numba.core.utils import MutableSortedSet, MutableSortedMap

from .regionpasses import (
    RegionVisitor,
    RegionTransformer,
)

if TYPE_CHECKING:
    from .regionrenderer import GraphBuilder


DEBUG_GRAPH = int(os.environ.get("DEBUG_GRAPH", "0"))


T = TypeVar("T")


def _just(v: Optional[T]) -> T:
    """Unpack optional type."""
    assert v is not None
    return v


@dataclass(frozen=True)
class ValueState:
    """For representing a RVSDG Value (State).

    For most compiler passes, Value and State can be treated as the same.
    """
    parent: Optional["Op"]
    """Optional. The parent Op that output this ValueState.
    """
    name: str
    """Name of the Value(State).
    """
    out_index: int
    """The output port index in the parent Op.
    """
    is_effect: bool = False
    """True if-and-only-if this is a state.
    """

    def short_identity(self) -> str:
        args = f"{id(self.parent):x}, {self.name}, {self.out_index}"
        return f"ValueState({args})"

    def __hash__(self):
        return id(self)


@dataclass(frozen=True)
class Op:
    """For representing a RVSDG operation.

    An Op has inputs and outputs ports that are ValueStates.
    An Op can optionally have a bytecode instruction associated with it.
    """
    opname: str
    """Operation name.
    """
    bc_inst: Optional[dis.Instruction]
    """Optional. The bytecode instruction.
    """
    _inputs: dict[str, ValueState] = field(default_factory=dict)
    """The input ports.
    """
    _outputs: dict[str, ValueState] = field(default_factory=dict)
    """The output ports.
    """

    def add_input(self, name, vs: ValueState):
        self._inputs[name] = vs

    def add_output(self, name: str, is_effect=False) -> ValueState:
        vs = ValueState(
            parent=self,
            name=name,
            out_index=len(self._outputs),
            is_effect=is_effect,
        )
        self._outputs[name] = vs
        return vs

    def short_identity(self) -> str:
        return f"Op({self.opname}, {id(self):x})"

    def summary(self) -> str:
        ins = ", ".join([k for k in self._inputs])
        outs = ", ".join([k for k in self._outputs])
        bc = "---"
        if self.bc_inst is not None:
            bc = f"{self.bc_inst.opname}({self.bc_inst.argrepr})"
        return f"Op\n{self.opname}\n{bc}\n({ins}) -> ({outs}) "

    @property
    def outputs(self) -> list[ValueState]:
        return list(self._outputs.values())

    @property
    def inputs(self) -> list[ValueState]:
        return list(self._inputs.values())

    @property
    def input_ports(self) -> Mapping[str, ValueState]:
        return self._inputs

    @property
    def output_ports(self) -> Mapping[str, ValueState]:
        return self._outputs

    def __hash__(self):
        return id(self)


@runtime_checkable
class DDGProtocol(Protocol):
    incoming_states: MutableSortedSet[str]
    outgoing_states: MutableSortedSet[str]


@dataclass(frozen=True)
class DDGRegion(RegionBlock):
    incoming_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )
    outgoing_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )


@dataclass(frozen=True)
class DDGBranch(SyntheticBranch):
    incoming_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )
    outgoing_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )


@dataclass(frozen=True)
class DDGControlVariable(SyntheticAssignment):
    incoming_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )
    outgoing_states: MutableSortedSet[str] = field(
        default_factory=MutableSortedSet
    )


@dataclass(frozen=True)
class DDGBlock(BasicBlock):
    in_effect: ValueState | None = None
    out_effect: ValueState | None = None
    in_stackvars: list[ValueState] = field(default_factory=list)
    out_stackvars: list[ValueState] = field(default_factory=list)
    in_vars: MutableSortedMap[str, ValueState] = field(
        default_factory=MutableSortedMap
    )
    out_vars: MutableSortedMap[str, ValueState] = field(
        default_factory=MutableSortedMap
    )

    exported_stackvars: MutableSortedMap[str, ValueState] = field(
        default_factory=MutableSortedMap
    )

    def __post_init__(self):
        assert isinstance(self.in_vars, MutableSortedMap)
        assert isinstance(self.out_vars, MutableSortedMap)

    def _gather_reachable(
        self, vs: ValueState, reached: set[ValueState]
    ) -> set[ValueState]:
        reached.add(vs)
        if vs.parent is not None:
            for ivs in vs.parent.inputs:
                if ivs not in reached:
                    self._gather_reachable(ivs, reached)
        return reached

    def render_graph(self, builder: "GraphBuilder"):
        reached_vs: set[ValueState] = set()
        for vs in [*self.out_vars.values(), _just(self.out_effect)]:
            self._gather_reachable(vs, reached_vs)
        reached_vs.add(_just(self.in_effect))
        reached_vs.update(self.in_vars.values())

        reached_op = {vs.parent for vs in reached_vs if vs.parent is not None}

        for vs in reached_vs:
            self._render_vs(builder, vs)

        for op in reached_op:
            self._render_op(builder, op)

        # Make outgoing node
        ports = []
        outgoing_nodename = f"outgoing_{self.name}"
        for k, vs in self.out_vars.items():
            ports.append(k)

            builder.graph.add_edge(
                vs.short_identity(), outgoing_nodename, dst_port=k
            )

        outgoing_node = builder.node_maker.make_node(
            kind="ports",
            ports=ports,
            data=dict(body="outgoing"),
        )
        builder.graph.add_node(outgoing_nodename, outgoing_node)

        # Make incoming node
        ports = []
        incoming_nodename = f"incoming_{self.name}"
        for k, vs in self.in_vars.items():
            ports.append(k)

            builder.graph.add_edge(
                incoming_nodename,
                _just(vs.parent).short_identity(),
                src_port=k,
            )

        incoming_node = builder.node_maker.make_node(
            kind="ports",
            ports=ports,
            data=dict(body="incoming"),
        )
        builder.graph.add_node(incoming_nodename, incoming_node)

        # Make jump target node for debugging
        jt_node = builder.node_maker.make_node(
            kind="meta",
            data=dict(body=f"jump-targets: {self._jump_targets}"),
        )
        builder.graph.add_node(f"jt_{self.name}", jt_node)

    def _render_vs(self, builder: "GraphBuilder", vs: ValueState):
        if vs.is_effect:
            node = builder.node_maker.make_node(
                kind="effect",
                data=dict(body=str(vs.name)),
            )
            builder.graph.add_node(vs.short_identity(), node)
        else:
            node = builder.node_maker.make_node(
                kind="valuestate",
                data=dict(body=str(vs.name)),
            )
            builder.graph.add_node(vs.short_identity(), node)

    def _render_op(self, builder, op: Op):
        op_anchor = op.short_identity()

        node = builder.node_maker.make_node(
            kind="op",
            data=dict(body=str(op.summary())),
        )
        builder.graph.add_node(op_anchor, node)

        # draw edge
        for edgename, vs in op._outputs.items():
            self._add_vs_edge(builder, op_anchor, vs, taillabel=f"{edgename}")
        for edgename, vs in op._inputs.items():
            self._add_vs_edge(builder, vs, op_anchor, headlabel=f"{edgename}")

    def _add_vs_edge(self, builder, src, dst, **attrs):
        is_effect = (isinstance(src, ValueState) and src.is_effect) or (
            isinstance(dst, ValueState) and dst.is_effect
        )
        if isinstance(src, ValueState):
            src = src.short_identity()
        if isinstance(dst, ValueState):
            dst = dst.short_identity()

        kwargs = attrs
        if is_effect:
            kwargs["kind"] = "effect"

        builder.graph.add_edge(src, dst, **kwargs)

    @property
    def incoming_states(self) -> MutableSortedSet:
        return MutableSortedSet(self.in_vars)

    @property
    def outgoing_states(self) -> MutableSortedSet:
        return MutableSortedSet(self.out_vars)

    def get_toposorted_ops(self) -> list[Op]:
        """Get a topologically sorted list of ``Op`` according
        to the data-dependence.

        Operations stored later in the list may depend on earlier operations,
        but the reverse can never be true.
        """
        res: list[Op] = []

        # Initially, the available states are the inputs
        avail: set[ValueState] = {*self.in_vars.values(), _just(self.in_effect)}
        # Seed the pending Operations using the outputs.
        pending: list[Op] = [
            vs.parent for vs in self.out_vars.values() if vs.parent is not None
        ]
        assert self.out_effect is not None  # for typing
        pending.append(_just(self.out_effect.parent))
        seen: set[Op] = set()

        while pending:
            op = pending[-1]
            if op in seen:
                pending.pop()
                continue
            # Get the set of incoming states that are not yet available
            # to the current operation.
            # NOTE: for stable ordering, change the following to a ordered-set.
            incomings = set()
            for vs in op._inputs.values():
                if vs not in avail and vs.parent is not None:
                    incomings.add(vs.parent)

            if not incomings:
                # All incoming states are already available
                avail |= set(op._outputs.values())
                pending.pop()
                res.append(op)
                seen.add(op)
            else:
                pending.extend(incomings)
        return res


def render_scfg(byteflow):
    bfr = ByteFlowRenderer()
    bfr.bcmap_from_bytecode(byteflow.bc)
    byteflow.scfg.view("scfg")


def _repl_jump_targets(block: BasicBlock, repl: dict[str, str]):
    if set(repl).intersection(set(block.jump_targets)):
        targets = [repl.get(k, k) for k in block.jump_targets]
        block = block.replace_jump_targets(tuple(targets))
    return block


def _canonicalize_scfg_switch(scfg: SCFG):
    """Introduce "switch" region to enclose "head", "branch", "tail" regions."""
    todos = set(scfg.graph)
    while todos:
        label = todos.pop()
        todos.discard(label)
        block = scfg[label]
        if isinstance(block, RegionBlock):
            if block.kind == "head":
                brlabels = block.jump_targets
                branches = [scfg[brlabel] for brlabel in brlabels]
                tail_label_candidates = set()
                for br in branches:
                    tail_label_candidates.update(br.jump_targets)
                [taillabel] = tail_label_candidates
                tail = scfg[taillabel]

                # Make new SCFG
                switch_labels = {label, taillabel, *brlabels}
                subregion_graph = {k: scfg[k] for k in switch_labels}
                scfg.remove_blocks(switch_labels)
                subregion_scfg = SCFG(
                    graph=subregion_graph, name_gen=scfg.name_gen
                )

                todos -= switch_labels

                # Make new region
                new_label = scfg.name_gen.new_region_name("switch")
                new_region = RegionBlock(
                    name=new_label,
                    _jump_targets=tail._jump_targets,
                    kind="switch",
                    parent_region=block.parent_region,
                    header=block.name,
                    exiting=taillabel,
                    subregion=subregion_scfg,
                )
                scfg.graph[new_label] = new_region
                # fixup branch to label
                for incoming_label, incoming_blk in scfg.graph.items():
                    if (
                        incoming_label != new_label
                        and label in incoming_blk.jump_targets
                    ):
                        repl = {label: new_label}
                        replblk = _repl_jump_targets(incoming_blk, repl)
                        scfg.graph[incoming_label] = replblk

                # fixup header
                if block.parent_region.header not in scfg.graph:
                    block.parent_region.replace_header(new_label)
                # fixup exiting
                if block.parent_region.exiting not in scfg.graph:
                    block.parent_region.replace_exiting(new_label)

                # recursively walk into the subregions
                _canonicalize_scfg_switch(subregion_graph[label].subregion)
                for br in brlabels:
                    _canonicalize_scfg_switch(subregion_graph[br].subregion)
                _canonicalize_scfg_switch(subregion_graph[taillabel].subregion)
            elif block.kind == "loop":
                _canonicalize_scfg_switch(block.subregion)


class CanonicalizeLoop(RegionTransformer[None]):
    """
    Make sure loops has non-region header.

    Preferably, the exiting block should be non-region as well but
    it's hard to do with the current numba_rvsdg API.

    Doing this so we don't have to fixup backedges as backedges will always
    point to a non-region node in ``_canonicalize_scfg_switch``.
    """

    def visit_loop(self, parent: SCFG, region: RegionBlock, data: None):
        # Fix header
        # Introduce a SyntheticFill block that just jump to the original header
        new_label = parent.name_gen.new_block_name(block_names.SYNTH_FILL)
        region.subregion.insert_SyntheticFill(new_label, {}, {region.header})
        region.replace_header(new_label)

        # Fix exiting block's backedges to point to the new header
        def get_inner_most_exiting(blk):
            while isinstance(blk, RegionBlock):
                parent, blk = blk, blk.subregion.graph[blk.exiting]
            return parent, blk

        tail_parent, tail_bb = get_inner_most_exiting(region)
        [backedge] = tail_bb.backedges
        repl = {backedge: new_label}
        new_tail_bb = replace(
            tail_bb,
            backedges=(new_label,),
            _jump_targets=tuple(
                [repl.get(x, x) for x in tail_bb._jump_targets]
            ),
        )
        tail_parent.subregion.graph[tail_bb.name] = new_tail_bb

        self.visit_linear(parent, region, data)

    def visit_block(self, parent: SCFG, block: BasicBlock, data: None):
        pass  # no-op

    def visit_switch(self, parent: SCFG, block: BasicBlock, data: None):
        pass  # no-op


def canonicalize_scfg(scfg: SCFG):
    CanonicalizeLoop().visit_graph(scfg, None)
    _canonicalize_scfg_switch(scfg)


class _ExtraBranch(NamedTuple):
    branch_instlists: tuple[tuple[str, ...], ...] = ()


@dataclass(frozen=True)
class ExtraBasicBlock(BasicBlock):
    inst_list: tuple[str, ...] = ()

    @classmethod
    def make(cls, label, jump_target, instlist):
        return ExtraBasicBlock(label, (jump_target,), inst_list=instlist)

    def __str__(self):
        args = "\n".join(f"{inst})" for inst in self.inst_list)
        return f"ExtraBasicBlock({args})"


class HandleConditionalPop:
    """Introduce pop-stack operations to the bytecode to correctly model
    operations that conditionally pop elements from the stack. Numba-rvsdg does
    not handle this. For example, FOR_ITER pop the stack when the iterator is
    exhausted.
    """

    def handle(self, inst: dis.Instruction) -> _ExtraBranch | None:
        fn = getattr(self, f"op_{inst.opname}", self._op_default)
        return fn(inst)

    def _op_default(self, inst: dis.Instruction) -> None:
        # Being defensive
        assert not inst.opname.endswith("OR_POP")
        return

    def op_FOR_ITER(self, inst: dis.Instruction) -> _ExtraBranch:
        # Bytecode semantic pop 1 for the iterator
        # but we need an extra pop because the indvar is pushed
        br0 = ("FOR_ITER_STORE_INDVAR",)
        br1 = ("POP",)
        return _ExtraBranch((br0, br1))

    def op_JUMP_IF_TRUE_OR_POP(self, inst: dis.Instruction) -> _ExtraBranch:
        br0 = ("POP",)
        br1 = ()
        return _ExtraBranch((br0, br1))

    def op_JUMP_IF_FALSE_OR_POP(self, inst: dis.Instruction) -> _ExtraBranch:
        br0 = ("POP",)
        br1 = ()
        return _ExtraBranch((br0, br1))


def _scfg_add_conditional_pop_stack(bcmap, scfg: SCFG):
    extra_records = {}
    for blk in scfg.graph.values():
        if isinstance(blk, PythonBytecodeBlock):
            # Check if last instruction has conditional pop
            last_inst = blk.get_instructions(bcmap)[-1]
            handler = HandleConditionalPop()
            res = handler.handle(last_inst)
            if res is not None:
                for br_index, instlist in enumerate(res.branch_instlists):
                    k = blk._jump_targets[br_index]
                    extra_records[k] = blk.name, (br_index, instlist)

    def _replace_jump_targets(blk, idx, repl):
        return replace(
            blk,
            _jump_targets=tuple(
                [
                    repl if i == idx else jt
                    for i, jt in enumerate(blk._jump_targets)
                ]
            ),
        )

    for label, (parent_label, (br_index, instlist)) in extra_records.items():
        newlabel = scfg.name_gen.new_block_name("python.extrabasicblock")
        scfg.graph[parent_label] = _replace_jump_targets(
            scfg.graph[parent_label], br_index, newlabel
        )
        ebb = ExtraBasicBlock.make(newlabel, label, instlist)
        scfg.graph[newlabel] = ebb


def build_rvsdg(code, argnames: tuple[str, ...]) -> SCFG:
    byteflow = ByteFlow.from_bytecode(code)
    bcmap = byteflow.scfg.bcmap_from_bytecode(byteflow.bc)
    _scfg_add_conditional_pop_stack(bcmap, byteflow.scfg)
    byteflow = byteflow.restructure()
    # if DEBUG_GRAPH:
    #     render_scfg(byteflow)
    canonicalize_scfg(byteflow.scfg)
    if DEBUG_GRAPH:
        render_scfg(byteflow)
    rvsdg = convert_to_dataflow(byteflow, argnames)
    rvsdg = propagate_states(rvsdg)
    if DEBUG_GRAPH:
        from .regionrenderer import RVSDGRenderer, to_graphviz

        to_graphviz(RVSDGRenderer().render(rvsdg)).view("rvsdg")

    return rvsdg


def _flatten_full_graph(scfg: SCFG):
    regions = [
        _flatten_full_graph(elem.subregion)
        for elem in scfg.graph.values()
        if isinstance(elem, RegionBlock)
    ]
    out = ChainMap(*regions, scfg.graph)
    for blk in out.values():
        assert not isinstance(blk, RegionBlock), type(blk)
    return out


DDGTypes = (DDGBlock, DDGControlVariable, DDGBranch)
_DDGTypeAnn = Union[DDGBlock, DDGControlVariable, DDGBranch]


def convert_to_dataflow(byteflow: ByteFlow, argnames: tuple[str, ...]) -> SCFG:
    bcmap = {inst.offset: inst for inst in byteflow.bc}
    rvsdg = convert_scfg_to_dataflow(byteflow.scfg, bcmap, argnames)
    return rvsdg


def propagate_states(rvsdg: SCFG) -> SCFG:
    propagate_stack(rvsdg)
    propagate_vars(rvsdg)
    connect_incoming_stack_vars(rvsdg)
    return rvsdg


def propagate_vars(rvsdg: SCFG):
    # Propagate variables
    visitor = PropagateVars()
    visitor.visit_graph(rvsdg, visitor.make_data())


_pvData = set[str]


class PropagateVars(RegionVisitor[_pvData]):
    """
    Depends on PropagateStack
    """

    def __init__(self, _debug: bool = False):
        super().__init__()
        self._debug = _debug

    def debug_print(self, *args, **kwargs):
        if self._debug:
            print(*args, **kwargs)

    def _apply(self, block: BasicBlock, data: _pvData) -> _pvData:
        assert isinstance(block, BasicBlock)
        if isinstance(block, DDGProtocol):
            if isinstance(block, DDGBlock):
                for k in data:
                    if k not in block.in_vars:
                        op = Op(opname="var.incoming", bc_inst=None)
                        vs = op.add_output(k)
                        block.in_vars[k] = vs
                        if k.startswith("tos."):
                            if k in block.exported_stackvars:
                                block.out_vars[k] = block.exported_stackvars[k]
                        elif k not in block.out_vars:
                            block.out_vars[k] = vs
            else:
                block.incoming_states.update(data)
                block.outgoing_states.update(data)

            data = set(block.outgoing_states)
            return data
        else:
            return data

    def visit_linear(self, region: RegionBlock, data: _pvData) -> _pvData:
        region.incoming_states.update(data)
        data = self.visit_graph(region.subregion, data)
        region.outgoing_states.update(data)
        return set(region.outgoing_states)

    def visit_block(self, block: BasicBlock, data: _pvData) -> _pvData:
        return self._apply(block, data)

    def visit_loop(self, region: RegionBlock, data: _pvData) -> _pvData:
        self.debug_print("---LOOP_ENTER", region.name, data)
        data = self.visit_linear(region, data)
        self.debug_print("---LOOP_END=", region.name, "vars", data)
        return data

    def visit_switch(self, region: RegionBlock, data: _pvData) -> _pvData:
        self.debug_print("---SWITCH_ENTER", region.name)
        region.incoming_states.update(data)
        header = region.header
        data_at_head = self.visit_linear(region.subregion[header], data)
        data_for_branches = []
        for blk in region.subregion.graph.values():
            if blk.kind == "branch":
                data_for_branches.append(self.visit_linear(blk, data_at_head))
        data_after_branches = reduce(operator.or_, data_for_branches)

        exiting = region.exiting

        data_at_tail = self.visit_linear(
            region.subregion[exiting], data_after_branches
        )

        self.debug_print("data_at_head", data_at_head)
        self.debug_print("data_for_branches", data_for_branches)
        self.debug_print("data_after_branches", data_after_branches)
        self.debug_print("data_at_tail", data_at_tail)
        self.debug_print("---SWITCH_END=", region.name, "vars", data_at_tail)
        region.outgoing_states.update(data_at_tail)
        return set(region.outgoing_states)

    def make_data(self) -> _pvData:
        return set()


_psData = tuple[str, ...]


class PropagateStack(RegionVisitor[_psData]):
    def __init__(self, _debug: bool = False):
        super(PropagateStack, self).__init__()
        self._debug = _debug

    def debug_print(self, *args, **kwargs):
        if self._debug:
            print(*args, **kwargs)

    def visit_block(self, block: BasicBlock, data: _psData) -> _psData:
        if isinstance(block, DDGBlock):
            nin = len(block.in_stackvars)
            inherited = data[: len(data) - nin]
            self.debug_print("--- stack", data)
            self.debug_print("--- inherited stack", inherited)
            out_stackvars = block.out_stackvars.copy()

            counts = reversed(list(range(len(inherited) + len(out_stackvars))))
            out_stack = block.out_stackvars[::-1]
            out_data = tuple([f"tos.{i}" for i in counts])

            unused_names = list(out_data)

            for vs in reversed(out_stack):
                k = unused_names.pop()
                op = Op("stack.export", bc_inst=None)
                op.add_input("0", vs)
                block.exported_stackvars[k] = vs = op.add_output(k)
                block.out_vars[k] = vs

            for orig in reversed(inherited):
                k = unused_names.pop()
                import_op = Op("var.incoming", bc_inst=None)
                block.in_vars[orig] = imported_vs = import_op.add_output(orig)

                op = Op("stack.export", bc_inst=None)
                op.add_input("0", imported_vs)
                vs = op.add_output(k)
                block.exported_stackvars[k] = vs
                block.out_vars[k] = vs

            self.debug_print("---=", block.name, "out stack", out_data)
            return out_data
        else:
            return data

    def visit_loop(self, region: RegionBlock, data: _psData) -> _psData:
        self.debug_print("---LOOP_ENTER", region.name)
        data = self.visit_linear(region, data)
        self.debug_print("---LOOP_END=", region.name, "stack", data)
        return data

    def visit_switch(self, region: RegionBlock, data: _psData) -> _psData:
        self.debug_print("---SWITCH_ENTER", region.name)
        header = region.header
        data_at_head = self.visit_linear(region.subregion[header], data)
        data_for_branches = []
        for blk in region.subregion.graph.values():
            if blk.kind == "branch":
                data_for_branches.append(self.visit_linear(blk, data_at_head))
        data_after_branches = max(data_for_branches, key=len)

        exiting = region.exiting

        data_at_tail = self.visit_linear(
            region.subregion[exiting], data_after_branches
        )

        self.debug_print("data_at_head", data_at_head)
        self.debug_print("data_for_branches", data_for_branches)
        self.debug_print("data_after_branches", data_after_branches)
        self.debug_print("data_at_tail", data_at_tail)
        self.debug_print("---SWITCH_END=", region.name, "stack", data_at_tail)
        return data_at_tail

    def make_data(self) -> _psData:
        # Stack data stored in a tuple
        return ()


class ConnectImportedStackVars(RegionVisitor[None]):
    def visit_block(self, block: BasicBlock, data: None):
        if isinstance(block, DDGBlock):
            # Connect stack.incoming node to the import stack variable.
            imported_stackvars = [
                var for k, var in block.in_vars.items() if k.startswith("tos.")
            ][::-1]
            n = len(block.in_stackvars)
            for vs, inc in zip(block.in_stackvars, imported_stackvars[-n:]):
                assert vs.parent is not None
                vs.parent.add_input("0", inc)

    def visit_loop(self, region: RegionBlock, data: None):
        self.visit_linear(region, data)

    def visit_switch(self, region: RegionBlock, data: None):
        header = region.header
        self.visit_linear(region.subregion[header], data)
        for blk in region.subregion.graph.values():
            if blk.kind == "branch":
                self.visit_linear(blk, None)
        exiting = region.exiting
        self.visit_linear(region.subregion[exiting], None)


def propagate_stack(rvsdg: SCFG):
    visitor = PropagateStack(_debug=False)
    visitor.visit_graph(rvsdg, visitor.make_data())


def connect_incoming_stack_vars(rvsdg: SCFG):
    ConnectImportedStackVars().visit_graph(rvsdg, None)


def _upgrade_dataclass(old, newcls, replacements=None):
    if replacements is None:
        replacements = {}
    fieldnames = [fd.name for fd in fields(old)]
    oldattrs = {k: getattr(old, k) for k in fieldnames if k not in replacements}
    return newcls(**oldattrs, **replacements)


def convert_scfg_to_dataflow(scfg, bcmap, argnames: tuple[str, ...]) -> SCFG:
    rvsdg = SCFG()
    for block in scfg.graph.values():
        # convert block
        if isinstance(block, PythonBytecodeBlock):
            ddg = convert_bc_to_ddg(block, bcmap, argnames)
            rvsdg.add_block(ddg)
        elif isinstance(block, RegionBlock):
            # Inside-out
            subregion = convert_scfg_to_dataflow(
                block.subregion, bcmap, argnames
            )
            rvsdg.add_block(
                _upgrade_dataclass(block, DDGRegion, dict(subregion=subregion))
            )
        elif isinstance(block, SyntheticBranch):
            rvsdg.add_block(_upgrade_dataclass(block, DDGBranch))
        elif isinstance(block, SyntheticAssignment):
            rvsdg.add_block(_upgrade_dataclass(block, DDGControlVariable))
        elif isinstance(block, ExtraBasicBlock):
            ddg = convert_extra_bb(block)
            rvsdg.add_block(ddg)
        elif isinstance(block, BasicBlock):
            start_env = Op("start", bc_inst=None)
            effect = start_env.add_output("env", is_effect=True)
            newblk = _upgrade_dataclass(
                block, DDGBlock, dict(in_effect=effect, out_effect=effect)
            )
            rvsdg.add_block(newblk)
        else:
            raise Exception("unreachable", type(block))

    return rvsdg


def _convert_bytecode(block, instlist, argnames: tuple[str, ...]) -> DDGBlock:
    converter = BC2DDG()
    if instlist[0].offset == 0:
        for arg in argnames:
            converter.load(f"var.{arg}")
    for inst in instlist:
        converter.convert(inst)
    return _converter_to_ddgblock(block, converter)


def _converter_to_ddgblock(block, converter) -> DDGBlock:
    blk = DDGBlock(
        name=block.name,
        _jump_targets=block._jump_targets,
        backedges=block.backedges,
        in_effect=converter.in_effect,
        out_effect=converter.effect,
        in_stackvars=list(converter.incoming_stackvars),
        out_stackvars=list(converter.stack),
        in_vars=MutableSortedMap(converter.incoming_vars),
        out_vars=MutableSortedMap(converter.varmap),
    )
    return blk


def convert_extra_bb(block: ExtraBasicBlock) -> DDGBlock:
    converter = BC2DDG()
    for opname in block.inst_list:
        if opname == "FOR_ITER_STORE_INDVAR":
            converter.push(converter.load("indvar"))
        elif opname == "POP":
            converter.pop()
        else:
            assert False, opname
    return _converter_to_ddgblock(block, converter)


def convert_bc_to_ddg(
    block: PythonBytecodeBlock,
    bcmap: dict[int, dis.Bytecode],
    argnames: tuple[str, ...],
) -> DDGBlock:
    instlist = block.get_instructions(bcmap)
    return _convert_bytecode(block, instlist, argnames)


class BC2DDG:
    stack: list[ValueState]
    effect: ValueState
    in_effect: ValueState
    varmap: dict[str, ValueState]
    incoming_vars: dict[str, ValueState]
    incoming_stackvars: list[ValueState]
    _kw_names: ValueState | None

    def __init__(self):
        self.stack = []
        start_env = Op("start", bc_inst=None)
        self.effect = start_env.add_output("env", is_effect=True)
        self.in_effect = self.effect
        self.varmap = {}
        self.incoming_vars = {}
        self.incoming_stackvars = []
        self._kw_names = None

    def push(self, val: ValueState):
        self.stack.append(val)

    def pop(self) -> ValueState:
        if not self.stack:
            op = Op(opname="stack.incoming", bc_inst=None)
            vs = op.add_output(f"stack.{len(self.incoming_stackvars)}")
            self.stack.append(vs)
            self.incoming_stackvars.append(vs)
        return self.stack.pop()

    def top(self) -> ValueState:
        tos = self.pop()
        self.push(tos)
        return tos

    def _decorate_varname(self, varname: str) -> str:
        return f"var.{varname}"

    def store(self, varname: str, value: ValueState):
        self.varmap[varname] = value

    def load(self, varname: str) -> ValueState:
        if varname not in self.varmap:
            op = Op(opname="var.incoming", bc_inst=None)
            vs = op.add_output(varname)
            self.incoming_vars[varname] = vs
            self.varmap[varname] = vs

        return self.varmap[varname]

    def replace_effect(self, env: ValueState):
        assert env.is_effect
        self.effect = env

    def convert(self, inst: dis.Instruction):
        fn = getattr(self, f"op_{inst.opname}")
        fn(inst)

    def set_kw_names(self, kw_vs: ValueState):
        assert self._kw_names is None
        self._kw_names = kw_vs

    def pop_kw_names(self):
        res = self._kw_names
        self._kw_names = None
        return res

    def op_POP_TOP(self, inst: dis.Instruction):
        self.pop()

    def op_RESUME(self, inst: dis.Instruction):
        pass  # no-op

    def op_COPY_FREE_VARS(self, inst: dis.Instruction):
        pass  # no-op

    def op_PUSH_NULL(self, inst: dis.Instruction):
        op = Op(opname="push_null", bc_inst=inst)
        null = op.add_output("null")
        self.push(null)

    def op_LOAD_GLOBAL(self, inst: dis.Instruction):
        assert isinstance(inst.arg, int)  # for typing
        load_null = inst.arg & 1
        op = Op(opname="global", bc_inst=inst)
        op.add_input("env", self.effect)
        null = op.add_output("null")
        if load_null:
            self.push(null)
        self.push(op.add_output(f"{inst.argval}"))

    def op_LOAD_CONST(self, inst: dis.Instruction):
        op = Op(opname="const", bc_inst=inst)
        self.push(op.add_output("out"))

    def op_STORE_FAST(self, inst: dis.Instruction):
        tos = self.pop()
        op = Op(opname="store", bc_inst=inst)
        op.add_input("value", tos)
        varname = self._decorate_varname(inst.argval)
        self.store(varname, op.add_output(varname))

    def op_LOAD_FAST(self, inst: dis.Instruction):
        varname = self._decorate_varname(inst.argval)
        self.push(self.load(varname))

    def op_LOAD_ATTR(self, inst: dis.Instruction):
        obj = self.pop()
        attr = inst.argval
        op = Op(opname=f"load_attr.{attr}", bc_inst=inst)
        op.add_input("obj", obj)
        self.push(op.add_output("out"))

    def op_LOAD_METHOD(self, inst: dis.Instruction):
        obj = self.pop()
        attr = inst.argval
        op = Op(opname=f"load_method.{attr}", bc_inst=inst)
        op.add_input("obj", obj)
        self.push(op.add_output("null"))
        self.push(op.add_output("out"))

    def op_LOAD_DEREF(self, inst: dis.Instruction):
        op = Op(opname="load_deref", bc_inst=inst)
        self.push(op.add_output("out"))

    def op_PRECALL(self, inst: dis.Instruction):
        pass  # no-op

    def op_KW_NAMES(self, inst: dis.Instruction):
        op = Op(opname="kw_names", bc_inst=inst)
        self.set_kw_names(op.add_output("out"))

    def op_CALL(self, inst: dis.Instruction):
        argc: int = inst.argval
        arg1plus = reversed([self.pop() for _ in range(argc)])
        arg0 = self.pop()  # TODO
        kw_names = self.pop_kw_names()

        args: list[ValueState] = [arg0, *arg1plus]
        callable = self.pop()  # TODO
        opname = "call" if kw_names is None else "call.kw"
        op = Op(opname=opname, bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("callee", callable)
        for i, arg in enumerate(args):
            op.add_input(f"arg.{i}", arg)

        if kw_names is not None:
            op.add_input("kw_names", kw_names)

        self.replace_effect(op.add_output("env", is_effect=True))
        self.push(op.add_output("ret"))

    def op_GET_ITER(self, inst: dis.Instruction):
        tos = self.pop()
        op = Op(opname="getiter", bc_inst=inst)
        op.add_input("obj", tos)
        self.push(op.add_output("iter"))

    def op_FOR_ITER(self, inst: dis.Instruction):
        tos = self.top()
        op = Op(opname="foriter", bc_inst=inst)
        op.add_input("iter", tos)
        # Store the indvar into an internal variable
        self.store("indvar", op.add_output("indvar"))

    def _binaryop(self, opname: str, inst: dis.Instruction):
        rhs = self.pop()
        lhs = self.pop()
        op = Op(opname=opname, bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("lhs", lhs)
        op.add_input("rhs", rhs)
        self.replace_effect(op.add_output("env", is_effect=True))
        self.push(op.add_output("out"))

    def op_BINARY_OP(self, inst: dis.Instruction):
        self._binaryop("binaryop", inst)

    def op_COMPARE_OP(self, inst: dis.Instruction):
        self._binaryop("compareop", inst)

    def op_IS_OP(self, inst: dis.Instruction):
        self._binaryop("is_op", inst)

    def _unaryop(self, opname: str, inst: dis.Instruction):
        op = Op(opname=opname, bc_inst=inst)
        op.add_input("val", self.pop())
        self.push(op.add_output("out"))

    def op_UNARY_NOT(self, inst: dis.Instruction):
        self._unaryop("not", inst)

    def op_BINARY_SUBSCR(self, inst: dis.Instruction):
        index = self.pop()
        target = self.pop()
        op = Op(opname="binary_subscr", bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("index", index)
        op.add_input("target", target)
        self.replace_effect(op.add_output("env", is_effect=True))
        self.push(op.add_output("out"))

    def op_STORE_SUBSCR(self, inst: dis.Instruction):
        index = self.pop()
        target = self.pop()
        value = self.pop()
        op = Op(opname="store_subscr", bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("index", index)
        op.add_input("target", target)
        op.add_input("value", value)
        self.replace_effect(op.add_output("env", is_effect=True))

    def op_BUILD_TUPLE(self, inst: dis.Instruction):
        count = inst.arg
        assert isinstance(count, int)
        items = list(reversed([self.pop() for _ in range(count)]))
        op = Op(opname="build_tuple", bc_inst=inst)
        for i, it in enumerate(items):
            op.add_input(str(i), it)
        self.push(op.add_output("out"))

    def op_BUILD_SLICE(self, inst: dis.Instruction):
        argc = inst.arg
        if argc == 2:
            tos = self.pop()
            tos1 = self.pop()
            start = tos1
            stop = tos
            step = None
        elif argc == 3:
            tos = self.pop()
            tos1 = self.pop()
            tos2 = self.pop()
            start = tos2
            stop = tos1
            step = tos
        else:
            raise Exception("unreachable")

        op = Op(opname="build_slice", bc_inst=inst)
        op.add_input("start", start)
        op.add_input("stop", stop)
        if step is not None:
            op.add_input("step", step)
        self.push(op.add_output("out"))

    def op_RETURN_VALUE(self, inst: dis.Instruction):
        tos = self.pop()
        op = Op(opname="ret", bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("retval", tos)
        self.replace_effect(op.add_output("env", is_effect=True))

    def op_RAISE_VARARGS(self, inst: dis.Instruction):
        if inst.arg == 0:
            exc = None
            # # No re-raising within a try-except block.
            # # But we allow bare reraise.
            # if state.has_active_try():
            #     raise UnsupportedError(
            #         "The re-raising of an exception is not yet supported.",
            #         loc=self.get_debug_loc(inst.lineno),
            #     )
            raise NotImplementedError
        elif inst.arg == 1:
            exc = self.pop()
        else:
            raise ValueError("Multiple argument raise is not supported.")
        op = Op(opname="raise_varargs", bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("exc", exc)
        self.replace_effect(op.add_output("env", is_effect=True))

    def op_JUMP_FORWARD(self, inst: dis.Instruction):
        pass  # no-op

    def op_JUMP_BACKWARD(self, inst: dis.Instruction):
        pass  # no-op

    def _POP_JUMP_X_IF_Y(self, inst: dis.Instruction, *, opname: str):
        tos = self.pop()
        op = Op(opname, bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("pred", tos)
        self.replace_effect(op.add_output("env", is_effect=True))

    def op_POP_JUMP_FORWARD_IF_TRUE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_true")

    def op_POP_JUMP_FORWARD_IF_FALSE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_false")

    def op_POP_JUMP_BACKWARD_IF_TRUE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_true")

    def op_POP_JUMP_BACKWARD_IF_FALSE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_false")

    def op_POP_JUMP_FORWARD_IF_NONE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_none")

    def op_POP_JUMP_FORWARD_IF_NOT_NONE(self, inst: dis.Instruction):
        self._POP_JUMP_X_IF_Y(inst, opname="jump.if_not_none")

    def _JUMP_IF_X_OR_POP(self, inst: dis.Instruction, *, opname):
        tos = self.top()
        op = Op(opname, bc_inst=inst)
        op.add_input("env", self.effect)
        op.add_input("pred", tos)
        self.replace_effect(op.add_output("env", is_effect=True))

    def op_JUMP_IF_TRUE_OR_POP(self, inst: dis.Instruction):
        self._JUMP_IF_X_OR_POP(inst, opname="jump.if_true")

    def op_JUMP_IF_FALSE_OR_POP(self, inst: dis.Instruction):
        self._JUMP_IF_X_OR_POP(inst, opname="jump.if_false")
