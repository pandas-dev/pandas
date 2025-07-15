"""
Implement transformation on Numba IR
"""


from collections import namedtuple, defaultdict
import logging
import operator

from numba.core.analysis import compute_cfg_from_blocks, find_top_level_loops
from numba.core import errors, ir, ir_utils
from numba.core.analysis import compute_use_defs, compute_cfg_from_blocks
from numba.core.utils import PYVERSION


_logger = logging.getLogger(__name__)


def _extract_loop_lifting_candidates(cfg, blocks):
    """
    Returns a list of loops that are candidate for loop lifting
    """
    # check well-formed-ness of the loop
    def same_exit_point(loop):
        "all exits must point to the same location"
        outedges = set()
        for k in loop.exits:
            succs = set(x for x, _ in cfg.successors(k))
            if not succs:
                # If the exit point has no successor, it contains an return
                # statement, which is not handled by the looplifting code.
                # Thus, this loop is not a candidate.
                _logger.debug("return-statement in loop.")
                return False
            outedges |= succs
        ok = len(outedges) == 1
        _logger.debug("same_exit_point=%s (%s)", ok, outedges)
        return ok

    def one_entry(loop):
        "there is one entry"
        ok = len(loop.entries) == 1
        _logger.debug("one_entry=%s", ok)
        return ok

    def cannot_yield(loop):
        "cannot have yield inside the loop"
        insiders = set(loop.body) | set(loop.entries) | set(loop.exits)
        for blk in map(blocks.__getitem__, insiders):
            for inst in blk.body:
                if isinstance(inst, ir.Assign):
                    if isinstance(inst.value, ir.Yield):
                        _logger.debug("has yield")
                        return False
        _logger.debug("no yield")
        return True

    _logger.info('finding looplift candidates')
    # the check for cfg.entry_point in the loop.entries is to prevent a bad
    # rewrite where a prelude for a lifted loop would get written into block -1
    # if a loop entry were in block 0
    candidates = []
    for loop in find_top_level_loops(cfg):
        _logger.debug("top-level loop: %s", loop)
        if (same_exit_point(loop) and one_entry(loop) and cannot_yield(loop) and
            cfg.entry_point() not in loop.entries):
            candidates.append(loop)
            _logger.debug("add candidate: %s", loop)
    return candidates


def find_region_inout_vars(blocks, livemap, callfrom, returnto, body_block_ids):
    """Find input and output variables to a block region.
    """
    inputs = livemap[callfrom]
    outputs = livemap[returnto]

    # ensure live variables are actually used in the blocks, else remove,
    # saves having to create something valid to run through postproc
    # to achieve similar
    loopblocks = {}
    for k in body_block_ids:
        loopblocks[k] = blocks[k]

    used_vars = set()
    def_vars = set()
    defs = compute_use_defs(loopblocks)
    for vs in defs.usemap.values():
        used_vars |= vs
    for vs in defs.defmap.values():
        def_vars |= vs
    used_or_defined = used_vars | def_vars

    # note: sorted for stable ordering
    inputs = sorted(set(inputs) & used_or_defined)
    outputs = sorted(set(outputs) & used_or_defined & def_vars)
    return inputs, outputs


_loop_lift_info = namedtuple('loop_lift_info',
                             'loop,inputs,outputs,callfrom,returnto')


def _loop_lift_get_candidate_infos(cfg, blocks, livemap):
    """
    Returns information on looplifting candidates.
    """
    loops = _extract_loop_lifting_candidates(cfg, blocks)
    loopinfos = []
    for loop in loops:

        [callfrom] = loop.entries   # requirement checked earlier
        an_exit = next(iter(loop.exits))  # anyone of the exit block
        if len(loop.exits) > 1:
            # has multiple exits
            [(returnto, _)] = cfg.successors(an_exit)  # requirement checked earlier
        else:
            # does not have multiple exits
            returnto = an_exit

        local_block_ids = set(loop.body) | set(loop.entries) | set(loop.exits)
        inputs, outputs = find_region_inout_vars(
            blocks=blocks,
            livemap=livemap,
            callfrom=callfrom,
            returnto=returnto,
            body_block_ids=local_block_ids,
        )

        lli = _loop_lift_info(loop=loop, inputs=inputs, outputs=outputs,
                              callfrom=callfrom, returnto=returnto)
        loopinfos.append(lli)

    return loopinfos


def _loop_lift_modify_call_block(liftedloop, block, inputs, outputs, returnto):
    """
    Transform calling block from top-level function to call the lifted loop.
    """
    scope = block.scope
    loc = block.loc
    blk = ir.Block(scope=scope, loc=loc)

    ir_utils.fill_block_with_call(
        newblock=blk,
        callee=liftedloop,
        label_next=returnto,
        inputs=inputs,
        outputs=outputs,
    )
    return blk


def _loop_lift_prepare_loop_func(loopinfo, blocks):
    """
    Inplace transform loop blocks for use as lifted loop.
    """
    entry_block = blocks[loopinfo.callfrom]
    scope = entry_block.scope
    loc = entry_block.loc

    # Lowering assumes the first block to be the one with the smallest offset
    firstblk = min(blocks) - 1
    blocks[firstblk] = ir_utils.fill_callee_prologue(
        block=ir.Block(scope=scope, loc=loc),
        inputs=loopinfo.inputs,
        label_next=loopinfo.callfrom,
    )
    blocks[loopinfo.returnto] = ir_utils.fill_callee_epilogue(
        block=ir.Block(scope=scope, loc=loc),
        outputs=loopinfo.outputs,
    )


def _loop_lift_modify_blocks(func_ir, loopinfo, blocks,
                             typingctx, targetctx, flags, locals):
    """
    Modify the block inplace to call to the lifted-loop.
    Returns a dictionary of blocks of the lifted-loop.
    """
    from numba.core.dispatcher import LiftedLoop

    # Copy loop blocks
    loop = loopinfo.loop

    loopblockkeys = set(loop.body) | set(loop.entries)
    if len(loop.exits) > 1:
        # has multiple exits
        loopblockkeys |= loop.exits
    loopblocks = dict((k, blocks[k].copy()) for k in loopblockkeys)
    # Modify the loop blocks
    _loop_lift_prepare_loop_func(loopinfo, loopblocks)
    # Since Python 3.13, [END_FOR, POP_TOP] sequence becomes the start of the
    # block causing the block to have line number of the start of previous loop.
    # Fix this using the loc of the first getiter.
    getiter_exprs = []
    for blk in loopblocks.values():
        getiter_exprs.extend(blk.find_exprs(op="getiter"))
    first_getiter = min(getiter_exprs, key=lambda x: x.loc.line)
    loop_loc = first_getiter.loc
    # Create a new IR for the lifted loop
    lifted_ir = func_ir.derive(blocks=loopblocks,
                               arg_names=tuple(loopinfo.inputs),
                               arg_count=len(loopinfo.inputs),
                               force_non_generator=True,
                               loc=loop_loc)
    liftedloop = LiftedLoop(lifted_ir,
                            typingctx, targetctx, flags, locals)

    # modify for calling into liftedloop
    callblock = _loop_lift_modify_call_block(liftedloop, blocks[loopinfo.callfrom],
                                             loopinfo.inputs, loopinfo.outputs,
                                             loopinfo.returnto)
    # remove blocks
    for k in loopblockkeys:
        del blocks[k]
    # update main interpreter callsite into the liftedloop
    blocks[loopinfo.callfrom] = callblock
    return liftedloop


def _has_multiple_loop_exits(cfg, lpinfo):
    """Returns True if there is more than one exit in the loop.

    NOTE: "common exits" refers to the situation where a loop exit has another
    loop exit as its successor. In that case, we do not need to alter it.
    """
    if len(lpinfo.exits) <= 1:
        return False
    exits = set(lpinfo.exits)
    pdom = cfg.post_dominators()

    # Eliminate blocks that have other blocks as post-dominators.
    processed = set()
    remain = set(exits) # create a copy to work on
    while remain:
        node = remain.pop()
        processed.add(node)
        exits -= pdom[node] - {node}
        remain = exits - processed

    return len(exits) > 1


def _pre_looplift_transform(func_ir):
    """Canonicalize loops for looplifting.
    """
    from numba.core.postproc import PostProcessor

    cfg = compute_cfg_from_blocks(func_ir.blocks)
    # For every loop that has multiple exits, combine the exits into one.
    for loop_info in cfg.loops().values():
        if _has_multiple_loop_exits(cfg, loop_info):
            func_ir, _common_key = _fix_multi_exit_blocks(
                func_ir, loop_info.exits
            )
    # Reset and reprocess the func_ir
    func_ir._reset_analysis_variables()
    PostProcessor(func_ir).run()
    return func_ir


def loop_lifting(func_ir, typingctx, targetctx, flags, locals):
    """
    Loop lifting transformation.

    Given a interpreter `func_ir` returns a 2 tuple of
    `(toplevel_interp, [loop0_interp, loop1_interp, ....])`
    """
    func_ir = _pre_looplift_transform(func_ir)
    blocks = func_ir.blocks.copy()
    cfg = compute_cfg_from_blocks(blocks)
    loopinfos = _loop_lift_get_candidate_infos(cfg, blocks,
                                               func_ir.variable_lifetime.livemap)
    loops = []
    if loopinfos:
        _logger.debug('loop lifting this IR with %d candidates:\n%s',
                      len(loopinfos), func_ir.dump_to_string())
    for loopinfo in loopinfos:
        lifted = _loop_lift_modify_blocks(func_ir, loopinfo, blocks,
                                          typingctx, targetctx, flags, locals)
        loops.append(lifted)

    # Make main IR
    main = func_ir.derive(blocks=blocks)

    return main, loops


def canonicalize_cfg_single_backedge(blocks):
    """
    Rewrite loops that have multiple backedges.
    """
    cfg = compute_cfg_from_blocks(blocks)
    newblocks = blocks.copy()

    def new_block_id():
        return max(newblocks.keys()) + 1

    def has_multiple_backedges(loop):
        count = 0
        for k in loop.body:
            blk = blocks[k]
            edges = blk.terminator.get_targets()
            # is a backedge?
            if loop.header in edges:
                count += 1
                if count > 1:
                    # early exit
                    return True
        return False

    def yield_loops_with_multiple_backedges():
        for lp in cfg.loops().values():
            if has_multiple_backedges(lp):
                yield lp

    def replace_target(term, src, dst):
        def replace(target):
            return (dst if target == src else target)

        if isinstance(term, ir.Branch):
            return ir.Branch(cond=term.cond,
                             truebr=replace(term.truebr),
                             falsebr=replace(term.falsebr),
                             loc=term.loc)
        elif isinstance(term, ir.Jump):
            return ir.Jump(target=replace(term.target), loc=term.loc)
        else:
            assert not term.get_targets()
            return term

    def rewrite_single_backedge(loop):
        """
        Add new tail block that gathers all the backedges
        """
        header = loop.header
        tailkey = new_block_id()
        for blkkey in loop.body:
            blk = newblocks[blkkey]
            if header in blk.terminator.get_targets():
                newblk = blk.copy()
                # rewrite backedge into jumps to new tail block
                newblk.body[-1] = replace_target(blk.terminator, header,
                                                 tailkey)
                newblocks[blkkey] = newblk
        # create new tail block
        entryblk = newblocks[header]
        tailblk = ir.Block(scope=entryblk.scope, loc=entryblk.loc)
        # add backedge
        tailblk.append(ir.Jump(target=header, loc=tailblk.loc))
        newblocks[tailkey] = tailblk

    for loop in yield_loops_with_multiple_backedges():
        rewrite_single_backedge(loop)

    return newblocks


def canonicalize_cfg(blocks):
    """
    Rewrite the given blocks to canonicalize the CFG.
    Returns a new dictionary of blocks.
    """
    return canonicalize_cfg_single_backedge(blocks)


def with_lifting(func_ir, typingctx, targetctx, flags, locals):
    """With-lifting transformation

    Rewrite the IR to extract all withs.
    Only the top-level withs are extracted.
    Returns the (the_new_ir, the_lifted_with_ir)
    """
    from numba.core import postproc

    def dispatcher_factory(func_ir, objectmode=False, **kwargs):
        from numba.core.dispatcher import LiftedWith, ObjModeLiftedWith

        myflags = flags.copy()
        if objectmode:
            # Lifted with-block cannot looplift
            myflags.enable_looplift = False
            # Lifted with-block uses object mode
            myflags.enable_pyobject = True
            myflags.force_pyobject = True
            myflags.no_cpython_wrapper = False
            cls = ObjModeLiftedWith
        else:
            cls = LiftedWith
        return cls(func_ir, typingctx, targetctx, myflags, locals, **kwargs)

    # find where with-contexts regions are
    withs, func_ir = find_setupwiths(func_ir)

    if not withs:
        return func_ir, []

    postproc.PostProcessor(func_ir).run()  # ensure we have variable lifetime
    assert func_ir.variable_lifetime
    vlt = func_ir.variable_lifetime
    blocks = func_ir.blocks.copy()
    cfg = vlt.cfg
    # For each with-regions, mutate them according to
    # the kind of contextmanager
    sub_irs = []
    for (blk_start, blk_end) in withs:
        body_blocks = []
        for node in _cfg_nodes_in_region(cfg, blk_start, blk_end):
            body_blocks.append(node)
        _legalize_with_head(blocks[blk_start])
        # Find the contextmanager
        cmkind, extra = _get_with_contextmanager(func_ir, blocks, blk_start)
        # Mutate the body and get new IR
        sub = cmkind.mutate_with_body(func_ir, blocks, blk_start, blk_end,
                                      body_blocks, dispatcher_factory,
                                      extra)
        sub_irs.append(sub)
    if not sub_irs:
        # Unchanged
        new_ir = func_ir
    else:
        new_ir = func_ir.derive(blocks)
    return new_ir, sub_irs


def _get_with_contextmanager(func_ir, blocks, blk_start):
    """Get the global object used for the context manager
    """
    _illegal_cm_msg = "Illegal use of context-manager."

    def get_var_dfn(var):
        """Get the definition given a variable"""
        return func_ir.get_definition(var)

    def get_ctxmgr_obj(var_ref):
        """Return the context-manager object and extra info.

        The extra contains the arguments if the context-manager is used
        as a call.
        """
        # If the contextmanager used as a Call
        dfn = func_ir.get_definition(var_ref)
        if isinstance(dfn, ir.Expr) and dfn.op == 'call':
            args = [get_var_dfn(x) for x in dfn.args]
            kws = {k: get_var_dfn(v) for k, v in dfn.kws}
            extra = {'args': args, 'kwargs': kws}
            var_ref = dfn.func
        else:
            extra = None

        ctxobj = ir_utils.guard(ir_utils.find_outer_value, func_ir, var_ref)

        # check the contextmanager object
        if ctxobj is ir.UNDEFINED:
            raise errors.CompilerError(
                "Undefined variable used as context manager",
                loc=blocks[blk_start].loc,
                )

        if ctxobj is None:
            raise errors.CompilerError(_illegal_cm_msg, loc=dfn.loc)

        return ctxobj, extra

    # Scan the start of the with-region for the contextmanager
    for stmt in blocks[blk_start].body:
        if isinstance(stmt, ir.EnterWith):
            var_ref = stmt.contextmanager
            ctxobj, extra = get_ctxmgr_obj(var_ref)
            if not hasattr(ctxobj, 'mutate_with_body'):
                raise errors.CompilerError(
                    "Unsupported context manager in use",
                    loc=blocks[blk_start].loc,
                    )
            return ctxobj, extra
    # No contextmanager found?
    raise errors.CompilerError(
        "malformed with-context usage",
        loc=blocks[blk_start].loc,
        )


def _legalize_with_head(blk):
    """Given *blk*, the head block of the with-context, check that it doesn't
    do anything else.
    """
    counters = defaultdict(int)
    for stmt in blk.body:
        counters[type(stmt)] += 1
    if counters.pop(ir.EnterWith) != 1:
        raise errors.CompilerError(
            "with's head-block must have exactly 1 ENTER_WITH",
            loc=blk.loc,
            )
    if counters.pop(ir.Jump, 0) != 1:
        raise errors.CompilerError(
            "with's head-block must have exactly 1 JUMP",
            loc=blk.loc,
            )
    # Can have any number of del
    counters.pop(ir.Del, None)
    # There MUST NOT be any other statements
    if counters:
        raise errors.CompilerError(
            "illegal statements in with's head-block",
            loc=blk.loc,
            )


def _cfg_nodes_in_region(cfg, region_begin, region_end):
    """Find the set of CFG nodes that are in the given region
    """
    region_nodes = set()
    stack = [region_begin]
    while stack:
        tos = stack.pop()
        succlist = list(cfg.successors(tos))
        # a single block function will have a empty successor list
        if succlist:
            succs, _ = zip(*succlist)
            nodes = set([node for node in succs
                        if node not in region_nodes and
                        node != region_end])
            stack.extend(nodes)
            region_nodes |= nodes

    return region_nodes


def find_setupwiths(func_ir):
    """Find all top-level with.

    Returns a list of ranges for the with-regions.
    """
    def find_ranges(blocks):

        cfg = compute_cfg_from_blocks(blocks)
        sus_setups, sus_pops = set(), set()
        # traverse the cfg and collect all suspected SETUP_WITH and POP_BLOCK
        # statements so that we can iterate over them
        for label, block in blocks.items():
            for stmt in block.body:
                if ir_utils.is_setup_with(stmt):
                    sus_setups.add(label)
                if ir_utils.is_pop_block(stmt):
                    sus_pops.add(label)

        # now that we do have the statements, iterate through them in reverse
        # topo order and from each start looking for pop_blocks
        setup_with_to_pop_blocks_map = defaultdict(set)
        for setup_block in cfg.topo_sort(sus_setups, reverse=True):
            # begin pop_block, search
            to_visit, seen = [], []
            to_visit.append(setup_block)
            while to_visit:
                # get whatever is next and record that we have seen it
                block = to_visit.pop()
                seen.append(block)
                # go through the body of the block, looking for statements
                for stmt in blocks[block].body:
                    # raise detected before pop_block
                    if ir_utils.is_raise(stmt):
                            raise errors.CompilerError(
                                'unsupported control flow due to raise '
                                'statements inside with block'
                                )
                    # if a pop_block, process it
                    if ir_utils.is_pop_block(stmt) and block in sus_pops:
                        # record the jump target of this block belonging to this setup
                        setup_with_to_pop_blocks_map[setup_block].add(block)
                        # remove the block from blocks to be matched
                        sus_pops.remove(block)
                        # stop looking, we have reached the frontier
                        break
                    # if we are still here, by the block terminator,
                    # add all its targets to the to_visit stack, unless we
                    # have seen them already
                    if ir_utils.is_terminator(stmt):
                        for t in stmt.get_targets():
                            if t not in seen:
                                to_visit.append(t)

        return setup_with_to_pop_blocks_map

    blocks = func_ir.blocks
    # initial find, will return a dictionary, mapping indices of blocks
    # containing SETUP_WITH statements to a set of indices of blocks containing
    # POP_BLOCK statements
    with_ranges_dict = find_ranges(blocks)
    # rewrite the CFG in case there are multiple POP_BLOCK statements for one
    # with
    func_ir = consolidate_multi_exit_withs(with_ranges_dict, blocks, func_ir)
    # here we need to turn the withs back into a list of tuples so that the
    # rest of the code can cope
    with_ranges_tuple = [(s, list(p)[0])
             for (s, p) in with_ranges_dict.items()]

    # check for POP_BLOCKS with multiple outgoing edges and reject
    for (_, p) in with_ranges_tuple:
        targets = blocks[p].terminator.get_targets()
        if len(targets) != 1:
            raise errors.CompilerError(
                "unsupported control flow: with-context contains branches "
                "(i.e. break/return/raise) that can leave the block "
            )
    # now we check for returns inside with and reject them
    for (_, p) in with_ranges_tuple:
        target_block = blocks[p]
        if ir_utils.is_return(func_ir.blocks[
                target_block.terminator.get_targets()[0]].terminator):
            _rewrite_return(func_ir, p)

    # now we need to rewrite the tuple such that we have SETUP_WITH matching the
    # successor of the block that contains the POP_BLOCK.
    with_ranges_tuple = [(s, func_ir.blocks[p].terminator.get_targets()[0])
                         for (s, p) in with_ranges_tuple]

    # finally we check for nested with statements and reject them
    with_ranges_tuple = _eliminate_nested_withs(with_ranges_tuple)

    return with_ranges_tuple, func_ir


def _rewrite_return(func_ir, target_block_label):
    """Rewrite a return block inside a with statement.

    Arguments
    ---------

    func_ir: Function IR
      the CFG to transform
    target_block_label: int
      the block index/label of the block containing the POP_BLOCK statement


    This implements a CFG transformation to insert a block between two other
    blocks.

    The input situation is:

    ┌───────────────┐
    │   top         │
    │   POP_BLOCK   │
    │   bottom      │
    └───────┬───────┘
            │
    ┌───────▼───────┐
    │               │
    │    RETURN     │
    │               │
    └───────────────┘

    If such a pattern is detected in IR, it means there is a `return` statement
    within a `with` context. The basic idea is to rewrite the CFG as follows:

    ┌───────────────┐
    │   top         │
    │   POP_BLOCK   │
    │               │
    └───────┬───────┘
            │
    ┌───────▼───────┐
    │               │
    │     bottom    │
    │               │
    └───────┬───────┘
            │
    ┌───────▼───────┐
    │               │
    │    RETURN     │
    │               │
    └───────────────┘

    We split the block that contains the `POP_BLOCK` statement into two blocks.
    Everything from the beginning of the block up to and including the
    `POP_BLOCK` statement is considered the 'top' and everything below is
    considered 'bottom'. Finally the jump statements are re-wired to make sure
    the CFG remains valid.

    """
    # the block itself from the index
    target_block = func_ir.blocks[target_block_label]
    # get the index of the block containing the return
    target_block_successor_label = target_block.terminator.get_targets()[0]
    # the return block
    target_block_successor = func_ir.blocks[target_block_successor_label]

    # create the new return block with an appropriate label
    max_label = ir_utils.find_max_label(func_ir.blocks)
    new_label = max_label + 1
    # create the new return block
    new_block_loc = target_block_successor.loc
    new_block_scope = ir.Scope(None, loc=new_block_loc)
    new_block = ir.Block(new_block_scope, loc=new_block_loc)

    # Split the block containing the POP_BLOCK into top and bottom
    # Block must be of the form:
    # -----------------
    # <some stmts>
    # POP_BLOCK
    # <some more stmts>
    # JUMP
    # -----------------
    top_body, bottom_body = [], []
    pop_blocks = [*target_block.find_insts(ir.PopBlock)]
    assert len(pop_blocks) == 1
    assert len([*target_block.find_insts(ir.Jump)]) == 1
    assert isinstance(target_block.body[-1], ir.Jump)
    pb_marker = pop_blocks[0]
    pb_is = target_block.body.index(pb_marker)
    top_body.extend(target_block.body[:pb_is])
    top_body.append(ir.Jump(target_block_successor_label, target_block.loc))
    bottom_body.extend(target_block.body[pb_is:-1])
    bottom_body.append(ir.Jump(new_label, target_block.loc))

    # get the contents of the return block
    return_body = func_ir.blocks[target_block_successor_label].body
    # finally, re-assign all blocks
    new_block.body.extend(return_body)
    target_block_successor.body.clear()
    target_block_successor.body.extend(bottom_body)
    target_block.body.clear()
    target_block.body.extend(top_body)

    # finally, append the new return block and rebuild the IR properties
    func_ir.blocks[new_label] = new_block
    func_ir._definitions = ir_utils.build_definitions(func_ir.blocks)
    return func_ir


def _eliminate_nested_withs(with_ranges):
    known_ranges = []
    def within_known_range(start, end, known_ranges):
        for a, b in known_ranges:
            # FIXME: this should be a comparison in topological order, right
            # now we are comparing the integers of the blocks, stuff probably
            # works by accident.
            if start > a and end < b:
                return True
        return False

    for s, e in sorted(with_ranges):
        if not within_known_range(s, e, known_ranges):
            known_ranges.append((s, e))

    return known_ranges

def consolidate_multi_exit_withs(withs: dict, blocks, func_ir):
    """Modify the FunctionIR to merge the exit blocks of with constructs.
    """
    for k in withs:
        vs : set = withs[k]
        if len(vs) > 1:
            func_ir, common = _fix_multi_exit_blocks(
                func_ir, vs, split_condition=ir_utils.is_pop_block,
            )
            withs[k] = {common}
    return func_ir


def _fix_multi_exit_blocks(func_ir, exit_nodes, *, split_condition=None):
    """Modify the FunctionIR to create a single common exit node given the
    original exit nodes.

    Parameters
    ----------
    func_ir :
        The FunctionIR. Mutated inplace.
    exit_nodes :
        The original exit nodes. A sequence of block keys.
    split_condition : callable or None
        If not None, it is a callable with the signature
        `split_condition(statement)` that determines if the `statement` is the
        splitting point (e.g. `POP_BLOCK`) in an exit node.
        If it's None, the exit node is not split.
    """

    # Convert the following:
    #
    #     |           |
    # +-------+   +-------+
    # | exit0 |   | exit1 |
    # +-------+   +-------+
    #     |           |
    # +-------+   +-------+
    # | after0|   | after1|
    # +-------+   +-------+
    #     |           |
    #
    # To roughly:
    #
    #     |           |
    # +-------+   +-------+
    # | exit0 |   | exit1 |
    # +-------+   +-------+
    #     |           |
    #     +-----+-----+
    #           |
    #      +---------+
    #      | common  |
    #      +---------+
    #           |
    #       +-------+
    #       | post  |
    #       +-------+
    #           |
    #     +-----+-----+
    #     |           |
    # +-------+   +-------+
    # | after0|   | after1|
    # +-------+   +-------+

    blocks = func_ir.blocks
    # Getting the scope
    any_blk = min(func_ir.blocks.values())
    scope = any_blk.scope
    # Getting the maximum block label
    max_label = max(func_ir.blocks) + 1
    # Define the new common block for the new exit.
    common_block = ir.Block(any_blk.scope, loc=ir.unknown_loc)
    common_label = max_label
    max_label += 1
    blocks[common_label] = common_block
    # Define the new block after the exit.
    post_block = ir.Block(any_blk.scope, loc=ir.unknown_loc)
    post_label = max_label
    max_label += 1
    blocks[post_label] = post_block

    # Adjust each exit node
    remainings = []
    for i, k in enumerate(exit_nodes):
        blk = blocks[k]

        # split the block if needed
        if split_condition is not None:
            for pt, stmt in enumerate(blk.body):
                if split_condition(stmt):
                    break
        else:
            # no splitting
            pt = -1

        before = blk.body[:pt]
        after = blk.body[pt:]
        remainings.append(after)

        # Add control-point variable to mark which exit block this is.
        blk.body = before
        loc = blk.loc
        blk.body.append(
            ir.Assign(value=ir.Const(i, loc=loc),
                      target=scope.get_or_define("$cp", loc=loc),
                      loc=loc)
        )
        # Replace terminator with a jump to the common block
        assert not blk.is_terminated
        blk.body.append(ir.Jump(common_label, loc=ir.unknown_loc))

    if split_condition is not None:
        # Move the splitting statement to the common block
        common_block.body.append(remainings[0][0])
    assert not common_block.is_terminated
    # Append jump from common block to post block
    common_block.body.append(ir.Jump(post_label, loc=loc))

    # Make if-else tree to jump to target
    remain_blocks = []
    for remain in remainings:
        remain_blocks.append(max_label)
        max_label += 1

    switch_block = post_block
    loc = ir.unknown_loc
    for i, remain in enumerate(remainings):
        match_expr = scope.redefine("$cp_check", loc=loc)
        match_rhs = scope.redefine("$cp_rhs", loc=loc)

        # Do comparison to match control-point variable to the exit block
        switch_block.body.append(
            ir.Assign(
                value=ir.Const(i, loc=loc),
                target=match_rhs,
                loc=loc
            ),
        )

        # Add assignment for the comparison
        switch_block.body.append(
            ir.Assign(
                value=ir.Expr.binop(
                    fn=operator.eq, lhs=scope.get("$cp"), rhs=match_rhs,
                    loc=loc,
                ),
                target=match_expr,
                loc=loc
            ),
        )

        # Insert jump to the next case
        [jump_target] = remain[-1].get_targets()
        switch_block.body.append(
            ir.Branch(match_expr, jump_target, remain_blocks[i], loc=loc),
        )
        switch_block = ir.Block(scope=scope, loc=loc)
        blocks[remain_blocks[i]] = switch_block

    # Add the final jump
    switch_block.body.append(ir.Jump(jump_target, loc=loc))

    return func_ir, common_label
