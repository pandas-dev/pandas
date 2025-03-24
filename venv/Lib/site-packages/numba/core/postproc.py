from functools import cached_property
from numba.core import ir, analysis, transforms, ir_utils


class YieldPoint(object):

    def __init__(self, block, inst):
        assert isinstance(block, ir.Block)
        assert isinstance(inst, ir.Yield)
        self.block = block
        self.inst = inst
        self.live_vars = None
        self.weak_live_vars = None


class GeneratorInfo(object):

    def __init__(self):
        # { index: YieldPoint }
        self.yield_points = {}
        # Ordered list of variable names
        self.state_vars = []

    def get_yield_points(self):
        """
        Return an iterable of YieldPoint instances.
        """
        return self.yield_points.values()


class VariableLifetime(object):
    """
    For lazily building information of variable lifetime
    """
    def __init__(self, blocks):
        self._blocks = blocks

    @cached_property
    def cfg(self):
        return analysis.compute_cfg_from_blocks(self._blocks)

    @cached_property
    def usedefs(self):
        return analysis.compute_use_defs(self._blocks)

    @cached_property
    def livemap(self):
        return analysis.compute_live_map(self.cfg, self._blocks,
                                         self.usedefs.usemap,
                                         self.usedefs.defmap)

    @cached_property
    def deadmaps(self):
        return analysis.compute_dead_maps(self.cfg, self._blocks, self.livemap,
                                          self.usedefs.defmap)


# other packages that define new nodes add calls for inserting dels
# format: {type:function}
ir_extension_insert_dels = {}


class PostProcessor(object):
    """
    A post-processor for Numba IR.
    """

    def __init__(self, func_ir):
        self.func_ir = func_ir

    def run(self, emit_dels: bool = False, extend_lifetimes: bool = False):
        """
        Run the following passes over Numba IR:
        - canonicalize the CFG
        - emit explicit `del` instructions for variables
        - compute lifetime of variables
        - compute generator info (if function is a generator function)
        """
        self.func_ir.blocks = transforms.canonicalize_cfg(self.func_ir.blocks)
        vlt = VariableLifetime(self.func_ir.blocks)
        self.func_ir.variable_lifetime = vlt

        bev = analysis.compute_live_variables(vlt.cfg, self.func_ir.blocks,
                                              vlt.usedefs.defmap,
                                              vlt.deadmaps.combined)
        for offset, ir_block in self.func_ir.blocks.items():
            self.func_ir.block_entry_vars[ir_block] = bev[offset]

        if self.func_ir.is_generator:
            self.func_ir.generator_info = GeneratorInfo()
            self._compute_generator_info()
        else:
            self.func_ir.generator_info = None

        # Emit del nodes, do this last as the generator info parsing generates
        # and then strips dels as part of its analysis.
        if emit_dels:
            self._insert_var_dels(extend_lifetimes=extend_lifetimes)

    def _populate_generator_info(self):
        """
        Fill `index` for the Yield instruction and create YieldPoints.
        """
        dct = self.func_ir.generator_info.yield_points
        assert not dct, 'rerunning _populate_generator_info'
        for block in self.func_ir.blocks.values():
            for inst in block.body:
                if isinstance(inst, ir.Assign):
                    yieldinst = inst.value
                    if isinstance(yieldinst, ir.Yield):
                        index = len(dct) + 1
                        yieldinst.index = index
                        yp = YieldPoint(block, yieldinst)
                        dct[yieldinst.index] = yp

    def _compute_generator_info(self):
        """
        Compute the generator's state variables as the union of live variables
        at all yield points.
        """
        # generate del info, it's used in analysis here, strip it out at the end
        self._insert_var_dels()
        self._populate_generator_info()
        gi = self.func_ir.generator_info
        for yp in gi.get_yield_points():
            live_vars = set(self.func_ir.get_block_entry_vars(yp.block))
            weak_live_vars = set()
            stmts = iter(yp.block.body)
            for stmt in stmts:
                if isinstance(stmt, ir.Assign):
                    if stmt.value is yp.inst:
                        break
                    live_vars.add(stmt.target.name)
                elif isinstance(stmt, ir.Del):
                    live_vars.remove(stmt.value)
            else:
                assert 0, "couldn't find yield point"
            # Try to optimize out any live vars that are deleted immediately
            # after the yield point.
            for stmt in stmts:
                if isinstance(stmt, ir.Del):
                    name = stmt.value
                    if name in live_vars:
                        live_vars.remove(name)
                        weak_live_vars.add(name)
                else:
                    break
            yp.live_vars = live_vars
            yp.weak_live_vars = weak_live_vars

        st = set()
        for yp in gi.get_yield_points():
            st |= yp.live_vars
            st |= yp.weak_live_vars
        gi.state_vars = sorted(st)
        self.remove_dels()

    def _insert_var_dels(self, extend_lifetimes=False):
        """
        Insert del statements for each variable.
        Returns a 2-tuple of (variable definition map, variable deletion map)
        which indicates variables defined and deleted in each block.

        The algorithm avoids relying on explicit knowledge on loops and
        distinguish between variables that are defined locally vs variables that
        come from incoming blocks.
        We start with simple usage (variable reference) and definition (variable
        creation) maps on each block. Propagate the liveness info to predecessor
        blocks until it stabilize, at which point we know which variables must
        exist before entering each block. Then, we compute the end of variable
        lives and insert del statements accordingly. Variables are deleted after
        the last use. Variable referenced by terminators (e.g. conditional
        branch and return) are deleted by the successors or the caller.
        """
        vlt = self.func_ir.variable_lifetime
        self._patch_var_dels(vlt.deadmaps.internal, vlt.deadmaps.escaping,
                             extend_lifetimes=extend_lifetimes)

    def _patch_var_dels(self, internal_dead_map, escaping_dead_map,
                        extend_lifetimes=False):
        """
        Insert delete in each block
        """
        for offset, ir_block in self.func_ir.blocks.items():
            # for each internal var, insert delete after the last use
            internal_dead_set = internal_dead_map[offset].copy()
            delete_pts = []
            # for each statement in reverse order
            for stmt in reversed(ir_block.body[:-1]):
                # internal vars that are used here
                live_set = set(v.name for v in stmt.list_vars())
                dead_set = live_set & internal_dead_set
                for T, def_func in ir_extension_insert_dels.items():
                    if isinstance(stmt, T):
                        done_dels = def_func(stmt, dead_set)
                        dead_set -= done_dels
                        internal_dead_set -= done_dels
                # used here but not afterwards
                delete_pts.append((stmt, dead_set))
                internal_dead_set -= dead_set

            # rewrite body and insert dels
            body = []
            lastloc = ir_block.loc
            del_store = []
            for stmt, delete_set in reversed(delete_pts):
                # If using extended lifetimes then the Dels are all put at the
                # block end just ahead of the terminator, so associate their
                # location with the terminator.
                if extend_lifetimes:
                    lastloc = ir_block.body[-1].loc
                else:
                    lastloc = stmt.loc
                # Ignore dels (assuming no user inserted deletes)
                if not isinstance(stmt, ir.Del):
                    body.append(stmt)
                # note: the reverse sort is not necessary for correctness
                #       it is just to minimize changes to test for now
                for var_name in sorted(delete_set, reverse=True):
                    delnode = ir.Del(var_name, loc=lastloc)
                    if extend_lifetimes:
                        del_store.append(delnode)
                    else:
                        body.append(delnode)
            if extend_lifetimes:
                body.extend(del_store)
            body.append(ir_block.body[-1])  # terminator
            ir_block.body = body

            # vars to delete at the start
            escape_dead_set = escaping_dead_map[offset]
            for var_name in sorted(escape_dead_set):
                ir_block.prepend(ir.Del(var_name, loc=ir_block.body[0].loc))

    def remove_dels(self):
        """
        Strips the IR of Del nodes
        """
        ir_utils.remove_dels(self.func_ir.blocks)
