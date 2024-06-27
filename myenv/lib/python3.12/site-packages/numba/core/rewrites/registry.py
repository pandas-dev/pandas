from collections import defaultdict

from numba.core import config


class Rewrite(object):
    '''Defines the abstract base class for Numba rewrites.
    '''

    def __init__(self, state=None):
        '''Constructor for the Rewrite class.
        '''
        pass

    def match(self, func_ir, block, typemap, calltypes):
        '''Overload this method to check an IR block for matching terms in the
        rewrite.
        '''
        return False

    def apply(self):
        '''Overload this method to return a rewritten IR basic block when a
        match has been found.
        '''
        raise NotImplementedError("Abstract Rewrite.apply() called!")


class RewriteRegistry(object):
    '''Defines a registry for Numba rewrites.
    '''
    _kinds = frozenset(['before-inference', 'after-inference'])

    def __init__(self):
        '''Constructor for the rewrite registry.  Initializes the rewrites
        member to an empty list.
        '''
        self.rewrites = defaultdict(list)

    def register(self, kind):
        """
        Decorator adding a subclass of Rewrite to the registry for
        the given *kind*.
        """
        if kind not in self._kinds:
            raise KeyError("invalid kind %r" % (kind,))
        def do_register(rewrite_cls):
            if not issubclass(rewrite_cls, Rewrite):
                raise TypeError('{0} is not a subclass of Rewrite'.format(
                    rewrite_cls))
            self.rewrites[kind].append(rewrite_cls)
            return rewrite_cls
        return do_register

    def apply(self, kind, state):
        '''Given a pipeline and a dictionary of basic blocks, exhaustively
        attempt to apply all registered rewrites to all basic blocks.
        '''
        assert kind in self._kinds
        blocks = state.func_ir.blocks
        old_blocks = blocks.copy()
        for rewrite_cls in self.rewrites[kind]:
            # Exhaustively apply a rewrite until it stops matching.
            rewrite = rewrite_cls(state)
            work_list = list(blocks.items())
            while work_list:
                key, block = work_list.pop()
                matches = rewrite.match(state.func_ir, block, state.typemap,
                                        state.calltypes)
                if matches:
                    if config.DEBUG or config.DUMP_IR:
                        print("_" * 70)
                        print("REWRITING (%s):" % rewrite_cls.__name__)
                        block.dump()
                        print("_" * 60)
                    new_block = rewrite.apply()
                    blocks[key] = new_block
                    work_list.append((key, new_block))
                    if config.DEBUG or config.DUMP_IR:
                        new_block.dump()
                        print("_" * 70)
        # If any blocks were changed, perform a sanity check.
        for key, block in blocks.items():
            if block != old_blocks[key]:
                block.verify()

        # Some passes, e.g. _inline_const_arraycall are known to occasionally
        # do invalid things WRT ir.Del, others, e.g. RewriteArrayExprs do valid
        # things with ir.Del, but the placement is not optimal. The lines below
        # fix-up the IR so that ref counts are valid and optimally placed,
        # see #4093 for context. This has to be run here opposed to in
        # apply() as the CFG needs computing so full IR is needed.
        from numba.core import postproc
        post_proc = postproc.PostProcessor(state.func_ir)
        post_proc.run()


rewrite_registry = RewriteRegistry()
register_rewrite = rewrite_registry.register
