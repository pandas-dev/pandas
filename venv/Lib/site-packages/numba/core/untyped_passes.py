from collections import defaultdict, namedtuple
from contextlib import contextmanager
from copy import deepcopy, copy
import warnings

from numba.core.compiler_machinery import (FunctionPass, AnalysisPass,
                                           SSACompliantMixin, register_pass)
from numba.core import (errors, types, ir, bytecode, postproc, rewrites, config,
                        transforms, consts)
from numba.misc.special import literal_unroll
from numba.core.analysis import (dead_branch_prune, rewrite_semantic_constants,
                                 find_literally_calls, compute_cfg_from_blocks,
                                 compute_use_defs)
from numba.core.ir_utils import (guard, resolve_func_from_module, simplify_CFG,
                                 GuardException, convert_code_obj_to_function,
                                 build_definitions,
                                 replace_var_names, get_name_var_table,
                                 compile_to_numba_ir, get_definition,
                                 find_max_label, rename_labels,
                                 transfer_scope, fixup_var_define_in_scope,
                                 )
from numba.core.ssa import reconstruct_ssa
from numba.core import interpreter


@contextmanager
def fallback_context(state, msg):
    """
    Wraps code that would signal a fallback to object mode
    """
    try:
        yield
    except Exception as e:
        if not state.status.can_fallback:
            raise
        else:
            # Clear all references attached to the traceback
            e = e.with_traceback(None)
            # this emits a warning containing the error message body in the
            # case of fallback from npm to objmode
            loop_lift = '' if state.flags.enable_looplift else 'OUT'
            msg_rewrite = ("\nCompilation is falling back to object mode "
                           "WITH%s looplifting enabled because %s"
                           % (loop_lift, msg))
            warnings.warn_explicit('%s due to: %s' % (msg_rewrite, e),
                                   errors.NumbaWarning,
                                   state.func_id.filename,
                                   state.func_id.firstlineno)
            raise


@register_pass(mutates_CFG=True, analysis_only=False)
class ExtractByteCode(FunctionPass):
    _name = "extract_bytecode"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        Extract bytecode from function
        """
        func_id = state['func_id']
        bc = bytecode.ByteCode(func_id)
        if config.DUMP_BYTECODE:
            print(bc.dump())

        state['bc'] = bc
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class TranslateByteCode(FunctionPass):
    _name = "translate_bytecode"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        Analyze bytecode and translating to Numba IR
        """
        func_id = state['func_id']
        bc = state['bc']
        interp = interpreter.Interpreter(func_id)
        func_ir = interp.interpret(bc)
        state["func_ir"] = func_ir
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class FixupArgs(FunctionPass):
    _name = "fixup_args"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        state['nargs'] = state['func_ir'].arg_count
        if not state['args'] and state['flags'].force_pyobject:
            # Allow an empty argument types specification when object mode
            # is explicitly requested.
            state['args'] = (types.pyobject,) * state['nargs']
        elif len(state['args']) != state['nargs']:
            raise TypeError("Signature mismatch: %d argument types given, "
                            "but function takes %d arguments"
                            % (len(state['args']), state['nargs']))
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class IRProcessing(FunctionPass):
    _name = "ir_processing"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        func_ir = state['func_ir']
        post_proc = postproc.PostProcessor(func_ir)
        post_proc.run()

        if config.DEBUG or config.DUMP_IR:
            name = func_ir.func_id.func_qualname
            print(("IR DUMP: %s" % name).center(80, "-"))
            func_ir.dump()
            if func_ir.is_generator:
                print(("GENERATOR INFO: %s" % name).center(80, "-"))
                func_ir.dump_generator_info()
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class RewriteSemanticConstants(FunctionPass):
    _name = "rewrite_semantic_constants"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        This prunes dead branches, a dead branch is one which is derivable as
        not taken at compile time purely based on const/literal evaluation.
        """
        assert state.func_ir
        msg = ('Internal error in pre-inference dead branch pruning '
               'pass encountered during compilation of '
               'function "%s"' % (state.func_id.func_name,))
        with fallback_context(state, msg):
            rewrite_semantic_constants(state.func_ir, state.args)

        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class DeadBranchPrune(SSACompliantMixin, FunctionPass):
    _name = "dead_branch_prune"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        This prunes dead branches, a dead branch is one which is derivable as
        not taken at compile time purely based on const/literal evaluation.
        """

        # purely for demonstration purposes, obtain the analysis from a pass
        # declare as a required dependent
        semantic_const_analysis = self.get_analysis(type(self))  # noqa

        assert state.func_ir
        msg = ('Internal error in pre-inference dead branch pruning '
               'pass encountered during compilation of '
               'function "%s"' % (state.func_id.func_name,))
        with fallback_context(state, msg):
            dead_branch_prune(state.func_ir, state.args)

        return True

    def get_analysis_usage(self, AU):
        AU.add_required(RewriteSemanticConstants)


@register_pass(mutates_CFG=True, analysis_only=False)
class InlineClosureLikes(FunctionPass):
    _name = "inline_closure_likes"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        # Ensure we have an IR and type information.
        assert state.func_ir

        # if the return type is a pyobject, there's no type info available and
        # no ability to resolve certain typed function calls in the array
        # inlining code, use this variable to indicate
        typed_pass = not isinstance(state.return_type, types.misc.PyObject)
        from numba.core.inline_closurecall import InlineClosureCallPass
        inline_pass = InlineClosureCallPass(
            state.func_ir,
            state.flags.auto_parallel,
            state.parfor_diagnostics.replaced_fns,
            typed_pass)
        inline_pass.run()

        # Remove all Dels, and re-run postproc
        post_proc = postproc.PostProcessor(state.func_ir)
        post_proc.run()

        fixup_var_define_in_scope(state.func_ir.blocks)

        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class GenericRewrites(FunctionPass):
    _name = "generic_rewrites"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        Perform any intermediate representation rewrites before type
        inference.
        """
        assert state.func_ir
        msg = ('Internal error in pre-inference rewriting '
               'pass encountered during compilation of '
               'function "%s"' % (state.func_id.func_name,))
        with fallback_context(state, msg):
            rewrites.rewrite_registry.apply('before-inference', state)
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class WithLifting(FunctionPass):
    _name = "with_lifting"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """
        Extract with-contexts
        """
        main, withs = transforms.with_lifting(
            func_ir=state.func_ir,
            typingctx=state.typingctx,
            targetctx=state.targetctx,
            flags=state.flags,
            locals=state.locals,
        )
        if withs:
            from numba.core.compiler import compile_ir, _EarlyPipelineCompletion
            cres = compile_ir(state.typingctx, state.targetctx, main,
                              state.args, state.return_type,
                              state.flags, state.locals,
                              lifted=tuple(withs), lifted_from=None,
                              pipeline_class=type(state.pipeline))
            raise _EarlyPipelineCompletion(cres)
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class InlineInlinables(FunctionPass):
    """
    This pass will inline a function wrapped by the numba.jit decorator directly
    into the site of its call depending on the value set in the 'inline' kwarg
    to the decorator.

    This is an untyped pass. CFG simplification is performed at the end of the
    pass but no block level clean up is performed on the mutated IR (typing
    information is not available to do so).
    """
    _name = "inline_inlinables"
    _DEBUG = False

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        """Run inlining of inlinables
        """
        if self._DEBUG:
            print('before inline'.center(80, '-'))
            print(state.func_ir.dump())
            print(''.center(80, '-'))

        from numba.core.inline_closurecall import (InlineWorker,
                                                   callee_ir_validator)
        inline_worker = InlineWorker(state.typingctx,
                                     state.targetctx,
                                     state.locals,
                                     state.pipeline,
                                     state.flags,
                                     validator=callee_ir_validator)

        modified = False
        # use a work list, look for call sites via `ir.Expr.op == call` and
        # then pass these to `self._do_work` to make decisions about inlining.
        work_list = list(state.func_ir.blocks.items())
        while work_list:
            label, block = work_list.pop()
            for i, instr in enumerate(block.body):
                if isinstance(instr, ir.Assign):
                    expr = instr.value
                    if isinstance(expr, ir.Expr) and expr.op == 'call':
                        if guard(self._do_work, state, work_list, block, i,
                                 expr, inline_worker):
                            modified = True
                            break  # because block structure changed

        if modified:
            # clean up unconditional branches that appear due to inlined
            # functions introducing blocks
            cfg = compute_cfg_from_blocks(state.func_ir.blocks)
            for dead in cfg.dead_nodes():
                del state.func_ir.blocks[dead]
            post_proc = postproc.PostProcessor(state.func_ir)
            post_proc.run()
            state.func_ir.blocks = simplify_CFG(state.func_ir.blocks)

        if self._DEBUG:
            print('after inline'.center(80, '-'))
            print(state.func_ir.dump())
            print(''.center(80, '-'))
        return True

    def _do_work(self, state, work_list, block, i, expr, inline_worker):
        from numba.core.compiler import run_frontend
        from numba.core.cpu import InlineOptions

        # try and get a definition for the call, this isn't always possible as
        # it might be a eval(str)/part generated awaiting update etc. (parfors)
        to_inline = None
        try:
            to_inline = state.func_ir.get_definition(expr.func)
        except Exception:
            if self._DEBUG:
                print("Cannot find definition for %s" % expr.func)
            return False
        # do not handle closure inlining here, another pass deals with that.
        if getattr(to_inline, 'op', False) == 'make_function':
            return False

        # see if the definition is a "getattr", in which case walk the IR to
        # try and find the python function via the module from which it's
        # imported, this should all be encoded in the IR.
        if getattr(to_inline, 'op', False) == 'getattr':
            val = resolve_func_from_module(state.func_ir, to_inline)
        else:
            # This is likely a freevar or global
            #
            # NOTE: getattr 'value' on a call may fail if it's an ir.Expr as
            # getattr is overloaded to look in _kws.
            try:
                val = getattr(to_inline, 'value', False)
            except Exception:
                raise GuardException

        # if something was found...
        if val:
            # check it's dispatcher-like, the targetoptions attr holds the
            # kwargs supplied in the jit decorator and is where 'inline' will
            # be if it is present.
            topt = getattr(val, 'targetoptions', False)
            if topt:
                inline_type = topt.get('inline', None)
                # has 'inline' been specified?
                if inline_type is not None:
                    inline_opt = InlineOptions(inline_type)
                    # Could this be inlinable?
                    if not inline_opt.is_never_inline:
                        # yes, it could be inlinable
                        do_inline = True
                        pyfunc = val.py_func
                        # Has it got an associated cost model?
                        if inline_opt.has_cost_model:
                            # yes, it has a cost model, use it to determine
                            # whether to do the inline
                            py_func_ir = run_frontend(pyfunc)
                            do_inline = inline_type(expr, state.func_ir,
                                                    py_func_ir)
                        # if do_inline is True then inline!
                        if do_inline:
                            _, _, _, new_blocks = \
                                inline_worker.inline_function(state.func_ir,
                                                              block, i, pyfunc,)
                            if work_list is not None:
                                for blk in new_blocks:
                                    work_list.append(blk)
                            return True
        return False


@register_pass(mutates_CFG=False, analysis_only=False)
class PreserveIR(AnalysisPass):
    """
    Preserves the IR in the metadata
    """

    _name = "preserve_ir"

    def __init__(self):
        AnalysisPass.__init__(self)

    def run_pass(self, state):
        state.metadata['preserved_ir'] = state.func_ir.copy()
        return False


@register_pass(mutates_CFG=False, analysis_only=True)
class FindLiterallyCalls(FunctionPass):
    """Find calls to `numba.literally()` and signal if its requirement is not
    satisfied.
    """
    _name = "find_literally"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        find_literally_calls(state.func_ir, state.args)
        return False


@register_pass(mutates_CFG=True, analysis_only=False)
class CanonicalizeLoopExit(FunctionPass):
    """A pass to canonicalize loop exit by splitting it from function exit.
    """
    _name = "canonicalize_loop_exit"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        fir = state.func_ir
        cfg = compute_cfg_from_blocks(fir.blocks)
        status = False
        for loop in cfg.loops().values():
            for exit_label in loop.exits:
                if exit_label in cfg.exit_points():
                    self._split_exit_block(fir, cfg, exit_label)
                    status = True

        fir._reset_analysis_variables()

        vlt = postproc.VariableLifetime(fir.blocks)
        fir.variable_lifetime = vlt
        return status

    def _split_exit_block(self, fir, cfg, exit_label):
        curblock = fir.blocks[exit_label]
        newlabel = exit_label + 1
        newlabel = find_max_label(fir.blocks) + 1
        fir.blocks[newlabel] = curblock
        newblock = ir.Block(scope=curblock.scope, loc=curblock.loc)
        newblock.append(ir.Jump(newlabel, loc=curblock.loc))
        fir.blocks[exit_label] = newblock
        # Rename all labels
        fir.blocks = rename_labels(fir.blocks)


@register_pass(mutates_CFG=True, analysis_only=False)
class CanonicalizeLoopEntry(FunctionPass):
    """A pass to canonicalize loop header by splitting it from function entry.

    This is needed for loop-lifting; esp in py3.8
    """
    _name = "canonicalize_loop_entry"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        fir = state.func_ir
        cfg = compute_cfg_from_blocks(fir.blocks)
        status = False
        for loop in cfg.loops().values():
            if len(loop.entries) == 1:
                [entry_label] = loop.entries
                if entry_label == cfg.entry_point():
                    self._split_entry_block(fir, cfg, loop, entry_label)
                    status = True
        fir._reset_analysis_variables()

        vlt = postproc.VariableLifetime(fir.blocks)
        fir.variable_lifetime = vlt
        return status

    def _split_entry_block(self, fir, cfg, loop, entry_label):
        # Find iterator inputs into the for-loop header
        header_block = fir.blocks[loop.header]
        deps = set()
        for expr in header_block.find_exprs(op="iternext"):
            deps.add(expr.value)
        # Find the getiter for each iterator
        entry_block = fir.blocks[entry_label]

        # Find the start of loop entry statement that needs to be included.
        startpt = None
        list_of_insts = list(entry_block.find_insts(ir.Assign))
        for assign in reversed(list_of_insts):
            if assign.target in deps:
                rhs = assign.value
                if isinstance(rhs, ir.Var):
                    if rhs.is_temp:
                        deps.add(rhs)
                elif isinstance(rhs, ir.Expr):
                    expr = rhs
                    if expr.op == 'getiter':
                        startpt = assign
                        if expr.value.is_temp:
                            deps.add(expr.value)
                    elif expr.op == 'call':
                        defn = guard(get_definition, fir, expr.func)
                        if isinstance(defn, ir.Global):
                            if expr.func.is_temp:
                                deps.add(expr.func)
                elif isinstance(rhs, ir.Global) and rhs.value is range:
                    startpt = assign

        if startpt is None:
            return

        splitpt = entry_block.body.index(startpt)
        new_block = entry_block.copy()
        new_block.body = new_block.body[splitpt:]
        new_block.loc = new_block.body[0].loc
        new_label = find_max_label(fir.blocks) + 1
        entry_block.body = entry_block.body[:splitpt]
        entry_block.append(ir.Jump(new_label, loc=new_block.loc))

        fir.blocks[new_label] = new_block
        # Rename all labels
        fir.blocks = rename_labels(fir.blocks)


@register_pass(mutates_CFG=False, analysis_only=True)
class PrintIRCFG(FunctionPass):
    _name = "print_ir_cfg"

    def __init__(self):
        FunctionPass.__init__(self)
        self._ver = 0

    def run_pass(self, state):
        fir = state.func_ir
        self._ver += 1
        fir.render_dot(filename_prefix='v{}'.format(self._ver)).render()
        return False


@register_pass(mutates_CFG=True, analysis_only=False)
class MakeFunctionToJitFunction(FunctionPass):
    """
    This swaps an ir.Expr.op == "make_function" i.e. a closure, for a compiled
    function containing the closure body and puts it in ir.Global. It's a 1:1
    statement value swap. `make_function` is already untyped
    """
    _name = "make_function_op_code_to_jit_function"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        from numba import njit
        func_ir = state.func_ir
        mutated = False
        for idx, blk in func_ir.blocks.items():
            for stmt in blk.body:
                if isinstance(stmt, ir.Assign):
                    if isinstance(stmt.value, ir.Expr):
                        if stmt.value.op == "make_function":
                            node = stmt.value
                            getdef = func_ir.get_definition
                            kw_default = getdef(node.defaults)
                            ok = False
                            if (kw_default is None or
                                    isinstance(kw_default, ir.Const)):
                                ok = True
                            elif isinstance(kw_default, tuple):
                                ok = all([isinstance(getdef(x), ir.Const)
                                          for x in kw_default])
                            elif isinstance(kw_default, ir.Expr):
                                if kw_default.op != "build_tuple":
                                    continue
                                ok = all([isinstance(getdef(x), ir.Const)
                                          for x in kw_default.items])
                            if not ok:
                                continue

                            pyfunc = convert_code_obj_to_function(node, func_ir)
                            func = njit()(pyfunc)
                            new_node = ir.Global(node.code.co_name, func,
                                                 stmt.loc)
                            stmt.value = new_node
                            mutated |= True

        # if a change was made the del ordering is probably wrong, patch up
        if mutated:
            post_proc = postproc.PostProcessor(func_ir)
            post_proc.run()

        return mutated


@register_pass(mutates_CFG=True, analysis_only=False)
class TransformLiteralUnrollConstListToTuple(FunctionPass):
    """ This pass spots a `literal_unroll([<constant values>])` and rewrites it
    as a `literal_unroll(tuple(<constant values>))`.
    """
    _name = "transform_literal_unroll_const_list_to_tuple"

    _accepted_types = (types.BaseTuple, types.LiteralList)

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        mutated = False
        func_ir = state.func_ir
        for label, blk in func_ir.blocks.items():
            calls = [_ for _ in blk.find_exprs('call')]
            for call in calls:
                glbl = guard(get_definition, func_ir, call.func)
                if glbl and isinstance(glbl, (ir.Global, ir.FreeVar)):
                    # find a literal_unroll
                    if glbl.value is literal_unroll:
                        if len(call.args) > 1:
                            msg = "literal_unroll takes one argument, found %s"
                            raise errors.UnsupportedError(msg % len(call.args),
                                                          call.loc)
                        # get the arg, make sure its a build_list
                        unroll_var = call.args[0]
                        to_unroll = guard(get_definition, func_ir, unroll_var)
                        if (isinstance(to_unroll, ir.Expr) and
                                to_unroll.op == "build_list"):
                            # make sure they are all const items in the list
                            for i, item in enumerate(to_unroll.items):
                                val = guard(get_definition, func_ir, item)
                                if not val:
                                    msg = ("multiple definitions for variable "
                                           "%s, cannot resolve constant")
                                    raise errors.UnsupportedError(msg % item,
                                                                  to_unroll.loc)
                                if not isinstance(val, ir.Const):
                                    msg = ("Found non-constant value at "
                                           "position %s in a list argument to "
                                           "literal_unroll" % i)
                                    raise errors.UnsupportedError(msg,
                                                                  to_unroll.loc)
                            # The above appears ok, now swap the build_list for
                            # a built tuple.

                            # find the assignment for the unroll target
                            to_unroll_lhs = guard(get_definition, func_ir,
                                                  unroll_var, lhs_only=True)

                            if to_unroll_lhs is None:
                                msg = ("multiple definitions for variable "
                                       "%s, cannot resolve constant")
                                raise errors.UnsupportedError(msg % unroll_var,
                                                              to_unroll.loc)
                            # scan all blocks looking for the LHS
                            for b in func_ir.blocks.values():
                                asgn = b.find_variable_assignment(
                                    to_unroll_lhs.name)
                                if asgn is not None:
                                    break
                            else:
                                msg = ("Cannot find assignment for known "
                                       "variable %s") % to_unroll_lhs.name
                                raise errors.CompilerError(msg, to_unroll.loc)

                            # Create a tuple with the list items as contents
                            tup = ir.Expr.build_tuple(to_unroll.items,
                                                      to_unroll.loc)

                            # swap the list for the tuple
                            asgn.value = tup
                            mutated = True
                        elif (isinstance(to_unroll, ir.Expr) and
                              to_unroll.op == "build_tuple"):
                            # this is fine, do nothing
                            pass
                        elif (isinstance(to_unroll, (ir.Global, ir.FreeVar)) and
                              isinstance(to_unroll.value, tuple)):
                            # this is fine, do nothing
                            pass
                        elif isinstance(to_unroll, ir.Arg):
                            # this is only fine if the arg is a tuple
                            ty = state.typemap[to_unroll.name]
                            if not isinstance(ty, self._accepted_types):
                                msg = ("Invalid use of literal_unroll with a "
                                       "function argument, only tuples are "
                                       "supported as function arguments, found "
                                       "%s") % ty
                                raise errors.UnsupportedError(msg,
                                                              to_unroll.loc)
                        else:
                            extra = None
                            if isinstance(to_unroll, ir.Expr):
                                # probably a slice
                                if to_unroll.op == "getitem":
                                    ty = state.typemap[to_unroll.value.name]
                                    # check if this is a tuple slice
                                    if not isinstance(ty, self._accepted_types):
                                        extra = "operation %s" % to_unroll.op
                                        loc = to_unroll.loc
                            elif isinstance(to_unroll, ir.Arg):
                                extra = "non-const argument %s" % to_unroll.name
                                loc = to_unroll.loc
                            else:
                                if to_unroll is None:
                                    extra = ('multiple definitions of '
                                             'variable "%s".' % unroll_var.name)
                                    loc = unroll_var.loc
                                else:
                                    loc = to_unroll.loc
                                    extra = "unknown problem"

                            if extra:
                                msg = ("Invalid use of literal_unroll, "
                                       "argument should be a tuple or a list "
                                       "of constant values. Failure reason: "
                                       "found %s" % extra)
                                raise errors.UnsupportedError(msg, loc)
        return mutated


@register_pass(mutates_CFG=True, analysis_only=False)
class MixedContainerUnroller(FunctionPass):
    _name = "mixed_container_unroller"

    _DEBUG = False

    _accepted_types = (types.BaseTuple, types.LiteralList)

    def __init__(self):
        FunctionPass.__init__(self)

    def analyse_tuple(self, tup):
        """
        Returns a map of type->list(indexes) for a typed tuple
        """
        d = defaultdict(list)
        for i, ty in enumerate(tup):
            d[ty].append(i)
        return d

    def add_offset_to_labels_w_ignore(self, blocks, offset, ignore=None):
        """add an offset to all block labels and jump/branch targets
        don't add an offset to anything in the ignore list
        """
        if ignore is None:
            ignore = set()

        new_blocks = {}
        for l, b in blocks.items():
            # some parfor last blocks might be empty
            term = None
            if b.body:
                term = b.body[-1]
            if isinstance(term, ir.Jump):
                if term.target not in ignore:
                    b.body[-1] = ir.Jump(term.target + offset, term.loc)
            if isinstance(term, ir.Branch):
                if term.truebr not in ignore:
                    new_true = term.truebr + offset
                else:
                    new_true = term.truebr

                if term.falsebr not in ignore:
                    new_false = term.falsebr + offset
                else:
                    new_false = term.falsebr
                b.body[-1] = ir.Branch(term.cond, new_true, new_false, term.loc)
            new_blocks[l + offset] = b
        return new_blocks

    def inject_loop_body(self, switch_ir, loop_ir, caller_max_label,
                         dont_replace, switch_data):
        """
        Injects the "loop body" held in `loop_ir` into `switch_ir` where ever
        there is a statement of the form `SENTINEL.<int> = RHS`. It also:
        * Finds and then deliberately does not relabel non-local jumps so as to
          make the switch table suitable for injection into the IR from which
          the loop body was derived.
        * Looks for `typed_getitem` and wires them up to loop body version
          specific variables or, if possible, directly writes in their constant
          value at their use site.

        Args:
        - switch_ir, the switch table with SENTINELS as generated by
          self.gen_switch
        - loop_ir, the IR of the loop blocks (derived from the original func_ir)
        - caller_max_label, the maximum label in the func_ir caller
        - dont_replace, variables that should not be renamed (to handle
          references to variables that are incoming at the loop head/escaping at
          the loop exit.
        - switch_data, the switch table data used to generated the switch_ir,
          can be generated by self.analyse_tuple.

        Returns:
        - A type specific switch table with each case containing a versioned
          loop body suitable for injection as a replacement for the loop_ir.
        """

        # Switch IR came from code gen, immediately relabel to prevent
        # collisions with IR derived from the user code (caller)
        switch_ir.blocks = self.add_offset_to_labels_w_ignore(
            switch_ir.blocks, caller_max_label + 1)

        # Find the sentinels and validate the form
        sentinel_exits = set()
        sentinel_blocks = []
        for lbl, blk in switch_ir.blocks.items():
            for i, stmt in enumerate(blk.body):
                if isinstance(stmt, ir.Assign):
                    if "SENTINEL" in stmt.target.name:
                        sentinel_blocks.append(lbl)
                        sentinel_exits.add(blk.body[-1].target)
                        break

        assert len(sentinel_exits) == 1  # should only be 1 exit
        switch_ir.blocks.pop(sentinel_exits.pop())  # kill the exit, it's dead

        # find jumps that are non-local, we won't relabel these
        ignore_set = set()
        local_lbl = [x for x in loop_ir.blocks.keys()]
        for lbl, blk in loop_ir.blocks.items():
            for i, stmt in enumerate(blk.body):
                if isinstance(stmt, ir.Jump):
                    if stmt.target not in local_lbl:
                        ignore_set.add(stmt.target)
                if isinstance(stmt, ir.Branch):
                    if stmt.truebr not in local_lbl:
                        ignore_set.add(stmt.truebr)
                    if stmt.falsebr not in local_lbl:
                        ignore_set.add(stmt.falsebr)

        # make sure the generated switch table matches the switch data
        assert len(sentinel_blocks) == len(switch_data)

        # replace the sentinel_blocks with the loop body
        for lbl, branch_ty in zip(sentinel_blocks, switch_data.keys()):
            loop_blocks = deepcopy(loop_ir.blocks)
            # relabel blocks WRT switch table, each block replacement will shift
            # the maximum label
            max_label = max(switch_ir.blocks.keys())
            loop_blocks = self.add_offset_to_labels_w_ignore(
                loop_blocks, max_label + 1, ignore_set)

            # start label
            loop_start_lbl = min(loop_blocks.keys())

            # fix the typed_getitem locations in the loop blocks
            for blk in loop_blocks.values():
                new_body = []
                for stmt in blk.body:
                    if isinstance(stmt, ir.Assign):
                        if (isinstance(stmt.value, ir.Expr) and
                                stmt.value.op == "typed_getitem"):
                            if isinstance(branch_ty, types.Literal):
                                scope = switch_ir.blocks[lbl].scope
                                new_const_name = scope.redefine(
                                    "branch_const", stmt.loc).name
                                new_const_var = ir.Var(
                                    blk.scope, new_const_name, stmt.loc)
                                new_const_val = ir.Const(
                                    branch_ty.literal_value, stmt.loc)
                                const_assign = ir.Assign(
                                    new_const_val, new_const_var, stmt.loc)
                                new_assign = ir.Assign(
                                    new_const_var, stmt.target, stmt.loc)
                                new_body.append(const_assign)
                                new_body.append(new_assign)
                                dont_replace.append(new_const_name)
                            else:
                                orig = stmt.value
                                new_typed_getitem = ir.Expr.typed_getitem(
                                    value=orig.value, dtype=branch_ty,
                                    index=orig.index, loc=orig.loc)
                                new_assign = ir.Assign(
                                    new_typed_getitem, stmt.target, stmt.loc)
                                new_body.append(new_assign)
                        else:
                            new_body.append(stmt)
                    else:
                        new_body.append(stmt)
                blk.body = new_body

            # rename
            var_table = get_name_var_table(loop_blocks)
            drop_keys = []
            for k, v in var_table.items():
                if v.name in dont_replace:
                    drop_keys.append(k)
            for k in drop_keys:
                var_table.pop(k)

            new_var_dict = {}
            for name, var in var_table.items():
                scope = switch_ir.blocks[lbl].scope
                try:
                    scope.get_exact(name)
                except errors.NotDefinedError:
                    # In case the scope doesn't have the variable, we need to
                    # define it prior creating new copies of it! This is
                    # because the scope of the function and the scope of the
                    # loop are different and the variable needs to be redefined
                    # within the scope of the loop.
                    scope.define(name, var.loc)
                new_var_dict[name] = scope.redefine(name, var.loc).name
            replace_var_names(loop_blocks, new_var_dict)

            # clobber the sentinel body and then stuff in the rest
            switch_ir.blocks[lbl] = deepcopy(loop_blocks[loop_start_lbl])
            remaining_keys = [y for y in loop_blocks.keys()]
            remaining_keys.remove(loop_start_lbl)
            for k in remaining_keys:
                switch_ir.blocks[k] = deepcopy(loop_blocks[k])

        if self._DEBUG:
            print("-" * 80 + "EXIT STUFFER")
            switch_ir.dump()
            print("-" * 80)

        return switch_ir

    def gen_switch(self, data, index):
        """
        Generates a function with a switch table like
        def foo():
            if PLACEHOLDER_INDEX in (<integers>):
                SENTINEL = None
            elif PLACEHOLDER_INDEX in (<integers>):
                SENTINEL = None
            ...
            else:
                raise RuntimeError

        The data is a map of (type : indexes) for example:
        (int64, int64, float64)
        might give:
        {int64: [0, 1], float64: [2]}

        The index is the index variable for the driving range loop over the
        mixed tuple.
        """
        elif_tplt = "\n\telif PLACEHOLDER_INDEX in (%s,):\n\t\tSENTINEL = None"

        # Note regarding the insertion of the garbage/defeat variables below:
        # These values have been designed and inserted to defeat a specific
        # behaviour of the cpython optimizer. The optimization was introduced
        # in Python 3.10.

        # The URL for the BPO is:
        # https://bugs.python.org/issue44626
        # The code for the optimization can be found at:
        # https://github.com/python/cpython/blob/d41abe8/Python/compile.c#L7533-L7557

        # Essentially the CPython optimizer will inline the exit block under
        # certain circumstances and thus replace the jump with a return if the
        # exit block is small enough.  This is an issue for unroller, as it
        # looks for a jump, not a return, when it inserts the generated switch
        # table.

        # Part of the condition for this optimization to be applied is that the
        # exit block not exceed a certain (4 at the time of writing) number of
        # bytecode instructions. We defeat the optimizer by inserting a
        # sufficient number of instructions so that the exit block is big
        # enough. We don't care about this garbage, because the generated exit
        # block is discarded anyway when we smash the switch table into the
        # original function and so all the inserted garbage is dropped again.

        # The final lines of the stacktrace w/o this will look like:
        #
        #  File "/numba/numba/core/untyped_passes.py", line 830, \
        #       in inject_loop_body
        #   sentinel_exits.add(blk.body[-1].target)
        # AttributeError: Failed in nopython mode pipeline \
        #       (step: handles literal_unroll)
        # Failed in literal_unroll_subpipeline mode pipeline \
        #       (step: performs mixed container unroll)
        # 'Return' object has no attribute 'target'
        #
        # Which indicates that a Return has been found instead of a Jump

        b = ('def foo():\n\tif PLACEHOLDER_INDEX in (%s,):\n\t\t'
             'SENTINEL = None\n%s\n\telse:\n\t\t'
             'raise RuntimeError("Unreachable")\n\t'
             'py310_defeat1 = 1\n\t'
             'py310_defeat2 = 2\n\t'
             'py310_defeat3 = 3\n\t'
             'py310_defeat4 = 4\n\t'
             )
        keys = [k for k in data.keys()]

        elifs = []
        for i in range(1, len(keys)):
            elifs.append(elif_tplt % ','.join(map(str, data[keys[i]])))
        src = b % (','.join(map(str, data[keys[0]])), ''.join(elifs))
        wstr = src
        l = {}
        exec(wstr, {}, l)
        bfunc = l['foo']
        branches = compile_to_numba_ir(bfunc, {})
        for lbl, blk in branches.blocks.items():
            for stmt in blk.body:
                if isinstance(stmt, ir.Assign):
                    if isinstance(stmt.value, ir.Global):
                        if stmt.value.name == "PLACEHOLDER_INDEX":
                            stmt.value = index
        return branches

    def apply_transform(self, state):
        # compute new CFG
        func_ir = state.func_ir
        cfg = compute_cfg_from_blocks(func_ir.blocks)
        # find loops
        loops = cfg.loops()

        # 0. Find the loops containing literal_unroll and store this
        #    information
        unroll_info = namedtuple(
            "unroll_info", [
                "loop", "call", "arg", "getitem"])

        def get_call_args(init_arg, want):
            # Chases the assignment of a called value back through a specific
            # call to a global function "want" and returns the arguments
            # supplied to that function's call
            some_call = get_definition(func_ir, init_arg)
            if not isinstance(some_call, ir.Expr):
                raise GuardException
            if not some_call.op == "call":
                raise GuardException
            the_global = get_definition(func_ir, some_call.func)
            if not isinstance(the_global, ir.Global):
                raise GuardException
            if the_global.value is not want:
                raise GuardException
            return some_call

        def find_unroll_loops(loops):
            """This finds loops which are compliant with the form:
            for i in range(len(literal_unroll(<something>>)))"""
            unroll_loops = {}
            for header_lbl, loop in loops.items():
                # TODO: check the loop head has literal_unroll, if it does but
                # does not conform to the following then raise

                # scan loop header
                iternexts = [_ for _ in
                             func_ir.blocks[loop.header].find_exprs('iternext')]
                # needs to be an single iternext driven loop
                if len(iternexts) != 1:
                    continue
                for iternext in iternexts:
                    # Walk the canonicalised loop structure and check it
                    # Check loop form range(literal_unroll(container)))
                    phi = guard(get_definition, func_ir,  iternext.value)
                    if phi is None:
                        continue

                    # check call global "range"
                    range_call = guard(get_call_args, phi.value, range)
                    if range_call is None:
                        continue
                    range_arg = range_call.args[0]

                    # check call global "len"
                    len_call = guard(get_call_args, range_arg, len)
                    if len_call is None:
                        continue
                    len_arg = len_call.args[0]

                    # check literal_unroll
                    literal_unroll_call = guard(get_definition, func_ir,
                                                len_arg)
                    if literal_unroll_call is None:
                        continue
                    if not isinstance(literal_unroll_call, ir.Expr):
                        continue
                    if literal_unroll_call.op != "call":
                        continue
                    literal_func = getattr(literal_unroll_call, 'func', None)
                    if not literal_func:
                        continue
                    call_func = guard(get_definition, func_ir,
                                      literal_unroll_call.func)
                    if call_func is None:
                        continue
                    call_func_value = call_func.value

                    if call_func_value is literal_unroll:
                        assert len(literal_unroll_call.args) == 1
                        unroll_loops[loop] = literal_unroll_call
            return unroll_loops

        def ensure_no_nested_unroll(unroll_loops):
            # Validate loop nests, nested literal_unroll loops are unsupported.
            # This doesn't check that there's a getitem or anything else
            # required for the transform to work, simply just that there's no
            # nesting.
            for test_loop in unroll_loops:
                for ref_loop in unroll_loops:
                    if test_loop == ref_loop:  # comparing to self! skip
                        continue
                    if test_loop.header in ref_loop.body:
                        msg = ("Nesting of literal_unroll is unsupported")
                        loc = func_ir.blocks[test_loop.header].loc
                        raise errors.UnsupportedError(msg, loc)

        def collect_literal_unroll_info(literal_unroll_loops):
            """Finds the loops induced by `literal_unroll`, returns a list of
            unroll_info namedtuples for use in the transform pass.
            """

            literal_unroll_info = []
            for loop, literal_unroll_call in literal_unroll_loops.items():
                arg = literal_unroll_call.args[0]
                typemap = state.typemap
                resolved_arg = guard(get_definition, func_ir, arg,
                                     lhs_only=True)
                ty = typemap[resolved_arg.name]
                assert isinstance(ty, self._accepted_types)
                # loop header is spelled ok, now make sure the body
                # actually contains a getitem

                # find a "getitem"... only looks in the blocks that belong
                # _solely_ to this literal_unroll (there should not be nested
                # literal_unroll loops, this is unsupported).
                tuple_getitem = None
                for lbli in loop.body:
                    blk = func_ir.blocks[lbli]
                    for stmt in blk.body:
                        if isinstance(stmt, ir.Assign):
                            if (isinstance(stmt.value, ir.Expr) and
                                    stmt.value.op == "getitem"):
                                # check for something like a[i]
                                if stmt.value.value != arg:
                                    # that failed, so check for the
                                    # definition
                                    dfn = guard(get_definition, func_ir,
                                                stmt.value.value)
                                    if dfn is None:
                                        continue
                                    try:
                                        args = getattr(dfn, 'args', False)
                                    except KeyError:
                                        continue
                                    if not args:
                                        continue
                                    if not args[0] == arg:
                                        continue
                                target_ty = state.typemap[arg.name]
                                if not isinstance(target_ty,
                                                  self._accepted_types):
                                    continue
                                tuple_getitem = stmt
                                break
                    if tuple_getitem:
                        break
                else:
                    continue  # no getitem in this loop

                ui = unroll_info(loop, literal_unroll_call, arg,
                                 tuple_getitem)
                literal_unroll_info.append(ui)
            return literal_unroll_info

        # 1. Collect info about the literal_unroll loops, ensure they are legal
        literal_unroll_loops = find_unroll_loops(loops)
        # validate
        ensure_no_nested_unroll(literal_unroll_loops)
        # assemble info
        literal_unroll_info = collect_literal_unroll_info(literal_unroll_loops)
        if not literal_unroll_info:
            return False

        # 2. Do the unroll, get a loop and process it!
        info = literal_unroll_info[0]
        self.unroll_loop(state, info)

        # 3. Rebuild the state, the IR has taken a hammering
        func_ir.blocks = simplify_CFG(func_ir.blocks)
        post_proc = postproc.PostProcessor(func_ir)
        post_proc.run()
        if self._DEBUG:
            print('-' * 80 + "END OF PASS, SIMPLIFY DONE")
            func_ir.dump()
        func_ir._definitions = build_definitions(func_ir.blocks)
        return True

    def unroll_loop(self, state, loop_info):
        # The general idea here is to:
        # 1. Find *a* getitem that conforms to the literal_unroll semantic,
        #    i.e. one that is targeting a tuple with a loop induced index
        # 2. Compute a structure from the tuple that describes which
        #    iterations of a loop will have which type
        # 3. Generate a switch table in IR form for the structure in 2
        # 4. Switch out getitems for the tuple for a `typed_getitem`
        # 5. Inject switch table as replacement loop body
        # 6. Patch up
        func_ir = state.func_ir
        getitem_target = loop_info.arg
        target_ty = state.typemap[getitem_target.name]
        assert isinstance(target_ty, self._accepted_types)

        # 1. find a "getitem" that conforms
        tuple_getitem = []
        for lbl in loop_info.loop.body:
            blk = func_ir.blocks[lbl]
            for stmt in blk.body:
                if isinstance(stmt, ir.Assign):
                    if isinstance(stmt.value,
                                  ir.Expr) and stmt.value.op == "getitem":
                        # try a couple of spellings... a[i] and ref(a)[i]
                        if stmt.value.value != getitem_target:
                            dfn = func_ir.get_definition(stmt.value.value)
                            try:
                                args = getattr(dfn, 'args', False)
                            except KeyError:
                                continue
                            if not args:
                                continue
                            if not args[0] == getitem_target:
                                continue
                        target_ty = state.typemap[getitem_target.name]
                        if not isinstance(target_ty, self._accepted_types):
                            continue
                        tuple_getitem.append(stmt)

        if not tuple_getitem:
            msg = ("Loop unrolling analysis has failed, there's no getitem "
                   "in loop body that conforms to literal_unroll "
                   "requirements.")
            LOC = func_ir.blocks[loop_info.loop.header].loc
            raise errors.CompilerError(msg, LOC)

        # 2. get switch data
        switch_data = self.analyse_tuple(target_ty)

        # 3. generate switch IR
        index = func_ir._definitions[tuple_getitem[0].value.index.name][0]
        branches = self.gen_switch(switch_data, index)

        # 4. swap getitems for a typed_getitem, these are actually just
        # placeholders at this point. When the loop is duplicated they can
        # be swapped for a typed_getitem of the correct type or if the item
        # is literal it can be shoved straight into the duplicated loop body
        for item in tuple_getitem:
            old = item.value
            new = ir.Expr.typed_getitem(
                old.value, types.void, old.index, old.loc)
            item.value = new

        # 5. Inject switch table

        # Find the actual loop without the header (that won't get replaced)
        # and derive some new IR for this set of blocks
        this_loop = loop_info.loop
        this_loop_body = this_loop.body - \
            set([this_loop.header])
        loop_blocks = {
            x: func_ir.blocks[x] for x in this_loop_body}
        new_ir = func_ir.derive(loop_blocks)

        # Work out what is live on entry and exit so as to prevent
        # replacement (defined vars can escape, used vars live at the header
        # need to remain as-is so their references are correct, they can
        # also escape).

        usedefs = compute_use_defs(func_ir.blocks)
        idx = this_loop.header
        keep = set()
        keep |= usedefs.usemap[idx] | usedefs.defmap[idx]
        keep |= func_ir.variable_lifetime.livemap[idx]
        dont_replace = [x for x in (keep)]

        # compute the unrolled body
        unrolled_body = self.inject_loop_body(
            branches, new_ir, max(func_ir.blocks.keys()) + 1,
            dont_replace, switch_data)

        # 6. Patch in the unrolled body and fix up
        blks = state.func_ir.blocks
        the_scope = next(iter(blks.values())).scope
        orig_lbl = tuple(this_loop_body)

        replace, *delete = orig_lbl
        unroll, header_block = unrolled_body, this_loop.header
        unroll_lbl = [x for x in sorted(unroll.blocks.keys())]
        blks[replace] = transfer_scope(unroll.blocks[unroll_lbl[0]], the_scope)
        [blks.pop(d) for d in delete]
        for k in unroll_lbl[1:]:
            blks[k] = transfer_scope(unroll.blocks[k], the_scope)
        # stitch up the loop predicate true -> new loop body jump
        blks[header_block].body[-1].truebr = replace

    def run_pass(self, state):
        mutated = False
        func_ir = state.func_ir
        # first limit the work by squashing the CFG if possible
        func_ir.blocks = simplify_CFG(func_ir.blocks)

        if self._DEBUG:
            print("-" * 80 + "PASS ENTRY")
            func_ir.dump()
            print("-" * 80)

        # limitations:
        # 1. No nested unrolls
        # 2. Opt in via `numba.literal_unroll`
        # 3. No multiple mix-tuple use

        # keep running the transform loop until it reports no more changes
        while (True):
            stat = self.apply_transform(state)
            mutated |= stat
            if not stat:
                break

        # reset type inference now we are done with the partial results
        state.typemap = {}
        state.calltypes = None

        return mutated


@register_pass(mutates_CFG=True, analysis_only=False)
class IterLoopCanonicalization(FunctionPass):
    """ Transforms loops that are induced by `getiter` into range() driven loops
    If the typemap is available this will only impact Tuple and UniTuple, if it
    is not available it will impact all matching loops.
    """
    _name = "iter_loop_canonicalisation"

    _DEBUG = False

    # if partial typing info is available it will only look at these types
    _accepted_types = (types.BaseTuple, types.LiteralList)
    _accepted_calls = (literal_unroll,)

    def __init__(self):
        FunctionPass.__init__(self)

    def assess_loop(self, loop, func_ir, partial_typemap=None):
        # it's a iter loop if:
        # - loop header is driven by an iternext
        # - the iternext value is a phi derived from getiter()

        # check header
        iternexts = [_ for _ in
                     func_ir.blocks[loop.header].find_exprs('iternext')]
        if len(iternexts) != 1:
            return False
        for iternext in iternexts:
            phi = guard(get_definition, func_ir,  iternext.value)
            if phi is None:
                return False
            if getattr(phi, 'op', False) == 'getiter':
                if partial_typemap:
                    # check that the call site is accepted, until we're
                    # confident that tuple unrolling is behaving require opt-in
                    # guard of `literal_unroll`, remove this later!
                    phi_val_defn = guard(get_definition, func_ir,  phi.value)
                    if not isinstance(phi_val_defn, ir.Expr):
                        return False
                    if not phi_val_defn.op == "call":
                        return False
                    call = guard(get_definition, func_ir,  phi_val_defn)
                    if call is None or len(call.args) != 1:
                        return False
                    func_var = guard(get_definition, func_ir,  call.func)
                    func = guard(get_definition, func_ir,  func_var)
                    if func is None or not isinstance(func,
                                                      (ir.Global, ir.FreeVar)):
                        return False
                    if (func.value is None or
                            func.value not in self._accepted_calls):
                        return False

                    # now check the type is supported
                    ty = partial_typemap.get(call.args[0].name, None)
                    if ty and isinstance(ty, self._accepted_types):
                        return len(loop.entries) == 1
                else:
                    return len(loop.entries) == 1

    def transform(self, loop, func_ir, cfg):
        def get_range(a):
            return range(len(a))

        iternext = [_ for _ in
                    func_ir.blocks[loop.header].find_exprs('iternext')][0]
        LOC = func_ir.blocks[loop.header].loc
        scope = func_ir.blocks[loop.header].scope
        get_range_var = scope.redefine("CANONICALISER_get_range_gbl", LOC)
        get_range_global = ir.Global('get_range', get_range, LOC)
        assgn = ir.Assign(get_range_global, get_range_var, LOC)

        loop_entry = tuple(loop.entries)[0]
        entry_block = func_ir.blocks[loop_entry]
        entry_block.body.insert(0, assgn)

        iterarg = guard(get_definition, func_ir,  iternext.value)
        if iterarg is not None:
            iterarg = iterarg.value

        # look for iternext
        idx = 0
        for stmt in entry_block.body:
            if isinstance(stmt, ir.Assign):
                if isinstance(stmt.value,
                              ir.Expr) and stmt.value.op == 'getiter':
                    break
            idx += 1
        else:
            raise ValueError("problem")

        # create a range(len(tup)) and inject it
        call_get_range_var = scope.redefine('CANONICALISER_call_get_range', LOC)
        make_call = ir.Expr.call(get_range_var, (stmt.value.value,), (), LOC)
        assgn_call = ir.Assign(make_call, call_get_range_var, LOC)
        entry_block.body.insert(idx, assgn_call)
        entry_block.body[idx + 1].value.value = call_get_range_var

        glbls = copy(func_ir.func_id.func.__globals__)
        from numba.core.inline_closurecall import inline_closure_call
        inline_closure_call(func_ir, glbls, entry_block, idx, get_range,)
        kill = entry_block.body.index(assgn)
        entry_block.body.pop(kill)

        # find the induction variable + references in the loop header
        # fixed point iter to do this, it's a bit clunky
        induction_vars = set()
        header_block = func_ir.blocks[loop.header]

        # find induction var
        ind = [x for x in header_block.find_exprs('pair_first')]
        for x in ind:
            induction_vars.add(func_ir.get_assignee(x, loop.header))
        # find aliases of the induction var
        tmp = set()
        for x in induction_vars:
            try:  # there's not always an alias, e.g. loop from inlined closure
                tmp.add(func_ir.get_assignee(x, loop.header))
            except ValueError:
                pass
        induction_vars |= tmp
        induction_var_names = set([x.name for x in induction_vars])

        # Find the downstream blocks that might reference the induction var
        succ = set()
        for lbl in loop.exits:
            succ |= set([x[0] for x in cfg.successors(lbl)])
        check_blocks = (loop.body | loop.exits | succ) ^ {loop.header}

        # replace RHS use of induction var with getitem
        for lbl in check_blocks:
            for stmt in func_ir.blocks[lbl].body:
                if isinstance(stmt, ir.Assign):
                    # check for aliases
                    try:
                        lookup = getattr(stmt.value, 'name', None)
                    except KeyError:
                        continue
                    if lookup and lookup in induction_var_names:
                        stmt.value = ir.Expr.getitem(
                            iterarg, stmt.value, stmt.loc)

        post_proc = postproc.PostProcessor(func_ir)
        post_proc.run()

    def run_pass(self, state):
        func_ir = state.func_ir
        cfg = compute_cfg_from_blocks(func_ir.blocks)
        loops = cfg.loops()

        mutated = False
        for header, loop in loops.items():
            stat = self.assess_loop(loop, func_ir, state.typemap)
            if stat:
                if self._DEBUG:
                    print("Canonicalising loop", loop)
                self.transform(loop, func_ir, cfg)
                mutated = True
            else:
                if self._DEBUG:
                    print("NOT Canonicalising loop", loop)

        func_ir.blocks = simplify_CFG(func_ir.blocks)
        return mutated


@register_pass(mutates_CFG=False, analysis_only=False)
class PropagateLiterals(FunctionPass):
    """Implement literal propagation based on partial type inference"""
    _name = "PropagateLiterals"

    def __init__(self):
        FunctionPass.__init__(self)

    def get_analysis_usage(self, AU):
        AU.add_required(ReconstructSSA)

    def run_pass(self, state):
        func_ir = state.func_ir
        typemap = state.typemap
        flags = state.flags

        accepted_functions = ('isinstance', 'hasattr')

        if not hasattr(func_ir, '_definitions') \
                and not flags.enable_ssa:
            func_ir._definitions = build_definitions(func_ir.blocks)

        changed = False

        for block in func_ir.blocks.values():
            for assign in block.find_insts(ir.Assign):
                value = assign.value
                if isinstance(value, (ir.Arg, ir.Const, ir.FreeVar, ir.Global)):
                    continue

                # 1) Don't change return stmt in the form
                # $return_xyz = cast(value=ABC)
                # 2) Don't propagate literal values that are not primitives
                if isinstance(value, ir.Expr) and \
                        value.op in ('cast', 'build_map', 'build_list',
                                     'build_tuple', 'build_set'):
                    continue

                target = assign.target
                if not flags.enable_ssa:
                    # SSA is disabled when doing inlining
                    if guard(get_definition, func_ir, target.name) is None:  # noqa: E501
                        continue

                # Numba cannot safely determine if an isinstance call
                # with a PHI node is True/False. For instance, in
                # the case below, the partial type inference step can coerce
                # '$z' to float, so any call to 'isinstance(z, int)' would fail.
                #
                #   def fn(x):
                #       if x > 4:
                #           z = 1
                #       else:
                #           z = 3.14
                #       if isinstance(z, int):
                #           print('int')
                #       else:
                #           print('float')
                #
                # At the moment, one avoid propagating the literal
                # value if the argument is a PHI node

                if isinstance(value, ir.Expr) and value.op == 'call':

                    fn = guard(get_definition, func_ir, value.func.name)
                    if fn is None:
                        continue

                    if not (isinstance(fn, ir.Global) and fn.name in
                            accepted_functions):
                        continue

                    for arg in value.args:
                        # check if any of the args to isinstance is a PHI node
                        iv = func_ir._definitions[arg.name]
                        assert len(iv) == 1  # SSA!
                        if isinstance(iv[0], ir.Expr) and iv[0].op == 'phi':
                            msg = (f'{fn.name}() cannot determine the '
                                   f'type of variable "{arg.unversioned_name}" '
                                   'due to a branch.')
                            raise errors.NumbaTypeError(msg, loc=assign.loc)

                # Only propagate a PHI node if all arguments are the same
                # constant
                if isinstance(value, ir.Expr) and value.op == 'phi':
                    # typemap will return None in case `inc.name` not in typemap
                    v = [typemap.get(inc.name) for inc in value.incoming_values]
                    # stop if the elements in `v` do not hold the same value
                    if v[0] is not None and any([v[0] != vi for vi in v]):
                        continue

                lit = typemap.get(target.name, None)
                if lit and isinstance(lit, types.Literal):
                    # replace assign instruction by ir.Const(lit) iff
                    # lit is a literal value
                    rhs = ir.Const(lit.literal_value, assign.loc)
                    new_assign = ir.Assign(rhs, target, assign.loc)

                    # replace instruction
                    block.insert_after(new_assign, assign)
                    block.remove(assign)

                    changed = True

        # reset type inference now we are done with the partial results
        state.typemap = None
        state.calltypes = None

        if changed:
            # Rebuild definitions
            func_ir._definitions = build_definitions(func_ir.blocks)

        return changed


@register_pass(mutates_CFG=True, analysis_only=False)
class LiteralPropagationSubPipelinePass(FunctionPass):
    """Implement literal propagation based on partial type inference"""
    _name = "LiteralPropagation"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        # Determine whether to even attempt this pass... if there's no
        # `isinstance` as a global or as a freevar then just skip.

        found = False
        func_ir = state.func_ir
        for blk in func_ir.blocks.values():
            for asgn in blk.find_insts(ir.Assign):
                if isinstance(asgn.value, (ir.Global, ir.FreeVar)):
                    value = asgn.value.value
                    if value is isinstance or value is hasattr:
                        found = True
                        break
            if found:
                break
        if not found:
            return False

        # run as subpipeline
        from numba.core.compiler_machinery import PassManager
        from numba.core.typed_passes import PartialTypeInference
        pm = PassManager("literal_propagation_subpipeline")

        pm.add_pass(PartialTypeInference, "performs partial type inference")
        pm.add_pass(PropagateLiterals, "performs propagation of literal values")

        # rewrite consts / dead branch pruning
        pm.add_pass(RewriteSemanticConstants, "rewrite semantic constants")
        pm.add_pass(DeadBranchPrune, "dead branch pruning")

        pm.finalize()
        pm.run(state)
        return True

    def get_analysis_usage(self, AU):
        AU.add_required(ReconstructSSA)


@register_pass(mutates_CFG=True, analysis_only=False)
class LiteralUnroll(FunctionPass):
    """Implement the literal_unroll semantics"""
    _name = "literal_unroll"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        # Determine whether to even attempt this pass... if there's no
        # `literal_unroll` as a global or as a freevar then just skip.
        found = False
        func_ir = state.func_ir
        for blk in func_ir.blocks.values():
            for asgn in blk.find_insts(ir.Assign):
                if isinstance(asgn.value, (ir.Global, ir.FreeVar)):
                    if asgn.value.value is literal_unroll:
                        found = True
                        break
            if found:
                break
        if not found:
            return False

        # run as subpipeline
        from numba.core.compiler_machinery import PassManager
        from numba.core.typed_passes import PartialTypeInference
        pm = PassManager("literal_unroll_subpipeline")
        # get types where possible to help with list->tuple change
        pm.add_pass(PartialTypeInference, "performs partial type inference")
        # make const lists tuples
        pm.add_pass(TransformLiteralUnrollConstListToTuple,
                    "switch const list for tuples")
        # recompute partial typemap following IR change
        pm.add_pass(PartialTypeInference, "performs partial type inference")
        # canonicalise loops
        pm.add_pass(IterLoopCanonicalization,
                    "switch iter loops for range driven loops")
        # rewrite consts
        pm.add_pass(RewriteSemanticConstants, "rewrite semantic constants")
        # do the unroll
        pm.add_pass(MixedContainerUnroller, "performs mixed container unroll")
        # rewrite dynamic getitem to static getitem as it's possible some more
        # getitems will now be statically resolvable
        pm.add_pass(GenericRewrites, "Generic Rewrites")
        pm.add_pass(RewriteSemanticConstants, "rewrite semantic constants")
        pm.finalize()
        pm.run(state)
        return True


@register_pass(mutates_CFG=True, analysis_only=False)
class SimplifyCFG(FunctionPass):
    """Perform CFG simplification"""
    _name = "simplify_cfg"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        blks = state.func_ir.blocks
        new_blks = simplify_CFG(blks)
        state.func_ir.blocks = new_blks
        mutated = blks != new_blks
        return mutated


@register_pass(mutates_CFG=False, analysis_only=False)
class ReconstructSSA(FunctionPass):
    """Perform SSA-reconstruction

    Produces minimal SSA.
    """
    _name = "reconstruct_ssa"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        state.func_ir = reconstruct_ssa(state.func_ir)
        self._patch_locals(state)

        # Rebuild definitions
        state.func_ir._definitions = build_definitions(state.func_ir.blocks)

        # Rerun postprocessor to update metadata
        # example generator_info
        post_proc = postproc.PostProcessor(state.func_ir)
        post_proc.run(emit_dels=False)

        if config.DEBUG or config.DUMP_SSA:
            name = state.func_ir.func_id.func_qualname
            print(f"SSA IR DUMP: {name}".center(80, "-"))
            state.func_ir.dump()

        return True      # XXX detect if it actually got changed

    def _patch_locals(self, state):
        # Fix dispatcher locals dictionary type annotation
        locals_dict = state.get('locals')
        if locals_dict is None:
            return

        first_blk, *_ = state.func_ir.blocks.values()
        scope = first_blk.scope
        for parent, redefs in scope.var_redefinitions.items():
            if parent in locals_dict:
                typ = locals_dict[parent]
                for derived in redefs:
                    locals_dict[derived] = typ


@register_pass(mutates_CFG=False, analysis_only=False)
class RewriteDynamicRaises(FunctionPass):
    """Replace existing raise statements by dynamic raises in Numba IR.
    """
    _name = "Rewrite dynamic raises"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        func_ir = state.func_ir
        changed = False

        for block in func_ir.blocks.values():
            for raise_ in block.find_insts((ir.Raise, ir.TryRaise)):
                call_inst = guard(get_definition, func_ir, raise_.exception)
                if call_inst is None:
                    continue
                exc_type = func_ir.infer_constant(call_inst.func.name)
                exc_args = []
                for exc_arg in call_inst.args:
                    try:
                        const = func_ir.infer_constant(exc_arg)
                        exc_args.append(const)
                    except consts.ConstantInferenceError:
                        exc_args.append(exc_arg)
                loc = raise_.loc

                cls = {
                    ir.TryRaise: ir.DynamicTryRaise,
                    ir.Raise: ir.DynamicRaise,
                }[type(raise_)]

                dyn_raise = cls(exc_type, tuple(exc_args), loc)
                block.insert_after(dyn_raise, raise_)
                block.remove(raise_)
                changed = True
        return changed
