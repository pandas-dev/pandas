from collections import namedtuple, defaultdict
import operator
import warnings
from functools import partial

import llvmlite.ir
from llvmlite.ir import Constant, IRBuilder

from numba.core import (typing, utils, types, ir, debuginfo, funcdesc,
                        generators, config, ir_utils, cgutils, removerefctpass,
                        targetconfig)
from numba.core.errors import (LoweringError, new_error_context, TypingError,
                               LiteralTypingError, UnsupportedError,
                               NumbaDebugInfoWarning)
from numba.core.funcdesc import default_mangler
from numba.core.environment import Environment
from numba.core.analysis import compute_use_defs, must_use_alloca
from numba.misc.firstlinefinder import get_func_body_first_lineno


_VarArgItem = namedtuple("_VarArgItem", ("vararg", "index"))


class BaseLower(object):
    """
    Lower IR to LLVM
    """

    def __init__(self, context, library, fndesc, func_ir, metadata=None):
        self.library = library
        self.fndesc = fndesc
        self.blocks = utils.SortedMap(func_ir.blocks.items())
        self.func_ir = func_ir
        self.generator_info = func_ir.generator_info
        self.metadata = metadata
        self.flags = targetconfig.ConfigStack.top_or_none()

        # Initialize LLVM
        self.module = self.library.create_ir_module(self.fndesc.unique_name)

        # Python execution environment (will be available to the compiled
        # function).
        self.env = Environment.from_fndesc(self.fndesc)

        # Internal states
        self.blkmap = {}
        self.pending_phis = {}
        self.varmap = {}
        self.firstblk = min(self.blocks.keys())
        self.loc = -1

        # Specializes the target context as seen inside the Lowerer
        # This adds:
        #  - environment: the python execution environment
        self.context = context.subtarget(environment=self.env,
                                         fndesc=self.fndesc)

        # Debuginfo
        dibuildercls = (self.context.DIBuilder
                        if self.context.enable_debuginfo
                        else debuginfo.DummyDIBuilder)

        # debuginfo def location
        self.defn_loc = self._compute_def_location()

        directives_only = self.flags.dbg_directives_only
        self.debuginfo = dibuildercls(module=self.module,
                                      filepath=func_ir.loc.filename,
                                      cgctx=context,
                                      directives_only=directives_only)

        # Subclass initialization
        self.init()

    @property
    def call_conv(self):
        return self.context.call_conv

    def init(self):
        pass

    def init_pyapi(self):
        """
        Init the Python API and Environment Manager for the function being
        lowered.
        """
        if self.pyapi is not None:
            return
        self.pyapi = self.context.get_python_api(self.builder)

        # Store environment argument for later use
        self.env_manager = self.context.get_env_manager(self.builder)
        self.env_body = self.env_manager.env_body
        self.envarg = self.env_manager.env_ptr

    def _compute_def_location(self):
        # Debuginfo requires source to be accurate. Find it and warn if not
        # found. If it's not found, use the func_ir line + 1, this assumes that
        # the function definition is decorated with a 1 line jit decorator.
        defn_loc = self.func_ir.loc.with_lineno(self.func_ir.loc.line + 1)
        if self.context.enable_debuginfo:
            fn = self.func_ir.func_id.func
            optional_lno = get_func_body_first_lineno(fn)
            if optional_lno is not None:
                # -1 as lines start at 1 and this is an offset.
                offset = optional_lno - 1
                defn_loc = self.func_ir.loc.with_lineno(offset)
            else:
                msg = ("Could not find source for function: "
                       f"{self.func_ir.func_id.func}. Debug line information "
                       "may be inaccurate.")
                warnings.warn(NumbaDebugInfoWarning(msg))
        return defn_loc

    def pre_lower(self):
        """
        Called before lowering all blocks.
        """
        # A given Lower object can be used for several LL functions
        # (for generators) and it's important to use a new API and
        # EnvironmentManager.
        self.pyapi = None
        self.debuginfo.mark_subprogram(function=self.builder.function,
                                       qualname=self.fndesc.qualname,
                                       argnames=self.fndesc.args,
                                       argtypes=self.fndesc.argtypes,
                                       line=self.defn_loc.line)

        # When full debug info is enabled, disable inlining where possible, to
        # improve the quality of the debug experience. 'alwaysinline' functions
        # cannot have inlining disabled.
        attributes = self.builder.function.attributes
        full_debug = self.flags.debuginfo and not self.flags.dbg_directives_only
        if full_debug and 'alwaysinline' not in attributes:
            attributes.add('noinline')

    def post_lower(self):
        """
        Called after all blocks are lowered
        """
        self.debuginfo.finalize()

    def pre_block(self, block):
        """
        Called before lowering a block.
        """

    def post_block(self, block):
        """
        Called after lowering a block.
        """

    def return_dynamic_exception(self, exc_class, exc_args, nb_types, loc=None):
        self.call_conv.return_dynamic_user_exc(
            self.builder, exc_class, exc_args, nb_types,
            loc=loc, func_name=self.func_ir.func_id.func_name,
        )

    def return_exception(self, exc_class, exc_args=None, loc=None):
        """Propagate exception to the caller.
        """
        self.call_conv.return_user_exc(
            self.builder, exc_class, exc_args,
            loc=loc, func_name=self.func_ir.func_id.func_name,
        )

    def set_exception(self, exc_class, exc_args=None, loc=None):
        """Set exception state in the current function.
        """
        self.call_conv.set_static_user_exc(
            self.builder, exc_class, exc_args,
            loc=loc, func_name=self.func_ir.func_id.func_name,
        )

    def emit_environment_object(self):
        """Emit a pointer to hold the Environment object.
        """
        # Define global for the environment and initialize it to NULL
        envname = self.context.get_env_name(self.fndesc)
        self.context.declare_env_global(self.module, envname)

    def lower(self):
        # Emit the Env into the module
        self.emit_environment_object()
        if self.generator_info is None:
            self.genlower = None
            self.lower_normal_function(self.fndesc)
        else:
            self.genlower = self.GeneratorLower(self)
            self.gentype = self.genlower.gentype

            self.genlower.lower_init_func(self)
            self.genlower.lower_next_func(self)
            if self.gentype.has_finalizer:
                self.genlower.lower_finalize_func(self)

        if config.DUMP_LLVM:
            utils.dump_llvm(self.fndesc, self.module)

        # Special optimization to remove NRT on functions that do not need it.
        if self.context.enable_nrt and self.generator_info is None:
            removerefctpass.remove_unnecessary_nrt_usage(self.function,
                                                         context=self.context,
                                                         fndesc=self.fndesc)

        # Run target specific post lowering transformation
        self.context.post_lowering(self.module, self.library)

        # Materialize LLVM Module
        self.library.add_ir_module(self.module)

    def extract_function_arguments(self):
        self.fnargs = self.call_conv.decode_arguments(self.builder,
                                                      self.fndesc.argtypes,
                                                      self.function)
        return self.fnargs

    def lower_normal_function(self, fndesc):
        """
        Lower non-generator *fndesc*.
        """
        self.setup_function(fndesc)

        # Init argument values
        self.extract_function_arguments()
        entry_block_tail = self.lower_function_body()

        # Close tail of entry block, do not emit debug metadata else the
        # unconditional jump gets associated with the metadata from the function
        # body end.
        with debuginfo.suspend_emission(self.builder):
            self.builder.position_at_end(entry_block_tail)
            self.builder.branch(self.blkmap[self.firstblk])

    def lower_function_body(self):
        """
        Lower the current function's body, and return the entry block.
        """
        # Init Python blocks
        for offset in self.blocks:
            bname = "B%s" % offset
            self.blkmap[offset] = self.function.append_basic_block(bname)

        self.pre_lower()
        # pre_lower() may have changed the current basic block
        entry_block_tail = self.builder.basic_block

        self.debug_print("# function begin: {0}".format(
            self.fndesc.unique_name))

        # Lower all blocks
        for offset, block in sorted(self.blocks.items()):
            bb = self.blkmap[offset]
            self.builder.position_at_end(bb)
            self.debug_print(f"# lower block: {offset}")
            self.lower_block(block)
        self.post_lower()
        return entry_block_tail

    def lower_block(self, block):
        """
        Lower the given block.
        """
        self.pre_block(block)
        for inst in block.body:
            self.loc = inst.loc
            defaulterrcls = partial(LoweringError, loc=self.loc)
            with new_error_context('lowering "{inst}" at {loc}', inst=inst,
                                   loc=self.loc, errcls_=defaulterrcls):
                self.lower_inst(inst)
        self.post_block(block)

    def create_cpython_wrapper(self, release_gil=False):
        """
        Create CPython wrapper(s) around this function (or generator).
        """
        if self.genlower:
            self.context.create_cpython_wrapper(self.library,
                                                self.genlower.gendesc,
                                                self.env, self.call_helper,
                                                release_gil=release_gil)
        self.context.create_cpython_wrapper(self.library, self.fndesc,
                                            self.env, self.call_helper,
                                            release_gil=release_gil)

    def create_cfunc_wrapper(self):
        """
        Create C wrapper around this function.
        """
        if self.genlower:
            raise UnsupportedError('generator as a first-class function type')
        self.context.create_cfunc_wrapper(self.library, self.fndesc,
                                          self.env, self.call_helper)

    def setup_function(self, fndesc):
        # Setup function
        self.function = self.context.declare_function(self.module, fndesc)
        if self.flags.dbg_optnone:
            attrset = self.function.attributes
            if "alwaysinline" not in attrset:
                attrset.add("optnone")
                attrset.add("noinline")
        self.entry_block = self.function.append_basic_block('entry')
        self.builder = IRBuilder(self.entry_block)
        self.call_helper = self.call_conv.init_call_helper(self.builder)

    def typeof(self, varname):
        return self.fndesc.typemap[varname]

    def debug_print(self, msg):
        if config.DEBUG_JIT:
            self.context.debug_print(
                self.builder, f"DEBUGJIT [{self.fndesc.qualname}]: {msg}")

    def print_variable(self, msg, varname):
        """Helper to emit ``print(msg, varname)`` for debugging.

        Parameters
        ----------
        msg : str
            Literal string to be printed.
        varname : str
            A variable name whose value will be printed.
        """
        argtys = (
            types.literal(msg),
            self.fndesc.typemap[varname]
        )
        args = (
            self.context.get_dummy_value(),
            self.loadvar(varname),
        )
        sig = typing.signature(types.none, *argtys)

        impl = self.context.get_function(print, sig)
        impl(self.builder, args)


class Lower(BaseLower):
    GeneratorLower = generators.GeneratorLower

    def init(self):
        super().init()
        # find all singly assigned variables
        self._find_singly_assigned_variable()

    @property
    def _disable_sroa_like_opt(self):
        """Flags that the SROA like optimisation that Numba performs (which
        prevent alloca and subsequent load/store for locals) should be disabled.
        Currently, this is conditional solely on the presence of a request for
        the emission of debug information."""
        if self.flags is None:
            return False

        return self.flags.debuginfo and not self.flags.dbg_directives_only

    def _find_singly_assigned_variable(self):
        func_ir = self.func_ir
        blocks = func_ir.blocks

        sav = set()

        if not self.func_ir.func_id.is_generator:
            use_defs = compute_use_defs(blocks)
            alloca_vars = must_use_alloca(blocks)

            # Compute where variables are defined
            var_assign_map = defaultdict(set)
            for blk, vl in use_defs.defmap.items():
                for var in vl:
                    var_assign_map[var].add(blk)

            # Compute where variables are used
            var_use_map = defaultdict(set)
            for blk, vl in use_defs.usemap.items():
                for var in vl:
                    var_use_map[var].add(blk)

            # Keep only variables that are defined locally and used locally
            for var in var_assign_map:
                if var not in alloca_vars and len(var_assign_map[var]) == 1:
                    # Usemap does not keep locally defined variables.
                    if len(var_use_map[var]) == 0:
                        # Ensure that the variable is not defined multiple times
                        # in the block
                        [defblk] = var_assign_map[var]
                        assign_stmts = self.blocks[defblk].find_insts(ir.Assign)
                        assigns = [stmt for stmt in assign_stmts
                                   if stmt.target.name == var]
                        if len(assigns) == 1:
                            sav.add(var)

        self._singly_assigned_vars = sav
        self._blk_local_varmap = {}

    def pre_block(self, block):
        from numba.core.unsafe import eh

        super(Lower, self).pre_block(block)
        self._cur_ir_block = block

        if block == self.firstblk:
            # create slots for all the vars, irrespective of whether they are
            # initialized, SSA will pick this up and warn users about using
            # uninitialized variables. Slots are added as alloca in the first
            # block
            bb = self.blkmap[self.firstblk]
            self.builder.position_at_end(bb)
            all_names = set()
            for block in self.blocks.values():
                for x in block.find_insts(ir.Del):
                    if x.value not in all_names:
                        all_names.add(x.value)
            for name in all_names:
                fetype = self.typeof(name)
                self._alloca_var(name, fetype)

        # Detect if we are in a TRY block by looking for a call to
        # `eh.exception_check`.
        for call in block.find_exprs(op='call'):
            defn = ir_utils.guard(
                ir_utils.get_definition, self.func_ir, call.func,
            )
            if defn is not None and isinstance(defn, ir.Global):
                if defn.value is eh.exception_check:
                    if isinstance(block.terminator, ir.Branch):
                        targetblk = self.blkmap[block.terminator.truebr]
                        # NOTE: This hacks in an attribute for call_conv to
                        #       pick up. This hack is no longer needed when
                        #       all old-style implementations are gone.
                        self.builder._in_try_block = {'target': targetblk}
                        break

    def post_block(self, block):
        # Clean-up
        try:
            del self.builder._in_try_block
        except AttributeError:
            pass

    def lower_inst(self, inst):
        # Set debug location for all subsequent LL instructions
        self.debuginfo.mark_location(self.builder, self.loc.line)
        self.debug_print(str(inst))
        if isinstance(inst, ir.Assign):
            ty = self.typeof(inst.target.name)
            val = self.lower_assign(ty, inst)
            argidx = None
            # If this is a store from an arg, like x = arg.x then tell debuginfo
            # that this is the arg
            if isinstance(inst.value, ir.Arg):
                # NOTE: debug location is the `def <func>` line
                self.debuginfo.mark_location(self.builder, self.defn_loc.line)
                argidx = inst.value.index + 1 # args start at 1
            self.storevar(val, inst.target.name, argidx=argidx)

        elif isinstance(inst, ir.Branch):
            cond = self.loadvar(inst.cond.name)
            tr = self.blkmap[inst.truebr]
            fl = self.blkmap[inst.falsebr]

            condty = self.typeof(inst.cond.name)
            pred = self.context.cast(self.builder, cond, condty, types.boolean)
            assert pred.type == llvmlite.ir.IntType(1),\
                ("cond is not i1: %s" % pred.type)
            self.builder.cbranch(pred, tr, fl)

        elif isinstance(inst, ir.Jump):
            target = self.blkmap[inst.target]
            self.builder.branch(target)

        elif isinstance(inst, ir.Return):
            if self.generator_info:
                # StopIteration
                self.genlower.return_from_generator(self)
                return
            val = self.loadvar(inst.value.name)
            oty = self.typeof(inst.value.name)
            ty = self.fndesc.restype
            if isinstance(ty, types.Optional):
                # If returning an optional type
                self.call_conv.return_optional_value(self.builder, ty, oty, val)
                return
            assert ty == oty, (
                "type '{}' does not match return type '{}'".format(oty, ty))
            retval = self.context.get_return_value(self.builder, ty, val)
            self.call_conv.return_value(self.builder, retval)

        elif isinstance(inst, ir.PopBlock):
            pass # this is just a marker

        elif isinstance(inst, ir.StaticSetItem):
            signature = self.fndesc.calltypes[inst]
            assert signature is not None
            try:
                impl = self.context.get_function('static_setitem', signature)
            except NotImplementedError:
                return self.lower_setitem(inst.target, inst.index_var,
                                          inst.value, signature)
            else:
                target = self.loadvar(inst.target.name)
                value = self.loadvar(inst.value.name)
                valuety = self.typeof(inst.value.name)
                value = self.context.cast(self.builder, value, valuety,
                                          signature.args[2])
                return impl(self.builder, (target, inst.index, value))

        elif isinstance(inst, ir.Print):
            self.lower_print(inst)

        elif isinstance(inst, ir.SetItem):
            signature = self.fndesc.calltypes[inst]
            assert signature is not None
            return self.lower_setitem(inst.target, inst.index, inst.value,
                                      signature)

        elif isinstance(inst, ir.StoreMap):
            signature = self.fndesc.calltypes[inst]
            assert signature is not None
            return self.lower_setitem(inst.dct, inst.key, inst.value, signature)

        elif isinstance(inst, ir.DelItem):
            target = self.loadvar(inst.target.name)
            index = self.loadvar(inst.index.name)

            targetty = self.typeof(inst.target.name)
            indexty = self.typeof(inst.index.name)

            signature = self.fndesc.calltypes[inst]
            assert signature is not None

            op = operator.delitem
            fnop = self.context.typing_context.resolve_value_type(op)
            callsig = fnop.get_call_type(
                self.context.typing_context, signature.args, {},
            )
            impl = self.context.get_function(fnop, callsig)

            assert targetty == signature.args[0]
            index = self.context.cast(self.builder, index, indexty,
                                      signature.args[1])

            return impl(self.builder, (target, index))

        elif isinstance(inst, ir.Del):
            self.delvar(inst.value)

        elif isinstance(inst, ir.SetAttr):
            target = self.loadvar(inst.target.name)
            value = self.loadvar(inst.value.name)
            signature = self.fndesc.calltypes[inst]

            targetty = self.typeof(inst.target.name)
            valuety = self.typeof(inst.value.name)
            assert signature is not None
            assert signature.args[0] == targetty
            impl = self.context.get_setattr(inst.attr, signature)

            # Convert argument to match
            value = self.context.cast(self.builder, value, valuety,
                                      signature.args[1])

            return impl(self.builder, (target, value))

        elif isinstance(inst, ir.DynamicRaise):
            self.lower_dynamic_raise(inst)

        elif isinstance(inst, ir.DynamicTryRaise):
            self.lower_try_dynamic_raise(inst)

        elif isinstance(inst, ir.StaticRaise):
            self.lower_static_raise(inst)

        elif isinstance(inst, ir.StaticTryRaise):
            self.lower_static_try_raise(inst)

        else:
            raise NotImplementedError(type(inst))

    def lower_setitem(self, target_var, index_var, value_var, signature):
        target = self.loadvar(target_var.name)
        value = self.loadvar(value_var.name)
        index = self.loadvar(index_var.name)

        targetty = self.typeof(target_var.name)
        valuety = self.typeof(value_var.name)
        indexty = self.typeof(index_var.name)

        op = operator.setitem
        fnop = self.context.typing_context.resolve_value_type(op)
        callsig = fnop.get_call_type(
            self.context.typing_context, signature.args, {},
        )
        impl = self.context.get_function(fnop, callsig)

        # Convert argument to match
        if isinstance(targetty, types.Optional):
            target = self.context.cast(self.builder, target, targetty,
                                       targetty.type)
        else:
            ul = types.unliteral
            assert ul(targetty) == ul(signature.args[0])

        index = self.context.cast(self.builder, index, indexty,
                                  signature.args[1])
        value = self.context.cast(self.builder, value, valuety,
                                  signature.args[2])

        return impl(self.builder, (target, index, value))

    def lower_try_dynamic_raise(self, inst):
        # Numba is a bit limited in what it can do with exceptions in a try
        # block. Thus, it is safe to use the same code as the static try raise.
        self.lower_static_try_raise(inst)

    def lower_dynamic_raise(self, inst):
        exc_args = inst.exc_args
        args = []
        nb_types = []
        for exc_arg in exc_args:
            if isinstance(exc_arg, ir.Var):
                # dynamic values
                typ = self.typeof(exc_arg.name)
                val = self.loadvar(exc_arg.name)
                self.incref(typ, val)
            else:
                typ = None
                val = exc_arg
            nb_types.append(typ)
            args.append(val)

        self.return_dynamic_exception(inst.exc_class, tuple(args),
                                      tuple(nb_types), loc=self.loc)

    def lower_static_raise(self, inst):
        if inst.exc_class is None:
            # Reraise
            self.return_exception(None, loc=self.loc)
        else:
            self.return_exception(inst.exc_class, inst.exc_args, loc=self.loc)

    def lower_static_try_raise(self, inst):
        if inst.exc_class is None:
            # Reraise
            self.set_exception(None, loc=self.loc)
        else:
            self.set_exception(inst.exc_class, inst.exc_args, loc=self.loc)

    def lower_assign(self, ty, inst):
        value = inst.value
        # In nopython mode, closure vars are frozen like globals
        if isinstance(value, (ir.Const, ir.Global, ir.FreeVar)):
            res = self.context.get_constant_generic(self.builder, ty,
                                                    value.value)
            self.incref(ty, res)
            return res

        elif isinstance(value, ir.Expr):
            return self.lower_expr(ty, value)

        elif isinstance(value, ir.Var):
            val = self.loadvar(value.name)
            oty = self.typeof(value.name)
            res = self.context.cast(self.builder, val, oty, ty)
            self.incref(ty, res)
            return res

        elif isinstance(value, ir.Arg):
            # Suspend debug info else all the arg repacking ends up being
            # associated with some line or other and it's actually just a detail
            # of Numba's CC.
            with debuginfo.suspend_emission(self.builder):
                # Cast from the argument type to the local variable type
                # (note the "arg.FOO" convention as used in typeinfer)
                argty = self.typeof("arg." + value.name)
                if isinstance(argty, types.Omitted):
                    pyval = argty.value
                    tyctx = self.context.typing_context
                    valty = tyctx.resolve_value_type_prefer_literal(pyval)
                    # use the type of the constant value
                    const = self.context.get_constant_generic(
                        self.builder, valty, pyval,
                    )
                    # cast it to the variable type
                    res = self.context.cast(self.builder, const, valty, ty)
                else:
                    val = self.fnargs[value.index]
                    res = self.context.cast(self.builder, val, argty, ty)
                self.incref(ty, res)
                return res

        elif isinstance(value, ir.Yield):
            res = self.lower_yield(ty, value)
            self.incref(ty, res)
            return res

        raise NotImplementedError(type(value), value)

    def lower_yield(self, retty, inst):
        yp = self.generator_info.yield_points[inst.index]
        assert yp.inst is inst
        y = generators.LowerYield(self, yp, yp.live_vars)
        y.lower_yield_suspend()
        # Yield to caller
        val = self.loadvar(inst.value.name)
        typ = self.typeof(inst.value.name)
        actual_rettyp = self.gentype.yield_type

        # cast the local val to the type yielded
        yret = self.context.cast(self.builder, val, typ, actual_rettyp)

        # get the return repr of yielded value
        retval = self.context.get_return_value(
            self.builder, actual_rettyp, yret,
        )

        # return
        self.call_conv.return_value(self.builder, retval)

        # Resumption point
        y.lower_yield_resume()
        # None is returned by the yield expression
        return self.context.get_constant_generic(self.builder, retty, None)

    def lower_binop(self, resty, expr, op):
        # if op in utils.OPERATORS_TO_BUILTINS:
        # map operator.the_op => the corresponding types.Function()
        # TODO: is this looks dodgy ...
        op = self.context.typing_context.resolve_value_type(op)

        lhs = expr.lhs
        rhs = expr.rhs
        static_lhs = expr.static_lhs
        static_rhs = expr.static_rhs
        lty = self.typeof(lhs.name)
        rty = self.typeof(rhs.name)
        lhs = self.loadvar(lhs.name)
        rhs = self.loadvar(rhs.name)

        # Convert argument to match
        signature = self.fndesc.calltypes[expr]
        lhs = self.context.cast(self.builder, lhs, lty, signature.args[0])
        rhs = self.context.cast(self.builder, rhs, rty, signature.args[1])

        def cast_result(res):
            return self.context.cast(self.builder, res,
                                     signature.return_type, resty)

        # First try with static operands, if known
        def try_static_impl(tys, args):
            if any(a is ir.UNDEFINED for a in args):
                return None
            try:
                if isinstance(op, types.Function):
                    static_sig = op.get_call_type(self.context.typing_context,
                                                  tys, {})
                else:
                    static_sig = typing.signature(signature.return_type, *tys)
            except TypingError:
                return None
            try:
                static_impl = self.context.get_function(op, static_sig)
                return static_impl(self.builder, args)
            except NotImplementedError:
                return None

        res = try_static_impl(
            (_lit_or_omitted(static_lhs), _lit_or_omitted(static_rhs)),
            (static_lhs, static_rhs),
        )
        if res is not None:
            return cast_result(res)

        res = try_static_impl(
            (_lit_or_omitted(static_lhs), rty),
            (static_lhs, rhs),
        )
        if res is not None:
            return cast_result(res)

        res = try_static_impl(
            (lty, _lit_or_omitted(static_rhs)),
            (lhs, static_rhs),
        )
        if res is not None:
            return cast_result(res)

        # Normal implementation for generic arguments

        sig = op.get_call_type(self.context.typing_context, signature.args, {})
        impl = self.context.get_function(op, sig)
        res = impl(self.builder, (lhs, rhs))
        return cast_result(res)

    def lower_getitem(self, resty, expr, value, index, signature):
        baseval = self.loadvar(value.name)
        indexval = self.loadvar(index.name)
        # Get implementation of getitem
        op = operator.getitem
        fnop = self.context.typing_context.resolve_value_type(op)
        callsig = fnop.get_call_type(
            self.context.typing_context, signature.args, {},
        )
        impl = self.context.get_function(fnop, callsig)

        argvals = (baseval, indexval)
        argtyps = (self.typeof(value.name),
                   self.typeof(index.name))
        castvals = [self.context.cast(self.builder, av, at, ft)
                    for av, at, ft in zip(argvals, argtyps,
                                          signature.args)]
        res = impl(self.builder, castvals)
        return self.context.cast(self.builder, res,
                                 signature.return_type,
                                 resty)

    def _cast_var(self, var, ty):
        """
        Cast a Numba IR variable to the given Numba type, returning a
        low-level value.
        """
        if isinstance(var, _VarArgItem):
            varty = self.typeof(var.vararg.name)[var.index]
            val = self.builder.extract_value(self.loadvar(var.vararg.name),
                                             var.index)
        else:
            varty = self.typeof(var.name)
            val = self.loadvar(var.name)
        return self.context.cast(self.builder, val, varty, ty)

    def fold_call_args(self, fnty, signature, pos_args, vararg, kw_args):
        if vararg:
            # Inject *args from function call
            # The lowering will be done in _cast_var() above.
            tp_vararg = self.typeof(vararg.name)
            assert isinstance(tp_vararg, types.BaseTuple)
            pos_args = pos_args + [_VarArgItem(vararg, i)
                                   for i in range(len(tp_vararg))]

        # Fold keyword arguments and resolve default argument values
        pysig = signature.pysig
        if pysig is None:
            if kw_args:
                raise NotImplementedError("unsupported keyword arguments "
                                          "when calling %s" % (fnty,))
            argvals = [self._cast_var(var, sigty)
                       for var, sigty in zip(pos_args, signature.args)]
        else:
            def normal_handler(index, param, var):
                return self._cast_var(var, signature.args[index])

            def default_handler(index, param, default):
                return self.context.get_constant_generic(
                    self.builder, signature.args[index], default)

            def stararg_handler(index, param, vars):
                stararg_ty = signature.args[index]
                assert isinstance(stararg_ty, types.BaseTuple), stararg_ty
                values = [self._cast_var(var, sigty)
                          for var, sigty in zip(vars, stararg_ty)]
                return cgutils.make_anonymous_struct(self.builder, values)

            argvals = typing.fold_arguments(pysig,
                                            pos_args, dict(kw_args),
                                            normal_handler,
                                            default_handler,
                                            stararg_handler)
        return argvals

    def lower_print(self, inst):
        """
        Lower a ir.Print()
        """
        # We handle this, as far as possible, as a normal call to built-in
        # print().  This will make it easy to undo the special ir.Print
        # rewrite when it becomes unnecessary (e.g. when we have native
        # strings).
        sig = self.fndesc.calltypes[inst]
        assert sig.return_type == types.none
        fnty = self.context.typing_context.resolve_value_type(print)

        # Fix the call signature to inject any constant-inferred
        # string argument
        pos_tys = list(sig.args)
        pos_args = list(inst.args)
        for i in range(len(pos_args)):
            if i in inst.consts:
                pyval = inst.consts[i]
                if isinstance(pyval, str):
                    pos_tys[i] = types.literal(pyval)

        fixed_sig = typing.signature(sig.return_type, *pos_tys)
        fixed_sig = fixed_sig.replace(pysig=sig.pysig)

        argvals = self.fold_call_args(fnty, sig, pos_args, inst.vararg, {})
        impl = self.context.get_function(print, fixed_sig)
        impl(self.builder, argvals)

    def lower_call(self, resty, expr):
        signature = self.fndesc.calltypes[expr]
        self.debug_print("# lower_call: expr = {0}".format(expr))
        if isinstance(signature.return_type, types.Phantom):
            return self.context.get_dummy_value()

        fnty = self.typeof(expr.func.name)

        if isinstance(fnty, types.ObjModeDispatcher):
            res = self._lower_call_ObjModeDispatcher(fnty, expr, signature)

        elif isinstance(fnty, types.ExternalFunction):
            res = self._lower_call_ExternalFunction(fnty, expr, signature)

        elif isinstance(fnty, types.ExternalFunctionPointer):
            res = self._lower_call_ExternalFunctionPointer(
                fnty, expr, signature)

        elif isinstance(fnty, types.RecursiveCall):
            res = self._lower_call_RecursiveCall(fnty, expr, signature)

        elif isinstance(fnty, types.FunctionType):
            res = self._lower_call_FunctionType(fnty, expr, signature)

        else:
            res = self._lower_call_normal(fnty, expr, signature)

        # If lowering the call returned None, interpret that as returning dummy
        # value if the return type of the function is void, otherwise there is
        # a problem
        if res is None:
            if signature.return_type == types.void:
                res = self.context.get_dummy_value()
            else:
                raise LoweringError(
                    msg="non-void function returns None from implementation",
                    loc=self.loc
                )

        return self.context.cast(self.builder, res, signature.return_type,
                                 resty)

    def _lower_call_ObjModeDispatcher(self, fnty, expr, signature):
        from numba.core.pythonapi import ObjModeUtils

        self.init_pyapi()
        # Acquire the GIL
        gil_state = self.pyapi.gil_ensure()
        # Fix types
        argnames = [a.name for a in expr.args]
        argtypes = [self.typeof(a) for a in argnames]
        argvalues = [self.loadvar(a) for a in argnames]
        for v, ty in zip(argvalues, argtypes):
            # Because .from_native_value steal the reference
            self.incref(ty, v)

        argobjs = [self.pyapi.from_native_value(atyp, aval,
                                                self.env_manager)
                   for atyp, aval in zip(argtypes, argvalues)]

        # Load objmode dispatcher
        callee = ObjModeUtils(self.pyapi).load_dispatcher(fnty, argtypes)
        # Make Call
        ret_obj = self.pyapi.call_function_objargs(callee, argobjs)
        has_exception = cgutils.is_null(self.builder, ret_obj)
        with self. builder.if_else(has_exception) as (then, orelse):
            # Handles exception
            # This branch must exit the function
            with then:
                # Clean arg
                for obj in argobjs:
                    self.pyapi.decref(obj)

                # Release the GIL
                self.pyapi.gil_release(gil_state)

                # Return and signal exception
                self.call_conv.return_exc(self.builder)

            # Handles normal return
            with orelse:
                # Fix output value
                native = self.pyapi.to_native_value(
                    fnty.dispatcher.output_types,
                    ret_obj,
                )
                output = native.value

                # Release objs
                self.pyapi.decref(ret_obj)
                for obj in argobjs:
                    self.pyapi.decref(obj)

                # cleanup output
                if callable(native.cleanup):
                    native.cleanup()

                # Release the GIL
                self.pyapi.gil_release(gil_state)

                # Error during unboxing
                with self.builder.if_then(native.is_error):
                    self.call_conv.return_exc(self.builder)

                return output

    def _lower_call_ExternalFunction(self, fnty, expr, signature):
        # Handle a named external function
        self.debug_print("# external function")
        argvals = self.fold_call_args(
            fnty, signature, expr.args, expr.vararg, expr.kws,
        )
        fndesc = funcdesc.ExternalFunctionDescriptor(
            fnty.symbol, fnty.sig.return_type, fnty.sig.args)
        func = self.context.declare_external_function(
            self.builder.module, fndesc)
        return self.context.call_external_function(
            self.builder, func, fndesc.argtypes, argvals,
        )

    def _lower_call_ExternalFunctionPointer(self, fnty, expr, signature):
        # Handle a C function pointer
        self.debug_print("# calling external function pointer")
        argvals = self.fold_call_args(
            fnty, signature, expr.args, expr.vararg, expr.kws,
        )
        pointer = self.loadvar(expr.func.name)
        # If the external function pointer uses libpython
        if fnty.requires_gil:
            self.init_pyapi()
            # Acquire the GIL
            gil_state = self.pyapi.gil_ensure()
            # Make PyObjects
            newargvals = []
            pyvals = []
            for exptyp, gottyp, aval in zip(fnty.sig.args, signature.args,
                                            argvals):
                # Adjust argument values to pyobjects
                if exptyp == types.ffi_forced_object:
                    self.incref(gottyp, aval)
                    obj = self.pyapi.from_native_value(
                        gottyp, aval, self.env_manager,
                    )
                    newargvals.append(obj)
                    pyvals.append(obj)
                else:
                    newargvals.append(aval)

            # Call external function
            res = self.context.call_function_pointer(
                self.builder, pointer, newargvals, fnty.cconv,
            )
            # Release PyObjects
            for obj in pyvals:
                self.pyapi.decref(obj)

            # Release the GIL
            self.pyapi.gil_release(gil_state)
        # If the external function pointer does NOT use libpython
        else:
            res = self.context.call_function_pointer(
                self.builder, pointer, argvals, fnty.cconv,
            )
        return res

    def _lower_call_RecursiveCall(self, fnty, expr, signature):
        # Recursive call
        argvals = self.fold_call_args(
            fnty, signature, expr.args, expr.vararg, expr.kws,
        )
        rec_ov = fnty.get_overloads(signature.args)
        mangler = self.context.mangler or default_mangler
        abi_tags = self.fndesc.abi_tags
        mangled_name = mangler(rec_ov.qualname, signature.args,
                               abi_tags=abi_tags, uid=rec_ov.uid)
        # special case self recursion
        if self.builder.function.name.startswith(mangled_name):
            res = self.context.call_internal(
                self.builder, self.fndesc, signature, argvals,
            )
        else:
            res = self.context.call_unresolved(
                self.builder, mangled_name, signature, argvals,
            )
        return res

    def _lower_call_FunctionType(self, fnty, expr, signature):
        self.debug_print("# calling first-class function type")
        sig = types.unliteral(signature)
        if not fnty.check_signature(signature):
            # value dependent polymorphism?
            raise UnsupportedError(
                f'mismatch of function types:'
                f' expected {fnty} but got {types.FunctionType(sig)}')
        ftype = fnty.ftype
        argvals = self.fold_call_args(
            fnty, sig, expr.args, expr.vararg, expr.kws,
        )
        func_ptr = self.__get_function_pointer(ftype, expr.func.name, sig=sig)
        res = self.builder.call(func_ptr, argvals, cconv=fnty.cconv)
        return res

    def __get_function_pointer(self, ftype, fname, sig=None):
        from numba.experimental.function_type import lower_get_wrapper_address

        llty = self.context.get_value_type(ftype)
        fstruct = self.loadvar(fname)
        addr = self.builder.extract_value(fstruct, 0,
                                          name='addr_of_%s' % (fname))

        fptr = cgutils.alloca_once(self.builder, llty,
                                   name="fptr_of_%s" % (fname))
        with self.builder.if_else(
                cgutils.is_null(self.builder, addr),
                likely=False) as (then, orelse):
            with then:
                self.init_pyapi()
                # Acquire the GIL
                gil_state = self.pyapi.gil_ensure()
                pyaddr = self.builder.extract_value(
                    fstruct, 1,
                    name='pyaddr_of_%s' % (fname))
                # try to recover the function address, see
                # test_zero_address BadToGood example in
                # test_function_type.py
                addr1 = lower_get_wrapper_address(
                    self.context, self.builder, pyaddr, sig,
                    failure_mode='ignore')
                with self.builder.if_then(
                        cgutils.is_null(self.builder, addr1), likely=False):
                    self.return_exception(
                        RuntimeError,
                        exc_args=(f"{ftype} function address is null",),
                        loc=self.loc)
                addr2 = self.pyapi.long_as_voidptr(addr1)
                self.builder.store(self.builder.bitcast(addr2, llty), fptr)
                self.pyapi.decref(addr1)
                self.pyapi.gil_release(gil_state)
            with orelse:
                self.builder.store(self.builder.bitcast(addr, llty), fptr)
        return self.builder.load(fptr)

    def _lower_call_normal(self, fnty, expr, signature):
        # Normal function resolution
        self.debug_print("# calling normal function: {0}".format(fnty))
        self.debug_print("# signature: {0}".format(signature))
        if isinstance(fnty, types.ObjModeDispatcher):
            argvals = expr.func.args
        else:
            argvals = self.fold_call_args(
                fnty, signature, expr.args, expr.vararg, expr.kws,
            )
        tname = expr.target
        if tname is not None:
            from numba.core.target_extension import resolve_dispatcher_from_str
            disp = resolve_dispatcher_from_str(tname)
            hw_ctx = disp.targetdescr.target_context
            impl = hw_ctx.get_function(fnty, signature)
        else:
            impl = self.context.get_function(fnty, signature)
        if signature.recvr:
            # The "self" object is passed as the function object
            # for bounded function
            the_self = self.loadvar(expr.func.name)
            # Prepend the self reference
            argvals = [the_self] + list(argvals)

        res = impl(self.builder, argvals, self.loc)
        return res

    def lower_expr(self, resty, expr):
        if expr.op == 'binop':
            return self.lower_binop(resty, expr, expr.fn)
        elif expr.op == 'inplace_binop':
            lty = self.typeof(expr.lhs.name)
            if lty.mutable:
                return self.lower_binop(resty, expr, expr.fn)
            else:
                # inplace operators on non-mutable types reuse the same
                # definition as the corresponding copying operators.)
                return self.lower_binop(resty, expr, expr.immutable_fn)
        elif expr.op == 'unary':
            val = self.loadvar(expr.value.name)
            typ = self.typeof(expr.value.name)
            func_ty = self.context.typing_context.resolve_value_type(expr.fn)
            # Get function
            signature = self.fndesc.calltypes[expr]
            impl = self.context.get_function(func_ty, signature)
            # Convert argument to match
            val = self.context.cast(self.builder, val, typ, signature.args[0])
            res = impl(self.builder, [val])
            res = self.context.cast(self.builder, res,
                                    signature.return_type, resty)
            return res

        elif expr.op == 'call':
            res = self.lower_call(resty, expr)
            return res

        elif expr.op == 'pair_first':
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)
            res = self.context.pair_first(self.builder, val, ty)
            self.incref(resty, res)
            return res

        elif expr.op == 'pair_second':
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)
            res = self.context.pair_second(self.builder, val, ty)
            self.incref(resty, res)
            return res

        elif expr.op in ('getiter', 'iternext'):
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)
            signature = self.fndesc.calltypes[expr]
            impl = self.context.get_function(expr.op, signature)
            [fty] = signature.args
            castval = self.context.cast(self.builder, val, ty, fty)
            res = impl(self.builder, (castval,))
            res = self.context.cast(self.builder, res, signature.return_type,
                                    resty)
            return res

        elif expr.op == 'exhaust_iter':
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)
            # Unpack optional
            if isinstance(ty, types.Optional):
                val = self.context.cast(self.builder, val, ty, ty.type)
                ty = ty.type

            # If we have a tuple, we needn't do anything
            # (and we can't iterate over the heterogeneous ones).
            if isinstance(ty, types.BaseTuple):
                assert ty == resty
                self.incref(ty, val)
                return val

            itemty = ty.iterator_type.yield_type
            tup = self.context.get_constant_undef(resty)
            pairty = types.Pair(itemty, types.boolean)
            getiter_sig = typing.signature(ty.iterator_type, ty)
            getiter_impl = self.context.get_function('getiter',
                                                     getiter_sig)
            iternext_sig = typing.signature(pairty, ty.iterator_type)
            iternext_impl = self.context.get_function('iternext',
                                                      iternext_sig)
            iterobj = getiter_impl(self.builder, (val,))
            # We call iternext() as many times as desired (`expr.count`).
            for i in range(expr.count):
                pair = iternext_impl(self.builder, (iterobj,))
                is_valid = self.context.pair_second(self.builder,
                                                    pair, pairty)
                with cgutils.if_unlikely(self.builder,
                                         self.builder.not_(is_valid)):
                    self.return_exception(ValueError, loc=self.loc)
                item = self.context.pair_first(self.builder,
                                               pair, pairty)
                tup = self.builder.insert_value(tup, item, i)

            # Call iternext() once more to check that the iterator
            # is exhausted.
            pair = iternext_impl(self.builder, (iterobj,))
            is_valid = self.context.pair_second(self.builder,
                                                pair, pairty)
            with cgutils.if_unlikely(self.builder, is_valid):
                self.return_exception(ValueError, loc=self.loc)

            self.decref(ty.iterator_type, iterobj)
            return tup

        elif expr.op == "getattr":
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)

            if isinstance(resty, types.BoundFunction):
                # if we are getting out a method, assume we have typed this
                # properly and just build a bound function object
                casted = self.context.cast(self.builder, val, ty, resty.this)
                res = self.context.get_bound_function(self.builder, casted,
                                                      resty.this)
                self.incref(resty, res)
                return res
            else:
                impl = self.context.get_getattr(ty, expr.attr)
                attrty = self.context.typing_context.resolve_getattr(ty,
                                                                     expr.attr)

                if impl is None:
                    # ignore the attribute
                    return self.context.get_dummy_value()
                else:
                    res = impl(self.context, self.builder, ty, val, expr.attr)

                # Cast the attribute type to the expected output type
                res = self.context.cast(self.builder, res, attrty, resty)
                return res

        elif expr.op == "static_getitem":
            signature = typing.signature(
                resty,
                self.typeof(expr.value.name),
                _lit_or_omitted(expr.index),
            )
            try:
                # Both get_function() and the returned implementation can
                # raise NotImplementedError if the types aren't supported
                impl = self.context.get_function("static_getitem", signature)
                return impl(self.builder,
                            (self.loadvar(expr.value.name), expr.index))
            except NotImplementedError:
                if expr.index_var is None:
                    raise
                # Fall back on the generic getitem() implementation
                # for this type.
                signature = self.fndesc.calltypes[expr]
                return self.lower_getitem(resty, expr, expr.value,
                                          expr.index_var, signature)
        elif expr.op == "typed_getitem":
            signature = typing.signature(
                resty,
                self.typeof(expr.value.name),
                self.typeof(expr.index.name),
            )
            impl = self.context.get_function("typed_getitem", signature)
            return impl(self.builder, (self.loadvar(expr.value.name),
                        self.loadvar(expr.index.name)))
        elif expr.op == "getitem":
            signature = self.fndesc.calltypes[expr]
            return self.lower_getitem(resty, expr, expr.value, expr.index,
                                      signature)

        elif expr.op == "build_tuple":
            itemvals = [self.loadvar(i.name) for i in expr.items]
            itemtys = [self.typeof(i.name) for i in expr.items]
            castvals = [self.context.cast(self.builder, val, fromty, toty)
                        for val, toty, fromty in zip(itemvals, resty, itemtys)]
            tup = self.context.make_tuple(self.builder, resty, castvals)
            self.incref(resty, tup)
            return tup

        elif expr.op == "build_list":
            itemvals = [self.loadvar(i.name) for i in expr.items]
            itemtys = [self.typeof(i.name) for i in expr.items]
            if isinstance(resty, types.LiteralList):
                castvals = [self.context.cast(self.builder, val, fromty, toty)
                            for val, toty, fromty in zip(itemvals, resty.types,
                                                         itemtys)]
                tup = self.context.make_tuple(self.builder,
                                              types.Tuple(resty.types),
                                              castvals)
                self.incref(resty, tup)
                return tup
            else:
                castvals = [self.context.cast(self.builder, val, fromty,
                                              resty.dtype)
                            for val, fromty in zip(itemvals, itemtys)]
                return self.context.build_list(self.builder, resty, castvals)

        elif expr.op == "build_set":
            # Insert in reverse order, as Python does
            items = expr.items[::-1]
            itemvals = [self.loadvar(i.name) for i in items]
            itemtys = [self.typeof(i.name) for i in items]
            castvals = [self.context.cast(self.builder, val, fromty,
                                          resty.dtype)
                        for val, fromty in zip(itemvals, itemtys)]
            return self.context.build_set(self.builder, resty, castvals)

        elif expr.op == "build_map":
            items = expr.items
            keys, values = [], []
            key_types, value_types = [], []
            for k, v in items:
                key = self.loadvar(k.name)
                keytype = self.typeof(k.name)
                val = self.loadvar(v.name)
                valtype = self.typeof(v.name)
                keys.append(key)
                values.append(val)
                key_types.append(keytype)
                value_types.append(valtype)
            return self.context.build_map(self.builder, resty,
                                          list(zip(key_types, value_types)),
                                          list(zip(keys, values)))

        elif expr.op == "cast":
            val = self.loadvar(expr.value.name)
            ty = self.typeof(expr.value.name)
            castval = self.context.cast(self.builder, val, ty, resty)
            self.incref(resty, castval)
            return castval

        elif expr.op == "phi":
            raise LoweringError("PHI not stripped")

        elif expr.op == 'null':
            return self.context.get_constant_null(resty)

        elif expr.op == 'undef':
            # Numba does not raise an UnboundLocalError for undefined variables.
            # The variable is set to zero.
            return self.context.get_constant_null(resty)

        elif expr.op in self.context.special_ops:
            res = self.context.special_ops[expr.op](self, expr)
            return res

        raise NotImplementedError(expr)

    def _alloca_var(self, name, fetype):
        """
        Ensure the given variable has an allocated stack slot (if needed).
        """
        if name in self.varmap:
            # quit early
            return

        # If the name is used in multiple blocks or lowering with debuginfo...
        if ((name not in self._singly_assigned_vars) or
                self._disable_sroa_like_opt):
            # If not already defined, allocate it
            ptr = self.alloca(name, fetype)
            # Remember the pointer
            self.varmap[name] = ptr

    def getvar(self, name):
        """
        Get a pointer to the given variable's slot.
        """
        if not self._disable_sroa_like_opt:
            assert name not in self._blk_local_varmap
            assert name not in self._singly_assigned_vars
        if name not in self.varmap:
            # Allocate undefined variable as needed.
            # NOTE: Py3.12 use of LOAD_FAST_AND_CLEAR will allow variable be
            # referenced before it is defined.
            self._alloca_var(name, self.typeof(name))
        return self.varmap[name]

    def loadvar(self, name):
        """
        Load the given variable's value.
        """
        if name in self._blk_local_varmap and not self._disable_sroa_like_opt:
            return self._blk_local_varmap[name]
        ptr = self.getvar(name)

        # Don't associate debuginfo with the load for a function arg else it
        # creates instructions ahead of the first source line of the
        # function which then causes problems with breaking on the function
        # symbol (it hits the symbol, not the first line).
        if name in self.func_ir.arg_names:
            with debuginfo.suspend_emission(self.builder):
                return self.builder.load(ptr)
        else:
            return self.builder.load(ptr)

    def storevar(self, value, name, argidx=None):
        """
        Store the value into the given variable.
        """
        fetype = self.typeof(name)
        # Define if not already
        self._alloca_var(name, fetype)

        # Store variable
        if (name in self._singly_assigned_vars and
                not self._disable_sroa_like_opt):
            self._blk_local_varmap[name] = value
        else:
            if argidx is None:
                # Clean up existing value stored in the variable, not needed
                # if it's an arg
                old = self.loadvar(name)
                self.decref(fetype, old)

            # stack stored variable
            ptr = self.getvar(name)
            if value.type != ptr.type.pointee:
                msg = ("Storing {value.type} to ptr of {ptr.type.pointee} "
                       "('{name}'). FE type {fetype}").format(value=value,
                                                              ptr=ptr,
                                                              fetype=fetype,
                                                              name=name)
                raise AssertionError(msg)

            # If this store is associated with an argument to the function (i.e.
            # store following reassemble from CC splatting structs as many args
            # to the function) then mark this variable as such.
            if argidx is not None:
                with debuginfo.suspend_emission(self.builder):
                    self.builder.store(value, ptr)
                loc = self.defn_loc # the line with `def <func>`
                lltype = self.context.get_value_type(fetype)
                sizeof = self.context.get_abi_sizeof(lltype)
                datamodel = self.context.data_model_manager[fetype]
                self.debuginfo.mark_variable(self.builder, ptr, name=name,
                                             lltype=lltype, size=sizeof,
                                             line=loc.line, datamodel=datamodel,
                                             argidx=argidx)
            else:
                self.builder.store(value, ptr)

    def delvar(self, name):
        """
        Delete the given variable.
        """
        fetype = self.typeof(name)

        # Out-of-order
        if (name not in self._blk_local_varmap and
                not self._disable_sroa_like_opt):
            if name in self._singly_assigned_vars:
                self._singly_assigned_vars.discard(name)

        # Define if not already (may happen if the variable is deleted
        # at the beginning of a loop, but only set later in the loop)
        self._alloca_var(name, fetype)

        if name in self._blk_local_varmap and not self._disable_sroa_like_opt:
            llval = self._blk_local_varmap[name]
            self.decref(fetype, llval)
        else:
            ptr = self.getvar(name)
            self.decref(fetype, self.builder.load(ptr))
            # Zero-fill variable to avoid double frees on subsequent dels
            self.builder.store(Constant(ptr.type.pointee, None), ptr)

    def alloca(self, name, type):
        lltype = self.context.get_value_type(type)
        datamodel = self.context.data_model_manager[type]
        return self.alloca_lltype(name, lltype, datamodel=datamodel)

    def alloca_lltype(self, name, lltype, datamodel=None):
        # Is user variable?
        is_uservar = not name.startswith('$')
        # Allocate space for variable
        aptr = cgutils.alloca_once(self.builder, lltype,
                                   name=name, zfill=False)

        # Emit debug info for user variable
        if is_uservar:
            # Don't associate debuginfo with the alloca for a function arg, this
            # is handled by the first store to the alloca so that repacking the
            # splatted args from the CC is dealt with.
            if name not in self.func_ir.arg_names:
                sizeof = self.context.get_abi_sizeof(lltype)
                self.debuginfo.mark_variable(self.builder, aptr, name=name,
                                             lltype=lltype, size=sizeof,
                                             line=self.loc.line,
                                             datamodel=datamodel,)
        return aptr

    def incref(self, typ, val):
        if not self.context.enable_nrt:
            return

        self.context.nrt.incref(self.builder, typ, val)

    def decref(self, typ, val):
        if not self.context.enable_nrt:
            return

        # do not associate decref with "use", it creates "jumpy" line info as
        # the decrefs are usually where the ir.Del nodes are, which is at the
        # end of the block.
        with debuginfo.suspend_emission(self.builder):
            self.context.nrt.decref(self.builder, typ, val)


def _lit_or_omitted(value):
    """Returns a Literal instance if the type of value is supported;
    otherwise, return `Omitted(value)`.
    """
    try:
        return types.literal(value)
    except LiteralTypingError:
        return types.Omitted(value)
