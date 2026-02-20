import numba
from numba.core import errors, ir, ir_utils, sigutils, types
from numba.core.ir_utils import build_definitions
from numba.core.transforms import find_region_inout_vars
from numba.core.typing.typeof import typeof_impl


class WithContext(object):
    """A dummy object for use as contextmanager.
    This can be used as a contextmanager.
    """

    is_callable = False

    def __enter__(self):
        pass

    def __exit__(self, typ, val, tb):
        pass

    def mutate_with_body(
        self,
        func_ir,
        blocks,
        blk_start,
        blk_end,
        body_blocks,
        dispatcher_factory,
        extra,
    ):
        """Mutate the *blocks* to implement this contextmanager.

        Parameters
        ----------
        func_ir : FunctionIR
        blocks : dict[ir.Block]
        blk_start, blk_end : int
            labels of the starting and ending block of the context-manager.
        body_block: sequence[int]
            A sequence of int's representing labels of the with-body
        dispatcher_factory : callable
            A callable that takes a `FunctionIR` and returns a `Dispatcher`.
        """
        raise NotImplementedError


@typeof_impl.register(WithContext)
def typeof_contextmanager(val, c):
    return types.ContextManager(val)


def _get_var_parent(name):
    """Get parent of the variable given its name"""
    # If not a temporary variable
    if not name.startswith("$"):
        # Return the base component of the name
        return name.split(
            ".",
        )[0]


def _clear_blocks(blocks, to_clear):
    """Remove keys in *to_clear* from *blocks*."""
    for b in to_clear:
        del blocks[b]


class _ByPassContextType(WithContext):
    """A simple context-manager that tells the compiler to bypass the body
    of the with-block.
    """

    def mutate_with_body(
        self,
        func_ir,
        blocks,
        blk_start,
        blk_end,
        body_blocks,
        dispatcher_factory,
        extra,
    ):
        assert extra is None
        # Determine variables that need forwarding
        vlt = func_ir.variable_lifetime
        inmap = {_get_var_parent(k): k for k in vlt.livemap[blk_start]}
        outmap = {_get_var_parent(k): k for k in vlt.livemap[blk_end]}
        forwardvars = {inmap[k]: outmap[k] for k in filter(bool, outmap)}
        # Transform the block
        _bypass_with_context(blocks, blk_start, blk_end, forwardvars)
        _clear_blocks(blocks, body_blocks)


bypass_context = _ByPassContextType()


class _CallContextType(WithContext):
    """A simple context-manager that tells the compiler to lift the body of the
    with-block as another function.
    """

    def mutate_with_body(
        self,
        func_ir,
        blocks,
        blk_start,
        blk_end,
        body_blocks,
        dispatcher_factory,
        extra,
    ):
        assert extra is None
        vlt = func_ir.variable_lifetime

        inputs, outputs = find_region_inout_vars(
            blocks=blocks,
            livemap=vlt.livemap,
            callfrom=blk_start,
            returnto=blk_end,
            body_block_ids=set(body_blocks),
        )

        lifted_blks = {k: blocks[k] for k in body_blocks}
        _mutate_with_block_callee(
            lifted_blks, blk_start, blk_end, inputs, outputs
        )

        # XXX: transform body-blocks to return the output variables
        lifted_ir = func_ir.derive(
            blocks=lifted_blks,
            arg_names=tuple(inputs),
            arg_count=len(inputs),
            force_non_generator=True,
        )

        dispatcher = dispatcher_factory(lifted_ir)

        newblk = _mutate_with_block_caller(
            dispatcher,
            blocks,
            blk_start,
            blk_end,
            inputs,
            outputs,
        )

        blocks[blk_start] = newblk
        _clear_blocks(blocks, body_blocks)
        return dispatcher


call_context = _CallContextType()


class _ObjModeContextType(WithContext):
    """Creates a contextmanager to be used inside jitted functions to enter
    *object-mode* for using interpreter features.  The body of the with-context
    is lifted into a function that is compiled in *object-mode*.  This
    transformation process is limited and cannot process all possible
    Python code.  However, users can wrap complicated logic in another
    Python function, which will then be executed by the interpreter.

    Use this as a function that takes keyword arguments only.
    The argument names must correspond to the output variables from the
    with-block.  Their respective values can be:

    1. strings representing the expected types; i.e. ``"float32"``.
    2. compile-time bound global or nonlocal variables referring to the
       expected type. The variables are read at compile time.

    When exiting the with-context, the output variables are converted
    to the expected nopython types according to the annotation.  This process
    is the same as passing Python objects into arguments of a nopython
    function.

    Example::

        import numpy as np
        from numba import njit, objmode, types

        def bar(x):
            # This code is executed by the interpreter.
            return np.asarray(list(reversed(x.tolist())))

        # Output type as global variable
        out_ty = types.intp[:]

        @njit
        def foo():
            x = np.arange(5)
            y = np.zeros_like(x)
            with objmode(y='intp[:]', z=out_ty):  # annotate return type
                # this region is executed by object-mode.
                y += bar(x)
                z = y
            return y, z

    .. note:: Known limitations:

        - with-block cannot use incoming list objects.
        - with-block cannot use incoming function objects.
        - with-block cannot ``yield``, ``break``, ``return`` or ``raise`` \
          such that the execution will leave the with-block immediately.
        - with-block cannot contain `with` statements.
        - random number generator states do not synchronize; i.e. \
          nopython-mode and object-mode uses different RNG states.

    .. note:: When used outside of no-python mode, the context-manager has no
        effect.

    .. warning:: This feature is experimental.  The supported features may
        change with or without notice.

    """

    is_callable = True

    def _legalize_args(
        self, func_ir, args, kwargs, loc, func_globals, func_closures
    ):
        """
        Legalize arguments to the context-manager

        Parameters
        ----------
        func_ir: FunctionIR
        args: tuple
            Positional arguments to the with-context call as IR nodes.
        kwargs: dict
            Keyword arguments to the with-context call as IR nodes.
        loc: numba.core.ir.Loc
            Source location of the with-context call.
        func_globals: dict
            The globals dictionary of the calling function.
        func_closures: dict
            The resolved closure variables of the calling function.
        """
        if args:
            raise errors.CompilerError(
                "objectmode context doesn't take any positional arguments",
            )
        typeanns = {}

        def report_error(varname, msg, loc):
            raise errors.CompilerError(
                f"Error handling objmode argument {varname!r}. {msg}",
                loc=loc,
            )

        for k, v in kwargs.items():
            if isinstance(v, ir.Const) and isinstance(v.value, str):
                typeanns[k] = sigutils._parse_signature_string(v.value)
            elif isinstance(v, ir.FreeVar):
                try:
                    v = func_closures[v.name]
                except KeyError:
                    report_error(
                        varname=k,
                        msg=f"Freevar {v.name!r} is not defined.",
                        loc=loc,
                    )
                typeanns[k] = v
            elif isinstance(v, ir.Global):
                try:
                    v = func_globals[v.name]
                except KeyError:
                    report_error(
                        varname=k,
                        msg=f"Global {v.name!r} is not defined.",
                        loc=loc,
                    )
                typeanns[k] = v
            elif isinstance(v, ir.Expr) and v.op == "getattr":
                try:
                    base_obj = func_ir.infer_constant(v.value)
                    typ = getattr(base_obj, v.attr)
                except (errors.ConstantInferenceError, AttributeError):
                    report_error(
                        varname=k,
                        msg="Getattr cannot be resolved at compile-time.",
                        loc=loc,
                    )
                else:
                    typeanns[k] = typ
            else:
                report_error(
                    varname=k,
                    msg=(
                        "The value must be a compile-time constant either as "
                        "a non-local variable or a getattr expression that "
                        "refers to a Numba type."
                    ),
                    loc=loc,
                )

        # Legalize the types for objmode
        for name, typ in typeanns.items():
            self._legalize_arg_type(name, typ, loc)

        return typeanns

    def _legalize_arg_type(self, name, typ, loc):
        """Legalize the argument type

        Parameters
        ----------
        name: str
            argument name.
        typ: numba.core.types.Type
            argument type.
        loc: numba.core.ir.Loc
            source location for error reporting.
        """
        if getattr(typ, "reflected", False):
            msgbuf = [
                "Objmode context failed.",
                f"Argument {name!r} is declared as "
                f"an unsupported type: {typ}.",
                "Reflected types are not supported.",
            ]
            raise errors.CompilerError(" ".join(msgbuf), loc=loc)

    def mutate_with_body(
        self,
        func_ir,
        blocks,
        blk_start,
        blk_end,
        body_blocks,
        dispatcher_factory,
        extra,
    ):
        cellnames = func_ir.func_id.func.__code__.co_freevars
        closures = func_ir.func_id.func.__closure__
        func_globals = func_ir.func_id.func.__globals__
        if closures is not None:
            # Resolve free variables
            func_closures = {}
            for cellname, closure in zip(cellnames, closures):
                try:
                    cellval = closure.cell_contents
                except ValueError as e:
                    # empty cell will raise
                    if str(e) != "Cell is empty":
                        raise
                else:
                    func_closures[cellname] = cellval
        else:
            # Missing closure object
            func_closures = {}
        args = extra["args"] if extra else ()
        kwargs = extra["kwargs"] if extra else {}

        typeanns = self._legalize_args(
            func_ir=func_ir,
            args=args,
            kwargs=kwargs,
            loc=blocks[blk_start].loc,
            func_globals=func_globals,
            func_closures=func_closures,
        )
        vlt = func_ir.variable_lifetime

        inputs, outputs = find_region_inout_vars(
            blocks=blocks,
            livemap=vlt.livemap,
            callfrom=blk_start,
            returnto=blk_end,
            body_block_ids=set(body_blocks),
        )

        # Determine types in the output tuple
        def strip_var_ver(x):
            return x.split(".", 1)[0]

        stripped_outs = list(map(strip_var_ver, outputs))

        # Verify that only outputs are annotated
        extra_annotated = set(typeanns) - set(stripped_outs)
        if extra_annotated:
            msg = (
                "Invalid type annotation on non-outgoing variables: {}."
                "Suggestion: remove annotation of the listed variables"
            )
            raise errors.TypingError(msg.format(extra_annotated))

        # Verify that all outputs are annotated

        # Note on "$cp" variable:
        # ``transforms.consolidate_multi_exit_withs()`` introduces the variable
        # for the control-point to determine the correct exit block. This
        # variable crosses the with-region boundary. Thus, it will be consider
        # an output variable leaving the lifted with-region.
        typeanns["$cp"] = types.int32
        not_annotated = set(stripped_outs) - set(typeanns)
        if not_annotated:
            msg = (
                "Missing type annotation on outgoing variable(s): {0}\n\n"
                "Example code: with objmode({1}='<"
                "add_type_as_string_here>')\n"
            )
            stable_ann = sorted(not_annotated)
            raise errors.TypingError(msg.format(stable_ann, stable_ann[0]))

        # Get output types
        outtup = types.Tuple([typeanns[v] for v in stripped_outs])

        lifted_blks = {k: blocks[k] for k in body_blocks}
        _mutate_with_block_callee(
            lifted_blks, blk_start, blk_end, inputs, outputs
        )

        lifted_ir = func_ir.derive(
            blocks=lifted_blks,
            arg_names=tuple(inputs),
            arg_count=len(inputs),
            force_non_generator=True,
        )

        dispatcher = dispatcher_factory(
            lifted_ir, objectmode=True, output_types=outtup
        )

        newblk = _mutate_with_block_caller(
            dispatcher,
            blocks,
            blk_start,
            blk_end,
            inputs,
            outputs,
        )

        blocks[blk_start] = newblk
        _clear_blocks(blocks, body_blocks)
        return dispatcher

    def __call__(self, *args, **kwargs):
        # No effect when used in pure-python
        return self


objmode_context = _ObjModeContextType()


def _bypass_with_context(blocks, blk_start, blk_end, forwardvars):
    """Given the starting and ending block of the with-context,
    replaces the head block with a new block that jumps to the end.

    *blocks* is modified inplace.
    """
    sblk = blocks[blk_start]
    scope = sblk.scope
    loc = sblk.loc
    newblk = ir.Block(scope=scope, loc=loc)
    for k, v in forwardvars.items():
        newblk.append(
            ir.Assign(
                value=scope.get_exact(k), target=scope.get_exact(v), loc=loc
            )
        )
    newblk.append(ir.Jump(target=blk_end, loc=loc))
    blocks[blk_start] = newblk


def _mutate_with_block_caller(
    dispatcher, blocks, blk_start, blk_end, inputs, outputs
):
    """Make a new block that calls into the lifeted with-context.

    Parameters
    ----------
    dispatcher : Dispatcher
    blocks : dict[ir.Block]
    blk_start, blk_end : int
        labels of the starting and ending block of the context-manager.
    inputs: sequence[str]
        Input variable names
    outputs: sequence[str]
        Output variable names
    """
    sblk = blocks[blk_start]
    scope = sblk.scope
    loc = sblk.loc
    newblock = ir.Block(scope=scope, loc=loc)

    ir_utils.fill_block_with_call(
        newblock=newblock,
        callee=dispatcher,
        label_next=blk_end,
        inputs=inputs,
        outputs=outputs,
    )
    return newblock


def _mutate_with_block_callee(blocks, blk_start, blk_end, inputs, outputs):
    """Mutate *blocks* for the callee of a with-context.

    Parameters
    ----------
    blocks : dict[ir.Block]
    blk_start, blk_end : int
        labels of the starting and ending block of the context-manager.
    inputs: sequence[str]
        Input variable names
    outputs: sequence[str]
        Output variable names
    """
    if not blocks:
        raise errors.NumbaValueError("No blocks in with-context block")
    head_blk = min(blocks)
    temp_blk = blocks[head_blk]
    scope = temp_blk.scope
    loc = temp_blk.loc

    blocks[blk_start] = ir_utils.fill_callee_prologue(
        block=ir.Block(scope=scope, loc=loc),
        inputs=inputs,
        label_next=head_blk,
    )
    blocks[blk_end] = ir_utils.fill_callee_epilogue(
        block=ir.Block(scope=scope, loc=loc),
        outputs=outputs,
    )


class _ParallelChunksize(WithContext):
    is_callable = True

    """A context-manager that on entry stores the current chunksize
    for the executing parfors and then changes the current chunksize
    to the programmer specified value. On exit the original
    chunksize is restored.
    """

    def mutate_with_body(
        self,
        func_ir,
        blocks,
        blk_start,
        blk_end,
        body_blocks,
        dispatcher_factory,
        extra,
    ):
        ir_utils.dprint_func_ir(func_ir, "Before with changes", blocks=blocks)
        assert extra is not None
        args = extra["args"]
        assert len(args) == 1
        arg = args[0]
        scope = blocks[blk_start].scope
        loc = blocks[blk_start].loc
        if isinstance(arg, ir.Arg):
            arg = ir.Var(scope, arg.name, loc)

        set_state = []
        restore_state = []

        # global for Numba itself
        gvar = scope.redefine("$ngvar", loc)
        set_state.append(ir.Assign(ir.Global("numba", numba, loc), gvar, loc))
        # getattr for set chunksize function in Numba
        spcattr = ir.Expr.getattr(gvar, "set_parallel_chunksize", loc)
        spcvar = scope.redefine("$spc", loc)
        set_state.append(ir.Assign(spcattr, spcvar, loc))
        # call set_parallel_chunksize
        orig_pc_var = scope.redefine("$save_pc", loc)
        cs_var = scope.redefine("$cs_var", loc)
        set_state.append(ir.Assign(arg, cs_var, loc))
        spc_call = ir.Expr.call(spcvar, [cs_var], (), loc)
        set_state.append(ir.Assign(spc_call, orig_pc_var, loc))

        restore_spc_call = ir.Expr.call(spcvar, [orig_pc_var], (), loc)
        restore_state.append(ir.Assign(restore_spc_call, orig_pc_var, loc))

        blocks[blk_start].body = (
            blocks[blk_start].body[1:-1]
            + set_state
            + [blocks[blk_start].body[-1]]
        )
        blocks[blk_end].body = restore_state + blocks[blk_end].body
        func_ir._definitions = build_definitions(blocks)
        ir_utils.dprint_func_ir(func_ir, "After with changes", blocks=blocks)

    def __call__(self, *args, **kwargs):
        """Act like a function and enforce the contract that
        setting the chunksize takes only one integer input.
        """
        if len(args) != 1 or kwargs or not isinstance(args[0], int):
            raise ValueError(
                "parallel_chunksize takes only a " "single integer argument."
            )

        self.chunksize = args[0]
        return self

    def __enter__(self):
        self.orig_chunksize = numba.get_parallel_chunksize()
        numba.set_parallel_chunksize(self.chunksize)

    def __exit__(self, typ, val, tb):
        numba.set_parallel_chunksize(self.orig_chunksize)


parallel_chunksize = _ParallelChunksize()
