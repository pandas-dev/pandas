from collections import namedtuple
import copy
import warnings
from numba.core.tracing import event

from numba.core import (errors, interpreter, bytecode, postproc, config,
                        callconv, cpu)
from numba.parfors.parfor import ParforDiagnostics
from numba.core.errors import CompilerError
from numba.core.environment import lookup_environment

from numba.core.compiler_machinery import PassManager

from numba.core.untyped_passes import (ExtractByteCode, TranslateByteCode,
                                       FixupArgs, IRProcessing, DeadBranchPrune,
                                       RewriteSemanticConstants,
                                       InlineClosureLikes, GenericRewrites,
                                       WithLifting, InlineInlinables,
                                       FindLiterallyCalls,
                                       MakeFunctionToJitFunction,
                                       CanonicalizeLoopExit,
                                       CanonicalizeLoopEntry, LiteralUnroll,
                                       ReconstructSSA, RewriteDynamicRaises,
                                       LiteralPropagationSubPipelinePass,
                                       )

from numba.core.typed_passes import (NopythonTypeInference, AnnotateTypes,
                                     NopythonRewrites, PreParforPass,
                                     ParforPass, DumpParforDiagnostics,
                                     IRLegalization, NoPythonBackend,
                                     InlineOverloads, PreLowerStripPhis,
                                     NativeLowering, NativeParforLowering,
                                     NoPythonSupportedFeatureValidation,
                                     ParforFusionPass, ParforPreLoweringPass
                                     )

from numba.core.object_mode_passes import (ObjectModeFrontEnd,
                                           ObjectModeBackEnd)
from numba.core.targetconfig import TargetConfig, Option, ConfigStack


class Flags(TargetConfig):
    __slots__ = ()

    enable_looplift = Option(
        type=bool,
        default=False,
        doc="Enable loop-lifting",
    )
    enable_pyobject = Option(
        type=bool,
        default=False,
        doc="Enable pyobject mode (in general)",
    )
    enable_pyobject_looplift = Option(
        type=bool,
        default=False,
        doc="Enable pyobject mode inside lifted loops",
    )
    enable_ssa = Option(
        type=bool,
        default=True,
        doc="Enable SSA",
    )
    force_pyobject = Option(
        type=bool,
        default=False,
        doc="Force pyobject mode inside the whole function",
    )
    release_gil = Option(
        type=bool,
        default=False,
        doc="Release GIL inside the native function",
    )
    no_compile = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    debuginfo = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    boundscheck = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    forceinline = Option(
        type=bool,
        default=False,
        doc="Force inlining of the function. Overrides _dbg_optnone.",
    )
    no_cpython_wrapper = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    no_cfunc_wrapper = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    auto_parallel = Option(
        type=cpu.ParallelOptions,
        default=cpu.ParallelOptions(False),
        doc="""Enable automatic parallel optimization, can be fine-tuned by
taking a dictionary of sub-options instead of a boolean, see parfor.py for
detail""",
    )
    nrt = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    no_rewrites = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    error_model = Option(
        type=str,
        default="python",
        doc="TODO",
    )
    fastmath = Option(
        type=cpu.FastMathOptions,
        default=cpu.FastMathOptions(False),
        doc="TODO",
    )
    noalias = Option(
        type=bool,
        default=False,
        doc="TODO",
    )
    inline = Option(
        type=cpu.InlineOptions,
        default=cpu.InlineOptions("never"),
        doc="TODO",
    )

    dbg_extend_lifetimes = Option(
        type=bool,
        default=False,
        doc=("Extend variable lifetime for debugging. "
             "This automatically turns on with debug=True."),
    )

    dbg_optnone = Option(
        type=bool,
        default=False,
        doc=("Disable optimization for debug. "
             "Equivalent to adding optnone attribute in the LLVM Function.")
    )

    dbg_directives_only = Option(
        type=bool,
        default=False,
        doc=("Make debug emissions directives-only. "
             "Used when generating lineinfo.")
    )


DEFAULT_FLAGS = Flags()
DEFAULT_FLAGS.nrt = True


CR_FIELDS = ["typing_context",
             "target_context",
             "entry_point",
             "typing_error",
             "type_annotation",
             "signature",
             "objectmode",
             "lifted",
             "fndesc",
             "library",
             "call_helper",
             "environment",
             "metadata",
             # List of functions to call to initialize on unserialization
             # (i.e cache load).
             "reload_init",
             "referenced_envs",
             ]


class CompileResult(namedtuple("_CompileResult", CR_FIELDS)):
    """
    A structure holding results from the compilation of a function.
    """

    __slots__ = ()

    def _reduce(self):
        """
        Reduce a CompileResult to picklable components.
        """
        libdata = self.library.serialize_using_object_code()
        # Make it (un)picklable efficiently
        typeann = str(self.type_annotation)
        fndesc = self.fndesc
        # Those don't need to be pickled and may fail
        fndesc.typemap = fndesc.calltypes = None
        # Include all referenced environments
        referenced_envs = self._find_referenced_environments()
        return (libdata, self.fndesc, self.environment, self.signature,
                self.objectmode, self.lifted, typeann, self.reload_init,
                tuple(referenced_envs))

    def _find_referenced_environments(self):
        """Returns a list of referenced environments
        """
        mod = self.library._final_module
        # Find environments
        referenced_envs = []
        for gv in mod.global_variables:
            gvn = gv.name
            if gvn.startswith("_ZN08NumbaEnv"):
                env = lookup_environment(gvn)
                if env is not None:
                    if env.can_cache():
                        referenced_envs.append(env)
        return referenced_envs

    @classmethod
    def _rebuild(cls, target_context, libdata, fndesc, env,
                 signature, objectmode, lifted, typeann,
                 reload_init, referenced_envs):
        if reload_init:
            # Re-run all
            for fn in reload_init:
                fn()

        library = target_context.codegen().unserialize_library(libdata)
        cfunc = target_context.get_executable(library, fndesc, env)
        cr = cls(target_context=target_context,
                 typing_context=target_context.typing_context,
                 library=library,
                 environment=env,
                 entry_point=cfunc,
                 fndesc=fndesc,
                 type_annotation=typeann,
                 signature=signature,
                 objectmode=objectmode,
                 lifted=lifted,
                 typing_error=None,
                 call_helper=None,
                 metadata=None,  # Do not store, arbitrary & potentially large!
                 reload_init=reload_init,
                 referenced_envs=referenced_envs,
                 )

        # Load Environments
        for env in referenced_envs:
            library.codegen.set_env(env.env_name, env)

        return cr

    @property
    def codegen(self):
        return self.target_context.codegen()

    def dump(self, tab=''):
        print(f'{tab}DUMP {type(self).__name__} {self.entry_point}')
        self.signature.dump(tab=tab + '  ')
        print(f'{tab}END DUMP')


_LowerResult = namedtuple("_LowerResult", [
    "fndesc",
    "call_helper",
    "cfunc",
    "env",
])


def sanitize_compile_result_entries(entries):
    keys = set(entries.keys())
    fieldset = set(CR_FIELDS)
    badnames = keys - fieldset
    if badnames:
        raise NameError(*badnames)
    missing = fieldset - keys
    for k in missing:
        entries[k] = None
    # Avoid keeping alive traceback variables
    err = entries['typing_error']
    if err is not None:
        entries['typing_error'] = err.with_traceback(None)
    return entries


def compile_result(**entries):
    entries = sanitize_compile_result_entries(entries)
    return CompileResult(**entries)


def run_frontend(func, inline_closures=False, emit_dels=False):
    """
    Run the compiler frontend over the given Python function, and return
    the function's canonical Numba IR.

    If inline_closures is Truthy then closure inlining will be run
    If emit_dels is Truthy the ir.Del nodes will be emitted appropriately
    """
    # XXX make this a dedicated Pipeline?
    func_id = bytecode.FunctionIdentity.from_function(func)
    interp = interpreter.Interpreter(func_id)
    bc = bytecode.ByteCode(func_id=func_id)
    func_ir = interp.interpret(bc)
    if inline_closures:
        from numba.core.inline_closurecall import InlineClosureCallPass
        inline_pass = InlineClosureCallPass(func_ir, cpu.ParallelOptions(False),
                                            {}, False)
        inline_pass.run()
    post_proc = postproc.PostProcessor(func_ir)
    post_proc.run(emit_dels)
    return func_ir


class _CompileStatus(object):
    """
    Describes the state of compilation. Used like a C record.
    """
    __slots__ = ['fail_reason', 'can_fallback']

    def __init__(self, can_fallback):
        self.fail_reason = None
        self.can_fallback = can_fallback

    def __repr__(self):
        vals = []
        for k in self.__slots__:
            vals.append("{k}={v}".format(k=k, v=getattr(self, k)))
        return ', '.join(vals)


class _EarlyPipelineCompletion(Exception):
    """
    Raised to indicate that a pipeline has completed early
    """

    def __init__(self, result):
        self.result = result


class StateDict(dict):
    """
    A dictionary that has an overloaded getattr and setattr to permit getting
    and setting key/values through the use of attributes.
    """

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        self[attr] = value


def _make_subtarget(targetctx, flags):
    """
    Make a new target context from the given target context and flags.
    """
    subtargetoptions = {}
    if flags.debuginfo:
        subtargetoptions['enable_debuginfo'] = True
    if flags.boundscheck:
        subtargetoptions['enable_boundscheck'] = True
    if flags.nrt:
        subtargetoptions['enable_nrt'] = True
    if flags.auto_parallel:
        subtargetoptions['auto_parallel'] = flags.auto_parallel
    if flags.fastmath:
        subtargetoptions['fastmath'] = flags.fastmath
    error_model = callconv.create_error_model(flags.error_model, targetctx)
    subtargetoptions['error_model'] = error_model

    return targetctx.subtarget(**subtargetoptions)


class CompilerBase(object):
    """
    Stores and manages states for the compiler
    """

    def __init__(self, typingctx, targetctx, library, args, return_type, flags,
                 locals):
        # Make sure the environment is reloaded
        config.reload_config()
        typingctx.refresh()
        targetctx.refresh()

        self.state = StateDict()

        self.state.typingctx = typingctx
        self.state.targetctx = _make_subtarget(targetctx, flags)
        self.state.library = library
        self.state.args = args
        self.state.return_type = return_type
        self.state.flags = flags
        self.state.locals = locals

        # Results of various steps of the compilation pipeline
        self.state.bc = None
        self.state.func_id = None
        self.state.func_ir = None
        self.state.lifted = None
        self.state.lifted_from = None
        self.state.typemap = None
        self.state.calltypes = None
        self.state.type_annotation = None
        # holds arbitrary inter-pipeline stage meta data
        self.state.metadata = {}
        self.state.reload_init = []
        # hold this for e.g. with_lifting, null out on exit
        self.state.pipeline = self

        # parfor diagnostics info, add to metadata
        self.state.parfor_diagnostics = ParforDiagnostics()
        self.state.metadata['parfor_diagnostics'] = \
            self.state.parfor_diagnostics
        self.state.metadata['parfors'] = {}

        self.state.status = _CompileStatus(
            can_fallback=self.state.flags.enable_pyobject
        )

    def compile_extra(self, func):
        self.state.func_id = bytecode.FunctionIdentity.from_function(func)
        ExtractByteCode().run_pass(self.state)

        self.state.lifted = ()
        self.state.lifted_from = None
        return self._compile_bytecode()

    def compile_ir(self, func_ir, lifted=(), lifted_from=None):
        self.state.func_id = func_ir.func_id
        self.state.lifted = lifted
        self.state.lifted_from = lifted_from
        self.state.func_ir = func_ir
        self.state.nargs = self.state.func_ir.arg_count

        FixupArgs().run_pass(self.state)
        return self._compile_ir()

    def define_pipelines(self):
        """Child classes override this to customize the pipelines in use.
        """
        raise NotImplementedError()

    def _compile_core(self):
        """
        Populate and run compiler pipeline
        """
        with ConfigStack().enter(self.state.flags.copy()):
            pms = self.define_pipelines()
            for pm in pms:
                pipeline_name = pm.pipeline_name
                func_name = "%s.%s" % (self.state.func_id.modname,
                                       self.state.func_id.func_qualname)

                event("Pipeline: %s for %s" % (pipeline_name, func_name))
                self.state.metadata['pipeline_times'] = {pipeline_name:
                                                         pm.exec_times}
                is_final_pipeline = pm == pms[-1]
                res = None
                try:
                    pm.run(self.state)
                    if self.state.cr is not None:
                        break
                except _EarlyPipelineCompletion as e:
                    res = e.result
                    break
                except Exception as e:
                    if not isinstance(e, errors.NumbaError):
                        raise e
                    self.state.status.fail_reason = e
                    if is_final_pipeline:
                        raise e
            else:
                raise CompilerError("All available pipelines exhausted")

            # Pipeline is done, remove self reference to release refs to user
            # code
            self.state.pipeline = None

            # organise a return
            if res is not None:
                # Early pipeline completion
                return res
            else:
                assert self.state.cr is not None
                return self.state.cr

    def _compile_bytecode(self):
        """
        Populate and run pipeline for bytecode input
        """
        assert self.state.func_ir is None
        return self._compile_core()

    def _compile_ir(self):
        """
        Populate and run pipeline for IR input
        """
        assert self.state.func_ir is not None
        return self._compile_core()


class Compiler(CompilerBase):
    """The default compiler
    """

    def define_pipelines(self):
        if self.state.flags.force_pyobject:
            # either object mode
            return [DefaultPassBuilder.define_objectmode_pipeline(self.state),]
        else:
            # or nopython mode
            return [DefaultPassBuilder.define_nopython_pipeline(self.state),]


class DefaultPassBuilder(object):
    """
    This is the default pass builder, it contains the "classic" default
    pipelines as pre-canned PassManager instances:
      - nopython
      - objectmode
      - interpreted
      - typed
      - untyped
      - nopython lowering
    """
    @staticmethod
    def define_nopython_pipeline(state, name='nopython'):
        """Returns an nopython mode pipeline based PassManager
        """
        # compose pipeline from untyped, typed and lowering parts
        dpb = DefaultPassBuilder
        pm = PassManager(name)
        untyped_passes = dpb.define_untyped_pipeline(state)
        pm.passes.extend(untyped_passes.passes)

        typed_passes = dpb.define_typed_pipeline(state)
        pm.passes.extend(typed_passes.passes)

        lowering_passes = dpb.define_nopython_lowering_pipeline(state)
        pm.passes.extend(lowering_passes.passes)

        pm.finalize()
        return pm

    @staticmethod
    def define_nopython_lowering_pipeline(state, name='nopython_lowering'):
        pm = PassManager(name)
        # legalise
        pm.add_pass(NoPythonSupportedFeatureValidation,
                    "ensure features that are in use are in a valid form")
        pm.add_pass(IRLegalization,
                    "ensure IR is legal prior to lowering")
        # Annotate only once legalized
        pm.add_pass(AnnotateTypes, "annotate types")
        # lower
        if state.flags.auto_parallel.enabled:
            pm.add_pass(NativeParforLowering, "native parfor lowering")
        else:
            pm.add_pass(NativeLowering, "native lowering")
        pm.add_pass(NoPythonBackend, "nopython mode backend")
        pm.add_pass(DumpParforDiagnostics, "dump parfor diagnostics")
        pm.finalize()
        return pm

    @staticmethod
    def define_parfor_gufunc_nopython_lowering_pipeline(
            state, name='parfor_gufunc_nopython_lowering'):
        pm = PassManager(name)
        # legalise
        pm.add_pass(NoPythonSupportedFeatureValidation,
                    "ensure features that are in use are in a valid form")
        pm.add_pass(IRLegalization,
                    "ensure IR is legal prior to lowering")
        # Annotate only once legalized
        pm.add_pass(AnnotateTypes, "annotate types")
        # lower
        if state.flags.auto_parallel.enabled:
            pm.add_pass(NativeParforLowering, "native parfor lowering")
        else:
            pm.add_pass(NativeLowering, "native lowering")
        pm.add_pass(NoPythonBackend, "nopython mode backend")
        pm.finalize()
        return pm

    @staticmethod
    def define_typed_pipeline(state, name="typed"):
        """Returns the typed part of the nopython pipeline"""
        pm = PassManager(name)
        # typing
        pm.add_pass(NopythonTypeInference, "nopython frontend")

        # strip phis
        pm.add_pass(PreLowerStripPhis, "remove phis nodes")

        # optimisation
        pm.add_pass(InlineOverloads, "inline overloaded functions")
        if state.flags.auto_parallel.enabled:
            pm.add_pass(PreParforPass, "Preprocessing for parfors")
        if not state.flags.no_rewrites:
            pm.add_pass(NopythonRewrites, "nopython rewrites")
        if state.flags.auto_parallel.enabled:
            pm.add_pass(ParforPass, "convert to parfors")
            pm.add_pass(ParforFusionPass, "fuse parfors")
            pm.add_pass(ParforPreLoweringPass, "parfor prelowering")

        pm.finalize()
        return pm

    @staticmethod
    def define_parfor_gufunc_pipeline(state, name="parfor_gufunc_typed"):
        """Returns the typed part of the nopython pipeline"""
        pm = PassManager(name)
        assert state.func_ir
        pm.add_pass(IRProcessing, "processing IR")
        pm.add_pass(NopythonTypeInference, "nopython frontend")
        pm.add_pass(ParforPreLoweringPass, "parfor prelowering")

        pm.finalize()
        return pm

    @staticmethod
    def define_untyped_pipeline(state, name='untyped'):
        """Returns an untyped part of the nopython pipeline"""
        pm = PassManager(name)
        if state.func_ir is None:
            pm.add_pass(TranslateByteCode, "analyzing bytecode")
            pm.add_pass(FixupArgs, "fix up args")
        pm.add_pass(IRProcessing, "processing IR")
        pm.add_pass(WithLifting, "Handle with contexts")

        # inline closures early in case they are using nonlocal's
        # see issue #6585.
        pm.add_pass(InlineClosureLikes,
                    "inline calls to locally defined closures")

        # pre typing
        if not state.flags.no_rewrites:
            pm.add_pass(RewriteSemanticConstants, "rewrite semantic constants")
            pm.add_pass(DeadBranchPrune, "dead branch pruning")
            pm.add_pass(GenericRewrites, "nopython rewrites")

        pm.add_pass(RewriteDynamicRaises, "rewrite dynamic raises")

        # convert any remaining closures into functions
        pm.add_pass(MakeFunctionToJitFunction,
                    "convert make_function into JIT functions")
        # inline functions that have been determined as inlinable and rerun
        # branch pruning, this needs to be run after closures are inlined as
        # the IR repr of a closure masks call sites if an inlinable is called
        # inside a closure
        pm.add_pass(InlineInlinables, "inline inlinable functions")
        if not state.flags.no_rewrites:
            pm.add_pass(DeadBranchPrune, "dead branch pruning")

        pm.add_pass(FindLiterallyCalls, "find literally calls")
        pm.add_pass(LiteralUnroll, "handles literal_unroll")

        if state.flags.enable_ssa:
            pm.add_pass(ReconstructSSA, "ssa")

        if not state.flags.no_rewrites:
            pm.add_pass(DeadBranchPrune, "dead branch pruning")

        pm.add_pass(LiteralPropagationSubPipelinePass, "Literal propagation")

        pm.finalize()
        return pm

    @staticmethod
    def define_objectmode_pipeline(state, name='object'):
        """Returns an object-mode pipeline based PassManager
        """
        pm = PassManager(name)
        if state.func_ir is None:
            pm.add_pass(TranslateByteCode, "analyzing bytecode")
            pm.add_pass(FixupArgs, "fix up args")
        else:
            # Reaches here if it's a fallback from nopython mode.
            # Strip the phi nodes.
            pm.add_pass(PreLowerStripPhis, "remove phis nodes")
        pm.add_pass(IRProcessing, "processing IR")

        # The following passes are needed to adjust for looplifting
        pm.add_pass(CanonicalizeLoopEntry, "canonicalize loop entry")
        pm.add_pass(CanonicalizeLoopExit, "canonicalize loop exit")

        pm.add_pass(ObjectModeFrontEnd, "object mode frontend")
        pm.add_pass(InlineClosureLikes,
                    "inline calls to locally defined closures")
        # convert any remaining closures into functions
        pm.add_pass(MakeFunctionToJitFunction,
                    "convert make_function into JIT functions")
        pm.add_pass(IRLegalization, "ensure IR is legal prior to lowering")
        pm.add_pass(AnnotateTypes, "annotate types")
        pm.add_pass(ObjectModeBackEnd, "object mode backend")
        pm.finalize()
        return pm


def compile_extra(typingctx, targetctx, func, args, return_type, flags,
                  locals, library=None, pipeline_class=Compiler):
    """Compiler entry point

    Parameter
    ---------
    typingctx :
        typing context
    targetctx :
        target context
    func : function
        the python function to be compiled
    args : tuple, list
        argument types
    return_type :
        Use ``None`` to indicate void return
    flags : numba.compiler.Flags
        compiler flags
    library : numba.codegen.CodeLibrary
        Used to store the compiled code.
        If it is ``None``, a new CodeLibrary is used.
    pipeline_class : type like numba.compiler.CompilerBase
        compiler pipeline
    """
    pipeline = pipeline_class(typingctx, targetctx, library,
                              args, return_type, flags, locals)
    return pipeline.compile_extra(func)


def compile_ir(typingctx, targetctx, func_ir, args, return_type, flags,
               locals, lifted=(), lifted_from=None, is_lifted_loop=False,
               library=None, pipeline_class=Compiler):
    """
    Compile a function with the given IR.

    For internal use only.
    """

    # This is a special branch that should only run on IR from a lifted loop
    if is_lifted_loop:
        # This code is pessimistic and costly, but it is a not often trodden
        # path and it will go away once IR is made immutable. The problem is
        # that the rewrite passes can mutate the IR into a state that makes
        # it possible for invalid tokens to be transmitted to lowering which
        # then trickle through into LLVM IR and causes RuntimeErrors as LLVM
        # cannot compile it. As a result the following approach is taken:
        # 1. Create some new flags that copy the original ones but switch
        #    off rewrites.
        # 2. Compile with 1. to get a compile result
        # 3. Try and compile another compile result but this time with the
        #    original flags (and IR being rewritten).
        # 4. If 3 was successful, use the result, else use 2.

        # create flags with no rewrites
        norw_flags = copy.deepcopy(flags)
        norw_flags.no_rewrites = True

        def compile_local(the_ir, the_flags):
            pipeline = pipeline_class(typingctx, targetctx, library,
                                      args, return_type, the_flags, locals)
            return pipeline.compile_ir(func_ir=the_ir, lifted=lifted,
                                       lifted_from=lifted_from)

        # compile with rewrites off, IR shouldn't be mutated irreparably
        norw_cres = compile_local(func_ir.copy(), norw_flags)

        # try and compile with rewrites on if no_rewrites was not set in the
        # original flags, IR might get broken but we've got a CompileResult
        # that's usable from above.
        rw_cres = None
        if not flags.no_rewrites:
            # Suppress warnings in compilation retry
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", errors.NumbaWarning)
                try:
                    rw_cres = compile_local(func_ir.copy(), flags)
                except Exception:
                    pass
        # if the rewrite variant of compilation worked, use it, else use
        # the norewrites backup
        if rw_cres is not None:
            cres = rw_cres
        else:
            cres = norw_cres
        return cres

    else:
        pipeline = pipeline_class(typingctx, targetctx, library,
                                  args, return_type, flags, locals)
        return pipeline.compile_ir(func_ir=func_ir, lifted=lifted,
                                   lifted_from=lifted_from)


def compile_internal(typingctx, targetctx, library,
                     func, args, return_type, flags, locals):
    """
    For internal use only.
    """
    pipeline = Compiler(typingctx, targetctx, library,
                        args, return_type, flags, locals)
    return pipeline.compile_extra(func)
