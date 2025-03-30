from ctypes import (c_bool, c_char_p, c_int, c_size_t, Structure, byref,
                    POINTER)
from collections import namedtuple
from enum import IntFlag
from llvmlite.binding import ffi
from llvmlite.binding.initfini import llvm_version_info
import os
from tempfile import mkstemp
from llvmlite.binding.common import _encode_string

llvm_version_major = llvm_version_info[0]

_prunestats = namedtuple('PruneStats',
                         ('basicblock diamond fanout fanout_raise'))


class PruneStats(_prunestats):
    """ Holds statistics from reference count pruning.
    """

    def __add__(self, other):
        if not isinstance(other, PruneStats):
            msg = 'PruneStats can only be added to another PruneStats, got {}.'
            raise TypeError(msg.format(type(other)))
        return PruneStats(self.basicblock + other.basicblock,
                          self.diamond + other.diamond,
                          self.fanout + other.fanout,
                          self.fanout_raise + other.fanout_raise)

    def __sub__(self, other):
        if not isinstance(other, PruneStats):
            msg = ('PruneStats can only be subtracted from another PruneStats, '
                   'got {}.')
            raise TypeError(msg.format(type(other)))
        return PruneStats(self.basicblock - other.basicblock,
                          self.diamond - other.diamond,
                          self.fanout - other.fanout,
                          self.fanout_raise - other.fanout_raise)


class _c_PruneStats(Structure):
    _fields_ = [
        ('basicblock', c_size_t),
        ('diamond', c_size_t),
        ('fanout', c_size_t),
        ('fanout_raise', c_size_t)]


def dump_refprune_stats(printout=False):
    """ Returns a namedtuple containing the current values for the refop pruning
    statistics. If kwarg `printout` is True the stats are printed to stderr,
    default is False.
    """

    stats = _c_PruneStats(0, 0, 0, 0)
    do_print = c_bool(printout)

    ffi.lib.LLVMPY_DumpRefPruneStats(byref(stats), do_print)
    return PruneStats(stats.basicblock, stats.diamond, stats.fanout,
                      stats.fanout_raise)


def set_time_passes(enable):
    """Enable or disable the pass timers.

    Parameters
    ----------
    enable : bool
        Set to True to enable the pass timers.
        Set to False to disable the pass timers.
    """
    ffi.lib.LLVMPY_SetTimePasses(c_bool(enable))


def report_and_reset_timings():
    """Returns the pass timings report and resets the LLVM internal timers.

    Pass timers are enabled by ``set_time_passes()``. If the timers are not
    enabled, this function will return an empty string.

    Returns
    -------
    res : str
        LLVM generated timing report.
    """
    with ffi.OutputString() as buf:
        ffi.lib.LLVMPY_ReportAndResetTimings(buf)
        return str(buf)


def create_module_pass_manager():
    return ModulePassManager()


def create_function_pass_manager(module):
    return FunctionPassManager(module)


class RefPruneSubpasses(IntFlag):
    PER_BB       = 0b0001    # noqa: E221
    DIAMOND      = 0b0010    # noqa: E221
    FANOUT       = 0b0100    # noqa: E221
    FANOUT_RAISE = 0b1000
    ALL = PER_BB | DIAMOND | FANOUT | FANOUT_RAISE


class PassManager(ffi.ObjectRef):
    """PassManager
    """

    def _dispose(self):
        self._capi.LLVMPY_DisposePassManager(self)

    def add_aa_eval_pass(self):
        """
        See https://llvm.org/docs/Passes.html#aa-eval-exhaustive-alias-analysis-precision-evaluator

        LLVM 14: `llvm::createAAEvalPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddAAEvalPass(self)

    def add_basic_aa_pass(self):
        """
        See https://llvm.org/docs/Passes.html#basic-aa-basic-alias-analysis-stateless-aa-impl

        LLVM 14: `llvm::createBasicAAWrapperPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddBasicAAWrapperPass(self)

    def add_constant_merge_pass(self):
        """
        See http://llvm.org/docs/Passes.html#constmerge-merge-duplicate-global-constants

        LLVM 14: `LLVMAddConstantMergePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddConstantMergePass(self)

    def add_dead_arg_elimination_pass(self):
        """
        See http://llvm.org/docs/Passes.html#deadargelim-dead-argument-elimination

        LLVM 14: `LLVMAddDeadArgEliminationPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDeadArgEliminationPass(self)

    def add_dependence_analysis_pass(self):
        """
        See https://llvm.org/docs/Passes.html#da-dependence-analysis

        LLVM 14: `llvm::createDependenceAnalysisWrapperPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDependenceAnalysisPass(self)

    def add_dot_call_graph_pass(self):
        """
        See https://llvm.org/docs/Passes.html#dot-callgraph-print-call-graph-to-dot-file

        LLVM 14: `llvm::createCallGraphDOTPrinterPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddCallGraphDOTPrinterPass(self)

    def add_dot_cfg_printer_pass(self):
        """
        See https://llvm.org/docs/Passes.html#dot-cfg-print-cfg-of-function-to-dot-file

        LLVM 14: `llvm::createCFGPrinterLegacyPassPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddCFGPrinterPass(self)

    def add_dot_dom_printer_pass(self, show_body=False):
        """
        See https://llvm.org/docs/Passes.html#dot-dom-print-dominance-tree-of-function-to-dot-file

        LLVM 14: `llvm::createDomPrinterPass` and `llvm::createDomOnlyPrinterPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDotDomPrinterPass(self, show_body)

    def add_dot_postdom_printer_pass(self, show_body=False):
        """
        See https://llvm.org/docs/Passes.html#dot-postdom-print-postdominance-tree-of-function-to-dot-file

        LLVM 14: `llvm::createPostDomPrinterPass` and `llvm::createPostDomOnlyPrinterPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDotPostDomPrinterPass(self, show_body)

    def add_globals_mod_ref_aa_pass(self):
        """
        See https://llvm.org/docs/Passes.html#globalsmodref-aa-simple-mod-ref-analysis-for-globals

        LLVM 14: `llvm::createGlobalsAAWrapperPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddGlobalsModRefAAPass(self)

    def add_iv_users_pass(self):
        """
        See https://llvm.org/docs/Passes.html#iv-users-induction-variable-users

        LLVM 14: `llvm::createIVUsersPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddIVUsersPass(self)

    def add_lint_pass(self):
        """
        See https://llvm.org/docs/Passes.html#lint-statically-lint-checks-llvm-ir

        LLVM 14: `llvm::createLintLegacyPassPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLintPass(self)

    def add_lazy_value_info_pass(self):
        """
        See https://llvm.org/docs/Passes.html#lazy-value-info-lazy-value-information-analysis

        LLVM 14: `llvm::createLazyValueInfoPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLazyValueInfoPass(self)

    def add_module_debug_info_pass(self):
        """
        See https://llvm.org/docs/Passes.html#module-debuginfo-decodes-module-level-debug-info

        LLVM 14: `llvm::createModuleDebugInfoPrinterPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddModuleDebugInfoPrinterPass(self)

    def add_region_info_pass(self):
        """
        See https://llvm.org/docs/Passes.html#regions-detect-single-entry-single-exit-regions

        LLVM 14: `llvm::createRegionInfoPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddRegionInfoPass(self)

    def add_scalar_evolution_aa_pass(self):
        """
        See https://llvm.org/docs/Passes.html#scev-aa-scalarevolution-based-alias-analysis

        LLVM 14: `llvm::createSCEVAAWrapperPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddScalarEvolutionAAPass(self)

    def add_aggressive_dead_code_elimination_pass(self):
        """
        See https://llvm.org/docs/Passes.html#adce-aggressive-dead-code-elimination

        LLVM 14: `llvm::createAggressiveDCEPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddAggressiveDCEPass(self)

    def add_always_inliner_pass(self, insert_lifetime=True):
        """
        See https://llvm.org/docs/Passes.html#always-inline-inliner-for-always-inline-functions

        LLVM 14: `llvm::createAlwaysInlinerLegacyPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddAlwaysInlinerPass(self, insert_lifetime)

    def add_arg_promotion_pass(self, max_elements=3):
        """
        See https://llvm.org/docs/Passes.html#argpromotion-promote-by-reference-arguments-to-scalars

        LLVM 14: `llvm::createArgumentPromotionPass`
        """  # noqa E501
        raise RuntimeError('ArgumentPromotionPass unavailable in LLVM > 14')

    def add_break_critical_edges_pass(self):
        """
        See https://llvm.org/docs/Passes.html#break-crit-edges-break-critical-edges-in-cfg

        LLVM 14: `llvm::createBreakCriticalEdgesPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddBreakCriticalEdgesPass(self)

    def add_dead_store_elimination_pass(self):
        """
        See https://llvm.org/docs/Passes.html#dse-dead-store-elimination

        LLVM 14: `llvm::createDeadStoreEliminationPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDeadStoreEliminationPass(self)

    def add_reverse_post_order_function_attrs_pass(self):
        """
        See https://llvm.org/docs/Passes.html#function-attrs-deduce-function-attributes

        LLVM 14: `llvm::createReversePostOrderFunctionAttrsPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddReversePostOrderFunctionAttrsPass(self)

    def add_function_attrs_pass(self):
        """
        See http://llvm.org/docs/Passes.html#functionattrs-deduce-function-attributes

        LLVM 14: `LLVMAddFunctionAttrsPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddFunctionAttrsPass(self)

    def add_function_inlining_pass(self, threshold):
        """
        See http://llvm.org/docs/Passes.html#inline-function-integration-inlining

        LLVM 14: `createFunctionInliningPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddFunctionInliningPass(self, threshold)

    def add_global_dce_pass(self):
        """
        See http://llvm.org/docs/Passes.html#globaldce-dead-global-elimination

        LLVM 14: `LLVMAddGlobalDCEPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddGlobalDCEPass(self)

    def add_global_optimizer_pass(self):
        """
        See http://llvm.org/docs/Passes.html#globalopt-global-variable-optimizer

        LLVM 14: `LLVMAddGlobalOptimizerPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddGlobalOptimizerPass(self)

    def add_ipsccp_pass(self):
        """
        See http://llvm.org/docs/Passes.html#ipsccp-interprocedural-sparse-conditional-constant-propagation

        LLVM 14: `LLVMAddIPSCCPPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddIPSCCPPass(self)

    def add_dead_code_elimination_pass(self):
        """
        See http://llvm.org/docs/Passes.html#dce-dead-code-elimination
        LLVM 14: `llvm::createDeadCodeEliminationPass`
        """
        ffi.lib.LLVMPY_AddDeadCodeEliminationPass(self)

    def add_aggressive_instruction_combining_pass(self):
        """
        See https://llvm.org/docs/Passes.html#aggressive-instcombine-combine-expression-patterns

        LLVM 14: `llvm::createAggressiveInstCombinerPass`
        """  # noqa E501
        if llvm_version_major > 15:
            msg = "AggressiveInstrCombinerPass unavailable in LLVM > 15"
            raise RuntimeError(msg)

        ffi.lib.LLVMPY_AddAggressiveInstructionCombiningPass(self)

    def add_internalize_pass(self):
        """
        See https://llvm.org/docs/Passes.html#internalize-internalize-global-symbols

        LLVM 14: `llvm::createInternalizePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddInternalizePass(self)

    def add_cfg_simplification_pass(self):
        """
        See http://llvm.org/docs/Passes.html#simplifycfg-simplify-the-cfg

        LLVM 14: `LLVMAddCFGSimplificationPass`
        """
        ffi.lib.LLVMPY_AddCFGSimplificationPass(self)

    def add_jump_threading_pass(self, threshold=-1):
        """
        See https://llvm.org/docs/Passes.html#jump-threading-jump-threading

        LLVM 14: `llvm::createJumpThreadingPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddJumpThreadingPass(self, threshold)

    def add_lcssa_pass(self):
        """
        See https://llvm.org/docs/Passes.html#lcssa-loop-closed-ssa-form-pass

        LLVM 14: `llvm::createLCSSAPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLCSSAPass(self)

    def add_gvn_pass(self):
        """
        See http://llvm.org/docs/Passes.html#gvn-global-value-numbering

        LLVM 14: `LLVMAddGVNPass`
        """
        ffi.lib.LLVMPY_AddGVNPass(self)

    def add_instruction_combining_pass(self):
        """
        See http://llvm.org/docs/Passes.html#passes-instcombine

        LLVM 14: `LLVMAddInstructionCombiningPass`
        """
        ffi.lib.LLVMPY_AddInstructionCombiningPass(self)

    def add_licm_pass(self):
        """
        See http://llvm.org/docs/Passes.html#licm-loop-invariant-code-motion

        LLVM 14: `LLVMAddLICMPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLICMPass(self)

    def add_loop_deletion_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-deletion-delete-dead-loops

        LLVM 14: `llvm::createLoopDeletionPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopDeletionPass(self)

    def add_loop_extractor_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-extract-extract-loops-into-new-functions

        LLVM 14: `llvm::createLoopExtractorPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopExtractorPass(self)

    def add_single_loop_extractor_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-extract-single-extract-at-most-one-loop-into-a-new-function

        LLVM 14: `llvm::createSingleLoopExtractorPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddSingleLoopExtractorPass(self)

    def add_sccp_pass(self):
        """
        See http://llvm.org/docs/Passes.html#sccp-sparse-conditional-constant-propagation

        LLVM 14: `LLVMAddSCCPPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddSCCPPass(self)

    def add_loop_strength_reduce_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-reduce-loop-strength-reduction

        LLVM 14: `llvm::createLoopStrengthReducePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopStrengthReducePass(self)

    def add_loop_simplification_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-simplify-canonicalize-natural-loops

        LLVM 14: `llvm::createLoopSimplifyPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopSimplificationPass(self)

    def add_loop_unroll_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-unroll-unroll-loops

        LLVM 14: `LLVMAddLoopUnrollPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopUnrollPass(self)

    def add_loop_unroll_and_jam_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loop-unroll-and-jam-unroll-and-jam-loops

        LLVM 14: `LLVMAddLoopUnrollAndJamPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopUnrollAndJamPass(self)

    def add_loop_unswitch_pass(self,
                               optimize_for_size=False,
                               has_branch_divergence=False):
        """
        See https://llvm.org/docs/Passes.html#loop-unswitch-unswitch-loops

        LLVM 14: `llvm::createLoopUnswitchPass`
        LLVM 15: `llvm::createSimpleLoopUnswitchLegacyPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLoopUnswitchPass(self, optimize_for_size,
                                           has_branch_divergence)

    def add_lower_atomic_pass(self):
        """
        See https://llvm.org/docs/Passes.html#loweratomic-lower-atomic-intrinsics-to-non-atomic-form

        LLVM 14: `llvm::createLowerAtomicPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLowerAtomicPass(self)

    def add_lower_invoke_pass(self):
        """
        See https://llvm.org/docs/Passes.html#lowerinvoke-lower-invokes-to-calls-for-unwindless-code-generators

        LLVM 14: `llvm::createLowerInvokePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLowerInvokePass(self)

    def add_lower_switch_pass(self):
        """
        See https://llvm.org/docs/Passes.html#lowerswitch-lower-switchinsts-to-branches

        LLVM 14: `llvm::createLowerSwitchPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddLowerSwitchPass(self)

    def add_memcpy_optimization_pass(self):
        """
        See https://llvm.org/docs/Passes.html#memcpyopt-memcpy-optimization

        LLVM 14: `llvm::createMemCpyOptPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddMemCpyOptimizationPass(self)

    def add_merge_functions_pass(self):
        """
        See https://llvm.org/docs/Passes.html#mergefunc-merge-functions

        LLVM 14: `llvm::createMergeFunctionsPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddMergeFunctionsPass(self)

    def add_merge_returns_pass(self):
        """
        See https://llvm.org/docs/Passes.html#mergereturn-unify-function-exit-nodes

        LLVM 14: `llvm::createUnifyFunctionExitNodesPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddMergeReturnsPass(self)

    def add_partial_inlining_pass(self):
        """
        See https://llvm.org/docs/Passes.html#partial-inliner-partial-inliner

        LLVM 14: `llvm::createPartialInliningPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddPartialInliningPass(self)

    def add_prune_exception_handling_pass(self):
        """
        See https://llvm.org/docs/Passes.html#prune-eh-remove-unused-exception-handling-info

        LLVM 14: `llvm::createPruneEHPass`
        """  # noqa E501
        if llvm_version_major > 15:
            raise RuntimeError("PruneEHPass unavailable in LLVM > 15")
        ffi.lib.LLVMPY_AddPruneExceptionHandlingPass(self)

    def add_reassociate_expressions_pass(self):
        """
        See https://llvm.org/docs/Passes.html#reassociate-reassociate-expressions

        LLVM 14: `llvm::createReassociatePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddReassociatePass(self)

    def add_demote_register_to_memory_pass(self):
        """
        See https://llvm.org/docs/Passes.html#rel-lookup-table-converter-relative-lookup-table-converter

        LLVM 14: `llvm::createDemoteRegisterToMemoryPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddDemoteRegisterToMemoryPass(self)

    def add_sroa_pass(self):
        """
        See http://llvm.org/docs/Passes.html#scalarrepl-scalar-replacement-of-aggregates-dt
        Note that this pass corresponds to the ``opt -sroa`` command-line option,
        despite the link above.

        LLVM 14: `llvm::createSROAPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddSROAPass(self)

    def add_sink_pass(self):
        """
        See https://llvm.org/docs/Passes.html#sink-code-sinking

        LLVM 14: `llvm::createSinkingPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddSinkPass(self)

    def add_strip_symbols_pass(self, only_debug=False):
        """
        See https://llvm.org/docs/Passes.html#strip-strip-all-symbols-from-a-module

        LLVM 14: `llvm::createStripSymbolsPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddStripSymbolsPass(self, only_debug)

    def add_strip_dead_debug_info_pass(self):
        """
        See https://llvm.org/docs/Passes.html#strip-dead-debug-info-strip-debug-info-for-unused-symbols

        LLVM 14: `llvm::createStripDeadDebugInfoPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddStripDeadDebugInfoPass(self)

    def add_strip_dead_prototypes_pass(self):
        """
        See https://llvm.org/docs/Passes.html#strip-dead-prototypes-strip-unused-function-prototypes

        LLVM 14: `llvm::createStripDeadPrototypesPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddStripDeadPrototypesPass(self)

    def add_strip_debug_declare_pass(self):
        """
        See https://llvm.org/docs/Passes.html#strip-debug-declare-strip-all-llvm-dbg-declare-intrinsics

        LLVM 14: `llvm::createStripDebugDeclarePass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddStripDebugDeclarePrototypesPass(self)

    def add_strip_nondebug_symbols_pass(self):
        """
        See https://llvm.org/docs/Passes.html#strip-nondebug-strip-all-symbols-except-dbg-symbols-from-a-module

        LLVM 14: `llvm::createStripNonDebugSymbolsPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddStripNondebugSymbolsPass(self)

    def add_tail_call_elimination_pass(self):
        """
        See https://llvm.org/docs/Passes.html#tailcallelim-tail-call-elimination

        LLVM 14: `llvm::createTailCallEliminationPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddTailCallEliminationPass(self)

    def add_type_based_alias_analysis_pass(self):
        """
        LLVM 14: `LLVMAddTypeBasedAliasAnalysisPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddTypeBasedAliasAnalysisPass(self)

    def add_basic_alias_analysis_pass(self):
        """
        See http://llvm.org/docs/AliasAnalysis.html#the-basicaa-pass

        LLVM 14: `LLVMAddBasicAliasAnalysisPass`
        """
        ffi.lib.LLVMPY_AddBasicAliasAnalysisPass(self)

    def add_loop_rotate_pass(self):
        """http://llvm.org/docs/Passes.html#loop-rotate-rotate-loops."""
        ffi.lib.LLVMPY_AddLoopRotatePass(self)

    def add_target_library_info(self, triple):
        ffi.lib.LLVMPY_AddTargetLibraryInfoPass(self, _encode_string(triple))

    def add_instruction_namer_pass(self):
        """
        See https://llvm.org/docs/Passes.html#instnamer-assign-names-to-anonymous-instructions.

        LLVM 14: `llvm::createInstructionNamerPass`
        """  # noqa E501
        ffi.lib.LLVMPY_AddInstructionNamerPass(self)

    # Non-standard LLVM passes

    def add_refprune_pass(self, subpasses_flags=RefPruneSubpasses.ALL,
                          subgraph_limit=1000):
        """Add Numba specific Reference count pruning pass.

        Parameters
        ----------
        subpasses_flags : RefPruneSubpasses
            A bitmask to control the subpasses to be enabled.
        subgraph_limit : int
            Limit the fanout pruners to working on a subgraph no bigger than
            this number of basic-blocks to avoid spending too much time in very
            large graphs. Default is 1000. Subject to change in future
            versions.
        """
        iflags = RefPruneSubpasses(subpasses_flags)
        ffi.lib.LLVMPY_AddLegacyRefPrunePass(self, iflags, subgraph_limit)


class ModulePassManager(PassManager):

    def __init__(self, ptr=None):
        if ptr is None:
            ptr = ffi.lib.LLVMPY_CreatePassManager()
        PassManager.__init__(self, ptr)

    def run(self, module, remarks_file=None, remarks_format='yaml',
            remarks_filter=''):
        """
        Run optimization passes on the given module.

        Parameters
        ----------
        module : llvmlite.binding.ModuleRef
            The module to be optimized inplace
        remarks_file : str; optional
            If not `None`, it is the file to store the optimization remarks.
        remarks_format : str; optional
            The format to write; YAML is default
        remarks_filter : str; optional
            The filter that should be applied to the remarks output.
        """
        if remarks_file is None:
            return ffi.lib.LLVMPY_RunPassManager(self, module)
        else:
            r = ffi.lib.LLVMPY_RunPassManagerWithRemarks(
                self, module, _encode_string(remarks_format),
                _encode_string(remarks_filter),
                _encode_string(remarks_file))
            if r == -1:
                raise IOError("Failed to initialize remarks file.")
            return r > 0

    def run_with_remarks(self, module, remarks_format='yaml',
                         remarks_filter=''):
        """
        Run optimization passes on the given module and returns the result and
        the remarks data.

        Parameters
        ----------
        module : llvmlite.binding.ModuleRef
            The module to be optimized
        remarks_format : str
            The remarks output; YAML is the default
        remarks_filter : str; optional
            The filter that should be applied to the remarks output.
        """
        remarkdesc, remarkfile = mkstemp()
        try:
            with os.fdopen(remarkdesc, 'r'):
                pass
            r = self.run(module, remarkfile, remarks_format, remarks_filter)
            if r == -1:
                raise IOError("Failed to initialize remarks file.")
            with open(remarkfile) as f:
                return bool(r), f.read()
        finally:
            os.unlink(remarkfile)


class FunctionPassManager(PassManager):

    def __init__(self, module):
        ptr = ffi.lib.LLVMPY_CreateFunctionPassManager(module)
        self._module = module
        module._owned = True
        PassManager.__init__(self, ptr)

    def initialize(self):
        """
        Initialize the FunctionPassManager.  Returns True if it produced
        any changes (?).
        """
        return ffi.lib.LLVMPY_InitializeFunctionPassManager(self)

    def finalize(self):
        """
        Finalize the FunctionPassManager.  Returns True if it produced
        any changes (?).
        """
        return ffi.lib.LLVMPY_FinalizeFunctionPassManager(self)

    def run(self, function, remarks_file=None, remarks_format='yaml',
            remarks_filter=''):
        """
        Run optimization passes on the given function.

        Parameters
        ----------
        function : llvmlite.binding.FunctionRef
            The function to be optimized inplace
        remarks_file : str; optional
            If not `None`, it is the file to store the optimization remarks.
        remarks_format : str; optional
            The format of the remarks file; the default is YAML
        remarks_filter : str; optional
            The filter that should be applied to the remarks output.
        """
        if remarks_file is None:
            return ffi.lib.LLVMPY_RunFunctionPassManager(self, function)
        else:
            r = ffi.lib.LLVMPY_RunFunctionPassManagerWithRemarks(
                self, function, _encode_string(remarks_format),
                _encode_string(remarks_filter),
                _encode_string(remarks_file))
            if r == -1:
                raise IOError("Failed to initialize remarks file.")
            return bool(r)

    def run_with_remarks(self, function, remarks_format='yaml',
                         remarks_filter=''):
        """
        Run optimization passes on the given function and returns the result
        and the remarks data.

        Parameters
        ----------
        function : llvmlite.binding.FunctionRef
            The function to be optimized inplace
        remarks_format : str; optional
            The format of the remarks file; the default is YAML
        remarks_filter : str; optional
            The filter that should be applied to the remarks output.
        """
        # LLVM is going to need to close this file and then reopen it, so we
        # can't use an unlinked temporary file.
        remarkdesc, remarkfile = mkstemp()
        try:
            # We get an open handle, but we need LLVM to write first, so close
            # it.
            with os.fdopen(remarkdesc, 'r'):
                pass
            r = self.run(function, remarkfile, remarks_format, remarks_filter)
            if r == -1:
                raise IOError("Failed to initialize remarks file.")
            with open(remarkfile) as f:
                return bool(r), f.read()
        finally:
            os.unlink(remarkfile)


# ============================================================================
# FFI

ffi.lib.LLVMPY_CreatePassManager.restype = ffi.LLVMPassManagerRef

ffi.lib.LLVMPY_CreateFunctionPassManager.argtypes = [ffi.LLVMModuleRef]
ffi.lib.LLVMPY_CreateFunctionPassManager.restype = ffi.LLVMPassManagerRef

ffi.lib.LLVMPY_DisposePassManager.argtypes = [ffi.LLVMPassManagerRef]

ffi.lib.LLVMPY_RunPassManager.argtypes = [ffi.LLVMPassManagerRef,
                                          ffi.LLVMModuleRef]
ffi.lib.LLVMPY_RunPassManager.restype = c_bool

ffi.lib.LLVMPY_RunPassManagerWithRemarks.argtypes = [ffi.LLVMPassManagerRef,
                                                     ffi.LLVMModuleRef,
                                                     c_char_p,
                                                     c_char_p,
                                                     c_char_p]
ffi.lib.LLVMPY_RunPassManagerWithRemarks.restype = c_int

ffi.lib.LLVMPY_InitializeFunctionPassManager.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_InitializeFunctionPassManager.restype = c_bool

ffi.lib.LLVMPY_FinalizeFunctionPassManager.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_FinalizeFunctionPassManager.restype = c_bool

ffi.lib.LLVMPY_RunFunctionPassManager.argtypes = [ffi.LLVMPassManagerRef,
                                                  ffi.LLVMValueRef]
ffi.lib.LLVMPY_RunFunctionPassManager.restype = c_bool

ffi.lib.LLVMPY_RunFunctionPassManagerWithRemarks.argtypes = [
    ffi.LLVMPassManagerRef, ffi.LLVMValueRef, c_char_p, c_char_p, c_char_p
]
ffi.lib.LLVMPY_RunFunctionPassManagerWithRemarks.restype = c_int

ffi.lib.LLVMPY_AddAAEvalPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddBasicAAWrapperPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddConstantMergePass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddDeadArgEliminationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddDependenceAnalysisPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddCallGraphDOTPrinterPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddCFGPrinterPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddDotDomPrinterPass.argtypes = [ffi.LLVMPassManagerRef, c_bool]
ffi.lib.LLVMPY_AddDotPostDomPrinterPass.argtypes = [
    ffi.LLVMPassManagerRef,
    c_bool]
ffi.lib.LLVMPY_AddGlobalsModRefAAPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddInstructionCountPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddIVUsersPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLazyValueInfoPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLintPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddModuleDebugInfoPrinterPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddRegionInfoPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddScalarEvolutionAAPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddAggressiveDCEPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddAlwaysInlinerPass.argtypes = [ffi.LLVMPassManagerRef, c_bool]
ffi.lib.LLVMPY_AddBreakCriticalEdgesPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddDeadStoreEliminationPass.argtypes = [
    ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddReversePostOrderFunctionAttrsPass.argtypes = [
    ffi.LLVMPassManagerRef]

if llvm_version_major < 16:
    ffi.lib.LLVMPY_AddAggressiveInstructionCombiningPass.argtypes = [
        ffi.LLVMPassManagerRef]

ffi.lib.LLVMPY_AddInternalizePass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLCSSAPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopDeletionPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopExtractorPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddSingleLoopExtractorPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopStrengthReducePass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopSimplificationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopUnrollPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopUnrollAndJamPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLoopUnswitchPass.argtypes = [ffi.LLVMPassManagerRef, c_bool,
                                               c_bool]
ffi.lib.LLVMPY_AddLowerAtomicPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLowerInvokePass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLowerSwitchPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddMemCpyOptimizationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddMergeFunctionsPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddMergeReturnsPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddPartialInliningPass.argtypes = [ffi.LLVMPassManagerRef]

if llvm_version_major < 16:
    ffi.lib.LLVMPY_AddPruneExceptionHandlingPass.argtypes = [
        ffi.LLVMPassManagerRef
    ]

ffi.lib.LLVMPY_AddReassociatePass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddDemoteRegisterToMemoryPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddSinkPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddStripSymbolsPass.argtypes = [ffi.LLVMPassManagerRef, c_bool]
ffi.lib.LLVMPY_AddStripDeadDebugInfoPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddStripDeadPrototypesPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddStripDebugDeclarePrototypesPass.argtypes = [
    ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddStripNondebugSymbolsPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddTailCallEliminationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddJumpThreadingPass.argtypes = [ffi.LLVMPassManagerRef, c_int]
ffi.lib.LLVMPY_AddFunctionAttrsPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddFunctionInliningPass.argtypes = [
    ffi.LLVMPassManagerRef, c_int]
ffi.lib.LLVMPY_AddGlobalDCEPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddGlobalOptimizerPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddIPSCCPPass.argtypes = [ffi.LLVMPassManagerRef]

ffi.lib.LLVMPY_AddDeadCodeEliminationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddCFGSimplificationPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddGVNPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddInstructionCombiningPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddLICMPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddSCCPPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddSROAPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddTypeBasedAliasAnalysisPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddBasicAliasAnalysisPass.argtypes = [ffi.LLVMPassManagerRef]
ffi.lib.LLVMPY_AddTargetLibraryInfoPass.argtypes = [ffi.LLVMPassManagerRef,
                                                    c_char_p]
ffi.lib.LLVMPY_AddInstructionNamerPass.argtypes = [ffi.LLVMPassManagerRef]

ffi.lib.LLVMPY_AddLegacyRefPrunePass.argtypes = [ffi.LLVMPassManagerRef, c_int,
                                                 c_size_t]

ffi.lib.LLVMPY_DumpRefPruneStats.argtypes = [POINTER(_c_PruneStats), c_bool]
