from llvmlite import ir
from numba.core.typing.templates import ConcreteTemplate
from numba.core import types, typing, funcdesc, config, compiler, sigutils
from numba.core.compiler import (sanitize_compile_result_entries, CompilerBase,
                                 DefaultPassBuilder, Flags, Option,
                                 CompileResult)
from numba.core.compiler_lock import global_compiler_lock
from numba.core.compiler_machinery import (LoweringPass,
                                           PassManager, register_pass)
from numba.core.errors import NumbaInvalidConfigWarning
from numba.core.typed_passes import (IRLegalization, NativeLowering,
                                     AnnotateTypes)
from warnings import warn
from numba.cuda.api import get_current_device
from numba.cuda.target import CUDACABICallConv


def _nvvm_options_type(x):
    if x is None:
        return None

    else:
        assert isinstance(x, dict)
        return x


class CUDAFlags(Flags):
    nvvm_options = Option(
        type=_nvvm_options_type,
        default=None,
        doc="NVVM options",
    )
    compute_capability = Option(
        type=tuple,
        default=None,
        doc="Compute Capability",
    )


# The CUDACompileResult (CCR) has a specially-defined entry point equal to its
# id.  This is because the entry point is used as a key into a dict of
# overloads by the base dispatcher. The id of the CCR is the only small and
# unique property of a CompileResult in the CUDA target (cf. the CPU target,
# which uses its entry_point, which is a pointer value).
#
# This does feel a little hackish, and there are two ways in which this could
# be improved:
#
# 1. We could change the core of Numba so that each CompileResult has its own
#    unique ID that can be used as a key - e.g. a count, similar to the way in
#    which types have unique counts.
# 2. At some future time when kernel launch uses a compiled function, the entry
#    point will no longer need to be a synthetic value, but will instead be a
#    pointer to the compiled function as in the CPU target.

class CUDACompileResult(CompileResult):
    @property
    def entry_point(self):
        return id(self)


def cuda_compile_result(**entries):
    entries = sanitize_compile_result_entries(entries)
    return CUDACompileResult(**entries)


@register_pass(mutates_CFG=True, analysis_only=False)
class CUDABackend(LoweringPass):

    _name = "cuda_backend"

    def __init__(self):
        LoweringPass.__init__(self)

    def run_pass(self, state):
        """
        Back-end: Packages lowering output in a compile result
        """
        lowered = state['cr']
        signature = typing.signature(state.return_type, *state.args)

        state.cr = cuda_compile_result(
            typing_context=state.typingctx,
            target_context=state.targetctx,
            typing_error=state.status.fail_reason,
            type_annotation=state.type_annotation,
            library=state.library,
            call_helper=lowered.call_helper,
            signature=signature,
            fndesc=lowered.fndesc,
        )
        return True


@register_pass(mutates_CFG=False, analysis_only=False)
class CreateLibrary(LoweringPass):
    """
    Create a CUDACodeLibrary for the NativeLowering pass to populate. The
    NativeLowering pass will create a code library if none exists, but we need
    to set it up with nvvm_options from the flags if they are present.
    """

    _name = "create_library"

    def __init__(self):
        LoweringPass.__init__(self)

    def run_pass(self, state):
        codegen = state.targetctx.codegen()
        name = state.func_id.func_qualname
        nvvm_options = state.flags.nvvm_options
        state.library = codegen.create_library(name, nvvm_options=nvvm_options)
        # Enable object caching upfront so that the library can be serialized.
        state.library.enable_object_caching()

        return True


class CUDACompiler(CompilerBase):
    def define_pipelines(self):
        dpb = DefaultPassBuilder
        pm = PassManager('cuda')

        untyped_passes = dpb.define_untyped_pipeline(self.state)
        pm.passes.extend(untyped_passes.passes)

        typed_passes = dpb.define_typed_pipeline(self.state)
        pm.passes.extend(typed_passes.passes)

        lowering_passes = self.define_cuda_lowering_pipeline(self.state)
        pm.passes.extend(lowering_passes.passes)

        pm.finalize()
        return [pm]

    def define_cuda_lowering_pipeline(self, state):
        pm = PassManager('cuda_lowering')
        # legalise
        pm.add_pass(IRLegalization,
                    "ensure IR is legal prior to lowering")
        pm.add_pass(AnnotateTypes, "annotate types")

        # lower
        pm.add_pass(CreateLibrary, "create library")
        pm.add_pass(NativeLowering, "native lowering")
        pm.add_pass(CUDABackend, "cuda backend")

        pm.finalize()
        return pm


@global_compiler_lock
def compile_cuda(pyfunc, return_type, args, debug=False, lineinfo=False,
                 inline=False, fastmath=False, nvvm_options=None,
                 cc=None):
    if cc is None:
        raise ValueError('Compute Capability must be supplied')

    from .descriptor import cuda_target
    typingctx = cuda_target.typing_context
    targetctx = cuda_target.target_context

    flags = CUDAFlags()
    # Do not compile (generate native code), just lower (to LLVM)
    flags.no_compile = True
    flags.no_cpython_wrapper = True
    flags.no_cfunc_wrapper = True

    # Both debug and lineinfo turn on debug information in the compiled code,
    # but we keep them separate arguments in case we later want to overload
    # some other behavior on the debug flag. In particular, -opt=3 is not
    # supported with debug enabled, and enabling only lineinfo should not
    # affect the error model.
    if debug or lineinfo:
        flags.debuginfo = True

    if lineinfo:
        flags.dbg_directives_only = True

    if debug:
        flags.error_model = 'python'
    else:
        flags.error_model = 'numpy'

    if inline:
        flags.forceinline = True
    if fastmath:
        flags.fastmath = True
    if nvvm_options:
        flags.nvvm_options = nvvm_options
    flags.compute_capability = cc

    # Run compilation pipeline
    from numba.core.target_extension import target_override
    with target_override('cuda'):
        cres = compiler.compile_extra(typingctx=typingctx,
                                      targetctx=targetctx,
                                      func=pyfunc,
                                      args=args,
                                      return_type=return_type,
                                      flags=flags,
                                      locals={},
                                      pipeline_class=CUDACompiler)

    library = cres.library
    library.finalize()

    return cres


def cabi_wrap_function(context, lib, fndesc, wrapper_function_name,
                       nvvm_options):
    """
    Wrap a Numba ABI function in a C ABI wrapper at the NVVM IR level.

    The C ABI wrapper will have the same name as the source Python function.
    """
    # The wrapper will be contained in a new library that links to the wrapped
    # function's library
    library = lib.codegen.create_library(f'{lib.name}_function_',
                                         entry_name=wrapper_function_name,
                                         nvvm_options=nvvm_options)
    library.add_linking_library(lib)

    # Determine the caller (C ABI) and wrapper (Numba ABI) function types
    argtypes = fndesc.argtypes
    restype = fndesc.restype
    c_call_conv = CUDACABICallConv(context)
    wrapfnty = c_call_conv.get_function_type(restype, argtypes)
    fnty = context.call_conv.get_function_type(fndesc.restype, argtypes)

    # Create a new module and declare the callee
    wrapper_module = context.create_module("cuda.cabi.wrapper")
    func = ir.Function(wrapper_module, fnty, fndesc.llvm_func_name)

    # Define the caller - populate it with a call to the callee and return
    # its return value

    wrapfn = ir.Function(wrapper_module, wrapfnty, wrapper_function_name)
    builder = ir.IRBuilder(wrapfn.append_basic_block(''))

    arginfo = context.get_arg_packer(argtypes)
    callargs = arginfo.from_arguments(builder, wrapfn.args)
    # We get (status, return_value), but we ignore the status since we
    # can't propagate it through the C ABI anyway
    _, return_value = context.call_conv.call_function(
        builder, func, restype, argtypes, callargs)
    builder.ret(return_value)

    library.add_ir_module(wrapper_module)
    library.finalize()
    return library


@global_compiler_lock
def compile(pyfunc, sig, debug=False, lineinfo=False, device=True,
            fastmath=False, cc=None, opt=True, abi="c", abi_info=None,
            output='ptx'):
    """Compile a Python function to PTX or LTO-IR for a given set of argument
    types.

    :param pyfunc: The Python function to compile.
    :param sig: The signature representing the function's input and output
                types. If this is a tuple of argument types without a return
                type, the inferred return type is returned by this function. If
                a signature including a return type is passed, the compiled code
                will include a cast from the inferred return type to the
                specified return type, and this function will return the
                specified return type.
    :param debug: Whether to include debug info in the compiled code.
    :type debug: bool
    :param lineinfo: Whether to include a line mapping from the compiled code
                     to the source code. Usually this is used with optimized
                     code (since debug mode would automatically include this),
                     so we want debug info in the LLVM IR but only the line
                     mapping in the final output.
    :type lineinfo: bool
    :param device: Whether to compile a device function.
    :type device: bool
    :param fastmath: Whether to enable fast math flags (ftz=1, prec_sqrt=0,
                     prec_div=, and fma=1)
    :type fastmath: bool
    :param cc: Compute capability to compile for, as a tuple
               ``(MAJOR, MINOR)``. Defaults to ``(5, 0)``.
    :type cc: tuple
    :param opt: Enable optimizations. Defaults to ``True``.
    :type opt: bool
    :param abi: The ABI for a compiled function - either ``"numba"`` or
                ``"c"``. Note that the Numba ABI is not considered stable.
                The C ABI is only supported for device functions at present.
    :type abi: str
    :param abi_info: A dict of ABI-specific options. The ``"c"`` ABI supports
                     one option, ``"abi_name"``, for providing the wrapper
                     function's name. The ``"numba"`` ABI has no options.
    :type abi_info: dict
    :param output: Type of output to generate, either ``"ptx"`` or ``"ltoir"``.
    :type output: str
    :return: (code, resty): The compiled code and inferred return type
    :rtype: tuple
    """
    if abi not in ("numba", "c"):
        raise NotImplementedError(f'Unsupported ABI: {abi}')

    if abi == 'c' and not device:
        raise NotImplementedError('The C ABI is not supported for kernels')

    if output not in ("ptx", "ltoir"):
        raise NotImplementedError(f'Unsupported output type: {output}')

    if debug and opt:
        msg = ("debug=True with opt=True (the default) "
               "is not supported by CUDA. This may result in a crash"
               " - set debug=False or opt=False.")
        warn(NumbaInvalidConfigWarning(msg))

    lto = (output == 'ltoir')
    abi_info = abi_info or dict()

    nvvm_options = {
        'fastmath': fastmath,
        'opt': 3 if opt else 0
    }

    if lto:
        nvvm_options['gen-lto'] = None

    args, return_type = sigutils.normalize_signature(sig)

    cc = cc or config.CUDA_DEFAULT_PTX_CC
    cres = compile_cuda(pyfunc, return_type, args, debug=debug,
                        lineinfo=lineinfo, fastmath=fastmath,
                        nvvm_options=nvvm_options, cc=cc)
    resty = cres.signature.return_type

    if resty and not device and resty != types.void:
        raise TypeError("CUDA kernel must have void return type.")

    tgt = cres.target_context

    if device:
        lib = cres.library
        if abi == "c":
            wrapper_name = abi_info.get('abi_name', pyfunc.__name__)
            lib = cabi_wrap_function(tgt, lib, cres.fndesc, wrapper_name,
                                     nvvm_options)
    else:
        code = pyfunc.__code__
        filename = code.co_filename
        linenum = code.co_firstlineno

        lib, kernel = tgt.prepare_cuda_kernel(cres.library, cres.fndesc, debug,
                                              lineinfo, nvvm_options, filename,
                                              linenum)

    if lto:
        code = lib.get_ltoir(cc=cc)
    else:
        code = lib.get_asm_str(cc=cc)
    return code, resty


def compile_for_current_device(pyfunc, sig, debug=False, lineinfo=False,
                               device=True, fastmath=False, opt=True,
                               abi="c", abi_info=None, output='ptx'):
    """Compile a Python function to PTX or LTO-IR for a given signature for the
    current device's compute capabilility. This calls :func:`compile` with an
    appropriate ``cc`` value for the current device."""
    cc = get_current_device().compute_capability
    return compile(pyfunc, sig, debug=debug, lineinfo=lineinfo, device=device,
                   fastmath=fastmath, cc=cc, opt=opt, abi=abi,
                   abi_info=abi_info, output=output)


def compile_ptx(pyfunc, sig, debug=False, lineinfo=False, device=False,
                fastmath=False, cc=None, opt=True, abi="numba", abi_info=None):
    """Compile a Python function to PTX for a given signature. See
    :func:`compile`. The defaults for this function are to compile a kernel
    with the Numba ABI, rather than :func:`compile`'s default of compiling a
    device function with the C ABI."""
    return compile(pyfunc, sig, debug=debug, lineinfo=lineinfo, device=device,
                   fastmath=fastmath, cc=cc, opt=opt, abi=abi,
                   abi_info=abi_info, output='ptx')


def compile_ptx_for_current_device(pyfunc, sig, debug=False, lineinfo=False,
                                   device=False, fastmath=False, opt=True,
                                   abi="numba", abi_info=None):
    """Compile a Python function to PTX for a given signature for the current
    device's compute capabilility. See :func:`compile_ptx`."""
    cc = get_current_device().compute_capability
    return compile_ptx(pyfunc, sig, debug=debug, lineinfo=lineinfo,
                       device=device, fastmath=fastmath, cc=cc, opt=opt,
                       abi=abi, abi_info=abi_info)


def declare_device_function(name, restype, argtypes):
    return declare_device_function_template(name, restype, argtypes).key


def declare_device_function_template(name, restype, argtypes):
    from .descriptor import cuda_target
    typingctx = cuda_target.typing_context
    targetctx = cuda_target.target_context
    sig = typing.signature(restype, *argtypes)
    extfn = ExternFunction(name, sig)

    class device_function_template(ConcreteTemplate):
        key = extfn
        cases = [sig]

    fndesc = funcdesc.ExternalFunctionDescriptor(
        name=name, restype=restype, argtypes=argtypes)
    typingctx.insert_user_function(extfn, device_function_template)
    targetctx.insert_user_function(extfn, fndesc)

    return device_function_template


class ExternFunction(object):
    def __init__(self, name, sig):
        self.name = name
        self.sig = sig
