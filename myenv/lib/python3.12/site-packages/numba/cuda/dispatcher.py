import numpy as np
import os
import sys
import ctypes
import functools

from numba.core import config, serialize, sigutils, types, typing, utils
from numba.core.caching import Cache, CacheImpl
from numba.core.compiler_lock import global_compiler_lock
from numba.core.dispatcher import Dispatcher
from numba.core.errors import NumbaPerformanceWarning
from numba.core.typing.typeof import Purpose, typeof

from numba.cuda.api import get_current_device
from numba.cuda.args import wrap_arg
from numba.cuda.compiler import compile_cuda, CUDACompiler
from numba.cuda.cudadrv import driver
from numba.cuda.cudadrv.devices import get_context
from numba.cuda.descriptor import cuda_target
from numba.cuda.errors import (missing_launch_config_msg,
                               normalize_kernel_dimensions)
from numba.cuda import types as cuda_types

from numba import cuda
from numba import _dispatcher

from warnings import warn

cuda_fp16_math_funcs = ['hsin', 'hcos',
                        'hlog', 'hlog10',
                        'hlog2',
                        'hexp', 'hexp10',
                        'hexp2',
                        'hsqrt', 'hrsqrt',
                        'hfloor', 'hceil',
                        'hrcp', 'hrint',
                        'htrunc', 'hdiv']


class _Kernel(serialize.ReduceMixin):
    '''
    CUDA Kernel specialized for a given set of argument types. When called, this
    object launches the kernel on the device.
    '''

    @global_compiler_lock
    def __init__(self, py_func, argtypes, link=None, debug=False,
                 lineinfo=False, inline=False, fastmath=False, extensions=None,
                 max_registers=None, opt=True, device=False):

        if device:
            raise RuntimeError('Cannot compile a device function as a kernel')

        super().__init__()

        # _DispatcherBase.nopython_signatures() expects this attribute to be
        # present, because it assumes an overload is a CompileResult. In the
        # CUDA target, _Kernel instances are stored instead, so we provide this
        # attribute here to avoid duplicating nopython_signatures() in the CUDA
        # target with slight modifications.
        self.objectmode = False

        # The finalizer constructed by _DispatcherBase._make_finalizer also
        # expects overloads to be a CompileResult. It uses the entry_point to
        # remove a CompileResult from a target context. However, since we never
        # insert kernels into a target context (there is no need because they
        # cannot be called by other functions, only through the dispatcher) it
        # suffices to pretend we have an entry point of None.
        self.entry_point = None

        self.py_func = py_func
        self.argtypes = argtypes
        self.debug = debug
        self.lineinfo = lineinfo
        self.extensions = extensions or []

        nvvm_options = {
            'fastmath': fastmath,
            'opt': 3 if opt else 0
        }

        cc = get_current_device().compute_capability
        cres = compile_cuda(self.py_func, types.void, self.argtypes,
                            debug=self.debug,
                            lineinfo=lineinfo,
                            inline=inline,
                            fastmath=fastmath,
                            nvvm_options=nvvm_options,
                            cc=cc)
        tgt_ctx = cres.target_context
        code = self.py_func.__code__
        filename = code.co_filename
        linenum = code.co_firstlineno
        lib, kernel = tgt_ctx.prepare_cuda_kernel(cres.library, cres.fndesc,
                                                  debug, lineinfo, nvvm_options,
                                                  filename, linenum,
                                                  max_registers)

        if not link:
            link = []

        # A kernel needs cooperative launch if grid_sync is being used.
        self.cooperative = 'cudaCGGetIntrinsicHandle' in lib.get_asm_str()
        # We need to link against cudadevrt if grid sync is being used.
        if self.cooperative:
            lib.needs_cudadevrt = True

        res = [fn for fn in cuda_fp16_math_funcs
               if (f'__numba_wrapper_{fn}' in lib.get_asm_str())]

        if res:
            # Path to the source containing the foreign function
            basedir = os.path.dirname(os.path.abspath(__file__))
            functions_cu_path = os.path.join(basedir,
                                             'cpp_function_wrappers.cu')
            link.append(functions_cu_path)

        for filepath in link:
            lib.add_linking_file(filepath)

        # populate members
        self.entry_name = kernel.name
        self.signature = cres.signature
        self._type_annotation = cres.type_annotation
        self._codelibrary = lib
        self.call_helper = cres.call_helper

        # The following are referred to by the cache implementation. Note:
        # - There are no referenced environments in CUDA.
        # - Kernels don't have lifted code.
        # - reload_init is only for parfors.
        self.target_context = tgt_ctx
        self.fndesc = cres.fndesc
        self.environment = cres.environment
        self._referenced_environments = []
        self.lifted = []
        self.reload_init = []

    @property
    def library(self):
        return self._codelibrary

    @property
    def type_annotation(self):
        return self._type_annotation

    def _find_referenced_environments(self):
        return self._referenced_environments

    @property
    def codegen(self):
        return self.target_context.codegen()

    @property
    def argument_types(self):
        return tuple(self.signature.args)

    @classmethod
    def _rebuild(cls, cooperative, name, signature, codelibrary,
                 debug, lineinfo, call_helper, extensions):
        """
        Rebuild an instance.
        """
        instance = cls.__new__(cls)
        # invoke parent constructor
        super(cls, instance).__init__()
        # populate members
        instance.entry_point = None
        instance.cooperative = cooperative
        instance.entry_name = name
        instance.signature = signature
        instance._type_annotation = None
        instance._codelibrary = codelibrary
        instance.debug = debug
        instance.lineinfo = lineinfo
        instance.call_helper = call_helper
        instance.extensions = extensions
        return instance

    def _reduce_states(self):
        """
        Reduce the instance for serialization.
        Compiled definitions are serialized in PTX form.
        Type annotation are discarded.
        Thread, block and shared memory configuration are serialized.
        Stream information is discarded.
        """
        return dict(cooperative=self.cooperative, name=self.entry_name,
                    signature=self.signature, codelibrary=self._codelibrary,
                    debug=self.debug, lineinfo=self.lineinfo,
                    call_helper=self.call_helper, extensions=self.extensions)

    def bind(self):
        """
        Force binding to current CUDA context
        """
        self._codelibrary.get_cufunc()

    @property
    def regs_per_thread(self):
        '''
        The number of registers used by each thread for this kernel.
        '''
        return self._codelibrary.get_cufunc().attrs.regs

    @property
    def const_mem_size(self):
        '''
        The amount of constant memory used by this kernel.
        '''
        return self._codelibrary.get_cufunc().attrs.const

    @property
    def shared_mem_per_block(self):
        '''
        The amount of shared memory used per block for this kernel.
        '''
        return self._codelibrary.get_cufunc().attrs.shared

    @property
    def max_threads_per_block(self):
        '''
        The maximum allowable threads per block.
        '''
        return self._codelibrary.get_cufunc().attrs.maxthreads

    @property
    def local_mem_per_thread(self):
        '''
        The amount of local memory used per thread for this kernel.
        '''
        return self._codelibrary.get_cufunc().attrs.local

    def inspect_llvm(self):
        '''
        Returns the LLVM IR for this kernel.
        '''
        return self._codelibrary.get_llvm_str()

    def inspect_asm(self, cc):
        '''
        Returns the PTX code for this kernel.
        '''
        return self._codelibrary.get_asm_str(cc=cc)

    def inspect_sass_cfg(self):
        '''
        Returns the CFG of the SASS for this kernel.

        Requires nvdisasm to be available on the PATH.
        '''
        return self._codelibrary.get_sass_cfg()

    def inspect_sass(self):
        '''
        Returns the SASS code for this kernel.

        Requires nvdisasm to be available on the PATH.
        '''
        return self._codelibrary.get_sass()

    def inspect_types(self, file=None):
        '''
        Produce a dump of the Python source of this function annotated with the
        corresponding Numba IR and type information. The dump is written to
        *file*, or *sys.stdout* if *file* is *None*.
        '''
        if self._type_annotation is None:
            raise ValueError("Type annotation is not available")

        if file is None:
            file = sys.stdout

        print("%s %s" % (self.entry_name, self.argument_types), file=file)
        print('-' * 80, file=file)
        print(self._type_annotation, file=file)
        print('=' * 80, file=file)

    def max_cooperative_grid_blocks(self, blockdim, dynsmemsize=0):
        '''
        Calculates the maximum number of blocks that can be launched for this
        kernel in a cooperative grid in the current context, for the given block
        and dynamic shared memory sizes.

        :param blockdim: Block dimensions, either as a scalar for a 1D block, or
                         a tuple for 2D or 3D blocks.
        :param dynsmemsize: Dynamic shared memory size in bytes.
        :return: The maximum number of blocks in the grid.
        '''
        ctx = get_context()
        cufunc = self._codelibrary.get_cufunc()

        if isinstance(blockdim, tuple):
            blockdim = functools.reduce(lambda x, y: x * y, blockdim)
        active_per_sm = ctx.get_active_blocks_per_multiprocessor(cufunc,
                                                                 blockdim,
                                                                 dynsmemsize)
        sm_count = ctx.device.MULTIPROCESSOR_COUNT
        return active_per_sm * sm_count

    def launch(self, args, griddim, blockdim, stream=0, sharedmem=0):
        # Prepare kernel
        cufunc = self._codelibrary.get_cufunc()

        if self.debug:
            excname = cufunc.name + "__errcode__"
            excmem, excsz = cufunc.module.get_global_symbol(excname)
            assert excsz == ctypes.sizeof(ctypes.c_int)
            excval = ctypes.c_int()
            excmem.memset(0, stream=stream)

        # Prepare arguments
        retr = []                       # hold functors for writeback

        kernelargs = []
        for t, v in zip(self.argument_types, args):
            self._prepare_args(t, v, stream, retr, kernelargs)

        if driver.USE_NV_BINDING:
            zero_stream = driver.binding.CUstream(0)
        else:
            zero_stream = None

        stream_handle = stream and stream.handle or zero_stream

        # Invoke kernel
        driver.launch_kernel(cufunc.handle,
                             *griddim,
                             *blockdim,
                             sharedmem,
                             stream_handle,
                             kernelargs,
                             cooperative=self.cooperative)

        if self.debug:
            driver.device_to_host(ctypes.addressof(excval), excmem, excsz)
            if excval.value != 0:
                # An error occurred
                def load_symbol(name):
                    mem, sz = cufunc.module.get_global_symbol("%s__%s__" %
                                                              (cufunc.name,
                                                               name))
                    val = ctypes.c_int()
                    driver.device_to_host(ctypes.addressof(val), mem, sz)
                    return val.value

                tid = [load_symbol("tid" + i) for i in 'zyx']
                ctaid = [load_symbol("ctaid" + i) for i in 'zyx']
                code = excval.value
                exccls, exc_args, loc = self.call_helper.get_exception(code)
                # Prefix the exception message with the source location
                if loc is None:
                    locinfo = ''
                else:
                    sym, filepath, lineno = loc
                    filepath = os.path.abspath(filepath)
                    locinfo = 'In function %r, file %s, line %s, ' % (sym,
                                                                      filepath,
                                                                      lineno,)
                # Prefix the exception message with the thread position
                prefix = "%stid=%s ctaid=%s" % (locinfo, tid, ctaid)
                if exc_args:
                    exc_args = ("%s: %s" % (prefix, exc_args[0]),) + \
                        exc_args[1:]
                else:
                    exc_args = prefix,
                raise exccls(*exc_args)

        # retrieve auto converted arrays
        for wb in retr:
            wb()

    def _prepare_args(self, ty, val, stream, retr, kernelargs):
        """
        Convert arguments to ctypes and append to kernelargs
        """

        # map the arguments using any extension you've registered
        for extension in reversed(self.extensions):
            ty, val = extension.prepare_args(
                ty,
                val,
                stream=stream,
                retr=retr)

        if isinstance(ty, types.Array):
            devary = wrap_arg(val).to_device(retr, stream)

            c_intp = ctypes.c_ssize_t

            meminfo = ctypes.c_void_p(0)
            parent = ctypes.c_void_p(0)
            nitems = c_intp(devary.size)
            itemsize = c_intp(devary.dtype.itemsize)

            ptr = driver.device_pointer(devary)

            if driver.USE_NV_BINDING:
                ptr = int(ptr)

            data = ctypes.c_void_p(ptr)

            kernelargs.append(meminfo)
            kernelargs.append(parent)
            kernelargs.append(nitems)
            kernelargs.append(itemsize)
            kernelargs.append(data)
            for ax in range(devary.ndim):
                kernelargs.append(c_intp(devary.shape[ax]))
            for ax in range(devary.ndim):
                kernelargs.append(c_intp(devary.strides[ax]))

        elif isinstance(ty, types.Integer):
            cval = getattr(ctypes, "c_%s" % ty)(val)
            kernelargs.append(cval)

        elif ty == types.float16:
            cval = ctypes.c_uint16(np.float16(val).view(np.uint16))
            kernelargs.append(cval)

        elif ty == types.float64:
            cval = ctypes.c_double(val)
            kernelargs.append(cval)

        elif ty == types.float32:
            cval = ctypes.c_float(val)
            kernelargs.append(cval)

        elif ty == types.boolean:
            cval = ctypes.c_uint8(int(val))
            kernelargs.append(cval)

        elif ty == types.complex64:
            kernelargs.append(ctypes.c_float(val.real))
            kernelargs.append(ctypes.c_float(val.imag))

        elif ty == types.complex128:
            kernelargs.append(ctypes.c_double(val.real))
            kernelargs.append(ctypes.c_double(val.imag))

        elif isinstance(ty, (types.NPDatetime, types.NPTimedelta)):
            kernelargs.append(ctypes.c_int64(val.view(np.int64)))

        elif isinstance(ty, types.Record):
            devrec = wrap_arg(val).to_device(retr, stream)
            ptr = devrec.device_ctypes_pointer
            if driver.USE_NV_BINDING:
                ptr = ctypes.c_void_p(int(ptr))
            kernelargs.append(ptr)

        elif isinstance(ty, types.BaseTuple):
            assert len(ty) == len(val)
            for t, v in zip(ty, val):
                self._prepare_args(t, v, stream, retr, kernelargs)

        elif isinstance(ty, types.EnumMember):
            try:
                self._prepare_args(
                    ty.dtype, val.value, stream, retr, kernelargs
                )
            except NotImplementedError:
                raise NotImplementedError(ty, val)

        else:
            raise NotImplementedError(ty, val)


class ForAll(object):
    def __init__(self, dispatcher, ntasks, tpb, stream, sharedmem):
        if ntasks < 0:
            raise ValueError("Can't create ForAll with negative task count: %s"
                             % ntasks)
        self.dispatcher = dispatcher
        self.ntasks = ntasks
        self.thread_per_block = tpb
        self.stream = stream
        self.sharedmem = sharedmem

    def __call__(self, *args):
        if self.ntasks == 0:
            return

        if self.dispatcher.specialized:
            specialized = self.dispatcher
        else:
            specialized = self.dispatcher.specialize(*args)
        blockdim = self._compute_thread_per_block(specialized)
        griddim = (self.ntasks + blockdim - 1) // blockdim

        return specialized[griddim, blockdim, self.stream,
                           self.sharedmem](*args)

    def _compute_thread_per_block(self, dispatcher):
        tpb = self.thread_per_block
        # Prefer user-specified config
        if tpb != 0:
            return tpb
        # Else, ask the driver to give a good config
        else:
            ctx = get_context()
            # Dispatcher is specialized, so there's only one definition - get
            # it so we can get the cufunc from the code library
            kernel = next(iter(dispatcher.overloads.values()))
            kwargs = dict(
                func=kernel._codelibrary.get_cufunc(),
                b2d_func=0,     # dynamic-shared memory is constant to blksz
                memsize=self.sharedmem,
                blocksizelimit=1024,
            )
            _, tpb = ctx.get_max_potential_block_size(**kwargs)
            return tpb


class _LaunchConfiguration:
    def __init__(self, dispatcher, griddim, blockdim, stream, sharedmem):
        self.dispatcher = dispatcher
        self.griddim = griddim
        self.blockdim = blockdim
        self.stream = stream
        self.sharedmem = sharedmem

        if config.CUDA_LOW_OCCUPANCY_WARNINGS:
            # Warn when the grid has fewer than 128 blocks. This number is
            # chosen somewhat heuristically - ideally the minimum is 2 times
            # the number of SMs, but the number of SMs varies between devices -
            # some very small GPUs might only have 4 SMs, but an H100-SXM5 has
            # 132. In general kernels should be launched with large grids
            # (hundreds or thousands of blocks), so warning when fewer than 128
            # blocks are used will likely catch most beginner errors, where the
            # grid tends to be very small (single-digit or low tens of blocks).
            min_grid_size = 128
            grid_size = griddim[0] * griddim[1] * griddim[2]
            if grid_size < min_grid_size:
                msg = (f"Grid size {grid_size} will likely result in GPU "
                       "under-utilization due to low occupancy.")
                warn(NumbaPerformanceWarning(msg))

    def __call__(self, *args):
        return self.dispatcher.call(args, self.griddim, self.blockdim,
                                    self.stream, self.sharedmem)


class CUDACacheImpl(CacheImpl):
    def reduce(self, kernel):
        return kernel._reduce_states()

    def rebuild(self, target_context, payload):
        return _Kernel._rebuild(**payload)

    def check_cachable(self, cres):
        # CUDA Kernels are always cachable - the reasons for an entity not to
        # be cachable are:
        #
        # - The presence of lifted loops, or
        # - The presence of dynamic globals.
        #
        # neither of which apply to CUDA kernels.
        return True


class CUDACache(Cache):
    """
    Implements a cache that saves and loads CUDA kernels and compile results.
    """
    _impl_class = CUDACacheImpl

    def load_overload(self, sig, target_context):
        # Loading an overload refreshes the context to ensure it is
        # initialized. To initialize the correct (i.e. CUDA) target, we need to
        # enforce that the current target is the CUDA target.
        from numba.core.target_extension import target_override
        with target_override('cuda'):
            return super().load_overload(sig, target_context)


class CUDADispatcher(Dispatcher, serialize.ReduceMixin):
    '''
    CUDA Dispatcher object. When configured and called, the dispatcher will
    specialize itself for the given arguments (if no suitable specialized
    version already exists) & compute capability, and launch on the device
    associated with the current context.

    Dispatcher objects are not to be constructed by the user, but instead are
    created using the :func:`numba.cuda.jit` decorator.
    '''

    # Whether to fold named arguments and default values. Default values are
    # presently unsupported on CUDA, so we can leave this as False in all
    # cases.
    _fold_args = False

    targetdescr = cuda_target

    def __init__(self, py_func, targetoptions, pipeline_class=CUDACompiler):
        super().__init__(py_func, targetoptions=targetoptions,
                         pipeline_class=pipeline_class)

        # The following properties are for specialization of CUDADispatchers. A
        # specialized CUDADispatcher is one that is compiled for exactly one
        # set of argument types, and bypasses some argument type checking for
        # faster kernel launches.

        # Is this a specialized dispatcher?
        self._specialized = False

        # If we produced specialized dispatchers, we cache them for each set of
        # argument types
        self.specializations = {}

    @property
    def _numba_type_(self):
        return cuda_types.CUDADispatcher(self)

    def enable_caching(self):
        self._cache = CUDACache(self.py_func)

    @functools.lru_cache(maxsize=128)
    def configure(self, griddim, blockdim, stream=0, sharedmem=0):
        griddim, blockdim = normalize_kernel_dimensions(griddim, blockdim)
        return _LaunchConfiguration(self, griddim, blockdim, stream, sharedmem)

    def __getitem__(self, args):
        if len(args) not in [2, 3, 4]:
            raise ValueError('must specify at least the griddim and blockdim')
        return self.configure(*args)

    def forall(self, ntasks, tpb=0, stream=0, sharedmem=0):
        """Returns a 1D-configured dispatcher for a given number of tasks.

        This assumes that:

        - the kernel maps the Global Thread ID ``cuda.grid(1)`` to tasks on a
          1-1 basis.
        - the kernel checks that the Global Thread ID is upper-bounded by
          ``ntasks``, and does nothing if it is not.

        :param ntasks: The number of tasks.
        :param tpb: The size of a block. An appropriate value is chosen if this
                    parameter is not supplied.
        :param stream: The stream on which the configured dispatcher will be
                       launched.
        :param sharedmem: The number of bytes of dynamic shared memory required
                          by the kernel.
        :return: A configured dispatcher, ready to launch on a set of
                 arguments."""

        return ForAll(self, ntasks, tpb=tpb, stream=stream, sharedmem=sharedmem)

    @property
    def extensions(self):
        '''
        A list of objects that must have a `prepare_args` function. When a
        specialized kernel is called, each argument will be passed through
        to the `prepare_args` (from the last object in this list to the
        first). The arguments to `prepare_args` are:

        - `ty` the numba type of the argument
        - `val` the argument value itself
        - `stream` the CUDA stream used for the current call to the kernel
        - `retr` a list of zero-arg functions that you may want to append
          post-call cleanup work to.

        The `prepare_args` function must return a tuple `(ty, val)`, which
        will be passed in turn to the next right-most `extension`. After all
        the extensions have been called, the resulting `(ty, val)` will be
        passed into Numba's default argument marshalling logic.
        '''
        return self.targetoptions.get('extensions')

    def __call__(self, *args, **kwargs):
        # An attempt to launch an unconfigured kernel
        raise ValueError(missing_launch_config_msg)

    def call(self, args, griddim, blockdim, stream, sharedmem):
        '''
        Compile if necessary and invoke this kernel with *args*.
        '''
        if self.specialized:
            kernel = next(iter(self.overloads.values()))
        else:
            kernel = _dispatcher.Dispatcher._cuda_call(self, *args)

        kernel.launch(args, griddim, blockdim, stream, sharedmem)

    def _compile_for_args(self, *args, **kws):
        # Based on _DispatcherBase._compile_for_args.
        assert not kws
        argtypes = [self.typeof_pyval(a) for a in args]
        return self.compile(tuple(argtypes))

    def typeof_pyval(self, val):
        # Based on _DispatcherBase.typeof_pyval, but differs from it to support
        # the CUDA Array Interface.
        try:
            return typeof(val, Purpose.argument)
        except ValueError:
            if cuda.is_cuda_array(val):
                # When typing, we don't need to synchronize on the array's
                # stream - this is done when the kernel is launched.
                return typeof(cuda.as_cuda_array(val, sync=False),
                              Purpose.argument)
            else:
                raise

    def specialize(self, *args):
        '''
        Create a new instance of this dispatcher specialized for the given
        *args*.
        '''
        cc = get_current_device().compute_capability
        argtypes = tuple(
            [self.typingctx.resolve_argument_type(a) for a in args])
        if self.specialized:
            raise RuntimeError('Dispatcher already specialized')

        specialization = self.specializations.get((cc, argtypes))
        if specialization:
            return specialization

        targetoptions = self.targetoptions
        specialization = CUDADispatcher(self.py_func,
                                        targetoptions=targetoptions)
        specialization.compile(argtypes)
        specialization.disable_compile()
        specialization._specialized = True
        self.specializations[cc, argtypes] = specialization
        return specialization

    @property
    def specialized(self):
        """
        True if the Dispatcher has been specialized.
        """
        return self._specialized

    def get_regs_per_thread(self, signature=None):
        '''
        Returns the number of registers used by each thread in this kernel for
        the device in the current context.

        :param signature: The signature of the compiled kernel to get register
                          usage for. This may be omitted for a specialized
                          kernel.
        :return: The number of registers used by the compiled variant of the
                 kernel for the given signature and current device.
        '''
        if signature is not None:
            return self.overloads[signature.args].regs_per_thread
        if self.specialized:
            return next(iter(self.overloads.values())).regs_per_thread
        else:
            return {sig: overload.regs_per_thread
                    for sig, overload in self.overloads.items()}

    def get_const_mem_size(self, signature=None):
        '''
        Returns the size in bytes of constant memory used by this kernel for
        the device in the current context.

        :param signature: The signature of the compiled kernel to get constant
                          memory usage for. This may be omitted for a
                          specialized kernel.
        :return: The size in bytes of constant memory allocated by the
                 compiled variant of the kernel for the given signature and
                 current device.
        '''
        if signature is not None:
            return self.overloads[signature.args].const_mem_size
        if self.specialized:
            return next(iter(self.overloads.values())).const_mem_size
        else:
            return {sig: overload.const_mem_size
                    for sig, overload in self.overloads.items()}

    def get_shared_mem_per_block(self, signature=None):
        '''
        Returns the size in bytes of statically allocated shared memory
        for this kernel.

        :param signature: The signature of the compiled kernel to get shared
                          memory usage for. This may be omitted for a
                          specialized kernel.
        :return: The amount of shared memory allocated by the compiled variant
                 of the kernel for the given signature and current device.
        '''
        if signature is not None:
            return self.overloads[signature.args].shared_mem_per_block
        if self.specialized:
            return next(iter(self.overloads.values())).shared_mem_per_block
        else:
            return {sig: overload.shared_mem_per_block
                    for sig, overload in self.overloads.items()}

    def get_max_threads_per_block(self, signature=None):
        '''
        Returns the maximum allowable number of threads per block
        for this kernel. Exceeding this threshold will result in
        the kernel failing to launch.

        :param signature: The signature of the compiled kernel to get the max
                          threads per block for. This may be omitted for a
                          specialized kernel.
        :return: The maximum allowable threads per block for the compiled
                 variant of the kernel for the given signature and current
                 device.
        '''
        if signature is not None:
            return self.overloads[signature.args].max_threads_per_block
        if self.specialized:
            return next(iter(self.overloads.values())).max_threads_per_block
        else:
            return {sig: overload.max_threads_per_block
                    for sig, overload in self.overloads.items()}

    def get_local_mem_per_thread(self, signature=None):
        '''
        Returns the size in bytes of local memory per thread
        for this kernel.

        :param signature: The signature of the compiled kernel to get local
                          memory usage for. This may be omitted for a
                          specialized kernel.
        :return: The amount of local memory allocated by the compiled variant
                 of the kernel for the given signature and current device.
        '''
        if signature is not None:
            return self.overloads[signature.args].local_mem_per_thread
        if self.specialized:
            return next(iter(self.overloads.values())).local_mem_per_thread
        else:
            return {sig: overload.local_mem_per_thread
                    for sig, overload in self.overloads.items()}

    def get_call_template(self, args, kws):
        # Originally copied from _DispatcherBase.get_call_template. This
        # version deviates slightly from the _DispatcherBase version in order
        # to force casts when calling device functions. See e.g.
        # TestDeviceFunc.test_device_casting, added in PR #7496.
        """
        Get a typing.ConcreteTemplate for this dispatcher and the given
        *args* and *kws* types.  This allows resolution of the return type.

        A (template, pysig, args, kws) tuple is returned.
        """
        # Ensure an exactly-matching overload is available if we can
        # compile. We proceed with the typing even if we can't compile
        # because we may be able to force a cast on the caller side.
        if self._can_compile:
            self.compile_device(tuple(args))

        # Create function type for typing
        func_name = self.py_func.__name__
        name = "CallTemplate({0})".format(func_name)

        call_template = typing.make_concrete_template(
            name, key=func_name, signatures=self.nopython_signatures)
        pysig = utils.pysignature(self.py_func)

        return call_template, pysig, args, kws

    def compile_device(self, args, return_type=None):
        """Compile the device function for the given argument types.

        Each signature is compiled once by caching the compiled function inside
        this object.

        Returns the `CompileResult`.
        """
        if args not in self.overloads:
            with self._compiling_counter:

                debug = self.targetoptions.get('debug')
                lineinfo = self.targetoptions.get('lineinfo')
                inline = self.targetoptions.get('inline')
                fastmath = self.targetoptions.get('fastmath')

                nvvm_options = {
                    'opt': 3 if self.targetoptions.get('opt') else 0,
                    'fastmath': fastmath
                }

                cc = get_current_device().compute_capability
                cres = compile_cuda(self.py_func, return_type, args,
                                    debug=debug,
                                    lineinfo=lineinfo,
                                    inline=inline,
                                    fastmath=fastmath,
                                    nvvm_options=nvvm_options,
                                    cc=cc)
                self.overloads[args] = cres

                cres.target_context.insert_user_function(cres.entry_point,
                                                         cres.fndesc,
                                                         [cres.library])
        else:
            cres = self.overloads[args]

        return cres

    def add_overload(self, kernel, argtypes):
        c_sig = [a._code for a in argtypes]
        self._insert(c_sig, kernel, cuda=True)
        self.overloads[argtypes] = kernel

    def compile(self, sig):
        '''
        Compile and bind to the current context a version of this kernel
        specialized for the given signature.
        '''
        argtypes, return_type = sigutils.normalize_signature(sig)
        assert return_type is None or return_type == types.none

        # Do we already have an in-memory compiled kernel?
        if self.specialized:
            return next(iter(self.overloads.values()))
        else:
            kernel = self.overloads.get(argtypes)
            if kernel is not None:
                return kernel

        # Can we load from the disk cache?
        kernel = self._cache.load_overload(sig, self.targetctx)

        if kernel is not None:
            self._cache_hits[sig] += 1
        else:
            # We need to compile a new kernel
            self._cache_misses[sig] += 1
            if not self._can_compile:
                raise RuntimeError("Compilation disabled")

            kernel = _Kernel(self.py_func, argtypes, **self.targetoptions)
            # We call bind to force codegen, so that there is a cubin to cache
            kernel.bind()
            self._cache.save_overload(sig, kernel)

        self.add_overload(kernel, argtypes)

        return kernel

    def inspect_llvm(self, signature=None):
        '''
        Return the LLVM IR for this kernel.

        :param signature: A tuple of argument types.
        :return: The LLVM IR for the given signature, or a dict of LLVM IR
                 for all previously-encountered signatures.

        '''
        device = self.targetoptions.get('device')
        if signature is not None:
            if device:
                return self.overloads[signature].library.get_llvm_str()
            else:
                return self.overloads[signature].inspect_llvm()
        else:
            if device:
                return {sig: overload.library.get_llvm_str()
                        for sig, overload in self.overloads.items()}
            else:
                return {sig: overload.inspect_llvm()
                        for sig, overload in self.overloads.items()}

    def inspect_asm(self, signature=None):
        '''
        Return this kernel's PTX assembly code for for the device in the
        current context.

        :param signature: A tuple of argument types.
        :return: The PTX code for the given signature, or a dict of PTX codes
                 for all previously-encountered signatures.
        '''
        cc = get_current_device().compute_capability
        device = self.targetoptions.get('device')
        if signature is not None:
            if device:
                return self.overloads[signature].library.get_asm_str(cc)
            else:
                return self.overloads[signature].inspect_asm(cc)
        else:
            if device:
                return {sig: overload.library.get_asm_str(cc)
                        for sig, overload in self.overloads.items()}
            else:
                return {sig: overload.inspect_asm(cc)
                        for sig, overload in self.overloads.items()}

    def inspect_sass_cfg(self, signature=None):
        '''
        Return this kernel's CFG for the device in the current context.

        :param signature: A tuple of argument types.
        :return: The CFG for the given signature, or a dict of CFGs
                 for all previously-encountered signatures.

        The CFG for the device in the current context is returned.

        Requires nvdisasm to be available on the PATH.
        '''
        if self.targetoptions.get('device'):
            raise RuntimeError('Cannot get the CFG of a device function')

        if signature is not None:
            return self.overloads[signature].inspect_sass_cfg()
        else:
            return {sig: defn.inspect_sass_cfg()
                    for sig, defn in self.overloads.items()}

    def inspect_sass(self, signature=None):
        '''
        Return this kernel's SASS assembly code for for the device in the
        current context.

        :param signature: A tuple of argument types.
        :return: The SASS code for the given signature, or a dict of SASS codes
                 for all previously-encountered signatures.

        SASS for the device in the current context is returned.

        Requires nvdisasm to be available on the PATH.
        '''
        if self.targetoptions.get('device'):
            raise RuntimeError('Cannot inspect SASS of a device function')

        if signature is not None:
            return self.overloads[signature].inspect_sass()
        else:
            return {sig: defn.inspect_sass()
                    for sig, defn in self.overloads.items()}

    def inspect_types(self, file=None):
        '''
        Produce a dump of the Python source of this function annotated with the
        corresponding Numba IR and type information. The dump is written to
        *file*, or *sys.stdout* if *file* is *None*.
        '''
        if file is None:
            file = sys.stdout

        for _, defn in self.overloads.items():
            defn.inspect_types(file=file)

    @classmethod
    def _rebuild(cls, py_func, targetoptions):
        """
        Rebuild an instance.
        """
        instance = cls(py_func, targetoptions)
        return instance

    def _reduce_states(self):
        """
        Reduce the instance for serialization.
        Compiled definitions are discarded.
        """
        return dict(py_func=self.py_func,
                    targetoptions=self.targetoptions)
