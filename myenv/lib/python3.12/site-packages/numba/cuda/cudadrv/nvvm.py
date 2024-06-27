"""
This is a direct translation of nvvm.h
"""
import logging
import re
import sys
import warnings
from ctypes import (c_void_p, c_int, POINTER, c_char_p, c_size_t, byref,
                    c_char)

import threading

from llvmlite import ir

from .error import NvvmError, NvvmSupportError, NvvmWarning
from .libs import get_libdevice, open_libdevice, open_cudalib
from numba.core import cgutils, config


logger = logging.getLogger(__name__)

ADDRSPACE_GENERIC = 0
ADDRSPACE_GLOBAL = 1
ADDRSPACE_SHARED = 3
ADDRSPACE_CONSTANT = 4
ADDRSPACE_LOCAL = 5

# Opaque handle for compilation unit
nvvm_program = c_void_p

# Result code
nvvm_result = c_int

RESULT_CODE_NAMES = '''
NVVM_SUCCESS
NVVM_ERROR_OUT_OF_MEMORY
NVVM_ERROR_PROGRAM_CREATION_FAILURE
NVVM_ERROR_IR_VERSION_MISMATCH
NVVM_ERROR_INVALID_INPUT
NVVM_ERROR_INVALID_PROGRAM
NVVM_ERROR_INVALID_IR
NVVM_ERROR_INVALID_OPTION
NVVM_ERROR_NO_MODULE_IN_PROGRAM
NVVM_ERROR_COMPILATION
'''.split()

for i, k in enumerate(RESULT_CODE_NAMES):
    setattr(sys.modules[__name__], k, i)

# Data layouts. NVVM IR 1.8 (CUDA 11.6) introduced 128-bit integer support.

_datalayout_original = ('e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-'
                        'i64:64:64-f32:32:32-f64:64:64-v16:16:16-v32:32:32-'
                        'v64:64:64-v128:128:128-n16:32:64')
_datalayout_i128 = ('e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-'
                    'i128:128:128-f32:32:32-f64:64:64-v16:16:16-v32:32:32-'
                    'v64:64:64-v128:128:128-n16:32:64')


def is_available():
    """
    Return if libNVVM is available
    """
    try:
        NVVM()
    except NvvmSupportError:
        return False
    else:
        return True


_nvvm_lock = threading.Lock()


class NVVM(object):
    '''Process-wide singleton.
    '''
    _PROTOTYPES = {

        # nvvmResult nvvmVersion(int *major, int *minor)
        'nvvmVersion': (nvvm_result, POINTER(c_int), POINTER(c_int)),

        # nvvmResult nvvmCreateProgram(nvvmProgram *cu)
        'nvvmCreateProgram': (nvvm_result, POINTER(nvvm_program)),

        # nvvmResult nvvmDestroyProgram(nvvmProgram *cu)
        'nvvmDestroyProgram': (nvvm_result, POINTER(nvvm_program)),

        # nvvmResult nvvmAddModuleToProgram(nvvmProgram cu, const char *buffer,
        #                                   size_t size, const char *name)
        'nvvmAddModuleToProgram': (
            nvvm_result, nvvm_program, c_char_p, c_size_t, c_char_p),

        # nvvmResult nvvmLazyAddModuleToProgram(nvvmProgram cu,
        #                                       const char* buffer,
        #                                       size_t size,
        #                                       const char *name)
        'nvvmLazyAddModuleToProgram': (
            nvvm_result, nvvm_program, c_char_p, c_size_t, c_char_p),

        # nvvmResult nvvmCompileProgram(nvvmProgram cu, int numOptions,
        #                          const char **options)
        'nvvmCompileProgram': (
            nvvm_result, nvvm_program, c_int, POINTER(c_char_p)),

        # nvvmResult nvvmGetCompiledResultSize(nvvmProgram cu,
        #                                      size_t *bufferSizeRet)
        'nvvmGetCompiledResultSize': (
            nvvm_result, nvvm_program, POINTER(c_size_t)),

        # nvvmResult nvvmGetCompiledResult(nvvmProgram cu, char *buffer)
        'nvvmGetCompiledResult': (nvvm_result, nvvm_program, c_char_p),

        # nvvmResult nvvmGetProgramLogSize(nvvmProgram cu,
        #                                      size_t *bufferSizeRet)
        'nvvmGetProgramLogSize': (nvvm_result, nvvm_program, POINTER(c_size_t)),

        # nvvmResult nvvmGetProgramLog(nvvmProgram cu, char *buffer)
        'nvvmGetProgramLog': (nvvm_result, nvvm_program, c_char_p),

        # nvvmResult nvvmIRVersion (int* majorIR, int* minorIR, int* majorDbg,
        #                           int* minorDbg )
        'nvvmIRVersion': (nvvm_result, POINTER(c_int), POINTER(c_int),
                          POINTER(c_int), POINTER(c_int)),
        # nvvmResult nvvmVerifyProgram (nvvmProgram prog, int numOptions,
        #                               const char** options)
        'nvvmVerifyProgram': (nvvm_result, nvvm_program, c_int,
                              POINTER(c_char_p))
    }

    # Singleton reference
    __INSTANCE = None

    def __new__(cls):
        with _nvvm_lock:
            if cls.__INSTANCE is None:
                cls.__INSTANCE = inst = object.__new__(cls)
                try:
                    inst.driver = open_cudalib('nvvm')
                except OSError as e:
                    cls.__INSTANCE = None
                    errmsg = ("libNVVM cannot be found. Do `conda install "
                              "cudatoolkit`:\n%s")
                    raise NvvmSupportError(errmsg % e)

                # Find & populate functions
                for name, proto in inst._PROTOTYPES.items():
                    func = getattr(inst.driver, name)
                    func.restype = proto[0]
                    func.argtypes = proto[1:]
                    setattr(inst, name, func)

        return cls.__INSTANCE

    def __init__(self):
        ir_versions = self.get_ir_version()
        self._majorIR = ir_versions[0]
        self._minorIR = ir_versions[1]
        self._majorDbg = ir_versions[2]
        self._minorDbg = ir_versions[3]
        self._supported_ccs = get_supported_ccs()

    @property
    def data_layout(self):
        if (self._majorIR, self._minorIR) < (1, 8):
            return _datalayout_original
        else:
            return _datalayout_i128

    @property
    def supported_ccs(self):
        return self._supported_ccs

    def get_version(self):
        major = c_int()
        minor = c_int()
        err = self.nvvmVersion(byref(major), byref(minor))
        self.check_error(err, 'Failed to get version.')
        return major.value, minor.value

    def get_ir_version(self):
        majorIR = c_int()
        minorIR = c_int()
        majorDbg = c_int()
        minorDbg = c_int()
        err = self.nvvmIRVersion(byref(majorIR), byref(minorIR),
                                 byref(majorDbg), byref(minorDbg))
        self.check_error(err, 'Failed to get IR version.')
        return majorIR.value, minorIR.value, majorDbg.value, minorDbg.value

    def check_error(self, error, msg, exit=False):
        if error:
            exc = NvvmError(msg, RESULT_CODE_NAMES[error])
            if exit:
                print(exc)
                sys.exit(1)
            else:
                raise exc


class CompilationUnit(object):
    def __init__(self):
        self.driver = NVVM()
        self._handle = nvvm_program()
        err = self.driver.nvvmCreateProgram(byref(self._handle))
        self.driver.check_error(err, 'Failed to create CU')

    def __del__(self):
        driver = NVVM()
        err = driver.nvvmDestroyProgram(byref(self._handle))
        driver.check_error(err, 'Failed to destroy CU', exit=True)

    def add_module(self, buffer):
        """
         Add a module level NVVM IR to a compilation unit.
         - The buffer should contain an NVVM module IR either in the bitcode
           representation (LLVM3.0) or in the text representation.
        """
        err = self.driver.nvvmAddModuleToProgram(self._handle, buffer,
                                                 len(buffer), None)
        self.driver.check_error(err, 'Failed to add module')

    def lazy_add_module(self, buffer):
        """
        Lazily add an NVVM IR module to a compilation unit.
        The buffer should contain NVVM module IR either in the bitcode
        representation or in the text representation.
        """
        err = self.driver.nvvmLazyAddModuleToProgram(self._handle, buffer,
                                                     len(buffer), None)
        self.driver.check_error(err, 'Failed to add module')

    def compile(self, **options):
        """Perform Compilation.

        Compilation options are accepted as keyword arguments, with the
        following considerations:

        - Underscores (`_`) in option names are converted to dashes (`-`), to
          match NVVM's option name format.
        - Options that take a value will be emitted in the form
          "-<name>=<value>".
        - Booleans passed as option values will be converted to integers.
        - Options which take no value (such as `-gen-lto`) should have a value
          of `None` passed in and will be emitted in the form "-<name>".

        For documentation on NVVM compilation options, see the CUDA Toolkit
        Documentation:

        https://docs.nvidia.com/cuda/libnvvm-api/index.html#_CPPv418nvvmCompileProgram11nvvmProgramiPPKc
        """

        def stringify_option(k, v):
            k = k.replace('_', '-')

            if v is None:
                return f'-{k}'

            if isinstance(v, bool):
                v = int(v)

            return f'-{k}={v}'

        options = [stringify_option(k, v) for k, v in options.items()]

        c_opts = (c_char_p * len(options))(*[c_char_p(x.encode('utf8'))
                                             for x in options])
        # verify
        err = self.driver.nvvmVerifyProgram(self._handle, len(options), c_opts)
        self._try_error(err, 'Failed to verify\n')

        # compile
        err = self.driver.nvvmCompileProgram(self._handle, len(options), c_opts)
        self._try_error(err, 'Failed to compile\n')

        # get result
        reslen = c_size_t()
        err = self.driver.nvvmGetCompiledResultSize(self._handle, byref(reslen))

        self._try_error(err, 'Failed to get size of compiled result.')

        output_buffer = (c_char * reslen.value)()
        err = self.driver.nvvmGetCompiledResult(self._handle, output_buffer)
        self._try_error(err, 'Failed to get compiled result.')

        # get log
        self.log = self.get_log()
        if self.log:
            warnings.warn(self.log, category=NvvmWarning)

        return output_buffer[:]

    def _try_error(self, err, msg):
        self.driver.check_error(err, "%s\n%s" % (msg, self.get_log()))

    def get_log(self):
        reslen = c_size_t()
        err = self.driver.nvvmGetProgramLogSize(self._handle, byref(reslen))
        self.driver.check_error(err, 'Failed to get compilation log size.')

        if reslen.value > 1:
            logbuf = (c_char * reslen.value)()
            err = self.driver.nvvmGetProgramLog(self._handle, logbuf)
            self.driver.check_error(err, 'Failed to get compilation log.')

            return logbuf.value.decode('utf8')  # populate log attribute

        return ''


COMPUTE_CAPABILITIES = (
    (3, 5), (3, 7),
    (5, 0), (5, 2), (5, 3),
    (6, 0), (6, 1), (6, 2),
    (7, 0), (7, 2), (7, 5),
    (8, 0), (8, 6), (8, 7), (8, 9),
    (9, 0)
)

# Maps CTK version -> (min supported cc, max supported cc) inclusive
CTK_SUPPORTED = {
    (11, 2): ((3, 5), (8, 6)),
    (11, 3): ((3, 5), (8, 6)),
    (11, 4): ((3, 5), (8, 7)),
    (11, 5): ((3, 5), (8, 7)),
    (11, 6): ((3, 5), (8, 7)),
    (11, 7): ((3, 5), (8, 7)),
    (11, 8): ((3, 5), (9, 0)),
    (12, 0): ((5, 0), (9, 0)),
    (12, 1): ((5, 0), (9, 0)),
    (12, 2): ((5, 0), (9, 0)),
    (12, 3): ((5, 0), (9, 0)),
    (12, 4): ((5, 0), (9, 0)),
}


def ccs_supported_by_ctk(ctk_version):
    try:
        # For supported versions, we look up the range of supported CCs
        min_cc, max_cc = CTK_SUPPORTED[ctk_version]
        return tuple([cc for cc in COMPUTE_CAPABILITIES
                      if min_cc <= cc <= max_cc])
    except KeyError:
        # For unsupported CUDA toolkit versions, all we can do is assume all
        # non-deprecated versions we are aware of are supported.
        return tuple([cc for cc in COMPUTE_CAPABILITIES
                      if cc >= config.CUDA_DEFAULT_PTX_CC])


def get_supported_ccs():
    try:
        from numba.cuda.cudadrv.runtime import runtime
        cudart_version = runtime.get_version()
    except: # noqa: E722
        # We can't support anything if there's an error getting the runtime
        # version (e.g. if it's not present or there's another issue)
        _supported_cc = ()
        return _supported_cc

    # Ensure the minimum CTK version requirement is met
    min_cudart = min(CTK_SUPPORTED)
    if cudart_version < min_cudart:
        _supported_cc = ()
        ctk_ver = f"{cudart_version[0]}.{cudart_version[1]}"
        unsupported_ver = (f"CUDA Toolkit {ctk_ver} is unsupported by Numba - "
                           f"{min_cudart[0]}.{min_cudart[1]} is the minimum "
                           "required version.")
        warnings.warn(unsupported_ver)
        return _supported_cc

    _supported_cc = ccs_supported_by_ctk(cudart_version)
    return _supported_cc


def find_closest_arch(mycc):
    """
    Given a compute capability, return the closest compute capability supported
    by the CUDA toolkit.

    :param mycc: Compute capability as a tuple ``(MAJOR, MINOR)``
    :return: Closest supported CC as a tuple ``(MAJOR, MINOR)``
    """
    supported_ccs = NVVM().supported_ccs

    if not supported_ccs:
        msg = "No supported GPU compute capabilities found. " \
              "Please check your cudatoolkit version matches your CUDA version."
        raise NvvmSupportError(msg)

    for i, cc in enumerate(supported_ccs):
        if cc == mycc:
            # Matches
            return cc
        elif cc > mycc:
            # Exceeded
            if i == 0:
                # CC lower than supported
                msg = "GPU compute capability %d.%d is not supported" \
                      "(requires >=%d.%d)" % (mycc + cc)
                raise NvvmSupportError(msg)
            else:
                # return the previous CC
                return supported_ccs[i - 1]

    # CC higher than supported
    return supported_ccs[-1]  # Choose the highest


def get_arch_option(major, minor):
    """Matches with the closest architecture option
    """
    if config.FORCE_CUDA_CC:
        arch = config.FORCE_CUDA_CC
    else:
        arch = find_closest_arch((major, minor))
    return 'compute_%d%d' % arch


MISSING_LIBDEVICE_FILE_MSG = '''Missing libdevice file.
Please ensure you have a CUDA Toolkit 11.2 or higher.
For CUDA 12, ``cuda-nvcc`` and ``cuda-nvrtc`` are required:

    $ conda install -c conda-forge cuda-nvcc cuda-nvrtc "cuda-version>=12.0"

For CUDA 11, ``cudatoolkit`` is required:

    $ conda install -c conda-forge cudatoolkit "cuda-version>=11.2,<12.0"
'''


class LibDevice(object):
    _cache_ = None

    def __init__(self):
        if self._cache_ is None:
            if get_libdevice() is None:
                raise RuntimeError(MISSING_LIBDEVICE_FILE_MSG)
            self._cache_ = open_libdevice()

        self.bc = self._cache_

    def get(self):
        return self.bc


cas_nvvm = """
    %cas_success = cmpxchg volatile {Ti}* %iptr, {Ti} %old, {Ti} %new monotonic monotonic
    %cas = extractvalue {{ {Ti}, i1 }} %cas_success, 0
""" # noqa: E501


# Translation of code from CUDA Programming Guide v6.5, section B.12
ir_numba_atomic_binary_template = """
define internal {T} @___numba_atomic_{T}_{FUNC}({T}* %ptr, {T} %val) alwaysinline {{
entry:
    %iptr = bitcast {T}* %ptr to {Ti}*
    %old2 = load volatile {Ti}, {Ti}* %iptr
    br label %attempt

attempt:
    %old = phi {Ti} [ %old2, %entry ], [ %cas, %attempt ]
    %dold = bitcast {Ti} %old to {T}
    %dnew = {OP} {T} %dold, %val
    %new = bitcast {T} %dnew to {Ti}
    {CAS}
    %repeat = icmp ne {Ti} %cas, %old
    br i1 %repeat, label %attempt, label %done

done:
    %result = bitcast {Ti} %old to {T}
    ret {T} %result
}}
""" # noqa: E501

ir_numba_atomic_inc_template = """
define internal {T} @___numba_atomic_{Tu}_inc({T}* %iptr, {T} %val) alwaysinline {{
entry:
    %old2 = load volatile {T}, {T}* %iptr
    br label %attempt

attempt:
    %old = phi {T} [ %old2, %entry ], [ %cas, %attempt ]
    %bndchk = icmp ult {T} %old, %val
    %inc = add {T} %old, 1
    %new = select i1 %bndchk, {T} %inc, {T} 0
    {CAS}
    %repeat = icmp ne {T} %cas, %old
    br i1 %repeat, label %attempt, label %done

done:
    ret {T} %old
}}
""" # noqa: E501

ir_numba_atomic_dec_template = """
define internal {T} @___numba_atomic_{Tu}_dec({T}* %iptr, {T} %val) alwaysinline {{
entry:
    %old2 = load volatile {T}, {T}* %iptr
    br label %attempt

attempt:
    %old = phi {T} [ %old2, %entry ], [ %cas, %attempt ]
    %dec = add {T} %old, -1
    %bndchk = icmp ult {T} %dec, %val
    %new = select i1 %bndchk, {T} %dec, {T} %val
    {CAS}
    %repeat = icmp ne {T} %cas, %old
    br i1 %repeat, label %attempt, label %done

done:
    ret {T} %old
}}
""" # noqa: E501

ir_numba_atomic_minmax_template = """
define internal {T} @___numba_atomic_{T}_{NAN}{FUNC}({T}* %ptr, {T} %val) alwaysinline {{
entry:
    %ptrval = load volatile {T}, {T}* %ptr
    ; Return early when:
    ; - For nanmin / nanmax when val is a NaN
    ; - For min / max when val or ptr is a NaN
    %early_return = fcmp uno {T} %val, %{PTR_OR_VAL}val
    br i1 %early_return, label %done, label %lt_check

lt_check:
    %dold = phi {T} [ %ptrval, %entry ], [ %dcas, %attempt ]
    ; Continue attempts if dold less or greater than val (depending on whether min or max)
    ; or if dold is NaN (for nanmin / nanmax)
    %cmp = fcmp {OP} {T} %dold, %val
    br i1 %cmp, label %attempt, label %done

attempt:
    ; Attempt to swap in the value
    %old = bitcast {T} %dold to {Ti}
    %iptr = bitcast {T}* %ptr to {Ti}*
    %new = bitcast {T} %val to {Ti}
    {CAS}
    %dcas = bitcast {Ti} %cas to {T}
    br label %lt_check

done:
    ret {T} %ptrval
}}
""" # noqa: E501


def ir_cas(Ti):
    return cas_nvvm.format(Ti=Ti)


def ir_numba_atomic_binary(T, Ti, OP, FUNC):
    params = dict(T=T, Ti=Ti, OP=OP, FUNC=FUNC, CAS=ir_cas(Ti))
    return ir_numba_atomic_binary_template.format(**params)


def ir_numba_atomic_minmax(T, Ti, NAN, OP, PTR_OR_VAL, FUNC):
    params = dict(T=T, Ti=Ti, NAN=NAN, OP=OP, PTR_OR_VAL=PTR_OR_VAL,
                  FUNC=FUNC, CAS=ir_cas(Ti))

    return ir_numba_atomic_minmax_template.format(**params)


def ir_numba_atomic_inc(T, Tu):
    return ir_numba_atomic_inc_template.format(T=T, Tu=Tu, CAS=ir_cas(T))


def ir_numba_atomic_dec(T, Tu):
    return ir_numba_atomic_dec_template.format(T=T, Tu=Tu, CAS=ir_cas(T))


def llvm_replace(llvmir):
    replacements = [
        ('declare double @"___numba_atomic_double_add"(double* %".1", double %".2")',     # noqa: E501
         ir_numba_atomic_binary(T='double', Ti='i64', OP='fadd', FUNC='add')),
        ('declare float @"___numba_atomic_float_sub"(float* %".1", float %".2")',         # noqa: E501
         ir_numba_atomic_binary(T='float', Ti='i32', OP='fsub', FUNC='sub')),
        ('declare double @"___numba_atomic_double_sub"(double* %".1", double %".2")',     # noqa: E501
         ir_numba_atomic_binary(T='double', Ti='i64', OP='fsub', FUNC='sub')),
        ('declare i64 @"___numba_atomic_u64_inc"(i64* %".1", i64 %".2")',
         ir_numba_atomic_inc(T='i64', Tu='u64')),
        ('declare i64 @"___numba_atomic_u64_dec"(i64* %".1", i64 %".2")',
         ir_numba_atomic_dec(T='i64', Tu='u64')),
        ('declare float @"___numba_atomic_float_max"(float* %".1", float %".2")',         # noqa: E501
         ir_numba_atomic_minmax(T='float', Ti='i32', NAN='', OP='nnan olt',
                                PTR_OR_VAL='ptr', FUNC='max')),
        ('declare double @"___numba_atomic_double_max"(double* %".1", double %".2")',     # noqa: E501
         ir_numba_atomic_minmax(T='double', Ti='i64', NAN='', OP='nnan olt',
                                PTR_OR_VAL='ptr', FUNC='max')),
        ('declare float @"___numba_atomic_float_min"(float* %".1", float %".2")',         # noqa: E501
         ir_numba_atomic_minmax(T='float', Ti='i32', NAN='', OP='nnan ogt',
                                PTR_OR_VAL='ptr', FUNC='min')),
        ('declare double @"___numba_atomic_double_min"(double* %".1", double %".2")',     # noqa: E501
         ir_numba_atomic_minmax(T='double', Ti='i64', NAN='', OP='nnan ogt',
                                PTR_OR_VAL='ptr', FUNC='min')),
        ('declare float @"___numba_atomic_float_nanmax"(float* %".1", float %".2")',      # noqa: E501
         ir_numba_atomic_minmax(T='float', Ti='i32', NAN='nan', OP='ult',
                                PTR_OR_VAL='', FUNC='max')),
        ('declare double @"___numba_atomic_double_nanmax"(double* %".1", double %".2")',  # noqa: E501
         ir_numba_atomic_minmax(T='double', Ti='i64', NAN='nan', OP='ult',
                                PTR_OR_VAL='', FUNC='max')),
        ('declare float @"___numba_atomic_float_nanmin"(float* %".1", float %".2")',      # noqa: E501
         ir_numba_atomic_minmax(T='float', Ti='i32', NAN='nan', OP='ugt',
                                PTR_OR_VAL='', FUNC='min')),
        ('declare double @"___numba_atomic_double_nanmin"(double* %".1", double %".2")',  # noqa: E501
         ir_numba_atomic_minmax(T='double', Ti='i64', NAN='nan', OP='ugt',
                                PTR_OR_VAL='', FUNC='min')),
        ('immarg', '')
    ]

    for decl, fn in replacements:
        llvmir = llvmir.replace(decl, fn)

    llvmir = llvm140_to_70_ir(llvmir)

    return llvmir


def compile_ir(llvmir, **opts):
    if isinstance(llvmir, str):
        llvmir = [llvmir]

    if opts.pop('fastmath', False):
        opts.update({
            'ftz': True,
            'fma': True,
            'prec_div': False,
            'prec_sqrt': False,
        })

    cu = CompilationUnit()
    libdevice = LibDevice()

    for mod in llvmir:
        mod = llvm_replace(mod)
        cu.add_module(mod.encode('utf8'))
    cu.lazy_add_module(libdevice.get())

    return cu.compile(**opts)


re_attributes_def = re.compile(r"^attributes #\d+ = \{ ([\w\s]+)\ }")


def llvm140_to_70_ir(ir):
    """
    Convert LLVM 14.0 IR for LLVM 7.0.
    """
    buf = []
    for line in ir.splitlines():
        if line.startswith('attributes #'):
            # Remove function attributes unsupported by LLVM 7.0
            m = re_attributes_def.match(line)
            attrs = m.group(1).split()
            attrs = ' '.join(a for a in attrs if a != 'willreturn')
            line = line.replace(m.group(1), attrs)

        buf.append(line)

    return '\n'.join(buf)


def set_cuda_kernel(function):
    """
    Mark a function as a CUDA kernel. Kernels have the following requirements:

    - Metadata that marks them as a kernel.
    - Addition to the @llvm.used list, so that they will not be discarded.
    - The noinline attribute is not permitted, because this causes NVVM to emit
      a warning, which counts as failing IR verification.

    Presently it is assumed that there is one kernel per module, which holds
    for Numba-jitted functions. If this changes in future or this function is
    to be used externally, this function may need modification to add to the
    @llvm.used list rather than creating it.
    """
    module = function.module

    # Add kernel metadata
    mdstr = ir.MetaDataString(module, "kernel")
    mdvalue = ir.Constant(ir.IntType(32), 1)
    md = module.add_metadata((function, mdstr, mdvalue))

    nmd = cgutils.get_or_insert_named_metadata(module, 'nvvm.annotations')
    nmd.add(md)

    # Create the used list
    ptrty = ir.IntType(8).as_pointer()
    usedty = ir.ArrayType(ptrty, 1)

    fnptr = function.bitcast(ptrty)

    llvm_used = ir.GlobalVariable(module, usedty, 'llvm.used')
    llvm_used.linkage = 'appending'
    llvm_used.section = 'llvm.metadata'
    llvm_used.initializer = ir.Constant(usedty, [fnptr])

    # Remove 'noinline' if it is present.
    function.attributes.discard('noinline')


def add_ir_version(mod):
    """Add NVVM IR version to module"""
    # We specify the IR version to match the current NVVM's IR version
    i32 = ir.IntType(32)
    ir_versions = [i32(v) for v in NVVM().get_ir_version()]
    md_ver = mod.add_metadata(ir_versions)
    mod.add_named_metadata('nvvmir.version', md_ver)
