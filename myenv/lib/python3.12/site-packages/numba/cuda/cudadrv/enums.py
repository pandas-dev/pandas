"""
Enum values for CUDA driver. Information about the values
can be found on the official NVIDIA documentation website.
ref: https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__TYPES.html
anchor: #group__CUDA__TYPES
"""


# Error codes

CUDA_SUCCESS = 0
CUDA_ERROR_INVALID_VALUE = 1
CUDA_ERROR_OUT_OF_MEMORY = 2
CUDA_ERROR_NOT_INITIALIZED = 3
CUDA_ERROR_DEINITIALIZED = 4
CUDA_ERROR_PROFILER_DISABLED = 5
CUDA_ERROR_PROFILER_NOT_INITIALIZED = 6
CUDA_ERROR_PROFILER_ALREADY_STARTED = 7
CUDA_ERROR_PROFILER_ALREADY_STOPPED = 8
CUDA_ERROR_STUB_LIBRARY = 34
CUDA_ERROR_DEVICE_UNAVAILABLE = 46
CUDA_ERROR_NO_DEVICE = 100
CUDA_ERROR_INVALID_DEVICE = 101
CUDA_ERROR_DEVICE_NOT_LICENSED = 102
CUDA_ERROR_INVALID_IMAGE = 200
CUDA_ERROR_INVALID_CONTEXT = 201
CUDA_ERROR_CONTEXT_ALREADY_CURRENT = 202
CUDA_ERROR_MAP_FAILED = 205
CUDA_ERROR_UNMAP_FAILED = 206
CUDA_ERROR_ARRAY_IS_MAPPED = 207
CUDA_ERROR_ALREADY_MAPPED = 208
CUDA_ERROR_NO_BINARY_FOR_GPU = 209
CUDA_ERROR_ALREADY_ACQUIRED = 210
CUDA_ERROR_NOT_MAPPED = 211
CUDA_ERROR_NOT_MAPPED_AS_ARRAY = 212
CUDA_ERROR_NOT_MAPPED_AS_POINTER = 213
CUDA_ERROR_ECC_UNCORRECTABLE = 214
CUDA_ERROR_UNSUPPORTED_LIMIT = 215
CUDA_ERROR_CONTEXT_ALREADY_IN_USE = 216
CUDA_ERROR_PEER_ACCESS_UNSUPPORTED = 217
CUDA_ERROR_INVALID_PTX = 218
CUDA_ERROR_INVALID_GRAPHICS_CONTEXT = 219
CUDA_ERROR_NVLINK_UNCORRECTABLE = 220
CUDA_ERROR_JIT_COMPILER_NOT_FOUND = 221
CUDA_ERROR_UNSUPPORTED_PTX_VERSION = 222
CUDA_ERROR_JIT_COMPILATION_DISABLED = 223
CUDA_ERROR_UNSUPPORTED_EXEC_AFFINITY = 224
CUDA_ERROR_UNSUPPORTED_DEVSIDE_SYNC = 225
CUDA_ERROR_INVALID_SOURCE = 300
CUDA_ERROR_FILE_NOT_FOUND = 301
CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND = 302
CUDA_ERROR_SHARED_OBJECT_INIT_FAILED = 303
CUDA_ERROR_OPERATING_SYSTEM = 304
CUDA_ERROR_INVALID_HANDLE = 400
CUDA_ERROR_ILLEGAL_STATE = 401
CUDA_ERROR_NOT_FOUND = 500
CUDA_ERROR_NOT_READY = 600
CUDA_ERROR_LAUNCH_FAILED = 700
CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES = 701
CUDA_ERROR_LAUNCH_TIMEOUT = 702
CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING = 703
CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED = 704
CUDA_ERROR_PEER_ACCESS_NOT_ENABLED = 705
CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE = 708
CUDA_ERROR_CONTEXT_IS_DESTROYED = 709
CUDA_ERROR_ASSERT = 710
CUDA_ERROR_TOO_MANY_PEERS = 711
CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED = 712
CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED = 713
CUDA_ERROR_HARDWARE_STACK_ERROR = 714
CUDA_ERROR_ILLEGAL_INSTRUCTION = 715
CUDA_ERROR_MISALIGNED_ADDRESS = 716
CUDA_ERROR_INVALID_ADDRESS_SPACE = 717
CUDA_ERROR_INVALID_PC = 718
CUDA_ERROR_LAUNCH_FAILED = 719
CUDA_ERROR_COOPERATIVE_LAUNCH_TOO_LARGE = 720
CUDA_ERROR_NOT_PERMITTED = 800
CUDA_ERROR_NOT_SUPPORTED = 801
CUDA_ERROR_SYSTEM_NOT_READY = 802
CUDA_ERROR_SYSTEM_DRIVER_MISMATCH = 803
CUDA_ERROR_COMPAT_NOT_SUPPORTED_ON_DEVICE = 804
CUDA_ERROR_MPS_CONNECTION_FAILED = 805
CUDA_ERROR_MPS_RPC_FAILURE = 806
CUDA_ERROR_MPS_SERVER_NOT_READY = 807
CUDA_ERROR_MPS_MAX_CLIENTS_REACHED = 808
CUDA_ERROR_MPS_MAX_CONNECTIONS_REACHED = 809
CUDA_ERROR_MPS_CLIENT_TERMINATED = 810
CUDA_ERROR_CDP_NOT_SUPPORTED = 811
CUDA_ERROR_CDP_VERSION_MISMATCH = 812
CUDA_ERROR_STREAM_CAPTURE_UNSUPPORTED = 900
CUDA_ERROR_STREAM_CAPTURE_INVALIDATED = 901
CUDA_ERROR_STREAM_CAPTURE_MERGE = 902
CUDA_ERROR_STREAM_CAPTURE_UNMATCHED = 903
CUDA_ERROR_STREAM_CAPTURE_UNJOINED = 904
CUDA_ERROR_STREAM_CAPTURE_ISOLATION = 905
CUDA_ERROR_STREAM_CAPTURE_IMPLICIT = 906
CUDA_ERROR_CAPTURED_EVENT = 907
CUDA_ERROR_STREAM_CAPTURE_WRONG_THREAD = 908
CUDA_ERROR_TIMEOUT = 909
CUDA_ERROR_GRAPH_EXEC_UPDATE_FAILURE = 910
CUDA_ERROR_EXTERNAL_DEVICE = 911
CUDA_ERROR_INVALID_CLUSTER_SIZE = 912
CUDA_ERROR_UNKNOWN = 999


# Function cache configurations

# no preference for shared memory or L1 (default)
CU_FUNC_CACHE_PREFER_NONE = 0x00
# prefer larger shared memory and smaller L1 cache
CU_FUNC_CACHE_PREFER_SHARED = 0x01
# prefer larger L1 cache and smaller shared memory
CU_FUNC_CACHE_PREFER_L1 = 0x02
# prefer equal sized L1 cache and shared memory
CU_FUNC_CACHE_PREFER_EQUAL = 0x03


# Context creation flags

# Automatic scheduling
CU_CTX_SCHED_AUTO = 0x00
# Set spin as default scheduling
CU_CTX_SCHED_SPIN = 0x01
# Set yield as default scheduling
CU_CTX_SCHED_YIELD = 0x02
# Set blocking synchronization as default scheduling
CU_CTX_SCHED_BLOCKING_SYNC = 0x04

CU_CTX_SCHED_MASK = 0x07
# Support mapped pinned allocations
#   This flag was deprecated as of CUDA 11.0 and it no longer has effect.
#   All contexts as of CUDA 3.2 behave as though the flag is enabled.
CU_CTX_MAP_HOST = 0x08
# Keep local memory allocation after launch
CU_CTX_LMEM_RESIZE_TO_MAX = 0x10
# Trigger coredumps from exceptions in this context
CU_CTX_COREDUMP_ENABLE = 0x20
# Enable user pipe to trigger coredumps in this context
CU_CTX_USER_COREDUMP_ENABLE = 0x40
# Force synchronous blocking on cudaMemcpy/cudaMemset
CU_CTX_SYNC_MEMOPS = 0x80

CU_CTX_FLAGS_MASK = 0xff


# DEFINES

# If set, host memory is portable between CUDA contexts.
# Flag for cuMemHostAlloc()
CU_MEMHOSTALLOC_PORTABLE = 0x01

# If set, host memory is mapped into CUDA address space and
# cuMemHostGetDevicePointer() may be called on the host pointer.
# Flag for cuMemHostAlloc()
CU_MEMHOSTALLOC_DEVICEMAP = 0x02

# If set, host memory is allocated as write-combined - fast to write,
# faster to DMA, slow to read except via SSE4 streaming load instruction
# (MOVNTDQA).
# Flag for cuMemHostAlloc()
CU_MEMHOSTALLOC_WRITECOMBINED = 0x04


# If set, host memory is portable between CUDA contexts.
# Flag for cuMemHostRegister()
CU_MEMHOSTREGISTER_PORTABLE = 0x01

# If set, host memory is mapped into CUDA address space and
# cuMemHostGetDevicePointer() may be called on the host pointer.
# Flag for cuMemHostRegister()
CU_MEMHOSTREGISTER_DEVICEMAP = 0x02

# If set, the passed memory pointer is treated as pointing to some
# memory-mapped I/O space, e.g. belonging to a third-party PCIe device.
# On Windows the flag is a no-op. On Linux that memory is marked
# as non cache-coherent for the GPU and is expected
# to be physically contiguous. It may return CUDA_ERROR_NOT_PERMITTED
# if run as an unprivileged user, CUDA_ERROR_NOT_SUPPORTED on older
# Linux kernel versions. On all other platforms, it is not supported
# and CUDA_ERROR_NOT_SUPPORTED is returned.
# Flag for cuMemHostRegister()
CU_MEMHOSTREGISTER_IOMEMORY = 0x04

# If set, the passed memory pointer is treated as pointing to memory
# that is considered read-only by the device. On platforms without
# CU_DEVICE_ATTRIBUTE_PAGEABLE_MEMORY_ACCESS_USES_HOST_PAGE_TABLES,
# this flag is required in order to register memory mapped
# to the CPU as read-only. Support for the use of this flag can be
# queried from the device attribute
# CU_DEVICE_ATTRIBUTE_READ_ONLY_HOST_REGISTER_SUPPORTED.
# Using this flag with a current context associated with a device
# that does not have this attribute set will cause cuMemHostRegister
# to error with CUDA_ERROR_NOT_SUPPORTED.
CU_MEMHOSTREGISTER_READ_ONLY = 0x08


# CUDA Mem Attach Flags

# If set, managed memory is accessible from all streams on all devices.
CU_MEM_ATTACH_GLOBAL = 0x01

# If set on a platform where the device attribute
# cudaDevAttrConcurrentManagedAccess is zero, then managed memory is
# only accessible on the host (unless explicitly attached to a stream
# with cudaStreamAttachMemAsync, in which case it can be used in kernels
# launched on that stream).
CU_MEM_ATTACH_HOST = 0x02

# If set on a platform where the device attribute
# cudaDevAttrConcurrentManagedAccess is zero, then managed memory accesses
# on the associated device must only be from a single stream.
CU_MEM_ATTACH_SINGLE = 0x04


# Event creation flags

# Default event flag
CU_EVENT_DEFAULT = 0x0
# Event uses blocking synchronization
CU_EVENT_BLOCKING_SYNC = 0x1
# Event will not record timing data
CU_EVENT_DISABLE_TIMING = 0x2
# Event is suitable for interprocess use. CU_EVENT_DISABLE_TIMING must be set
CU_EVENT_INTERPROCESS = 0x4


# Pointer information

# The CUcontext on which a pointer was allocated or registered
CU_POINTER_ATTRIBUTE_CONTEXT = 1
# The CUmemorytype describing the physical location of a pointer
CU_POINTER_ATTRIBUTE_MEMORY_TYPE = 2
# The address at which a pointer's memory may be accessed on the device
CU_POINTER_ATTRIBUTE_DEVICE_POINTER = 3
# The address at which a pointer's memory may be accessed on the host
CU_POINTER_ATTRIBUTE_HOST_POINTER = 4
# A pair of tokens for use with the nv-p2p.h Linux kernel interface
CU_POINTER_ATTRIBUTE_P2P_TOKENS = 5
# Synchronize every synchronous memory operation initiated on this region
CU_POINTER_ATTRIBUTE_SYNC_MEMOPS = 6
# A process-wide unique ID for an allocated memory region
CU_POINTER_ATTRIBUTE_BUFFER_ID = 7
# Indicates if the pointer points to managed memory
CU_POINTER_ATTRIBUTE_IS_MANAGED = 8
# A device ordinal of a device on which a pointer was allocated or registered
CU_POINTER_ATTRIBUTE_DEVICE_ORDINAL = 9
# 1 if this pointer maps to an allocation
# that is suitable for cudaIpcGetMemHandle, 0 otherwise
CU_POINTER_ATTRIBUTE_IS_LEGACY_CUDA_IPC_CAPABLE = 10
# Starting address for this requested pointer
CU_POINTER_ATTRIBUTE_RANGE_START_ADDR = 11
# Size of the address range for this requested pointer
CU_POINTER_ATTRIBUTE_RANGE_SIZE = 12
# 1 if this pointer is in a valid address range
# that is mapped to a backing allocation, 0 otherwise
CU_POINTER_ATTRIBUTE_MAPPED = 13
# Bitmask of allowed CUmemAllocationHandleType for this allocation
CU_POINTER_ATTRIBUTE_ALLOWED_HANDLE_TYPES = 14
# 1 if the memory this pointer is referencing
# can be used with the GPUDirect RDMA API
CU_POINTER_ATTRIBUTE_IS_GPU_DIRECT_RDMA_CAPABLE = 15
# Returns the access flags the device associated
# with the current context has on the corresponding
# memory referenced by the pointer given
CU_POINTER_ATTRIBUTE_ACCESS_FLAGS = 16
# Returns the mempool handle for the allocation
# if it was allocated from a mempool. Otherwise returns NULL
CU_POINTER_ATTRIBUTE_MEMPOOL_HANDLE = 17
# Size of the actual underlying mapping that the pointer belongs to
CU_POINTER_ATTRIBUTE_MAPPING_SIZE = 18
# The start address of the mapping that the pointer belongs to
CU_POINTER_ATTRIBUTE_MAPPING_BASE_ADDR = 19
# A process-wide unique id corresponding to the
# physical allocation the pointer belongs to
CU_POINTER_ATTRIBUTE_MEMORY_BLOCK_ID = 20


# Memory types

# Host memory
CU_MEMORYTYPE_HOST = 0x01
# Device memory
CU_MEMORYTYPE_DEVICE = 0x02
# Array memory
CU_MEMORYTYPE_ARRAY = 0x03
# Unified device or host memory
CU_MEMORYTYPE_UNIFIED = 0x04


# Device code formats

# Compiled device-class-specific device code
# Applicable options: none
CU_JIT_INPUT_CUBIN = 0

# PTX source code
# Applicable options: PTX compiler options
CU_JIT_INPUT_PTX = 1

# Bundle of multiple cubins and/or PTX of some device code
# Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
CU_JIT_INPUT_FATBINARY = 2

# Host object with embedded device code
# Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
CU_JIT_INPUT_OBJECT = 3

# Archive of host objects with embedded device code
# Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
CU_JIT_INPUT_LIBRARY = 4

CU_JIT_NUM_INPUT_TYPES = 6


# Online compiler and linker options

# Max number of registers that a thread may use.
# Option type: unsigned int
# Applies to: compiler only
CU_JIT_MAX_REGISTERS = 0

# IN: Specifies minimum number of threads per block to target compilation
# for
# OUT: Returns the number of threads the compiler actually targeted.
# This restricts the resource utilization fo the compiler (e.g. max
# registers) such that a block with the given number of threads should be
# able to launch based on register limitations. Note, this option does not
# currently take into account any other resource limitations, such as
# shared memory utilization.
# Cannot be combined with ::CU_JIT_TARGET.
# Option type: unsigned int
# Applies to: compiler only
CU_JIT_THREADS_PER_BLOCK = 1

# Overwrites the option value with the total wall clock time, in
# milliseconds, spent in the compiler and linker
# Option type: float
# Applies to: compiler and linker
CU_JIT_WALL_TIME = 2

# Pointer to a buffer in which to print any log messages
# that are informational in nature (the buffer size is specified via
# option ::CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES)
# Option type: char *
# Applies to: compiler and linker
CU_JIT_INFO_LOG_BUFFER = 3

# IN: Log buffer size in bytes.  Log messages will be capped at this size
# (including null terminator)
# OUT: Amount of log buffer filled with messages
# Option type: unsigned int
# Applies to: compiler and linker
CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES = 4

# Pointer to a buffer in which to print any log messages that
# reflect errors (the buffer size is specified via option
# ::CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES)
# Option type: char *
# Applies to: compiler and linker
CU_JIT_ERROR_LOG_BUFFER = 5

# IN: Log buffer size in bytes.  Log messages will be capped at this size
# (including null terminator)
# OUT: Amount of log buffer filled with messages
# Option type: unsigned int
# Applies to: compiler and linker
CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES = 6

# Level of optimizations to apply to generated code (0 - 4), with 4
# being the default and highest level of optimizations.
# Option type: unsigned int
# Applies to: compiler only
CU_JIT_OPTIMIZATION_LEVEL = 7

# No option value required. Determines the target based on the current
# attached context (default)
# Option type: No option value needed
# Applies to: compiler and linker
CU_JIT_TARGET_FROM_CUCONTEXT = 8

# Target is chosen based on supplied ::CUjit_target.  Cannot be
# combined with ::CU_JIT_THREADS_PER_BLOCK.
# Option type: unsigned int for enumerated type ::CUjit_target
# Applies to: compiler and linker
CU_JIT_TARGET = 9

# Specifies choice of fallback strategy if matching cubin is not found.
# Choice is based on supplied ::CUjit_fallback.
# Option type: unsigned int for enumerated type ::CUjit_fallback
# Applies to: compiler only
CU_JIT_FALLBACK_STRATEGY = 10

# Specifies whether to create debug information in output (-g)
# (0: false, default)
# Option type: int
# Applies to: compiler and linker
CU_JIT_GENERATE_DEBUG_INFO = 11

# Generate verbose log messages (0: false, default)
# Option type: int
# Applies to: compiler and linker
CU_JIT_LOG_VERBOSE = 12

# Generate line number information (-lineinfo) (0: false, default)
# Option type: int
# Applies to: compiler only
CU_JIT_GENERATE_LINE_INFO = 13

# Specifies whether to enable caching explicitly (-dlcm)
# Choice is based on supplied ::CUjit_cacheMode_enum.
# Option type: unsigned int for enumerated type ::CUjit_cacheMode_enum
# Applies to: compiler only
CU_JIT_CACHE_MODE = 14


# CUfunction_attribute

# The maximum number of threads per block, beyond which a launch of the
# function would fail. This number depends on both the function and the
# device on which the function is currently loaded.
CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 0

# The size in bytes of statically-allocated shared memory required by
# this function. This does not include dynamically-allocated shared
# memory requested by the user at runtime.
CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES = 1

# The size in bytes of user-allocated constant memory required by this
# function.
CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES = 2

# The size in bytes of local memory used by each thread of this function.
CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES = 3

# The number of registers used by each thread of this function.
CU_FUNC_ATTRIBUTE_NUM_REGS = 4

# The PTX virtual architecture version for which the function was
# compiled. This value is the major PTX version * 10 + the minor PTX
# version, so a PTX version 1.3 function would return the value 13.
# Note that this may return the undefined value of 0 for cubins
# compiled prior to CUDA 3.0.
CU_FUNC_ATTRIBUTE_PTX_VERSION = 5

# The binary architecture version for which the function was compiled.
# This value is the major binary version * 10 + the minor binary version,
# so a binary version 1.3 function would return the value 13. Note that
# this will return a value of 10 for legacy cubins that do not have a
# properly-encoded binary architecture version.
CU_FUNC_ATTRIBUTE_BINARY_VERSION = 6

# The attribute to indicate whether the function has been compiled
# with user specified option "-Xptxas --dlcm=ca" set
CU_FUNC_ATTRIBUTE_CACHE_MODE_CA = 7

# The maximum size in bytes of dynamically-allocated shared memory
# that can be used by this function. If the user-specified
# dynamic shared memory size is larger than this value,
# the launch will fail. See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES = 8

# On devices where the L1 cache and shared memory use the same
# hardware resources, this sets the shared memory carveout preference,
# in percent of the total shared memory. Refer to
# CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR.
# This is only a hint, and the driver can choose a different ratio
# if required to execute the function.
# See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_PREFERRED_SHARED_MEMORY_CARVEOUT = 9

# If this attribute is set, the kernel must launch with a valid cluster
# size specified. See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_CLUSTER_SIZE_MUST_BE_SET = 10

# The required cluster width in blocks. The values must either all be 0
# or all be positive. The validity of the cluster dimensions
# is otherwise checked at launch time. If the value is set during
# compile time, it cannot be set at runtime.
# Setting it at runtime will return CUDA_ERROR_NOT_PERMITTED.
# See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_REQUIRED_CLUSTER_WIDTH = 11

# The required cluster height in blocks. The values must either all be 0
# or all be positive. The validity of the cluster dimensions
# is otherwise checked at launch time.If the value is set during
# compile time, it cannot be set at runtime.
# Setting it at runtime should return CUDA_ERROR_NOT_PERMITTED.
# See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_REQUIRED_CLUSTER_HEIGHT = 12

# The required cluster depth in blocks. The values must either all be 0
# or all be positive. The validity of the cluster dimensions
# is otherwise checked at launch time.If the value is set during
# compile time, it cannot be set at runtime.
# Setting it at runtime should return CUDA_ERROR_NOT_PERMITTED.
# See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_REQUIRED_CLUSTER_DEPTH = 13

# Whether the function can be launched with non-portable cluster size.
# 1 is allowed, 0 is disallowed. A non-portable cluster size may only
# function on the specific SKUs the program is tested on.
# The launch might fail if the program is run on a different hardware platform.
# For more details refer to link :
# https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__TYPES.html#group__CUDA__TYPES
CU_FUNC_ATTRIBUTE_NON_PORTABLE_CLUSTER_SIZE_ALLOWED = 14

# The block scheduling policy of a function.
# The value type is CUclusterSchedulingPolicy / cudaClusterSchedulingPolicy.
# See cuFuncSetAttribute, cuKernelSetAttribute
CU_FUNC_ATTRIBUTE_CLUSTER_SCHEDULING_POLICY_PREFERENCE = 15


# Device attributes

CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1
CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X = 2
CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y = 3
CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z = 4
CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X = 5
CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y = 6
CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z = 7
CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8
CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY = 9
CU_DEVICE_ATTRIBUTE_WARP_SIZE = 10
CU_DEVICE_ATTRIBUTE_MAX_PITCH = 11
CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK = 12
CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13
CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT = 14
CU_DEVICE_ATTRIBUTE_GPU_OVERLAP = 15
CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16
CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT = 17
CU_DEVICE_ATTRIBUTE_INTEGRATED = 18
CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY = 19
CU_DEVICE_ATTRIBUTE_COMPUTE_MODE = 20
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_1D_WIDTH = 21
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_WIDTH = 22
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_HEIGHT = 23
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_WIDTH = 24
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_HEIGHT = 25
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_DEPTH = 26
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LAYERED_WIDTH = 27
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LAYERED_HEIGHT = 28
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LAYERED_LAYERS = 29
CU_DEVICE_ATTRIBUTE_SURFACE_ALIGNMENT = 30
CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS = 31
CU_DEVICE_ATTRIBUTE_ECC_ENABLED = 32
CU_DEVICE_ATTRIBUTE_PCI_BUS_ID = 33
CU_DEVICE_ATTRIBUTE_PCI_DEVICE_ID = 34
CU_DEVICE_ATTRIBUTE_TCC_DRIVER = 35
CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE = 36
CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH = 37
CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE = 38
CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTI_PROCESSOR = 39
CU_DEVICE_ATTRIBUTE_ASYNC_ENGINE_COUNT = 40
CU_DEVICE_ATTRIBUTE_UNIFIED_ADDRESSING = 41
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_1D_LAYERED_WIDTH = 42
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_1D_LAYERED_LAYERS = 43
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_GATHER_WIDTH = 45
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_GATHER_HEIGHT = 46
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_WIDTH_ALT = 47
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_HEIGHT_ALT = 48
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_3D_DEPTH_ALT = 49
CU_DEVICE_ATTRIBUTE_PCI_DOMAIN_ID = 50
CU_DEVICE_ATTRIBUTE_TEXTURE_PITCH_ALIGNMENT = 51
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_CUBEMAP_WIDTH = 52
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_CUBEMAP_LAYERED_WIDTH = 53
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_CUBEMAP_LAYERED_LAYERS = 54
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_1D_WIDTH = 55
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_2D_WIDTH = 56
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_2D_HEIGHT = 57
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_3D_WIDTH = 58
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_3D_HEIGHT = 59
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_3D_DEPTH = 60
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_1D_LAYERED_WIDTH = 61
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_1D_LAYERED_LAYERS = 62
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_2D_LAYERED_WIDTH = 63
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_2D_LAYERED_HEIGHT = 64
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_2D_LAYERED_LAYERS = 65
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_CUBEMAP_WIDTH = 66
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_CUBEMAP_LAYERED_WIDTH = 67
CU_DEVICE_ATTRIBUTE_MAX_SURFACE_CUBEMAP_LAYERED_LAYERS = 68
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_1D_LINEAR_WIDTH = 69
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LINEAR_WIDTH = 70
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LINEAR_HEIGHT = 71
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_LINEAR_PITCH = 72
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_2D_MIPMAPPED_WIDTH = 73
CU_DEVICE_ATTRIBUTE_MAX_MAX_TEXTURE_2D_MIPMAPPED_HEIGHT = 74
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76
CU_DEVICE_ATTRIBUTE_MAX_TEXTURE_1D_MIPMAPPED_WIDTH = 77
CU_DEVICE_ATTRIBUTE_STREAM_PRIORITIES_SUPPORTED = 78
CU_DEVICE_ATTRIBUTE_GLOBAL_L1_CACHE_SUPPORTED = 79
CU_DEVICE_ATTRIBUTE_LOCAL_L1_CACHE_SUPPORTED = 80
CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR = 81
CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_MULTIPROCESSOR = 82
CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY = 83
CU_DEVICE_ATTRIBUTE_IS_MULTI_GPU_BOARD = 84
CU_DEVICE_ATTRIBUTE_MULTI_GPU_BOARD_GROUP_ID = 85
CU_DEVICE_ATTRIBUTE_HOST_NATIVE_ATOMIC_SUPPORTED = 86
CU_DEVICE_ATTRIBUTE_SINGLE_TO_DOUBLE_PRECISION_PERF_RATIO = 87
CU_DEVICE_ATTRIBUTE_PAGEABLE_MEMORY_ACCESS = 88
CU_DEVICE_ATTRIBUTE_CONCURRENT_MANAGED_ACCESS = 89
CU_DEVICE_ATTRIBUTE_COMPUTE_PREEMPTION_SUPPORTED = 90
CU_DEVICE_ATTRIBUTE_CAN_USE_HOST_POINTER_FOR_REGISTERED_MEM = 91
CU_DEVICE_ATTRIBUTE_COOPERATIVE_LAUNCH = 95
CU_DEVICE_ATTRIBUTE_COOPERATIVE_MULTI_DEVICE_LAUNCH = 96
CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN = 97
