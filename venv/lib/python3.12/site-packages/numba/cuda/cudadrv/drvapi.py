from ctypes import (c_byte, c_char_p, c_float, c_int, c_size_t, c_uint,
                    c_uint8, c_void_p, py_object, CFUNCTYPE, POINTER)

from numba.cuda.cudadrv import _extras

cu_device = c_int
cu_device_attribute = c_int     # enum
cu_context = c_void_p           # an opaque handle
cu_module = c_void_p            # an opaque handle
cu_jit_option = c_int           # enum
cu_jit_input_type = c_int       # enum
cu_function = c_void_p          # an opaque handle
cu_device_ptr = c_size_t        # defined as unsigned long long
cu_stream = c_void_p            # an opaque handle
cu_event = c_void_p
cu_link_state = c_void_p
cu_function_attribute = c_int
cu_ipc_mem_handle = (c_byte * _extras.CUDA_IPC_HANDLE_SIZE)   # 64 bytes wide
cu_uuid = (c_byte * 16)         # Device UUID

cu_stream_callback_pyobj = CFUNCTYPE(None, cu_stream, c_int, py_object)

cu_occupancy_b2d_size = CFUNCTYPE(c_size_t, c_int)

# See https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__TYPES.html
CU_STREAM_DEFAULT = 0
CU_STREAM_LEGACY = 1
CU_STREAM_PER_THREAD = 2

API_PROTOTYPES = {
    # CUresult cuInit(unsigned int Flags);
    'cuInit' : (c_int, c_uint),

    # CUresult cuDriverGetVersion (int* driverVersion )
    'cuDriverGetVersion': (c_int, POINTER(c_int)),

    # CUresult cuDeviceGetCount(int *count);
    'cuDeviceGetCount': (c_int, POINTER(c_int)),

    # CUresult cuDeviceGet(CUdevice *device, int ordinal);
    'cuDeviceGet': (c_int, POINTER(cu_device), c_int),

    # CUresult cuDeviceGetName ( char* name, int  len, CUdevice dev )
    'cuDeviceGetName': (c_int, c_char_p, c_int, cu_device),

    # CUresult cuDeviceGetAttribute(int *pi, CUdevice_attribute attrib,
    #                               CUdevice dev);
    'cuDeviceGetAttribute': (c_int, POINTER(c_int), cu_device_attribute,
                             cu_device),

    # CUresult cuDeviceComputeCapability(int *major, int *minor,
    #                                    CUdevice dev);
    'cuDeviceComputeCapability': (c_int, POINTER(c_int), POINTER(c_int),
                                  cu_device),

    # CUresult cuDevicePrimaryCtxGetState(
    #              CUdevice dev,
    #              unsigned int* flags,
    #              int* active)
    'cuDevicePrimaryCtxGetState': (c_int,
                                   cu_device, POINTER(c_uint), POINTER(c_int)),

    # CUresult cuDevicePrimaryCtxRelease ( CUdevice dev )
    'cuDevicePrimaryCtxRelease': (c_int, cu_device),

    # CUresult cuDevicePrimaryCtxReset ( CUdevice dev )
    'cuDevicePrimaryCtxReset': (c_int, cu_device),

    # CUresult cuDevicePrimaryCtxRetain ( CUcontext* pctx, CUdevice dev )
    'cuDevicePrimaryCtxRetain': (c_int, POINTER(cu_context), cu_device),

    # CUresult cuDevicePrimaryCtxSetFlags ( CUdevice dev, unsigned int  flags )
    'cuDevicePrimaryCtxSetFlags': (c_int, cu_device, c_uint),

    # CUresult cuCtxCreate(CUcontext *pctx, unsigned int flags,
    #                      CUdevice dev);
    'cuCtxCreate': (c_int, POINTER(cu_context), c_uint, cu_device),

    # CUresult cuCtxGetDevice (	CUdevice * 	device	 )
    'cuCtxGetDevice': (c_int, POINTER(cu_device)),

    # CUresult cuCtxGetCurrent (CUcontext *pctx);
    'cuCtxGetCurrent': (c_int, POINTER(cu_context)),

    # CUresult cuCtxPushCurrent (CUcontext pctx);
    'cuCtxPushCurrent': (c_int, cu_context),

    # CUresult cuCtxPopCurrent (CUcontext *pctx);
    'cuCtxPopCurrent': (c_int, POINTER(cu_context)),

    # CUresult cuCtxDestroy(CUcontext pctx);
    'cuCtxDestroy': (c_int, cu_context),

    # CUresult cuModuleLoadDataEx(CUmodule *module, const void *image,
    #                             unsigned int numOptions,
    #                             CUjit_option *options,
    #                             void **optionValues);
    'cuModuleLoadDataEx': (c_int, cu_module, c_void_p, c_uint,
                           POINTER(cu_jit_option), POINTER(c_void_p)),

    # CUresult cuModuleUnload(CUmodule hmod);
    'cuModuleUnload': (c_int, cu_module),

    # CUresult cuModuleGetFunction(CUfunction *hfunc, CUmodule hmod,
    #                              const char *name);
    'cuModuleGetFunction': (c_int, cu_function, cu_module, c_char_p),

    # CUresult cuModuleGetGlobal ( CUdeviceptr* dptr, size_t* bytes, CUmodule
    #                              hmod, const char* name )
    'cuModuleGetGlobal': (c_int, POINTER(cu_device_ptr), POINTER(c_size_t),
                          cu_module, c_char_p),

    # CUresult CUDAAPI cuFuncSetCacheConfig(CUfunction hfunc,
    #                                       CUfunc_cache config);
    'cuFuncSetCacheConfig': (c_int, cu_function, c_uint),

    # CUresult cuMemAlloc(CUdeviceptr *dptr, size_t bytesize);
    'cuMemAlloc': (c_int, POINTER(cu_device_ptr), c_size_t),

    # CUresult cuMemAllocManaged(CUdeviceptr *dptr, size_t bytesize,
    #                            unsigned int flags);
    'cuMemAllocManaged': (c_int, c_void_p, c_size_t, c_uint),

    # CUresult cuMemsetD8(CUdeviceptr dstDevice, unsigned char uc, size_t N)
    'cuMemsetD8': (c_int, cu_device_ptr, c_uint8, c_size_t),

    # CUresult cuMemsetD8Async(CUdeviceptr dstDevice, unsigned char uc,
    #                          size_t N, CUstream hStream);
    'cuMemsetD8Async': (c_int,
                        cu_device_ptr, c_uint8, c_size_t, cu_stream),

    # CUresult cuMemcpyHtoD(CUdeviceptr dstDevice, const void *srcHost,
    #                       size_t ByteCount);
    'cuMemcpyHtoD': (c_int, cu_device_ptr, c_void_p, c_size_t),

    # CUresult cuMemcpyHtoDAsync(CUdeviceptr dstDevice, const void *srcHost,
    #                            size_t ByteCount, CUstream hStream);
    'cuMemcpyHtoDAsync': (c_int, cu_device_ptr, c_void_p, c_size_t,
                          cu_stream),

    # CUresult cuMemcpyDtoD(CUdeviceptr dstDevice, const void *srcDevice,
    #                       size_t ByteCount);
    'cuMemcpyDtoD': (c_int, cu_device_ptr, cu_device_ptr, c_size_t),

    # CUresult cuMemcpyDtoDAsync(CUdeviceptr dstDevice, const void *srcDevice,
    #                            size_t ByteCount, CUstream hStream);
    'cuMemcpyDtoDAsync': (c_int, cu_device_ptr, cu_device_ptr, c_size_t,
                          cu_stream),


    # CUresult cuMemcpyDtoH(void *dstHost, CUdeviceptr srcDevice,
    #                       size_t ByteCount);
    'cuMemcpyDtoH': (c_int, c_void_p, cu_device_ptr, c_size_t),

    # CUresult cuMemcpyDtoHAsync(void *dstHost, CUdeviceptr srcDevice,
    #                            size_t ByteCount, CUstream hStream);
    'cuMemcpyDtoHAsync': (c_int, c_void_p, cu_device_ptr, c_size_t,
                          cu_stream),

    # CUresult cuMemFree(CUdeviceptr dptr);
    'cuMemFree': (c_int, cu_device_ptr),

    # CUresult cuStreamCreate(CUstream *phStream, unsigned int Flags);
    'cuStreamCreate': (c_int, POINTER(cu_stream), c_uint),

    # CUresult cuStreamDestroy(CUstream hStream);
    'cuStreamDestroy': (c_int, cu_stream),

    # CUresult cuStreamSynchronize(CUstream hStream);
    'cuStreamSynchronize': (c_int, cu_stream),

    # CUresult cuStreamAddCallback(
    #              CUstream hStream,
    #              CUstreamCallback callback,
    #              void* userData,
    #              unsigned int flags)
    'cuStreamAddCallback': (c_int, cu_stream, cu_stream_callback_pyobj,
                            py_object, c_uint),

    # CUresult cuLaunchKernel(CUfunction f, unsigned int gridDimX,
    #                        unsigned int gridDimY,
    #                        unsigned int gridDimZ,
    #                        unsigned int blockDimX,
    #                        unsigned int blockDimY,
    #                        unsigned int blockDimZ,
    #                        unsigned int sharedMemBytes,
    #                        CUstream hStream, void **kernelParams,
    #                        void ** extra)
    'cuLaunchKernel': (c_int, cu_function, c_uint, c_uint, c_uint,
                       c_uint, c_uint, c_uint, c_uint, cu_stream,
                       POINTER(c_void_p), POINTER(c_void_p)),

    # CUresult cuLaunchCooperativeKernel(CUfunction f, unsigned int gridDimX,
    #                                   unsigned int gridDimY,
    #                                   unsigned int gridDimZ,
    #                                   unsigned int blockDimX,
    #                                   unsigned int blockDimY,
    #                                   unsigned int blockDimZ,
    #                                   unsigned int sharedMemBytes,
    #                                   CUstream hStream, void **kernelParams)
    'cuLaunchCooperativeKernel': (c_int, cu_function, c_uint, c_uint, c_uint,
                                  c_uint, c_uint, c_uint, c_uint, cu_stream,
                                  POINTER(c_void_p)),

    #  CUresult cuMemHostAlloc (	void ** 	pp,
    #                               size_t 	bytesize,
    #                               unsigned int 	Flags
    #                           )
    'cuMemHostAlloc': (c_int, c_void_p, c_size_t, c_uint),

    #  CUresult cuMemFreeHost (	void * 	p	 )
    'cuMemFreeHost': (c_int, c_void_p),

    # CUresult cuMemHostRegister(void * 	p,
    #                            size_t 	bytesize,
    #                            unsigned int 	Flags)
    'cuMemHostRegister': (c_int, c_void_p, c_size_t, c_uint),

    # CUresult cuMemHostUnregister(void * 	p)
    'cuMemHostUnregister': (c_int, c_void_p),

    # CUresult cuMemHostGetDevicePointer(CUdeviceptr * pdptr,
    #                                    void *        p,
    #                                    unsigned int  Flags)
    'cuMemHostGetDevicePointer': (c_int, POINTER(cu_device_ptr),
                                  c_void_p, c_uint),

    # CUresult cuMemGetInfo(size_t * free, size_t * total)
    'cuMemGetInfo' : (c_int, POINTER(c_size_t), POINTER(c_size_t)),

    # CUresult cuEventCreate (	CUevent * 	phEvent,
    #                               unsigned int 	Flags )
    'cuEventCreate': (c_int, POINTER(cu_event), c_uint),

    # CUresult cuEventDestroy (	CUevent 	hEvent	 )
    'cuEventDestroy': (c_int, cu_event),

    # CUresult cuEventElapsedTime (	float * 	pMilliseconds,
    #                                   CUevent 	hStart,
    #                                   CUevent 	hEnd )
    'cuEventElapsedTime': (c_int, POINTER(c_float), cu_event, cu_event),

    # CUresult cuEventQuery (	CUevent 	hEvent	 )
    'cuEventQuery': (c_int, cu_event),

    # CUresult cuEventRecord (	CUevent 	hEvent,
    #                               CUstream 	hStream )
    'cuEventRecord': (c_int, cu_event, cu_stream),

    # CUresult cuEventSynchronize (	CUevent 	hEvent	 )
    'cuEventSynchronize': (c_int, cu_event),


    # CUresult cuStreamWaitEvent (	CUstream        hStream,
    #                                   CUevent         hEvent,
    #                                	unsigned int 	Flags )
    'cuStreamWaitEvent': (c_int, cu_stream, cu_event, c_uint),

    # CUresult 	cuPointerGetAttribute (
    #               void *data,
    #               CUpointer_attribute attribute,
    #               CUdeviceptr ptr)
    'cuPointerGetAttribute': (c_int, c_void_p, c_uint, cu_device_ptr),

    #    CUresult cuMemGetAddressRange (	CUdeviceptr * 	pbase,
    #                                        size_t * 	psize,
    #                                        CUdeviceptr 	dptr
    #                                        )
    'cuMemGetAddressRange': (c_int,
                             POINTER(cu_device_ptr),
                             POINTER(c_size_t),
                             cu_device_ptr),

    #    CUresult cuMemHostGetFlags (	unsigned int * 	pFlags,
    #                                   void * 	p )
    'cuMemHostGetFlags': (c_int,
                          POINTER(c_uint),
                          c_void_p),

    #   CUresult cuCtxSynchronize ( void )
    'cuCtxSynchronize' : (c_int,),

    #    CUresult
    #    cuLinkCreate(unsigned int numOptions, CUjit_option *options,
    #                 void **optionValues, CUlinkState *stateOut);
    'cuLinkCreate': (c_int,
                     c_uint, POINTER(cu_jit_option),
                     POINTER(c_void_p), POINTER(cu_link_state)),

    #    CUresult
    #    cuLinkAddData(CUlinkState state, CUjitInputType type, void *data,
    #                  size_t size, const char *name, unsigned
    #                  int numOptions, CUjit_option *options,
    #                  void **optionValues);
    'cuLinkAddData': (c_int,
                      cu_link_state, cu_jit_input_type, c_void_p,
                      c_size_t, c_char_p, c_uint, POINTER(cu_jit_option),
                      POINTER(c_void_p)),

    #    CUresult
    #    cuLinkAddFile(CUlinkState state, CUjitInputType type,
    #                  const char *path, unsigned int numOptions,
    #                  CUjit_option *options, void **optionValues);

    'cuLinkAddFile': (c_int,
                      cu_link_state, cu_jit_input_type, c_char_p, c_uint,
                      POINTER(cu_jit_option), POINTER(c_void_p)),

    #    CUresult CUDAAPI
    #    cuLinkComplete(CUlinkState state, void **cubinOut, size_t *sizeOut)
    'cuLinkComplete': (c_int,
                       cu_link_state, POINTER(c_void_p), POINTER(c_size_t)),

    #    CUresult CUDAAPI
    #    cuLinkDestroy(CUlinkState state)
    'cuLinkDestroy': (c_int, cu_link_state),

    # cuProfilerStart ( void )
    'cuProfilerStart': (c_int,),

    # cuProfilerStop ( void )
    'cuProfilerStop': (c_int,),

    # CUresult cuFuncGetAttribute ( int* pi, CUfunction_attribute attrib,
    #                              CUfunction hfunc )
    'cuFuncGetAttribute': (c_int,
                           POINTER(c_int), cu_function_attribute, cu_function),

    # CUresult CUDAAPI cuOccupancyMaxActiveBlocksPerMultiprocessor(
    #                      int *numBlocks,
    #                      CUfunction func,
    #                      int blockSize,
    #                      size_t dynamicSMemSize);
    'cuOccupancyMaxActiveBlocksPerMultiprocessor': (c_int, POINTER(c_int),
                                                    cu_function, c_size_t,
                                                    c_uint),

    # CUresult CUDAAPI cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(
    #                      int *numBlocks,
    #                      CUfunction func,
    #                      int blockSize,
    #                      size_t dynamicSMemSize,
    #                      unsigned int flags);
    'cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags': (c_int,
                                                             POINTER(c_int),
                                                             cu_function,
                                                             c_size_t, c_uint),

    # CUresult CUDAAPI cuOccupancyMaxPotentialBlockSize(
    #                      int *minGridSize, int *blockSize,
    #                      CUfunction func,
    #                      CUoccupancyB2DSize blockSizeToDynamicSMemSize,
    #                      size_t dynamicSMemSize, int blockSizeLimit);
    'cuOccupancyMaxPotentialBlockSize': (c_int, POINTER(c_int), POINTER(c_int),
                                         cu_function, cu_occupancy_b2d_size,
                                         c_size_t, c_int),

    # CUresult CUDAAPI cuOccupancyMaxPotentialBlockSizeWithFlags(
    #                      int *minGridSize, int *blockSize,
    #                      CUfunction func,
    #                      CUoccupancyB2DSize blockSizeToDynamicSMemSize,
    #                      size_t dynamicSMemSize, int blockSizeLimit,
    #                      unsigned int flags);
    'cuOccupancyMaxPotentialBlockSizeWithFlags': (c_int, POINTER(c_int),
                                                  POINTER(c_int), cu_function,
                                                  cu_occupancy_b2d_size,
                                                  c_size_t, c_int, c_uint),

    # CUresult cuIpcGetMemHandle ( CUipcMemHandle* pHandle, CUdeviceptr dptr )
    'cuIpcGetMemHandle': (c_int,
                          POINTER(cu_ipc_mem_handle), cu_device_ptr),

    # CUresult cuIpcOpenMemHandle(
    #              CUdeviceptr* pdptr,
    #              CUipcMemHandle handle,
    #              unsigned int Flags)
    'cuIpcOpenMemHandle': (c_int, POINTER(cu_device_ptr), cu_ipc_mem_handle,
                           c_uint),

    # CUresult cuIpcCloseMemHandle ( CUdeviceptr dptr )

    'cuIpcCloseMemHandle': (c_int, cu_device_ptr),

    # CUresult cuCtxEnablePeerAccess (CUcontext peerContext, unsigned int Flags)
    'cuCtxEnablePeerAccess': (c_int, cu_context, c_int),

    # CUresult cuDeviceCanAccessPeer ( int* canAccessPeer,
    #                                  CUdevice dev, CUdevice peerDev )
    'cuDeviceCanAccessPeer': (c_int,
                              POINTER(c_int), cu_device, cu_device),

    # CUresult cuDeviceGetUuid ( CUuuid* uuid, CUdevice dev )
    'cuDeviceGetUuid': (c_int, POINTER(cu_uuid), cu_device),
}
