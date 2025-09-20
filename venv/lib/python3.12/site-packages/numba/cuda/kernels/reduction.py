"""
A library written in CUDA Python for generating reduction kernels
"""

from numba.np.numpy_support import from_dtype


_WARPSIZE = 32
_NUMWARPS = 4


def _gpu_reduce_factory(fn, nbtype):
    from numba import cuda

    reduce_op = cuda.jit(device=True)(fn)
    inner_sm_size = _WARPSIZE + 1   # plus one to avoid SM collision
    max_blocksize = _NUMWARPS * _WARPSIZE

    @cuda.jit(device=True)
    def inner_warp_reduction(sm_partials, init):
        """
        Compute reduction within a single warp
        """
        tid = cuda.threadIdx.x
        warpid = tid // _WARPSIZE
        laneid = tid % _WARPSIZE

        sm_this = sm_partials[warpid, :]
        sm_this[laneid] = init
        cuda.syncwarp()

        width = _WARPSIZE // 2
        while width:
            if laneid < width:
                old = sm_this[laneid]
                sm_this[laneid] = reduce_op(old, sm_this[laneid + width])
            cuda.syncwarp()
            width //= 2

    @cuda.jit(device=True)
    def device_reduce_full_block(arr, partials, sm_partials):
        """
        Partially reduce `arr` into `partials` using `sm_partials` as working
        space.  The algorithm goes like:

            array chunks of 128:  |   0 | 128 | 256 | 384 | 512 |
                        block-0:  |   x |     |     |   x |     |
                        block-1:  |     |   x |     |     |   x |
                        block-2:  |     |     |   x |     |     |

        The array is divided into chunks of 128 (size of a threadblock).
        The threadblocks consumes the chunks in roundrobin scheduling.
        First, a threadblock loads a chunk into temp memory.  Then, all
        subsequent chunks are combined into the temp memory.

        Once all chunks are processed.  Inner-block reduction is performed
        on the temp memory.  So that, there will just be one scalar result
        per block.  The result from each block is stored to `partials` at
        the dedicated slot.
        """
        tid = cuda.threadIdx.x
        blkid = cuda.blockIdx.x
        blksz = cuda.blockDim.x
        gridsz = cuda.gridDim.x

        # block strided loop to compute the reduction
        start = tid + blksz * blkid
        stop = arr.size
        step = blksz * gridsz

        # load first value
        tmp = arr[start]
        # loop over all values in block-stride
        for i in range(start + step, stop, step):
            tmp = reduce_op(tmp, arr[i])

        cuda.syncthreads()
        # inner-warp reduction
        inner_warp_reduction(sm_partials, tmp)

        cuda.syncthreads()
        # at this point, only the first slot for each warp in tsm_partials
        # is valid.

        # finish up block reduction
        # warning: this is assuming 4 warps.
        # assert numwarps == 4
        if tid < 2:
            sm_partials[tid, 0] = reduce_op(sm_partials[tid, 0],
                                            sm_partials[tid + 2, 0])
            cuda.syncwarp()
        if tid == 0:
            partials[blkid] = reduce_op(sm_partials[0, 0], sm_partials[1, 0])

    @cuda.jit(device=True)
    def device_reduce_partial_block(arr, partials, sm_partials):
        """
        This computes reduction on `arr`.
        This device function must be used by 1 threadblock only.
        The blocksize must match `arr.size` and must not be greater than 128.
        """
        tid = cuda.threadIdx.x
        blkid = cuda.blockIdx.x
        blksz = cuda.blockDim.x
        warpid = tid // _WARPSIZE
        laneid = tid % _WARPSIZE

        size = arr.size
        # load first value
        tid = cuda.threadIdx.x
        value = arr[tid]
        sm_partials[warpid, laneid] = value

        cuda.syncthreads()

        if (warpid + 1) * _WARPSIZE < size:
            # fully populated warps
            inner_warp_reduction(sm_partials, value)
        else:
            # partially populated warps
            # NOTE: this uses a very inefficient sequential algorithm
            if laneid == 0:
                sm_this = sm_partials[warpid, :]
                base = warpid * _WARPSIZE
                for i in range(1, size - base):
                    sm_this[0] = reduce_op(sm_this[0], sm_this[i])

        cuda.syncthreads()
        # finish up
        if tid == 0:
            num_active_warps = (blksz + _WARPSIZE - 1) // _WARPSIZE

            result = sm_partials[0, 0]
            for i in range(1, num_active_warps):
                result = reduce_op(result, sm_partials[i, 0])

            partials[blkid] = result

    def gpu_reduce_block_strided(arr, partials, init, use_init):
        """
        Perform reductions on *arr* and writing out partial reduction result
        into *partials*.  The length of *partials* is determined by the
        number of threadblocks. The initial value is set with *init*.

        Launch config:

        Blocksize must be multiple of warpsize and it is limited to 4 warps.
        """
        tid = cuda.threadIdx.x

        sm_partials = cuda.shared.array((_NUMWARPS, inner_sm_size),
                                        dtype=nbtype)
        if cuda.blockDim.x == max_blocksize:
            device_reduce_full_block(arr, partials, sm_partials)
        else:
            device_reduce_partial_block(arr, partials, sm_partials)
        # deal with the initializer
        if use_init and tid == 0 and cuda.blockIdx.x == 0:
            partials[0] = reduce_op(partials[0], init)

    return cuda.jit(gpu_reduce_block_strided)


class Reduce(object):
    """Create a reduction object that reduces values using a given binary
    function. The binary function is compiled once and cached inside this
    object. Keeping this object alive will prevent re-compilation.
    """

    _cache = {}

    def __init__(self, functor):
        """
        :param functor: A function implementing a binary operation for
                        reduction. It will be compiled as a CUDA device
                        function using ``cuda.jit(device=True)``.
        """
        self._functor = functor

    def _compile(self, dtype):
        key = self._functor, dtype
        if key in self._cache:
            kernel = self._cache[key]
        else:
            kernel = _gpu_reduce_factory(self._functor, from_dtype(dtype))
            self._cache[key] = kernel
        return kernel

    def __call__(self, arr, size=None, res=None, init=0, stream=0):
        """Performs a full reduction.

        :param arr: A host or device array.
        :param size: Optional integer specifying the number of elements in
                    ``arr`` to reduce. If this parameter is not specified, the
                    entire array is reduced.
        :param res: Optional device array into which to write the reduction
                    result to. The result is written into the first element of
                    this array. If this parameter is specified, then no
                    communication of the reduction output takes place from the
                    device to the host.
        :param init: Optional initial value for the reduction, the type of which
                    must match ``arr.dtype``.
        :param stream: Optional CUDA stream in which to perform the reduction.
                    If no stream is specified, the default stream of 0 is
                    used.
        :return: If ``res`` is specified, ``None`` is returned. Otherwise, the
                result of the reduction is returned.
        """
        from numba import cuda

        # ensure 1d array
        if arr.ndim != 1:
            raise TypeError("only support 1D array")

        # adjust array size
        if size is not None:
            arr = arr[:size]

        init = arr.dtype.type(init)  # ensure the right type

        # return `init` if `arr` is empty
        if arr.size < 1:
            return init

        kernel = self._compile(arr.dtype)

        # Perform the reduction on the GPU
        blocksize = _NUMWARPS * _WARPSIZE
        size_full = (arr.size // blocksize) * blocksize
        size_partial = arr.size - size_full
        full_blockct = min(size_full // blocksize, _WARPSIZE * 2)

        # allocate size of partials array
        partials_size = full_blockct
        if size_partial:
            partials_size += 1
        partials = cuda.device_array(shape=partials_size, dtype=arr.dtype)

        if size_full:
            # kernel for the fully populated threadblocks
            kernel[full_blockct, blocksize, stream](arr[:size_full],
                                                    partials[:full_blockct],
                                                    init,
                                                    True)

        if size_partial:
            # kernel for partially populated threadblocks
            kernel[1, size_partial, stream](arr[size_full:],
                                            partials[full_blockct:],
                                            init,
                                            not full_blockct)

        if partials.size > 1:
            # finish up
            kernel[1, partials_size, stream](partials, partials, init, False)

        # handle return value
        if res is not None:
            res[:1].copy_to_device(partials[:1], stream=stream)
            return
        else:
            return partials[0]
