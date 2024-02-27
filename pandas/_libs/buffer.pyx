from cpython cimport Py_buffer
from libc.stdint cimport (
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)

ctypedef fused supported_buffer_t:
    float
    double
    int8_t
    int16_t
    int32_t
    int64_t
    uint8_t
    uint16_t
    uint32_t
    uint64_t


cdef class PandasBuffer:
    """
    Data in the buffer is guaranteed to be contiguous in memory.

    Note that there is no dtype attribute present, a buffer can be thought of
    as simply a block of memory. However, if the column that the buffer is
    attached to has a dtype that's supported by DLPack and ``__dlpack__`` is
    implemented, then that dtype information will be contained in the return
    value from ``__dlpack__``.

    This distinction is useful to support both data exchange via DLPack on a
    buffer and (b) dtypes like variable-length strings which do not have a
    fixed number of bytes per element.
    """

    def __init__(self, supported_buffer_t[:] buf, allow_copy: bool = True) -> None:
        """
        Handle only regular columns (= numpy arrays) for now.
        """
        if buf.strides[0] and not buf.strides == (buf.dtype.itemsize,):
            # The protocol does not support strided buffers, so a copy is
            # necessary. If that's not allowed, we need to raise an exception.
            if allow_copy:
                buf = buf.copy()
            else:
                raise RuntimeError(
                    "Exports cannot be zero-copy in the case "
                    "of a non-contiguous buffer"
                )

        # Store the numpy array in which the data resides as a private
        # attribute, so we can use it to retrieve the public attributes
        self.buf = buf

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        buffer.buf = &self.buf[0]
        # assumes sizeof(unsigned char) == sizeof(uint8_t)
        # TODO: use C11 static_assert macro in Cython
        buffer.format = self.buf.format
        buffer.itemsize = self.buf.itemsize
        buffer.len = len(self.buf)
        buffer.ndim = self.buf.ndim
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = self.buf.shape
        buffer.strides = self.buf.strides
        buffer.suboffsets = self.buf.suboffsets

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    @property
    def bufsize(self) -> int:
        """
        Buffer size in bytes.
        """
        return self.buf.size * self.buf.itemsize

    @property
    def ptr(self) -> int:
        """
        Pointer to start of the buffer as an integer.
        """
        return &self.buf[0]

    def __dlpack__(self):
        """
        Represent this structure as DLPack interface.
        """
        raise NotImplementedError

    def __dlpack_device__(self):
        """
        Device type and device ID for where the data in the buffer resides.
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        return (
            "PandasBuffer("
            + str(
                {
                    "bufsize": self.bufsize,
                    "ptr": self.ptr,
                }
            )
            + ")"
        )
