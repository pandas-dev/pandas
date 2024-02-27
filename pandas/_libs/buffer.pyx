from cpython cimport (
    Py_buffer,
    PyLong_FromVoidPtr,
)
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

    # we cannot use a fused type as a class attribute, so we instead
    # unpack the items we need for the buffer protocol in __init__

    cdef:
        void *ptr_
        Py_ssize_t len_
        Py_ssize_t itemsize

        int readonly
        int ndim
        bytes format
        Py_ssize_t *shape
        Py_ssize_t *strides
        Py_ssize_t *suboffsets

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
        self.ptr_ = &buf[0]
        self.len_ = len(buf)
        self.itemsize = buf.itemsize
        self.readonly = buf.readonly
        self.ndim = buf.ndim
        self.format = buf.format
        self.shape = buf.shape
        self.strides = buf.strides
        self.suboffsets = buf.suboffsets

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        buffer.buf = self.ptr_
        # assumes sizeof(unsigned char) == sizeof(uint8_t)
        # TODO: use C11 static_assert macro in Cython
        buffer.format = self.format
        buffer.itemsize = self.itemsize
        buffer.len = self.len_
        buffer.ndim = self.ndim
        buffer.obj = self.obj
        buffer.readonly = self.readonly
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = self.suboffsets

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
        return PyLong_FromVoidPtr(self.ptr_)

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
