from cpython cimport (
    Py_buffer,
    PyLong_AsVoidPtr,
)
from libc.stdint cimport uint8_t


cdef class CBuffer:
    cdef object _x_handle

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(uint8_t)
        cdef Py_ssize_t[1] shape = tuple((self.bufsize // itemsize,))
        cdef Py_ssize_t[1] strides = tuple((itemsize,))
        buffer.buf = PyLong_AsVoidPtr(self.ptr)
        # assumes sizeof(unsigned char) == sizeof(uint8_t)
        # TODO: use C11 static_assert macro in Cython
        buffer.format = "B"
        buffer.itemsize = itemsize
        buffer.len = self.bufsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL
        self._x_handle = self._x

    def __releasebuffer__(self, Py_buffer *buffer):
        pass
