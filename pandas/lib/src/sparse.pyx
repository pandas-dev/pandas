from numpy cimport ndarray, int32_t, float64_t
cimport numpy as np

cimport cython

import numpy as np
import sys

ctypedef Py_ssize_t pyst

#-------------------------------------------------------------------------------
# Preamble stuff

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

cdef float64_t __add(float64_t a, float64_t b):
    return a + b
cdef float64_t __sub(float64_t a, float64_t b):
    return a - b
cdef float64_t __div(float64_t a, float64_t b):
    return a / b
cdef float64_t __mul(float64_t a, float64_t b):
    return a * b
cdef float64_t __eq(float64_t a, float64_t b):
    return a == b
cdef float64_t __ne(float64_t a, float64_t b):
    return a != b
cdef float64_t __lt(float64_t a, float64_t b):
    return a < b
cdef float64_t __gt(float64_t a, float64_t b):
    return a > b
cdef float64_t __pow(float64_t a, float64_t b):
    return a ** b

ctypedef float64_t (* double_func)(float64_t a, float64_t b)

#-------------------------------------------------------------------------------

cdef class SparseIndex:
    pass

cdef class DenseIndex(SparseIndex):

    cdef readonly:
        pyst length
        ndarray indices

    cdef:
        int32_t* indbuf

    def __init__(self, pyst length, indices):
        self.length = length
        self.indices = np.ascontiguousarray(indices, dtype=np.int32)

        self.indbuf = <int32_t*> self.indices.data

    def __repr__(self):
        output = 'sparse.DenseIndex\n'
        output += 'Indices: %s\n' % repr(self.indices)
        return output

    def to_block(self):
        pass

    cpdef intersect(self, SparseIndex):
        pass

cdef class BlockIndex(SparseIndex):
    '''

    '''

    cdef readonly:
        pyst nblocks, npoints, length
        ndarray blocs, blengths

    cdef:
        int32_t* locbuf, *lenbuf

    def __init__(self, length, blocs, blengths):

        self.blocs = np.ascontiguousarray(blocs, dtype=np.int32)
        self.blengths = np.ascontiguousarray(blengths, dtype=np.int32)

        # in case we need
        self.locbuf = <int32_t*> self.blocs.data
        self.lenbuf = <int32_t*> self.blengths.data

        self.length = length
        self.nblocks = len(self.blocs)
        self.npoints = self.blengths.sum()

        # self.block_start = blocs
        # self.block_end = blocs + blengths

        self.check_integrity()

    def __repr__(self):
        output = 'sparse.BlockIndex\n'
        output += 'Block locations: %s\n' % repr(self.blocs)
        output += 'Block lengths: %s\n' % repr(self.blengths)

        return output

    cpdef check_integrity(self):
        '''
        Check:
        - Locations are in ascending order
        - No overlapping blocks
        - Blocks to not start after end of index, nor extend beyond end
        '''
        cdef:
            pyst i
            ndarray[int32_t, ndim=1] blocs, blengths

        blocs = self.blocs
        blengths = self.blengths

        if len(blocs) != len(blengths):
            raise ValueError('block bound arrays must be same length')

        for i from 0 <= i < self.nblocks:
            if i > 0:
                if blocs[i] <= blocs[i-1]:
                    raise ValueError('Locations not in ascending order')

            if i < self.nblocks - 1:
                if blocs[i] + blengths[i] > blocs[i + 1]:
                    raise ValueError('Block %d overlaps' % i)
            else:
                if blocs[i] + blengths[i] > self.length:
                    raise ValueError('Block %d extends beyond end' % i)

            # no zero-length blocks
            if self.blengths[i] == 0:
                raise ValueError('Zero-length block %d' % i)

    def to_dense(self):
        cdef:
            pyst i = 0, j, b
            int32_t offset
            ndarray[int32_t, ndim=1] indices

        indices = np.empty(self.npoints, dtype=np.int32)

        for b from 0 <= b < self.nblocks:
            offset = self.locbuf[b]

            for j from 0 <= j < self.lenbuf[b]:
                indices[i] = offset + j
                i += 1

        return DenseIndex(self.length, indices)

    cpdef BlockIndex intersect(self, BlockIndex other):
        cdef:
            pyst out_length
            ndarray[int32_t, ndim=1] xloc, xlen, yloc, ylen

            list out_blocs = []
            list out_blengths = []

        # unwise? should enforce same length?
        out_length = int_max(self.length, other.length)

        xloc = self.blocs
        xlen = self.blengths
        yloc = other.blocs
        ylen = other.blengths

        cdef pyst xi = 0, yi = 0
        cdef int32_t cur_loc, cur_length, xend, yend, diff
        while True:
            # we are done (or possibly never began)
            if xi >= self.nblocks or yi >= other.nblocks:
                break

            # completely symmetric...would like to avoid code dup but oh well
            if xloc[xi] >= yloc[yi]:
                cur_loc = xloc[xi]
                diff = xloc[xi] - yloc[yi]

                if ylen[yi] - diff <= 0:
                    # have to skip this block
                    yi += 1
                    continue

                if ylen[yi] - diff < xlen[xi]:
                    # take end of y block, move onward
                    cur_length = ylen[yi] - diff
                    yi += 1
                else:
                    # take end of x block
                    cur_length = xlen[xi]
                    xi += 1

            else: # xloc[xi] < yloc[yi]
                cur_loc = yloc[yi]
                diff = yloc[yi] - xloc[xi]

                if xlen[xi] - diff <= 0:
                    # have to skip this block
                    xi += 1
                    continue

                if xlen[xi] - diff < ylen[yi]:
                    # take end of x block, move onward
                    cur_length = xlen[xi] - diff
                    xi += 1
                else:
                    # take end of y block
                    cur_length = ylen[yi]
                    yi += 1

            out_blocs.append(cur_loc)
            out_blengths.append(cur_length)

        return BlockIndex(self.length, out_blocs, out_blengths)

cdef class SparseVector:
    '''
    Data structure for storing sparse representation of floating point data
    '''

    cdef public:
        ndarray values
        SparseIndex index

    cdef:
        float64_t* vbuf

    def __init__(self, ndarray values, SparseIndex index):
        self.values = np.ascontiguousarray(values, dtype=np.float64)
        self.vbuf = <float64_t*> self.values.data

    def __repr__(self):
        pass

    cpdef reindex(self):
        pass

    cpdef add(self, SparseVector other):
        cdef double_func op = __add

    cpdef sub(self, SparseVector other):
        pass

    cpdef mul(self, SparseVector other):
        pass

    cpdef div(self, SparseVector other):
        pass

    cdef ndarray _combine(self, SparseVector other, double_func op):
        cdef SparseIndex out_index = self.index.intersect(other.index)
        cdef ndarray out = np.empty(out_index.npoints, dtype=np.float64)

        if isinstance(out_index, BlockIndex):
            block_op(self.vbuf, other.vbuf, <float64_t*> out.data,
                     op, self.index, other.index)
        elif isinstance(out_index, DenseIndex):
            pass

        return SparseVector(out, out_index)

cdef block_op(float64_t* xbuf, float64_t* ybuf, float64_t* out,
              double_func op,
              BlockIndex xindex, BlockIndex yindex) except -1:
    pass
