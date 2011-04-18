from numpy cimport ndarray, int32_t, float64_t
cimport numpy as np

cimport cython

import numpy as np
import operator
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

    def __init__(self):
        raise NotImplementedError

cdef class IntIndex(SparseIndex):

    cdef readonly:
        pyst length, npoints
        ndarray indices

    cdef:
        int32_t* indp

    def __init__(self, pyst length, indices):
        self.length = length
        self.indices = np.ascontiguousarray(indices, dtype=np.int32)

        self.npoints = len(self.indices)

        self.indp = <int32_t*> self.indices.data

    def __repr__(self):
        output = 'sparse.IntIndex\n'
        output += 'Indices: %s\n' % repr(self.indices)
        return output

    def equals(self, other):
        if not isinstance(other, IntIndex):
            raise Exception('Can only compare with like object')

        same_length = self.length == other.length
        same_indices = np.array_equal(self.indices, other.indices)
        return same_length and same_indices

    def to_int_index(self):
        return self

    def to_block_index(self):
        locs, lens = get_blocks(self.indices)
        return BlockIndex(self.length, locs, lens)

    cpdef intersect(self, SparseIndex y_):
        cdef:
            pyst i, xi, yi = 0
            int32_t xind
            list new_list = []
            IntIndex y

        # if is one already, returns self
        y = y_.to_int_index()

        for xi from 0 <= xi < self.npoints:
            xind = self.indp[xi]

            while yi < y.npoints and y.indp[yi] < xind:
                yi += 1

            if yi >= y.npoints:
                break

            # TODO: would a two-pass algorithm be faster?
            if y.indp[yi] == xind:
                new_list.append(xind)

        return IntIndex(self.length, new_list)


cpdef get_blocks(ndarray[int32_t, ndim=1] indices):
    cdef:
        pyst i, npoints
        int32_t block, length = 1, cur, prev
        list locs = [], lens = []

    npoints = len(indices)

    # just handle the special empty case separately
    if npoints == 0:
        return [], []

    # TODO: two-pass algorithm faster?
    prev = block = indices[0]
    for i from 1 <= i < npoints:
        cur = indices[i]
        if cur - prev > 1:
            # new block
            locs.append(block)
            lens.append(length)
            block = cur
            length = 1
        else:
            # same block, increment length
            length += 1

        prev = cur

    locs.append(block)
    lens.append(length)
    return locs, lens

#-------------------------------------------------------------------------------
# BlockIndex

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
        output += 'Block lengths: %s' % repr(self.blengths)

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

    def equals(self, other):
        if not isinstance(other, BlockIndex):
            raise Exception('Can only compare with like object')

        same_length = self.length == other.length
        same_blocks = (np.array_equal(self.blocs, other.blocs) and
                       np.array_equal(self.blengths, other.blengths))
        return same_length and same_blocks

    def to_block_index(self):
        return self

    def to_int_index(self):
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

        return IntIndex(self.length, indices)

    cpdef BlockIndex intersect(self, SparseIndex other):
        '''
        Intersect two BlockIndex objects

        Parameters
        ----------

        Returns
        -------
        intersection : BlockIndex
        '''
        cdef:
            BlockIndex y
            pyst out_length
            ndarray[int32_t, ndim=1] xloc, xlen, yloc, ylen

            list out_blocs = []
            list out_blengths = []

            pyst xi = 0, yi = 0
            int32_t cur_loc, cur_length, xend, yend, diff


        y = other.to_block_index()

        # unwise? should enforce same length?
        out_length = int_max(self.length, y.length)

        xloc = self.blocs
        xlen = self.blengths
        yloc = y.blocs
        ylen = y.blengths

        while True:
            # we are done (or possibly never began)
            if xi >= self.nblocks or yi >= y.nblocks:
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

    cpdef BlockIndex make_union(self, SparseIndex other):
        '''
        Combine together two BlockIndex objects, accepting indices if contained
        in one or the other

        Parameters
        ----------
        other : SparseIndex

        Notes
        -----
        union is a protected keyword in Cython, hence make_union

        Returns
        -------
        union : BlockIndex
        '''
        pass

#-------------------------------------------------------------------------------
# Sparse operations

cpdef sparse_nanadd(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __add)

cpdef sparse_nansub(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __sub)

cpdef sparse_nanmul(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __mul)

cpdef sparse_nandiv(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __div)

cdef tuple sparse_nancombine(ndarray x, SparseIndex xindex,
                             ndarray y, SparseIndex yindex, double_func op):
    if isinstance(xindex, BlockIndex):
        return block_nanop(x, xindex.to_block_index(),
                           y, yindex.to_block_index(), op)
    elif isinstance(xindex, IntIndex):
        return int_nanop(x, xindex.to_int_index(),
                         y, yindex.to_int_index(), op)

# NaN-based arithmetic operation-- no handling of fill values
# TODO: faster to convert everything to dense?

cdef tuple block_nanop(ndarray[float64_t, ndim=1] x, BlockIndex xindex,
                       ndarray[float64_t, ndim=1] y, BlockIndex yindex,
                       double_func op):
    cdef:
        BlockIndex out_index
        int xi = 0, yi = 0, out_i = 0 # fp buf indices
        int xbp = 0, ybp = 0, obp = 0 # block positions
        pyst xblock = 0, yblock = 0, outblock = 0 # block numbers

        ndarray[float64_t, ndim=1] out

    out_index = xindex.intersect(yindex)
    out = np.empty(out_index.npoints, dtype=np.float64)

    # walk the two SparseVectors, adding matched locations...
    for out_i from 0 <= out_i < out_index.npoints:

        # I have a feeling this is inefficient

        # walk x
        while xindex.locbuf[xblock] + xbp < out_index.locbuf[outblock] + obp:
            xbp += 1
            xi += 1
            if xbp == xindex.lenbuf[xblock]:
                xblock += 1
                xbp = 0

        # walk y
        while yindex.locbuf[yblock] + ybp < out_index.locbuf[outblock] + obp:
            ybp += 1
            yi += 1
            if ybp == yindex.lenbuf[yblock]:
                yblock += 1
                ybp = 0

        out[out_i] = op(x[xi], y[yi])

        # advance. strikes me as too complicated
        xi += 1
        yi += 1

        xbp += 1
        if xbp == xindex.lenbuf[xblock]:
            xblock += 1
            xbp = 0

        ybp += 1
        if ybp == yindex.lenbuf[yblock]:
            yblock += 1
            ybp = 0

        obp += 1
        if obp == out_index.lenbuf[outblock]:
            outblock += 1
            obp = 0

    return out, out_index

cdef tuple int_nanop(ndarray[float64_t, ndim=1] x, IntIndex xindex,
                     ndarray[float64_t, ndim=1] y, IntIndex yindex,
                     double_func op):
    cdef:
        IntIndex out_index
        int xi = 0, yi = 0, out_i = 0 # fp buf indices
        ndarray[float64_t, ndim=1] out

    # need to do this first to know size of result array
    out_index = xindex.intersect(yindex).to_int_index()
    out = np.empty(out_index.npoints, dtype=np.float64)

    # walk the two SparseVectors, adding matched locations...
    for out_i from 0 <= out_i < out_index.npoints:

        # walk x
        while xindex.indp[xi] < out_index.indp[out_i]:
            xi += 1

        # walk y
        while yindex.indp[yi] < out_index.indp[out_i]:
            yi += 1

        out[out_i] = op(x[xi], y[yi])

        # advance
        xi += 1
        yi += 1

    return out, out_index
