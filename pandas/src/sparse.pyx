from numpy cimport ndarray, int32_t, float64_t
cimport numpy as np

cimport cython

import numpy as np
import operator
import sys

np.import_array()
np.import_ufunc()

#-------------------------------------------------------------------------------
# Preamble stuff

cdef float64_t NaN = <float64_t> np.NaN
cdef float64_t INF = <float64_t> np.inf

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

#-------------------------------------------------------------------------------


cdef class SparseIndex:
    '''
    Abstract superclass for sparse index types
    '''
    def __init__(self):
        raise NotImplementedError


cdef class IntIndex(SparseIndex):
    '''
    Object for holding exact integer sparse indexing information

    Parameters
    ----------
    length : integer
    indices : array-like
        Contains integers corresponding to
    '''
    cdef readonly:
        Py_ssize_t length, npoints
        ndarray indices

    def __init__(self, Py_ssize_t length, indices):
        self.length = length
        self.indices = np.ascontiguousarray(indices, dtype=np.int32)
        self.npoints = len(self.indices)

    def __reduce__(self):
        args = (self.length, self.indices)
        return (IntIndex, args)

    def __repr__(self):
        output = 'IntIndex\n'
        output += 'Indices: %s\n' % repr(self.indices)
        return output

    def check_integrity(self):
        '''
        Only need be strictly ascending and nothing less than 0 or greater than
        totall ength
        '''
        pass

    def equals(self, other):
        if not isinstance(other, IntIndex):
            return False

        if self is other:
            return True

        same_length = self.length == other.length
        same_indices = np.array_equal(self.indices, other.indices)
        return same_length and same_indices

    def to_int_index(self):
        return self

    def to_block_index(self):
        locs, lens = get_blocks(self.indices)
        return BlockIndex(self.length, locs, lens)

    cpdef IntIndex intersect(self, SparseIndex y_):
        cdef:
            Py_ssize_t out_length, xi, yi = 0
            int32_t xind
            ndarray[int32_t, ndim=1] xindices, yindices
            list new_list = []
            IntIndex y

        # if is one already, returns self
        y = y_.to_int_index()

        if self.length != y.length:
            raise Exception('Indices must reference same underlying length')

        xindices = self.indices
        yindices = y.indices

        for xi from 0 <= xi < self.npoints:
            xind = xindices[xi]

            while yi < y.npoints and yindices[yi] < xind:
                yi += 1

            if yi >= y.npoints:
                break

            # TODO: would a two-pass algorithm be faster?
            if yindices[yi] == xind:
                new_list.append(xind)

        return IntIndex(self.length, new_list)

    cpdef IntIndex make_union(self, SparseIndex y_):
        cdef:
            Py_ssize_t out_length, i, xi, yi
            int32_t xind
            ndarray[int32_t, ndim=1] xindices, yindices
            list new_list = []
            IntIndex x, y

        x = self

        # if is one already, returns self
        y = y_.to_int_index()

        if self.length != y.length:
            raise Exception('Indices must reference same underlying length')

        xindices = self.indices
        yindices = y.indices

        xi = yi = 0
        while True:
            if xi == x.npoints:
                while yi < y.npoints:
                    new_list.append(yindices[yi])
                    yi += 1
                break
            elif yi == y.npoints:
                while xi < x.npoints:
                    new_list.append(xindices[xi])
                    xi += 1
                break

            xind = xindices[xi]
            yind = yindices[yi]

            if xind == yind:
                new_list.append(xind)
                xi += 1
                yi += 1
            elif xind < yind:
                new_list.append(xind)
                xi += 1
            else:
                new_list.append(yind)
                yi += 1

        return IntIndex(x.length, new_list)

    @cython.wraparound(False)
    cpdef lookup(self, Py_ssize_t index):
        cdef:
            Py_ssize_t res, n, cum_len = 0
            ndarray[int32_t, ndim=1] inds

        inds = self.indices
        res = inds.searchsorted(index)
        if res == self.npoints:
            return -1
        elif inds[res] == index:
            return res
        else:
            return -1

    cpdef ndarray reindex(self, ndarray[float64_t, ndim=1] values,
                          float64_t fill_value, SparseIndex other_):
        cdef:
            Py_ssize_t i = 0, j = 0
            IntIndex other
            ndarray[float64_t, ndim=1] result
            ndarray[int32_t, ndim=1] sinds, oinds

        other = other_.to_int_index()

        oinds = other.indices
        sinds = self.indices

        result = np.empty(other.npoints, dtype=np.float64)
        result.fill(fill_value)

        for 0 <= i < other.npoints:
            while oinds[i] > sinds[j] and j < self.npoints:
                j += 1

            if j == self.npoints:
                break

            if oinds[i] < sinds[j]:
                continue
            elif oinds[i] == sinds[j]:
                result[i] = values[j]
                j += 1

        return result

    cpdef put(self, ndarray[float64_t, ndim=1] values,
              ndarray[int32_t, ndim=1] indices, object to_put):
        pass

    cpdef take(self, ndarray[float64_t, ndim=1] values,
               ndarray[int32_t, ndim=1] indices):
        pass

cpdef get_blocks(ndarray[int32_t, ndim=1] indices):
    cdef:
        Py_ssize_t i, npoints
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
    Object for holding block-based sparse indexing information

    Parameters
    ----------
    '''
    cdef readonly:
        Py_ssize_t nblocks, npoints, length
        ndarray blocs, blengths

    cdef:
        object __weakref__ # need to be picklable
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

    def __reduce__(self):
        args = (self.length, self.blocs, self.blengths)
        return (BlockIndex, args)

    def __repr__(self):
        output = 'BlockIndex\n'
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
            Py_ssize_t i
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
            if blengths[i] == 0:
                raise ValueError('Zero-length block %d' % i)

    def equals(self, other):
        if not isinstance(other, BlockIndex):
            return False

        if self is other:
            return True

        same_length = self.length == other.length
        same_blocks = (np.array_equal(self.blocs, other.blocs) and
                       np.array_equal(self.blengths, other.blengths))
        return same_length and same_blocks

    def to_block_index(self):
        return self

    def to_int_index(self):
        cdef:
            Py_ssize_t i = 0, j, b
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
            ndarray[int32_t, ndim=1] xloc, xlen, yloc, ylen

            list out_blocs = []
            list out_blengths = []

            Py_ssize_t xi = 0, yi = 0
            int32_t cur_loc, cur_length, diff

        y = other.to_block_index()

        if self.length != y.length:
            raise Exception('Indices must reference same underlying length')

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

                if ylen[yi] <= diff:
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

                if xlen[xi] <= diff:
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

    cpdef BlockIndex make_union(self, SparseIndex y):
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
        return BlockUnion(self, y.to_block_index()).result

    cpdef lookup(self, Py_ssize_t index):
        '''

        Returns -1 if not found
        '''
        cdef:
            Py_ssize_t i, cum_len
            ndarray[int32_t, ndim=1] locs, lens

        locs = self.blocs
        lens = self.blengths

        if self.nblocks == 0:
            return -1
        elif index < locs[0]:
            return -1

        cum_len = 0
        for i from 0 <= i < self.nblocks:
            if index >= locs[i] and index < locs[i] + lens[i]:
                return cum_len + index - locs[i]
            cum_len += lens[i]

        return -1

    cpdef ndarray reindex(self, ndarray[float64_t, ndim=1] values,
                          float64_t fill_value, SparseIndex other_):
        cdef:
            Py_ssize_t i = 0, j = 0, ocur, ocurlen
            BlockIndex other
            ndarray[float64_t, ndim=1] result
            ndarray[int32_t, ndim=1] slocs, slens, olocs, olens

        other = other_.to_block_index()

        olocs = other.blocs
        olens = other.blengths
        slocs = self.blocs
        slens = self.blengths

        result = np.empty(other.npoints, dtype=np.float64)

        for 0 <= i < other.nblocks:
            ocur = olocs[i]
            ocurlen = olens[i]

            while slocs[j] + slens[j] < ocur:
                j += 1

    cpdef put(self, ndarray[float64_t, ndim=1] values,
              ndarray[int32_t, ndim=1] indices, object to_put):
        pass

    cpdef take(self, ndarray[float64_t, ndim=1] values,
               ndarray[int32_t, ndim=1] indices):
        pass


cdef class BlockMerge(object):
    '''
    Object-oriented approach makes sharing state between recursive functions a
    lot easier and reduces code duplication
    '''
    cdef:
        BlockIndex x, y, result
        ndarray xstart, xlen, xend, ystart, ylen, yend
        int32_t xi, yi # block indices

    def __init__(self, BlockIndex x, BlockIndex y):
        self.x = x
        self.y = y

        if x.length != y.length:
            raise Exception('Indices must reference same underlying length')

        self.xstart = self.x.blocs
        self.ystart = self.y.blocs

        self.xend = self.x.blocs + self.x.blengths
        self.yend = self.y.blocs + self.y.blengths

        # self.xlen = self.x.blengths
        # self.ylen = self.y.blengths

        self.xi = 0
        self.yi = 0

        self.result = self._make_merged_blocks()

    cdef _make_merged_blocks(self):
        raise NotImplementedError

    cdef _set_current_indices(self, int32_t xi, int32_t yi, bint mode):
        if mode == 0:
            self.xi = xi
            self.yi = yi
        else:
            self.xi = yi
            self.yi = xi

cdef class BlockIntersection(BlockMerge):
    '''
    not done yet
    '''
    pass

cdef class BlockUnion(BlockMerge):
    '''
    Object-oriented approach makes sharing state between recursive functions a
    lot easier and reduces code duplication
    '''

    cdef _make_merged_blocks(self):
        cdef:
            ndarray[int32_t, ndim=1] xstart, xend, ystart, yend
            int32_t nstart, nend, diff
            list out_blocs = [], out_blengths = []

        xstart = self.xstart
        xend = self.xend
        ystart = self.ystart
        yend = self.yend

        while True:
            # we are done (or possibly never began)
            if self.xi >= self.x.nblocks and self.yi >= self.y.nblocks:
                break
            elif self.yi >= self.y.nblocks:
                # through with y, just pass through x blocks
                nstart = xstart[self.xi]
                nend = xend[self.xi]
                self.xi += 1
            elif self.xi >= self.x.nblocks:
                # through with x, just pass through y blocks
                nstart = ystart[self.yi]
                nend = yend[self.yi]
                self.yi += 1
            else:
                # find end of new block
                if xstart[self.xi] < ystart[self.yi]:
                    nstart = xstart[self.xi]
                    nend = self._find_next_block_end(0)
                else:
                    nstart = ystart[self.yi]
                    nend = self._find_next_block_end(1)

            out_blocs.append(nstart)
            out_blengths.append(nend - nstart)

        return BlockIndex(self.x.length, out_blocs, out_blengths)

    cdef int32_t _find_next_block_end(self, bint mode) except -1:
        '''
        Wow, this got complicated in a hurry

        mode 0: block started in index x
        mode 1: block started in index y
        '''
        cdef:
            ndarray[int32_t, ndim=1] xstart, xend, ystart, yend
            int32_t xi, yi, xnblocks, ynblocks, nend

        if mode != 0 and mode != 1:
            raise Exception('Mode must be 0 or 1')

        # so symmetric code will work
        if mode == 0:
            xstart = self.xstart
            xend = self.xend
            xi = self.xi

            ystart = self.ystart
            yend = self.yend
            yi = self.yi
            ynblocks = self.y.nblocks
        else:
            xstart = self.ystart
            xend = self.yend
            xi = self.yi

            ystart = self.xstart
            yend = self.xend
            yi = self.xi
            ynblocks = self.x.nblocks

        nend = xend[xi]

        # print 'here xi=%d, yi=%d, mode=%d, nend=%d' % (self.xi, self.yi,
        #                                                mode, nend)

        # done with y?
        if yi == ynblocks:
            self._set_current_indices(xi + 1, yi, mode)
            return nend
        elif nend < ystart[yi]:
            # block ends before y block
            self._set_current_indices(xi + 1, yi, mode)
            return nend
        else:
            while yi < ynblocks and nend > yend[yi]:
                yi += 1

            self._set_current_indices(xi + 1, yi, mode)

            if yi == ynblocks:
                return nend

            if nend < ystart[yi]:
                # we're done, return the block end
                return nend
            else:
                # merge blocks, continue searching
                # this also catches the case where blocks
                return self._find_next_block_end(1 - mode)


#-------------------------------------------------------------------------------
# Sparse arithmetic

ctypedef float64_t (* double_func)(float64_t a, float64_t b)

cdef inline tuple sparse_nancombine(ndarray x, SparseIndex xindex,
                                    ndarray y, SparseIndex yindex,
                                    double_func op):
    # faster to convert to IntIndex
    return int_nanop(x, xindex.to_int_index(),
                     y, yindex.to_int_index(), op)

    # if isinstance(xindex, BlockIndex):
    #     return block_nanop(x, xindex.to_block_index(),
    #                        y, yindex.to_block_index(), op)
    # elif isinstance(xindex, IntIndex):
    #     return int_nanop(x, xindex.to_int_index(),
    #                      y, yindex.to_int_index(), op)


cdef inline tuple sparse_combine(ndarray x, SparseIndex xindex, float64_t xfill,
                                 ndarray y, SparseIndex yindex, float64_t yfill,
                                 double_func op):
    if isinstance(xindex, BlockIndex):
        return block_op(x, xindex.to_block_index(), xfill,
                        y, yindex.to_block_index(), yfill, op)
    elif isinstance(xindex, IntIndex):
        return int_op(x, xindex.to_int_index(), xfill,
                      y, yindex.to_int_index(), yfill, op)

# NaN-based arithmetic operation-- no handling of fill values
# TODO: faster to convert everything to dense?

@cython.boundscheck(False)
cdef inline tuple block_nanop(ndarray x_, BlockIndex xindex,
                              ndarray y_, BlockIndex yindex,
                              double_func op):
    cdef:
        BlockIndex out_index
        Py_ssize_t xi = 0, yi = 0, out_i = 0 # fp buf indices
        Py_ssize_t xbp = 0, ybp = 0, obp = 0 # block positions
        Py_ssize_t xblock = 0, yblock = 0, outblock = 0 # block numbers

        ndarray[float64_t, ndim=1] x, y
        ndarray[float64_t, ndim=1] out

    # suppress Cython compiler warnings due to inlining
    x = x_
    y = y_

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

@cython.boundscheck(False)
cdef inline tuple int_nanop(ndarray x_, IntIndex xindex,
                            ndarray y_, IntIndex yindex,
                            double_func op):
    cdef:
        IntIndex out_index
        Py_ssize_t xi = 0, yi = 0, out_i = 0 # fp buf indices
        ndarray[int32_t, ndim=1] xindices, yindices, out_indices
        ndarray[float64_t, ndim=1] x, y
        ndarray[float64_t, ndim=1] out

    # suppress Cython compiler warnings due to inlining
    x = x_
    y = y_

    # need to do this first to know size of result array
    out_index = xindex.intersect(yindex)
    out = np.empty(out_index.npoints, dtype=np.float64)

    xindices = xindex.indices
    yindices = yindex.indices
    out_indices = out_index.indices

    # walk the two SparseVectors, adding matched locations...
    for out_i from 0 <= out_i < out_index.npoints:

        # walk x
        while xindices[xi] < out_indices[out_i]:
            xi += 1

        # walk y
        while yindices[yi] < out_indices[out_i]:
            yi += 1

        out[out_i] = op(x[xi], y[yi])

        # advance
        xi += 1
        yi += 1

    return out, out_index


@cython.boundscheck(False)
cdef inline tuple block_op(ndarray x_, BlockIndex xindex, float64_t xfill,
                           ndarray y_, BlockIndex yindex, float64_t yfill,
                           double_func op):
    '''
    Binary operator on BlockIndex objects with fill values
    '''

    cdef:
        BlockIndex out_index
        Py_ssize_t xi = 0, yi = 0, out_i = 0 # fp buf indices
        Py_ssize_t xbp = 0, ybp = 0 # block positions
        int32_t xloc, yloc
        Py_ssize_t xblock = 0, yblock = 0 # block numbers

        ndarray[float64_t, ndim=1] x, y
        ndarray[float64_t, ndim=1] out

    # to suppress Cython warning
    x = x_
    y = y_

    out_index = xindex.make_union(yindex)
    out = np.empty(out_index.npoints, dtype=np.float64)

    # Wow, what a hack job. Need to do something about this

    # walk the two SparseVectors, adding matched locations...
    for out_i from 0 <= out_i < out_index.npoints:
        if yblock == yindex.nblocks:
            # use y fill value
            out[out_i] = op(x[xi], yfill)
            xi += 1

            # advance x location
            xbp += 1
            if xbp == xindex.lenbuf[xblock]:
                xblock += 1
                xbp = 0
            continue

        if xblock == xindex.nblocks:
            # use x fill value
            out[out_i] = op(xfill, y[yi])
            yi += 1

            # advance y location
            ybp += 1
            if ybp == yindex.lenbuf[yblock]:
                yblock += 1
                ybp = 0
            continue

        yloc = yindex.locbuf[yblock] + ybp
        xloc = xindex.locbuf[xblock] + xbp

        # each index in the out_index had to come from either x, y, or both
        if xloc == yloc:
            out[out_i] = op(x[xi], y[yi])
            xi += 1
            yi += 1

            # advance both locations
            xbp += 1
            if xbp == xindex.lenbuf[xblock]:
                xblock += 1
                xbp = 0

            ybp += 1
            if ybp == yindex.lenbuf[yblock]:
                yblock += 1
                ybp = 0

        elif xloc < yloc:
            # use y fill value
            out[out_i] = op(x[xi], yfill)
            xi += 1

            # advance x location
            xbp += 1
            if xbp == xindex.lenbuf[xblock]:
                xblock += 1
                xbp = 0
        else:
            # use x fill value
            out[out_i] = op(xfill, y[yi])
            yi += 1

            # advance y location
            ybp += 1
            if ybp == yindex.lenbuf[yblock]:
                yblock += 1
                ybp = 0

    return out, out_index


@cython.boundscheck(False)
cdef inline tuple int_op(ndarray x_, IntIndex xindex, float64_t xfill,
                         ndarray y_, IntIndex yindex, float64_t yfill,
                         double_func op):
    cdef:
        IntIndex out_index
        Py_ssize_t xi = 0, yi = 0, out_i = 0 # fp buf indices
        int32_t xloc, yloc
        ndarray[int32_t, ndim=1] xindices, yindices, out_indices
        ndarray[float64_t, ndim=1] x, y
        ndarray[float64_t, ndim=1] out

    # suppress Cython compiler warnings due to inlining
    x = x_
    y = y_

    # need to do this first to know size of result array
    out_index = xindex.make_union(yindex)
    out = np.empty(out_index.npoints, dtype=np.float64)

    xindices = xindex.indices
    yindices = yindex.indices
    out_indices = out_index.indices

    # walk the two SparseVectors, adding matched locations...
    for out_i from 0 <= out_i < out_index.npoints:
        if xi == xindex.npoints:
            # use x fill value
            out[out_i] = op(xfill, y[yi])
            yi += 1
            continue

        if yi == yindex.npoints:
            # use y fill value
            out[out_i] = op(x[xi], yfill)
            xi += 1
            continue

        xloc = xindices[xi]
        yloc = yindices[yi]

        # each index in the out_index had to come from either x, y, or both
        if xloc == yloc:
            out[out_i] = op(x[xi], y[yi])
            xi += 1
            yi += 1
        elif xloc < yloc:
            # use y fill value
            out[out_i] = op(x[xi], yfill)
            xi += 1
        else:
            # use x fill value
            out[out_i] = op(xfill, y[yi])
            yi += 1

    return out, out_index

cdef inline float64_t __add(float64_t a, float64_t b):
    return a + b

cdef inline float64_t __sub(float64_t a, float64_t b):
    return a - b

cdef inline float64_t __rsub(float64_t a, float64_t b):
    return b - a

cdef inline float64_t __div(float64_t a, float64_t b):
    if b == 0:
        if a >= 0:
            return INF
        else:
            return -INF
    else:
        return a / b

cdef inline float64_t __rdiv(float64_t a, float64_t b):
    return __div(b, a)

cdef inline float64_t __floordiv(float64_t a, float64_t b):
    if b == 0:
        if a >= 0:
            return INF
        else:
            return -INF
    else:
        return a // b

cdef inline float64_t __rfloordiv(float64_t a, float64_t b):
    return __floordiv(b, a)

cdef inline float64_t __mul(float64_t a, float64_t b):
    return a * b
cdef inline float64_t __eq(float64_t a, float64_t b):
    return a == b
cdef inline float64_t __ne(float64_t a, float64_t b):
    return a != b
cdef inline float64_t __lt(float64_t a, float64_t b):
    return a < b
cdef inline float64_t __gt(float64_t a, float64_t b):
    return a > b

cdef inline float64_t __pow(float64_t a, float64_t b):
    # NaN
    if a != a or b != b:
        return NaN
    return a ** b

cdef inline float64_t __rpow(float64_t a, float64_t b):
    return __pow(b, a)


# This probably needs to be "templated" to achieve maximum performance.
# TODO: quantify performance boost to "templating"

cpdef sparse_nanadd(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __add)

cpdef sparse_nansub(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __sub)

cpdef sparse_nanrsub(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __rsub)

cpdef sparse_nanmul(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __mul)

cpdef sparse_nandiv(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __div)

cpdef sparse_nanrdiv(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __rdiv)

sparse_nantruediv = sparse_nandiv
sparse_nanrtruediv = sparse_nanrdiv

cpdef sparse_nanfloordiv(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __floordiv)

cpdef sparse_nanrfloordiv(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __rfloordiv)

cpdef sparse_nanpow(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __pow)

cpdef sparse_nanrpow(ndarray x, SparseIndex xindex,
                    ndarray y, SparseIndex yindex):
    return sparse_nancombine(x, xindex, y, yindex, __rpow)

cpdef sparse_add(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __add)

cpdef sparse_sub(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __sub)

cpdef sparse_rsub(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __rsub)

cpdef sparse_mul(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __mul)

cpdef sparse_div(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __div)

cpdef sparse_rdiv(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __rdiv)

sparse_truediv = sparse_div
sparse_rtruediv = sparse_rdiv

cpdef sparse_floordiv(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __floordiv)

cpdef sparse_rfloordiv(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __rfloordiv)

cpdef sparse_pow(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __pow)

cpdef sparse_rpow(ndarray x, SparseIndex xindex, float64_t xfill,
                 ndarray y, SparseIndex yindex, float64_t yfill):
    return sparse_combine(x, xindex, xfill,
                             y, yindex, yfill, __rpow)


#-------------------------------------------------------------------------------
# Indexing operations

def get_reindexer(ndarray[object, ndim=1] values, dict index_map):
    cdef object idx
    cdef Py_ssize_t i
    cdef Py_ssize_t new_length = len(values)
    cdef ndarray[int32_t, ndim=1] indexer

    indexer = np.empty(new_length, dtype=np.int32)

    for i in range(new_length):
        idx = values[i]
        if idx in index_map:
            indexer[i] = index_map[idx]
        else:
            indexer[i] = -1

    return indexer

# def reindex_block(ndarray[float64_t, ndim=1] values,
#                   BlockIndex sparse_index,
#                   ndarray[int32_t, ndim=1] indexer):
#     cdef:
#         Py_ssize_t i, length
#         ndarray[float64_t, ndim=1] out

#     out = np.empty(length, dtype=np.float64)

#     for i from 0 <= i < length:
#         if indexer[i] == -1:
#             pass


# cdef class SparseCruncher(object):
#     '''
#     Class to acquire float pointer for convenient operations on sparse data
#     structures
#     '''
#     cdef:
#         SparseIndex index
#         float64_t* buf

#     def __init__(self, ndarray[float64_t, ndim=1, mode='c'] values,
#                  SparseIndex index):

#         self.index = index
#         self.buf = <float64_t*> values.data


def reindex_integer(ndarray[float64_t, ndim=1] values,
                    IntIndex sparse_index,
                    ndarray[int32_t, ndim=1] indexer):
    pass
