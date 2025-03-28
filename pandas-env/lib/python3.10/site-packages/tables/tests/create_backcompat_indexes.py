# Script for creating different kind of indexes in a small space as possible.
# This is intended for testing purposes.

import tables as tb


class Descr(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4, shape=(), dflt='', pos=0)
    var2 = tb.BoolCol(shape=(), dflt=False, pos=1)
    var3 = tb.Int32Col(shape=(), dflt=0, pos=2)
    var4 = tb.Float64Col(shape=(), dflt=0.0, pos=3)


# Parameters for the table and index creation
small_chunkshape = (2,)
small_blocksizes = (64, 32, 16, 8)
nrows = 43

# Create the new file
h5fname = 'indexes_2_1.h5'
h5file = tb.open_file(h5fname, 'w')
t1 = h5file.create_table(h5file.root, 'table1', Descr)
row = t1.row
for i in range(nrows):
    row['var1'] = i
    row['var2'] = i
    row['var3'] = i
    row['var4'] = i
    row.append()
t1.flush()

# Do a copy of table1
t1.copy(h5file.root, 'table2')

# Create indexes of all kinds
t1.cols.var1.create_index(0, 'ultralight', _blocksizes=small_blocksizes)
t1.cols.var2.create_index(3, 'light', _blocksizes=small_blocksizes)
t1.cols.var3.create_index(6, 'medium', _blocksizes=small_blocksizes)
t1.cols.var4.create_index(9, 'full', _blocksizes=small_blocksizes)

h5file.close()
