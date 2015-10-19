import nose
from nose import with_setup
import pandas as pd
import numpy as np
import os, sys

def create_test_file():
    global xbed, xstore, xgroup
    xbed = "testtable.tab"
    xstore = 'tempstore.h5'
    xgroup = "x"

    col_nums = [0]
    df = pd.DataFrame({"V1":["a","b","c","d","e", "aaaah!!!"], 
                              "W":["c","d","c","d","c","c"],
                              "ZZZ":np.arange(6)})
    df.set_index(["V1","W"], inplace = True)
    df.to_csv( xbed, sep = "\t")


def clear_files():
    os.remove(xbed)
    os.remove(xstore)

def write_hdf5_11364(indexcols):
    sep = "\t"
    chunksize=5
    try:
        os.remove(xstore)
    except OSError:
        pass
    # create a store
    with pd.HDFStore(xstore) as store:
        for nn, chunk in enumerate(pd.read_table(xbed, chunksize=chunksize, sep = sep, index_col= indexcols if not  indexcols==["index"] else 0)):
            #print(chunk.index.names)
            store.append(xgroup, chunk, format = "table", min_itemsize =  \
                         #{"index":32} if len(indexcols)==1 else \
                         dict(zip(chunk.index.names, [32]*len(chunk.index.names))))
            print("chunk #" , nn, file = sys.stderr)

    print("index columns:", indexcols, file = sys.stderr)
    assert True

def read_hdf5_11364(indexcols):
    with pd.HDFStore(xstore) as store:
        df = store.get(xgroup)
        print(df.shape)
    assert (df.shape==(6,3 - len(indexcols))), "wrong shape"

@with_setup(create_test_file, clear_files )
def test_write_read_hdf5_11364_indexcol():
    indexcols = ["index"]
    write_hdf5_11364(indexcols)
    read_hdf5_11364(indexcols)
    return

@with_setup(create_test_file, clear_files )
def test_write_read_hdf5_11364_1col():
    indexcols =[0]
    write_hdf5_11364(indexcols)
    read_hdf5_11364(indexcols)
    return

@with_setup(create_test_file, clear_files )
def test_write_read_hdf5_11364_2col():
    indexcols =[0,1]
    write_hdf5_11364(indexcols)
    read_hdf5_11364(indexcols)
    return


