import pandas as pd
import os

def create_test_file():
    col_nums = [0]
    df = pd.DataFrame({"V1":["a","b","c","d","e", "aaaah!!!"], 
                              "W":["c","d","c","d","c","c"],
                              "ZZZ":np.arange(6)})
    df.set_index(["V1","W"], inplace = True)
    df.to_csv("testtable.tab",sep = "\t")


def test_write_hdf5_11364():
    sep = "\t"
    indexcols =[0]
    chunksize=5

    xbed = "testtable.tab"
    os.remove(xbed)
    # create a store
    with pd.HDFStore('tempstore.h5') as store:
        for nn, chunk in enumerate(pd.read_table(xbed, chunksize=chunksize, sep = sep, index_col= indexcols)):
            group = "x"
            #print(chunk.index.names)
            store.append(group, chunk, format = "table", min_itemsize =  \
                         {"index":32} if len(indexcols)==1 else \
                         dict(zip(chunk.index.names, [32]*len(chunk.index.names))))
            print("chunk #" , nn, file = sys.stderr)

    os.remove(xbed)
    assert True

def test_read_hdf5_11364():
    with pd.HDFStore('tempstore.h5') as store:
        df = store.get(group)
        print(df.shape)
    assert (df.shape==(6,3 - len(indexcols))), "wrong shape"
