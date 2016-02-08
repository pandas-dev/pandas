import numpy as np
from pandas import *


def panda_test():

    # generate some data
    data = np.random.rand(50, 5)
    # generate some dates
    dates = DatetimeIndex('1/1/1969', periods=50)
    # generate column headings
    cols = ['A', 'B', 'C', 'D', 'E']

    df = DataFrame(data, index=dates, columns=cols)

    # save to HDF5Store
    store = HDFStore('bugzilla.h5', mode='w')
    store['df'] = df  # This gives: OverflowError: mktime argument out of range
    store.close()


if __name__ == '__main__':
    panda_test()
