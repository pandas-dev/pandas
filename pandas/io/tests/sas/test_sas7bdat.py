import pandas as pd
from pandas.compat import PY2
import pandas.util.testing as tm
import os
import io
import numpy as np


class TestSAS7BDAT(tm.TestCase):

    def setUp(self):
        self.dirpath = tm.get_data_path()
        self.data = []
        self.test_ix = [list(range(1, 16)), [16]]
        for j in 1, 2:
            fname = os.path.join(self.dirpath, "test_sas7bdat_%d.csv" % j)
            df = pd.read_csv(fname)
            epoch = pd.datetime(1960, 1, 1)
            t1 = pd.to_timedelta(df["Column4"], unit='d')
            df["Column4"] = epoch + t1
            t2 = pd.to_timedelta(df["Column12"], unit='d')
            df["Column12"] = epoch + t2
            for k in range(df.shape[1]):
                col = df.iloc[:, k]
                if col.dtype == np.int64:
                    df.iloc[:, k] = df.iloc[:, k].astype(np.float64)
                elif col.dtype == np.dtype('O'):
                    if PY2:
                        f = lambda x: (x.decode('utf-8') if
                                       isinstance(x, str) else x)
                        df.iloc[:, k] = df.iloc[:, k].apply(f)
            self.data.append(df)

    def test_from_file(self):
        for j in 0, 1:
            df0 = self.data[j]
            for k in self.test_ix[j]:
                fname = os.path.join(self.dirpath, "test%d.sas7bdat" % k)
                df = pd.read_sas(fname, encoding='utf-8')
                tm.assert_frame_equal(df, df0)

    def test_from_buffer(self):
        for j in 0, 1:
            df0 = self.data[j]
            for k in self.test_ix[j]:
                fname = os.path.join(self.dirpath, "test%d.sas7bdat" % k)
                byts = open(fname, 'rb').read()
                buf = io.BytesIO(byts)
                df = pd.read_sas(buf, format="sas7bdat", encoding='utf-8')
                tm.assert_frame_equal(df, df0)

    def test_from_iterator(self):
        for j in 0, 1:
            df0 = self.data[j]
            for k in self.test_ix[j]:
                fname = os.path.join(self.dirpath, "test%d.sas7bdat" % k)
                byts = open(fname, 'rb').read()
                buf = io.BytesIO(byts)
                rdr = pd.read_sas(buf, format="sas7bdat",
                                  iterator=True, encoding='utf-8')
                df = rdr.read(2)
                tm.assert_frame_equal(df, df0.iloc[0:2, :])
                df = rdr.read(3)
                tm.assert_frame_equal(df, df0.iloc[2:5, :])
