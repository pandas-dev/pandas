import numpy as np

from pandas import (
    DataFrame,
    date_range,
    read_stata,
)

from ..pandas_vb_common import (
    BaseIO,
    tm,
)


class Stata(BaseIO):

    params = ["tc", "td", "tm", "tw", "th", "tq", "ty"]
    param_names = ["convert_dates"]

    def setup(self, convert_dates):
        self.fname = "__test__.dta"
        N = self.N = 100000
        C = self.C = 5
        self.df = DataFrame(
            np.random.randn(N, C),
            columns=[f"float{i}" for i in range(C)],
            index=date_range("20000101", periods=N, freq="H"),
        )
        self.df["object"] = tm.makeStringIndex(self.N)
        self.df["int8_"] = np.random.randint(
            np.iinfo(np.int8).min, np.iinfo(np.int8).max - 27, N
        )
        self.df["int16_"] = np.random.randint(
            np.iinfo(np.int16).min, np.iinfo(np.int16).max - 27, N
        )
        self.df["int32_"] = np.random.randint(
            np.iinfo(np.int32).min, np.iinfo(np.int32).max - 27, N
        )
        self.df["float32_"] = np.array(np.random.randn(N), dtype=np.float32)
        self.convert_dates = {"index": convert_dates}
        self.df.to_stata(self.fname, self.convert_dates)

    def time_read_stata(self, convert_dates):
        read_stata(self.fname)

    def time_write_stata(self, convert_dates):
        self.df.to_stata(self.fname, self.convert_dates)


class StataMissing(Stata):
    def setup(self, convert_dates):
        super().setup(convert_dates)
        for i in range(10):
            missing_data = np.random.randn(self.N)
            missing_data[missing_data < 0] = np.nan
            self.df[f"missing_{i}"] = missing_data
        self.df.to_stata(self.fname, self.convert_dates)


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
