import sys

import numpy as np

from pandas import (
    DataFrame,
    concat,
    date_range,
    json_normalize,
    read_json,
    timedelta_range,
)

from ..pandas_vb_common import (
    BaseIO,
    tm,
)


class ReadJSON(BaseIO):

    fname = "__test__.json"
    params = (["split", "index", "records"], ["int", "datetime"])
    param_names = ["orient", "index"]

    def setup(self, orient, index):
        N = 100000
        indexes = {
            "int": np.arange(N),
            "datetime": date_range("20000101", periods=N, freq="H"),
        }
        df = DataFrame(
            np.random.randn(N, 5),
            columns=[f"float_{i}" for i in range(5)],
            index=indexes[index],
        )
        df.to_json(self.fname, orient=orient)

    def time_read_json(self, orient, index):
        read_json(self.fname, orient=orient)


class ReadJSONLines(BaseIO):

    fname = "__test_lines__.json"
    params = ["int", "datetime"]
    param_names = ["index"]

    def setup(self, index):
        N = 100000
        indexes = {
            "int": np.arange(N),
            "datetime": date_range("20000101", periods=N, freq="H"),
        }
        df = DataFrame(
            np.random.randn(N, 5),
            columns=[f"float_{i}" for i in range(5)],
            index=indexes[index],
        )
        df.to_json(self.fname, orient="records", lines=True)

    def time_read_json_lines(self, index):
        read_json(self.fname, orient="records", lines=True)

    def time_read_json_lines_concat(self, index):
        concat(read_json(self.fname, orient="records", lines=True, chunksize=25000))

    def time_read_json_lines_nrows(self, index):
        read_json(self.fname, orient="records", lines=True, nrows=25000)

    def peakmem_read_json_lines(self, index):
        read_json(self.fname, orient="records", lines=True)

    def peakmem_read_json_lines_concat(self, index):
        concat(read_json(self.fname, orient="records", lines=True, chunksize=25000))

    def peakmem_read_json_lines_nrows(self, index):
        read_json(self.fname, orient="records", lines=True, nrows=15000)


class NormalizeJSON(BaseIO):
    fname = "__test__.json"
    params = [
        ["split", "columns", "index", "values", "records"],
        ["df", "df_date_idx", "df_td_int_ts", "df_int_floats", "df_int_float_str"],
    ]
    param_names = ["orient", "frame"]

    def setup(self, orient, frame):
        data = {
            "hello": ["thisisatest", 999898, "mixed types"],
            "nest1": {"nest2": {"nest3": "nest3_value", "nest3_int": 3445}},
            "nest1_list": {"nest2": ["blah", 32423, 546456.876, 92030234]},
            "hello2": "string",
        }
        self.data = [data for i in range(10000)]

    def time_normalize_json(self, orient, frame):
        json_normalize(self.data)


class ToJSON(BaseIO):

    fname = "__test__.json"
    params = [
        ["split", "columns", "index", "values", "records"],
        ["df", "df_date_idx", "df_td_int_ts", "df_int_floats", "df_int_float_str"],
    ]
    param_names = ["orient", "frame"]

    def setup(self, orient, frame):
        N = 10 ** 5
        ncols = 5
        index = date_range("20000101", periods=N, freq="H")
        timedeltas = timedelta_range(start=1, periods=N, freq="s")
        datetimes = date_range(start=1, periods=N, freq="s")
        ints = np.random.randint(100000000, size=N)
        longints = sys.maxsize * np.random.randint(100000000, size=N)
        floats = np.random.randn(N)
        strings = tm.makeStringIndex(N)
        self.df = DataFrame(np.random.randn(N, ncols), index=np.arange(N))
        self.df_date_idx = DataFrame(np.random.randn(N, ncols), index=index)
        self.df_td_int_ts = DataFrame(
            {
                "td_1": timedeltas,
                "td_2": timedeltas,
                "int_1": ints,
                "int_2": ints,
                "ts_1": datetimes,
                "ts_2": datetimes,
            },
            index=index,
        )
        self.df_int_floats = DataFrame(
            {
                "int_1": ints,
                "int_2": ints,
                "int_3": ints,
                "float_1": floats,
                "float_2": floats,
                "float_3": floats,
            },
            index=index,
        )
        self.df_int_float_str = DataFrame(
            {
                "int_1": ints,
                "int_2": ints,
                "float_1": floats,
                "float_2": floats,
                "str_1": strings,
                "str_2": strings,
            },
            index=index,
        )

        self.df_longint_float_str = DataFrame(
            {
                "longint_1": longints,
                "longint_2": longints,
                "float_1": floats,
                "float_2": floats,
                "str_1": strings,
                "str_2": strings,
            },
            index=index,
        )

    def time_to_json(self, orient, frame):
        getattr(self, frame).to_json(self.fname, orient=orient)

    def peakmem_to_json(self, orient, frame):
        getattr(self, frame).to_json(self.fname, orient=orient)

    def time_to_json_wide(self, orient, frame):
        base_df = getattr(self, frame).copy()
        df = concat([base_df.iloc[:100]] * 1000, ignore_index=True, axis=1)
        df.to_json(self.fname, orient=orient)

    def peakmem_to_json_wide(self, orient, frame):
        base_df = getattr(self, frame).copy()
        df = concat([base_df.iloc[:100]] * 1000, ignore_index=True, axis=1)
        df.to_json(self.fname, orient=orient)


class ToJSONISO(BaseIO):
    fname = "__test__.json"
    params = [["split", "columns", "index", "values", "records"]]
    param_names = ["orient"]

    def setup(self, orient):
        N = 10 ** 5
        index = date_range("20000101", periods=N, freq="H")
        timedeltas = timedelta_range(start=1, periods=N, freq="s")
        datetimes = date_range(start=1, periods=N, freq="s")
        self.df = DataFrame(
            {
                "td_1": timedeltas,
                "td_2": timedeltas,
                "ts_1": datetimes,
                "ts_2": datetimes,
            },
            index=index,
        )

    def time_iso_format(self, orient):
        self.df.to_json(orient=orient, date_format="iso")


class ToJSONLines(BaseIO):

    fname = "__test__.json"

    def setup(self):
        N = 10 ** 5
        ncols = 5
        index = date_range("20000101", periods=N, freq="H")
        timedeltas = timedelta_range(start=1, periods=N, freq="s")
        datetimes = date_range(start=1, periods=N, freq="s")
        ints = np.random.randint(100000000, size=N)
        longints = sys.maxsize * np.random.randint(100000000, size=N)
        floats = np.random.randn(N)
        strings = tm.makeStringIndex(N)
        self.df = DataFrame(np.random.randn(N, ncols), index=np.arange(N))
        self.df_date_idx = DataFrame(np.random.randn(N, ncols), index=index)
        self.df_td_int_ts = DataFrame(
            {
                "td_1": timedeltas,
                "td_2": timedeltas,
                "int_1": ints,
                "int_2": ints,
                "ts_1": datetimes,
                "ts_2": datetimes,
            },
            index=index,
        )
        self.df_int_floats = DataFrame(
            {
                "int_1": ints,
                "int_2": ints,
                "int_3": ints,
                "float_1": floats,
                "float_2": floats,
                "float_3": floats,
            },
            index=index,
        )
        self.df_int_float_str = DataFrame(
            {
                "int_1": ints,
                "int_2": ints,
                "float_1": floats,
                "float_2": floats,
                "str_1": strings,
                "str_2": strings,
            },
            index=index,
        )
        self.df_longint_float_str = DataFrame(
            {
                "longint_1": longints,
                "longint_2": longints,
                "float_1": floats,
                "float_2": floats,
                "str_1": strings,
                "str_2": strings,
            },
            index=index,
        )

    def time_floats_with_int_idex_lines(self):
        self.df.to_json(self.fname, orient="records", lines=True)

    def time_floats_with_dt_index_lines(self):
        self.df_date_idx.to_json(self.fname, orient="records", lines=True)

    def time_delta_int_tstamp_lines(self):
        self.df_td_int_ts.to_json(self.fname, orient="records", lines=True)

    def time_float_int_lines(self):
        self.df_int_floats.to_json(self.fname, orient="records", lines=True)

    def time_float_int_str_lines(self):
        self.df_int_float_str.to_json(self.fname, orient="records", lines=True)

    def time_float_longint_str_lines(self):
        self.df_longint_float_str.to_json(self.fname, orient="records", lines=True)


class ToJSONMem:
    def setup_cache(self):
        df = DataFrame([[1]])
        frames = {"int": df, "float": df.astype(float)}

        return frames

    def peakmem_int(self, frames):
        df = frames["int"]
        for _ in range(100_000):
            df.to_json()

    def peakmem_float(self, frames):
        df = frames["float"]
        for _ in range(100_000):
            df.to_json()


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
