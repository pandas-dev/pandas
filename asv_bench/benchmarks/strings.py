import warnings

import numpy as np

from pandas import (
    NA,
    Categorical,
    DataFrame,
    Index,
    Series,
)
from pandas.arrays import StringArray


class Dtypes:
    params = ["str", "string[python]", "string[pyarrow]"]
    param_names = ["dtype"]

    def setup(self, dtype):
        try:
            self.s = Series(
                Index([f"i-{i}" for i in range(10000)], dtype=object), dtype=dtype
            )
        except ImportError:
            raise NotImplementedError


class Construction:
    params = (
        ["series", "frame", "categorical_series"],
        ["str", "string[python]", "string[pyarrow]"],
    )
    param_names = ["pd_type", "dtype"]
    pd_mapping = {"series": Series, "frame": DataFrame, "categorical_series": Series}
    dtype_mapping = {"str": "str", "string[python]": object, "string[pyarrow]": object}

    def setup(self, pd_type, dtype):
        series_arr = np.array(
            [str(i) * 10 for i in range(100_000)], dtype=self.dtype_mapping[dtype]
        )
        if pd_type == "series":
            self.arr = series_arr
        elif pd_type == "frame":
            self.arr = series_arr.reshape((50_000, 2)).copy()
        elif pd_type == "categorical_series":
            # GH37371. Testing construction of string series/frames from ExtensionArrays
            self.arr = Categorical(series_arr)

    def time_construction(self, pd_type, dtype):
        self.pd_mapping[pd_type](self.arr, dtype=dtype)

    def peakmem_construction(self, pd_type, dtype):
        self.pd_mapping[pd_type](self.arr, dtype=dtype)


class Methods(Dtypes):
    def time_center(self, dtype):
        self.s.str.center(100)

    def time_count(self, dtype):
        self.s.str.count("A")

    def time_endswith(self, dtype):
        self.s.str.endswith("A")

    def time_extract(self, dtype):
        with warnings.catch_warnings(record=True):
            self.s.str.extract("(\\w*)A(\\w*)")

    def time_findall(self, dtype):
        self.s.str.findall("[A-Z]+")

    def time_find(self, dtype):
        self.s.str.find("[A-Z]+")

    def time_rfind(self, dtype):
        self.s.str.rfind("[A-Z]+")

    def time_fullmatch(self, dtype):
        self.s.str.fullmatch("A")

    def time_get(self, dtype):
        self.s.str.get(0)

    def time_len(self, dtype):
        self.s.str.len()

    def time_join(self, dtype):
        self.s.str.join(" ")

    def time_match(self, dtype):
        self.s.str.match("A")

    def time_normalize(self, dtype):
        self.s.str.normalize("NFC")

    def time_pad(self, dtype):
        self.s.str.pad(100, side="both")

    def time_partition(self, dtype):
        self.s.str.partition("A")

    def time_rpartition(self, dtype):
        self.s.str.rpartition("A")

    def time_replace(self, dtype):
        self.s.str.replace("A", "\x01\x01")

    def time_translate(self, dtype):
        self.s.str.translate({"A": "\x01\x01"})

    def time_slice(self, dtype):
        self.s.str.slice(5, 15, 2)

    def time_startswith(self, dtype):
        self.s.str.startswith("A")

    def time_strip(self, dtype):
        self.s.str.strip("A")

    def time_rstrip(self, dtype):
        self.s.str.rstrip("A")

    def time_lstrip(self, dtype):
        self.s.str.lstrip("A")

    def time_title(self, dtype):
        self.s.str.title()

    def time_upper(self, dtype):
        self.s.str.upper()

    def time_lower(self, dtype):
        self.s.str.lower()

    def time_wrap(self, dtype):
        self.s.str.wrap(10)

    def time_zfill(self, dtype):
        self.s.str.zfill(10)

    def time_isalnum(self, dtype):
        self.s.str.isalnum()

    def time_isalpha(self, dtype):
        self.s.str.isalpha()

    def time_isdecimal(self, dtype):
        self.s.str.isdecimal()

    def time_isdigit(self, dtype):
        self.s.str.isdigit()

    def time_islower(self, dtype):
        self.s.str.islower()

    def time_isnumeric(self, dtype):
        self.s.str.isnumeric()

    def time_isspace(self, dtype):
        self.s.str.isspace()

    def time_istitle(self, dtype):
        self.s.str.istitle()

    def time_isupper(self, dtype):
        self.s.str.isupper()


class Repeat:
    params = ["int", "array"]
    param_names = ["repeats"]

    def setup(self, repeats):
        N = 10**5
        self.s = Series(Index([f"i-{i}" for i in range(N)], dtype=object))
        repeat = {"int": 1, "array": np.random.randint(1, 3, N)}
        self.values = repeat[repeats]

    def time_repeat(self, repeats):
        self.s.str.repeat(self.values)


class Cat:
    params = ([0, 3], [None, ","], [None, "-"], [0.0, 0.001, 0.15])
    param_names = ["other_cols", "sep", "na_rep", "na_frac"]

    def setup(self, other_cols, sep, na_rep, na_frac):
        N = 10**5
        mask_gen = lambda: np.random.choice([True, False], N, p=[1 - na_frac, na_frac])
        self.s = Series(Index([f"i-{i}" for i in range(N)], dtype=object)).where(
            mask_gen()
        )
        if other_cols == 0:
            # str.cat self-concatenates only for others=None
            self.others = None
        else:
            self.others = DataFrame(
                {
                    i: Index([f"i-{i}" for i in range(N)], dtype=object).where(
                        mask_gen()
                    )
                    for i in range(other_cols)
                }
            )

    def time_cat(self, other_cols, sep, na_rep, na_frac):
        # before the concatenation (one caller + other_cols columns), the total
        # expected fraction of rows containing any NaN is:
        # reduce(lambda t, _: t + (1 - t) * na_frac, range(other_cols + 1), 0)
        # for other_cols=3 and na_frac=0.15, this works out to ~48%
        self.s.str.cat(others=self.others, sep=sep, na_rep=na_rep)


class Contains(Dtypes):
    params = (Dtypes.params, [True, False])
    param_names = ["dtype", "regex"]

    def setup(self, dtype, regex):
        super().setup(dtype)

    def time_contains(self, dtype, regex):
        self.s.str.contains("A", regex=regex)


class Split(Dtypes):
    params = (Dtypes.params, [True, False])
    param_names = ["dtype", "expand"]

    def setup(self, dtype, expand):
        super().setup(dtype)
        self.s = self.s.str.join("--")

    def time_split(self, dtype, expand):
        self.s.str.split("--", expand=expand)

    def time_rsplit(self, dtype, expand):
        self.s.str.rsplit("--", expand=expand)


class Extract(Dtypes):
    params = (Dtypes.params, [True, False])
    param_names = ["dtype", "expand"]

    def setup(self, dtype, expand):
        super().setup(dtype)

    def time_extract_single_group(self, dtype, expand):
        with warnings.catch_warnings(record=True):
            self.s.str.extract("(\\w*)A", expand=expand)


class Dummies(Dtypes):
    def setup(self, dtype):
        super().setup(dtype)
        N = len(self.s) // 5
        self.s = self.s[:N].str.join("|")

    def time_get_dummies(self, dtype):
        self.s.str.get_dummies("|")


class Encode:
    def setup(self):
        self.ser = Series(Index([f"i-{i}" for i in range(10_000)], dtype=object))

    def time_encode_decode(self):
        self.ser.str.encode("utf-8").str.decode("utf-8")


class Slice:
    def setup(self):
        self.s = Series(["abcdefg", np.nan] * 500000)

    def time_vector_slice(self):
        # GH 2602
        self.s.str[:5]


class Iter(Dtypes):
    def time_iter(self, dtype):
        for i in self.s:
            pass


class StringArrayConstruction:
    def setup(self):
        self.series_arr = np.array([str(i) * 10 for i in range(10**5)], dtype=object)
        self.series_arr_nan = np.concatenate([self.series_arr, np.array([NA] * 1000)])

    def time_string_array_construction(self):
        StringArray(self.series_arr)

    def time_string_array_with_nan_construction(self):
        StringArray(self.series_arr_nan)

    def peakmem_stringarray_construction(self):
        StringArray(self.series_arr)
