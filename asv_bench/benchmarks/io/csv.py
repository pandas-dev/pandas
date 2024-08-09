from io import (
    BytesIO,
    StringIO,
)
import random
import string

import numpy as np

from pandas import (
    Categorical,
    DataFrame,
    Index,
    concat,
    date_range,
    period_range,
    read_csv,
    to_datetime,
)

from ..pandas_vb_common import BaseIO


class ToCSV(BaseIO):
    fname = "__test__.csv"
    params = ["wide", "long", "mixed"]
    param_names = ["kind"]

    def setup(self, kind):
        wide_frame = DataFrame(np.random.randn(3000, 30))
        long_frame = DataFrame(
            {
                "A": np.arange(50000),
                "B": np.arange(50000) + 1.0,
                "C": np.arange(50000) + 2.0,
                "D": np.arange(50000) + 3.0,
            }
        )
        mixed_frame = DataFrame(
            {
                "float": np.random.randn(5000),
                "int": np.random.randn(5000).astype(int),
                "bool": (np.arange(5000) % 2) == 0,
                "datetime": date_range("2001", freq="s", periods=5000),
                "object": ["foo"] * 5000,
            }
        )
        mixed_frame.loc[30:500, "float"] = np.nan
        data = {"wide": wide_frame, "long": long_frame, "mixed": mixed_frame}
        self.df = data[kind]

    def time_frame(self, kind):
        self.df.to_csv(self.fname)


class ToCSVMultiIndexUnusedLevels(BaseIO):
    fname = "__test__.csv"

    def setup(self):
        df = DataFrame({"a": np.random.randn(100_000), "b": 1, "c": 1})
        self.df = df.set_index(["a", "b"])
        self.df_unused_levels = self.df.iloc[:10_000]
        self.df_single_index = df.set_index(["a"]).iloc[:10_000]

    def time_full_frame(self):
        self.df.to_csv(self.fname)

    def time_sliced_frame(self):
        self.df_unused_levels.to_csv(self.fname)

    def time_single_index_frame(self):
        self.df_single_index.to_csv(self.fname)


class ToCSVDatetime(BaseIO):
    fname = "__test__.csv"

    def setup(self):
        rng = date_range("1/1/2000", periods=1000)
        self.data = DataFrame(rng, index=rng)

    def time_frame_date_formatting_default(self):
        self.data.to_csv(self.fname)

    def time_frame_date_formatting_default_explicit(self):
        self.data.to_csv(self.fname, date_format="%Y-%m-%d")

    def time_frame_date_formatting_custom(self):
        self.data.to_csv(self.fname, date_format="%Y%m%d")


class ToCSVDatetimeIndex(BaseIO):
    fname = "__test__.csv"

    def setup(self):
        rng = date_range("2000", periods=100_000, freq="s")
        self.data = DataFrame({"a": 1}, index=rng)

    def time_frame_date_formatting_index_default(self):
        self.data.to_csv(self.fname)

    def time_frame_date_formatting_index_default_explicit(self):
        self.data.to_csv(self.fname, date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_date_formatting_index_custom(self):
        self.data.to_csv(self.fname, date_format="%Y-%m-%d__%H:%M:%S")


class ToCSVPeriod(BaseIO):
    fname = "__test__.csv"

    params = ([1000, 10000], ["D", "h"])
    param_names = ["nobs", "freq"]

    def setup(self, nobs, freq):
        rng = period_range(start="2000-01-01", periods=nobs, freq=freq)
        self.data = DataFrame(rng)
        if freq == "D":
            self.default_fmt = "%Y-%m-%d"
        elif freq == "h":
            self.default_fmt = "%Y-%m-%d %H:00"

    def time_frame_period_formatting_default(self, nobs, freq):
        self.data.to_csv(self.fname)

    def time_frame_period_formatting_default_explicit(self, nobs, freq):
        self.data.to_csv(self.fname, date_format=self.default_fmt)

    def time_frame_period_formatting_custom(self, nobs, freq):
        # Nb: `date_format` is not actually taken into account here today, so the
        # performance is currently identical to `time_frame_period_formatting_default`
        # above. This timer is therefore expected to degrade when GH#51621 is fixed.
        # (Remove this comment when GH#51621 is fixed.)
        self.data.to_csv(self.fname, date_format="%Y-%m-%d___%H:%M:%S")


class ToCSVPeriodIndex(BaseIO):
    fname = "__test__.csv"

    params = ([1000, 10000], ["D", "h"])
    param_names = ["nobs", "freq"]

    def setup(self, nobs, freq):
        rng = period_range(start="2000-01-01", periods=nobs, freq=freq)
        self.data = DataFrame({"a": 1}, index=rng)
        if freq == "D":
            self.default_fmt = "%Y-%m-%d"
        elif freq == "h":
            self.default_fmt = "%Y-%m-%d %H:00"

    def time_frame_period_formatting_index_default(self, nobs, freq):
        self.data.to_csv(self.fname)

    def time_frame_period_formatting_index_default_explicit(self, nobs, freq):
        self.data.to_csv(self.fname, date_format=self.default_fmt)

    def time_frame_period_formatting_index_custom(self, nobs, freq):
        self.data.to_csv(self.fname, date_format="%Y-%m-%d___%H:%M:%S")


class ToCSVDatetimeBig(BaseIO):
    fname = "__test__.csv"
    timeout = 1500
    params = [1000, 10000, 100000]
    param_names = ["nobs"]

    def setup(self, nobs):
        d = "2018-11-29"
        dt = "2018-11-26 11:18:27.0"
        self.data = DataFrame(
            {
                "dt": [np.datetime64(dt)] * nobs,
                "d": [np.datetime64(d)] * nobs,
                "r": [np.random.uniform()] * nobs,
            }
        )

    def time_frame_formatting_default(self, nobs):
        self.data.to_csv(self.fname)

    def time_frame_date_formatting_custom(self, nobs):
        self.data.to_csv(self.fname, date_format="%Y%m%d__%H%M%S")


class ToCSVIndexes(BaseIO):
    fname = "__test__.csv"

    @staticmethod
    def _create_df(rows, cols):
        index_cols = {
            "index1": np.random.randint(0, rows, rows),
            "index2": np.full(rows, 1, dtype=int),
            "index3": np.full(rows, 1, dtype=int),
        }
        data_cols = {
            f"col{i}": np.random.uniform(0, 100000.0, rows) for i in range(cols)
        }
        df = DataFrame({**index_cols, **data_cols})
        return df

    def setup(self):
        ROWS = 100000
        COLS = 5
        # For tests using .head(), create an initial dataframe with this many times
        # more rows
        HEAD_ROW_MULTIPLIER = 10

        self.df_standard_index = self._create_df(ROWS, COLS)

        self.df_custom_index_then_head = (
            self._create_df(ROWS * HEAD_ROW_MULTIPLIER, COLS)
            .set_index(["index1", "index2", "index3"])
            .head(ROWS)
        )

        self.df_head_then_custom_index = (
            self._create_df(ROWS * HEAD_ROW_MULTIPLIER, COLS)
            .head(ROWS)
            .set_index(["index1", "index2", "index3"])
        )

    def time_standard_index(self):
        self.df_standard_index.to_csv(self.fname)

    def time_multiindex(self):
        self.df_head_then_custom_index.to_csv(self.fname)

    def time_head_of_multiindex(self):
        self.df_custom_index_then_head.to_csv(self.fname)


class StringIORewind:
    def data(self, stringio_object):
        stringio_object.seek(0)
        return stringio_object


class ReadCSVDInferDatetimeFormat(StringIORewind):
    params = [None, "custom", "iso8601", "ymd"]
    param_names = ["format"]

    def setup(self, format):
        rng = date_range("1/1/2000", periods=1000)
        formats = {
            None: None,
            "custom": "%m/%d/%Y %H:%M:%S.%f",
            "iso8601": "%Y-%m-%d %H:%M:%S",
            "ymd": "%Y%m%d",
        }
        dt_format = formats[format]
        self.StringIO_input = StringIO("\n".join(rng.strftime(dt_format).tolist()))

    def time_read_csv(self, format):
        read_csv(
            self.data(self.StringIO_input),
            header=None,
            names=["foo"],
            parse_dates=["foo"],
        )


class ReadCSVConcatDatetime(StringIORewind):
    iso8601 = "%Y-%m-%d %H:%M:%S"

    def setup(self):
        rng = date_range("1/1/2000", periods=50000, freq="s")
        self.StringIO_input = StringIO("\n".join(rng.strftime(self.iso8601).tolist()))

    def time_read_csv(self):
        read_csv(
            self.data(self.StringIO_input),
            header=None,
            names=["foo"],
            parse_dates=["foo"],
        )


class ReadCSVConcatDatetimeBadDateValue(StringIORewind):
    params = (["nan", "0", ""],)
    param_names = ["bad_date_value"]

    def setup(self, bad_date_value):
        self.StringIO_input = StringIO((f"{bad_date_value},\n") * 50000)

    def time_read_csv(self, bad_date_value):
        read_csv(
            self.data(self.StringIO_input),
            header=None,
            names=["foo", "bar"],
            parse_dates=["foo"],
        )


class ReadCSVSkipRows(BaseIO):
    fname = "__test__.csv"
    params = ([None, 10000], ["c", "python", "pyarrow"])
    param_names = ["skiprows", "engine"]

    def setup(self, skiprows, engine):
        N = 20000
        index = Index([f"i-{i}" for i in range(N)], dtype=object)
        df = DataFrame(
            {
                "float1": np.random.randn(N),
                "float2": np.random.randn(N),
                "string1": ["foo"] * N,
                "bool1": [True] * N,
                "int1": np.random.randint(0, N, size=N),
            },
            index=index,
        )
        df.to_csv(self.fname)

    def time_skipprows(self, skiprows, engine):
        read_csv(self.fname, skiprows=skiprows, engine=engine)


class ReadUint64Integers(StringIORewind):
    def setup(self):
        self.na_values = [2**63 + 500]
        arr = np.arange(10000).astype("uint64") + 2**63
        self.data1 = StringIO("\n".join(arr.astype(str).tolist()))
        arr = arr.astype(object)
        arr[500] = -1
        self.data2 = StringIO("\n".join(arr.astype(str).tolist()))

    def time_read_uint64(self):
        read_csv(self.data(self.data1), header=None, names=["foo"])

    def time_read_uint64_neg_values(self):
        read_csv(self.data(self.data2), header=None, names=["foo"])

    def time_read_uint64_na_values(self):
        read_csv(
            self.data(self.data1), header=None, names=["foo"], na_values=self.na_values
        )


class ReadCSVThousands(BaseIO):
    fname = "__test__.csv"
    params = ([",", "|"], [None, ","], ["c", "python"])
    param_names = ["sep", "thousands", "engine"]

    def setup(self, sep, thousands, engine):
        N = 10000
        K = 8
        data = np.random.randn(N, K) * np.random.randint(100, 10000, (N, K))
        df = DataFrame(data)
        if thousands is not None:
            fmt = f":{thousands}"
            fmt = "{" + fmt + "}"
            df = df.map(lambda x: fmt.format(x))
        df.to_csv(self.fname, sep=sep)

    def time_thousands(self, sep, thousands, engine):
        read_csv(self.fname, sep=sep, thousands=thousands, engine=engine)


class ReadCSVComment(StringIORewind):
    params = ["c", "python"]
    param_names = ["engine"]

    def setup(self, engine):
        data = ["A,B,C"] + (["1,2,3 # comment"] * 100000)
        self.StringIO_input = StringIO("\n".join(data))

    def time_comment(self, engine):
        read_csv(
            self.data(self.StringIO_input), comment="#", header=None, names=list("abc")
        )


class ReadCSVFloatPrecision(StringIORewind):
    params = ([",", ";"], [".", "_"], [None, "high", "round_trip"])
    param_names = ["sep", "decimal", "float_precision"]

    def setup(self, sep, decimal, float_precision):
        floats = [
            "".join([random.choice(string.digits) for _ in range(28)])
            for _ in range(15)
        ]
        rows = sep.join([f"0{decimal}{{}}"] * 3) + "\n"
        data = rows * 5
        data = data.format(*floats) * 200  # 1000 x 3 strings csv
        self.StringIO_input = StringIO(data)

    def time_read_csv(self, sep, decimal, float_precision):
        read_csv(
            self.data(self.StringIO_input),
            sep=sep,
            header=None,
            names=list("abc"),
            float_precision=float_precision,
        )

    def time_read_csv_python_engine(self, sep, decimal, float_precision):
        read_csv(
            self.data(self.StringIO_input),
            sep=sep,
            header=None,
            engine="python",
            float_precision=None,
            names=list("abc"),
        )


class ReadCSVEngine(StringIORewind):
    params = ["c", "python", "pyarrow"]
    param_names = ["engine"]

    def setup(self, engine):
        data = ["A,B,C,D,E"] + (["1,2,3,4,5"] * 100000)
        self.StringIO_input = StringIO("\n".join(data))
        # simulate reading from file
        self.BytesIO_input = BytesIO(self.StringIO_input.read().encode("utf-8"))

    def time_read_stringcsv(self, engine):
        read_csv(self.data(self.StringIO_input), engine=engine)

    def time_read_bytescsv(self, engine):
        read_csv(self.data(self.BytesIO_input), engine=engine)

    def peakmem_read_csv(self, engine):
        read_csv(self.data(self.BytesIO_input), engine=engine)


class ReadCSVCategorical(BaseIO):
    fname = "__test__.csv"
    params = ["c", "python"]
    param_names = ["engine"]

    def setup(self, engine):
        N = 100000
        group1 = ["aaaaaaaa", "bbbbbbb", "cccccccc", "dddddddd", "eeeeeeee"]
        df = DataFrame(np.random.choice(group1, (N, 3)), columns=list("abc"))
        df.to_csv(self.fname, index=False)

    def time_convert_post(self, engine):
        read_csv(self.fname, engine=engine).apply(Categorical)

    def time_convert_direct(self, engine):
        read_csv(self.fname, engine=engine, dtype="category")


class ReadCSVParseDates(StringIORewind):
    params = ["c", "python"]
    param_names = ["engine"]

    def setup(self, engine):
        data = """{},19:00:00,18:56:00,0.8100,2.8100,7.2000,0.0000,280.0000\n
                  {},20:00:00,19:56:00,0.0100,2.2100,7.2000,0.0000,260.0000\n
                  {},21:00:00,20:56:00,-0.5900,2.2100,5.7000,0.0000,280.0000\n
                  {},21:00:00,21:18:00,-0.9900,2.0100,3.6000,0.0000,270.0000\n
                  {},22:00:00,21:56:00,-0.5900,1.7100,5.1000,0.0000,290.0000\n
               """
        two_cols = ["KORD,19990127"] * 5
        data = data.format(*two_cols)
        self.StringIO_input = StringIO(data)

    def time_baseline(self, engine):
        read_csv(
            self.data(self.StringIO_input),
            engine=engine,
            sep=",",
            header=None,
            parse_dates=[1],
            names=list(string.digits[:9]),
        )


class ReadCSVCachedParseDates(StringIORewind):
    params = ([True, False], ["c", "python"])
    param_names = ["do_cache", "engine"]

    def setup(self, do_cache, engine):
        data = ("\n".join([f"10/{year}" for year in range(2000, 2100)]) + "\n") * 10
        self.StringIO_input = StringIO(data)

    def time_read_csv_cached(self, do_cache, engine):
        try:
            read_csv(
                self.data(self.StringIO_input),
                engine=engine,
                header=None,
                parse_dates=[0],
                cache_dates=do_cache,
            )
        except TypeError:
            # cache_dates is a new keyword in 0.25
            pass


class ReadCSVMemoryGrowth(BaseIO):
    chunksize = 20
    num_rows = 1000
    fname = "__test__.csv"
    params = ["c", "python"]
    param_names = ["engine"]

    def setup(self, engine):
        with open(self.fname, "w", encoding="utf-8") as f:
            for i in range(self.num_rows):
                f.write(f"{i}\n")

    def mem_parser_chunks(self, engine):
        # see gh-24805.
        result = read_csv(self.fname, chunksize=self.chunksize, engine=engine)

        for _ in result:
            pass


class ReadCSVParseSpecialDate(StringIORewind):
    params = (["mY", "mdY", "hm"], ["c", "python"])
    param_names = ["value", "engine"]
    objects = {
        "mY": "01-2019\n10-2019\n02/2000\n",
        "mdY": "12/02/2010\n",
        "hm": "21:34\n",
    }

    def setup(self, value, engine):
        count_elem = 10000
        data = self.objects[value] * count_elem
        self.StringIO_input = StringIO(data)

    def time_read_special_date(self, value, engine):
        read_csv(
            self.data(self.StringIO_input),
            engine=engine,
            sep=",",
            header=None,
            names=["Date"],
            parse_dates=["Date"],
        )


class ReadCSVMemMapUTF8:
    fname = "__test__.csv"
    number = 5

    def setup(self):
        lines = []
        line_length = 128
        start_char = " "
        end_char = "\U00010080"
        # This for loop creates a list of 128-char strings
        # consisting of consecutive Unicode chars
        for lnum in range(ord(start_char), ord(end_char), line_length):
            line = "".join([chr(c) for c in range(lnum, lnum + 0x80)]) + "\n"
            try:
                line.encode("utf-8")
            except UnicodeEncodeError:
                # Some 16-bit words are not valid Unicode chars and must be skipped
                continue
            lines.append(line)
        df = DataFrame(lines)
        df = concat([df for n in range(100)], ignore_index=True)
        df.to_csv(self.fname, index=False, header=False, encoding="utf-8")

    def time_read_memmapped_utf8(self):
        read_csv(self.fname, header=None, memory_map=True, encoding="utf-8", engine="c")


class ParseDateComparison(StringIORewind):
    params = ([False, True],)
    param_names = ["cache_dates"]

    def setup(self, cache_dates):
        count_elem = 10000
        data = "12-02-2010\n" * count_elem
        self.StringIO_input = StringIO(data)

    def time_read_csv_dayfirst(self, cache_dates):
        try:
            read_csv(
                self.data(self.StringIO_input),
                sep=",",
                header=None,
                names=["Date"],
                parse_dates=["Date"],
                cache_dates=cache_dates,
                dayfirst=True,
            )
        except TypeError:
            # cache_dates is a new keyword in 0.25
            pass

    def time_to_datetime_dayfirst(self, cache_dates):
        df = read_csv(
            self.data(self.StringIO_input), dtype={"date": str}, names=["date"]
        )
        to_datetime(df["date"], cache=cache_dates, dayfirst=True)

    def time_to_datetime_format_DD_MM_YYYY(self, cache_dates):
        df = read_csv(
            self.data(self.StringIO_input), dtype={"date": str}, names=["date"]
        )
        to_datetime(df["date"], cache=cache_dates, format="%d-%m-%Y")


class ReadCSVIndexCol(StringIORewind):
    def setup(self):
        count_elem = 100_000
        data = "a,b\n" + "1,2\n" * count_elem
        self.StringIO_input = StringIO(data)

    def time_read_csv_index_col(self):
        read_csv(self.StringIO_input, index_col="a")


class ReadCSVDatePyarrowEngine(StringIORewind):
    def setup(self):
        count_elem = 100_000
        data = "a\n" + "2019-12-31\n" * count_elem
        self.StringIO_input = StringIO(data)

    def time_read_csv_index_col(self):
        read_csv(
            self.StringIO_input,
            parse_dates=["a"],
            engine="pyarrow",
            dtype_backend="pyarrow",
        )


class ReadCSVCParserLowMemory:
    # GH 16798
    def setup(self):
        self.csv = StringIO(
            "strings\n" + "\n".join(["x" * (1 << 20) for _ in range(2100)])
        )

    def peakmem_over_2gb_input(self):
        read_csv(self.csv, engine="c", low_memory=False)


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
