from datetime import datetime, timedelta
from io import StringIO
import warnings

import numpy as np
import pytest

from pandas import (
    Categorical,
    DataFrame,
    Series,
    date_range,
    option_context,
    period_range,
)
import pandas._testing as tm

import pandas.io.formats.format as fmt

# Segregated collection of methods that require the BlockManager internal data
# structure


class TestDataFrameReprInfoEtc:
    def test_repr_empty(self):
        # empty
        repr(DataFrame())

        # empty with index
        frame = DataFrame(index=np.arange(1000))
        repr(frame)

    def test_repr_mixed(self, float_string_frame):
        buf = StringIO()

        # mixed
        repr(float_string_frame)
        float_string_frame.info(verbose=False, buf=buf)

    @pytest.mark.slow
    def test_repr_mixed_big(self):
        # big mixed
        biggie = DataFrame(
            {"A": np.random.randn(200), "B": tm.makeStringIndex(200)}, index=range(200)
        )
        biggie.loc[:20, "A"] = np.nan
        biggie.loc[:20, "B"] = np.nan

        repr(biggie)

    def test_repr(self, float_frame):
        buf = StringIO()

        # small one
        repr(float_frame)
        float_frame.info(verbose=False, buf=buf)

        # even smaller
        float_frame.reindex(columns=["A"]).info(verbose=False, buf=buf)
        float_frame.reindex(columns=["A", "B"]).info(verbose=False, buf=buf)

        # exhausting cases in DataFrame.info

        # columns but no index
        no_index = DataFrame(columns=[0, 1, 3])
        repr(no_index)

        # no columns or index
        DataFrame().info(buf=buf)

        df = DataFrame(["a\n\r\tb"], columns=["a\n\r\td"], index=["a\n\r\tf"])
        assert "\t" not in repr(df)
        assert "\r" not in repr(df)
        assert "a\n" not in repr(df)

    def test_repr_dimensions(self):
        df = DataFrame([[1, 2], [3, 4]])
        with option_context("display.show_dimensions", True):
            assert "2 rows x 2 columns" in repr(df)

        with option_context("display.show_dimensions", False):
            assert "2 rows x 2 columns" not in repr(df)

        with option_context("display.show_dimensions", "truncate"):
            assert "2 rows x 2 columns" not in repr(df)

    @pytest.mark.slow
    def test_repr_big(self):
        # big one
        biggie = DataFrame(np.zeros((200, 4)), columns=range(4), index=range(200))
        repr(biggie)

    def test_repr_unsortable(self, float_frame):
        # columns are not sortable

        warn_filters = warnings.filters
        warnings.filterwarnings("ignore", category=FutureWarning, module=".*format")

        unsortable = DataFrame(
            {
                "foo": [1] * 50,
                datetime.today(): [1] * 50,
                "bar": ["bar"] * 50,
                datetime.today() + timedelta(1): ["bar"] * 50,
            },
            index=np.arange(50),
        )
        repr(unsortable)

        fmt.set_option("display.precision", 3, "display.column_space", 10)
        repr(float_frame)

        fmt.set_option("display.max_rows", 10, "display.max_columns", 2)
        repr(float_frame)

        fmt.set_option("display.max_rows", 1000, "display.max_columns", 1000)
        repr(float_frame)

        tm.reset_display_options()

        warnings.filters = warn_filters

    def test_repr_unicode(self):
        uval = "\u03c3\u03c3\u03c3\u03c3"

        df = DataFrame({"A": [uval, uval]})

        result = repr(df)
        ex_top = "      A"
        assert result.split("\n")[0].rstrip() == ex_top

        df = DataFrame({"A": [uval, uval]})
        result = repr(df)
        assert result.split("\n")[0].rstrip() == ex_top

    def test_unicode_string_with_unicode(self):
        df = DataFrame({"A": ["\u05d0"]})
        str(df)

    def test_str_to_bytes_raises(self):
        # GH 26447
        df = DataFrame({"A": ["abc"]})
        msg = "^'str' object cannot be interpreted as an integer$"
        with pytest.raises(TypeError, match=msg):
            bytes(df)

    def test_very_wide_info_repr(self):
        df = DataFrame(np.random.randn(10, 20), columns=tm.rands_array(10, 20))
        repr(df)

    def test_repr_column_name_unicode_truncation_bug(self):
        # #1906
        df = DataFrame(
            {
                "Id": [7117434],
                "StringCol": (
                    "Is it possible to modify drop plot code"
                    "so that the output graph is displayed "
                    "in iphone simulator, Is it possible to "
                    "modify drop plot code so that the "
                    "output graph is \xe2\x80\xa8displayed "
                    "in iphone simulator.Now we are adding "
                    "the CSV file externally. I want to Call "
                    "the File through the code.."
                ),
            }
        )

        with option_context("display.max_columns", 20):
            assert "StringCol" in repr(df)

    def test_latex_repr(self):
        result = r"""\begin{tabular}{llll}
\toprule
{} &         0 &  1 &  2 \\
\midrule
0 &  $\alpha$ &  b &  c \\
1 &         1 &  2 &  3 \\
\bottomrule
\end{tabular}
"""
        with option_context("display.latex.escape", False, "display.latex.repr", True):
            df = DataFrame([[r"$\alpha$", "b", "c"], [1, 2, 3]])
            assert result == df._repr_latex_()

        # GH 12182
        assert df._repr_latex_() is None

    def test_repr_categorical_dates_periods(self):
        # normal DataFrame
        dt = date_range("2011-01-01 09:00", freq="H", periods=5, tz="US/Eastern")
        p = period_range("2011-01", freq="M", periods=5)
        df = DataFrame({"dt": dt, "p": p})
        exp = """                         dt        p
0 2011-01-01 09:00:00-05:00  2011-01
1 2011-01-01 10:00:00-05:00  2011-02
2 2011-01-01 11:00:00-05:00  2011-03
3 2011-01-01 12:00:00-05:00  2011-04
4 2011-01-01 13:00:00-05:00  2011-05"""

        assert repr(df) == exp

        df2 = DataFrame({"dt": Categorical(dt), "p": Categorical(p)})
        assert repr(df2) == exp

    @pytest.mark.parametrize("arg", [np.datetime64, np.timedelta64])
    @pytest.mark.parametrize(
        "box, expected",
        [[Series, "0    NaT\ndtype: object"], [DataFrame, "     0\n0  NaT"]],
    )
    def test_repr_np_nat_with_object(self, arg, box, expected):
        # GH 25445
        result = repr(box([arg("NaT")], dtype=object))
        assert result == expected

    def test_frame_datetime64_pre1900_repr(self):
        df = DataFrame({"year": date_range("1/1/1700", periods=50, freq="A-DEC")})
        # it works!
        repr(df)
