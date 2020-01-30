from pandas import DataFrame
import pandas._testing as tm


def test_sanitize_column_names():
    df = DataFrame(
        {
            "A  ": (1, 1),
            "b": (1, 1),
            "  485 5(468a)44 44   4 ?  $@e3   *   C cc    c D  ": (1, 1),
        }
    )
    df.sanitize_column_names()
    expected_df = DataFrame(
        {"a": (1, 1), "b": (1, 1), "485_5468a44_44_4_e3_c_cc_c_d": (1, 1)}
    )
    tm.assert_frame_equal(df, expected_df)
