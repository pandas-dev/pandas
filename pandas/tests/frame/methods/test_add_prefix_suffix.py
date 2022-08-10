from pandas import Index
import pandas._testing as tm


def test_add_prefix_suffix(float_frame):
    with_prefix = float_frame.add_prefix("foo#")
    expected = Index([f"foo#{c}" for c in float_frame.columns])
    tm.assert_index_equal(with_prefix.columns, expected)

    with_suffix = float_frame.add_suffix("#foo")
    expected = Index([f"{c}#foo" for c in float_frame.columns])
    tm.assert_index_equal(with_suffix.columns, expected)

    with_pct_prefix = float_frame.add_prefix("%")
    expected = Index([f"%{c}" for c in float_frame.columns])
    tm.assert_index_equal(with_pct_prefix.columns, expected)

    with_pct_suffix = float_frame.add_suffix("%")
    expected = Index([f"{c}%" for c in float_frame.columns])
    tm.assert_index_equal(with_pct_suffix.columns, expected)


def test_add_prefix_suffix_copy(float_frame):
    # GH#47934
    ser = float_frame.iloc[0]

    with_prefix = float_frame.add_prefix("foo#", copy=True)
    expected = Index([f"foo#{c}" for c in float_frame.columns])
    tm.assert_index_equal(with_prefix.columns, expected)
    assert not any(
        tm.shares_memory(float_frame.iloc[:, i], with_prefix.iloc[:, i])
        for i in range(float_frame.shape[1])
    )

    ser_with_prefix = ser.add_prefix("foo#", copy=True)
    tm.assert_index_equal(ser_with_prefix.index, expected)
    assert not tm.shares_memory(ser_with_prefix, ser)

    with_prefix = float_frame.add_prefix("foo#", copy=False)
    expected = Index([f"foo#{c}" for c in float_frame.columns])
    tm.assert_index_equal(with_prefix.columns, expected)
    assert all(
        tm.shares_memory(float_frame.iloc[:, i], with_prefix.iloc[:, i])
        for i in range(float_frame.shape[1])
    )

    ser_with_prefix = ser.add_prefix("foo#", copy=False)
    tm.assert_index_equal(ser_with_prefix.index, expected)
    assert tm.shares_memory(ser_with_prefix, ser)

    with_suffix = float_frame.add_suffix("#foo", copy=True)
    expected = Index([f"{c}#foo" for c in float_frame.columns])
    tm.assert_index_equal(with_suffix.columns, expected)
    assert not any(
        tm.shares_memory(float_frame.iloc[:, i], with_suffix.iloc[:, i])
        for i in range(float_frame.shape[1])
    )

    ser_with_suffix = ser.add_suffix("#foo", copy=True)
    tm.assert_index_equal(ser_with_suffix.index, expected)
    assert not tm.shares_memory(ser_with_suffix, ser)

    with_suffix = float_frame.add_suffix("#foo", copy=False)
    expected = Index([f"{c}#foo" for c in float_frame.columns])
    tm.assert_index_equal(with_suffix.columns, expected)
    assert all(
        tm.shares_memory(float_frame.iloc[:, i], with_suffix.iloc[:, i])
        for i in range(float_frame.shape[1])
    )

    ser_with_suffix = ser.add_suffix("#foo", copy=False)
    tm.assert_index_equal(ser_with_suffix.index, expected)
    assert tm.shares_memory(ser_with_suffix, ser)
