import numpy as np
import pytest

from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    date_range,
)
import pandas._testing as tm


def test_doc_string():
    df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    df.ewm(com=0.5).mean()


def test_constructor(frame_or_series):
    c = frame_or_series(range(5)).ewm

    # valid
    c(com=0.5)
    c(span=1.5)
    c(alpha=0.5)
    c(halflife=0.75)
    c(com=0.5, span=None)
    c(alpha=0.5, com=None)
    c(halflife=0.75, alpha=None)

    # not valid: mutually exclusive
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        c(com=0.5, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        c(span=1.5, halflife=0.75)
    with pytest.raises(ValueError, match=msg):
        c(alpha=0.5, span=1.5)

    # not valid: com < 0
    msg = "comass must satisfy: comass >= 0"
    with pytest.raises(ValueError, match=msg):
        c(com=-0.5)

    # not valid: span < 1
    msg = "span must satisfy: span >= 1"
    with pytest.raises(ValueError, match=msg):
        c(span=0.5)

    # not valid: halflife <= 0
    msg = "halflife must satisfy: halflife > 0"
    with pytest.raises(ValueError, match=msg):
        c(halflife=0)

    # not valid: alpha <= 0 or alpha > 1
    msg = "alpha must satisfy: 0 < alpha <= 1"
    for alpha in (-0.5, 1.5):
        with pytest.raises(ValueError, match=msg):
            c(alpha=alpha)


def test_ewma_times_not_datetime_type():
    msg = r"times must be datetime64 dtype."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(times=np.arange(5))


def test_ewma_times_not_same_length():
    msg = "times must be the same length as the object."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(times=np.arange(4).astype("datetime64[ns]"))


def test_ewma_halflife_not_correct_type():
    msg = "halflife must be a timedelta convertible object"
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(halflife=1, times=np.arange(5).astype("datetime64[ns]"))


def test_ewma_halflife_without_times(halflife_with_times):
    msg = "halflife can only be a timedelta convertible argument if times is not None."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(halflife=halflife_with_times)


@pytest.mark.parametrize(
    "times",
    [
        np.arange(10).astype("datetime64[D]").astype("datetime64[ns]"),
        date_range("2000", freq="D", periods=10),
        date_range("2000", freq="D", periods=10).tz_localize("UTC"),
    ],
)
@pytest.mark.parametrize("min_periods", [0, 2])
def test_ewma_with_times_equal_spacing(
    halflife_with_times, times, min_periods, adjust, ignore_na
):
    # GH#66523 equally spaced times must match no-times for every
    # adjust/ignore_na combination, including on NaN data.
    halflife = halflife_with_times
    data = np.arange(10.0)
    data[::2] = np.nan
    df = DataFrame({"A": data})
    result = df.ewm(
        halflife=halflife,
        min_periods=min_periods,
        times=times,
        adjust=adjust,
        ignore_na=ignore_na,
    ).mean()
    expected = df.ewm(
        halflife=1.0, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
    ).mean()
    tm.assert_frame_equal(result, expected)


def test_ewma_with_times_variable_spacing(tz_aware_fixture, unit, adjust):
    # GH 54328
    tz = tz_aware_fixture
    halflife = "23 days"
    times = (
        DatetimeIndex(["2020-01-01", "2020-01-10T00:04:05", "2020-02-23T05:00:23"])
        .tz_localize(tz)
        .as_unit(unit)
    )
    data = np.arange(3)
    df = DataFrame(data)
    result = df.ewm(halflife=halflife, times=times, adjust=adjust).mean()
    if adjust:
        expected = DataFrame([0.0, 0.5674161888241773, 1.545239952073459])
    else:
        expected = DataFrame([0.0, 0.23762518642226227, 1.534926369128742])
    tm.assert_frame_equal(result, expected)


def test_ewm_with_nat_raises(halflife_with_times):
    # GH#38535
    ser = Series(range(1))
    times = DatetimeIndex(["NaT"])
    with pytest.raises(ValueError, match="Cannot convert NaT values to integer"):
        ser.ewm(com=0.1, halflife=halflife_with_times, times=times)


def test_ewm_with_times_getitem(halflife_with_times):
    # GH 40164
    halflife = halflife_with_times
    data = np.arange(10.0)
    data[::2] = np.nan
    times = date_range("2000", freq="D", periods=10)
    df = DataFrame({"A": data, "B": data})
    result = df.ewm(halflife=halflife, times=times)["A"].mean()
    expected = df.ewm(halflife=1.0)["A"].mean()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("arg", ["com", "halflife", "span", "alpha"])
def test_ewm_getitem_attributes_retained(arg, adjust, ignore_na):
    # GH 40164
    kwargs = {arg: 1, "adjust": adjust, "ignore_na": ignore_na}
    ewm = DataFrame({"A": range(1), "B": range(1)}).ewm(**kwargs)
    expected = {attr: getattr(ewm, attr) for attr in ewm._attributes}
    ewm_slice = ewm["A"]
    result = {attr: getattr(ewm, attr) for attr in ewm_slice._attributes}
    assert result == expected


def test_ewma_times_adjust_false_with_disallowed_com():
    # GH 54328
    with pytest.raises(
        NotImplementedError,
        match=(
            "None of com, span, or alpha can be specified "
            "if times is provided and adjust=False"
        ),
    ):
        Series(range(1)).ewm(
            0.1,
            adjust=False,
            times=date_range("2000", freq="D", periods=1),
            halflife="1D",
        )


def test_ewma_times_adjust_false_with_disallowed_alpha():
    # GH 54328
    with pytest.raises(
        NotImplementedError,
        match=(
            "None of com, span, or alpha can be specified "
            "if times is provided and adjust=False"
        ),
    ):
        Series(range(1)).ewm(
            0.1,
            adjust=False,
            times=date_range("2000", freq="D", periods=1),
            alpha=0.5,
            halflife="1D",
        )


def test_ewma_times_adjust_false_with_disallowed_span():
    # GH 54328
    with pytest.raises(
        NotImplementedError,
        match=(
            "None of com, span, or alpha can be specified "
            "if times is provided and adjust=False"
        ),
    ):
        Series(range(1)).ewm(
            0.1,
            adjust=False,
            times=date_range("2000", freq="D", periods=1),
            span=10,
            halflife="1D",
        )


def test_times_string_col_raises():
    # GH 43265
    df = DataFrame(
        {"A": np.arange(10.0), "time_col": date_range("2000", freq="D", periods=10)}
    )
    with pytest.raises(ValueError, match="times must be datetime64"):
        df.ewm(halflife="1 day", min_periods=0, times="time_col")


def test_ewm_sum_adjust_false_notimplemented():
    data = Series(range(1)).ewm(com=1, adjust=False)
    with pytest.raises(NotImplementedError, match="sum is not"):
        data.sum()


def test_ewm_aggregate_callable_gh41700(frame_or_series):
    # GH#41700
    obj = frame_or_series(range(5)).ewm(alpha=0.1)
    msg = "aggregate does not support arbitrary callables"
    with pytest.raises(NotImplementedError, match=msg):
        obj.agg(lambda x: x)
    with pytest.raises(NotImplementedError, match=msg):
        obj.agg(np.sum)
    # supported string aggregations still work
    tm.assert_equal(obj.agg("mean"), obj.mean())


@pytest.mark.parametrize("method", ["sum", "std", "var", "cov", "corr"])
def test_times_only_mean_implemented(frame_or_series, method):
    # GH 51695
    halflife = "1 day"
    times = date_range("2000", freq="D", periods=10)
    ewm = frame_or_series(range(10)).ewm(halflife=halflife, times=times)
    with pytest.raises(
        NotImplementedError, match=f"{method} is not implemented with times"
    ):
        getattr(ewm, method)()


@pytest.mark.parametrize(
    "expected_data, ignore",
    [[[10.0, 5.0, 2.5, 11.25], False], [[10.0, 5.0, 5.0, 12.5], True]],
)
def test_ewm_sum(expected_data, ignore):
    # xref from Numbagg tests
    # https://github.com/numbagg/numbagg/blob/v0.2.1/numbagg/test/test_moving.py#L50
    data = Series([10, 0, np.nan, 10])
    result = data.ewm(alpha=0.5, ignore_na=ignore).sum()
    expected = Series(expected_data)
    tm.assert_series_equal(result, expected)


def test_ewma_adjust():
    vals = Series(np.zeros(1000))
    vals[5] = 1
    result = vals.ewm(span=100, adjust=False).mean().sum()
    assert np.abs(result - 1) < 1e-2


def test_ewma_cases(adjust, ignore_na):
    # try adjust/ignore_na args matrix

    s = Series([1.0, 2.0, 4.0, 8.0])

    if adjust:
        expected = Series([1.0, 1.6, 2.736842, 4.923077])
    else:
        expected = Series([1.0, 1.333333, 2.222222, 4.148148])

    result = s.ewm(com=2.0, adjust=adjust, ignore_na=ignore_na).mean()
    tm.assert_series_equal(result, expected)


def test_ewma_nan_handling():
    s = Series([1.0] + [np.nan] * 5 + [1.0])
    result = s.ewm(com=5).mean()
    tm.assert_series_equal(result, Series([1.0] * len(s)))

    s = Series([np.nan] * 2 + [1.0] + [np.nan] * 2 + [1.0])
    result = s.ewm(com=5).mean()
    tm.assert_series_equal(result, Series([np.nan] * 2 + [1.0] * 4))


@pytest.mark.parametrize(
    "s, adjust, ignore_na, w",
    [
        (
            [np.nan, 1.0, 101.0],
            True,
            False,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))), 1.0],
        ),
        (
            [np.nan, 1.0, 101.0],
            True,
            True,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))), 1.0],
        ),
        (
            [np.nan, 1.0, 101.0],
            False,
            False,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))), (1.0 / (1.0 + 2.0))],
        ),
        (
            [np.nan, 1.0, 101.0],
            False,
            True,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))), (1.0 / (1.0 + 2.0))],
        ),
        (
            [1.0, np.nan, 101.0],
            True,
            False,
            [(1.0 - (1.0 / (1.0 + 2.0))) ** 2, np.nan, 1.0],
        ),
        (
            [1.0, np.nan, 101.0],
            True,
            True,
            [(1.0 - (1.0 / (1.0 + 2.0))), np.nan, 1.0],
        ),
        (
            [1.0, np.nan, 101.0],
            False,
            False,
            [(1.0 - (1.0 / (1.0 + 2.0))) ** 2, np.nan, (1.0 / (1.0 + 2.0))],
        ),
        (
            [1.0, np.nan, 101.0],
            False,
            True,
            [(1.0 - (1.0 / (1.0 + 2.0))), np.nan, (1.0 / (1.0 + 2.0))],
        ),
        (
            [np.nan, 1.0, np.nan, np.nan, 101.0, np.nan],
            True,
            False,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))) ** 3, np.nan, np.nan, 1.0, np.nan],
        ),
        (
            [np.nan, 1.0, np.nan, np.nan, 101.0, np.nan],
            True,
            True,
            [np.nan, (1.0 - (1.0 / (1.0 + 2.0))), np.nan, np.nan, 1.0, np.nan],
        ),
        (
            [np.nan, 1.0, np.nan, np.nan, 101.0, np.nan],
            False,
            False,
            [
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))) ** 3,
                np.nan,
                np.nan,
                (1.0 / (1.0 + 2.0)),
                np.nan,
            ],
        ),
        (
            [np.nan, 1.0, np.nan, np.nan, 101.0, np.nan],
            False,
            True,
            [
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))),
                np.nan,
                np.nan,
                (1.0 / (1.0 + 2.0)),
                np.nan,
            ],
        ),
        (
            [1.0, np.nan, 101.0, 50.0],
            True,
            False,
            [
                (1.0 - (1.0 / (1.0 + 2.0))) ** 3,
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))),
                1.0,
            ],
        ),
        (
            [1.0, np.nan, 101.0, 50.0],
            True,
            True,
            [
                (1.0 - (1.0 / (1.0 + 2.0))) ** 2,
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))),
                1.0,
            ],
        ),
        (
            [1.0, np.nan, 101.0, 50.0],
            False,
            False,
            [
                (1.0 - (1.0 / (1.0 + 2.0))) ** 3,
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))) * (1.0 / (1.0 + 2.0)),
                (1.0 / (1.0 + 2.0))
                * ((1.0 - (1.0 / (1.0 + 2.0))) ** 2 + (1.0 / (1.0 + 2.0))),
            ],
        ),
        (
            [1.0, np.nan, 101.0, 50.0],
            False,
            True,
            [
                (1.0 - (1.0 / (1.0 + 2.0))) ** 2,
                np.nan,
                (1.0 - (1.0 / (1.0 + 2.0))) * (1.0 / (1.0 + 2.0)),
                (1.0 / (1.0 + 2.0)),
            ],
        ),
    ],
)
def test_ewma_nan_handling_cases(s, adjust, ignore_na, w):
    # GH 7603
    s = Series(s)
    expected = (s.multiply(w).cumsum() / Series(w).cumsum()).ffill()
    result = s.ewm(com=2.0, adjust=adjust, ignore_na=ignore_na).mean()

    tm.assert_series_equal(result, expected)
    if ignore_na is False:
        # check that ignore_na defaults to False
        result = s.ewm(com=2.0, adjust=adjust).mean()
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs", [{"com": 1}, {"span": 3}, {"alpha": 0.5}, {"halflife": 1}]
)
def test_ewma_nan_handling_adjust_false_com1(kwargs):
    # GH#31178 / GH#66523: com=1 (and equivalent spellings) must use the
    # documented ignore_na=False recursion, not a convex combination over
    # the accumulated gap. alpha = 0.5, so
    # y[2] = ((1-alpha)**2 * 1 + alpha * 5) / ((1-alpha)**2 + alpha) = 11/3
    ser = Series([1.0, np.nan, 5.0])
    result = ser.ewm(adjust=False, ignore_na=False, **kwargs).mean()
    expected = Series([1.0, 1.0, 11 / 3])
    tm.assert_series_equal(result, expected)


def test_ewma_nan_handling_adjust_false_com1_continuous():
    # GH#31178 the com=1 result must be continuous in alpha
    ser = Series([1.0, np.nan, 5.0])
    result = ser.ewm(com=1, adjust=False, ignore_na=False).mean()
    for com in (1 - 1e-9, 1 + 1e-9):
        neighbor = ser.ewm(com=com, adjust=False, ignore_na=False).mean()
        tm.assert_series_equal(result, neighbor, atol=1e-8, rtol=0)


def test_ewma_times_equal_spacing_nan_adjust_false():
    # GH#66523 equally spaced times with a NaN gap use the documented
    # ignore_na=False recursion (11/3), not the convex-over-gap value 4.0.
    ser = Series([1.0, np.nan, 5.0])
    times = date_range("2000", freq="D", periods=3)
    result = ser.ewm(halflife="1D", times=times, adjust=False, ignore_na=False).mean()
    expected = ser.ewm(com=1, adjust=False, ignore_na=False).mean()
    tm.assert_series_equal(result, expected)
    tm.assert_almost_equal(result.iloc[2], 11 / 3)


def test_ewma_times_variable_spacing_with_nan():
    # GH#66523 elapsed time of the *current* step is convex; decay
    # accumulated across an ignore_na=False NaN row stays in old_wt and is
    # renormalized. times are 0, 1, 3 days with values [1, NaN, 5]:
    #   old_wt = 0.5 * 0.5**2 = 0.125
    #   new_wt = 1 - 0.5**2 = 0.75
    #   y = (0.125 * 1 + 0.75 * 5) / (0.125 + 0.75) = 31/7
    # The convex-over-accumulated-gap form would give 4.5; treating the
    # whole 0-to-3 day span as one normalized delta=3 step would give 4.2.
    times = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-04"])
    ser = Series([1.0, np.nan, 5.0])
    result = ser.ewm(halflife="1D", times=times, adjust=False, ignore_na=False).mean()
    expected = Series([1.0, 1.0, 31 / 7])
    tm.assert_series_equal(result, expected)


def test_ewma_ignore_na_true_is_skipped_row_step_not_dropna():
    # GH#66523 ignore_na=True skips the row-to-row step that lands on a
    # NaN. It is not "drop NaNs and keep the calendar gap".
    # Equally spaced [1, nan, 5] must stay 3.0 (one unit step after the
    # skip) to match no-times; dropna-and-keep-times would use a 2-day
    # gap and give 4.0, breaking that invariant.
    ser = Series([1.0, np.nan, 5.0])
    daily = date_range("2000", freq="D", periods=3)
    result = ser.ewm(halflife="1D", times=daily, adjust=False, ignore_na=True).mean()
    expected = ser.ewm(com=1, adjust=False, ignore_na=True).mean()
    tm.assert_series_equal(result, expected)
    tm.assert_almost_equal(result.iloc[2], 3.0)

    irregular = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-04"])
    keep = ser.ewm(halflife="1D", times=irregular, adjust=False, ignore_na=True).mean()
    # skip the NaN row, then apply only the last hop (2 days): 0.25*1 + 0.75*5
    tm.assert_series_equal(keep, Series([1.0, 1.0, 4.0]))
    dropped = (
        Series([1.0, 5.0])
        .ewm(
            halflife="1D",
            times=DatetimeIndex(["2000-01-01", "2000-01-04"]),
            adjust=False,
        )
        .mean()
    )
    # dropna-and-keep-times would use the full 3-day gap -> 4.5
    tm.assert_almost_equal(dropped.iloc[-1], 4.5)
    assert keep.iloc[-1] != dropped.iloc[-1]


def test_ewma_times_zero_length_interval_has_no_new_weight():
    # GH#66523 a zero-length step is 1 - (1-alpha)**0 = 0 new weight, so a
    # later row at the same timestamp does not move the mean. That is the
    # same rule with or without a preceding NaN row.
    times = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-02"])
    no_nan = (
        Series([1.0, 2.0, 9.0])
        .ewm(halflife="1D", times=times, adjust=False, ignore_na=False)
        .mean()
    )
    tm.assert_series_equal(no_nan, Series([1.0, 1.5, 1.5]))

    after_nan = (
        Series([1.0, np.nan, 9.0])
        .ewm(halflife="1D", times=times, adjust=False, ignore_na=False)
        .mean()
    )
    tm.assert_series_equal(after_nan, Series([1.0, 1.0, 1.0]))


def test_ewma_times_equal_spacing_other_step():
    # A uniform 2-day grid with halflife="2D" is one half-life per row,
    # same as no-times with numeric halflife=1.
    times = date_range("2000", freq="2D", periods=6)
    vals = [1.0, np.nan, 5.0, 6.0, np.nan, 9.0]
    ser = Series(vals)
    for adjust in (True, False):
        for ignore_na in (True, False):
            result = ser.ewm(
                halflife="2D", times=times, adjust=adjust, ignore_na=ignore_na
            ).mean()
            expected = ser.ewm(halflife=1.0, adjust=adjust, ignore_na=ignore_na).mean()
            tm.assert_series_equal(result, expected)


def test_ewma_times_long_nan_then_short_hop():
    # After a long ignore_na=False NaN gap, leftover old_wt is tiny, so
    # renormalizing a one-day new_wt still moves almost all the way to
    # the new observation. Convex-over-the-whole-gap would be even closer
    # to 5; the point is the leftover mass, not the last hop alone (which
    # would give 3.0).
    times = DatetimeIndex(["2000-01-01", "2000-01-11", "2000-01-12"])
    ser = Series([1.0, np.nan, 5.0])
    result = ser.ewm(halflife="1D", times=times, adjust=False, ignore_na=False).mean()
    old_wt = (0.5**10) * (0.5**1)
    new_wt = 1.0 - 0.5**1
    expected = Series([1.0, 1.0, (old_wt * 1.0 + new_wt * 5.0) / (old_wt + new_wt)])
    tm.assert_series_equal(result, expected)
    assert result.iloc[2] > 4.99


def test_ewma_groupby_times_uses_per_group_calendar():
    # Interleaved groups have a 2-day calendar gap. The times path must
    # use that gap (so it differs from no-times unit steps) and must
    # match ewm on each group alone.
    df = DataFrame(
        {"A": ["a", "b", "a", "b", "a", "b"], "B": [1.0, 1.0, np.nan, np.nan, 5.0, 5.0]}
    )
    times = date_range("2000", freq="D", periods=6)
    result = (
        df.groupby("A")
        .ewm(halflife="1D", times=times, adjust=False, ignore_na=False)
        .mean()
    )
    notimes = df.groupby("A").ewm(com=1, adjust=False, ignore_na=False).mean()
    # 2-day NaN then 2-day observation: (0.5**2 * 0.5**2 * 1 + (1-0.5**2)*5)
    # / (0.5**4 + (1-0.5**2)) = 4.6923..., not the unit-step 11/3.
    tm.assert_almost_equal(result.loc[("a", 4), "B"], 4.6923076923076925)
    assert not np.isclose(result.loc[("a", 4), "B"], notimes.loc[("a", 4), "B"])
    for key in ("a", "b"):
        mask = df["A"] == key
        direct = (
            Series(df.loc[mask, "B"].to_numpy())
            .ewm(halflife="1D", times=times[mask], adjust=False, ignore_na=False)
            .mean()
        )
        tm.assert_almost_equal(result.loc[key, "B"].to_numpy(), direct.to_numpy())


def _reference_ewm_mean(values, deltas, com, adjust, ignore_na):
    """Running-mass mean from the documented row-to-row split (GH#66523).

    Each walked row decays the mass by ``(1-alpha)**delta``. A new
    observation is mixed with weight ``1`` (adjust=True) or
    ``1 - (1-alpha)**delta`` (adjust=False, this row-to-row step only).
    ``ignore_na=True`` skips the step that lands on a NaN.
    """
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    if deltas is None:
        deltas = np.ones(max(n - 1, 0), dtype=np.float64)
    fade = 1.0 - 1.0 / (1.0 + com)
    out = np.full(n, np.nan)
    y = np.nan
    mass = 1.0
    for i in range(n):
        x = values[i]
        observed = x == x
        if y != y:
            if observed:
                y = x
                mass = 1.0
        elif observed or not ignore_na:
            step = deltas[i - 1] if i else 0.0
            decay = fade**step
            mass *= decay
            if observed:
                new_wt = 1.0 if adjust else (1.0 - decay)
                if y != x:
                    y = (mass * y + new_wt * x) / (mass + new_wt)
                mass = mass + new_wt if adjust else 1.0
        out[i] = y
    return out


@pytest.mark.parametrize("adjust", [True, False])
@pytest.mark.parametrize("ignore_na", [True, False])
@pytest.mark.parametrize("com", [0.0, 1.0, 2.5])
def test_ewm_mean_matches_row_step_reference(adjust, ignore_na, com):
    # Independent of aggregations.pyx: equal-spacing / no-times path.
    rng = np.random.default_rng(66523)
    values = rng.normal(size=12)
    values[rng.random(12) < 0.3] = np.nan
    result = Series(values).ewm(com=com, adjust=adjust, ignore_na=ignore_na).mean()
    expected = Series(_reference_ewm_mean(values, None, com, adjust, ignore_na))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("adjust", [True, False])
@pytest.mark.parametrize("ignore_na", [True, False])
def test_ewm_times_matches_row_step_reference(adjust, ignore_na):
    # Irregular gaps, a zero-length interval, and NaNs.
    rng = np.random.default_rng(665230)
    values = rng.normal(size=8)
    values[[1, 4]] = np.nan
    deltas = np.array([0.5, 1.0, 0.0, 2.0, 1.0, 3.0, 0.25])
    offsets = np.concatenate([[0.0], np.cumsum(deltas)])
    times = DatetimeIndex(
        np.datetime64("2000-01-01")
        + (offsets * 86_400 * 1_000_000_000).astype("timedelta64[ns]")
    )
    result = (
        Series(values)
        .ewm(halflife="1D", times=times, adjust=adjust, ignore_na=ignore_na)
        .mean()
    )
    expected = Series(_reference_ewm_mean(values, deltas, 1.0, adjust, ignore_na))
    tm.assert_series_equal(result, expected)


def test_ewm_alpha():
    # GH 10789
    arr = np.random.default_rng(2).standard_normal(100)
    locs = np.arange(20, 40)
    arr[locs] = np.nan

    s = Series(arr)
    a = s.ewm(alpha=0.61722699889169674).mean()
    b = s.ewm(com=0.62014947789973052).mean()
    c = s.ewm(span=2.240298955799461).mean()
    d = s.ewm(halflife=0.721792864318).mean()
    tm.assert_series_equal(a, b)
    tm.assert_series_equal(a, c)
    tm.assert_series_equal(a, d)


def test_ewm_domain_checks():
    # GH 12492
    arr = np.random.default_rng(2).standard_normal(100)
    locs = np.arange(20, 40)
    arr[locs] = np.nan

    s = Series(arr)
    msg = "comass must satisfy: comass >= 0"
    with pytest.raises(ValueError, match=msg):
        s.ewm(com=-0.1)
    s.ewm(com=0.0)
    s.ewm(com=0.1)

    msg = "span must satisfy: span >= 1"
    with pytest.raises(ValueError, match=msg):
        s.ewm(span=-0.1)
    with pytest.raises(ValueError, match=msg):
        s.ewm(span=0.0)
    with pytest.raises(ValueError, match=msg):
        s.ewm(span=0.9)
    s.ewm(span=1.0)
    s.ewm(span=1.1)

    msg = "halflife must satisfy: halflife > 0"
    with pytest.raises(ValueError, match=msg):
        s.ewm(halflife=-0.1)
    with pytest.raises(ValueError, match=msg):
        s.ewm(halflife=0.0)
    s.ewm(halflife=0.1)

    msg = "alpha must satisfy: 0 < alpha <= 1"
    with pytest.raises(ValueError, match=msg):
        s.ewm(alpha=-0.1)
    with pytest.raises(ValueError, match=msg):
        s.ewm(alpha=0.0)
    s.ewm(alpha=0.1)
    s.ewm(alpha=1.0)
    with pytest.raises(ValueError, match=msg):
        s.ewm(alpha=1.1)


@pytest.mark.parametrize("method", ["mean", "std", "var"])
def test_ew_empty_series(method):
    vals = Series([], dtype=np.float64)

    ewm = vals.ewm(3)
    result = getattr(ewm, method)()
    tm.assert_almost_equal(result, vals)


@pytest.mark.parametrize("min_periods", [0, 1])
@pytest.mark.parametrize("name", ["mean", "var", "std"])
def test_ew_min_periods(min_periods, name):
    # excluding NaNs correctly
    arr = np.random.default_rng(2).standard_normal(50)
    arr[:10] = np.nan
    arr[-10:] = np.nan
    s = Series(arr)

    # check min_periods
    # GH 7898
    result = getattr(s.ewm(com=50, min_periods=2), name)()
    assert result[:11].isna().all()
    assert not result[11:].isna().any()

    result = getattr(s.ewm(com=50, min_periods=min_periods), name)()
    if name == "mean":
        assert result[:10].isna().all()
        assert not result[10:].isna().any()
    else:
        # ewm.std, ewm.var (with bias=False) require at least
        # two values
        assert result[:11].isna().all()
        assert not result[11:].isna().any()

    # check series of length 0
    result = getattr(Series(dtype=object).ewm(com=50, min_periods=min_periods), name)()
    tm.assert_series_equal(result, Series(dtype="float64"))

    # check series of length 1
    result = getattr(Series([1.0]).ewm(50, min_periods=min_periods), name)()
    if name == "mean":
        tm.assert_series_equal(result, Series([1.0]))
    else:
        # ewm.std, ewm.var with bias=False require at least
        # two values
        tm.assert_series_equal(result, Series([np.nan]))

    # pass in ints
    result2 = getattr(Series(np.arange(50)).ewm(span=10), name)()
    assert result2.dtype == np.float64


@pytest.mark.parametrize("name", ["cov", "corr"])
def test_ewm_corr_cov(name):
    A = Series(np.random.default_rng(2).standard_normal(50), index=range(50))
    B = A[2:] + np.random.default_rng(2).standard_normal(48)

    A[:10] = np.nan
    B.iloc[-10:] = np.nan

    result = getattr(A.ewm(com=20, min_periods=5), name)(B)
    assert np.isnan(result.values[:14]).all()
    assert not np.isnan(result.values[14:]).any()


@pytest.mark.parametrize("min_periods", [0, 1, 2])
@pytest.mark.parametrize("name", ["cov", "corr"])
def test_ewm_corr_cov_min_periods(name, min_periods):
    # GH 7898
    A = Series(np.random.default_rng(2).standard_normal(50), index=range(50))
    B = A[2:] + np.random.default_rng(2).standard_normal(48)

    A[:10] = np.nan
    B.iloc[-10:] = np.nan

    result = getattr(A.ewm(com=20, min_periods=min_periods), name)(B)
    # binary functions (ewmcov, ewmcorr) with bias=False require at
    # least two values
    assert np.isnan(result.values[:11]).all()
    assert not np.isnan(result.values[11:]).any()

    # check series of length 0
    empty = Series([], dtype=np.float64)
    result = getattr(empty.ewm(com=50, min_periods=min_periods), name)(empty)
    tm.assert_series_equal(result, empty)

    # check series of length 1
    result = getattr(Series([1.0]).ewm(com=50, min_periods=min_periods), name)(
        Series([1.0])
    )
    tm.assert_series_equal(result, Series([np.nan]))


@pytest.mark.parametrize("name", ["cov", "corr"])
def test_different_input_array_raise_exception(name):
    A = Series(np.random.default_rng(2).standard_normal(50), index=range(50))
    A[:10] = np.nan

    msg = "other must be a DataFrame or Series"
    # exception raised is Exception
    with pytest.raises(ValueError, match=msg):
        getattr(A.ewm(com=20, min_periods=5), name)(
            np.random.default_rng(2).standard_normal(50)
        )


@pytest.mark.parametrize("name", ["var", "std", "mean"])
def test_ewma_series(series, name):
    series_result = getattr(series.ewm(com=10), name)()
    assert isinstance(series_result, Series)


@pytest.mark.parametrize("name", ["var", "std", "mean"])
def test_ewma_frame(frame, name):
    frame_result = getattr(frame.ewm(com=10), name)()
    assert isinstance(frame_result, DataFrame)


def test_ewma_span_com_args(series):
    A = series.ewm(com=9.5).mean()
    B = series.ewm(span=20).mean()
    tm.assert_almost_equal(A, B)
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, span=20)

    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        series.ewm().mean()


def test_ewma_halflife_arg(series):
    A = series.ewm(com=13.932726172912965).mean()
    B = series.ewm(halflife=10.0).mean()
    tm.assert_almost_equal(A, B)
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        series.ewm(span=20, halflife=50)
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, halflife=50)
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, span=20, halflife=50)
    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        series.ewm()


def test_ewm_alpha_arg(series):
    # GH 10789
    s = series
    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        s.ewm()

    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        s.ewm(com=10.0, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        s.ewm(span=10.0, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        s.ewm(halflife=10.0, alpha=0.5)


@pytest.mark.parametrize("func", ["cov", "corr"])
def test_ewm_pairwise_cov_corr(func, frame):
    result = getattr(frame.ewm(span=10, min_periods=5), func)()
    result = result.loc[(slice(None), 1), 5]
    result.index = result.index.droplevel(1)
    expected = getattr(frame[1].ewm(span=10, min_periods=5), func)(frame[5])
    tm.assert_series_equal(result, expected, check_names=False)


def test_numeric_only_frame(arithmetic_win_operators, numeric_only):
    # GH#46560
    kernel = arithmetic_win_operators
    df = DataFrame({"a": [1], "b": 2, "c": 3})
    df["c"] = df["c"].astype(object)
    ewm = df.ewm(span=2, min_periods=1)
    op = getattr(ewm, kernel, None)
    if op is not None:
        result = op(numeric_only=numeric_only)

        columns = ["a", "b"] if numeric_only else ["a", "b", "c"]
        expected = df[columns].agg([kernel]).reset_index(drop=True).astype(float)
        assert list(expected.columns) == columns

        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kernel", ["corr", "cov"])
@pytest.mark.parametrize("use_arg", [True, False])
def test_numeric_only_corr_cov_frame(kernel, numeric_only, use_arg):
    # GH#46560
    df = DataFrame({"a": [1, 2, 3], "b": 2, "c": 3})
    df["c"] = df["c"].astype(object)
    arg = (df,) if use_arg else ()
    ewm = df.ewm(span=2, min_periods=1)
    op = getattr(ewm, kernel)
    result = op(*arg, numeric_only=numeric_only)

    # Compare result to op using float dtypes, dropping c when numeric_only is True
    columns = ["a", "b"] if numeric_only else ["a", "b", "c"]
    df2 = df[columns].astype(float)
    arg2 = (df2,) if use_arg else ()
    ewm2 = df2.ewm(span=2, min_periods=1)
    op2 = getattr(ewm2, kernel)
    expected = op2(*arg2, numeric_only=numeric_only)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dtype", [int, object])
def test_numeric_only_series(arithmetic_win_operators, numeric_only, dtype):
    # GH#46560
    kernel = arithmetic_win_operators
    ser = Series([1], dtype=dtype)
    ewm = ser.ewm(span=2, min_periods=1)
    op = getattr(ewm, kernel, None)
    if op is None:
        # Nothing to test
        pytest.skip("No op to test")
    if numeric_only and dtype is object:
        msg = f"ExponentialMovingWindow.{kernel} does not implement numeric_only"
        with pytest.raises(NotImplementedError, match=msg):
            op(numeric_only=numeric_only)
    else:
        result = op(numeric_only=numeric_only)
        expected = ser.agg([kernel]).reset_index(drop=True).astype(float)
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("kernel", ["corr", "cov"])
@pytest.mark.parametrize("use_arg", [True, False])
@pytest.mark.parametrize("dtype", [int, object])
def test_numeric_only_corr_cov_series(kernel, use_arg, numeric_only, dtype):
    # GH#46560
    ser = Series([1, 2, 3], dtype=dtype)
    arg = (ser,) if use_arg else ()
    ewm = ser.ewm(span=2, min_periods=1)
    op = getattr(ewm, kernel)
    if numeric_only and dtype is object:
        msg = f"ExponentialMovingWindow.{kernel} does not implement numeric_only"
        with pytest.raises(NotImplementedError, match=msg):
            op(*arg, numeric_only=numeric_only)
    else:
        result = op(*arg, numeric_only=numeric_only)

        ser2 = ser.astype(float)
        arg2 = (ser2,) if use_arg else ()
        ewm2 = ser2.ewm(span=2, min_periods=1)
        op2 = getattr(ewm2, kernel)
        expected = op2(*arg2, numeric_only=numeric_only)
        tm.assert_series_equal(result, expected)
