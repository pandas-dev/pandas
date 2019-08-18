import pytest

from pandas import DataFrame, Timestamp
from pandas.util.testing import assert_frame_equal


@pytest.mark.parametrize(
    "data, expected_shift, expected_bfill, expected_ffill",
    [
        (
            {
                "id": ["A", "B", "A", "B", "A", "B"],
                "time": [
                    Timestamp("2019-01-01 12:00:00"),
                    Timestamp("2019-01-01 12:30:00"),
                    None,
                    None,
                    Timestamp("2019-01-01 14:00:00"),
                    Timestamp("2019-01-01 14:30:00"),
                ],
            },
            {
                "time": [
                    None,
                    None,
                    Timestamp("2019-01-01 12:00:00"),
                    Timestamp("2019-01-01 12:30:00"),
                    None,
                    None,
                ]
            },
            {
                "time": [
                    Timestamp("2019-01-01 12:00:00"),
                    Timestamp("2019-01-01 12:30:00"),
                    Timestamp("2019-01-01 14:00:00"),
                    Timestamp("2019-01-01 14:30:00"),
                    Timestamp("2019-01-01 14:00:00"),
                    Timestamp("2019-01-01 14:30:00"),
                ]
            },
            {
                "time": [
                    Timestamp("2019-01-01 12:00:00"),
                    Timestamp("2019-01-01 12:30:00"),
                    Timestamp("2019-01-01 12:00:00"),
                    Timestamp("2019-01-01 12:30:00"),
                    Timestamp("2019-01-01 14:00:00"),
                    Timestamp("2019-01-01 14:30:00"),
                ]
            },
        ),
        (
            {
                "id": ["A", "B", "A", "B", "A", "B"],
                "time": [
                    Timestamp("2019-01-01 12:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:30:00", tz="Asia/Tokyo"),
                    None,
                    None,
                    Timestamp("2019-01-01 14:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:30:00", tz="Asia/Tokyo"),
                ],
            },
            {
                "time": [
                    None,
                    None,
                    Timestamp("2019-01-01 12:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:30:00", tz="Asia/Tokyo"),
                    None,
                    None,
                ]
            },
            {
                "time": [
                    Timestamp("2019-01-01 12:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:30:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:30:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:30:00", tz="Asia/Tokyo"),
                ]
            },
            {
                "time": [
                    Timestamp("2019-01-01 12:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:30:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 12:30:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:00:00", tz="Asia/Tokyo"),
                    Timestamp("2019-01-01 14:30:00", tz="Asia/Tokyo"),
                ]
            },
        ),
    ],
)
def test_shift_bfill_ffill_tz(data, expected_shift, expected_bfill, expected_ffill):
    # GH27992: Check that timezone does not drop in shift, bfill, and ffill
    df = DataFrame(data)

    result = df.groupby("id").shift()
    expected = DataFrame(expected_shift)
    assert_frame_equal(result, expected)

    result = df.groupby("id").bfill()
    expected = DataFrame(expected_bfill)
    assert_frame_equal(result, expected)

    result = df.groupby("id").ffill()
    expected = DataFrame(expected_ffill)
    assert_frame_equal(result, expected)
