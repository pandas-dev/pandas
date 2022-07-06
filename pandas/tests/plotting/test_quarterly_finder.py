import numpy as np
import pytest
import pandas as pd
from pandas.plotting._matplotlib import converter


@pytest.mark.parametrize("year_span", [11.25, 30, 80, 150, 300, 584.5])
# valid ranges are from 11.25 to 584.5 years, limited at the bottom by
# if statements in the _quarterly_finder() function and at the top end
# by the Timestamp format
def test_quarterly_finder(year_span):
    # earliest start date given pd.Timestamp.min
    start_date = pd.to_datetime("1677Q4")
    daterange = pd.period_range(
        start_date, periods=year_span * 4, freq="Q").asi8
    vmin = daterange[0]
    vmax = daterange[-1]
    span = vmax - vmin + 1
    if span < 45:  # the quarterly finder is only invoked if the span is >= 45
        return
    nyears = span / 4
    (min_anndef, maj_anndef) = converter._get_default_annual_spacing(nyears)
    result = converter._quarterly_finder(vmin, vmax, "Q")
    quarters = pd.PeriodIndex(
        pd.arrays.PeriodArray(np.array([x[0] for x in result]), freq="Q")
    )
    majors = np.array([x[1] for x in result])
    minors = np.array([x[2] for x in result])
    major_quarters = quarters[majors]
    minor_quarters = quarters[minors]
    check_major_years = major_quarters.year % maj_anndef == 0
    check_minor_years = minor_quarters.year % min_anndef == 0
    check_major_quarters = major_quarters.quarter == 1
    check_minor_quarters = minor_quarters.quarter == 1
    assert (
        np.all(check_major_years)
        and np.all(check_minor_years)
        and np.all(check_major_quarters)
        and np.all(check_minor_quarters)
    )
