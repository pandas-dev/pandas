import pytest

import pandas.tseries.offsets as offsets

DATE_OFFSETS = {
    "Day",
    "MonthBegin",
    "MonthEnd",
    "SemiMonthBegin",
    "SemiMonthEnd",
    "YearBegin",
    "YearEnd",
    "QuarterBegin",
    "QuarterEnd",
    "LastWeekOfMonth",
    "Week",
    "WeekOfMonth",
    "Easter",
    "Hour",
    "Minute",
    "Second",
    "Milli",
    "Micro",
    "Nano",
    "DateOffset",
}

# BUSINESS_OFFSETS = {
#     "BusinessDay",
#     "BDay",
#     "CustomBusinessDay",
#     "CBMonthBegin",
#     "CBMonthEnd",
#     "BMonthBegin",
#     "BMonthEnd",
#     "BusinessHour",
#     "CustomBusinessHour",
#     "BYearBegin",
#     "BYearEnd",
#     "CDay",
#     "BQuarterBegin",
#     "BQuarterEnd",
#     "FY5253Quarter",
#     "FY5253",
# }


@pytest.fixture(params=[getattr(offsets, o) for o in DATE_OFFSETS])
def date_offset_types(request):
    """
    Fixture for all the datetime offsets available for a time series.
    """
    return request.param


# @pytest.fixture(params=[getattr(offsets, o) for o in BUSINESS_OFFSETS])
# def business_offset_types(request):
#     """
#     Fixture for all the datetime offsets available for a time series.
#     """
#     return request.param


@pytest.fixture(
    params=[
        getattr(offsets, o)
        for o in DATE_OFFSETS
        if issubclass(getattr(offsets, o), offsets.MonthOffset) and o != "MonthOffset"
    ]
)
def date_month_classes(request):
    """
    Fixture for month based datetime offsets available for a time series.
    """
    return request.param


# @pytest.fixture(
#     params=[
#         getattr(offsets, o)
#         for o in BUSINESS_OFFSETS
#         if issubclass(getattr(offsets, o), offsets.MonthOffset) and o != "MonthOffset"
#     ]
# )
# def business_month_classes(request):
#     """
#     Fixture for month based datetime offsets available for a time series.
#     """
#     return request.param