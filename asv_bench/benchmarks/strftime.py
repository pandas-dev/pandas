import numpy as np

import pandas as pd
from pandas import offsets


class DatetimeStrftime:
    timeout = 1500
    params = [1000, 10000]
    param_names = ["nobs"]

    def setup(self, nobs):
        d = "2018-11-29"
        dt = "2018-11-26 11:18:27.0"
        self.data = pd.DataFrame(
            {
                "dt": [np.datetime64(dt)] * nobs,
                "d": [np.datetime64(d)] * nobs,
                "r": [np.random.uniform()] * nobs,
            }
        )

    def time_frame_date_to_str(self, nobs):
        self.data["d"].astype(str)

    def time_frame_date_formatting_default(self, nobs):
        self.data["d"].dt.strftime(date_format=None)

    def time_frame_date_formatting_default_explicit(self, nobs):
        self.data["d"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_date_formatting_custom(self, nobs):
        self.data["d"].dt.strftime(date_format="%Y---%m---%d")

    def time_frame_datetime_to_str(self, nobs):
        self.data["dt"].astype(str)

    def time_frame_datetime_formatting_default(self, nobs):
        self.data["dt"].dt.strftime(date_format=None)

    def time_frame_datetime_formatting_default_explicit_date_only(self, nobs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_datetime_formatting_default_explicit(self, nobs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_datetime_formatting_default_with_float(self, nobs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_datetime_formatting_custom(self, nobs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")


class PeriodStrftime:
    timeout = 1500
    params = ([1000, 10000], ["D", "h"])
    param_names = ["nobs", "freq"]

    def setup(self, nobs, freq):
        self.data = pd.DataFrame(
            {
                "p": pd.period_range(start="2000-01-01", periods=nobs, freq=freq),
                "r": [np.random.uniform()] * nobs,
            }
        )
        self.data["i"] = self.data["p"]
        self.data.set_index("i", inplace=True)
        if freq == "D":
            self.default_fmt = "%Y-%m-%d"
        elif freq == "h":
            self.default_fmt = "%Y-%m-%d %H:00"

    def time_frame_period_to_str(self, nobs, freq):
        self.data["p"].astype(str)

    def time_frame_period_formatting_default(self, nobs, freq):
        self.data["p"].dt.strftime(date_format=None)

    def time_frame_period_formatting_default_explicit(self, nobs, freq):
        self.data["p"].dt.strftime(date_format=self.default_fmt)

    def time_frame_period_formatting_custom(self, nobs, freq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")

    def time_frame_period_formatting_iso8601_strftime_Z(self, nobs, freq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%SZ")

    def time_frame_period_formatting_iso8601_strftime_offset(self, nobs, freq):
        """Not optimized yet as %z is not supported by `convert_strftime_format`"""
        self.data["p"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%S%z")


class BusinessHourStrftime:
    timeout = 1500
    params = [1000, 10000]
    param_names = ["nobs"]

    def setup(self, nobs):
        self.data = pd.DataFrame(
            {
                "off": [offsets.BusinessHour()] * nobs,
            }
        )

    def time_frame_offset_str(self, nobs):
        self.data["off"].apply(str)

    def time_frame_offset_repr(self, nobs):
        self.data["off"].apply(repr)
