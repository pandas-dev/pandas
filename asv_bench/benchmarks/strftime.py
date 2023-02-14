import numpy as np

import pandas as pd
from pandas import offsets


class DatetimeStrftime:
    timeout = 1500
    params = ([1000, 10000], [False, True])
    param_names = ["obs", "tz_aware"]

    def setup(self, obs, tz_aware):
        d = "2018-11-29"
        dt = "2018-11-26 11:18:27.0"
        self.data = pd.DataFrame(
            {
                "dt": [np.datetime64(dt)] * obs,
                "d": [np.datetime64(d)] * obs,
                "r": [np.random.uniform()] * obs,
            }
        )
        if tz_aware:
            self.data["dt"] = self.data["dt"].dt.tz_localize("UTC")
            self.data["d"] = self.data["d"].dt.tz_localize("UTC")

    def time_frame_date_to_str(self, obs, tz_aware):
        self.data["d"].astype(str)

    def time_frame_date_formatting_default(self, obs, tz_aware):
        self.data["d"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_date_formatting_index_default(self, obs, tz_aware):
        self.data.set_index("d").index.format()

    def time_frame_date_formatting_custom(self, obs, tz_aware):
        self.data["d"].dt.strftime(date_format="%Y---%m---%d")

    def time_frame_date_formatting_index_custom(self, obs, tz_aware):
        self.data.set_index("d").index.format(date_format="%Y---%m---%d")

    def time_frame_datetime_to_str(self, obs, tz_aware):
        self.data["dt"].astype(str)

    def time_frame_datetime_formatting_default_date_only(self, obs, tz_aware):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_datetime_formatting_default(self, obs, tz_aware):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_datetime_formatting_default_with_float(self, obs, tz_aware):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_datetime_formatting_index_default(self, obs, tz_aware):
        self.data.set_index("dt").index.format()

    def time_frame_datetime_formatting_custom(self, obs, tz_aware):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")

    def time_frame_datetime_formatting_index_custom(self, obs, tz_aware):
        self.data.set_index("dt").index.format(date_format="%Y-%m-%d --- %H:%M:%S")

    def time_frame_datetime_formatting_iso8601_map(self, obs, tz_aware):
        self.data["dt"].map(lambda timestamp: timestamp.isoformat())

    def time_frame_datetime_formatting_iso8601_strftime_Z(self, obs, tz_aware):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%SZ")

    def time_frame_datetime_formatting_iso8601_strftime_offset(self, obs, tz_aware):
        """Not optimized yet as %z is not supported by `convert_strftime_format`"""
        self.data["dt"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%S%z")

    # def time_frame_datetime_formatting_iso8601_isoformat(self, obs, tz_aware):
    #     TODO this PR is probably a good opportunity to add this too, or maybe
    #      another PR
    #     self.data["dt"].dt.isoformat()


class PeriodStrftime:
    timeout = 1500
    params = ([1000, 10000], ["D", "H"])
    param_names = ["obs", "fq"]

    def setup(self, obs, fq):
        self.data = pd.DataFrame(
            {
                "p": pd.period_range(start="2000-01-01", periods=obs, freq=fq),
                "r": [np.random.uniform()] * obs,
            }
        )

    def time_frame_period_to_str(self, obs, fq):
        self.data["p"].astype(str)

    def time_frame_period_formatting_default_date_only(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_period_formatting_default(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_period_formatting_index_default(self, obs, fq):
        self.data.set_index("p").index.format()

    def time_frame_period_formatting_default_with_float(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_period_formatting_custom(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")

    # def time_frame_period_formatting_iso8601_map(self, obs, fq):
    #     self.data["p"].map(lambda p: p.isoformat())

    def time_frame_period_formatting_iso8601_strftime_Z(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%SZ")

    def time_frame_period_formatting_iso8601_strftime_offset(self, obs, fq):
        """Not optimized yet as %z is not supported by `convert_strftime_format`"""
        self.data["p"].dt.strftime(date_format="%Y-%m-%dT%H:%M:%S%z")


class BusinessHourStrftime:
    timeout = 1500
    params = [1000, 10000]
    param_names = ["obs"]

    def setup(self, obs):
        self.data = pd.DataFrame(
            {
                "off": [offsets.BusinessHour()] * obs,
            }
        )

    def time_frame_offset_str(self, obs):
        self.data["off"].apply(str)

    def time_frame_offset_repr(self, obs):
        self.data["off"].apply(repr)
