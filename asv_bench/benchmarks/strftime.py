import numpy as np

import pandas as pd
from pandas import offsets


class DatetimeStrftime:
    timeout = 1500
    params = [1000, 10000]
    param_names = ["obs"]

    def setup(self, obs):
        d = "2018-11-29"
        dt = "2018-11-26 11:18:27.0"
        self.data = pd.DataFrame(
            {
                "dt": [np.datetime64(dt)] * obs,
                "d": [np.datetime64(d)] * obs,
                "r": [np.random.uniform()] * obs,
            }
        )

    def time_frame_date_to_str(self, obs):
        self.data["d"].astype(str)

    def time_frame_date_formatting_default(self, obs):
        self.data["d"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_date_formatting_custom(self, obs):
        self.data["d"].dt.strftime(date_format="%Y---%m---%d")

    def time_frame_datetime_to_str(self, obs):
        self.data["dt"].astype(str)

    def time_frame_datetime_formatting_default_date_only(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_datetime_formatting_default(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_datetime_formatting_default_with_float(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_datetime_formatting_custom(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")


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
