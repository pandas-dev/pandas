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

    def time_frame_datetime_formatting_default(self, obs):
        self.data["dt"].dt.strftime(date_format=None)

    def time_frame_datetime_formatting_default_explicit_date_only(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_datetime_formatting_default_explicit(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_datetime_formatting_default_with_float(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_datetime_formatting_custom(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")


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

    def time_frame_period_formatting_default(self, obs, fq):
        """Note that as opposed to datetimes, the default format of periods are
        many and depend from the period characteristics, so we have almost no chance
        to reach the same level of performance if a 'default' format string is
        explicitly provided by the user. See
        time_frame_datetime_formatting_default_explicit above."""
        self.data["p"].dt.strftime(date_format=None)

    def time_frame_period_formatting_index_default(self, obs, fq):
        self.data.set_index("p").index.format()

    def time_frame_period_formatting_custom(self, obs, fq):
        self.data["p"].dt.strftime(date_format="%Y-%m-%d --- %H:%M:%S")

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


if __name__ == "__main__":
    # A __main__ to easily debug this script
    for cls in (DatetimeStrftime, PeriodStrftime, BusinessHourStrftime):
        all_params = dict()
        all_p_values = cls.params
        if len(cls.param_names) == 1:
            all_p_values = (all_p_values,)
        for p_name, p_values in zip(cls.param_names, all_p_values):
            all_params[p_name] = p_values

        from itertools import product

        for case in product(*all_params.values()):
            p_dict = {p_name: p_val for p_name, p_val in zip(all_params.keys(), case)}
            print(f"{cls.__name__} - {p_dict}")
            o = cls()
            o.setup(**p_dict)
            for m_name, m in cls.__dict__.items():
                if callable(m):
                    print(m_name)
                    m(o, **p_dict)
