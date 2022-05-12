from inspect import getmembers, ismethod

import numpy as np
import pandas as pd

from pandas import offsets


class DatetimeStrftime:
    fname = "__test__.csv"
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

    def time_frame_date_no_formatting(self, obs):
        self.data["d"].astype(str)

    def time_frame_date_formatting(self, obs):
        self.data["d"].dt.strftime(date_format="%Y-%m-%d")

    def time_frame_datetime_no_formatting(self, obs):
        self.data["dt"].astype(str)

    def time_frame_datetime_formatting(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S")

    def time_frame_datetime_formatting_with_float(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d %H:%M:%S.%f")

    def time_frame_datetime_formatting_date_only(self, obs):
        self.data["dt"].dt.strftime(date_format="%Y-%m-%d")


class BusinessHourStrftime:
    fname = "__test__.csv"
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


if __name__ == '__main__':
    # To debug this ASV benchmark, simply run this file with python
    from itertools import product
    for cls in (DatetimeStrftime, BusinessHourStrftime):
        if len(cls.param_names) == 1:
            all_params = [{cls.param_names[0]: p} for p in cls.params]
        else:
            all_params = [{n: p for n, p in zip(cls.param_names, ps)}
                          for ps in product(*cls.params)]
        for kwargs in all_params:
            kwargs_str = ','.join([f'{k}={v}' for k, v in kwargs.items()])
            print(f"Executing {cls} with {kwargs_str}")
            o = cls()
            o.setup(**kwargs)
            for k, v in getmembers(o, predicate=ismethod):
                if k != "setup":
                    print(f"Executing {v.__name__}({kwargs_str})")
                    v(**kwargs)
                    print(f"Executing {v.__name__}({kwargs_str}): DONE")
            print(f"Executing {cls} with {kwargs_str}: DONE")
