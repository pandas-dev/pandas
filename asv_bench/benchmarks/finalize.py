import pandas as pd


class Finalize:
    param_names = ["series", "frame"]
    params = [pd.Series, pd.DataFrame]

    def setup(self, param):
        N = 1000
        df = param(dtype=float)
        for i in range(N):
            df.attrs[i] = i
        self.df = df

    def time_finalize_micro(self, param):
        self.df.__finalize__(self.df, method="__finalize__")
