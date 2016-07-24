import numpy as np
import pandas as pd


class algorithm(object):
    goal_time = 0.2

    def setup(self):
        N = 100000
        self.int = pd.Int64Index(np.arange(N).repeat(5))
        self.float = pd.Float64Index(np.random.randn(N).repeat(5))

    def time_int_factorize(self):
        self.int.factorize()

    def time_float_factorize(self):
        self.int.factorize()
