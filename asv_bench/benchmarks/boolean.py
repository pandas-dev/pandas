import numpy as np

import pandas as pd


class TimeLogicalOps:
    def setup(self):
        N = 10_000
        left, right, lmask, rmask = np.random.randint(0, 2, size=(4, N)).astype("bool")
        self.left = pd.arrays.BooleanArray(left, lmask)
        self.right = pd.arrays.BooleanArray(right, rmask)

    def time_or_scalar(self):
        self.left | True
        self.left | False

    def time_or_array(self):
        self.left | self.right

    def time_and_scalar(self):
        self.left & True
        self.left & False

    def time_and_array(self):
        self.left & self.right

    def time_xor_scalar(self):
        self.left ^ True
        self.left ^ False

    def time_xor_array(self):
        self.left ^ self.right
