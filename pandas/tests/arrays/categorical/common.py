from pandas import Categorical


class TestCategorical:
    def setup_method(self):
        self.factor = Categorical(
            ["a", "b", "b", "a", "a", "c", "c", "c"], ordered=True
        )
