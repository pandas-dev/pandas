import pandas as pd


class MyDtype(pd.api.extensions.ExtensionDtype):
    name = "mydtype"
    type = list


class MyEA(pd.api.extensions.ExtensionArray):
    def __init__(self, data):
        self.data = data
        self._dtype = MyDtype()

    @property
    def dtype(self):
        return self._dtype

    def __len__(self):
        return 1

    def __array__(self, dtype=None):
        raise ValueError("Cannot be converted to an array!")

    def _format_array(
        self,
        formatter: None,
        float_format: None,
        na_rep="NaN",
        digits=None,
        space=None,
        justify="right",
        decimal=".",
        leading_space=True,
        quoting=None,
    ):
        return ["<MyEA>([1])"]


def test_no_conversion():
    s = pd.Series(MyEA([1]))
    repr(s)  # OK!

    df = pd.DataFrame({"A": MyEA([1])}, copy=False)
    repr(df)  # OK!
