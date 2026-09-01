import pandas as pd
import pandas._testing as tm


class SubclassedCategorical(pd.Categorical):
    pass


class TestCategoricalSubclassing:
    def test_constructor(self):
        sc = SubclassedCategorical(["a", "b", "c"])
        assert isinstance(sc, SubclassedCategorical)
        tm.assert_categorical_equal(sc, pd.Categorical(["a", "b", "c"]))

    def test_from_codes(self):
        sc = SubclassedCategorical.from_codes([1, 0, 2], ["a", "b", "c"])
        assert isinstance(sc, SubclassedCategorical)
        exp = pd.Categorical.from_codes([1, 0, 2], ["a", "b", "c"])
        tm.assert_categorical_equal(sc, exp)

    def test_map(self):
        sc = SubclassedCategorical(["a", "b", "c"])
        res = sc.map(lambda x: x.upper(), na_action=None)
        assert isinstance(res, SubclassedCategorical)
        exp = pd.Categorical(["A", "B", "C"])
        tm.assert_categorical_equal(res, exp)
